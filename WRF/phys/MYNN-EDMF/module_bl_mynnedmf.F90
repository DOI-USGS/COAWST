!>\file module_bl_mynnedmf.F90
!! This file contains the entity of MYNN-EDMF PBL scheme.
! **********************************************************************
! *   An improved Mellor-Yamada turbulence closure model               *
! *                                                                    *
! *      Original author: M. Nakanishi (N.D.A), naka@nda.ac.jp         *
! *      Translated into F90 and implemented in WRF-ARW by:            *
! *                       Mariusz Pagowski (NOAA-GSL)                  *
! *      Subsequently developed by:                                    *
! *                 Joseph Olson, Jaymes Kenyon (NOAA/GSL),            *
! *                 Wayne Angevine (NOAA/CSL), Kay Suselj (NASA/JPL),  *
! *                 Franciano Puhales (UFSM), Laura Fowler (NCAR),     *
! *                 Elynn Wu (UCSD), and Jordan Schnell (NOAA/GSL),    *
! *                 Tiziani Cherubin (MKWC, Univ. Hawaii),             *
! *                 Xia Sun (NOAA/GSL/CIRES)                           *
! *                                                                    *
! *   Contents:                                                        *
! *                                                                    *
! *   mynnedmf - the main subroutine which calls all other routines    *
! *   --------------                                                   *
! *     1. mym_initialize  (to be called once initially)               *
! *        gives the closure constants and initializes the turbulent   *
! *        quantities.                                                 *
! *     2. get_pblh                                                    *
! *        Calculates the boundary layer height                        *
! *     3. scale_aware                                                 *
! *        Calculates scale-adaptive tapering functions                *
! *     4. mym_condensation                                            *
! *        determines the liquid water content and the cloud fraction  *
! *        diagnostically.                                             *
! *     5. dmp_mf                                                      *
! *        Calls the (nonlocal) mass-flux component                    *
! *     6. ddmp_mf                                                     *
! *        Calls the downdraft mass-flux component                     *
! *    (-) mym_level2      (called in the other subroutines)           *
! *        calculates the stability functions at Level 2.              *
! *    (-) mym_length      (called in the other subroutines)           *
! *        calculates the master length scale.                         *
! *     7. mym_turbulence                                              *
! *        calculates the vertical diffusivity coefficients and the    *
! *        production terms for the turbulent quantities.              *
! *     8. mym_predict                                                 *
! *        predicts the turbulent quantities at the next step.         *
! *                                                                    *
! *             call mym_initialize                                    *
! *                  |                                                 *
! *                  |<----------------+                               *
! *                  |                 |                               *
! *             call get_pblh          |                               *
! *             call scale_aware       |                               *
! *             call mym_condensation  |                               *
! *             call dmp_mf            |                               *
! *             call ddmp_mf           |                               *
! *             call mym_turbulence    |                               *
! *             call mym_predict       |                               *
! *                  |                 |                               *
! *                  |-----------------+                               *
! *                  |                                                 *
! *                 end                                                *
! *                                                                    *
! *   Variables worthy of special mention:                             *
! *     tref   : Reference temperature                                 *
! *     thl    : Liquid water potential temperature                    *
! *     qw     : Total water (water vapor+liquid water) content        *
! *     ql     : Liquid water content                                  *
! *     vt, vq : Functions for computing the buoyancy flux             *
! *     qke    : 2 * TKE                                               *
! *     el     : mixing length                                         *
! *                                                                    *
! *     If the water contents are unnecessary, e.g., in the case of    *
! *     ocean models, thl is the potential temperature and qw, ql, vt  *
! *     and vq are all zero.                                           *
! *                                                                    *
! *   Grid arrangement:                                                *
! *             k+1 +---------+                                        *
! *                 |         |     i = 1 - nx                         *
! *             (k) |    *    |     k = 1 - nz                         *
! *                 |         |                                        *
! *              k  +---------+                                        *
! *                 i   (i)  i+1                                       *
! *                                                                    *
! *     All the predicted variables are defined at the center (*) of   *
! *     the grid boxes. The diffusivity coefficients and two of their  *
! *     components (el and stability functions sh & sm) are, however,  *
! *     defined on the walls of the grid boxes.                        *
! *     # Upper boundary values are given at k=nz.                     *
! *                                                                    *
! *   References:                                                      *
! *     1. Nakanishi, M., 2001:                                        *
! *        Boundary-Layer Meteor., 99, 349-378.                        *
! *     2. Nakanishi, M. and H. Niino, 2004:                           *
! *        Boundary-Layer Meteor., 112, 1-31.                          *
! *     3. Nakanishi, M. and H. Niino, 2006:                           *
! *        Boundary-Layer Meteor., 119, 397-407.                       *
! *     4. Nakanishi, M. and H. Niino, 2009:                           *
! *        Jour. Meteor. Soc. Japan, 87, 895-912.                      *
! *     5. Olson J. and coauthors, 2019: A description of the          *
! *        MYNN-EDMF scheme and coupling to other components in        *
! *        WRF-ARW. NOAA Tech. Memo. OAR GSD, 61, 37 pp.,              *
! *        https://doi.org/10.25923/n9wm-be49.                         * 
! *     6. Puhales, Franciano S. and coauthors, 2020: Turbulent        *
! *        Kinetic Energy Budget for MYNN-EDMF PBL Scheme in WRF model.*
! *        Universidade Federal de Santa Maria Technical Note. 9 pp.   *
! *        https://repositorio.ufsm.br/handle/1/28327                  *
! *     7. Olson, J. B., W. M. Angevine, D. D. Turner, X. Sun,         *
! *        J. M. Simonson, C. Evans, J. S. Kenyon, H. Li, J. Schnell,  *
! *        F. S. Puhales, T. Cherubini, W. Li, and M. Zhang, 2026:     *
! *        A Description of the MYNN-EDMF Turbulence Scheme. NOAA Tech.*
! *        Memo. OAR GSL-77. 60 pp. https://doi.org/10.25923/rahr-sj70 *
!***********************************************************************
! notes on the versioning of the myhnn-edmf:
!
! Version 1.0.0: Described by references 1-4 (above). Implemented into WRF
!    WRFv3.0-v3.4.1 for HRRRv1/RAPv2.
!
! Version 2.0.0: Still best described by references 1-4 (above). Approximately pplies to
!    WRFv3.5-v3.6.1 and used in HRRRv2/RAPv3.
!
! Version 3.0.0:  Still best described by references 1-4 (above), despite significant
!    code departures from the original scheme. This version pproximately pplies to
!    WRFv3.7-3.9 and used for HRRRv3/RAPv4.
!
! Version 4.0.0: Best described by reference 5 (above). Approximately applies to
!    WRFv4.0-v4.4 (but not well synced during this period) and used for HRRRv4/RAPv5.
!
! Version 5.0.0: Mods stemming from both RRFSv1 development and community development in WRF.
!    Approximately applies to WRFv4.5-v4.7 / CCPP
!
! Version 6.0.0: Mods stemming from both RRFSv2+ development and community development in WRF.
!     This version marks the begining of the MYNN-EDMF submodule repository, where all
!     development is now tracked in a public-facing github repository (https://github.com/NCAR/MYNN-EDMF.git).
!     WRFv4.8 / CCPP / MPAS. Many updates captured in the Olson et al. (2026) MYNN-EDMF
!     tech note (see listed references above).
!====================================================================

module module_bl_mynnedmf

 implicit none

!--------------------------------------------------------------------
! global debugging configuration options:
!--------------------------------------------------------------------
 logical,parameter:: debug_code = .false.
 integer,parameter:: idbg       = 452 !specific i-point to write out
 integer,parameter:: jdbg       = 272 !specific j-point to write out

contains

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine is the MYNN-EDMF PBL driver routine, which
!! encompassed the majority of the subroutines that comprise the 
!! procedures that ultimately solve for tendencies of 
!! \f$U, V, \theta, q_v, q_c, and q_i\f$.
!!\section gen_mynnedmf_driver mynnedmf_driver General Algorithm
!> @{
  subroutine mynnedmf(            i                 , j                 , &
             initflag           , restart           , cycling           , &
             delt               , dz1               , dx                , &
             !3d state variables
             u1                 , v1                , w1                , &
             th1                , sqv1              , sqc1              , &
             sqi1               , sqs1              , qnc1              , &
             qni1               , qnwfa1            , qnifa1            , &
             qnbca1             , ozone1            , pres1             , &
             ex1                , rho1              , tk1               , &
             !2d surface fields
             xland              , ts                , qsfc              , &
             ps                 , ust               , ch                , &
             hfx                , qfx               , znt               , &
             wspd               , uoce              , voce              , &
             !eddy-diffusivity
             qke1               , qke_adv1          , el1               , &
             sh1                , sm1               , kh1               , &
             km1                ,                                         &
             !smoke variables
             nchem              , ndvel             ,                     &
             chem1              , vdep              , frp               , &
             emis_ant_no        , mix_chem,enh_mix  , settle1           , &
             !generic scalar array
             scalars            , nscalars          ,                     &
             !higher-order moments
             tsq1               , qsq1              , cov1              , &
             !tendencies
             du1                , dv1               , dth1              , &
             dqv1               , dqc1              , dqi1              , &
             dqnc1              , dqni1             , dqs1              , &
             dqnwfa1            , dqnifa1           , dqnbca1           , &
             dozone1            , rthraten1         ,                     &
             !2d output
             pblh               , kpbl              , maxwidth_dd       , &
             maxwidth           , maxmf             , ztop_plume        , &
             excess_h           , excess_q          , maxmf_dd          , &
             maxtkeprod         , cldtop_cooling    , ent_eff           , &
             !tke budget arrays
             dqke1              , qwt1              , qshear1           , &
             qbuoy1             , qdiss1            ,                     &
             !subgrid clouds
             qc_bl1             , qi_bl1            , cldfra_bl1        , &
             !namelist configurations option
             bl_mynn_tkeadvect  , tke_budget        , bl_mynn_cloudpdf  , &
             bl_mynn_mixlength  , bl_mynn_closure   , bl_mynn_ess       , &
             bl_mynn_edmf       , bl_mynn_edmf_mom  , bl_mynn_edmf_tke  , &
             bl_mynn_mixscalars , bl_mynn_mixaerosols,bl_mynn_mixnumcon , &
             bl_mynn_output     , bl_mynn_cloudmix  , bl_mynn_mixqt     , &
             bl_mynn_edmf_dd    ,                                         &
             !3d emdf output
             edmf_a1            , edmf_w1           , edmf_qt1          , &
             edmf_thl1          , edmf_ent1         , edmf_qc1          , &
             sub_thl1           , sub_sqv1          , det_thl1          , &
             det_sqv1           ,                                         &
             !spp
             spp_pbl            , pattern_spp_pbl1  ,                     &
             !check flags
             FLAG_QC            , FLAG_QI           , FLAG_QNC          , &
             FLAG_QNI           , FLAG_QS           , FLAG_QNWFA        , &
             FLAG_QNIFA         , FLAG_QNBCA        , FLAG_OZONE        , &
             KTS                , KTE               , errmsg, errflg      )

!-------------------------------------------------------------------
 use module_bl_mynnedmf_common, only: cp,xlvcp,xlscp,p608,p1000mb,rcp,&
      karman,gtr,cphh_st,cphm_unst,cphh_unst,b1,zero,one,p01,p5,      &
      five,kind_phys
    
 integer, intent(in) :: initflag
 logical, intent(in) :: restart,cycling
 integer, intent(in) :: tke_budget
 integer, intent(in) :: bl_mynn_cloudpdf
 integer, intent(in) :: bl_mynn_mixlength
 integer, intent(in) :: bl_mynn_edmf
 integer, intent(in) :: bl_mynn_edmf_dd
 logical, intent(in) :: bl_mynn_tkeadvect
 integer, intent(in) :: bl_mynn_edmf_mom
 integer, intent(in) :: bl_mynn_edmf_tke
 integer, intent(in) :: bl_mynn_mixscalars
 integer, intent(in) :: bl_mynn_mixaerosols
 integer, intent(in) :: bl_mynn_mixnumcon
 integer, intent(in) :: bl_mynn_output
 integer, intent(in) :: bl_mynn_cloudmix
 integer, intent(in) :: bl_mynn_mixqt
 integer, intent(in) :: bl_mynn_ess
 real(kind_phys), intent(in) :: bl_mynn_closure

 logical, intent(in) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC,&
                        FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA, &
                        FLAG_OZONE,FLAG_QS

 logical, intent(in) :: mix_chem,enh_mix

 integer, intent(in) :: KTS,KTE

 character(len=*),intent(out):: &
    errmsg        ! output error message (-).
 integer,intent(out):: &
    errflg        ! output error flag (-).

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

!scalar variables
 real(kind_phys), intent(in)    :: delt
 real(kind_phys), intent(in)    :: dx
 real(kind_phys), intent(in)    :: ust,ch,qsfc,ps,wspd,xland
 real(kind_phys), intent(in)    :: ts,znt,hfx,qfx,uoce,voce
 real(kind_phys), intent(inout) :: pblh,excess_h,excess_q
 real(kind_phys), intent(inout) :: maxmf,maxwidth,ztop_plume
 real(kind_phys), intent(inout) :: maxmf_dd,maxwidth_dd,maxtkeprod,cldtop_cooling,ent_eff
 integer,         intent(in)    :: i,j
 integer,         intent(inout) :: kpbl
!local
 real(kind_phys) :: psig_bl,psig_shcu,rmol,rmolh
 integer :: imd,jmd
 integer :: k,kproblem

!column variables (all end with a "1")
 real(kind_phys), dimension(kts:kte), intent(in)      ::            &
       dz1,u1,v1,w1,th1,pres1,ex1,rho1,tk1,rthraten1
 real(kind_phys), dimension(kts:kte), intent(inout)   ::            &
       sqv1,sqc1,sqi1,sqs1,qni1,qnc1,qnwfa1,qnifa1,qnbca1,ozone1,   &
       qke1,tsq1,qsq1,cov1,qke_adv1,                                &
       sh1,sm1,el1,                                                 & !interface, but kte+1 not included
       du1,dv1,dth1,dqv1,dqc1,dqi1,dqs1,                            &
       dqni1,dqnc1,dqnwfa1,dqnifa1,dqnbca1,dozone1,                 &
       qc_bl1,qi_bl1,cldfra_bl1,edmf_a1,edmf_w1,                    &
       edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1,                       &
       sub_thl1,sub_sqv1,det_thl1,det_sqv1
 real(kind_phys), dimension(kts:kte), intent(inout), optional ::    &
       qwt1,qshear1,qbuoy1,qdiss1,dqke1
 real(kind_phys), dimension(kts:kte), intent(out)     ::            & !interface
       kh1,km1
!local
 real(kind_phys), dimension(kts:kte)                  ::            &
       qc_bl1_old,qi_bl1_old,cldfra_bl1_old,dummy1,dummy2,          &
       diss_heat1,                                                  &
       thl1,thv1,thlv1,qv1,qc1,qi1,qs1,sqw1,                        &
       thl_tot1,qc_tot1,qi_tot1,                                    &
       dfm1, dfh1, dfq1, tcd1, qcd1,                                &
       pdk1, pdt1, pdq1, pdc1,                                      &
       vt1, vq1, sgm1, kzero1
 real(kind_phys), dimension(kts:kte)                  ::            &
       edmf_u1,edmf_v1,edmf_qv1
 
!smoke/chemical arrays
 integer, intent(in) ::   nchem, ndvel
 real(kind_phys), dimension(kts:kte,nchem), intent(inout) :: chem1
 real(kind_phys), dimension(kts:kte,nchem), intent(in   ) :: settle1
 real(kind_phys), dimension(ndvel), intent(in)    :: vdep
 real(kind_phys),                   intent(in)    :: frp,emis_ant_no
 real(kind_phys), dimension(kts:kte+1,nchem)      :: s_awchem1
 integer :: ic
!scalar array
 integer, intent(in) ::   nscalars 
 real(kind_phys), dimension(kts:kte,nscalars), intent(inout) :: scalars
 real(kind_phys), dimension(kts:kte+1,nscalars)      :: s_awscalars1

!local vars
!mass-flux variables
 real(kind_phys), dimension(kts:kte) ::                             &
       dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf
 real(kind_phys), dimension(kts:kte) ::                             &
       edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1,edmf_thl_dd1,              &
       edmf_ent_dd1,edmf_qc_dd1
 real(kind_phys), dimension(kts:kte) ::                             &
       sub_u1,sub_v1,det_sqc1,det_u1,det_v1
 real(kind_phys), dimension(kts:kte+1) ::                           & !interface
       s_aw1,s_awthl1,s_awqt1,                                      &
       s_awqv1,s_awqc1,s_awu1,s_awv1,s_awqke1,s_awqsq1,             &
       s_awqnc1,s_awqni1,s_awqnwfa1,s_awqnifa1,                     &
       s_awqnbca1
 real(kind_phys), dimension(kts:kte+1) ::                           & !interface
       sd_aw1,sd_awthl1,sd_awqt1,                                   &
       sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,sd_awqke1,                 &
       sd_awqi1,sd_awqnc1,sd_awqni1,                                &
       sd_awqnwfa1,sd_awqnifa1
 integer:: ktop_plume
 real(kind_phys), dimension(kts:kte+1) :: zw1              !interface
!mass flux (or analytical) tke production
 real(kind_phys), dimension(kts:kte) :: TKEprod_dn,TKEprod_up

!-------------------------------------------------- 
!-------- configuration options:
!--------------------------------------------------
!Option to switch flux-profile relationship for surface (from Puhales et al. 2020)
!0: use original Dyer-Hicks, 1: use Cheng-Brustaert and Blended COARE
 integer,parameter:: bl_mynn_stfunc = 0

!>Use Canuto/Kitamura mod (remove Ric and negative TKE) (real; 1.0:yes, 0.0:no)
!!For more info, see Canuto et al. (2008 JAS) and Kitamura (Journal of the
!!Meteorological Society of Japan, Vol. 88, No. 5, pp. 857-864, 2010).
!!Note that this change required further modification of other parameters
!!above (c2, c3). If you want to remove this option, set c2 and c3 constants
!!(above) back to NN2009 values (see commented out lines next to the
!!parameters above). This only removes the negative TKE problem
!!but does not necessarily improve performance - neutral impact. Also,
!!it is a real variable because it is used within an equation to
!!activate/deactive this option
 real(kind_phys),parameter:: CKmod   = 1.0

!>Parameter used for the conversion of tsq to qpe when using the TTE closure
 real(kind_phys),parameter:: taue    = 80.

!>Option to control the upwind/centerd finite differencing of explicit mass-flux solver when
!>bl_mynn_edmf = 2  (0: mass flux inactive, 1: implicit, 2: explicit)
 real(kind_phys),parameter:: upwind  = 1.0 ! upwind=1.0: use upwind approximation for mass-flux calculation
                                           ! upwind=0.5: use centered difference for mass-flux calculation
                                           ! explicit mass-flux can use either upwind or centered-difference
                                           ! implicit mass-flux only uses the centered differencing method.
!---------------------------------------------------
 
 real(kind_phys) :: cpm,sqcg,flt,fltv,flq,flqv,flqc,                &
       pmz,phh,exnerg,zet,phi_m,                                    &
       afk,abk,th_sfc,wsp

 logical :: INITIALIZE_QKE,problem

!stochastic fields
 integer,  intent(in)                             :: spp_pbl
 real(kind_phys), dimension(kts:kte), intent(in)  :: pattern_spp_pbl1

!substepping TKE
 integer :: nsub
 real(kind_phys) :: delt2

    errmsg = " "
    errflg = 0

    if (debug_code) then !check incoming values
       problem = .false.
       do k=kts,kte
          wsp  = sqrt(u1(k)**2 + v1(k)**2)
          if (abs(hfx) > 1200. .or. abs(qfx) > 0.001 .or.           &
              wsp > 270. .or. tk1(k) > 380. .or. tk1(k) < 160. .or. &
              sqv1(k) < zero .or. sqc1(k) < zero .or.               &
              (ps-pres1(kts)) < zero) then
             kproblem = k
             problem = .true.
             print*,"Incoming problem at: i=",i," j=",j," k=",k
             print*," QFX=",qfx," HFX=",hfx
             print*," wsp=",wsp," T=",tk1(k)
             print*," qv=",sqv1(k)," qc=",sqc1(k)
             print*," u*=",ust," wspd=",wspd
             print*," xland=",xland," ts=",ts
             print*," ps=",ps,"delp1=",ps-pres1(kts)
             print*," znt=",znt," dx=",dx," dz(1)=",dz1(1)
          endif
       enddo
       if (problem) then
          print*,"===tk:",tk1(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qv:",sqv1(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qc:",sqc1(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qi:",sqi1(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====u:",u1(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====v:",v1(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====p:",pres1(max(kproblem-3,1):min(kproblem+3,kte))
       endif
    endif

!***  Begin debugging
    IMD=10
    JMD=10
!***  End debugging

    !initialize 3d output arrays
    edmf_a1      =zero
    edmf_w1      =zero
    edmf_qt1     =zero
    edmf_thl1    =zero
    edmf_ent1    =zero
    edmf_qc1     =zero
    edmf_u1      =zero
    edmf_v1      =zero
    edmf_qv1     =zero
    sub_thl1     =zero
    sub_sqv1     =zero
    det_thl1     =zero
    det_sqv1     =zero
    edmf_a_dd1   =zero
    edmf_w_dd1   =zero
    edmf_qc_dd1  =zero
    !initialize 2d fields
    ktop_plume   =0    !int
    ztop_plume   =zero
    maxwidth     =zero
    maxmf        =zero
    excess_h     =zero
    excess_q     =zero
    kzero1       =zero

    ! DH* CHECK HOW MUCH OF THIS INIT IF-BLOCK IS ACTUALLY NEEDED FOR RESTARTS
!> - Within the MYNN-EDMF, there is a dependecy check for the first time step,
!! If true, a three-dimensional initialization loop is entered. Within this loop,
!! several arrays are initialized and k-oriented (vertical) subroutines are called 
!! at every i and j point, corresponding to the x- and y- directions, respectively.  
    if (initflag > 0 .and. .not.restart) then

       !test to see if we want to initialize qke
       if ( (restart .or. cycling)) then
          if (maxval(qke1(:)) < 0.0002) then
             initialize_qke = .true.
             !print*,"qke is too small, must initialize"
          else
             initialize_qke = .false.
             !print*,"using background qke, will not initialize"
          endif
       else ! not cycling or restarting:
          initialize_qke = .true.
          !print*,"not restart nor cycling, must initialize qke"
       endif
 
       if (.not.restart .or. .not.cycling) THEN
          sh1         =zero
          sm1         =zero
          el1         =zero
          tsq1        =zero
          qsq1        =zero
          cov1        =zero
          cldfra_bl1  =zero
          qc_bl1      =zero
          qi_bl1      =zero
          qke1        =zero
          qke_adv1    =zero
       end if
       dqc1           =zero
       dqi1           =zero
       dqni1          =zero
       dqnc1          =zero
       dqnwfa1        =zero
       dqnifa1        =zero
       dqnbca1        =zero
       dozone1        =zero
       qc_bl1_old     =zero
       cldfra_bl1_old =zero
       sgm1           =zero
       vt1            =zero
       vq1            =zero
       km1            =zero
       kh1            =zero

       if (tke_budget .eq. 1) then
          qwt1        =zero
          qshear1     =zero
          qbuoy1      =zero
          qdiss1      =zero
          dqke1       =zero
       endif

       zw1(kts)=zero
       do k=kts,kte
          !keep snow out for now - increases ceiling bias
          sqw1(k)=sqv1(k)+sqc1(k)+sqi1(k)!+sqs1(k)
          thl1(k)=th1(k) - xlvcp/ex1(k)*sqc1(k) &
               &         - xlscp/ex1(k)*(sqi1(k))!+sqs1(k))
          !Use form from Tripoli and Cotton (1981) with their
          !suggested min temperature to improve accuracy.
          !thl1(k)=th1(k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc1(k) &
          !    &             - xlscp/MAX(tk1(k),TKmin)*sqi1(k))
          thlv1(k)=thl1(k)*(one+p608*sqv1(k))
          thv1(k)=th1(k)*(one+p608*sqv1(k) - (sqc1(k)+sqi1(k)))
          zw1(k+1)=zw1(k)+dz1(k)
       enddo

       if (INITIALIZE_QKE) then
          !Initialize tke for initial PBLH calc only - using
          !simple PBLH form of Koracin and Berkowicz (1988, BLM)
          !to linearly taper off tke towards top of PBL.
          do k=kts,kte
             qke1(k)=five * ust * MAX((ust*700._kind_phys - zw1(k))/(MAX(ust,p01)*700.), p01)
          enddo
       endif


!>  - Call get_pblh() to calculate hybrid (\f$\theta_{v}-TKE\f$) PBL height.
       call get_pblh(kts,kte,pblh,thv1,qke1,ust,zw1,dz1,xland,kpbl)
             
!>  - Call scale_aware() to calculate similarity functions for scale-adaptive control
!! (\f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$).
       call scale_aware(dx,pblh,psig_bl,psig_shcu)

       ! DH* CHECK IF WE CAN DO WITHOUT CALLING THIS ROUTINE FOR RESTARTS
!>  - Call mym_initialize() to initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$. These variables are calculated after 
!! obtaining prerequisite variables by calling the following subroutines from 
!! within mym_initialize(): mym_level2() and mym_length().
       rmol = zero
       rmolh= zero
       call mym_initialize (                & 
            &kts,kte,ckmod,taue,xland,pblh, &
            &dz1, dx, zw1, pres1, ex1,      &
            &u1, v1, thl1, sqv1,            &
            &th1, thv1, thlv1,              &
            &sh1, sm1,                      &
            &ust, rmol, rmolh,              &
            &el1, qke1,                     &
            &tsq1, qsq1, cov1,              &
            &psig_bl, cldfra_bl1,           &
            &bl_mynn_mixlength,             &
            &bl_mynn_ess,                   &
            &edmf_w1,edmf_a1,               &
            &edmf_w_dd1,edmf_a_dd1,         &
            &initialize_qke,                &
            &spp_pbl,pattern_spp_pbl1       )

       if (.not.restart) then
          !initialize qke_adv array if using advection
          IF (bl_mynn_tkeadvect) THEN
             DO k=KTS,KTE
                qke_adv1(k)=qke1(k)
             ENDDO
          ENDIF
       ENDIF

!***  Begin debugging
!       IF (I==IMD .AND. J==JMD) THEN
!          PRINT*,"MYNN DRIVER INIT: k=",1," sh=",sh1(k)
!          PRINT*," sqw=",sqw1(k)," thl=",thl1(k)," k_m=",km1(k)
!          PRINT*," xland=",xland," rmol=",rmol," ust=",ust
!          PRINT*," qke=",qke1(k)," el=",el1(k)," tsq=",tsq1(k)
!          PRINT*," PBLH=",PBLH," u=",u1(k)," v=",v1(k)
!       ENDIF
!***  End debugging

    ENDIF ! end initflag

!> - After initializing all required variables, the regular procedures 
!! performed at every time step are ready for execution.
    IF (bl_mynn_tkeadvect) THEN
       qke1(kts:kte)=qke_adv1(kts:kte)
    ENDIF

    !Initialize some arrays
    if (tke_budget .eq. 1) then
       dqke1(kts:kte)=qke1(kts:kte)
    endif
    cldfra_bl1_old(kts:kte)=cldfra_bl1(kts:kte)
    qc_bl1_old(kts:kte)=qc_bl1(kts:kte)
    qi_bl1_old(kts:kte)=qi_bl1(kts:kte)
    cldfra_bl1 =zero
    qc_bl1     =zero
    qi_bl1     =zero
    dqc1       =zero
    dqi1       =zero
    dqs1       =zero
    dqni1      =zero
    dqnc1      =zero
    dqnwfa1    =zero
    dqnifa1    =zero
    dqnbca1    =zero
    dozone1    =zero
    !edmf
    edmf_a1    =zero
    edmf_w1    =zero
    edmf_qc1   =zero
    edmf_qv1   =zero
    edmf_thl1  =zero
    edmf_u1    =zero
    edmf_v1    =zero
    s_aw1      =zero
    s_awthl1   =zero
    s_awqt1    =zero
    s_awqv1    =zero
    s_awqc1    =zero
    s_awu1     =zero
    s_awv1     =zero
    s_awqke1   =zero
    s_awqsq1   =zero
    s_awqnc1   =zero
    s_awqni1   =zero
    s_awqnwfa1 =zero
    s_awqnifa1 =zero
    s_awqnbca1 =zero
    s_awchem1(kts:kte+1,1:nchem) = zero
    s_awscalars1(kts:kte+1,1:nscalars) = zero
    !edmf-downdraft
    edmf_a_dd1 =zero
    edmf_w_dd1 =zero
    edmf_qc_dd1=zero
    sd_aw1     =zero
    sd_awthl1  =zero
    sd_awqt1   =zero
    sd_awqv1   =zero
    sd_awqc1   =zero
    sd_awqi1   =zero
    sd_awqnc1  =zero
    sd_awqni1  =zero
    sd_awqnwfa1=zero
    sd_awqnifa1=zero
    sd_awu1    =zero
    sd_awv1    =zero
    sd_awqke1  =zero
    sub_thl1   =zero
    sub_sqv1   =zero
    sub_u1     =zero
    sub_v1     =zero
    det_thl1   =zero
    det_sqv1   =zero
    det_sqc1   =zero
    det_u1     =zero
    det_v1     =zero

    zw1(kts)=zero
    do k = kts,kte
       zw1(k+1)=zw1(k)+dz1(k)
       qv1(k) = sqv1(k)/(one-sqv1(k))
       qc1(k) =	sqc1(k)/(one-sqv1(k))
       qi1(k) =	sqi1(k)/(one-sqv1(k))
       qs1(k) =	sqs1(k)/(one-sqv1(k))
       !keep snow out for now - increases ceiling bias
       sqw1(k)= sqv1(k)+sqc1(k)+sqi1(k)!+sqs1(k)
       thl1(k)= th1(k) - xlvcp/ex1(k)*sqc1(k) &
             &         - xlscp/ex1(k)*(sqi1(k))!+sqs1(k))
       !Use form from Tripoli and Cotton (1981) with their
       !suggested min temperature to improve accuracy.
       !thl1(k)=th1(i,k)*(one- xlvcp/MAX(tk1(k),TKmin)*sqc1(k) &
       !    &               - xlscp/MAX(tk1(k),TKmin)*sqi1(k))
       thv1(k)=th1(k)*(one+p608*sqv1(k) - (sqc1(k)+sqi1(k)))
    enddo ! end k

!>  - Call get_pblh() to calculate the hybrid \f$\theta_{v}-TKE\f$
!! PBL height diagnostic.
    call get_pblh(kts,kte,pblh,thv1,qke1,ust,zw1,dz1,xland,kpbl)

!>  - Call scale_aware() to calculate the similarity functions,
!! \f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$, to control 
!! the scale-adaptive behaviour for the local and nonlocal 
!! components, respectively.
    call scale_aware(dx,pblh,psig_bl,psig_shcu)

    sqcg= zero   !ill-defined variable; qcg has been removed
    cpm=cp*(one + 0.84_kind_phys*max(sqv1(kts),1e-8_kind_phys))
    exnerg=(ps/p1000mb)**rcp

    !-----------------------------------------------------
    !ORIGINAL CODE
    !flt = hfx/( rho(i,kts)*cpm ) &
    ! +xlvcp*ch(i)*(sqc(kts)/exner(i,kts) -sqcg/exnerg)
    !flq = qfx/  rho(i,kts)       &
    !    -ch(i)*(sqc(kts)   -sqcg )
    !-----------------------------------------------------
    flqv = qfx/rho1(kts)
    flqc = zero !currently no sea-spray fluxes, fog settling handled elsewhere
    th_sfc = ts/ex1(kts)

    ! TURBULENT FLUX FOR TKE BOUNDARY CONDITIONS
    flq =flqv+flqc                   !! LATENT
    flt =hfx/(rho1(kts)*cpm )-xlvcp*flqc/ex1(kts)  !! Temperature flux
    fltv=flt + flqv*p608*th_sfc      !! Virtual temperature flux

    ! Update 1/L using updated sfc heat flux and friction velocity
    rmol= -karman*gtr*fltv/max(ust**3,1.0e-6_kind_phys)
    rmolh=-karman*gtr*flt /max(ust**3,1.0e-6_kind_phys)
    zet = p5*dz1(kts)*rmol
    zet = max(zet, -10._kind_phys)
    zet = min(zet,  10._kind_phys)
    !if(i.eq.idbg)print*,"updated z/L=",zet
    if (bl_mynn_stfunc == 0) then
       !Original Kansas-type stability functions
       if ((xland-1.5) .ge. zero) then       ! WATER
          if ( zet >= zero ) then
             !pmz = one + (cphm_st-one) * zet
             pmz = one + 3.0_kind_phys * zet
             phh = one +  cphh_st * zet
          else
             pmz = one/    (one-cphm_unst*zet)**0.25 - zet
             phh = one/sqrt(one-cphh_unst*zet)
          end if
       else                                  ! LAND
          if ( zet >= zero ) then
             !pmz = one + (cphm_st-one) * zet
             pmz = 0.95_kind_phys + 3.5_kind_phys * zet
             phh = one + cphh_st       * zet
          else
             !pmz = one/(one-cphm_unst*zet)**0.25 - zet
             pmz = 0.95_kind_phys/(one-cphm_unst*zet)**0.25 - zet
             phh = one       /sqrt(one-cphh_unst*zet)
          end if
       endif
       phi_m = phim(zet)
    else
       !Updated stability functions (Puhales, 2020)
       phi_m = phim(zet)
       pmz   = phi_m - zet
       phh   = phih(zet)
    end if

!>  - Call mym_condensation() to calculate the nonconvective component
!! of the subgrid cloud fraction and mixing ratio as well as the functions
!! used to calculate the buoyancy flux. Different cloud PDFs can be
!! selected by use of the namelist parameter \p bl_mynn_cloudpdf.

    call mym_condensation (kts,kte,                   &
         &dx,dz1,zw1,xland,                           &
         &thl1,sqw1,sqv1,sqc1,sqi1,sqs1,              &
         &pres1,ex1,tsq1,qsq1,cov1,                   &
         &sh1,el1,bl_mynn_cloudpdf,                   &
         &qc_bl1,qi_bl1,cldfra_bl1,                   &
         &pblh,hfx,                                   &
         &vt1, vq1, th1, sgm1,                        &
         &bl_mynn_closure,                            &
         &spp_pbl, pattern_spp_pbl1                   )


    TKEprod_up = zero
    if (bl_mynn_edmf > 0) then
       !PRINT*,"Calling DMP Mass-Flux"
       call DMP_mf(i,j,                               &
            &kts,kte,delt,zw1,dz1,pres1,rho1,         &
            &bl_mynn_edmf_mom,                        &
            &bl_mynn_edmf_tke,                        &
            &bl_mynn_mixscalars,                      &
            &bl_mynn_mixaerosols,                     &
            &bl_mynn_mixnumcon,                       &
            &bl_mynn_closure,                         &
            &u1,v1,w1,th1,thl1,thv1,tk1,              &
            &sqw1,sqv1,sqc1,qke1,qsq1,                &
            &qnc1,qni1,qnwfa1,qnifa1,qnbca1,          &
            &ex1,vt1,vq1,sgm1,                        &
            &ust,flt,fltv,flq,flqv,                   &
            &pblh,kpbl,dx,                            &
            &xland,th_sfc,                            &
            ! now outputs - tendencies
            ! &,dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf      &
            ! outputs - updraft properties
            &edmf_a1,edmf_w1,edmf_qt1,                &
            &edmf_thl1,edmf_ent1,edmf_qc1,            &
            &edmf_qv1,edmf_u1,edmf_v1,                &
            ! for the solver
            &s_aw1,s_awthl1,s_awqt1,                  &
            &s_awqv1,s_awqc1,                         &
            &s_awu1,s_awv1,s_awqke1,s_awqsq1,         &
            &s_awqnc1,s_awqni1,                       &
            &s_awqnwfa1,s_awqnifa1,s_awqnbca1,        &
            &sub_thl1,sub_sqv1,                       &
            &sub_u1,sub_v1,                           &
            &det_thl1,det_sqv1,det_sqc1,              &
            &det_u1,det_v1,                           &
            ! chem/smoke mixing
            &nchem,chem1,s_awchem1,                   &
            &mix_chem,                                &
            &nscalars,scalars,s_awscalars1,           &
            &qc_bl1,cldfra_bl1,                       &
            &qc_bl1_old,cldfra_bl1_old,               &
            &FLAG_QC,FLAG_QI,                         &
            &FLAG_QNC,FLAG_QNI,                       &
            &FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,        &
            &Psig_shcu,                               &
            &maxwidth,ktop_plume,                     &
            &maxmf,ztop_plume,excess_h,excess_q,      &
            &spp_pbl,pattern_spp_pbl1,                &
            &TKEprod_up,el1                           )
    endif

    !> Options for downdrafts or analytical top-down method:
    !!  - Add TKE source driven by cloud top cooling and
    !!  calculate the buoyancy production of TKE from cloud-top cooling
    !!  for downdraft or analytic options. 
    tkeprod_dn     = zero
    maxtkeprod     = zero
    cldtop_cooling = zero
    ent_eff        = zero
    maxmf_dd       = zero
    maxwidth_dd    = zero
    if (bl_mynn_edmf_dd > 0) then
       call ddmp_mf(kts,kte,delt,dx,zw1,dz1,pres1,    &
            &u1,v1,th1,thl1,thv1,tk1,                 &
            &sqw1,sqv1,sqc1,sqi1,qnc1,qni1,           &
            &qnwfa1,qnifa1,                           &
            &qke1,rho1,ex1,                           &
            &qc_bl1,qi_bl1,cldfra_bl1,                &
            &ust,flt,flq,fltv,                        &
            &pblh,kpbl,                               &
            &edmf_a_dd1,edmf_w_dd1,edmf_qt_dd1,       &
            &edmf_thl_dd1,edmf_ent_dd1,               &
            &edmf_qc_dd1,                             &
            &sd_aw1,sd_awthl1,sd_awqt1,               &
            &sd_awqv1,sd_awqc1,sd_awqi1,              &
            &sd_awqnc1,sd_awqni1,                     &
            &sd_awqnwfa1,sd_awqnifa1,                 &
            &sd_awu1,sd_awv1,                         &
            &sd_awqke1,                               &
            &tkeprod_dn,el1,                          &
            &rthraten1,psig_bl,                       &
            &maxmf_dd,maxwidth_dd                     )
       
       !make sure there is some tke production for shallow fog,
       !when the nonlocal (mf) approach is no longer appropriate.
       call topdown_cloudrad(kts,kte,                 &
            &dz1,zw1,fltv,u1(kts),v1(kts),            &
            &xland,kpbl,pblh,                         &
            &sqc1,sqi1,sqw1,thl1,th1,                 &
            &ex1,pres1,rho1,thv1,                     &
            &cldfra_bl1,qc_bl1,qi_bl1,rthraten1,      &
            &tkeprod_dn,psig_bl,                      &
            maxtkeprod,cldtop_cooling,ent_eff         )
    endif

    !Capability to substep the eddy-diffusivity portion
    !do nsub = 1,2
    delt2 = delt !*0.5    !only works if topdown=0

    !thlv with updated subgrid clouds
    do k=kts,kte
       qc_tot1(k) =max(qc_bl1(k),sqc1(k))
       qi_tot1(k) =max(qi_bl1(k),sqi1(k)+sqs1(k))
       thlv1(k)   =(th1(k) - xlvcp/ex1(k)*qc_tot1(k)  &
            &              - xlscp/ex1(k)*qi_tot1(k)) &
            &              * (one+p608*sqv1(k))
       thl_tot1(k)=th1(k)  - xlvcp/ex1(k)*qc_tot1(k)  &
            &              - xlscp/ex1(k)*qi_tot1(k)
       thv1(k)    =th1(k)*(one+p608*sqv1(k) - (qc_tot1(k)+qi_tot1(k)))
    enddo

    call mym_turbulence(                                 &
            &kts,kte,xland,bl_mynn_closure,ckmod,taue,   &
            &dz1, dx, zw1, pres1, ex1,                   &
            &u1, v1, thl1, thv1, thlv1,                  &
            &sqc1, sqw1,                                 &
            &qke1, tsq1, qsq1, cov1,                     &
            &vt1, vq1,                                   &
            &rmol, rmolh, flt, fltv, flq,                &
            &pblh, th1,                                  &
            &sh1,sm1,el1,                                &
            &Dfm1,Dfh1,Dfq1,                             &
            &Tcd1,Qcd1,Pdk1,                             &
            &Pdt1,Pdq1,Pdc1,                             &
            &qWT1,qSHEAR1,qBUOY1,qDISS1,                 &
            &tke_budget,                                 &
            &Psig_bl,Psig_shcu,                          &
            &cldfra_bl1,bl_mynn_mixlength,               &
            &bl_mynn_ess,                                &
            &edmf_w1,edmf_a1,                            &
            &edmf_w_dd1,edmf_a_dd1,                      &
            &TKEprod_dn,TKEprod_up,                      &
            &spp_pbl,pattern_spp_pbl1                    )

!>  - Call mym_predict() to solve TKE and 
!! \f$\theta^{'2}, q^{'2}, and \theta^{'}q^{'}\f$
!! for the following time step.
    call mym_predict(kts,kte,bl_mynn_closure,            &
            &delt2, dz1,                                 &
            &ust, flt, flq, pmz, phh,                    &
            &el1, dfq1, rho1, pdk1, pdt1, pdq1, pdc1,    &
            &qke1, tsq1, qsq1, cov1,                     &
            &s_aw1, s_awqke1, s_awqsq1,                  &
            &bl_mynn_edmf, bl_mynn_edmf_tke,             &
            &qWT1, qDISS1, tke_budget, upwind            )

!>  - Calculate the heating due to dissipation of TKE
!! Set max dissipative heating rate to 7.2 K per hour and limit heating above 100 mb.
    do k=kts,kte-1
       diss_heat1(k) = min(max((qke1(k)**1.5_kind_phys)/(b1*max(p5*(el1(k)+el1(k+1)),one))/cp, zero),0.002_kind_phys)
       diss_heat1(k) = diss_heat1(k) * exp(-10000._kind_phys/max(pres1(k),one)) 
    enddo
    diss_heat1(kte) = zero

!>  - Call mynn_tendencies() to solve for tendencies of 
!! \f$U, V, \theta, q_{v}, q_{c}, and q_{i}\f$.
    call mynn_tendencies(kts,kte,i,                      &
            &delt, dz1, zw1, xland, pblh, rho1,          &
            &u1, v1, th1, tk1, qv1,                      &
            &qc1, qi1, kzero1, qnc1, qni1,               & !kzero replaces qs1 - not mixing snow
            &ps, pres1, ex1, thl1,                       &
            &sqv1, sqc1, sqi1, kzero1, sqw1,             & !kzero replaces sqs - not mixing snow
            &thl_tot1, qc_tot1, qi_tot1,                 &
            &qnwfa1, qnifa1, qnbca1, ozone1,             &
            &ust,flt,flq,flqv,flqc,                      &
            &wspd,uoce,voce,                             &
            &tsq1, qsq1, cov1,                           &
            &tcd1, qcd1,                                 &
            &dfm1, dfh1, dfq1,                           &
            &Du1, Dv1, Dth1, Dqv1,                       &
            &Dqc1, Dqi1, Dqs1, Dqnc1, Dqni1,             &
            &Dqnwfa1, Dqnifa1, Dqnbca1,                  &
            &Dozone1,                                    &
            &diss_heat1,                                 &
            ! mass flux components
            &s_aw1,s_awthl1,s_awqt1,                     &
            &s_awqv1,s_awqc1,s_awu1,s_awv1,              &
            &s_awqnc1,s_awqni1,                          &
            &s_awqnwfa1,s_awqnifa1,s_awqnbca1,           &
            &sd_aw1,sd_awthl1,sd_awqt1,                  &
            &sd_awqv1,sd_awqc1,sd_awqi1,                 &
            &sd_awqnc1,sd_awqni1,                        &
            &sd_awqnwfa1,sd_awqnifa1,                    &
            &sd_awu1,sd_awv1,                            &
            &sub_thl1,sub_sqv1,                          &
            &sub_u1,sub_v1,                              &
            &det_thl1,det_sqv1,det_sqc1,                 &
            &det_u1,det_v1,                              &
            &FLAG_QC,FLAG_QI,FLAG_QNC,                   &
            &FLAG_QNI,FLAG_QS,                           &
            &FLAG_QNWFA,FLAG_QNIFA,                      &
            &FLAG_QNBCA,FLAG_OZONE,                      &
            &cldfra_bl1,                                 &
            &bl_mynn_cloudmix,                           &
            &bl_mynn_mixqt,                              &
            &bl_mynn_edmf,                               &
            &bl_mynn_edmf_mom,                           &
            &bl_mynn_mixscalars,                         &
            &bl_mynn_mixaerosols,                        &
            &bl_mynn_mixnumcon,                          &
            &upwind                                      )

    if ( mix_chem ) then
          call mynn_mix_chem(kts,kte,i,                  &
               &delt, dz1, pblh,                         &
               &nchem, ndvel,                            &
               &chem1, vdep, settle1,                    &
               &rho1, flt,                               &
               &tcd1, qcd1,                              &
               &dfh1,                                    &
               &s_aw1,s_awchem1,                         &
               &emis_ant_no,                             &
               &frp,                                     &
               &enh_mix,                                 &
               &bl_mynn_edmf, upwind                     )
       do ic = 1,nchem
          do k = kts,kte
             chem1(k,ic) = max(1.e-12, chem1(k,ic))
          enddo
       enddo
    endif

    if ( bl_mynn_mixscalars .eq. 1) then
       call mynn_mix_scalars(kts,kte,i,               &
            &delt, dz1,                               &
            &nscalars, scalars,                       &
            &rho1, flt,                               &
            &tcd1, qcd1,                              &
            &dfh1,                                    &
            &s_aw1,s_awscalars1,                      &
            &bl_mynn_edmf, upwind                     )
    endif
       
    
    call retrieve_exchange_coeffs(kts,kte,               &
         dfm1, dfh1, dz1, km1, kh1                       )

    !update the TKE budget
    if (tke_budget .eq. 1) then
       !! TKE budget is now given in m**2/s**-3 (Puhales, 2020)
       !! Lower boundary condtions (using similarity relationships such as the prognostic equation for Qke)
       k=kts
       qSHEAR1(k)   = 4.*(ust**3*phi_m/(karman*dz1(k)))-qSHEAR1(k+1) !! staggered
       qBUOY1(k)    = 4.*(-ust**3*zet/(karman*dz1(k)))-qBUOY1(k+1) !! staggered
       !! unstaggering SHEAR and BUOY and trasfering all TKE budget to 3D array               
       do k = kts,kte-1
          dummy1(k) = p5*(qSHEAR1(k)+qSHEAR1(k+1)) !!! unstaggering in z
          dummy2(k) = p5*(qBUOY1(k)+qBUOY1(k+1)) !!! unstaggering in z
          dqke1(k)  = p5*(qke1(k)-dqke1(k))/delt
       enddo
       qSHEAR1      = dummy1
       qBUOY1       = dummy2
       !! Upper boundary conditions
       k=kte
       qSHEAR1(k)   = zero
       qBUOY1(k)    = zero
       qWT1(k)      = zero
       qDISS1(k)    = zero
       dqke1(k)     = zero
    endif

    !update updraft/downdraft properties
    if (bl_mynn_output > 0) then !research mode == 1 or 2
       !if mode 2, then overwrite updrafts with downdrafts
       if (bl_mynn_output == 2 .and. bl_mynn_edmf_dd == 1) then
          edmf_a1(kts:kte)   =edmf_a_dd1(kts:kte)
          edmf_w1(kts:kte)   =edmf_w_dd1(kts:kte)
          edmf_qt1(kts:kte)  =edmf_qt_dd1(kts:kte)
          edmf_thl1(kts:kte) =edmf_thl_dd1(kts:kte)
          edmf_ent1(kts:kte) =edmf_ent_dd1(kts:kte)
          edmf_qc1(kts:kte)  =edmf_qc_dd1(kts:kte)
       endif
    endif

    !***  Begin debug prints
    if ( debug_code .and. (i .eq. idbg)) THEN
       if ( ABS(QFX)>.001)print*,&
          "SUSPICIOUS VALUES AT: i=",i," QFX=",QFX
       if ( ABS(HFX)>1100.)print*,&
          "SUSPICIOUS VALUES AT: i=",i," HFX=",HFX
       do k = kts,kte
          IF ( sh1(k) < 0. .OR. sh1(k)> 200.)print*,&
             "SUSPICIOUS VALUES AT: i,k=",i,k," sh=",sh1(k)
          IF ( ABS(vt1(k)) > 2.0 )print*,&
             "SUSPICIOUS VALUES AT: i,k=",i,k," vt=",vt1(k)
          IF ( ABS(vq1(k)) > 7000.)print*,&
             "SUSPICIOUS VALUES AT: i,k=",i,k," vq=",vq1(k)
          IF ( qke1(k) < -1. .OR. qke1(k)> 200.)print*,&
             "SUSPICIOUS VALUES AT: i,k=",i,k," qke=",qke1(k)
          IF ( el1(k) < 0. .OR. el1(k)> 1500.)print*,&
             "SUSPICIOUS VALUES AT: i,k=",i,k," el1=",el1(k)
          IF ( km1(k) < 0. .OR. km1(k)> 2000.)print*,&
             "SUSPICIOUS VALUES AT: i,k=",i,k," km=",km1(k)
          IF (cldfra_bl1(k) < zero .OR. cldfra_bl1(k)> one)THEN
             PRINT*,"SUSPICIOUS VALUES: CLDFRA_BL=",cldfra_bl1(k),&
                                          " qc_bl=",qc_bl1(k)
          endif
       enddo !end-k
    endif

    !ACF copy qke into qke_adv if using advection
    qke_adv1(kts:kte)=qke1(kts:kte)

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mynnedmf
!> @}

!=======================================================================
!     SUBROUTINE  mym_initialize:
!
!     Input variables:
!       iniflag         : <>0; turbulent quantities will be initialized
!                         = 0; turbulent quantities have been already
!                              given, i.e., they will not be initialized
!       nx, nz          : Dimension sizes of the
!                         x and z directions, respectively
!       tref            : Reference temperature                      (K)
!       dz(nz)          : Vertical grid spacings                     (m)
!                         # dz(nz)=dz(nz-1)
!       zw(nz+1)        : Heights of the walls of the grid boxes     (m)
!                         # zw(1)=0.0 and zw(k)=zw(k-1)+dz(k-1)
!       exner(nx,nz)    : Exner function at zw*h+zg             (J/kg K)
!                         defined by c_p*( p_basic/1000hPa )^kappa
!                         This is usually computed by integrating
!                         d(pi0)/dz = -h*g/tref.
!       rmol(nx)        : Inverse of the Obukhov length         (m^(-1))
!       flt, flq(nx)    : Turbulent fluxes of potential temperature and
!                         total water, respectively:
!                                    flt=-u_*Theta_*             (K m/s)
!                                    flq=-u_*qw_*            (kg/kg m/s)
!       ust(nx)         : Friction velocity                        (m/s)
!       pmz(nx)         : phi_m-zeta at z1*h+z0, where z1 (=0.5*dz(1))
!                         is the first grid point above the surafce, z0
!                         the roughness length and zeta=(z1*h+z0)*rmo
!       phh(nx)         : phi_h at z1*h+z0
!       u, v(nx,nz)     : Components of the horizontal wind        (m/s)
!       thl(nx,nz)      : Liquid water potential temperature
!                                                                    (K)
!       qw(nx,nz)       : Total water content Q_w                (kg/kg)
!
!     Output variables:
!       ql(nx,nz)       : Liquid water content                   (kg/kg)
!       vt, vq(nx,nz)   : Functions for computing the buoyancy flux
!       qke(nx,nz)      : Twice the turbulent kinetic energy q^2
!                                                              (m^2/s^2)
!       tsq(nx,nz)      : Variance of Theta_l                      (K^2)
!       qsq(nx,nz)      : Variance of Q_w
!       cov(nx,nz)      : Covariance of Theta_l and Q_w              (K)
!       el(nx,nz)       : Master length scale L                      (m)
!                         defined on the walls of the grid boxes
!
!     Work arrays:        see subroutine mym_level2
!       pd?(nx,nz,ny) : Half of the production terms at Level 2
!                         defined on the walls of the grid boxes
!       qkw(nx,nz,ny) : q on the walls of the grid boxes         (m/s)
!
!     # As to dtl, ...gh, see subroutine mym_turbulence.
!
!-------------------------------------------------------------------

!>\ingroup gsd_mynn_edmf
!! This subroutine initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$.
!!\section gen_mym_ini GSD MYNN-EDMF mym_initialize General Algorithm 
!> @{
  SUBROUTINE  mym_initialize (                                & 
       &            kts,kte,ckmod,taue,xland,pblh,            &
       &            dz, dx, zw, p, exner,                     &
       &            u, v, thl, qw,                            &
!       &            ust, rmo, pmz, phh, flt, flq,             &
       &            theta, thv, thlv, sh, sm,                 &
       &            ust, rmol, rmolh, el,                     &
       &            Qke, Tsq, Qsq, Cov,                       &
       &            Psig_bl, cldfra_bl1,                      &
       &            bl_mynn_mixlength,                        &
       &            bl_mynn_ess,                              &
       &            edmf_w1,edmf_a1,                          &
       &            edmf_w_dd1,edmf_a_dd1,                    &
       &            INITIALIZE_QKE,                           &
       &            spp_pbl,pattern_spp_pbl1                  )
!
!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: b1,b2,gtr,karman,    &
         qkemin,zero,one,two,three,hundred,p01,p333,p5,p666,  &
         kind_phys
    
    integer, intent(in)           :: kts,kte
    integer, intent(in)           :: bl_mynn_mixlength
    integer, intent(in)           :: bl_mynn_ess
    logical, intent(in)           :: INITIALIZE_QKE
!    real(kind_phys), intent(in)   :: ust, rmol, pmz, phh, flt, flq
    real(kind_phys), intent(in)   :: rmol, rmolh, Psig_bl, xland, pblh
    real(kind_phys), intent(in)   :: dx, ust, ckmod, taue
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), dimension(kts:kte),   intent(in) :: u,v,thl,&
         &thlv,qw,cldfra_bl1,edmf_w1,edmf_a1,edmf_w_dd1,edmf_a_dd1,&
         &p,exner
    real(kind_phys), dimension(kts:kte),   intent(inout) :: tsq,qsq,cov
    real(kind_phys), dimension(kts:kte),   intent(inout) :: el,qke
    real(kind_phys), dimension(kts:kte) ::                       &
         &ql,pdk,pdt,pdq,pdc,dtl,dqw,dtv,                        &
         &gm,gh,sm,sh,qkw,qtw,vt,vq,qpe
    integer :: k,l,lmax
    real(kind_phys):: phm,vkz,elq,elv,b1l,b2l,pmz=1.,phh=1.,     &
         &flt=0.,fltv=0.,flq=0.,tmpq
    real(kind_phys):: ustm
    real(kind_phys), dimension(kts:kte) :: theta,thv
    real(kind_phys), dimension(kts:kte) :: pattern_spp_pbl1
    integer ::spp_pbl

!> - At first ql, vt and vq are set to zero.
    DO k = kts,kte
       ql(k) = zero
       vt(k) = zero
       vq(k) = zero
       qkw(k)= zero
       qtw(k)= zero
    END DO
    ustm = max(p01, ust)
!
!> - Call mym_level2() to calculate the stability functions at level 2.
    CALL mym_level2 ( kts,kte,                      &
         &            ckmod, bl_mynn_ess,           &
         &            zw, dz, xland, pblh,          &
         &            u, v, thl, thv, thlv,         &
         &            theta, p, exner,              &
         &            qke, cldfra_bl1,              &
         &            edmf_a1, edmf_w1,             &
         &            qw, ql, vt, vq,               &
         &            dtl, dqw, dtv, gm, gh, sm, sh )
!
!   **  Preliminary setting  **

    el(kts) = zero
    IF (INITIALIZE_QKE) THEN
       !qke(kts) = ust**2 * ( b1*pmz )**(2.0/3.0)
       qke(kts) = 1.5_kind_phys * ustm**2 * ( b1*pmz )**p666
       DO k = kts+1,kte
          !qke(k) = 0.0
          !linearly taper off towards top of pbl
          qke(k)=qke(kts)*MAX((ust*700._kind_phys - zw(k))/(ustm*700._kind_phys), p01)
       ENDDO
    ENDIF
!
    phm      = phh*b2 / ( b1*pmz )**p333
    tsq(kts) = min(three,          phm*( max(zero,flt)/ustm )**2)
    qsq(kts) = min(5e-6_kind_phys, phm*( max(zero,flq)/ustm )**2)
    cov(kts) = phm*( flt/ustm )*( flq/ustm )
    qpe(kts) = min(three, max(zero, two * tsq(kts) * gtr**2 * taue**2))
!
    DO k = kts+1,kte
       vkz = karman*zw(k)
       el(k) = vkz/( one + vkz/hundred )
!       qke(k) = 0.0
!
       tsq(k) = zero
       qsq(k) = zero
       cov(k) = zero
       qpe(k) = zero
    END DO
!
!   **  Initialization with an iterative manner          **
!   **  lmax is the iteration count. This is arbitrary.  **
    lmax = 5
!
    DO l = 1,lmax
!
!> - call mym_length() to calculate the master length scale.
       CALL mym_length (                          &
            &            kts,kte,xland,           &
            &            dz, dx, zw,              &
            &            rmol, rmolh,             &
            &            flt, fltv, flq,          &
            &            vt, vq,                  &
            &            u, v, qke, qpe,          &
            &            dtv,                     &
            &            el,                      &
            &            pblh,theta,              &
            &            qkw, qtw,                &
            &            Psig_bl,cldfra_bl1,      &
            &            bl_mynn_mixlength,       &
            &            edmf_w1,edmf_a1,         &
            &            edmf_w_dd1,edmf_a_dd1    )
!
       DO k = kts+1,kte
          elq    = el(k)*qkw(k)
          pdk(k) = elq*( sm(k)*gm(k) + &
               &         sh(k)*gh(k) )
          pdt(k) = elq*  sh(k)*dtl(k)**2
          pdq(k) = elq*  sh(k)*dqw(k)**2
          pdc(k) = elq*  sh(k)*dtl(k)*dqw(k)
       END DO
!
!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
       vkz = karman*p5*dz(kts)
       elv = p5*( el(kts+1)+el(kts) ) /  vkz
       IF (INITIALIZE_QKE)THEN 
          !qke(kts) = ust**2 * ( b1*pmz*elv    )**(2.0/3.0)
          qke(kts) = one * ustm**2 * ( b1*pmz*elv    )**p666 
       ENDIF

       phm      = phh*b2 / ( b1*pmz/elv**2 )**p333
       tsq(kts) = min(three,          phm*( max(zero,flt)/ustm )**2)
       qsq(kts) = min(5e-6_kind_phys, phm*( max(zero,flq)/ustm )**2)
       cov(kts) = phm*( flt/ustm )*( flq/ustm )
       
       DO k = kts+1,kte-1
          b1l = b1*0.25*( el(k+1)+el(k) )
          !tmpq=MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin)
          !add MIN to limit unreasonable QKE
          tmpq=MIN(MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin),125._kind_phys)
          !print*,'tmpq=',tmpq,pdk(k+1),pdk(k)
          IF (INITIALIZE_QKE)THEN
             qke(k) = tmpq**p666
          ENDIF

          IF ( qke(k) .LE. zero ) THEN
             b2l = zero
          ELSE
             b2l = b2*( b1l/b1 ) / SQRT( qke(k) )
          END IF

          tsq(k) = min(three,          max(zero, b2l*( pdt(k+1)+pdt(k) )))
          qsq(k) = min(5e-6_kind_phys, max(zero, b2l*( pdq(k+1)+pdq(k) )))
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO

    END DO

    IF (INITIALIZE_QKE)THEN
       qke(kts)=p5*(qke(kts)+qke(kts+1))
       qke(kte)=qke(kte-1)
    ENDIF
    tsq(kte)=tsq(kte-1)
    qsq(kte)=qsq(kte-1)
    cov(kte)=cov(kte-1)

  END SUBROUTINE mym_initialize
!> @}
  
!
! ==================================================================
!     SUBROUTINE  mym_level2:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:
!       dtl(nx,nz,ny) : Vertical gradient of Theta_l             (K/m)
!       dqw(nx,nz,ny) : Vertical gradient of Q_w
!       dtv(nx,nz,ny) : Vertical gradient of Theta_V             (K/m)
!       gm (nx,nz,ny) : G_M divided by L^2/q^2                (s^(-2))
!       gh (nx,nz,ny) : G_H divided by L^2/q^2                (s^(-2))
!       sm (nx,nz,ny) : Stability function for momentum, at Level 2
!       sh (nx,nz,ny) : Stability function for heat, at Level 2
!
!       These are defined on the walls of the grid boxes.
!

!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the level 2, non-dimensional wind shear
!! \f$G_M\f$ and vertical temperature gradient \f$G_H\f$ as well as 
!! the level 2 stability funcitons \f$S_h\f$ and \f$S_m\f$.
!!\param kts    horizontal dimension
!!\param kte    vertical dimension
!!\param dz     vertical grid spacings (\f$m\f$)
!!\param u      west-east component of the horizontal wind (\f$m s^{-1}\f$)
!!\param v      south-north component of the horizontal wind (\f$m s^{-1}\f$)
!!\param thl    liquid water potential temperature
!!\param qw     total water content \f$Q_w\f$
!!\param ql     liquid water content (\f$kg kg^{-1}\f$)
!!\param vt
!!\param vq
!!\param dtl     vertical gradient of \f$\theta_l\f$ (\f$K m^{-1}\f$)
!!\param dqw     vertical gradient of \f$Q_w\f$
!!\param dtv     vertical gradient of \f$\theta_V\f$ (\f$K m^{-1}\f$)
!!\param gm      \f$G_M\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
!!\param gh      \f$G_H\f$ divided by \f$L^{2}/q^{2}\f$ (\f$s^{-2}\f$)
!!\param sm      stability function for momentum, at Level 2
!!\param sh      stability function for heat, at Level 2
!!\section gen_mym_level2 GSD MYNN-EDMF mym_level2 General Algorithm
!! @ {
  SUBROUTINE  mym_level2 (kts,kte,                &
       &            ckmod, bl_mynn_ess,           &
       &            zw, dz, xland, pblh,          &
       &            u, v, thl, thv, thlv,         &
       &            th, p, exner,                 &
       &            qke, cldfra,                  &
       &            edmf_a, edmf_w,               &
       &            qw, ql, vt, vq,               &
       &            dtl, dqw, dtv, gm, gh, sm, sh )
!
!-------------------------------------------------------------------
 use module_bl_mynnedmf_common, only: tv0,gtr,cp, &
      qkemin,a1,a2,b1,c1,c2,c5,g1,g2,zero,one,two,&
      three,four,ten,twenty,hundred,p1,p2,p25,p3, &
      p5,p95,kind_phys
    
 integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

 real(kind_phys), intent(in)::xland,pblh,ckmod
 integer, intent(in)        ::bl_mynn_ess
 real(kind_phys), dimension(kts:kte), intent(in)  :: dz
 real(kind_phys), dimension(kts:kte), intent(in)  :: u,v, &
      &thl,qw,ql,vt,vq,thv,thlv,th,p,exner,edmf_a,edmf_w, &
      &qke,cldfra
 real(kind_phys), dimension(kts:kte+1), intent(in):: zw
 real(kind_phys), dimension(kts:kte), intent(out) ::      &
      &dtl,dqw,dtv,gm,gh,sm,sh
 real(kind_phys), dimension(kts:kte):: ri
 integer :: k
 real(kind_phys):: rfc,f1,f2,rf1,rf2,smc,shc,             &
      &ri1,ri2,ri3,ri4,duz,dtz,dqz,vtt,vqq,dtq,dzk,       &
      &afk,abk,rf
 real(kind_phys):: a2fac,thv1,thv0,eth0,eth1,lambda,qsat, &
      &xl,tabs,clam0,clam,wt,mfi,qkei,cldfrai,ugrid,      &
      &cfac,uonset,taper
 real(kind_phys), parameter:: thvp = 0.0 !0.25 !percentage of thv in blend

!    ev  = 2.5e6
!    tv0 = 0.61*tref
!    tv1 = 1.61*tref
!    gtr = 9.81/tref
!
!intialize output
 dtl(kts)=zero
 dqw(kts)=zero
 dtv(kts)=zero
 gm(kts)=zero
 gh(kts)=zero
 sm(kts)=zero
 sh(kts)=zero
 ri(kts)=zero

 select case (bl_mynn_ess)  !moist static stability

    case (0) !use buoyancy flux functions to calculate effective static stability (ess)

       do k = kts+1,kte
          dzk = p5 *( dz(k)+dz(k-1) )
          afk = dz(k)/( dz(k)+dz(k-1) )
          abk = one -afk
          duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
          duz =   duz                    /dzk**2
          dtz = ( thl(k)-thl(k-1) )/( dzk )
          dqz = ( qw(k)-qw(k-1) )/( dzk )

          vtt =  one +vt(k)*abk +vt(k-1)*afk  ! Beta-theta in NN09, Eq. 39
          vqq =  tv0 +vq(k)*abk +vq(k-1)*afk  ! Beta-q
          dtq =  vtt*dtz +vqq*dqz

          dtl(k) =  dtz
          dqw(k) =  dqz
          dtv(k) =  dtq

          gm(k)  =  duz
          gh(k)  = -dtq*gtr

          !   **  Gradient Richardson number  **
          ri(k) = -gh(k)/MAX( duz, 1.0e-10_kind_phys )
       enddo
       
    case (1) !form of effective static stability (ess) similar to O'Gorman (2011, JAS)
       !but modified to better fit within the MYNN-EDMF. This approach
       !destabilizes the lapse rates in grid cells with tke/mf/clouds by adding
       !a negative eqiv potential temperature gradient proportional to the
       !turbulence and cloud fractions. The coefficient for lambda (clam) has
       !units of s/m (i.e., eddy timescale / pbl depth ~ 1.8), but needs to be
       !limited over land as noted below.
       !--------------------------------
       !taper ess for hurricane conditions
       ugrid  = sqrt(u(kts)**2 + v(kts)**2)
       uonset = 15._kind_phys
       taper  = one - p95*min(one, max(zero, ugrid - uonset)/twenty) !reduce to 0.1
       if ((xland-1.5).GE.zero) then !water
          clam0 = 1.8_kind_phys * taper
       else                          !land
          !since the MF and TKE is much large over land and pblhs are much deeper,
          !we need a smaller tuning constant.
          clam0 = 0.08_kind_phys
       endif

       do k = kts+1,kte
          dzk   = p5 *  ( dz(k)+dz(k-1) )
          afk   = dz(k)/( dz(k)+dz(k-1) )
          abk   = one -afk
          duz   = ( u(k)-u(k-1) )**2 + ( v(k)-v(k-1) )**2
          duz   =   duz               / dzk**2
          dtz   = ( thl(k)-thl(k-1) ) / dzk
          dqz   = ( qw(k) -qw(k-1)  ) / dzk

          !blend clam to free tropospheric value (0.2) above the pblh
          wt    = min(one, max(zero, zw(k)-(pblh+hundred))/max(400._kind_phys, p3*pblh)) !0 below pblh, 1 above pblh
          clam  = clam0*(one-wt) + p2*wt
          !choke of clam near the surface
          clam  = clam * min(one, max(zero, zw(k)-twenty)/200.0_kind_phys)
          
          !Use a blended subgrid-cloud-included theta-v and theta-l-v for background
          !thermodynamic profile:
          thv1  = thvp*thv(k)   + (one-thvp)*thlv(k)
          thv0  = thvp*thv(k-1) + (one-thvp)*thlv(k-1)
          !Then use the equivalent potential pemerature profile to adjust the sbability in partially
          !saturated conditions with an upwards vertical component.
          !equivalent potential temperature at k
          tabs  = th(k)*exner(k)
          xl    = xl_blend(tabs)            ! blended latent heat
          qsat  = qsat_blend(tabs,p(k))     ! saturation water vapor mixing ratio
          eth1  = th(k)*exp(xl*qsat/(cp*tabs))
          !equivalent potential temperature at k-1
          tabs  = th(k-1)*exner(k-1)
          xl    = xl_blend(tabs)            ! blended latent heat
          qsat  = qsat_blend(tabs,p(k-1))   ! saturation water vapor mixing ratio
          eth0  = th(k-1)*exp(xl*qsat/(cp*tabs))

          !mass flux at interface (0 to 0.2 m/s):
          mfi   = max(zero, min(p2, (p5*(edmf_a(k)+edmf_a(k-1))) * (p5*(edmf_w(k)+edmf_w(k-1)))))

          !normalize tke to be comparable in magnitude to marine mass flux (0 to 0.1 m/s).
          !Roughly half the tke is ascending, so 0.5*tke = 0.25*qke.
          qkei  = max(qkemin, p5*(qke(k)+qke(k-1)))  !avg qke at interface
          qkei  = min(p5, sqrt(p25*qkei))*p2
          
          !control factor for cloudy grid cells and/or have nonzero mass flux
          cldfrai =  p5*(cldfra(k)+cldfra(k-1)) !avg cloud fraction at interface
          !introduce factor for limit reduction of ess in fully resolved clouds where latent heating is explicit
          cfac    = max( 0.01_kind_phys, min(2.7_kind_phys * (cldfrai - one)**2, one))
          !cfac    = max( 0.01_kind_phys, min(4.5_kind_phys * (min(0.85_kind_phys,cldfrai) - 0.85_kind_phys)**2, one))
          cldfrai = min(p1, cldfrai*cfac)*ten
          !cldfrai  = max(ncld, min(p1,mfi)) !TEST: always allow some destabilization in grid cells with plumes.
          
          !lambda significantly departs from OGorman (2011) by using MYNN-EDMF-specific information: 
          lambda= clam * max(mfi, qkei) * cldfrai
          !lambda= clam * (mfi + qkei) * cldfrai
          dtq   = ( thv1-thv0 )/dzk + lambda*min(zero, ( eth1-eth0 )/dzk )
          !dtq = max(zero, dtq)

          dtl(k)=  dtz
          dqw(k)=  dqz
          dtv(k)=  dtq

          gm(k) =  duz
          gh(k) = -dtq*gtr

          !   **  Gradient Richardson number  **
          ri(k) = -gh(k)/MAX( duz, 1.0e-10_kind_phys )
       enddo
       
    end select
    
    do k = kts+1,kte
       !a2fac is needed for the Canuto/Kitamura mod
       a2fac = one/(one + ckmod*MAX(ri(k),zero))

       rfc = g1/( g1+g2 )
       f1  = b1*( g1-c1 ) +three*a2*a2fac *( one    -c2 )*( one-c5 ) &
           &              +two  *a1*( three-two*c2 )
       f2  = b1*( g1+g2 ) -three*a1*( one      -c2 )
       rf1 = b1*( g1-c1 )/f1
       rf2 = b1*  g1     /f2
       smc = a1   /(a2*a2fac)*  f1/f2
       shc = three*(a2*a2fac)*( g1+g2 )

       ri1 = p5/smc
       ri2 = rf1*smc
       ri3 = four*rf2*smc -two*ri2
       ri4 = ri2**2

       !   **  Flux Richardson number  **
       rf = MIN( ri1*( ri(k) + ri2-SQRT(ri(k)**2 - ri3*ri(k) + ri4) ), rfc )

       sh(k) = shc*( rfc-rf )/( one-rf )
       sm(k) = smc*( rf1-rf )/( rf2-rf ) * sh(k)
    enddo

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_level2
!! @}

! ==================================================================
!     SUBROUTINE  mym_length:
!
!     Input variables:    see subroutine mym_initialize
!
!     Output variables:   see subroutine mym_initialize
!
!     Work arrays:
!       elt(single point): Length scale depending on the PBL depth    (m)
!       vsc(single point): Velocity scale q_c                       (m/s)
!                          at first, used for computing elt
!
!     NOTE: the mixing lengths are meant to be calculated at the full-
!           sigmal levels (or interfaces beween the model layers).
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the mixing lengths.
  SUBROUTINE  mym_length (                     & 
    &            kts,kte,xland,                &
    &            dz, dx, zw,                   &
    &            rmol, rmolh, flt, fltv, flq,  &
    &            vt, vq,                       &
    &            u1, v1, qke, qpe,             &
    &            dtv,                          &
    &            el,                           &
    &            pblh, theta, qkw, qtw,        &
    &            Psig_bl, cldfra_bl1,          &
    &            bl_mynn_mixlength,            &
    &            edmf_w1,edmf_a1,              &
    &            edmf_w_dd1,edmf_a_dd1         )
    
!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only:       &
         qkemin,qmin,tv0,gtr,karman,grav,zmax, &
         zero,one,two,four,five,ten,twenty,    &
         thirty,forty,fifty,hundred,           &
         p1,p2,p25,p3,p333,p4,p5,kind_phys
    
    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in)   :: bl_mynn_mixlength
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in) :: rmol,rmolh,flt,fltv,flq,Psig_bl,xland
    real(kind_phys), intent(in) :: dx,pblh
    real(kind_phys), dimension(kts:kte), intent(in)   :: u1,v1,  &
         &qke,qpe,vt,vq,cldfra_bl1,edmf_w1,edmf_a1,edmf_w_dd1,edmf_a_dd1
    real(kind_phys), dimension(kts:kte), intent(out)  :: el
    real(kind_phys), dimension(kts:kte), intent(inout):: qkw,qtw
    real(kind_phys), dimension(kts:kte), intent(in)   :: dtv
    real(kind_phys):: elt,vsc
    real(kind_phys), dimension(kts:kte), intent(in) :: theta
    real(kind_phys), dimension(kts:kte) :: qtke,elBLmin,elBLavg,thetaw
    real(kind_phys):: wt,wt2,pblh2,h1,h2,hs,elBLmin0,elBLavg0,cldavg

    ! THE FOLLOWING CONSTANTS ARE IMPORTANT FOR REGULATING THE
    ! MIXING LENGTHS:
    real(kind_phys):: cns,   &   !< for surface layer (els) in stable conditions
            alp1,            &   !< for turbulent length scale (elt)
            alp2,            &   !< for buoyancy length scale (elb)
            alp3,            &   !< for buoyancy enhancement factor of elb
            alp4,            &   !< for surface layer (els) in unstable conditions
            alp5,            &   !< for BouLac mixing length or above PBLH
            alp6                 !< for mass-flux/

    !THE FOLLOWING LIMITS DO NOT DIRECTLY AFFECT THE ACTUAL PBLH.
    !THEY ONLY IMPOSE LIMITS ON THE CALCULATION OF THE MIXING LENGTH 
    !SCALES SO THAT THE BOULAC MIXING LENGTH (IN FREE ATMOS) DOES
    !NOT ENCROACH UPON THE BOUNDARY LAYER MIXING LENGTH (els, elb & elt).
    real(kind_phys), parameter :: minpblh     = 250.  !< min mixed-layer height
    real(kind_phys), parameter :: maxdz       = 750.  !< max (half) transition layer depth
                                     !! =0.3*2500 m PBLH, so the transition
                                     !! layer stops growing for PBLHs > 2.5 km.
    real(kind_phys), parameter :: mindz       = 250.  !< min (half) transition layer depth

    !SURFACE LAYER LENGTH SCALE MODS TO REDUCE IMPACT IN UPPER BOUNDARY LAYER
    real(kind_phys), parameter :: ZSLH        = 100.  !< Max height correlated to surface conditions (m)
    real(kind_phys), parameter :: CSL         = 2.    !< CSL = constant of proportionality to L O(1)
    real(kind_phys), parameter :: qkw_elb_min = 0.05  !this assumes some minumum TKE/TPE is present 

    integer :: i,j,k
    real(kind_phys):: afk,abk,zwk,zwk1,dzk,qdz,vflx,bv,tau_cloud,      &
           & wstar,elb,els,elf,el_stab,el_mf,el_stab_mf,elb_mf,elt_max,&
           & pblh_plus_ent,uonset,ugrid,wt_u1,wt_u2,el_les,qkw_mf,     &
           & z_m,el_unstab,els1,alp3z,cpblh,wt_dx
    real(kind_phys), parameter :: ctau = 1000. !constant for tau_cloud

!    tv0 = 0.61*tref
!    gtr = 9.81/tref

    SELECT CASE(bl_mynn_mixlength)

      CASE (0) ! ORIGINAL MYNN MIXING LENGTH + BouLac

        cns  = 2.7_kind_phys
        alp1 = 0.23_kind_phys
        alp2 = one
        alp3 = five
        alp4 = hundred
        alp5 = p3

        ! Impose limits on the height integration for elt and the transition layer depth
        pblh2= pblh+dz(kts)          !originally integrated to model top, not just pblh.
        h1   = max(p3*pblh2,mindz)
        h1   = min(h1,maxdz)         ! 1/2 transition layer depth
        h2   = h1*p5                ! 1/4 transition layer depth

        qkw(kts) = SQRT(MAX(qke(kts), qkemin))
        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = one -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk, qkemin))
        END DO
        qtw = qkw

        elt = 1.0e-5_kind_phys
        vsc = 1.0e-5_kind_phys

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k   = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. pblh2+h1)
           dzk = p5*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt  = alp1*elt/vsc
        vflx = ( vt(kts)+one )*flt +( vq(kts)+tv0 )*flq
        vsc  = ( gtr*elt*MAX( vflx, zero ) )**p333

        !   **  Strictly, el(i,k=1) is not zero.  **
        el(kts) = zero
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. zero ) THEN
              bv  = SQRT( gtr*dtv(k) )
              elb = alp2*qkw(k) / bv &
                  &       *( one + alp3/alp2*&
                  &SQRT( vsc/( bv*elt ) ) )
              elf = alp2 * qkw(k)/bv

           ELSE
              elb = 1.0e10_kind_phys
              elf = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmol .GT. 0.0 ) THEN
              els  = karman*zwk/(one+cns*MIN( zwk*rmol, zmax ))
           ELSE
              els  = karman*zwk*( one - alp4* zwk*rmol )**p2
           END IF

           !   ** HARMONC AVERGING OF MIXING LENGTH SCALES:
           !       el(k) =    MIN(elb/( elb/elt+elb/els+one ),elf)
           !       el(k) =    elb/( elb/elt+elb/els+one )

           wt=p5*TANH((zwk - (pblh2+h1))/h2) + p5

           el(k) = MIN(elb/( elb/elt+elb/els+one ),elf)

        END DO

      CASE (1) !NONLOCAL (using BouLac) FORM OF MIXING LENGTH

        !wt_u* are for hurricane tuning, meant to reduce diffusion in hurricanes
        ugrid = sqrt(u1(kts)**2 + v1(kts)**2)
        uonset= 15._kind_phys
        wt_u1 = one - p2*min(one, max(zero, ugrid - uonset)/forty) !reduce to 0.8
        wt_u2 = one - p5*min(one, max(zero, ugrid - uonset)/forty) !reduce to 0.5
        !scale-awareness for the mesoscale greyzone (4-16 km)
        wt_dx = one - min(one, (max(zero, dx-4000._kind_phys)/12000._kind_phys))
        
        cns   = 3.5_kind_phys
        alp1  = 0.23_kind_phys*wt_dx + (one-wt_dx)*0.24_kind_phys
        alp2  = p3*wt_dx             + (one-wt_dx)*0.40_kind_phys
        alp3  = five * wt_u2 !taper off bouyancy enhancement in shear-driven pbls
        if ((xland-1.5).GE.zero) then !hurricane tuning, over water only
           alp4  = twenty * wt_u2
        else
           alp4  = 15.0_kind_phys
        endif
        alp5  = p3
        alp6  = fifty

        ! Impose limits on the height integration for elt and the transition layer depth
        pblh2 = max(pblh,200._kind_phys)
        h1    = max(p3*pblh2,300._kind_phys)
        h1    = min(h1,600._kind_phys)   ! 1/2 transition layer depth
        h2    = h1*p5                    ! 1/4 transition layer depth

        qtke(kts)  =max(p5*qke(kts), p5*qkemin) !tke at full sigma levels
        thetaw(kts)=theta(kts)            !theta at full-sigma levels
        qkw(kts)   =sqrt(max(qke(kts)         , qkemin))
        qtw(kts)   =sqrt(max(qke(kts)+qpe(kts), qkemin))
        DO k = kts+1,kte
           afk      = dz(k)/( dz(k)+dz(k-1) )
           abk      = one -afk
           qkw(k)   = sqrt(max(qke(k)*abk+qke(k-1)*afk, qkemin))
           qtw(k)   = sqrt(max((qke(k)+qpe(k))*abk+(qke(k-1)+qpe(k-1))*afk, qkemin))
           qtke(k)  = max(p5*(qkw(k)**2), 0.005_kind_phys) ! q -> TKE
           thetaw(k)= theta(k)*abk + theta(k-1)*afk
        END DO

        elt = 1.0e-5_kind_phys
        vsc = 1.0e-5_kind_phys

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k   = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. pblh2+max(200.,min(p3*pblh2,600._kind_phys)))
           dzk = p5*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k), 0.01_kind_phys ), 30.0_kind_phys)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        if ((xland-1.5).GE.zero) then !hurricane tuning, over water only
           elt_max=350.+100.*min(one, max(zero, ugrid - twenty)/thirty)
        else
           elt_max=400._kind_phys
        endif
        elt = MIN( MAX( alp1*elt/vsc, 8._kind_phys), elt_max)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+one )*flt + ( vq(kts)+tv0 )*flq
        vflx= fltv-0.005_kind_phys
        vsc = ( gtr*elt*MAX( vflx, zero ) )**p333

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = zero
        zwk1    = zw(kts+1)         !full-sigma levels

        ! COMPUTE BouLac mixing length
        CALL boulac_length(kts,kte,zw,dz,qtke,thetaw,elBLmin,elBLavg)

        DO k = kts+1,kte
           zwk    = zw(k)          !full-sigma levels
           qkw_mf = max((p5*((edmf_a1(k)+edmf_a1(k-1))))*(p5*((edmf_w1(k)+edmf_w1(k-1)))), &
                  & abs(edmf_a_dd1(k-1)*edmf_w_dd1(k-1)))
           !pblh-dependent modifier
           cpblh  = min((zwk + p25*pblh2)/pblh2, one)
           qkw_mf = cpblh * qkw_mf
           alp3z  = cpblh * alp3
           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv     = max( sqrt( gtr*dtv(k) ), 0.001_kind_phys)
              elb_mf = alp6*qkw_mf/bv
              elb_mf = elb_mf / (one + (elb_mf/100._kind_phys))
              elb    = alp2*(max(qkw(k),qkw_elb_min))/bv    &
                   &  *( one + alp3z*SQRT( vsc/(bv*elt) ) )
              elb    = max(elb, elb_mf)
              elb    = MIN(elb, zwk)
              elf    = one * max(qkw(k), qkw_elb_min)/bv
              elf    = max(elf, elb_mf)
              elBLavg(k) = MAX(alp5*elBLavg(k), elb_mf)
           ELSE
              elb        = 1.0e10_kind_phys
              elf        = elb
              elBLavg(k) = p4*elBLavg(k)
           ENDIF

           z_m = MAX(karman, zwk - four)
           !   **  Length scale in the surface layer  **
           IF ( rmol .GT. 0.0 ) THEN
              els  = karman*zwk/(one+cns*MIN( zwk*rmol, zmax ))
              els1 = karman*z_m/(one+cns*MIN( zwk*rmol, zmax ))
           ELSE
              els  = karman*zwk*( one - alp4* zwk*rmol)**p2
              els1 = karman*z_m*( one - alp4* zwk*rmol)**p2
           END IF

           !   ** NOW BLEND THE MIXING LENGTH SCALES:
           wt =p5*TANH((zwk - (pblh2+h1))/h2) + p5
           !add blending to use BouLac mixing length in free atmos;
           !defined relative to the PBLH (pblh) + transition layer (h1)
           !el(k) = MIN(elb/( elb/elt+elb/els+one ),elf)
           !squared blending - but take out elb (makes it underdiffusive)
           !el_unstab = SQRT( els**2/(one + (els**2/elt**2) +(els**2/elb**2)))
           !el_unstab = sqrt( els**2/(one + (els**2/elt**2)))
           !non-squared blending:
           el_unstab = els/(one + (els1/elt))
           el(k)     = min(el_unstab, elb)
           el(k)     = min(el(k), elf)  !elf can be smaller than elb in upper pbl
           if ((xland-1.5).GE.zero) then !hurricane tuning, over water only
              el(k)  = el(k)*wt_u1
           endif
           el(k)     = el(k)*(one-wt) + elBLavg(k)*wt

           !if (el(k) > 1000.) then
           !   print*,"big ML at k=",k," el=",el(k),"elb=",elb," elBL=",elBLavg(k), &
           !   " elt=",elt," els=",els," N=",bv," qtke=",qtke(k), &
           !   " pblh=",pblh2," zwk=",zwk," wt=",wt
           !endif

        END DO

     CASE (2) !Local (mostly) mixing length formulation

        Uonset = 3.5_kind_phys + dz(kts)*p1
        Ugrid  = sqrt(u1(kts)**2 + v1(kts)**2)
        cns    = 3.5_kind_phys !  * (one - MIN(MAX(Ugrid - Uonset, 0.0)/10.0, 1.0))
        alp1   = 0.22_kind_phys
        alp2   = p3
        alp3   = 2.5_kind_phys
        alp4   = five
        alp5   = alp2  !like alp2, but for free atmosphere
        alp6   = fifty !used for MF mixing length

        ! Impose limits on the height integration for elt and the transition layer depth
        !pblh2=MAX(pblh,minpblh)
        pblh2=MAX(pblh, 200._kind_phys)
        !h1=MAX(0.3*pblh2,mindz)
        !h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h1=MAX(p3*pblh2,300._kind_phys)
        h1=MIN(h1,600._kind_phys)
        h2=h1*p5                ! 1/4 transition layer depth

        qtke(kts)=MAX(p5*qke(kts), p5*qkemin) !tke at full sigma levels
        qkw(kts) = SQRT(MAX(qke(kts)         , qkemin))
        qtw(kts) = SQRT(MAX(qke(kts)+qpe(kts), qkemin))
        DO k = kts+1,kte
           afk    = dz(k)/( dz(k)+dz(k-1) )
           abk    = one -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk                    , qkemin))
           qtw(k) = SQRT(MAX((qke(k)+qpe(k))*abk+(qke(k-1)+qpe(k-1))*afk, qkemin))
           qtke(k)= p5*qkw(k)**2  ! qkw -> TKE
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        PBLH_PLUS_ENT = MAX(pblh+h1, 100.)
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. PBLH_PLUS_ENT)
           dzk = p5*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.03 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt = MIN( MAX(alp1*elt/vsc, ten), 400.)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+one )*flt +( vq(kts)+tv0 )*flq
        vflx = fltv
        vsc  = ( gtr*elt*MAX( vflx, zero ) )**p333

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = zero
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels
           dzk = p5*( dz(k)+dz(k-1) )
           cldavg = p5*(cldfra_bl1(k-1)+cldfra_bl1(k))
           qkw_mf = max((p5*(edmf_a1(k)+edmf_a1(k-1)))*(p5*(edmf_w1(k)+edmf_w1(k-1))), &
                  & abs(edmf_a_dd1(k-1)*edmf_w_dd1(k-1)))

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              !impose min value on bv
              bv  = MAX( SQRT( gtr*dtv(k) ), 0.001)  
              !elb_mf = alp2*qkw(k) / bv  &
              elb_mf = MAX(alp2*qkw(k),                    &
                  &        alp6*qkw_mf) / bv               &
                  &  *( one + alp3*SQRT( vsc/( bv*elt ) ) )
              elb = MIN(MAX(alp5*qkw(k), alp6*qkw_mf)/bv, zwk)

              !tau_cloud = MIN(MAX(0.5*pblh/((gtr*pblh*MAX(vflx,1.0e-4))**p333),30.),150.)
              wstar = 1.25*(gtr*pblh*MAX(vflx,1.0e-4))**p333
              tau_cloud = MIN(MAX(ctau * wstar/grav, 30.), 150.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=p5*TANH((zwk - (pblh2+h1))/h2) + p5
              tau_cloud = tau_cloud*(1.-wt) + fifty*wt
              elf = MIN(MAX(tau_cloud*SQRT(MIN(qtke(k),forty)), &
                  &         alp6*qkw_mf/bv), zwk)

              !IF (zwk > pblh .AND. elf > 400.) THEN
              !   ! COMPUTE BouLac mixing length
              !   !CALL boulac_length0(k,kts,kte,zw,dz,qtke,thetaw,elBLmin0,elBLavg0)
              !   !elf = alp5*elBLavg0
              !   elf = MIN(MAX(50.*SQRT(qtke(k)), 400.), zwk)
              !ENDIF

           ELSE
              ! use version in development for RAP/HRRR 2016
              ! JAYMES-
              ! tau_cloud is an eddy turnover timescale;
              ! see Teixeira and Cheinet (2004), Eq. 1, and
              ! Cheinet and Teixeira (2003), Eq. 7.  The
              ! coefficient 0.5 is tuneable. Expression in
              ! denominator is identical to vsc (a convective
              ! velocity scale), except that elt is relpaced
              ! by pblh, and zero is replaced by 1.0e-4 to
              ! prevent division by zero.
              !tau_cloud = MIN(MAX(0.5*pblh/((gtr*pblh*MAX(vflx,1.0e-4))**p333),50.),150.)
              wstar     = 1.25*(gtr*pblh*MAX(vflx,1.0e-4))**p333
              tau_cloud = MIN(MAX(ctau * wstar/grav, fifty), 200.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt        = p5*TANH((zwk - (pblh2+h1))/h2) + p5
              !tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              tau_cloud = tau_cloud*(one-wt) + MAX(hundred,dzk*p25)*wt

              elb       = MIN(tau_cloud*SQRT(MIN(qtke(k),40.)), zwk)
              !elf = elb
              elf       = elb !/(1. + (elb/800.))  !bound free-atmos mixing length to < 800 m.
              elb_mf    = elb
         END IF
         elf    = elf/(one + (elf/800.))  !bound free-atmos mixing length to < 800 m.
         elb_mf = MAX(elb_mf, 0.01) !to avoid divide-by-zero below

         !   **  Length scale in the surface layer  **
         IF ( rmol .GT. 0.0 ) THEN
            els  = karman*zwk/(one+cns*MIN( zwk*rmol, zmax ))
         ELSE
            els  =  karman*zwk*( one - alp4* zwk*rmol)**p2
         END IF

         !   ** NOW BLEND THE MIXING LENGTH SCALES:
         wt=p5*TANH((zwk - (pblh2+h1))/h2) + p5

         !try squared-blending
         el(k) = SQRT( els**2/(one + (els**2/elt**2) +(els**2/elb_mf**2)))
         el(k) = el(k)*(one-wt) + elf*wt
       END DO

    END SELECT

    ! include scale-awareness. limit el to be < 0.25*dz)
    DO k = kts+1,kte
       el_les = p25*p5*( dz(k)+dz(k-1) )
       el(k)  = el(k)*Psig_bl + (one-Psig_bl)*min(el_les,el(k))
    ENDDO
      
#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_length

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine was taken from the BouLac scheme in WRF-ARW and modified for
!! integration into the MYNN PBL scheme. WHILE loops were added to reduce the
!! computational expense. This subroutine computes the length scales up and down
!! and then computes the min, average of the up/down length scales, and also
!! considers the distance to the surface.
!\param dlu  the distance a parcel can be lifted upwards give a finite
!  amount of TKE.
!\param dld  the distance a parcel can be displaced downwards given a
!  finite amount of TKE.
!\param lb1  the minimum of the length up and length down
!\param lb2  the average of the length up and length down
  SUBROUTINE boulac_length0(k,kts,kte,zw,dz,qtke,theta,lb1,lb2)
!
!    NOTE: This subroutine was taken from the BouLac scheme in WRF-ARW
!          and modified for integration into the MYNN PBL scheme.
!          WHILE loops were added to reduce the computational expense.
!          This subroutine computes the length scales up and down
!          and then computes the min, average of the up/down
!          length scales, and also considers the distance to the
!          surface.
!
!      dlu = the distance a parcel can be lifted upwards give a finite
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: gtr,zero,two,p1,p5,kind_phys
    
     integer, intent(in) :: k,kts,kte
     real(kind_phys), dimension(kts:kte), intent(in) :: qtke,dz,theta
     real(kind_phys), intent(out) :: lb1,lb2
     real(kind_phys), dimension(kts:kte+1), intent(in) :: zw

     !LOCAL VARS
     integer :: izz, found
     real(kind_phys):: dlu,dld
     real(kind_phys):: dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz


     !----------------------------------
     ! FIND DISTANCE UPWARD             
     !----------------------------------
     zup=zero
     dlu=zw(kte+1)-zw(k)-dz(k)*p5
     zzz=zero
     zup_inf=zero
     beta=gtr           !Buoyancy coefficient (g/tref)

     !print*,"FINDING Dup, k=",k," zw=",zw(k)

     if (k .lt. kte) then      !cant integrate upwards from highest level
        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .lt. kte) then
              dzt=dz(izz)                   ! layer depth above
              zup=zup-beta*theta(k)*dzt     ! initial PE the parcel has at k
              !print*,"  ",k,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt*p5 ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer k to izz+1
              !print*,"  PE=",zup," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zup .and. qtke(k).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. zero) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(k)) + &
                      & sqrt( max(zero,(beta*(theta(izz)-theta(k)))**2 + &
                      &       two*bbb*beta*(qtke(k)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k))then
                       tl=(qtke(k)-zup_inf)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=zero
                    endif
                 endif
                 dlu=zzz-dzt+tl
                 !print*,"  FOUND Dup:",dlu," z=",zw(izz)," tl=",tl
                 found =1
              endif
              zup_inf=zup
              izz=izz+1
           ELSE
              found = 1
           ENDIF

        ENDDO

     endif

     !----------------------------------
     ! FIND DISTANCE DOWN               
     !----------------------------------
     zdo=zero
     zdo_sup=zero
     dld=zw(k)
     zzz=zero

     !print*,"FINDING Ddown, k=",k," zwk=",zw(k)
     if (k .gt. kts) then  !cant integrate downwards from lowest level

        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(k)*dzt
              !print*,"  ",k,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt*p5
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zdo .and. qtke(k).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. zero) then
                    tl=(beta*(theta(izz)-theta(k))+ &
                      & sqrt( max(zero,(beta*(theta(izz)-theta(k)))**2 + &
                      &       two*bbb*beta*(qtke(k)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k)) then
                       tl=(qtke(k)-zdo_sup)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=zero
                    endif
                 endif
                 dld=zzz-dzt+tl
                 !print*,"  FOUND Ddown:",dld," z=",zw(izz)," tl=",tl
                 found = 1
              endif
              zdo_sup=zdo
              izz=izz-1
           ELSE
              found = 1
           ENDIF
        ENDDO

     endif

     !----------------------------------
     ! GET MINIMUM (OR AVERAGE)         
     !----------------------------------
     !The surface layer length scale can exceed z for large z/L,
     !so keep maximum distance down > z.
     dld = min(dld,zw(k+1))!not used in PBL anyway, only free atmos
     lb1 = min(dlu,dld)     !minimum
     !JOE-fight floating point errors
     dlu=MAX(p1,MIN(dlu,1000._kind_phys))
     dld=MAX(p1,MIN(dld,1000._kind_phys))
     lb2 = sqrt(dlu*dld)    !average - biased towards smallest
     !lb2 = 0.5*(dlu+dld)   !average

     if (k .eq. kte) then
        lb1 = zero
        lb2 = zero
     endif
     !print*,"IN MYNN-BouLac",k,lb1
     !print*,"IN MYNN-BouLac",k,dld,dlu

  END SUBROUTINE boulac_length0

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine was taken from the BouLac scheme in WRF-ARW
!! and modified for integration into the MYNN PBL scheme.
!! WHILE loops were added to reduce the computational expense.
!! This subroutine computes the length scales up and down
!! and then computes the min, average of the up/down
!! length scales, and also considers the distance to the
!! surface.
  SUBROUTINE boulac_length(kts,kte,zw,dz,qtke,theta,lb1,lb2)
!      dlu = the distance a parcel can be lifted upwards give a finite 
!            amount of TKE.
!      dld = the distance a parcel can be displaced downwards given a
!            finite amount of TKE.
!      lb1 = the minimum of the length up and length down
!      lb2 = the average of the length up and length down
!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: gtr,zero,one,two,p1,p5,kind_phys
    
     integer, intent(in) :: kts,kte
     real(kind_phys), dimension(kts:kte),   intent(in) :: qtke,dz,theta
     real(kind_phys), dimension(kts:kte),   intent(out):: lb1,lb2
     real(kind_phys), dimension(kts:kte+1), intent(in) :: zw

     !LOCAL VARS
     integer :: iz, izz, found
     real(kind_phys), dimension(kts:kte) :: dlu,dld
     real(kind_phys), parameter :: Lmax=1500.  !soft limit
     real(kind_phys):: dzt, zup, beta, zup_inf, bbb, tl, zdo, zdo_sup, zzz

     !print*,"IN MYNN-BouLac",kts, kte

     do iz=kts,kte

        !----------------------------------
        ! FIND DISTANCE UPWARD
        !----------------------------------
        zup    =zero
        dlu(iz)=zw(kte+1)-zw(iz)-dz(iz)*p5
        zzz    =zero
        zup_inf=zero
        beta   =gtr           !Buoyancy coefficient (g/tref)

        !print*,"FINDING Dup, k=",iz," zw=",zw(iz)

        if (iz .lt. kte) then      !cant integrate upwards from highest level

          found = 0
          izz=iz
          DO WHILE (found .EQ. 0)

            if (izz .lt. kte) then
              dzt=dz(izz)                    ! layer depth above
              zup=zup-beta*theta(iz)*dzt     ! initial PE the parcel has at iz
              !print*,"  ",iz,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt*p5 ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer iz to izz+1
              !print*,"  PE=",zup," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zup .and. qtke(iz).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(iz)) + &
                      & sqrt( max(zero,(beta*(theta(izz)-theta(iz)))**2 + &
                      &       two*bbb*beta*(qtke(iz)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz))then
                       tl=(qtke(iz)-zup_inf)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=zero
                    endif
                 endif            
                 dlu(iz)=zzz-dzt+tl
                 !print*,"  FOUND Dup:",dlu(iz)," z=",zw(izz)," tl=",tl
                 found =1
              endif
              zup_inf=zup
              izz=izz+1
             ELSE
              found = 1
            ENDIF

          ENDDO

        endif
                   
        !----------------------------------
        ! FIND DISTANCE DOWN
        !----------------------------------
        zdo    =zero
        zdo_sup=zero
        dld(iz)=zw(iz)
        zzz    =zero

        !print*,"FINDING Ddown, k=",iz," zwk=",zw(iz)
        if (iz .gt. kts) then  !cant integrate downwards from lowest level

          found = 0
          izz=iz       
          DO WHILE (found .EQ. 0) 

            if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(iz)*dzt
              !print*,"  ",iz,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt*p5
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zdo .and. qtke(iz).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(iz))+ &
                      & sqrt( max(zero,(beta*(theta(izz)-theta(iz)))**2 + &
                      &       two*bbb*beta*(qtke(iz)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz)) then
                       tl=(qtke(iz)-zdo_sup)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=zero
                    endif
                 endif            
                 dld(iz)=zzz-dzt+tl
                 !print*,"  FOUND Ddown:",dld(iz)," z=",zw(izz)," tl=",tl
                 found = 1
              endif
              zdo_sup=zdo
              izz=izz-1
            ELSE
              found = 1
            ENDIF
          ENDDO

        endif

        !----------------------------------
        ! GET MINIMUM (OR AVERAGE)
        !----------------------------------
        !The length scale down should not exceed z.
        dld(iz) = min(dld(iz),zw(iz+1))
        !apply soft upper limit (only impacts very large lb; lb=100 by 5%, lb=500 by 20%)
        dlu(iz)=MAX(p1, dlu(iz)/(one + (dlu(iz)/Lmax)) )
        dld(iz)=MAX(p1, dld(iz)/(one + (dld(iz)/Lmax)) )

        !minimum of up/down length scales
        lb1(iz) = min(dlu(iz),dld(iz))
        !average - biased towards smallest
        lb2(iz) = sqrt(dlu(iz)*dld(iz))

        !print*,"IN MYNN-BouLac",kts, kte,lb1(iz)
        !print*,"IN MYNN-BouLac",iz,dld(iz),dlu(iz)

     ENDDO

     lb1(kte) = lb1(kte-1)
     lb2(kte) = lb2(kte-1)
        
  END SUBROUTINE boulac_length
!
! ==================================================================
!     SUBROUTINE  mym_turbulence:
!
!     Input variables:    see subroutine mym_initialize
!       closure        : closure level (2.5, 2.6, or 3.0)
!
!     # ql, vt, vq, qke, tsq, qsq and cov are changed to input variables.
!
!     Output variables:   see subroutine mym_initialize
!       dfm(nx,nz,ny) : Diffusivity coefficient for momentum,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfh(nx,nz,ny) : Diffusivity coefficient for heat,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       dfq(nx,nz,ny) : Diffusivity coefficient for q^2,
!                         divided by dz (not dz*h(i,j))            (m/s)
!       tcd(nx,nz,ny)   : Countergradient diffusion term for Theta_l
!                                                                  (K/s)
!       qcd(nx,nz,ny)   : Countergradient diffusion term for Q_w
!                                                              (kg/kg s)
!       pd?(nx,nz,ny) : Half of the production terms
!
!       Only tcd and qcd are defined at the center of the grid boxes
!
!     # DO NOT forget that tcd and qcd are added on the right-hand side
!       of the equations for Theta_l and Q_w, respectively.
!
!     Work arrays:        see subroutine mym_initialize and level2
!
!     # dtl, dqw, dtv, gm and gh are allowed to share storage units with
!       dfm, dfh, dfq, tcd and qcd, respectively, for saving memory.
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the vertical diffusivity coefficients and the 
!! production terms for the turbulent quantities.      
!>\section gen_mym_turbulence GSD mym_turbulence General Algorithm
!! Two subroutines mym_level2() and mym_length() are called within this
!!subrouine to collect variable to carry out successive calculations:
!! - mym_level2() calculates the level 2 nondimensional wind shear \f$G_M\f$
!! and vertical temperature gradient \f$G_H\f$ as well as the level 2 stability
!! functions \f$S_h\f$ and \f$S_m\f$.
!! - mym_length() calculates the mixing lengths.
!! - The stability criteria from Helfand and Labraga (1989) are applied.
!! - The stability functions for level 2.5 or level 3.0 are calculated.
!! - If level 3.0 is used, counter-gradient terms are calculated.
!! - Production terms of TKE,\f$\theta^{'2}\f$,\f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$
!! are calculated.
!! - Eddy diffusivity \f$K_h\f$ and eddy viscosity \f$K_m\f$ are calculated.
!! - TKE budget terms are calculated (if the namelist parameter \p tke_budget 
!! is set to True)
  SUBROUTINE  mym_turbulence (                                &
    &            kts,kte,                                     &
    &            xland,closure,ckmod,taue,                    &
    &            dz, dx, zw, p, exner,                        &
    &            u, v, thl, thv, thlv, ql, qw,                &
    &            qke, tsq, qsq, cov,                          &
    &            vt, vq,                                      &
    &            rmol, rmolh, flt, fltv, flq,                 &
    &            pblh,theta,                                  &
    &            sh, sm,                                      &
    &            El,                                          &
    &            Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
    &		 qWT1,qSHEAR1,qBUOY1,qDISS1,                  &
    &            tke_budget,                                  &
    &            Psig_bl,Psig_shcu,cldfra_bl1,                &
    &            bl_mynn_mixlength,                           &
    &            bl_mynn_ess,                                 &
    &            edmf_w1,edmf_a1,                             &
    &            edmf_w_dd1,edmf_a_dd1,                       &
    &            TKEprod_dn,TKEprod_up,                       &
    &            spp_pbl,pattern_spp_pbl1                     )

!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: gtr,tv0,             &
         a1,a2,b1,b2,c1,e1c,e2c,e3c,e4c,e5c,cc3,              &
         zero,one,two,three,four,seven,nine,p1,p2,p5,         &
         kind_phys
    
    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in)               :: bl_mynn_mixlength
    integer, intent(in)               :: tke_budget
    integer, intent(in)               :: bl_mynn_ess
    real(kind_phys), intent(in)       :: closure,ckmod,taue
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in)       :: rmol,rmolh,flt,fltv,flq,          &
         &Psig_bl,Psig_shcu,xland,dx,pblh
    real(kind_phys), dimension(kts:kte), intent(in) :: u,v,thl,thv,        &
         &thlv,qw,ql,vt,vq,qke,tsq,qsq,cov,cldfra_bl1,edmf_w1,edmf_a1,     &
         &edmf_w_dd1,edmf_a_dd1,TKEprod_dn,TKEprod_up,p,exner

    real(kind_phys), dimension(kts:kte), intent(out) :: dfm,dfh,dfq,       &
         &pdk,pdt,pdq,pdc,tcd,qcd,el

    real(kind_phys), dimension(kts:kte), intent(inout), optional ::        &
         qWT1,qSHEAR1,qBUOY1,qDISS1
    real(kind_phys):: q3sq_old,qWTP_old,qWTP_new
    real(kind_phys):: dudz,dvdz,dTdz,upwp,vpwp,Tpwp

    real(kind_phys), dimension(kts:kte) :: qkw,qtw,dtl,dqw,dtv,gm,gh,sm,sh,qpe

    integer :: k
!    real(kind_phys):: cc2,cc3,e1c,e2c,e3c,e4c,e5c
    real(kind_phys):: e6c,dzk,afk,abk,vtt,vqq,                             &
         &cw25,clow,cupp,gamt,gamq,smd,gamv,elq,elh,elqt

    real(kind_phys):: cldavg
    real(kind_phys), dimension(kts:kte), intent(in) :: theta

    real(kind_phys)::  a2fac, duz, ri !JOE-Canuto/Kitamura mod

    real(kind_phys):: auh,aum,adh,adm,aeh,aem,Req,Rsl,Rsl2,                &
           gmelq,sm20,sh20,sm25max,sh25max,sm25min,sh25min,                &
           sm_pbl,sh_pbl,pblh2,wt,slht,mfmax,cpblh

    DOUBLE PRECISION  q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel
    DOUBLE PRECISION  q3sq, t3sq, r3sq, c3sq, dlsq, qdiv
    DOUBLE PRECISION  e1, e2, e3, e4, enum, eden, wden

!   Stochastic
    integer,         intent(in)                   :: spp_pbl
    real(kind_phys), dimension(kts:kte)           :: pattern_spp_pbl1
    real(kind_phys):: shb, Prlim
    real(kind_phys), parameter :: Prlimit_fre = 6.0 !Pr limit in free troposphere
    real(kind_phys), parameter :: Prlimit_sfc = 4.0 !Pr limit at the surface

!
!    tv0 = 0.61*tref
!    gtr = 9.81/tref
!
!    cc2 =  1.0-c2
!    cc3 =  1.0-c3
!    e1c =  3.0*a2*b2*cc3
!    e2c =  9.0*a1*a2*cc2
!    e3c =  9.0*a2*a2*cc2*( 1.0-c5 )
!    e4c = 12.0*a1*a2*cc2
!    e5c =  6.0*a1*a1

    qkw  = zero
    qtw  = zero
    pblh2= max(pblh, 200._kind_phys)
    if (closure .lt. 2.7) then
       qpe = zero
    else
       do k = kts,kte
          !Calculate qpe (twice tpe) as in Machulskaya and Mironov (2020, BLM), but with magnitudes
          !tapered at low-levels to limit high 10-m wind speed biases:
          cpblh  = min((zw(k) + p1*pblh2)/(0.8_kind_phys*pblh2), one)
          !qpe(k) = two * tsq(k) * gtr**2 * (p5*(el(k+1)+el(k))/10.)**2/(p5*max(0.01,qke(k)))
          qpe(k) = min(three, max(zero, cpblh * two * tsq(k) * gtr**2 * taue**2))
       enddo
    endif

    CALL mym_level2 (kts,kte,                   &
    &            ckmod, bl_mynn_ess,            &
    &            zw, dz, xland, pblh,           &
    &            u, v, thl, thv, thlv,          &
    &            theta, p, exner,               &
    &            qke, cldfra_bl1,               &
    &            edmf_a1, edmf_w1,   	      	&
    &            qw, ql, vt, vq,                &
    &            dtl, dqw, dtv, gm, gh, sm, sh  )
!
    CALL mym_length (                           &
    &            kts,kte,xland,                 &
    &            dz, dx, zw,                    &
    &            rmol, rmolh, flt, fltv, flq,   &
    &            vt, vq,                        &
    &            u, v, qke, qpe,                &
    &            dtv,                           &
    &            el,                            &
    &            pblh, theta,                   &
    &            qkw, qtw,                      &
    &            Psig_bl, cldfra_bl1,           &
    &            bl_mynn_mixlength,             &
    &            edmf_w1,edmf_a1,               &
    &            edmf_w_dd1,edmf_a_dd1          )
!

    DO k = kts+1,kte
       dzk  = p5 *( dz(k)+dz(k-1) )
       afk  = dz(k)/( dz(k)+dz(k-1) )
       abk  = one -afk
       elsq = el(k)**2
       q3sq = qkw(k)**2
       q2sq = b1*elsq*( sm(k)*gm(k)+sh(k)*gh(k) )

       sh20 = MAX(sh(k), 1e-5_kind_phys)
       sm20 = MAX(sm(k), 1e-5_kind_phys)
       sh(k)= MAX(sh(k), 1e-5_kind_phys)

       !Canuto/Kitamura mod
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       !   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10_kind_phys )
       !a2fac is needed for the Canuto-Kitamura modification (CKmod)
       a2fac = one/(one + ckmod*max(ri,zero))

       !The form of Zilitinkevich et al. (2006) but modified
       !half-way towards Esau and Grachev (2007, Wind Eng)
       !Prlim = MIN(0.76_kind_phys + 3.0_kind_phys*MAX(ri,zero), Prlimit_fre)
       !Prlim = MIN(0.76_kind_phys + 4.0_kind_phys*MAX(ri,zero), Prlimit_fre)

       !TZC - Kondo Correction
       if (ri >= one) then
          ! Kh/Km = 1/(7*Ri)
          Prlim = min(seven*ri, Prlimit_fre)
       elseif (ri >= 0.01 .and. ri <= one) then
          ! Kh/Km(i,k) = 1/(6.873*Ri + 1/(6.873*Ri))
          Prlim = min(6.873_kind_phys*ri + one/(6.873_kind_phys*ri), Prlimit_fre)
       else
          ! no Pr limit required?
          Prlim = Prlimit_fre
       end if
       wt    = min(one, zw(k)/max(pblh, 200._kind_phys))  !0 at sfc, 1 at/above pblh
       Prlim = wt*Prlim + (one-wt)*min(Prlimit_sfc, Prlim)
!     
!  Modified: Dec/22/2005, from here, (dlsq -> elsq)
       gmel = gm(k)*elsq
       ghel = gh(k)*elsq
!  Modified: Dec/22/2005, up to here

       ! Level 2.0 debug prints
       IF ( debug_code ) THEN
         IF (sh(k)<0.0 .OR. sm(k)<0.0) THEN
           print*,"MYNN; mym_turbulence 2.0; sh=",sh(k)," k=",k
           print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
           print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)," ri=",ri
           print*," PBLH=",pblh," u=",u(k)," v=",v(k)
         ENDIF
       ENDIF

!     **  Since qkw is set to more than 0.0, q3sq > 0.0.  **

!     new stability criteria in level 2.5 (as well as level 3) - little/no impact
!     **  Limitation on q, instead of L/q  **
       dlsq =  elsq
       IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)

       IF ( q3sq .LT. q2sq ) THEN
          !Apply Helfand & Labraga mod
          qdiv = SQRT( q3sq/q2sq )   !HL89: (1-alfa)
!
          !Use level 2.5 stability functions
          !e1   = q3sq - e1c*ghel*a2fac
          !e2   = q3sq - e2c*ghel*a2fac
          !e3   = e1   + e3c*ghel*a2fac**2
          !e4   = e1   - e4c*ghel*a2fac
          !eden = e2*e4 + e3*e5c*gmel
          !eden = MAX( eden, 1.0d-20 )
          !sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
          !sm(k) = sm(k) * qdiv

          !Use level 2.0 functions as in original MYNN
          sh(k) = sh(k) * qdiv
          sm(k) = sm(k) * qdiv

          !Recalculate terms for later use
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel * qdiv**2
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = e1   + e3c*ghel * qdiv**2
          !e4   = e1   - e4c*ghel * qdiv**2
          e1   = q3sq - e1c*ghel*a2fac * qdiv**2
          e2   = q3sq - e2c*ghel*a2fac * qdiv**2
          e3   = e1   + e3c*ghel*a2fac**2 * qdiv**2
          e4   = e1   - e4c*ghel*a2fac * qdiv**2
          eden = e2*e4 + e3*e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
       ELSE
          !JOE-Canuto/Kitamura mod
          !e1   = q3sq - e1c*ghel
          !e2   = q3sq - e2c*ghel
          !e3   = e1   + e3c*ghel
          !e4   = e1   - e4c*ghel
          e1   = q3sq - e1c*ghel*a2fac
          e2   = q3sq - e2c*ghel*a2fac
          e3   = e1   + e3c*ghel*a2fac**2
          e4   = e1   - e4c*ghel*a2fac
          eden = e2*e4 + e3*e5c*gmel
          eden = MAX( eden, 1.0d-20 )

          qdiv = one
          !Use level 2.5 stability functions
          sm(k) = q3sq*a1*( e3-three*c1*e4       )/eden
          !  sm_pbl = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          sh(k) = q3sq*(a2*a2fac)*( e2+three*c1*e5c*gmel )/eden
       END IF !end Helfand & Labraga check

       !Impose broad limits on Sh and Sm:
       gmelq    = max(real(gmel/q3sq, kind_phys), 1e-8_kind_phys)
       sm25max  = four  !MIN(sm20*3.0, SQRT(.1936/gmelq))
       sh25max  = four  !MIN(sh20*3.0, 0.76*b2)
       sm25min  = zero !MAX(sm20*0.1, 1e-6)
       sh25min  = zero !MAX(sh20*0.1, 1e-6)

       !JOE: Level 2.5 debug prints
       ! HL88 , lev2.5 criteria from eqs. 3.17, 3.19, & 3.20
       IF ( debug_code ) THEN
         IF ((sh(k)<sh25min .OR. sm(k)<sm25min .OR. &
              sh(k)>sh25max .OR. sm(k)>sm25max) ) THEN
           print*,"In mym_turbulence 2.5: k=",k
           print*," sm=",sm(k)," sh=",sh(k)
           print*," ri=",ri," Pr=",sm(k)/max(sh(k),1e-8_kind_phys)
           print*," gm=",gm(k)," gh=",gh(k)
           print*," q2sq=",q2sq," q3sq=",q3sq, q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)
           print*," PBLH=",pblh," u=",u(k)," v=",v(k)
           print*," SMnum=",q3sq*a1*( e3-three*c1*e4)," SMdenom=",eden
           print*," SHnum=",q3sq*(a2*a2fac)*( e2+three*c1*e5c*gmel ),&
                  " SHdenom=",eden
         ENDIF
       ENDIF

       !Enforce constraints for level 2.5 functions
       IF ( sh(k) > sh25max ) sh(k) = sh25max
       IF ( sh(k) < sh25min ) sh(k) = sh25min
       !IF ( sm(k) > sm25max ) sm(k) = sm25max
       !IF ( sm(k) < sm25min ) sm(k) = sm25min

       shb   = max(sh(k), 0.004_kind_phys)
       sm(k) = min(sm(k), Prlim*shb)

!   **  Level 3 : start  **
       IF ( closure .GE. 3.0 ) THEN
          t2sq = qdiv*b2*elsq*sh(k)*dtl(k)**2
          r2sq = qdiv*b2*elsq*sh(k)*dqw(k)**2
          c2sq = qdiv*b2*elsq*sh(k)*dtl(k)*dqw(k)
          t3sq = MAX( tsq(k)*abk+tsq(k-1)*afk, zero )
          r3sq = MAX( qsq(k)*abk+qsq(k-1)*afk, zero )
          c3sq =      cov(k)*abk+cov(k-1)*afk

!  Modified: Dec/22/2005, from here
          c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
!
          vtt  = one +vt(k)*abk +vt(k-1)*afk
          vqq  = tv0 +vq(k)*abk +vq(k-1)*afk

          t2sq = vtt*t2sq +vqq*c2sq
          r2sq = vtt*c2sq +vqq*r2sq
          c2sq = MAX( vtt*t2sq+vqq*r2sq, 0.0d0 )
          t3sq = vtt*t3sq +vqq*c3sq
          r3sq = vtt*c3sq +vqq*r3sq
          c3sq = MAX( vtt*t3sq+vqq*r3sq, 0.0d0 )
!
          cw25 = e1*( e2 + three*c1*e5c*gmel*qdiv**2 )/( three*eden )
!
!     **  Limitation on q, instead of L/q  **
          dlsq =  elsq
          IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)
!
!     **  Limitation on c3sq (0.12 =< cw =< 0.76) **
          ! Use Janjic's (2001; p 13-17) methodology (eqs 4.11-414 and 5.7-5.10)
          ! to calculate an exact limit for c3sq:
          auh = 27._kind_phys*a1*((a2*a2fac)**2)*b2*(gtr)**2
          aum = 54._kind_phys*(a1**2)*(a2*a2fac)*b2*c1*(gtr)
          adh = nine*a1*((a2*a2fac)**2)*(12._kind_phys*a1 + three*b2)*(gtr)**2
          adm = 18._kind_phys*(a1**2)*(a2*a2fac)*(b2 - three*(a2*a2fac))*(gtr)

          aeh = (nine*a1*((a2*a2fac)**2)*b1 +nine*a1*((a2*a2fac)**2)*            &
                (12._kind_phys*a1 + three*b2))*(gtr)
          aem = three*a1*(a2*a2fac)*b1*(three*(a2*a2fac) + three*b2*c1 +         &
                (18._kind_phys*a1*c1 - b2)) +                                    &
                (18._kind_phys)*(a1**2)*(a2*a2fac)*(b2 - three*(a2*a2fac))

          Req = -aeh/aem
          Rsl = (auh + aum*Req)/(three*adh + three*adm*Req)
          !For now, use default values, since tests showed little/no sensitivity
          Rsl = 0.12_kind_phys   !lower limit
          Rsl2= one - two*Rsl    !upper limit
          !IF (k==2)print*,"Dynamic limit RSL=",Rsl
          !IF (Rsl < 0.10 .OR. Rsl > 0.18) THEN
          !   print*,'--- ERROR: MYNN: Dynamic Cw '// &
          !        'limit exceeds reasonable limits'
          !   print*," MYNN: Dynamic Cw limit needs attention=",Rsl
          !ENDIF

          !JOE-Canuto/Kitamura mod
          !e2   = q3sq - e2c*ghel * qdiv**2
          !e3   = q3sq + e3c*ghel * qdiv**2
          !e4   = q3sq - e4c*ghel * qdiv**2
          e2   = q3sq - e2c*ghel*a2fac * qdiv**2
          e3   = q3sq + e3c*ghel*a2fac**2 * qdiv**2
          e4   = q3sq - e4c*ghel*a2fac * qdiv**2
          eden = e2*e4  + e3 *e5c*gmel * qdiv**2

          !JOE-Canuto/Kitamura mod
          !wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
          !     &        *( e2*e4c - e3c*e5c*gmel * qdiv**2 )
          wden = cc3*gtr**2 * dlsq**2/elsq * qdiv**2 &
               &        *( e2*e4c*a2fac - e3c*e5c*gmel*a2fac**2 * qdiv**2 )

          IF ( wden .NE. zero ) THEN
             !set dynamic limits
             clow = q3sq*( 0.12-cw25 )*eden/wden
             cupp = q3sq*( 0.76-cw25 )*eden/wden
             !clow = q3sq*( Rsl -cw25 )*eden/wden
             !cupp = q3sq*( Rsl2-cw25 )*eden/wden
!
             IF ( wden .GT. zero ) THEN
                c3sq  = MIN( MAX( c3sq, c2sq+clow ), c2sq+cupp )
             ELSE
                c3sq  = MAX( MIN( c3sq, c2sq+clow ), c2sq+cupp )
             END IF
          END IF
!
          e1   = e2 + e5c*gmel * qdiv**2
          eden = MAX( eden, 1.0d-20 )
!  Modified: Dec/22/2005, up to here

          !JOE-Canuto/Kitamura mod
          !e6c  = 3.0*a2*cc3*gtr * dlsq/elsq
          e6c  = 3.0*(a2*a2fac)*cc3*gtr * dlsq/elsq

          !============================
          !     **  for Gamma_theta  **
          !!          enum = qdiv*e6c*( t3sq-t2sq )
          IF ( t2sq .GE. zero ) THEN
             enum = MAX( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( t3sq-t2sq ), 0.0d0 )
          ENDIF
          gamt =-e1  *enum    /eden

          !============================
          !     **  for Gamma_q  **
          !!          enum = qdiv*e6c*( r3sq-r2sq )
          IF ( r2sq .GE. zero ) THEN
             enum = MAX( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ELSE
             enum = MIN( qdiv*e6c*( r3sq-r2sq ), 0.0d0 )
          ENDIF
          gamq =-e1  *enum    /eden

          !============================
          !     **  for Sm' and Sh'd(Theta_V)/dz  **
          !!          enum = qdiv*e6c*( c3sq-c2sq )
          enum = MAX( qdiv*e6c*( c3sq-c2sq ), 0.0d0)

          !JOE-Canuto/Kitamura mod
          !smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c+e4c)*a1/a2
          smd  = dlsq*enum*gtr/eden * qdiv**2 * (e3c*a2fac**2 + &
               & e4c*a2fac)*a1/(a2*a2fac)

          gamv = e1  *enum*gtr/eden
          sm(k) = sm(k) +smd

          !============================
          !     **  For elh (see below), qdiv at Level 3 is reset to 1.0.  **
          qdiv = one

          ! Level 3 debug prints
          IF ( debug_code ) THEN
            IF (sh(k)<-0.3 .OR. sm(k)<-0.3 .OR. &
              qke(k) < -0.1 .or. ABS(smd) .gt. 2.0) THEN
              print*," MYNN; mym_turbulence3.0; sh=",sh(k)," k=",k
              print*," gm=",gm(k)," gh=",gh(k)," sm=",sm(k)
              print*," q2sq=",q2sq," q3sq=",q3sq," q3/q2=",q3sq/q2sq
              print*," qke=",qke(k)," el=",el(k)," ri=",ri
              print*," PBLH=",pblh," u=",u(k)," v=",v(k)
            ENDIF
          ENDIF

!   **  Level 3 : end  **

       ELSE
!     **  At Level 2.5, qdiv is not reset.  **
          gamt = zero
          gamq = zero
          gamv = zero
       END IF
!
!      Add min background stability function (diffusivity) within model levels
!      with active plumes and clouds.
       cldavg = p5*(cldfra_bl1(k-1) + cldfra_bl1(k))
       mfmax  = max(p5*(edmf_a1(k-1)+edmf_a1(k))*p5*(edmf_w1(k-1)+edmf_w1(k)), abs(edmf_a_dd1(k)*edmf_w_dd1(k)))
       ! impose minimum for mass-flux columns
       sm(k) = max(sm(k), 0.04_kind_phys*min(10._kind_phys*mfmax,one) )
       sh(k) = max(sh(k), 0.04_kind_phys*min(10._kind_phys*mfmax,one) )
       ! impose minimum for clouds
       sm(k) = max(sm(k), 0.04_kind_phys*min(cldavg, p5) )
       sh(k) = max(sh(k), 0.04_kind_phys*min(cldavg, p5) )
       ! impose minimum for tkeprod_down
       if (tkeprod_dn(k)>1e-6) sm(k) = max(sm(k), 0.01_kind_phys)
       if (tkeprod_dn(k)>1e-6) sh(k) = max(sh(k), 0.01_kind_phys)
       ! impose minimum sm for tte configurations. this may overide Pr limitations above.
       if (closure .ge. 2.7 .and. closure .lt. 3.0) then
          sm(k) = max(sm(k),min(p2, max(zero,three*qpe(k))))
       endif
       !
       elq   = el(k)*qkw(k)
       elqt  = el(k)*qtw(k)  !for km with total turbulent energy; elq = elqt when closure < 2.7
       elh   = elq*qdiv

       ! Production of TKE (pdk), T-variance (pdt),
       ! q-variance (pdq), and covariance (pdc)
       pdk(k) = elq *(sm(k)*gm(k)                    &
            & +       sh(k)*gh(k)+gamv ) +           &
            &   p5*TKEprod_dn(k)         +           & ! xmchen
            &   p5*TKEprod_up(k)
       pdt(k) = elh *( sh(k)*dtl(k)+gamt )*dtl(k)
       pdq(k) = elh *( sh(k)*dqw(k)+gamq )*dqw(k)
       pdc(k) = elh *( sh(k)*dtl(k)+gamt )*dqw(k)*p5 &
            & + elh *( sh(k)*dqw(k)+gamq )*dtl(k)*p5

       ! Contergradient terms
       tcd(k) = elq*gamt
       qcd(k) = elq*gamq

       ! Eddy Diffusivity/Viscosity divided by dz
       dfm(k) = elqt*sm(k) / dzk
       dfh(k) = elq *sh(k) / dzk
!  Modified: Dec/22/2005, from here
!   **  In sub.mym_predict, dfq for the TKE and scalar variance **
!   **  are set to 3.0*dfm and 1.0*dfm, respectively. (Sqfac)   **
       dfq(k) =     dfm(k)
!  Modified: Dec/22/2005, up to here

   IF (tke_budget .eq. 1) THEN
       !TKE BUDGET
!       dudz = ( u(k)-u(k-1) )/dzk
!       dvdz = ( v(k)-v(k-1) )/dzk
!       dTdz = ( thl(k)-thl(k-1) )/dzk

!       upwp = -elq*sm(k)*dudz
!       vpwp = -elq*sm(k)*dvdz
!       Tpwp = -elq*sh(k)*dTdz
!       Tpwp = SIGN(MAX(ABS(Tpwp),1.E-6),Tpwp)

       
!!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB   

       !!!Shear Term
       !!!qSHEAR1(k)=-(upwp*dudz + vpwp*dvdz)
       qSHEAR1(k) = elq*sm(k)*gm(k) !staggered

       !!!Buoyancy Term    
       !!!qBUOY1(k)=grav*Tpwp/thl(k)
       !qBUOY1(k)= elq*(sh(k)*gh(k) + gamv)
       !qBUOY1(k) = elq*(sh(k)*(-dTdz*grav/thl(k)) + gamv) !! ORIGINAL CODE
       
       !! Buoyncy term takes the TKEprodTD(k) production now
       qBUOY1(k) = elq*(sh(k)*gh(k)+gamv)  +         &
       &           p5*TKEprod_dn(k)       +          & ! xmchen
       &           p5*TKEprod_up(k) 

       !!!Dissipation Term (now it evaluated in mym_predict)
       !qDISS1(k) = (q3sq**(3./2.))/(b1*MAX(el(k),1.)) !! ORIGINAL CODE
       
       !! >> EOB
    ENDIF

    END DO
!
    !make sure all inent(out) vars are initialized
    dfm(kts) = zero
    dfh(kts) = zero
    dfq(kts) = zero
    pdk(kts) = zero
    pdt(kts) = zero
    pdq(kts) = zero
    pdc(kts) = zero
    tcd(kts) = zero
    qcd(kts) = zero

    tcd(kte) = zero
    qcd(kte) = zero

!
    DO k = kts,kte-1
       dzk = dz(k)
       tcd(k) = ( tcd(k+1)-tcd(k) )/( dzk )
       qcd(k) = ( qcd(k+1)-qcd(k) )/( dzk )
    END DO
!
    if (spp_pbl==1) then
       DO k = kts,kte
          dfm(k)= dfm(k) + dfm(k)* pattern_spp_pbl1(k) * 1.5 * MAX(exp(-MAX(zw(k)-8000.,zero)/2000.),0.001)
          dfh(k)= dfh(k) + dfh(k)* pattern_spp_pbl1(k) * 1.5 * MAX(exp(-MAX(zw(k)-8000.,zero)/2000.),0.001)
       END DO
    endif

!    RETURN
#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_turbulence

! ==================================================================
!     SUBROUTINE  mym_predict:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       qke(nx,nz,ny) : qke at (n)th time level
!       tsq, ...cov     : ditto
!
!     Output variables:
!       qke(nx,nz,ny) : qke at (n+1)th time level
!       tsq, ...cov     : ditto
!
!     Work arrays:
!       qkw(nx,nz,ny)   : q at the center of the grid boxes        (m/s)
!       bp (nx,nz,ny)   : = 1/2*F,     see below
!       rp (nx,nz,ny)   : = P-1/2*F*Q, see below
!
!     # The equation for a turbulent quantity Q can be expressed as
!          dQ/dt + Ah + Av = Dh + Dv + P - F*Q,                      (1)
!       where A is the advection, D the diffusion, P the production,
!       F*Q the dissipation and h and v denote horizontal and vertical,
!       respectively. If Q is q^2, F is 2q/B_1L.
!       Using the Crank-Nicholson scheme for Av, Dv and F*Q, a finite
!       difference equation is written as
!          Q{n+1} - Q{n} = dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   )
!                        + dt/2*( Dv{n+1} - Av{n+1} - F*Q{n+1} ),    (2)
!       where n denotes the time level.
!       When the advection and diffusion terms are discretized as
!          dt/2*( Dv - Av ) = a(k)Q(k+1) - b(k)Q(k) + c(k)Q(k-1),    (3)
!       Eq.(2) can be rewritten as
!          - a(k)Q(k+1) + [ 1 + b(k) + dt/2*F ]Q(k) - c(k)Q(k-1)
!                 = Q{n} + dt  *( Dh{n}   - Ah{n}   + P{n} )
!                        + dt/2*( Dv{n}   - Av{n}   - F*Q{n}   ),    (4)
!       where Q on the left-hand side is at (n+1)th time level.
!
!       In this subroutine, a(k), b(k) and c(k) are obtained from
!       subprogram coefvu and are passed to subprogram tinteg via
!       common. 1/2*F and P-1/2*F*Q are stored in bp and rp,
!       respectively. Subprogram tinteg solves Eq.(4).
!
!       Modify this subroutine according to your numerical integration
!       scheme (program).
!
!-------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf
!! This subroutine predicts the turbulent quantities at the next step.
SUBROUTINE  mym_predict (kts,kte,                                     &
     &            closure,                                            &
     &            delt, dz,                                           &
     &            ust, flt, flq, pmz, phh,                            &
     &            el,  dfq, rho,                                      &
     &            pdk, pdt, pdq, pdc,                                 &
     &            qke, tsq, qsq, cov,                                 &
     &            s_aw1, s_awqke1, s_awqsq1,                          &
     &            bl_mynn_edmf, bl_mynn_edmf_tke,                     &
     &            qWT1, qDISS1,tke_budget, upwind                     )

!-------------------------------------------------------------------
  use module_bl_mynnedmf_common, only: qkemin,karman,sqfac,b1,b2,     &
       zero,one,two,three,p5,kind_phys
  
  integer, intent(in) :: kts,kte    

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

real(kind_phys), intent(in)    :: closure,upwind
integer, intent(in) :: bl_mynn_edmf,bl_mynn_edmf_tke,tke_budget
real(kind_phys), dimension(kts:kte), intent(in) :: dz, dfq, el, rho
real(kind_phys), dimension(kts:kte), intent(inout) :: pdk, pdt, pdq, pdc
real(kind_phys), intent(in)    :: flt, flq, pmz, phh
real(kind_phys), intent(in)    :: ust, delt
real(kind_phys), dimension(kts:kte), intent(inout) :: qke,tsq, qsq, cov
real(kind_phys), dimension(kts:kte+1), intent(inout) :: s_awqke1,s_awqsq1,s_aw1
    
!!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB 
real(kind_phys), dimension(kts:kte), intent(out) :: qWT1, qDISS1
real(kind_phys), dimension(kts:kte) :: tke_up,dzinv  
!! >> EOB
    
integer :: k
real(kind_phys), dimension(kts:kte) :: qkw, bp, rp, df3q
real(kind_phys):: vkz,pdk1,phm,pdt1,pdq1,pdc1,b1l,b2l,onoff
real(kind_phys), dimension(kts:kte) :: dtz,upcont,dncont
real(kind_phys), dimension(kts:kte) :: a,b,c,d,x

real(kind_phys), dimension(kts:kte) :: rhoinv
real(kind_phys), dimension(kts:kte+1) :: rhoz,kqdz,kmdz
!------------------------------------------------------------------
! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
IF (bl_mynn_edmf_tke == 0) THEN
   onoff=zero
ELSE
   onoff=one
ENDIF

!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
vkz = karman*p5*dz(kts)
!
!   **  dfq for the TKE is 3.0*dfm.  **
!
DO k = kts,kte
!   qke(k) = MAX(qke(k), zero)
   qkw(k) = SQRT( MAX( qke(k), qkemin ) )
   df3q(k)=Sqfac*dfq(k)
   dtz(k)=delt/dz(k)
END DO

!add stability criteria
!Prepare "constants" for diffusion equation.
!khdz = rho*Kh/dz = rho*dfh
rhoz(kts)  =rho(kts)
rhoinv(kts)=one/rho(kts)
kqdz(kts)  =rhoz(kts)*df3q(kts)
kmdz(kts)  =rhoz(kts)*dfq(kts)
DO k=kts+1,kte
   rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
   rhoz(k)  =    MAX(rhoz(k),1E-4_kind_phys)
   rhoinv(k)=one/MAX(rho(k) ,1E-4_kind_phys)
   kqdz(k)  = rhoz(k)*df3q(k) ! for TKE
   kmdz(k)  = rhoz(k)*dfq(k)  ! for T'2, q'2, and T'q'
ENDDO
rhoz(kte+1)=rhoz(kte)
kqdz(kte+1)=rhoz(kte+1)*df3q(kte)
kmdz(kte+1)=rhoz(kte+1)*dfq(kte)

if (bl_mynn_edmf == 1) then
   !stability criteria for mf
   DO k=kts+1,kte-1
      kqdz(k) = MAX(kqdz(k),  p5* s_aw1(k))
      kqdz(k) = MAX(kqdz(k), -p5*(s_aw1(k)-s_aw1(k+1)))
      kmdz(k) = MAX(kmdz(k),  p5* s_aw1(k))
      kmdz(k) = MAX(kmdz(k), -p5*(s_aw1(k)-s_aw1(k+1)))
!      kqdz(k) = max(kqdz(k),  0.5*(s_aw1(k)+sd_aw1(k)))
!      kqdz(k) = max(kqdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)) -0.5*(sd_aw1(k)-sd_aw1(k+1)) )
!      kmdz(k) = max(kmdz(k),  0.5*(s_aw1(k)+sd_aw1(k)))
!      kmdz(k) = max(kmdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)) -0.5*(sd_aw1(k)-sd_aw1(k+1)) )
   ENDDO
endif

!-----------------------------------------------------------------------
!   **  Compute production terms at the surface  **
!-----------------------------------------------------------------------
pdk1 = two*ust**3*pmz/( vkz )
phm  = two/ust   *phh/( vkz )
pdt1 = phm*flt**2
pdq1 = phm*flq**2
pdc1 = phm*flt*flq
!
! **  pdk(1)+pdk(2) corresponds to pdk1.  **
pdk(kts) = pdk1 - pdk(kts+1)

!!  pdt(kts) = pdt1 -pdt(kts+1)
!!  pdq(kts) = pdq1 -pdq(kts+1)
!!  pdc(kts) = pdc1 -pdc(kts+1)
pdt(kts) = pdt(kts+1)
pdq(kts) = pdq(kts+1)
pdc(kts) = pdc(kts+1)

!----------------------------------------------------------------------
!   **  Prediction of twice the turbulent kinetic energy  **
!----------------------------------------------------------------------
!! DO k = kts+1,kte-1
DO k = kts,kte-1
   b1l = b1*p5*( el(k+1)+el(k) )
   bp(k) = two*qkw(k) / b1l
   rp(k) = pdk(k+1) + pdk(k)
ENDDO
bp(kte) = zero

! Since df3q(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*df3q(k+1)+bp(k)*delt.
if (bl_mynn_edmf > 1) then

    DO k=kts+1,kte-1
       upcont(k)= s_awqke1(k)- s_aw1(k)*(qke(k)*upwind+qke(k-1)*(one-upwind))
       dncont(k)=zero !sd_awqke1(k)-sd_aw1(k)*(qke(k)*upwind+qke(k-1)*(one-upwind))                                                                                           
    ENDDO
    ! no flux at the top of the atmosphere
    upcont(kte)=zero
    dncont(kte)=zero

    k=kts
    a(1)=zero
    b(1)=one + dtz(k)*kqdz(k+1)*rhoinv(k) + bp(k)*delt
    c(1)=    - dtz(k)*kqdz(k+1)*rhoinv(k)
    d(1)=max(qkemin, qke(k)) + rp(k)*delt                        &
        &    - dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*kqdz(k)*rhoinv(k)
       b(k)=one+ dtz(k)*(kqdz(k)+kqdz(k+1))*rhoinv(k) + bp(k)*delt
       c(k)=   - dtz(k)*kqdz(k+1)*rhoinv(k)
       d(k)=max(qkemin, qke(k)) + rp(k)*delt                     &
          &    - dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
    ENDDO

else !implicit
   DO k=kts,kte-1
!       a(k-kts+1)=-dtz(k)*df3q(k)
!       b(k-kts+1)=1.+dtz(k)*(df3q(k)+df3q(k+1))+bp(k)*delt
!       c(k-kts+1)=-dtz(k)*df3q(k+1)
!       d(k-kts+1)=rp(k)*delt + qke(k)
!JOE 8/22/20 improve conservation
      a(k)=   - dtz(k)*kqdz(k)*rhoinv(k)                         &
          &   + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)*onoff
      b(k)=1. + dtz(k)*(kqdz(k)+kqdz(k+1))*rhoinv(k)             &
          &   + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*onoff &
          &   + bp(k)*delt
      c(k)=   - dtz(k)*kqdz(k+1)*rhoinv(k)                       &
          &   - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff
      d(k)=rp(k)*delt + qke(k)                                   &
          &   + dtz(k)*rhoinv(k)*(s_awqke1(k)-s_awqke1(k+1))*onoff
   ENDDO
endif

!!"no flux at top"
!  a(kte)=-1. !0.
!  b(kte)=1.
!  c(kte)=0.
!  d(kte)=0.
!! "prescribed value"
a(kte)=0.
b(kte)=1.
c(kte)=0.
d(kte)=qke(kte)

!CALL tridiag(kte,a,b,c,d)
CALL tridiag2(kte,a,b,c,d,x)

DO k=kts,kte
!   qke(k)=max(d(k-kts+1), qkemin)
   qke(k)=max(x(k), qkemin)
   qke(k)=min(qke(k), 150.)
ENDDO
      
   
!! TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB 
IF (tke_budget .eq. 1) THEN
   !! TKE Vertical transport << EOBvt
   tke_up=p5*qke
   dzinv=one/dz
   k=kts
   qWT1(k)=dzinv(k)*(                                                &
        &  (kqdz(k+1)*(tke_up(k+1)-tke_up(k))-kqdz(k)*tke_up(k))     &
        &  + p5*rhoinv(k)*(s_aw1(k+1)*tke_up(k+1)                    &
        &  +      (s_aw1(k+1)-s_aw1(k))*tke_up(k)                    &
        &  +      (s_awqke1(k)-s_awqke1(k+1)))*onoff) !unstaggered
   DO k=kts+1,kte-1
      qWT1(k)=dzinv(k)*(                                             &
            & (kqdz(k+1)*(tke_up(k+1)-tke_up(k))-kqdz(k)*(tke_up(k)-tke_up(k-1))) &
            &  + p5*rhoinv(k)*(s_aw1(k+1)*tke_up(k+1)                &
            &  +      (s_aw1(k+1)-s_aw1(k))*tke_up(k)                &
            &  -                  s_aw1(k)*tke_up(k-1)               &
            &  +      (s_awqke1(k)-s_awqke1(k+1)))*onoff) !unstaggered
   ENDDO
   k=kte
   qWT1(k)=dzinv(k)*(-kqdz(k)*(tke_up(k)-tke_up(k-1))                &
       &  + p5*rhoinv(k)*(-s_aw1(k)*tke_up(k)-s_aw1(k)*tke_up(k-1)+s_awqke1(k))*onoff) !unstaggered
   !!  >> EOBvt
   qDISS1=bp*tke_up !! TKE dissipation rate !unstaggered
END IF
!! >> EOB 

if ( closure > 2.5 ) then
   !-----------------------------------------------------------------
   !   **  Prediction of the moisture variance  **
   !-----------------------------------------------------------------
   do k = kts,kte-1
      b2l   = b2*p5*( el(k+1)+el(k) )
      bp(k) = two*qkw(k) / b2l
      rp(k) = pdq(k+1) + pdq(k)
   enddo
   bp(kte) = zero


   if (bl_mynn_edmf > 1) then
      do k=kts+1,kte-1
         upcont(k)= s_awqsq1(k)- s_aw1(k)*(qsq(k)*upwind+qsq(k-1)*(one-upwind))
         dncont(k)=zero !sd_awqsq1(k)-sd_aw1(k)*(qsq(k)*upwind+qsq(k-1)*(one-upwind))
      enddo
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*kmdz(k+1)*rhoinv(k) + bp(k)*delt
      c(1)=    - dtz(k)*kmdz(k+1)*rhoinv(k)
      d(1)=max(zero, qsq(k)) + rp(k)*delt                        &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      do k=kts+1,kte-1
         a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
         c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
         d(k)=max(zero, qsq(k)) + rp(k)*delt                     &
            &    - dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      enddo
       
   else !implicit

      do k=kts,kte-1
         a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)                         &
             &   + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)*onoff
         b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k)             &
             &   + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*onoff &
             &   + bp(k)*delt
         c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)                       &
             &   - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff
         d(k)=rp(k)*delt + qsq(k)                                   &
             &   + dtz(k)*rhoinv(k)*(s_awqsq1(k)-s_awqsq1(k+1))*onoff
      enddo
   endif

   a(kte)=-one
   b(kte)=one
   c(kte)=zero
   d(kte)=zero !qsq(kte)

   !call tridiag(kte,a,b,c,d)
   call tridiag2(kte,a,b,c,d,x)

   do k=kts,kte
      !qke(k)=max(d(k-kts+1), qkemin)
      qsq(k)=max(x(k), 1e-17_kind_phys)
      qsq(k)=min(qsq(k),6e-6_kind_phys)
   enddo

else
   !level 2.5 - use level 2 diagnostic
   do k = kts,kte-1
      if ( qkw(k) .le. zero ) then
         b2l = zero
      else
         b2l = b2*0.25_kind_phys*( el(k+1)+el(k) )/qkw(k)
      end if
      qsq(k) = min(6e-6_kind_phys, max(1e-17_kind_phys, b2l*( pdq(k+1)+pdq(k) )))
   enddo
   qsq(kte)=qsq(kte-1)
endif
    
!!!!!!!!!!!!!!!!!!!!!!end level 2.6   

    IF ( closure .GE. 2.7 ) THEN
!
!   **  dfq for the scalar variance is 1.0*dfm.  **
!
!   **  Prediction of the temperature variance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l   = b2*p5*( el(k+1)+el(k) )
          bp(k) = two*qkw(k) / b2l
          rp(k) = pdt(k+1) + pdt(k) 
       END DO
       
!zero gradient for tsq at bottom and top
       
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          !a(k-kts+1)=-dtz(k)*dfq(k)
          !b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          !c(k-kts+1)=-dtz(k)*dfq(k+1)
          !d(k-kts+1)=rp(k)*delt + tsq(k)
!JOE 8/22/20 improve conservation
          a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
          b(k)=one+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + max(zero, tsq(k))
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + tsq(k) - tsq(k)*bp(k)*delt
!!       ENDDO

       a(kte)=-1. !0.
       b(kte)=one
       c(kte)=zero
       d(kte)=zero
       
!       CALL tridiag(kte,a,b,c,d)
       CALL tridiag2(kte,a,b,c,d,x)

       DO k=kts,kte
          !tsq(k)=d(k-kts+1)
          tsq(k)=min(three, max(zero,x(k)))
       ENDDO
    ELSE
       !less than level 2.7 - default to level 2 diagnostic
       DO k = kts,kte-1
          IF ( qkw(k) .LE. zero ) THEN
             b2l = zero
          ELSE
             b2l = b2*0.25*( el(k+1)+el(k) )/qkw(k)
          END IF
          tsq(k) = min(three, max(zero, b2l*( pdt(k+1)+pdt(k) )))
       END DO
       tsq(kte)=tsq(kte-1)
    ENDIF

    !Full level 3 closure
    IF ( closure .GE. 3.0 ) THEN
!   **  Prediction of the temperature-moisture covariance  **
       DO k = kts,kte-1
          b2l = b2*p5*( el(k+1)+el(k) )
          bp(k) = two*qkw(k) / b2l
          rp(k) = pdc(k+1) + pdc(k) 
       END DO
       
!zero gradient for tqcov at bottom and top
!!       a(1)=0.
!!       b(1)=1.
!!       c(1)=-1.
!!       d(1)=0.

! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          !a(k-kts+1)=-dtz(k)*dfq(k)
          !b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))+bp(k)*delt
          !c(k-kts+1)=-dtz(k)*dfq(k+1)
          !d(k-kts+1)=rp(k)*delt + cov(k)
!JOE 8/22/20 improve conservation
          a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
          b(k)=one+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + cov(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + cov(k) - cov(k)*bp(k)*delt
!!       ENDDO

       a(kte)=-1. !0.
       b(kte)=one
       c(kte)=zero
       d(kte)=zero

!       CALL tridiag(kte,a,b,c,d)
       CALL tridiag2(kte,a,b,c,d,x)
       
       DO k=kts,kte
!          cov(k)=d(k-kts+1)
          cov(k)=x(k)
       ENDDO
       
    ELSE

       !Not level 3 - default to level 2 diagnostic
       DO k = kts,kte-1
          IF ( qkw(k) .LE. zero ) THEN
             b2l = zero
          ELSE
             b2l = b2*0.25*( el(k+1)+el(k) )/qkw(k)
          END IF
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO
       cov(kte)=cov(kte-1)
      
    END IF

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_predict
  
! ==================================================================
!     SUBROUTINE  mym_condensation:
!
!     Input variables:    see subroutine mym_initialize and turbulence
!       exner(nz)    : Perturbation of the Exner function    (J/kg K)
!                         defined on the walls of the grid boxes
!                         This is usually computed by integrating
!                         d(pi)/dz = h*g*tv/tref**2
!                         from the upper boundary, where tv is the
!                         virtual potential temperature minus tref.
!
!     Output variables:   see subroutine mym_initialize
!       cld(nx,nz,ny)   : Cloud fraction
!
!     Work arrays/variables:
!       qmq             : Q_w-Q_{sl}, where Q_{sl} is the saturation
!                         specific humidity at T=Tl
!       alp(nx,nz,ny)   : Functions in the condensation process
!       bet(nx,nz,ny)   : ditto
!       sgm(nx,nz,ny)   : Combined standard deviation sigma_s
!                         multiplied by 2/alp
!
!     # qmq, alp, bet and sgm are allowed to share storage units with
!       any four of other work arrays for saving memory.
!
!     # Results are sensitive particularly to values of cp and r_d.
!       Set these values to those adopted by you.
!
!-------------------------------------------------------------------
!>\ingroup gsd_mynn_edmf 
!! This subroutine calculates the nonconvective component of the 
!! subgrid cloud fraction and mixing ratio as well as the functions used to 
!! calculate the buoyancy flux. Different cloud PDFs can be selected by
!! use of the namelist parameter \p bl_mynn_cloudpdf .
  SUBROUTINE  mym_condensation (kts,kte,   &
    &            dx, dz, zw, xland,        &
    &            thl, qw, qv, qc, qi, qs,  &
    &            p,exner,                  &
    &            tsq, qsq, cov,            &
    &            Sh, el, bl_mynn_cloudpdf, &
    &            qc_bl1, qi_bl1,           &
    &            cldfra_bl1,               &
    &            PBLH1,HFX1,               &
    &            Vt, Vq, th, sgm,          &
    &            closure,                  &
    &            spp_pbl, pattern_spp_pbl1 )

!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: cp,xlv,r_d,r_v,xlvcp,ep_2,   &
         ep_3,p608,cpv,rr2,rrp,tv0,tice,tliq,b2,                      &
         zero,one,two,three,four,five,ten,hundred,p1,p2,p5,           &
         kind_phys
    
    integer, intent(in)   :: kts,kte, bl_mynn_cloudpdf

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), intent(in)      :: HFX1,xland
    real(kind_phys), intent(in)      :: dx,pblh1,closure
    real(kind_phys), dimension(kts:kte), intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), dimension(kts:kte), intent(in) :: p,exner,thl,qw,   &
         &qv,qc,qi,qs,tsq,qsq,cov,th

    real(kind_phys), dimension(kts:kte), intent(inout) :: vt,vq,sgm

    real(kind_phys), dimension(kts:kte) :: alp,a,bet,b,ql,q1,RH
    real(kind_phys), dimension(kts:kte), intent(out) :: qc_bl1,qi_bl1,   &
         &cldfra_bl1
    DOUBLE PRECISION :: t3sq, r3sq, c3sq

    real(kind_phys):: qsl,esat,qsat,dqsl,cld0,q1k,qlk,eq1,qll,           &
         &q2p,pt,rac,qt,t,xl,rsl,cpm,Fng,qww,alpha,beta,bb,sgmq,sgmc,    &
         &ls,wt,wt2,cld_factor,fac_damp,liq_frac,ql_ice,ql_water,        &
         &qmq,qsat_tk,q1_rh,rh_adj,zsl,maxqc,cldfra_rh,cldfra_qsq,       &
         &cldfra_rh0,cldfra_rh1,cldfra_qsq0,cldfra_qsq1,clim,qlim
    !lower limits for sgm (for mixing ratio estimates) in case sgm falls out (% of qw)
    real(kind_phys), parameter :: qlim_sfc =0.008
    real(kind_phys), parameter :: qlim_pbl =0.020
    real(kind_phys), parameter :: qlim_trp =0.025
    !lower limits for sqm (for cloud fraction) in case sgm falls out (% of qw)
    real(kind_phys), parameter :: clim_sfc =0.010
    real(kind_phys), parameter :: clim_pbl =0.020
    real(kind_phys), parameter :: clim_trp =0.030
    real(kind_phys), parameter :: rhcrit   =0.83 !for cloudpdf = 2
    real(kind_phys), parameter :: rhmax    =1.10 !for cloudpdf = 2
    integer :: i,j,k

    real(kind_phys):: erf

    !VARIABLES FOR ALTERNATIVE SIGMA
    real:: dth,dtl,dqw,dzk,els
    real(kind_phys), dimension(kts:kte), intent(in) :: Sh,el

    !variables for SGS BL clouds
    real(kind_phys)           :: zagl,damp,PBLH2
    real(kind_phys)           :: cfmax

    !JAYMES:  variables for tropopause-height estimation
    real(kind_phys)           :: theta1, theta2, ht1, ht2
    integer                   :: k_tropo

!   Stochastic
    integer,  intent(in)      :: spp_pbl
    real(kind_phys), dimension(kts:kte) :: pattern_spp_pbl1
    real(kind_phys)           :: qw_pert

! First, obtain an estimate for the tropopause height (k), using the method employed in the
! Thompson subgrid-cloud scheme.  This height will be a consideration later when determining 
! the "final" subgrid-cloud properties.
! JAYMES:  added 3 Nov 2016, adapted from G. Thompson

    DO k = kte-3, kts, -1
       theta1 = th(k)
       theta2 = th(k+2)
       ht1 = 44307.692 * (one - (p(k)/101325.)**0.190)
       ht2 = 44307.692 * (one - (p(k+2)/101325.)**0.190)
       if ( (((theta2-theta1)/(ht2-ht1)) .lt. 10./1500. ) .AND.       &
     &                       (ht1.lt.19000.) .and. (ht1.gt.4000.) ) then 
          goto 86
       endif
    ENDDO
 86   continue
    k_tropo = MAX(kts+2, k+2)

    zagl = 0.

    SELECT CASE(bl_mynn_cloudpdf)

      CASE (0) ! ORIGINAL MYNN PARTIAL-CONDENSATION SCHEME

        DO k = kts,kte-1
           t  = th(k)*exner(k)

!x      if ( ct .gt. zero ) then
!       a  =  17.27
!       b  = 237.3
!x      else
!x        a  =  21.87
!x        b  = 265.5
!x      end if
!
!   **  3.8 = 0.622*6.11 (hPa)  **

           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           !qsl=ep_2*esat/(p(k)-ep_3*esat)
           qsl=ep_2*esat/max(1.e-4,(p(k)-ep_3*esat))
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*xlv/( r_d*t**2 )

           alp(k) = one/( one+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           !Sommeria and Deardorff (1977) scheme, as implemented
           !in Nakanishi and Niino (2009), Appendix B
           t3sq = MAX( tsq(k), zero )
           r3sq = MAX( qsq(k), zero )
           c3sq =      cov(k)
           c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )
           r3sq = r3sq +bet(k)**2*t3sq -2.0*bet(k)*c3sq
           !DEFICIT/EXCESS WATER CONTENT
           qmq  = qw(k) -qsl
           !ORIGINAL STANDARD DEVIATION
           sgm(k) = SQRT( MAX( r3sq, 1.0d-10 ))
           !NORMALIZED DEPARTURE FROM SATURATION
           q1(k)   = qmq / sgm(k)
           !CLOUD FRACTION. rr2 = 1/SQRT(2) = 0.707
           cldfra_bl1(k) = 0.5*( one+erf( q1(k)*rr2 ) )

           q1k  = q1(k)
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cldfra_bl1(k)*q1k + eq1, zero )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql(k) = alp(k)*sgm(k)*qll
           !LIMIT SPECIES TO TEMPERATURE RANGES
           liq_frac = min(one, max(zero,(t-240.0)/29.0))
           qc_bl1(k) = liq_frac*ql(k)
           qi_bl1(k) = (one - liq_frac)*ql(k)

           !Now estimate the buoyancy flux functions
           q2p = xlvcp/exner(k)
           pt = thl(k) +q2p*ql(k) ! potential temp

           !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
           qt   = one +p608*qw(k) -(1.+p608)*(qc_bl1(k)+qi_bl1(k))*cldfra_bl1(k)
           rac  = alp(k)*( cldfra_bl1(K)-qll*eq1 )*( q2p*qt-(1.+p608)*pt )

           !BUOYANCY FACTORS: wherever vt and vq are used, there is a
           !"+1" and "+tv0", respectively, so these are subtracted out here.
           !vt is unitless and vq has units of K.
           vt(k) =      qt-one -rac*bet(k)
           vq(k) = p608*pt-tv0 +rac

        END DO

      CASE (1, -1) !ALTERNATIVE FORM (Nakanishi & Niino 2004 BLM, eq. B6, and
                       !Kuwano-Yoshida et al. 2010 QJRMS, eq. 7):
        DO k = kts,kte-1
           t  = th(k)*exner(k)
           !SATURATED VAPOR PRESSURE
           esat = esat_blend(t)
           !SATURATED SPECIFIC HUMIDITY
           !qsl=ep_2*esat/(p(k)-ep_3*esat)
           qsl=ep_2*esat/max(1.e-4,(p(k)-ep_3*esat))
           !dqw/dT: Clausius-Clapeyron
           dqsl = qsl*ep_2*xlv/( r_d*t**2 )

           alp(k) = one/( one+dqsl*xlvcp )
           bet(k) = dqsl*exner(k)

           if (k .eq. kts) then 
             dzk = 0.5*dz(k)
           else
             dzk = dz(k)
           end if
           dth = 0.5*(thl(k+1)+thl(k)) - 0.5*(thl(k)+thl(MAX(k-1,kts)))
           dqw = 0.5*(qw(k+1) + qw(k)) - 0.5*(qw(k) + qw(MAX(k-1,kts)))
           sgm(k) = SQRT( MAX( (alp(k)**2 * MAX(el(k)**2, p1) * &
                             b2 * MAX(Sh(k),0.03))/four * &
                      (dqw/dzk - bet(k)*(dth/dzk ))**2 , 1.0e-10) )
           qmq   = qw(k) -qsl
           q1(k) = qmq / sgm(k)
           cldfra_bl1(K) = p5*( one+erf( q1(k)*rr2 ) )

           !now compute estimated lwc for PBL scheme's use 
           !qll IS THE NORMALIZED LIQUID WATER CONTENT (Sommeria and
           !Deardorff (1977, eq 29a). rrp = 1/(sqrt(2*pi)) = 0.3989
           q1k  = q1(k)
           eq1  = rrp*EXP( -0.5*q1k*q1k )
           qll  = MAX( cldfra_bl1(K)*q1k + eq1, zero )
           !ESTIMATED LIQUID WATER CONTENT (UNNORMALIZED)
           ql (k) = alp(k)*sgm(k)*qll
           liq_frac = min(one, max(zero,(t-240.0)/29.0))
           qc_bl1(k) = liq_frac*ql(k)
           qi_bl1(k) = (one - liq_frac)*ql(k)

           !Now estimate the buoyancy flux functions
           q2p = xlvcp/exner(k)
           pt = thl(k) +q2p*ql(k) ! potential temp

           !qt is a THETA-V CONVERSION FOR TOTAL WATER (i.e., THETA-V = qt*THETA)
           qt   = one +p608*qw(k) -(1.+p608)*(qc_bl1(k)+qi_bl1(k))*cldfra_bl1(k)
           rac  = alp(k)*( cldfra_bl1(K)-qll*eq1 )*( q2p*qt-(1.+p608)*pt )

           !BUOYANCY FACTORS: wherever vt and vq are used, there is a
           !"+1" and "+tv0", respectively, so these are subtracted out here.
           !vt is unitless and vq has units of K.
           vt(k) =      qt-one -rac*bet(k)
           vq(k) = p608*pt-tv0 +rac

        END DO

      CASE (2, -2)

        !Diagnostic statistical scheme of Chaboureau and Bechtold (2002), JAS
        !but with use of higher-order moments to estimate sigma
        pblh2=MAX(ten,pblh1)
        DO k = kts,kte-1
           zagl   = zw(k) + p5*dz(k)

           t      = th(k)*exner(k)
           xl     = xl_blend(t)              ! obtain latent heat
           qsat_tk= qsat_blend(t,  p(k))     ! saturation water vapor mixing ratio at tk and p
           rh(k)  = MAX(MIN(rhmax, qw(k)/MAX(1E-10_kind_phys,qsat_tk)),0.001_kind_phys)

           !dqw/dT: Clausius-Clapeyron
           dqsl   = qsat_tk*ep_2*xlv/( r_d*t**2 )
           alp(k) = one/(one + dqsl*xlvcp )
           bet(k) = dqsl*exner(k)
 
           rsl    = xl*qsat_tk / (r_v*t**2)  ! slope of C-C curve at t (=abs temperature)
                                             ! CB02, Eqn. 4
           cpm    = cp + qw(k)*cpv           ! CB02, sec. 2, para. 1
           a(k)   = one/(one + xl*rsl/cpm)   ! CB02 variable "a"
           b(k)   = a(k)*rsl                 ! CB02 variable "b"

           !SPP
           qw_pert= qw(k) + qw(k)*p5*pattern_spp_pbl1(k)*real(spp_pbl,kind=kind_phys)

           !This form of qmq (the numerator of Q1) no longer uses the a(k) factor
           qmq    = qw_pert - qsat_tk          ! saturation deficit/excess;

           !Use the form of Eq. (6) in Chaboureau and Bechtold (2002)
           !except neglect all but the first term for sig_r
           r3sq   = max( qsq(k), zero )
           t3sq   = max( tsq(k), zero )
           c3sq   =      cov(k)
           c3sq   = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )

           if (closure .lt. 2.7) then
              r3sq = r3sq
           elseif ((closure .ge. 2.7) .and. (closure .lt. 3.0)) then
              r3sq = r3sq +bet(k)**2*t3sq
           else
              r3sq = r3sq +bet(k)**2*t3sq -two*bet(k)*c3sq
           endif
           
           !Calculate sigma using higher-order moments:
           sgm(k) = max(1e-13, sqrt( real(r3sq,kind_phys) ))
           !Set constraints on sigma relative to total water
           sgm(k) = min( sgm(k), qw(k)*p5 )
           
           !introduce vertical grid spacing dependence on min sgm
           wt     = min(one, max(zero, dz(k)-100.)/500.) !=0 for dz < 100 m, =1 for dz > 600 m
           sgm(k) = sgm(k) + sgm(k)*p2*wt !inflate sgm for coarse dz
           !save sgm for mixing ratio and cloud fraction estimates
           sgmq   = sgm(k)
           sgmc   = sgm(k)
           
           !allow minimum sgm to vary with z.
           wt     = min(one, max(zero, (zagl - (pblh2+10.)))/250.) !0 in pbl, 1 aloft
           clim   = clim_pbl*(one-wt) + clim_trp*wt
           zsl    = min(150., max(50., p1*pblh2))        !crude ekman layer
           wt     = min(one, max(zero, zagl - zsl)/150.)  !0 near sfc, 1 above 
           clim   = clim_sfc*(one-wt) + clim*wt
           sgmc   = max( sgmc, qw(k)*clim )
           !apply absolute lower limit in case qw = 0.
           sgmc   = max(1e-13, sgmc)
           !For cloud fractions in saturated conditions, apply lower limit on sgmc
           if (qmq .ge. zero) sgmc = max(0.02*qw(k), sgmc)
           
           q1(k)  = qmq  / sgmc  ! Q1, the normalized saturation

           !Add condition for falling/settling into low-RH layers, so at least
           !some cloud fraction is applied for all qc, qs, and qi.
           wt2    = min(one, max(zero, zagl - pblh2)/300.) !0 in pbl, 1 aloft
           !ensure adequate RH & q1 when qi is at least 1e-9 (above the PBLH)
           if ((qi(k)+qs(k))>1.e-10 .and. (zagl .gt. pblh2)) then
              rh_adj  =min(rhmax, rhcrit + wt2*0.037_kind_phys*(max(zero, ten + log10(qi(k)+qs(k))) ))
              rh_adj  =max(rh(k), rh_adj)
              !add rh-based q1
              q1_rh   =-3. + three*(rh_adj-rhcrit)/(one-rhcrit)
              q1(k)   =max(q1_rh, q1(k) )
           endif
           !ensure adequate rh & q1 when qc is at least 1e-5 (above the PBLH)
           if (qc(k)>1.e-5 .and. (zagl .gt. pblh2)) then
              rh_adj  =min(rhmax, rhcrit + wt2*0.07_kind_phys*(max(zero, five + log10(qc(k))) ))
              rh_adj  =max(rh(k), rh_adj)
              !add rh-based q1
              q1_rh   =-3. + three*(rh_adj-rhcrit)/(one-rhcrit)
              q1(k)   =max(q1_rh, q1(k) )
           endif

           q1k   = q1(k)          ! backup Q1 for later modification

           ! Specify cloud fraction
           !Original C-B cloud fraction, allows cloud fractions out to q1 = -3.5
           !cldfra_bl1(K) = max(0., min(1., 0.5+0.36*atan(1.55*q1(k)))) ! Eq. 7 in CB02
           !Waynes LES fit  - over-diffuse, when limits removed from vt & vq & fng
           !cldfra_bl1(K) = max(0., min(1., 0.5+0.36*atan(1.2*(q1(k)+0.4))))
           !Best compromise: Improves marine stratus without adding much cold bias.
           !RRFSv1
           !cldfra_bl1(k) = max(zero, min(one, 0.5+0.36*atan(1.8*(q1k+0.2))))
           !exp1
           !cldfra_bl1(k) = max(zero, min(one, p5+0.36*atan(1.65*q1k)))

           !Use specialized forms within and outside the pbl.
           wt2           = min(one, max(zero, (zagl - (pblh1-100.))/200.)) !0 in pbl, 1 aloft

           !cldfra_qsq0   = max(zero, min(one, p5+0.35*atan(4.1*(q1k))))
           cldfra_qsq0   = max(zero, min(one, p5+0.35*atan(3.6*(q1k+0.05))))
           !cldfra_qsq1   = max(zero, min(one, p5+0.37*atan(2.1*(q1k+0.4))))
           cldfra_qsq1   = max(zero, min(one, p5+0.39*atan(1.6*(q1k+0.55))))
           cldfra_qsq    = cldfra_qsq0*(one-wt2) + cldfra_qsq1*wt2

           !For ceiling detection, apply minimum rh-based cloud fraction
           cldfra_rh0    = min(one, max(zero, 0.56*tanh((rh(k)-0.976)/0.030)+p5))
           cldfra_rh1    = min(one, max(zero, 0.55*tanh((rh(k)-0.955)/0.055)+p5))
           cldfra_rh     = zero !cldfra_rh0*(one-wt2) + cldfra_rh1*wt2

           cldfra_bl1(k) = max(cldfra_qsq, cldfra_rh)
           
           !Specify hydrometeors (grid mean = in-cloud * cloud fraction)
           !allow minimum sgmq (lower limit of mixing ratios) to vary with z.
           wt     = min(one, max(zero, (zagl - (pblh2+ten)))/300.) !0 in pbl, 1 aloft
           qlim   = qlim_pbl*(one-wt) + qlim_trp*wt
           zsl    = min(150., max(50., p1*pblh2)) !height of surface layer
           wt     = min(one, max(zero, zagl - zsl)/200.)  !0 near sfc, 1 above
           qlim   = qlim_sfc*(one-wt) + qlim*wt
           sgmq   = max(sgmq, qw(k)*qlim)
           
           ql_water = min(sgmq, 0.025*qw(k))*cldfra_bl1(k)
           ql_ice   = min(sgmq, 0.025*qw(k))*cldfra_bl1(k)

           ! The cloud water formulations are taken from CB02, Eq. 8.
!           maxqc = max(qw(k) - qsat_tk, zero)
!           if (q1k < zero) then        !unsaturated
!              !orig: ql_water = sgm(k)*exp(1.2*q1k-one)
!              !orig: ql_ice   = sgm(k)*exp(1.2*q1k-one)
!              ql_water = min(sgm(k),0.03*qw(k))*cldfra_bl1(k)
!              ql_ice   = min(sgm(k),0.03*qw(k))*cldfra_bl1(k)
!           elseif (q1k > 2.) then !supersaturated
!              ql_water = min(sgm(k)*q1k, maxqc)
!              ql_ice   =     sgm(k)*q1k
!           else                      !slightly saturated (0 > q1 < 2)
!              !orig: ql_water = min(sgm(k)*(exp(-1.) + 0.66*q1k + 0.086*q1k**2), maxqc)
!              !orig: ql_ice   =     sgm(k)*(exp(-1.) + 0.66*q1k + 0.086*q1k**2)
!              ql_water = min(sgm(k),0.02*qw(k))*(exp(-1.) + 0.66*q1k + 0.086*q1k**2)
!              ql_ice   = min(sgm(k),0.02*qw(k))*(exp(-1.) + 0.66*q1k + 0.086*q1k**2)
!           endif

           !In saturated grid cells, use average of SGS and resolved values
           !if ( qc(k) > 1.e-6 ) ql_water = 0.5 * ( ql_water + qc(k) ) 
           !ql_ice is actually the total frozen condensate (snow+ice),
           !if ( (qi(k)+qs(k)) > 1.e-9 ) ql_ice = 0.5 * ( ql_ice + (qi(k)+qs(k)) )

           if (cldfra_bl1(k) < 0.001) then
              ql_ice        = zero
              ql_water      = zero
              cldfra_bl1(k) = zero
           endif

           liq_frac = MIN(one, MAX(zero, (t-tice)/(tliq-tice)))
           qc_bl1(k) = liq_frac*ql_water       ! apply liq_frac to ql_water and ql_ice
           qi_bl1(k) = (one-liq_frac)*ql_ice

           !Above tropopause:  eliminate subgrid clouds from CB scheme. Note that this was
           !"k_tropo - 1" as of 20 Feb 2023. Changed to allow more high-level clouds.
           if (k .ge. k_tropo) then
              cldfra_bl1(K) = zero
              qc_bl1(k)     = zero
              qi_bl1(k)     = zero
           endif

           !Buoyancy-flux-related calculations follow...
           !limiting Q1 to avoid too much diffusion in cloud layers
           !when using the buoyancy-flux functions.
           !q1k=max(Q1(k),-2.0)
           if ((xland-1.5).GE.zero) then  ! water
              q1k=max(Q1(k),-2.5)
           else                           ! land
              q1k=max(Q1(k),-2.0)
           endif
           ! "Fng" represents the non-Gaussian transport factor
           ! (non-dimensional) from Bechtold and Siebesma (1998, JAS)
           if (q1k .ge. one) then
              Fng = one
           elseif (q1k .ge. -1.7 .and. q1k .lt. one) then
              Fng = exp(-0.4*(q1k-one))
           elseif (q1k .ge. -2.5 .and. q1k .lt. -1.7) then
              Fng = three + exp(-3.8*(q1k+1.7))
           else
              Fng = min(23.9 + exp(-1.6*(q1k+2.5)), 60._kind_phys)
           endif

           cfmax = min(cldfra_bl1(k), 0.6_kind_phys)
           !Further limit the cf going into vt & vq near the surface
           zsl   = min(max(25., p1*pblh2), hundred)
           wt    = min(zagl/zsl, one) !=0 at z=0 m, =1 above ekman layer
           cfmax = cfmax*wt

           bb = b(k)*t/th(k) ! bb is "b" in BCMT95.  Their "b" differs from
                             ! "b" in CB02 (i.e., b(k) above) by a factor
                             ! of T/theta.  Strictly, b(k) above is formulated in
                             ! terms of sat. mixing ratio, but bb in BCMT95 is
                             ! cast in terms of sat. specific humidity.  The
                             ! conversion is neglected here.
           qww   = one + p608*qw(k)
           alpha = p608*th(k)
           beta  = (th(k)/t)*(xl/cp) - 1.61*th(k)
           vt(k) = qww   - cfmax*beta*bb*Fng   - one
           vq(k) = alpha + cfmax*beta*a(k)*Fng - tv0
           ! vt and vq correspond to beta-theta and beta-q, respectively,  
           ! in NN09, Eq. B8.  They also correspond to the bracketed
           ! expressions in BCMT95, Eq. 15, since (s*ql/sigma^2) = cldfra*Fng
           ! The "-1" and "-tv0" terms are included for consistency with 
           ! the legacy vt and vq formulations (above).
        enddo

      END SELECT !end cloudPDF option

      !For testing purposes only, option for isolating on the mass-flux clouds.
      IF (bl_mynn_cloudpdf .LT. 0) THEN
         DO k = kts,kte-1
            cldfra_bl1(k) = zero
            qc_bl1(k)     = zero
            qi_bl1(k)     = zero
         END DO
      ENDIF

      ql(kte)        = ql(kte-1)
      vt(kte)        = vt(kte-1)
      vq(kte)        = vq(kte-1)
      qc_bl1(kte)    = zero
      qi_bl1(kte)    = zero
      cldfra_bl1(kte)= zero
    RETURN

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mym_condensation

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine solves for tendencies of U, V, \f$\theta\f$, qv,
!! qc, and qi
  SUBROUTINE mynn_tendencies(kts,kte,i,       &
       &delt,dz,zw,xland,pblh,rho,            &
       &u,v,th,tk,qv,qc,qi,qs,qnc,qni,        &
       &psfc,p,exner,                         &
       &thl,sqv,sqc,sqi,sqs,sqw,              &
       &thl_tot1,qc_tot1,qi_tot1,             &
       &qnwfa,qnifa,qnbca,ozone,              &
       &ust,flt,flq,flqv,flqc,wspd,           &
       &uoce,voce,                            &
       &tsq,qsq,cov,                          &
       &tcd,qcd,                              &
       &dfm,dfh,dfq,                          &
       &Du,Dv,Dth,Dqv,Dqc,Dqi,Dqs,Dqnc,Dqni,  &
       &Dqnwfa,Dqnifa,Dqnbca,Dozone,          &
       &diss_heat,                            &
       &s_aw1,s_awthl1,                       &
       &s_awqt1,s_awqv1,s_awqc1,              &
       &s_awu1,s_awv1,                        &
       &s_awqnc1,s_awqni1,                    &
       &s_awqnwfa1,s_awqnifa1,s_awqnbca1,     &
       &sd_aw1,sd_awthl1,sd_awqt1,sd_awqv1,   &
       &sd_awqc1,sd_awqi1,                    &
       &sd_awqnc1,sd_awqni1,                  &
       &sd_awqnwfa1,sd_awqnifa1,              &
       &sd_awu1,sd_awv1,                      &
       &sub_thl,sub_sqv,                      &
       &sub_u,sub_v,                          &
       &det_thl,det_sqv,det_sqc,              &
       &det_u,det_v,                          &
       &FLAG_QC,FLAG_QI,FLAG_QNC,FLAG_QNI,    &
       &FLAG_QS,                              &
       &FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,     &
       &FLAG_OZONE,                           &
       &cldfra_bl1,                           &
       &bl_mynn_cloudmix,                     &
       &bl_mynn_mixqt,                        &
       &bl_mynn_edmf,                         &
       &bl_mynn_edmf_mom,                     &
       &bl_mynn_mixscalars,                   &
       &bl_mynn_mixaerosols,                  &
       &bl_mynn_mixnumcon,                    &
       &upwind                                )

!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: r_d,p608,grav,xlvcp,xlscp,       &
         wfa_max,wfa_min,ifa_max,ifa_min,wfa_ht,ifa_ht,zero,one,          &
         p333,p5,kind_phys
    
    integer, intent(in) :: kts,kte,i

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in) :: bl_mynn_cloudmix,bl_mynn_mixqt,                &
                           bl_mynn_edmf,bl_mynn_edmf_mom,                 &
                           bl_mynn_mixscalars,bl_mynn_mixaerosols,        &
                           bl_mynn_mixnumcon
    logical, intent(in) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QS,              &
         &FLAG_QNC,FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,FLAG_OZONE

! thl - liquid water potential temperature
! qw - total water
! dfm,dfh,dfq - diffusivities i.e., dfh(k) = elq*sh(k) / dzk
! flt - surface flux of thl
! flq - surface flux of qw

! mass-flux plumes
    real(kind_phys), dimension(kts:kte+1), intent(in) :: s_aw1,            &
         &s_awthl1,s_awqt1,s_awqnc1,s_awqni1,s_awqv1,s_awqc1,s_awu1,s_awv1,&
         &s_awqnwfa1,s_awqnifa1,s_awqnbca1,                                &
         &sd_aw1,sd_awthl1,sd_awqt1,sd_awqv1,sd_awqc1,sd_awqi1,            &
         &sd_awqnc1,sd_awqni1,sd_awqnwfa1,sd_awqnifa1,sd_awu1,sd_awv1
! tendencies from mass-flux environmental subsidence and detrainment
    real(kind_phys), dimension(kts:kte), intent(in) :: sub_thl,sub_sqv,   &
         &sub_u,sub_v,det_thl,det_sqv,det_sqc,det_u,det_v
    real(kind_phys), dimension(kts:kte), intent(in) :: u,v,th,tk,qv,qc,qi,&
         &qs,qni,qnc,rho,p,exner,dfq,dz,zw,tsq,qsq,cov,tcd,qcd,           &
         &cldfra_bl1,diss_heat,thl_tot1,qc_tot1,qi_tot1
    real(kind_phys), dimension(kts:kte), intent(inout) :: thl,sqw,sqv,sqc,&
         &sqi,sqs,qnwfa,qnifa,qnbca,ozone,dfm,dfh
    real(kind_phys), dimension(kts:kte), intent(inout) :: du,dv,dth,dqv,  &
         &dqc,dqi,dqs,dqni,dqnc,dqnwfa,dqnifa,dqnbca,dozone
    real(kind_phys), intent(in) :: flt,flq,flqv,flqc,uoce,voce
    real(kind_phys), intent(in) :: ust,delt,psfc,wspd,xland,pblh,upwind
    !debugging
    real(kind_phys):: wsp,wsp2,tk2,th2
    logical :: problem
    integer :: kproblem

!    real(kind_phys), intent(in) :: gradu_top,gradv_top,gradth_top,gradqv_top

!local vars

    real(kind_phys), dimension(kts:kte) :: dtz,dfhc,dfmc,delp
    real(kind_phys), dimension(kts:kte) :: sqv2,sqc2,sqi2,sqs2,sqw2,      &
          &qni2,qnc2,qnwfa2,qnifa2,qnbca2,ozone2
    real(kind_phys), dimension(kts:kte) :: zfac,plumeKh,rhoinv
    real(kind_phys), dimension(kts:kte) :: upcont,dncont ! updraft/downdraft contribution to fluxes for explicit calculation
    real(kind_phys), dimension(kts:kte) :: a,b,c,d,x
    real(kind_phys), dimension(kts:kte+1) :: rhoz,                        & !rho on model interface
          &khdz,kmdz
    real(kind_phys):: rhs,gfluxm,gfluxp,dztop,maxdfh,mindfh,maxcf,maxKh
    real(kind_phys):: t,esat,qsl,onoff,kh,km,dzk,rhosfc,ustovwsp
    real(kind_phys):: qvflux,aero_min,aero_max
    real(kind_phys):: th_new,portion_qc,portion_qi,condensate,qsat
    integer :: k,kk

    !Activate nonlocal mixing from the mass-flux scheme for
    !number concentrations and aerosols (0.0 = no; 1.0 = yes)
    real(kind_phys), parameter :: nonloc  = 1.0
    real(kind_phys), parameter :: nc_min  = 100.0
    real(kind_phys), parameter :: ni_min  = 1e-6
    real(kind_phys) :: wfa_max2, wt

    dztop=p5*(dz(kte)+dz(kte-1))

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    ! Note that s_awu and s_awv already come in as 0.0 if bl_mynn_edmf_mom == 0, so
    ! we only need to zero-out the MF term
    IF (bl_mynn_edmf_mom == 0) THEN
       onoff=zero
    ELSE
       onoff=one
    ENDIF

    ! USTAR/WSPD make sure it does not blow up when WSPD = 0.
    IF (wspd .le. 1e-6) THEN
       ustovwsp=zero
    ELSE 
       ustovwsp=ust/wspd
    ENDIF
    
    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhosfc     = psfc/(R_d*(tk(kts)+p608*qv(kts)))
    dtz(kts)   = delt/dz(kts)
    rhoz(kts)  = rho(kts)
    rhoinv(kts)= one/rho(kts)
    khdz(kts)  = rhoz(kts)*dfh(kts)
    kmdz(kts)  = rhoz(kts)*dfm(kts)
    DO k=kts+1,kte
       dtz(k)   = delt/dz(k)
       rhoz(k)  = (rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4_kind_phys)
       rhoinv(k)= one/MAX(rho(k),1E-4_kind_phys)
       dzk      = p5 *( dz(k)+dz(k-1) )
       khdz(k)  = rhoz(k)*dfh(k)
       kmdz(k)  = rhoz(k)*dfm(k)
    ENDDO
    rhoz(kte+1)= rhoz(kte)
    khdz(kte+1)= rhoz(kte+1)*dfh(kte)
    kmdz(kte+1)= rhoz(kte+1)*dfm(kte)

    !delta-p for the moisture check
    !delp(kts)  = psfc - (p(kts+1)*dz(kts) + p(kts)*dz(kts+1))/(dz(kts)+dz(kts+1))
    DO k=kts,kte !kts+1,kte-1
       !delp(k)  = (p(k)*dz(k-1) + p(k-1)*dz(k))/(dz(k)+dz(k-1)) - &
       !           (p(k+1)*dz(k) + p(k)*dz(k+1))/(dz(k)+dz(k+1))
       delp(k) = rho(k)*grav*dz(k)
    ENDDO
    !delp(kte)  =delp(kte-1)
    if ( delp(kts) < p333*delp(kts+1) )delp(kts)=p333*delp(kts+1)

    !stability criteria for implicit mf
    if (bl_mynn_edmf == 1) then
       DO k=kts+1,kte-1
          khdz(k) = max(khdz(k),  p5*(s_aw1(k) +sd_aw1(k)))
          khdz(k) = max(khdz(k), -p5*(s_aw1(k) -s_aw1(k+1))  &
                                 -p5*(sd_aw1(k)-sd_aw1(k+1)) )
          kmdz(k) = max(kmdz(k),  p5*(s_aw1(k) +sd_aw1(k)))
          kmdz(k) = max(kmdz(k), -p5*(s_aw1(k) -s_aw1(k+1))  &
                                 -p5*(sd_aw1(k)-sd_aw1(k+1)) )
       ENDDO
    endif

    dth(kts:kte) = zero  ! must initialize for moisture_check routine

!!============================================
!! u
!!============================================
 if (bl_mynn_edmf > 1) then

    do k=kts+1,kte-1
       upcont(k)=onoff*(s_awu1(k) - s_aw1(k)*(u(k)*upwind+u(k-1)*(one-upwind)))
       dncont(k)=onoff*(sd_awu1(k)-sd_aw1(k)*(u(k)*upwind+u(k-1)*(one-upwind)))
    enddo
    ! no flux at the top of the atmosphere
    upcont(kte)=zero
    dncont(kte)=zero
    ! upcont(1) and dncont(1) are not used so they don't need to be set

    k=kts
    a(1)=zero
    b(1)=one + dtz(k)*(kmdz(k+1)+rhosfc*ust*ustovwsp)*rhoinv(k) 
    c(1)=    - dtz(k)*kmdz(k+1)*rhoinv(k)
    d(1)=u(k)+ dtz(k)*uoce*ust**2/wspd*rho(k)                      &
             - dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
       b(k)=one+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) 
       c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
       d(k)=u(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k)) &
           &  + sub_u(k)*delt + det_u(k)*delt
    ENDDO

 else !implicit
    
    k=kts
    !rho-weighted (drag in b-vector):
    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
    b(k)=one+dtz(k)*(kmdz(k+1)+rhosfc*ust*ustovwsp)*rhoinv(k)      &
           & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff                &
           & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)                             &
           & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff                &
           & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    d(k)=u(k)  + dtz(k)*uoce*ust**2/wspd                           &
           & - dtz(k)*rhoinv(k)*s_awu1(k+1)*onoff                  &
           & + dtz(k)*rhoinv(k)*sd_awu1(k+1)*onoff                 &
           & + sub_u(k)*delt + det_u(k)*delt

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)                            &
           &  + p5*dtz(k)*rhoinv(k)*s_aw1(k)*onoff                 &
           &  + p5*dtz(k)*rhoinv(k)*sd_aw1(k)*onoff
       b(k)=one+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k)              &
           &  + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*onoff    &
           &  + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))*onoff
       c(k)=  - dtz(k)*kmdz(k+1)*rhoinv(k)                         &
           &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff               &
           &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
       d(k)=u(k) + dtz(k)*rhoinv(k)*(s_awu1(k)-s_awu1(k+1))*onoff  &
           &  - dtz(k)*rhoinv(k)*(sd_awu1(k)-sd_awu1(k+1))*onoff   &
           &  + sub_u(k)*delt + det_u(k)*delt
    enddo

 endif
 
!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

!! specified gradient at the top 
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradu_top*dztop

!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=u(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
!       du(k)=(d(k-kts+1)-u(k))/delt
       du(k)=(x(k)-u(k))/delt
    ENDDO

!!============================================
!! v
!!============================================

 if (bl_mynn_edmf > 1) then

    DO k=kts+1,kte-1
       upcont(k)=onoff*(s_awv1(k)-s_aw1(k)*(v(k)*upwind+v(k-1)*(one-upwind)))
       dncont(k)=onoff*(sd_awv1(k)-sd_aw1(k)*(v(k)*upwind+v(k-1)*(one-upwind)))
    ENDDO
    ! no flux at the top of the atmosphere
    upcont(kte)=0.
    dncont(kte)=0.
    ! upcont(1) and dncont(1) are not used so they don't need to be set

    k=kts
    a(1)=zero
    b(1)=one + dtz(k)*(kmdz(k+1)+rhosfc*ust*ustovwsp)*rhoinv(k)
    c(1)=    - dtz(k)*kmdz(k+1)*rhoinv(k)
    d(1)=v(k)+ dtz(k)*voce*ust**2/wspd*rho(k)                      &
             - dtz(k)*(upcont(k+1)+dncont(k+1))

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
       b(k)=one+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k)
       c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
       d(k)=v(k) -dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k)) &
           &  + sub_v(k)*delt + det_v(k)*delt
    ENDDO

 else !implicit

    k=kts
    !rho-weighted (drag in b-vector):
    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
    b(k)=one+dtz(k)*(kmdz(k+1) + rhosfc*ust**2/wspd)*rhoinv(k)    &
        &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff                 &
        &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)                            &
        &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff                 &
        &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    d(k)=v(k)  + dtz(k)*voce*ust**2/wspd                          &
        &  - dtz(k)*rhoinv(k)*s_awv1(k+1)*onoff                   &
        &  + dtz(k)*rhoinv(k)*sd_awv1(k+1)*onoff                  &
        &  + sub_v(k)*delt + det_v(k)*delt

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)                           &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)*onoff                   &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)*onoff
       b(k)=one+dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k)              &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*onoff      &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))*onoff
       c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)                         &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff                 &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
       d(k)=v(k) + dtz(k)*rhoinv(k)*(s_awv1(k)-s_awv1(k+1))*onoff &
         & - dtz(k)*rhoinv(k)*(sd_awv1(k)-sd_awv1(k+1))*onoff     &
         & + sub_v(k)*delt + det_v(k)*delt
    enddo
 endif
!! no flux at the top
!    a(kte)=-1.
!    b(kte)=one
!    c(kte)=zero
!    d(kte)=zero

!! specified gradient at the top
!    a(kte)=-1.
!    b(kte)=one
!    c(kte)=zero
!    d(kte)=gradv_top*dztop

!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=v(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
!       dv(k)=(d(k-kts+1)-v(k))/delt
       dv(k)=(x(k)-v(k))/delt
    ENDDO

!!============================================
!! thl tendency
!!============================================
 if (bl_mynn_edmf > 1) then

    DO k=kts+1,kte-1
       upcont(k)= s_awthl1(k)- s_aw1(k)*(thl(k)*upwind+thl(k-1)*(one-upwind))
       dncont(k)=sd_awthl1(k)-sd_aw1(k)*(thl(k)*upwind+thl(k-1)*(one-upwind))
    ENDDO
    ! no flux at the top of the atmosphere
    upcont(kte)=zero
    dncont(kte)=zero
    ! upcont(1) and dncont(1) are not used so they don't need to be set

    k=kts
    a(1)=zero
    b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
    c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
    d(1)=thl(k)+dtz(k)*rhosfc*flt*rhoinv(k) + tcd(k)*delt              &
         - dtz(k)*(upcont(k+1)+dncont(k+1))                            &
         + diss_heat(k)*delt + sub_thl(k)*delt + det_thl(k)*delt

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
       b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=thl(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k)) &
           &  + diss_heat(k)*delt + sub_thl(k)*delt + det_thl(k)*delt
    ENDDO

 else !implicit
    
    k=kts
!rho-weighted: rhosfc*x*rhoinv(k)
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)             &
       &   - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       &   - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                       &
       &   - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       &   - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=thl(k) + dtz(k)*rhosfc*flt*rhoinv(k) + tcd(k)*delt  &
!    d(k)=thl_tot1(k) + dtz(k)*rhosfc*flt*rhoinv(k) + tcd(k)*delt  &
       &   - dtz(k)*rhoinv(k)*s_awthl1(k+1)                  &
       &   + dtz(k)*rhoinv(k)*sd_awthl1(k+1)                 &
       &   + diss_heat(k)*delt + sub_thl(k)*delt + det_thl(k)*delt

    do k=kts+1,kte-1
       a(k)= -dtz(k)*khdz(k)*rhoinv(k)                       &
       &   + p5*dtz(k)*rhoinv(k)*s_aw1(k)                    &
       &   + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)          &
       &   + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))       &
       &   + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)= -dtz(k)*khdz(k+1)*rhoinv(k)                     &
       &   - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       &   - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=thl(k) + tcd(k)*delt                             &
!       d(k)=thl_tot1(k) + tcd(k)*delt                             &
       &   + dtz(k)*rhoinv(k)*(s_awthl1(k)-s_awthl1(k+1))    &
       &   - dtz(k)*rhoinv(k)*(sd_awthl1(k)-sd_awthl1(k+1))  &
       &   +     diss_heat(k)*delt                           &
       &   +     sub_thl(k)*delt + det_thl(k)*delt
    enddo
 endif
!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

!! specified gradient at the top
!assume gradthl_top=gradth_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradth_top*dztop

!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=thl(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !thl(k)=d(k-kts+1)
       thl(k)=x(k)
    ENDDO

IF (bl_mynn_mixqt > 0) THEN
 !============================================
 ! MIX total water (sqw = sqc + sqv + sqi)
 ! NOTE: no total water tendency is output; instead, we must calculate
 !       the saturation specific humidity and then 
 !       subtract out the moisture excess (sqc & sqi)
 !============================================
 if (bl_mynn_edmf > 1) then

    DO k=kts+1,kte-1
       upcont(k)= s_awqt1(k)- s_aw1(k)*(sqw(k)*upwind+sqw(k-1)*(one-upwind))
       dncont(k)=sd_awqt1(k)-sd_aw1(k)*(sqw(k)*upwind+sqw(k-1)*(one-upwind))
    ENDDO
    ! no flux at the top of the atmosphere
    upcont(kte)=zero
    dncont(kte)=zero
    ! upcont(1) and dncont(1) are not used so they don't need to be set

    k=kts
    a(1)=zero
    b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
    c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
    d(1)=sqw(k)+dtz(k)*rhosfc*flq*rhoinv(k) + qcd(k)*delt              &
         - dtz(k)*(upcont(k+1)+dncont(k+1))                            &
         + sub_sqv(k)*delt + det_sqv(k)*delt

    DO k=kts+1,kte-1
       a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
       b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=sqw(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k)) &
           &  + sub_sqv(k)*delt + det_sqv(k)*delt
    ENDDO

 else !implicit

    k=kts
    !rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)         &
       & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                 &
       & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                    &
       & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                 &
       & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=sqw(k)  + dtz(k)*rhosfc*flq*rhoinv(k) + qcd(k)*delt &
       &  - dtz(k)*rhoinv(k)*s_awqt1(k+1)                 &
       &  + dtz(k)*rhoinv(k)*sd_awqt1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                   &
       & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                   &
       & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)      &
       & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))      &
       & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                 &
       & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                 &
       & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=sqw(k) + qcd(k)*delt                          &
       & + dtz(k)*rhoinv(k)*(s_awqt1(k)-s_awqt1(k+1))     &
       & - dtz(k)*rhoinv(k)*(sd_awqt1(k)-sd_awqt1(k+1))
    enddo
endif
!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.
!! specified gradient at the top
!assume gradqw_top=gradqv_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradqv_top*dztop
!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=sqw(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqw2)
!    CALL tridiag3(kte,a,b,c,d,sqw2)

!    DO k=kts,kte
!       sqw2(k)=d(k-kts+1)
!    ENDDO
ELSE
    sqw2=sqw
ENDIF

IF (bl_mynn_mixqt == 0) THEN
!============================================
! cloud water ( sqc ). If mixing total water (bl_mynn_mixqt > 0),
! then sqc will be backed out of saturation check (below).
!============================================
 IF (bl_mynn_cloudmix > 0 .AND. FLAG_QC) THEN

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqc1(k)- s_aw1(k)*(sqc(k)*upwind+sqc(k-1)*(one-upwind))
         dncont(k)=sd_awqc1(k)-sd_aw1(k)*(sqc(k)*upwind+sqc(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero
      ! upcont(1) and dncont(1) are not used so they don't need to be set

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=sqc(k)+dtz(k)*rhosfc*flqc*rhoinv(k) + qcd(k)*delt             &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))                        &
          &    + det_sqc(k)*delt
    
      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=sqc(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k)) &
             &  + det_sqc(k)*delt
      ENDDO

    else !implicit

       k=kts
       !rho-weighted:
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
       &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
       &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=sqc(k)  + dtz(k)*rhosfc*flqc*rhoinv(k) + qcd(k)*delt &
       &  - dtz(k)*rhoinv(k)*s_awqc1(k+1)                   &
       &  + dtz(k)*rhoinv(k)*sd_awqc1(k+1)                  &
       &  + det_sqc(k)*delt

       do k=kts+1,kte-1
          a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                  &
          & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                  &
          & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
          b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
          & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))     &
          & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
          c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                &
          & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
          & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
          d(k)=sqc(k) + qcd(k)*delt                         &
          & + dtz(k)*rhoinv(k)*(s_awqc1(k)-s_awqc1(k+1))    &
          & - dtz(k)*rhoinv(k)*(sd_awqc1(k)-sd_awqc1(k+1))  &
          & + det_sqc(k)*delt
       enddo
    endif
    
! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=sqc(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqc2)
!    CALL tridiag3(kte,a,b,c,d,sqc2)

!    DO k=kts,kte
!       sqc2(k)=d(k-kts+1)
!    ENDDO
  ELSE
    !If not mixing clouds, set "updated" array equal to original array
    sqc2=sqc
  ENDIF
ENDIF

IF (bl_mynn_mixqt == 0) THEN
  !============================================
  ! MIX WATER VAPOR ONLY ( sqv ). If mixing total water (bl_mynn_mixqt > 0),
  ! then sqv will be backed out of saturation check (below).
  !============================================

   k=kts
   !limit unreasonably large negative fluxes:
   qvflux = flqv
   if (qvflux < zero) then
      !do not allow specified surface flux to reduce qv below 1e-8 kg/kg
      qvflux = max(qvflux, (min(0.9*sqv(kts) - 1e-8, zero)/dtz(kts)))
   endif

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqv1(k)- s_aw1(k)*(sqv(k)*upwind+sqv(k-1)*(one-upwind))
         dncont(k)=sd_awqv1(k)-sd_aw1(k)*(sqv(k)*upwind+sqv(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero
      ! upcont(1) and dncont(1) are not used so they don't need to be set

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=sqv(k)+dtz(k)*rhosfc*qvflux*rhoinv(k) + qcd(k)*delt             &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))			   &    
          &    + sub_sqv(k)*delt + det_sqv(k)*delt

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=sqv(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k)) &
             &   + sub_sqv(k)*delt + det_sqv(k)*delt
      ENDDO

   else !implicit 

      k=kts
      !rho-weighted:
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
      & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                   &
      & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
      & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                   &
      & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      d(k)=sqv(k)  + dtz(k)*rhosfc*qvflux*rhoinv(k) + qcd(k)*delt &
      & - dtz(k)*rhoinv(k)*s_awqv1(k+1)                    &
      & + dtz(k)*rhoinv(k)*sd_awqv1(k+1)                   &
      & + sub_sqv(k)*delt + det_sqv(k)*delt

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                  &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                  &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))     &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
         d(k)=sqv(k) + qcd(k)*delt                         &
         & + dtz(k)*rhoinv(k)*(s_awqv1(k)-s_awqv1(k+1))    &
         & - dtz(k)*rhoinv(k)*(sd_awqv1(k)-sd_awqv1(k+1))  &
         & + sub_sqv(k)*delt + det_sqv(k)*delt
      enddo
   endif
   
! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

! specified gradient at the top
! assume gradqw_top=gradqv_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradqv_top*dztop

! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=sqv(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqv2)
!    CALL tridiag3(kte,a,b,c,d,sqv2)

!    DO k=kts,kte
!       sqv2(k)=d(k-kts+1)
!    ENDDO
ELSE
    sqv2=sqv
ENDIF

!============================================
! MIX CLOUD ICE ( sqi )                      
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QI) THEN

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= zero !s_awqi1(k)- s_aw1(k)*(sqi(k)*upwind+sqi(k-1)*(one-upwind))
         dncont(k)=sd_awqi1(k)-sd_aw1(k)*(sqi(k)*upwind+sqi(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero
      ! upcont(1) and dncont(1) are not used so they don't need to be set

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=sqi(k)                                         &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=sqi(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
      
      k=kts
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)      &
!      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                 &
!      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      d(k)=sqi(k)                                        &
!      &  - dtz(k)*rhoinv(k)*s_awqi1(k+1)                &
      &  + dtz(k)*rhoinv(k)*sd_awqi1(k+1)

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                &
!         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)               &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)   &
!         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))  &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)              &
!         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)             &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
         d(k)=sqi(k)                                     &
!         & + dtz(k)*rhoinv(k)*(s_awqi1(k)-s_awqi1(k+1)) &
         & - dtz(k)*rhoinv(k)*(sd_awqi1(k)-sd_awqi1(k+1))
      enddo
   endif

!! no flux at the top
!    a(kte)=-1.       
!    b(kte)=1.        
!    c(kte)=0.        
!    d(kte)=0.        

!! specified gradient at the top
!assume gradqw_top=gradqv_top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradqv_top*dztop

!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=sqi(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqi2)
!    CALL tridiag3(kte,a,b,c,d,sqi2)

!    DO k=kts,kte
!       sqi2(k)=d(k-kts+1)
!    ENDDO
ELSE
   sqi2=sqi
ENDIF

!============================================
! MIX SNOW ( sqs )
!============================================
!hard-code to not mix snow
IF (bl_mynn_cloudmix > 0 .AND. .false.) THEN

    k=kts
!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
    d(k)=sqs(k)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=sqs(k)
    ENDDO

!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=sqs(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,sqs2)
!    CALL tridiag3(kte,a,b,c,d,sqs2)

!    DO k=kts,kte
!       sqs2(k)=d(k-kts+1)
!    ENDDO
ELSE
   sqs2=sqs
ENDIF

!!============================================
!! cloud ice number concentration (qni)
!!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNI .AND. &
      bl_mynn_mixnumcon > 0) THEN

   DO k=kts,kte
      qni2(k)=max(qni(k), zero)
      !enforce minimum number concentration
      if (sqi(k) .gt. 1e-12)qni2(k)=max(qni2(k), ni_min)
   ENDDO

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqni1(k)- s_aw1(k)*(qni2(k)*upwind+qni2(k-1)*(one-upwind))
         dncont(k)=sd_awqni1(k)-sd_aw1(k)*(qni2(k)*upwind+qni2(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero
      ! upcont(1) and dncont(1) are not used so they don't need to be set

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=qni2(k)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=qni2(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
      
      k=kts
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      d(k)=qni2(k)                                         &
      &  - dtz(k)*rhoinv(k)*s_awqni1(k+1)                  &
      &  + dtz(k)*rhoinv(k)*sd_awqni1(k+1)

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                  &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                  &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))     &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
         d(k)=qni2(k)                                      &
         & + dtz(k)*rhoinv(k)*(s_awqni1(k)-s_awqni1(k+1))  &
         & - dtz(k)*rhoinv(k)*(sd_awqni1(k)-sd_awqni1(k+1))
      enddo
   endif
   
!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=qni2(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qni2(k)=d(k-kts+1)
       qni2(k)=x(k)
       !enforce minimum number concentration
       if (sqi2(k) .gt. 1e-12)qni2(k)=max(qni2(k), ni_min)
    ENDDO

ELSE
    qni2=qni
ENDIF

!!============================================
!! cloud water number concentration (qnc)     
!! include non-local transport                
!!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNC .AND. &
    bl_mynn_mixnumcon > 0) THEN

   DO k=kts,kte
      qnc2(k)=max(qnc(k),zero)
      !enforce minimum number concentration
      if (sqc(k) .gt. 1e-12)qnc2(k)=max(qnc2(k), nc_min)
   ENDDO

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqnc1(k)- s_aw1(k)*(qnc2(k)*upwind+qnc2(k-1)*(one-upwind))
         dncont(k)=sd_awqnc1(k)-sd_aw1(k)*(qnc2(k)*upwind+qnc2(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero
      ! upcont(1) and dncont(1) are not used so they don't need to be set

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=qnc2(k)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=qnc2(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
    
      k=kts
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      d(k)=qnc2(k)                                         &
      &  - dtz(k)*rhoinv(k)*s_awqnc1(k+1)                  &
      &  + dtz(k)*rhoinv(k)*sd_awqnc1(k+1)

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                  &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                  &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))     &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
         d(k)=qnc2(k)                                      &
         & + dtz(k)*rhoinv(k)*(s_awqnc1(k)-s_awqnc1(k+1))  &
         & - dtz(k)*rhoinv(k)*(sd_awqnc1(k)-sd_awqnc1(k+1))
      enddo
   endif

!! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=qnc2(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnc2(k)=d(k-kts+1)
       qnc2(k)=x(k)
       !enforce a minimum number concentration
       if (sqc2(k) .gt. 1e-12)qnc2(k)=max(qnc2(k), nc_min)
    ENDDO

ELSE
    qnc2=qnc
ENDIF

!============================================
! Water-friendly aerosols ( qnwfa ).
!============================================
IF (FLAG_QNWFA .AND. bl_mynn_mixaerosols > 0) THEN

   do k=kts,kte
      qnwfa2(k)=max(qnwfa(k),zero)
   enddo

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqnwfa1(k)- s_aw1(k)*(qnwfa2(k)*upwind+qnwfa2(k-1)*(one-upwind))
         dncont(k)=sd_awqnwfa1(k)-sd_aw1(k)*(qnwfa2(k)*upwind+qnwfa2(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero
      ! upcont(1) and dncont(1) are not used so they don't need to be set

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=qnwfa2(k)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=qnwfa2(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
      
      k=kts
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)            &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                      &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                       &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                      &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      d(k)=qnwfa2(k)                                           &
      &  - dtz(k)*rhoinv(k)*s_awqnwfa1(k+1)                    &
      &  + dtz(k)*rhoinv(k)*sd_awqnwfa1(k+1)

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                      &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                      &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)         &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))         &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                    &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                    &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
         d(k)=qnwfa2(k)                                        &
         & + dtz(k)*rhoinv(k)*(s_awqnwfa1(k)-s_awqnwfa1(k+1))  &
         & - dtz(k)*rhoinv(k)*(sd_awqnwfa1(k)-sd_awqnwfa1(k+1))
      enddo
   endif
      
! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=qnwfa2(kte)

!    CALL tridiag(kte,a,b,c,d)
    call tridiag2(kte,a,b,c,d,x)
!    call tridiag3(kte,a,b,c,d,x)

    do k=kts,kte
       !qnwfa2(k)=d(k)
       qnwfa2(k)=max(x(k),zero)
    enddo

    !apply liberal bounds
    do k=kts,kte
       !patch for high wfa over water (useful for cold starts)
       if ((xland-1.5).ge.0)then   ! water
          wt       = min(one, max(zero, zw(k)-100._kind_phys)/max(200._kind_phys, pblh)) !0 below 100 m, 1 above 300 m
          wfa_max2 = (one - wt)*6.e8_kind_phys + wt*wfa_max
       else                        ! land
          wfa_max2 = wfa_max
       endif
       aero_min = wfa_min * exp(-zw(k)/wfa_ht)
       aero_max = wfa_max2* exp(-zw(k)/wfa_ht)
       qnwfa2(k)= min(max(aero_min, qnwfa2(k)), aero_max)
    enddo
else
    !if not mixing aerosols, set "updated" array equal to original array
    qnwfa2=qnwfa
endif

!============================================
! Ice-friendly aerosols ( qnifa ).
!============================================
IF (FLAG_QNIFA .AND. bl_mynn_mixaerosols > 0) THEN

   do k=kts,kte
      qnifa2(k)=max(qnifa(k),zero)
   enddo

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqnifa1(k)- s_aw1(k)*(qnifa2(k)*upwind+qnifa2(k-1)*(one-upwind))
         dncont(k)=sd_awqnifa1(k)-sd_aw1(k)*(qnifa2(k)*upwind+qnifa2(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=qnifa2(k)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=qnifa2(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
    
      k=kts
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)            &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                      &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                       &
      &  - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                      &
      &  - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
      d(k)=qnifa2(k)                                           &
      &  - dtz(k)*rhoinv(k)*s_awqnifa1(k+1)                    &
      &  + dtz(k)*rhoinv(k)*sd_awqnifa1(k+1)

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                      &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)                      &
         & + p5*dtz(k)*rhoinv(k)*sd_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)         &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))         &
         & + p5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                    &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)                    &
         & - p5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
         d(k)=qnifa2(k)                                        &
         & + dtz(k)*rhoinv(k)*(s_awqnifa1(k)-s_awqnifa1(k+1))  &
         & - dtz(k)*rhoinv(k)*(sd_awqnifa1(k)-sd_awqnifa1(k+1))
      enddo
   endif
   
! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=qnifa2(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnifa2(k)=d(k-kts+1)
       qnifa2(k)=max(x(k),zero)
    ENDDO

    !apply bounds
    DO k=kts,kte
       aero_min = ifa_min * exp(-zw(k)/ifa_ht)
       aero_max = ifa_max * exp(-zw(k)/ifa_ht)
       qnifa2(k)= min(max(aero_min, qnifa2(k)), aero_max)
    ENDDO
    
ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnifa2=qnifa
ENDIF

!============================================
! Black-carbon aerosols ( qnbca ).           
!============================================
IF (FLAG_QNBCA .AND. bl_mynn_mixaerosols > 0) THEN

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awqnbca1(k)- s_aw1(k)*(qnbca(k)*upwind+qnbca(k-1)*(one-upwind))
         dncont(k)=zero !sd_awqnbca1(k)-sd_aw1(k)*(qnbca(k)*upwind+qnbca(k-1)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=qnbca(k)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=qnbca(k)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
   
      k=kts
      a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k)          &
      & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*nonloc
      c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                       &
      & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*nonloc
      d(k)=qnbca(k)  - dtz(k)*rhoinv(k)*s_awqnbca1(k+1)*nonloc

      do k=kts+1,kte-1
         a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                      &
         & + p5*dtz(k)*rhoinv(k)*s_aw1(k)*nonloc
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)         &
         & + p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*nonloc
         c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                    &
         & - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)*nonloc
         d(k)=qnbca(k) + dtz(k)*rhoinv(k)*(s_awqnbca1(k)-s_awqnbca1(k+1))*nonloc
      enddo
   endif
   
! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=qnbca(kte)

!    CALL tridiag(kte,a,b,c,d)
   CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnbca2(k)=d(k-kts+1)
       qnbca2(k)=max(zero, x(k))
    ENDDO

ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnbca2=qnbca
ENDIF

!============================================
! Ozone - local mixing only
!============================================
IF (FLAG_OZONE) THEN
    k=kts

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
    d(k)=ozone(k)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=ozone(k)
    ENDDO

! prescribed value
    a(kte)=zero
    b(kte)=one
    c(kte)=zero
    d(kte)=ozone(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !ozone2(k)=d(k-kts+1)
       dozone(k)=(max(x(k),zero) - ozone(k))/delt
    ENDDO
ELSE
    dozone(:)=zero
ENDIF

!!============================================
!! Compute tendencies and convert to mixing ratios for WRF.
!! Note that the momentum tendencies are calculated above.
!!============================================

   IF (bl_mynn_mixqt > 0) THEN 
      DO k=kts,kte
         !compute updated theta using updated thl and old condensate
         th_new = thl(k) + xlvcp/exner(k)*sqc(k) &
           &             + xlscp/exner(k)*sqi(k)

         t  = th_new*exner(k)
         qsat = qsat_blend(t,p(k)) 
         !SATURATED VAPOR PRESSURE
         !esat=esat_blend(t)
         !SATURATED SPECIFIC HUMIDITY
         !qsl=ep_2*esat/(p(k)-ep_3*esat)
         !qsl=ep_2*esat/max(1.e-4,(p(k)-ep_3*esat))

         IF (sqc(k) > zero .or. sqi(k) > zero) THEN !initially saturated
            sqv2(k) = MIN(sqw2(k),qsat)
            portion_qc = sqc(k)/(sqc(k) + sqi(k))
            portion_qi = sqi(k)/(sqc(k) + sqi(k))
            condensate = MAX(sqw2(k) - qsat, zero)
            sqc2(k) = condensate*portion_qc
            sqi2(k) = condensate*portion_qi
         ELSE                     ! initially unsaturated -----
            sqv2(k) = sqw2(k)     ! let microphys decide what to do
            sqi2(k) = zero         ! if sqw2 > qsat 
            sqc2(k) = zero
         ENDIF
      ENDDO
   ENDIF


    !=====================
    ! WATER VAPOR TENDENCY
    !=====================
    DO k=kts,kte
       Dqv(k)=(sqv2(k) - sqv(k))/delt
       !if (sqv2(k) < zero)print*,"neg qv:",sqv2(k),k
    ENDDO

    IF (bl_mynn_cloudmix > 0) THEN
      !=====================
      ! CLOUD WATER TENDENCY
      !=====================
      !print*,"FLAG_QC:",FLAG_QC
      IF (FLAG_QC) THEN
         DO k=kts,kte
            Dqc(k)=(sqc2(k) - sqc(k))/delt
            !if (sqc2(k) < zero)print*,"neg qc:",sqc2(k),k
         ENDDO
      ELSE
         DO k=kts,kte
           Dqc(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD WATER NUM CONC TENDENCY
      !===================
      IF (FLAG_QNC .AND. bl_mynn_mixnumcon > 0) THEN
         DO k=kts,kte
           Dqnc(k) = (qnc2(k)-qnc(k))/delt
           !IF(Dqnc(k)*delt + qnc(k) < 0.)Dqnc(k)=-qnc(k)/delt
         ENDDO 
      ELSE
         DO k=kts,kte
           Dqnc(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD ICE TENDENCY
      !===================
      IF (FLAG_QI) THEN
         DO k=kts,kte
           Dqi(k)=(sqi2(k) - sqi(k))/delt
           !if (sqi2(k) < zero)print*,"neg qi:",sqi2(k),k
         ENDDO
      ELSE
         DO k=kts,kte
           Dqi(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD SNOW TENDENCY
      !===================
      IF (.false.) THEN !disabled
         DO k=kts,kte
           Dqs(k)=(sqs2(k) - sqs(k))/delt
         ENDDO
      ELSE
         DO k=kts,kte
           Dqs(k) = 0.
         ENDDO
      ENDIF

      !===================
      ! CLOUD ICE NUM CONC TENDENCY
      !===================
      IF (FLAG_QNI .AND. bl_mynn_mixnumcon > 0) THEN
         DO k=kts,kte
           Dqni(k)=(qni2(k)-qni(k))/delt
           !IF(Dqni(k)*delt + qni(k) < 0.)Dqni(k)=-qni(k)/delt
         ENDDO
      ELSE
         DO k=kts,kte
           Dqni(k)=0.
         ENDDO
      ENDIF
    ELSE !-MIX CLOUD SPECIES?
      !CLOUDS ARE NOT NIXED (when bl_mynn_cloudmix == 0)
      DO k=kts,kte
         Dqc(k) =0.
         Dqnc(k)=0.
         Dqi(k) =0.
         Dqni(k)=0.
         Dqs(k) =0.
      ENDDO
    ENDIF

    !ensure non-negative moist species
    CALL moisture_check(kte, delt, delp, exner,        &
                        sqv2, sqc2, sqi2, sqs2, thl,   &
                        dqv, dqc, dqi, dqs, dth        )

    !=====================
    ! OZONE TENDENCY CHECK
    !=====================
    DO k=kts,kte
       IF(Dozone(k)*delt + ozone(k) < 0.) THEN
         Dozone(k)=-ozone(k)*0.99/delt
       ENDIF
    ENDDO

    !===================
    ! THETA TENDENCY
    !===================
    IF (FLAG_QI) THEN
      DO k=kts,kte
         Dth(k)=(thl(k) + xlvcp/exner(k)*sqc(k)          &
           &            + xlscp/exner(k)*(sqi(k))        & !+sqs2(k)) &
           &            - th(k))/delt
         !Dth(k)=(thl(k) + xlvcp/exner(k)*qc_tot1(k)       &
         !  &            + xlscp/exner(k)*qi_tot1(k)       &
         !  &            - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy:
         !Dth(k)=(thl(k)*(one+xlvcp/MAX(tk(k),TKmin)*sqc(k)  &
         !  &               + xlscp/MAX(tk(k),TKmin)*sqi(k)) &
         !  &               - th(k))/delt
      ENDDO
    ELSE
      DO k=kts,kte
         Dth(k)=(thl(k)+xlvcp/exner(k)*sqc2(k) - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy.
         !Dth(k)=(thl(k)*(one+xlvcp/MAX(tk(k),TKmin)*sqc(k))  &
         !&               - th(k))/delt
      ENDDO
    ENDIF

    !===================
    ! AEROSOL TENDENCIES
    !===================
    IF (FLAG_QNWFA .AND. FLAG_QNIFA .AND. &
        bl_mynn_mixaerosols > 0) THEN
       DO k=kts,kte
          !=====================
          ! WATER-friendly aerosols
          !=====================
          Dqnwfa(k)=(qnwfa2(k) - qnwfa(k))/delt
          !=====================
          ! Ice-friendly aerosols
          !=====================
          Dqnifa(k)=(qnifa2(k) - qnifa(k))/delt
       ENDDO
    ELSE
       DO k=kts,kte
          Dqnwfa(k)=0.
          Dqnifa(k)=0.
       ENDDO
    ENDIF

    !========================
    ! BLACK-CARBON TENDENCIES
    !========================
    IF (FLAG_QNBCA .AND. bl_mynn_mixaerosols > 0) THEN
       DO k=kts,kte
          Dqnbca(k)=(qnbca2(k) - qnbca(k))/delt
       ENDDO
    ELSE
       DO k=kts,kte
          Dqnbca(k)=0.
       ENDDO
    ENDIF

    !ensure non-negative moist species
    !note: if called down here, dth needs to be updated, but
    !      if called before the theta-tendency calculation, do not compute dth
    !CALL moisture_check(kte, delt, delp, exner,     &
    !                    sqv, sqc, sqi, thl,         &
    !                    dqv, dqc, dqi, dth )

    if (debug_code) then
       problem = .false.
       do k=kts,kte
          wsp  = sqrt(u(k)**2 + v(k)**2)
          wsp2 = sqrt((u(k)+du(k)*delt)**2 + (v(k)+du(k)*delt)**2)
          th2  = th(k) + Dth(k)*delt
          tk2  = th2*exner(k)
          if (wsp2 > 200. .or. tk2 > 360. .or. tk2 < 160.) then
             problem = .true.
             print*,"After tendencies problem at: i=",i," k=",k
             print*," wsp=",wsp," updated wsp=",wsp2
             print*," T=",th(k)*exner(k)," updated T:",tk2
             print*," du=",du(k)*delt," dv=",dv(k)*delt," dth=",dth(k)*delt
             print*," km=",kmdz(k)*dz(k)," kh=",khdz(k)*dz(k)
             print*," u*=",ust," wspd=",wspd
             print*," rhosfc=",rhosfc," delp=",delp(k)
             print*," LH=",flq*rhosfc*1004.," HFX=",flt*rhosfc*1004.
             print*," flq=",flq," flt=",flt," exner=",exner(k)
             print*," drag term=",ust**2/wspd*dtz(k)*rhosfc/rho(kts)
             kproblem = k
          endif
       enddo
       if (problem) then
          print*,"==thl:",thl(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qv:",sqv2(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qc:",sqc2(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"===qi:",sqi2(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====u:",u(max(kproblem-3,1):min(kproblem+3,kte))
          print*,"====v:",v(max(kproblem-3,1):min(kproblem+3,kte))
       endif
    endif

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

  END SUBROUTINE mynn_tendencies

! ==================================================================
  SUBROUTINE moisture_check(kte, delt, dp, exner,   &
                            qv, qc, qi, qs, th,     &
                            dqv, dqc, dqi, dqs, dth )

  ! This subroutine was adopted from the CAM-UW ShCu scheme and 
  ! adapted for use here.
  !
  ! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  ! force them to be larger than minimum value by (1) condensating 
  ! water vapor into liquid or ice, and (2) by transporting water vapor 
  ! from the very lower layer.
  ! 
  ! We then update the final state variables and tendencies associated
  ! with this correction. If any condensation happens, update theta too.
  ! Note that (qv,qc,qi,th) are the final state variables after
  ! applying corresponding input tendencies and corrective tendencies.
    use module_bl_mynnedmf_common, only: xlvcp,xlscp,zero,two,p5,kind_phys
    
    implicit none
    integer,         intent(in)     :: kte
    real(kind_phys), intent(in)     :: delt
    real(kind_phys), dimension(kte), intent(in)     :: dp, exner
    real(kind_phys), dimension(kte), intent(inout)  :: qv, qc, qi, qs, th
    real(kind_phys), dimension(kte), intent(inout)  :: dqv, dqc, dqi, dqs, dth
    integer:: k
    real(kind_phys)::  dqc2, dqi2, dqs2, dqv2, sum, aa, dum
    real(kind_phys), parameter :: qvmin = 1e-20,       &
                                  qcmin = zero,         &
                                  qimin = zero

    do k = kte, 1, -1  ! From the top to the surface
       dqc2 = max(zero, qcmin-qc(k)) !qc deficit (>=0)
       dqi2 = max(zero, qimin-qi(k)) !qi deficit (>=0)
       dqs2 = max(zero, qimin-qs(k)) !qs deficit (>=0)

       !fix tendencies
       dqc(k) = dqc(k) +  dqc2/delt
       dqi(k) = dqi(k) +  dqi2/delt
       dqs(k) = dqs(k) +  dqs2/delt
       dqv(k) = dqv(k) - (dqc2+dqi2+dqs2)/delt
       dth(k) = dth(k) + xlvcp/exner(k)*(dqc2/delt) + &
                         xlscp/exner(k)*((dqi2+dqs2)/delt)
       !update species
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qs(k)  = qs(k)  +  dqs2
       qv(k)  = qv(k)  -  dqc2 - dqi2 - dqs2
       th(k)  = th(k)  +  xlvcp/exner(k)*dqc2 + &
                          xlscp/exner(k)*(dqi2+dqs2)

       !then fix qv
       dqv2   = max(zero, qvmin-qv(k)) !qv deficit (>=0)
       dqv(k) = dqv(k) + dqv2/delt
       qv(k)  = qv(k)  + dqv2
       if ( k .ne. 1 ) then
          !print*,k,"dqv2=",dqv2," dp(k)/dp(k-1)=",dp(k)/dp(k-1)
          qv(k-1)   = qv(k-1)  - dqv2*dp(k)/dp(k-1)
          dqv(k-1)  = dqv(k-1) - dqv2*dp(k)/dp(k-1)/delt
       endif
       qv(k) = max(qv(k),qvmin)
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
       qs(k) = max(qs(k),qimin)
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
                   dqv(k) = dqv(k) - dum/delt
               endif
            enddo
        else
        ! For testing purposes only (not yet found in any output):
        !    write(*,*) 'Full moisture conservation is impossible'
        endif
    endif

    return

  END SUBROUTINE moisture_check

! ==================================================================

SUBROUTINE mynn_mix_chem(kts,kte,i,     &
     delt,dz,pblh,                      &
     nchem, ndvel,                      &
     chem1, vd1, settle1,               &
     rho,                               &
     flt, tcd, qcd,                     &
     dfh,                               &
     s_aw1, s_awchem1,                  &
     emis_ant_no, frp,                  &
     enh_mix,                           &
     bl_mynn_edmf, upwind               )

!-------------------------------------------------------------------
  use module_bl_mynnedmf_common, only: zero,one,two,p5,kind_phys
  
integer, intent(in) :: kts,kte,i,bl_mynn_edmf
real(kind_phys), dimension(kts:kte), intent(in) :: dfh,dz,tcd,qcd
real(kind_phys), dimension(kts:kte), intent(in) :: rho
real(kind_phys), intent(in)    :: flt
real(kind_phys), intent(in)    :: delt,pblh,upwind
integer, intent(in) :: nchem, ndvel
real(kind_phys), dimension( kts:kte+1), intent(in) :: s_aw1
real(kind_phys), dimension( kts:kte, nchem ), intent(inout) :: chem1
real(kind_phys), dimension( kts:kte, nchem ), intent(in) :: settle1
real(kind_phys), dimension( kts:kte+1,nchem), intent(in) :: s_awchem1
real(kind_phys), dimension( ndvel ), intent(in) :: vd1
real(kind_phys), intent(in) :: emis_ant_no,frp
logical, intent(in) :: enh_mix
!local vars
real(kind_phys), dimension(kts:kte) :: dtz,upcont,dncont
real(kind_phys), dimension(kts:kte) :: a,b,c,d,x
real(kind_phys):: rhs,dztop
real(kind_phys):: t,dzk
real(kind_phys):: hght 
real(kind_phys):: khdz_old, khdz_back
integer :: k,kk,kmaxfire                         ! JLS 12/21/21
integer :: ic  ! Chemical array loop index
    
integer, SAVE :: icall

real(kind_phys), dimension(kts:kte) :: rhoinv
real(kind_phys), dimension(kts:kte+1) :: rhoz,khdz
real(kind_phys), parameter :: NO_threshold    = 0.1      ! For anthropogenic sources
real(kind_phys), parameter :: frp_threshold   = 10.0     ! RAR 02/11/22: I increased the frp threshold to enhance mixing over big fires
real(kind_phys), parameter :: pblh_threshold  = 100.0

dztop=.5*(dz(kte)+dz(kte-1))

DO k=kts,kte
   dtz(k)=delt/dz(k)
ENDDO

!Prepare "constants" for diffusion equation.
!khdz = rho*Kh/dz = rho*dfh
rhoz(kts)  =rho(kts)
rhoinv(kts)=one/rho(kts)
khdz(kts)  =rhoz(kts)*dfh(kts)

DO k=kts+1,kte
   rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
   rhoz(k)  =  MAX(rhoz(k),1E-4_kind_phys)
   rhoinv(k)=one/MAX(rho(k),1E-4_kind_phys)
   dzk      = p5  *( dz(k)+dz(k-1) )
   khdz(k)  = rhoz(k)*dfh(k)
ENDDO
rhoz(kte+1)=rhoz(kte)
khdz(kte+1)=rhoz(kte+1)*dfh(kte)

if (bl_mynn_edmf == 1) then
   !stability criteria for mf
   DO k=kts+1,kte-1
      khdz(k) = MAX(khdz(k),  p5*s_aw1(k))
      khdz(k) = MAX(khdz(k), -p5*(s_aw1(k)-s_aw1(k+1)))
   ENDDO
endif

!Enhanced mixing over fires
IF ( enh_mix ) THEN
   DO k=kts+1,kte-1
      khdz_old  = khdz(k)
      khdz_back = pblh * 0.15_kind_phys / dz(k)
      !Modify based on anthropogenic emissions of NO and FRP
      IF ( pblh < pblh_threshold ) THEN
         IF ( emis_ant_no > NO_threshold ) THEN
            khdz(k) = MAX(1.1*khdz(k),sqrt((emis_ant_no / NO_threshold)) / dz(k) * rhoz(k)) ! JLS 12/21/21
!            khdz(k) = MAX(khdz(k),khdz_back)
         ENDIF
         IF ( frp > frp_threshold ) THEN
            kmaxfire = ceiling(log(frp))
            khdz(k) = MAX(1.1*khdz(k), (one - k/(kmaxfire*two)) * ((log(frp))**2 - two*log(frp)) / dz(k)*rhoz(k)) ! JLS 12/21/21
!            khdz(k) = MAX(khdz(k),khdz_back)
         ENDIF
      ENDIF
   ENDDO
ENDIF

!============================================
! Patterned after mixing of water vapor in mynn_tendencies.
!============================================

DO ic = 1,nchem

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awchem1(k,ic)- s_aw1(k)*(chem1(k,ic)*upwind+chem1(k-1,ic)*(one-upwind))
         dncont(k)= zero !sd_awchem1(k,ic)-sd_aw1(k)*(chem1(k,ic)*upwind+chem1(k-1,ic)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=chem1(k,ic)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=chem1(k,ic)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit
      
      k=kts
      a(k)=   -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)
      c(k)=   -dtz(k)*khdz(k+1)*rhoinv(k)           - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)
      d(k)=chem1(k,ic) & !dtz(k)*flt  !neglecting surface sources 
           & - dtz(k)*vd1(ic)*chem1(k,ic) &
           & - dtz(k)*rhoinv(k)*s_awchem1(k+1,ic)

      DO k=kts+1,kte-1
         a(k)=   -dtz(k)*khdz(k)*rhoinv(k)     + p5*dtz(k)*rhoinv(k)*s_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
            &    p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))
         c(k)=   -dtz(k)*khdz(k+1)*rhoinv(k) - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)
         d(k)=chem1(k,ic) + dtz(k)*rhoinv(k)*(s_awchem1(k,ic)-s_awchem1(k+1,ic))
      ENDDO
   endif

   ! prescribed value at top
   a(kte)=zero
   b(kte)=one
   c(kte)=zero
   d(kte)=chem1(kte,ic)

   CALL tridiag3(kte,a,b,c,d,x)

   DO k=kts,kte
      chem1(k,ic)=x(k) + settle1(k,ic)
   ENDDO

ENDDO

END SUBROUTINE mynn_mix_chem

!==================================================================

SUBROUTINE mynn_mix_scalars(            &
     kts,kte,i,                         &
     delt,dz,                           &
     nscalars, scalars,                 &
     rho,                               &
     flt, tcd, qcd,                     &
     dfh,                               &
     s_aw1, s_awscalars1,               &
     bl_mynn_edmf, upwind               )

!-------------------------------------------------------------------
  use module_bl_mynnedmf_common, only: zero,one,p5,kind_phys
  
integer, intent(in) :: kts,kte,i,bl_mynn_edmf
real(kind_phys), dimension(kts:kte), intent(in) :: dfh,dz,tcd,qcd
real(kind_phys), dimension(kts:kte), intent(in) :: rho
real(kind_phys), intent(in)    :: flt
real(kind_phys), intent(in)    :: delt,upwind
integer, intent(in) :: nscalars
real(kind_phys), dimension( kts:kte+1), intent(in) :: s_aw1
real(kind_phys), dimension( kts:kte, nscalars ), intent(inout) :: scalars
real(kind_phys), dimension( kts:kte+1,nscalars), intent(in) :: s_awscalars1
!local vars
real(kind_phys), dimension(kts:kte) :: dtz,upcont,dncont
real(kind_phys), dimension(kts:kte) :: a,b,c,d,x
real(kind_phys):: rhs,dztop
real(kind_phys):: dzk 
integer :: k,ns  !loop indecies

real(kind_phys), dimension(kts:kte) :: rhoinv
real(kind_phys), dimension(kts:kte+1) :: rhoz,khdz

dztop=p5*(dz(kte)+dz(kte-1))
    
DO k=kts,kte
   dtz(k)=delt/dz(k)
ENDDO

!Prepare "constants" for diffusion equation.
!khdz = rho*Kh/dz = rho*dfh
rhoz(kts)  =rho(kts)
rhoinv(kts)=one/rho(kts)
khdz(kts)  =rhoz(kts)*dfh(kts)

DO k=kts+1,kte
   rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
   rhoz(k)  =   MAX(rhoz(k),1E-4_kind_phys)
   rhoinv(k)=one/MAX(rho(k),1E-4_kind_phys)
   dzk      = p5*( dz(k)+dz(k-1) )
   khdz(k)  = rhoz(k)*dfh(k)
ENDDO
rhoz(kte+1)=rhoz(kte)
khdz(kte+1)=rhoz(kte+1)*dfh(kte)

if (bl_mynn_edmf == 1) then
   !stability criteria for mf
   DO k=kts+1,kte-1
      khdz(k) = MAX(khdz(k),  p5*s_aw1(k))
      khdz(k) = MAX(khdz(k), -p5*(s_aw1(k)-s_aw1(k+1)))
   ENDDO
endif

!============================================
! Patterned after mixing of water vapor in mynn_tendencies.
!============================================

DO ns = 1,nscalars

   if (bl_mynn_edmf > 1) then

      DO k=kts+1,kte-1
         upcont(k)= s_awscalars1(k,ns)- s_aw1(k)*(scalars(k,ns)*upwind+scalars(k-1,ns)*(one-upwind))
         dncont(k)=zero !sd_awscalars1(k,ns)-sd_aw1(k)*(scalars(k,ns)*upwind+scalars(k-1,ns)*(one-upwind))
      ENDDO
      ! no flux at the top of the atmosphere
      upcont(kte)=zero
      dncont(kte)=zero

      k=kts
      a(1)=zero
      b(1)=one + dtz(k)*khdz(k+1)*rhoinv(k)
      c(1)=    - dtz(k)*khdz(k+1)*rhoinv(k)
      d(1)=scalars(k,ns)                                       &
          &    - dtz(k)*(upcont(k+1)+dncont(k+1))

      DO k=kts+1,kte-1
         a(k)=   - dtz(k)*khdz(k)*rhoinv(k)
         b(k)=one+ dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
         c(k)=   - dtz(k)*khdz(k+1)*rhoinv(k)
         d(k)=scalars(k,ns)-dtz(k)*(upcont(k+1)-upcont(k)+dncont(k+1)-dncont(k))
      ENDDO

   else !implicit 

      k=kts
      a(k)=   -dtz(k)*khdz(k)*rhoinv(k)
      b(k)=one+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)
      c(k)=   -dtz(k)*khdz(k+1)*rhoinv(k)           - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)
      d(k)=scalars(k,ns) - dtz(k)*rhoinv(k)*s_awscalars1(k+1,ns)

      DO k=kts+1,kte-1
         a(k)=   -dtz(k)*khdz(k)*rhoinv(k)     + p5*dtz(k)*rhoinv(k)*s_aw1(k)
         b(k)=one+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
             &    p5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))
         c(k)=   -dtz(k)*khdz(k+1)*rhoinv(k) - p5*dtz(k)*rhoinv(k)*s_aw1(k+1)
         d(k)=scalars(k,ns) + dtz(k)*rhoinv(k)*(s_awscalars1(k,ns)-s_awscalars1(k+1,ns))
      ENDDO
   endif

   ! prescribed value at top
   a(kte)=zero
   b(kte)=one
   c(kte)=zero
   d(kte)=scalars(kte,ns)

   CALL tridiag3(kte,a,b,c,d,x)

   DO k=kts,kte
      scalars(k,ns)=max(1e-12_kind_phys, x(k))
   ENDDO
ENDDO

END SUBROUTINE mynn_mix_scalars
  
! ==================================================================
!>\ingroup gsd_mynn_edmf
  SUBROUTINE retrieve_exchange_coeffs(kts,kte,dfm,dfh,dz,km1,kh1)

!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: zero,p5,kind_phys
    
    integer , intent(in) :: kts,kte

    real(kind_phys), dimension(kts:kte), intent(in)  :: dz,dfm,dfh
    real(kind_phys), dimension(kts:kte), intent(out) :: km1, kh1


    integer :: k
    real(kind_phys):: dzk

    km1(kts)=zero
    kh1(kts)=zero

    do k=kts+1,kte
       dzk   = p5 *( dz(k)+dz(k-1) )
       km1(k)=dfm(k)*dzk
       kh1(k)=dfh(k)*dzk
    enddo

!    km1(kte+1)=0.
!    kh1(kte+1)=0.

  END SUBROUTINE retrieve_exchange_coeffs

! ==================================================================
!>\ingroup gsd_mynn_edmf
  SUBROUTINE tridiag(n,a,b,c,d)

!! to solve system of linear eqs on tridiagonal matrix n times n
!! after Peaceman and Rachford, 1955
!! a,b,c,d - are vectors of order n 
!! a,b,c - are coefficients on the LHS
!! d - is initially RHS on the output becomes a solution vector
!-------------------------------------------------------------------
    use module_bl_mynnedmf_common, only: zero,one,kind_phys
    
    integer, intent(in):: n
    real(kind_phys), dimension(n), intent(in) :: a,b
    real(kind_phys), dimension(n), intent(inout) :: c,d
    
    integer :: i
    real(kind_phys):: p
    real(kind_phys), dimension(n) :: q
    
    c(n)=zero
    q(1)=-c(1)/b(1)
    d(1)=d(1)/b(1)
    
    DO i=2,n
       p=one/(b(i)+a(i)*q(i-1))
       q(i)=-c(i)*p
       d(i)=(d(i)-a(i)*d(i-1))*p
    ENDDO
    
    DO i=n-1,1,-1
       d(i)=d(i)+q(i)*d(i+1)
    ENDDO

  END SUBROUTINE tridiag

! ==================================================================
!>\ingroup gsd_mynn_edmf
      subroutine tridiag2(kte,a,b,c,d,x)

!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!      b - the main diagonal
!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!      d - right part
!      x - the answer
!      kte - number of unknowns (levels)
        use module_bl_mynnedmf_common, only: zero,one,kind_phys
        
        integer,intent(in) :: kte
        real(kind_phys), dimension(kte), intent(in) :: a,b,c,d
        real(kind_phys), dimension(kte), intent(out):: x
        real(kind_phys), dimension(kte)  :: cp1,dp1
        real(kind_phys):: m
        integer :: k

        ! initialize c-prime and d-prime
        cp1(1) = c(1)/b(1)
        dp1(1) = d(1)/b(1)
        ! solve for vectors c-prime and d-prime
        do k = 2,kte
           m = b(k)-cp1(k-1)*a(k)
           cp1(k) = c(k)/m
           dp1(k) = (d(k)-dp1(k-1)*a(k))/m
        enddo
        ! initialize x
        x(kte) = dp1(kte)
        ! solve for x from the vectors c-prime and d-prime
        do k = kte-1, 1, -1
           x(k) = dp1(k)-cp1(k)*x(k+1)
        end do

    end subroutine tridiag2
! ==================================================================
!>\ingroup gsd_mynn_edmf
       subroutine tridiag3(kte,a,b,c,d,x)

!ccccccccccccccccccccccccccccccc                                                                   
! Aim: Inversion and resolution of a tridiagonal matrix                                            
!          A X = D                                                                                 
! Input:                                                                                           
!  a(*) lower diagonal (Ai,i-1)                                                                  
!  b(*) principal diagonal (Ai,i)                                                                
!  c(*) upper diagonal (Ai,i+1)                                                                  
!  d                                                                                               
! Output                                                                                           
!  x     results                                                                                   
!ccccccccccccccccccccccccccccccc                                                                   
        use module_bl_mynnedmf_common, only: zero,one,kind_phys
         
        integer,intent(in)   :: kte
        integer, parameter   :: kts=1
        real(kind_phys), dimension(kte) :: a,b,c,d
        real(kind_phys), dimension(kte), intent(out) :: x
        integer :: k

        do k=kte-1,kts,-1
           d(k)=d(k)-c(k)*d(k+1)/b(k+1)
           b(k)=b(k)-c(k)*a(k+1)/b(k+1)
        enddo

        do k=kts+1,kte
           d(k)=d(k)-a(k)*d(k-1)/b(k-1)
        enddo

        do k=kts,kte
           x(k)=d(k)/b(k)
        enddo

        return
        end subroutine tridiag3

! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates hybrid diagnotic boundary-layer height (PBLH).
!!
!! NOTES ON THE PBLH FORMULATION:
!!Two different diagnostics are calculated, each of them specialized for a
!!distinct regime (stbale and unstable). The two diagnostics are then blended
!!according to low-level stability, where the magnitude of the of the unstable
!!pbl height is used as a proxy for the low-level stability. The unstable pblh uses
!!the delta-theta-increase method, which defines the pbl height as the level at
!!which the potential temperature first exceeds the minimum potential
!!!temperature within the boundary layer by some specified delta (about 1.5 K).
!!When applied to observed temperatures, this method has been shown to produce PBL-
!!height estimates that are unbiased relative to profiler-based
!!estimates (Nielsen-Gammon et al. 2008\cite Nielsen_Gammon_2008),
!!but it can be biased in stable conditions. Two options 
!!are now available for the stable pblh diagnostic, selectable by the
!!stable_method internal parameter. Option 0 uses a TKE-based PBL height,
!!which has been shown by Banta and Pichugina (2008) \cite Pichugina_2008 to be a good estimate of
!!the PBL height in stable LLJ conditions. Option 1 selects a friction-velocity method of
!!Koracin and Berkowicz (1988\cite Koracin_Berkowicz_1988), which is a very simple one-line
!!diagnostic that has been shown to be very reliable strictly within the stable regime
!!(Steeneveld et al. 2007\cite Steeneveld_et_al_2007). Therefore, a blending (hybrid)
!!method is implemented that uses both methods, weighting each higher in their
!!respective specialized regime. 
!>\section gen_get_pblh  GSD get_pblh General Algorithm
!> @{
SUBROUTINE GET_PBLH(KTS,KTE,pblh,thv1,qke1,ust,zw1,dz1,landsea,kpbl)

!---------------------------------------------------------------
  use module_bl_mynnedmf_common, only: zero,one,ten,hundred,p5,kind_phys
  
integer,intent(in) :: KTS,KTE

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

real(kind_phys), intent(out) :: pblh
real(kind_phys), intent(in) :: landsea,ust
real(kind_phys), dimension(kts:kte), intent(in) :: thv1, qke1, dz1
real(kind_phys), dimension(kts:kte+1), intent(in) :: zw1
!local variables
real(kind_phys):: pblh_stable,qtke,qtkem1,wt,maxqke,TKEeps,minthv,cpblh
real(kind_phys):: delt_thv   !delta theta-v; dependent on land/sea point
real(kind_phys), parameter :: sbl_lim  = 200. !below this height (m), the stable PBLH method will dominate the blending.
real(kind_phys), parameter :: sbl_damp = 170. !transition length for blending (m).
integer :: i,j,k,kthv,ktke,kpbl
integer, parameter :: stable_method=1 !0: TKE-based PBLH, 1: ust*700

!Initialize kpbl
kpbl = 2

!> - FIND MIN THETAV IN THE LOWEST 200 M AGL
k = kts+1
kthv = 1
minthv = 9.E9_kind_phys
do while (zw1(k) .le. 200.)
   if (minthv > thv1(k)) then
      minthv = thv1(k)
      kthv = k
   endif
   k = k+1
enddo

!> - FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
pblh=zero
k = kthv+1
if ((landsea-1.5) .ge. zero) then
   ! WATER
   delt_thv = 0.75_kind_phys
else
   ! LAND
   delt_thv = 1.25_kind_phys
endif

pblh=zero
k = kthv+1
do k=kts+1,kte-1
   if (thv1(k) .ge. (minthv + delt_thv)) then
      pblh = zw1(k) - dz1(k-1)* &
      & min((thv1(k)-(minthv + delt_thv))/ &
      & max(thv1(k)-thv1(k-1),1E-6_kind_phys),one)
   endif
   if (k .EQ. kte-1) pblh = zw1(kts+1) !EXIT SAFEGUARD
   if (pblh .gt. zero) exit
enddo
!print*,"IN GET_PBLH:",thsfc,pblh

if (stable_method == 0) then
   !> - for stable boundary layers, use tke method to complement the
   !! thetav-based definition (when the theta-v based pblh is below ~0.5 km).
   !!the tanh weighting function will make the tke-based definition negligible 
   !!when the theta-v-based definition is above ~1 km.
   ktke   = 1
   maxqke = max(qke1(kts),zero)
   !use 5% of tke max (kosovic and curry, 2000; jas)
   !tkeeps = maxtke/20. = maxqke/40.
   tkeeps = maxqke/40._kind_phys
   tkeeps = max(tkeeps,0.01_kind_phys) !0.025) 
   pblh_stable=zero

   k = ktke+1
   do k=kts+1,kte-1
      qtke  =max(p5*qke1(k)  , zero)
      qtkem1=max(p5*qke1(k-1), zero)
      if (qtke .le. tkeeps) then
         pblh_stable = zw1(k) - dz1(k-1)* &
         & min((tkeeps-qtke)/max(qtkem1-qtke, 1e-6_kind_phys), one)
         !in case of near zero tke, set pblh = lowest level.
         pblh_stable = max(pblh_stable,zw1(kts+1))
         !print *,"pblh_stable:",i,pblh_stable, qke1(k)/2., zw1(kts+1)
      endif
      !k = k+1
      if (k .eq. kte-1) pblh_stable = zw1(kts+1) !exit safeguard
      if (pblh_stable .ne. zero) exit
   enddo
   !> - the tke-based pblh can (rarely) become very large 
   !! in grid points with deep convection (> 8 km!),
   !! so a limit is imposed to not let pblh_stable exceed the
   !! theta_v-based pbl height +/- 350 m.
   !! this has no impact on 99% of the domain.
   pblh_stable = min(pblh_stable,pblh+350._kind_phys)
   pblh_stable = max(pblh_stable,max(pblh-350._kind_phys, ten))
   !if in old pool situation, default to theta_v-based def
   if (maxqke <= tkeeps) pblh_stable = pblh
else
   !> - Method based on Koracin and Berkowicz (1988):
   wt          = p5*TANH((pblh - hundred)/hundred) + p5
   cpblh       = 400._kind_phys*(one-wt) + 700._kind_phys*wt
   pblh_stable = max(ust*cpblh, ten)
   !Even though this estimate will dominate in stable conditions, it should
   !be liberally bounded:
   pblh_stable = max(pblh_stable, pblh-sbl_lim)
   pblh_stable = min(pblh_stable, pblh+sbl_lim)
endif

!blend the stable and unstable pblh heights
wt=p5*TANH((pblh - sbl_lim)/sbl_damp) + p5
pblh=pblh_stable*(one-wt) + pblh*wt

!Compute kpbl
do k=kts+1,kte-1
   if ( zw1(k) >= pblh) then
      kpbl = k-1
      exit
   endif
enddo

#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

END SUBROUTINE GET_PBLH
!> @}
  
! ==================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine is the Dynamic Multi-Plume (DMP) Mass-Flux Scheme.
!! 
!! dmp_mf() calculates the nonlocal turbulent transport from the dynamic
!! multiplume mass-flux scheme as well as the shallow-cumulus component of 
!! the subgrid clouds. Note that this mass-flux scheme is called when the
!! namelist paramter \p bl_mynn_edmf is set to 1 (recommended).
!!
!! Much thanks to Kay Suslj of NASA-JPL for contributing the original version
!! of this mass-flux scheme. Considerable changes have been made from it's
!! original form. Some additions include:
!!  -# scale-aware tapering as dx -> 0
!!  -# transport of TKE (extra namelist option)
!!  -# Chaboureau-Bechtold cloud fraction & coupling to radiation
!!  -# some extra limits for numerical stability
!!
!! This scheme remains under development, so consider it experimental code. 
!!
  SUBROUTINE DMP_mf(i,j,                           &
                 & kts,kte,dt,zw1,dz1,pres1,rho1,  &
                 & momentum_opt,                   &
                 & tke_opt,                        &
                 & scalar_opt,                     &
                 & aerosol_opt, numcon_opt,        &
                 & bl_mynn_closure,                &
                 & u1,v1,w1,th1,thl1,thv1,tk1,     &
                 & qt1,qv1,qc1,qke1,qsq1,          &
                 & qnc1,qni1,qnwfa1,qnifa1,qnbca1, &
                 & ex1,vt1,vq1,sgm1,               &
                 & ust,flt,fltv,flq,flqv,          &
                 & pblh,kpbl,dx,landsea,ts,        &
            ! outputs - updraft properties   
                 & edmf_a1,edmf_w1,                &
                 & edmf_qt1,edmf_thl1,             &
                 & edmf_ent1,edmf_qc1,             &
                 & edmf_qv1,edmf_u1,edmf_v1,       &
            ! outputs - variables needed for solver 
                 & s_aw1,s_awthl1,s_awqt1,         &
                 & s_awqv1,s_awqc1,                &
                 & s_awu1,s_awv1,s_awqke1,s_awqsq1,&
                 & s_awqnc1,s_awqni1,              &
                 & s_awqnwfa1,s_awqnifa1,          &
                 & s_awqnbca1,                     &
                 & sub_thl1,sub_sqv1,              &
                 & sub_u1,sub_v1,                  &
                 & det_thl1,det_sqv1,det_sqc1,     &
                 & det_u1,det_v1,                  &
            ! chem/smoke
                 & nchem,chem1,s_awchem1,          &
                 & mix_chem,                       &
            ! generic scalar array
                 & nscalars,scalars,s_awscalars1,  &                 
            ! in/outputs - subgrid scale clouds
                 & qc_bl1,cldfra_bl1,              &
                 & qc_bl1_old,cldfra_bl1_old,      &
            ! inputs - flags for moist arrays
                 & F_QC,F_QI,                      &
                 & F_QNC,F_QNI,                    &
                 & F_QNWFA,F_QNIFA,F_QNBCA,        &
                 & Psig_shcu,                      &
            ! output info
                 & maxwidth,ktop,maxmf,ztop,       &
                 & excess_h,excess_q,              &
            ! inputs for stochastic perturbations
                 & spp_pbl,pattern_spp_pbl1,       &
                 & tkeprod_up,el1                  )

    use module_bl_mynnedmf_common, only: gtr,p608,grav,xlvcp,ep_2,  &
         ep_3,r_v,cp,cpv,tv0,b1,                                    &
         zero,one,two,three,four,five,thirty,hundred,               &
         p2,p3,p333,p4,p5,kind_phys

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

! configuration inputs:
 integer, intent(in) ::                                            &
      kts,                kte,                  kpbl,              &
      momentum_opt,       tke_opt,              scalar_opt,        &
      aerosol_opt,        numcon_opt,                              &
      spp_pbl,            i,                    j
 real(kind_phys), intent(in):: bl_mynn_closure
 real(kind_phys), dimension(kts:kte), intent(in)  ::               &
      pattern_spp_pbl1
! state variables
 real(kind_phys), dimension(kts:kte), intent(in)  ::               &
      &u1,v1,w1,th1,thl1,tk1,qt1,qv1,qc1,                          &
      &ex1,dz1,thv1,pres1,rho1,qke1,qsq1,qnc1,qni1,                &
      &qnwfa1,qnifa1,qnbca1,el1
 real(kind_phys),dimension(kts:kte+1), intent(in) ::               &
      &zw1
 real(kind_phys), intent(in)                      ::               &
      &flt,fltv,flq,flqv,Psig_shcu,                                &
      &landsea,ts,dx,dt,ust,pblh
 logical, optional :: F_QC,F_QI,F_QNC,F_QNI,F_QNWFA,F_QNIFA,F_QNBCA

 ! outputs - updraft properties
 real(kind_phys),dimension(kts:kte), intent(inout) ::              &
      &edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1,      &
      &edmf_u1,edmf_v1,edmf_qv1
 ! add one local edmf variable:
  real(kind_phys),dimension(kts:kte) :: edmf_th1
 ! output
 integer,         intent(inout) :: ktop
 real(kind_phys), intent(inout) :: maxmf,ztop,maxwidth,            &
      &excess_h,excess_q
 ! outputs - variables needed for solver: sum ai*rho*wis_awphi
 real(kind_phys),dimension(kts:kte+1), intent(inout) ::            &
      &s_aw1,s_awthl1,s_awqt1,s_awqv1,s_awqc1,s_awqnc1,s_awqni1,   &
      &s_awqnwfa1,s_awqnifa1,s_awqnbca1,s_awu1,s_awv1,             &
      &s_awqke1,s_awqsq1
 real(kind_phys),dimension(kts:kte+1) :: s_aw2

 real(kind_phys),dimension(kts:kte), intent(inout) ::              &
      &qc_bl1,cldfra_bl1,qc_bl1_old,cldfra_bl1_old,tkeprod_up
 integer, parameter :: nup=8, debug_mf=0
 real(kind_phys)    :: nup2

 !------------- local variables -------------------
 ! updraft properties defined on interfaces (k=1 is the top of the
 ! first model layer
 real(kind_phys),dimension(kts:kte+1,1:NUP) ::                     &
      &UPW,UPTHL,UPQT,UPQC,UPQV,                                   &
      &UPA,UPU,UPV,UPTHV,UPQKE,UPQSQ,UPQNC,                        &
      &UPQNI,UPQNWFA,UPQNIFA,UPQNBCA
 ! entrainment defined as the mass-layer mean
 real(kind_phys),dimension(kts:kte,1:NUP) :: ENT
 ! internal variables
 integer :: k,ip,k50
 real(kind_phys):: fltv2,wstar,qstar,thstar,sigmaW,sigmaQT,        &
      &sigmaTH,z0,pwmin,pwmax,wmin,wmax,wlv,Psig_w,maxw,maxqc,wpbl
 real(kind_phys):: B,QTn,THLn,THVn,QCn,Un,Vn,QKEn,QSQn,QNCn,QNIn,  &
      &  QNWFAn,QNIFAn,QNBCAn,upak,                                &
      &  Wn2,Wn,EntEXP,EntEXM,EntW,BCOEFF,THVkm1,THVk,Pk,rho_int

 ! w parameters-used in the original entrainment form
 real(kind_phys), parameter ::                                     &
      &Wa=2./3.,      Wb=0.002,      Wc=1.5 

 ! Parameters/variables for regulating plumes:
 real(kind_phys), parameter :: Atot = 0.10 ! Maximum total fractional area of all updrafts
 real(kind_phys), parameter :: lmax = 1000.! diameter of largest plume (absolute maximum, can be smaller)
 real(kind_phys), parameter :: lmin = 300. ! diameter of smallest plume (absolute minimum, can be larger)
 real(kind_phys), parameter :: dlmin = 0.  ! delta increase in the diameter of smallest plume (large fltv) 
 real(kind_phys)            :: minwidth    ! actual width of smallest plume
 real(kind_phys)            :: dl          ! variable increment of plume size
 real(kind_phys), parameter :: dcut = 1.2  ! max diameter of plume to parameterize relative to dx (km)
 real(kind_phys)::  d     != -2.3 to -1.7  ;=-1.9 in Neggers paper; power law exponent for number density (N=Cl^d).
      ! Note that changing d to -2.0 makes each size plume equally contribute to the total coverage of all plumes.
      ! Note that changing d to -1.7 doubles the area coverage of the largest plumes relative to the smallest plumes.
 real(kind_phys):: cn,c,l,n,an2,hux,wspd_pbl,cloud_base,           &
      maxwidth_dx,maxwidth_pbl,maxwidth_cld,maxwidth_flx

 ! chem/smoke
 integer, intent(in) :: nchem
 real(kind_phys),dimension(kts:kte,   nchem)       :: chem1
 real(kind_phys),dimension(kts:kte+1, nchem)       :: s_awchem1
 real(kind_phys),dimension(nchem)                  :: chemn
 real(kind_phys),dimension(kts:kte+1,1:NUP, nchem) :: UPCHEM
 integer :: ic,cb_check
 real(kind_phys),dimension(kts:kte,   nchem)       :: edmf_chem
 logical, intent(in) :: mix_chem
 !generic scalars
 integer, intent(in) :: nscalars
 real(kind_phys),dimension(kts:kte,   nscalars)    :: scalars
 real(kind_phys),dimension(kts:kte+1, nscalars)    :: s_awscalars1
 real(kind_phys),dimension(nscalars)               :: scalarsn
 real(kind_phys),dimension(kts:kte+1,1:NUP,nscalars) :: upscalars
 !local ktop for each plume
 integer,dimension(1:NUP) :: ktop_plume
 logical :: superadiabatic

 ! Varaibles for mass flux cloud fraction
 real(kind_phys),dimension(kts:kte), intent(inout) :: vt1, vq1, sgm1
 real(kind_phys):: sigq,xl,rsl,cpm,a,qmq,Aup,Q1,diffqt,qsat_tk,     &
         Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid,entfac,   &
         cf_strat,qc_strat,cf_mf,qc_mf,qc_mf_min,pct_mf,wt2,        &
         dqwdz,tauc,mf_at_cb,qcfac,sigqfac,cffac
 real(kind_phys), parameter :: cf_thresh = 0.5 ! only overwrite stratus CF less than this value

 ! Variables interpolated to interface levels
 real(kind_phys),dimension(kts:kte) :: exneri,dzi,rhoi,qti
 ! Variables for plume interpolation/saturation check 
 real(kind_phys):: THp, QTp, QCp, QCs, esat, qsl
 real(kind_phys):: csigma,acfac,ac_wsp

 !plume overshoot
 integer :: overshoot
 real(kind_phys):: bvf, Frz, dzp

 !Flux limiter: not let mass-flux of heat between k=1&2 exceed (fluxportion)*(surface heat flux).
 !This limiter makes adjustments to the entire column.
 real(kind_phys):: adjustment, flx1, flt2
 real(kind_phys), parameter :: fluxportion=0.75 ! set liberally, so has minimal impact. Note that
                                     ! 0.5 starts to have a noticeable impact
                                     ! over land (decrease maxMF by 10-20%), but no impact over water.
 !----- environmental subsidence--------------------------------------
 !Option to activate environmental subsidence in mass-flux scheme
 logical,parameter:: env_subs   = .false.
 real(kind_phys),dimension(kts:kte) :: sub_thl1,sub_sqv1,         &
         sub_u1,sub_v1,det_thl1,det_sqv1,det_sqc1,det_u1,det_v1,  &  !tendencied due to detrainment
         envm_a,envm_w,envm_thl,envm_sqv,envm_sqc,                &
         envm_u,envm_v                  !environmental variables defined at middle of layer
 real(kind_phys),dimension(kts:kte+1) ::  envi_a,envi_w !environmental variables defined at model interface
 real(kind_phys):: temp,sublim,qc_ent,qv_ent,qt_ent,thl_ent,detrate, &
         detrateUV,oow,exc_fac,aratio,detturb,qc_grid,qc_sgs,        &
         qc_plume,exc_heat,exc_moist,tk_int,tvs,dthvdz,zagl
 real(kind_phys), parameter :: Cdet   = 1./45.
 real(kind_phys), parameter :: dzpmax = 300. !limit dz used in detrainment - can be excessing in thick layers
 !parameter "Csub" determines the propotion of upward vertical velocity that contributes to
 !environmenatal subsidence. Some portion is expected to be compensated by downdrafts instead of
 !gentle environmental subsidence. 1.0 assumes all upward vertical velocity in the mass-flux scheme
 !is compensated by "gentle" environmental subsidence.
 real(kind_phys), parameter :: Csub=0.25

 !Factor for the pressure gradient effects on momentum transport
 real(kind_phys), parameter :: pgfac = 0.00  ! Zhang and Wu showed 0.4 is more appropriate for lower troposphere
 real(kind_phys):: Uk,Ukm1,Vk,Vkm1,dxsa

 ! WA TEST 12/23/24 Vars for PBL-average QKE
 real(kind_phys):: qkebl

 !set for debugging at specific point
 !integer, parameter::idbg = 452, jdbg = 272
      
! Inititialize 2d ouput
 ktop      =0    !integer
 ztop      =zero
 maxmf     =zero
 maxwidth  =zero
 excess_h  =zero
 excess_q  =zero
 ! Initialize individual updraft properties
 upw       =zero
 upthl     =zero
 upthv     =zero
 upqt      =zero
 upa       =zero
 upu       =zero
 upv       =zero
 upqc      =zero
 upqv      =zero
 upqke     =zero
 upqsq     =zero
 upqnc     =zero
 upqni     =zero
 upqnwfa   =zero
 upqnifa   =zero
 upqnbca   =zero
 if ( mix_chem ) then
    upchem(kts:kte+1,1:nup,1:nchem)=zero
 endif
 if ( scalar_opt > 0 ) then
    upscalars(kts:kte+1,1:nup,1:nscalars)=zero
 endif
 ent       =0.001_kind_phys
 ! Initialize mean updraft properties
 edmf_a1   =zero
 edmf_w1   =zero
 edmf_qt1  =zero
 edmf_thl1 =zero
 edmf_ent1 =zero
 edmf_qc1  =zero
 edmf_qv1  =zero
 edmf_u1   =zero
 edmf_v1   =zero
 if ( mix_chem ) then
    edmf_chem(kts:kte,1:nchem) = zero
 endif
 ! Initialize the variables needed for implicit solver
 s_aw1     =zero
 s_awthl1  =zero
 s_awqt1   =zero
 s_awqv1   =zero
 s_awqc1   =zero
 s_awu1    =zero
 s_awv1    =zero
 s_awqke1  =zero
 s_awqsq1  =zero
 s_awqnc1  =zero
 s_awqni1  =zero
 s_awqnwfa1=zero
 s_awqnifa1=zero
 s_awqnbca1=zero
 if ( mix_chem ) then
    s_awchem1(kts:kte+1,1:nchem) = zero
 endif
 if ( scalar_opt > 0 ) then
    s_awscalars1(kts:kte+1,1:nscalars) = zero
 endif
 ! Initialize explicit tendencies for subsidence & detrainment
 sub_thl1 = zero
 sub_sqv1 = zero
 sub_u1   = zero
 sub_v1   = zero
 det_thl1 = zero
 det_sqv1 = zero
 det_sqc1 = zero
 det_u1   = zero
 det_v1   = zero
 nup2     = nup !start with nup, but set to zero if activation criteria fails
 tkeprod_up = zero

 if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
    print*,"===incoming forcing in mf component:"
    print*,"  flt=",flt," fltv=",fltv
    print*,"  Psig_shcu=",Psig_shcu," wspd=",sqrt(max(u1(kts)**2 + v1(kts)**2, 0.01_kind_phys))
 endif
      
 ! Taper off MF scheme when significant resolved-scale motions
 ! are present This function needs to be asymetric...
 maxw        = zero
 cloud_base  = 9000.0_kind_phys
 qkebl       = zero
 do k=1,kte-1
    zagl = zw1(k) + p5*dz1(k)
    if (zagl > (pblh + 500.)) exit

    wpbl = w1(k)
    if (w1(k) < zero)wpbl = two*w1(k)
    maxw = max(maxw,abs(wpbl))

    !Find highest k-level below 50m AGL
    if (zagl <= 50.)k50=k

    !Establish mean tke in pbl for entrainment
    if (k <= kpbl) then
       qkebl = qkebl + qke1(k)
    endif

    !Search for cloud base
    qc_sgs = max(qc1(k), qc_bl1(k))
    if ((qc_sgs > 1E-5) .and. (cldfra_bl1(k) .ge. p5) .and. cloud_base == 9000.0) then
       cloud_base = zw1(k) !height at interface below cloud mass level
    endif
 enddo
 qkebl = qkebl / real(kpbl, kind=kind_phys)
 
 !do nothing for small w (< 1 m/s), but linearly taper off for w > 1.0 m/s
 maxw = max(zero, maxw - one)
 Psig_w = max(zero, one - maxw)
 Psig_w = min(Psig_w, Psig_shcu)

 !Completely shut off MF scheme for strong resolved-scale vertical velocities.
 fltv2 = fltv
 if (Psig_w < 1e-2_kind_phys .and. fltv > zero) fltv2 = -1.*fltv

 if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
    print*,"===criteria for small w in pbl:"
    print*,"maxw=",maxw," Psig_w=",Psig_w," k50=",k50
 endif
      
 ! If surface buoyancy is positive we do integration, otherwise no.
 ! Also, ensure that it is at least slightly superadiabatic up through 50 m
 superadiabatic = .false.
 if ((landsea-1.5).ge.zero) then
    hux = -0.0024  ! WATER   ! dTHv/dz must be < -0.24 K per 100 m.
 else
    hux = -0.0030  ! LAND    ! dTHv/dz must be < -0.30 K per 100 m.
 endif
 tvs    = ts*(one+p608*qv1(kts))
 if ((thv1(kts)-tvs)/(p5*dz1(kts)) < hux) then
    !if superadiabatic at the surface, continue checking all layers below 50 m:
    superadiabatic = .true.
    do k=2,max(2,k50)
       hux = -0.0018  !allow for smaller superadiabatic layers above the surface
       dthvdz = (thv1(k)-thv1(k-1))/(p5*(dz1(k)+dz1(k-1)))
       if (dthvdz < hux) then
          superadiabatic = .true.
       else
          superadiabatic = .false.
          exit
       endif
    enddo
 endif

 if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
    print*," superadiatic=",superadiabatic," hux=",hux
 endif

 ! Determine the numer of updrafts/plumes in the grid column:
 ! Some of these criteria may be a little redundant but useful for bullet-proofing.
 !   (1) largest plume = 1.2 * dx.
 !   (2) Apply a scale-break, assuming no plumes with diameter larger than 1.1*PBLH can exist.
 !   (3) max plume size beneath clouds deck approx = 0.5 * cloud_base.
 !   (4) add wspd-dependent limit, when plume model breaks down. (hurricanes)
 !   (5) limit to reduce max plume sizes in weakly forced conditions. This is only
 !       meant to "soften" the activation of the mass-flux scheme.
 ! Criteria (1)
 maxwidth_dx = min(dx*dcut, lmax)
 !Criteria (2)
 maxwidth_pbl = min(1.1_kind_phys*pblh, lmax) 
 ! Criteria (3)
 if ((landsea-1.5) .lt. zero) then  !land
    maxwidth_cld = min(lmax, max(0.5_kind_phys*cloud_base, 300._kind_phys))
 else                               !water
    maxwidth_cld = min(lmax, max(0.8_kind_phys*cloud_base, 300._kind_phys))
 endif
 ! Criteria (4)
 wspd_pbl=sqrt(max(u1(kts)**2 + v1(kts)**2, 0.01_kind_phys))
 !Note: area fraction (acfac) is modified below
 ! Criteria (5) - function of fltv
 if ((landsea-1.5) .lt. zero) then  !land
    maxwidth_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.040)/0.04) + p5),1000._kind_phys), zero)
 else                             !water
    maxwidth_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.007)/0.02) + p5),1000._kind_phys), zero)
    !width_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.010)/0.025) + .5),1000._kind_phys), zero)
 endif
 maxwidth = MIN(maxwidth_dx, maxwidth_pbl)
 maxwidth = MIN(maxwidth,    maxwidth_cld)
 maxwidth = MIN(maxwidth,    maxwidth_flx)      
 minwidth = lmin

 if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then  !land
    print*,"===limiting factors on plume width:"
    print*,"maxwidth_dx= ",maxwidth_dx," dx=",dx
    print*,"maxwidth_pbl=",maxwidth_pbl," pblh=",pblh
    print*,"maxwidth_cld=",maxwidth_cld," cldbase=",cloud_base
    print*,"maxwidth_flx=",maxwidth_flx," fltv=",fltv
    print*,"===final check for activation criteria:"
    print*,"fltv2=",fltv2," superadiatic=",superadiabatic
    print*,"maxwidth=",maxwidth," minwidth=",minwidth
 endif

 if (maxwidth .le. minwidth) then ! deactivate MF component
    nup2 = 0
    maxwidth = zero
 endif

 !Begin plume processing if passes criteria
 if ( flt > zero .AND. fltv2 > 0.008 .AND. (maxwidth > minwidth) .AND. superadiabatic) then

    ! Find coef C for number size density N
    cn = zero
    d  =-1.9_kind_phys  !set d to value suggested by Neggers 2015 (JAMES).
    dl = (maxwidth - minwidth)/real(nup-1,kind=kind_phys)
    do ip=1,NUP
       ! diameter of plume
       l = minwidth + dl*real(ip-1,kind=kind_phys)
       cn = cn + l**d * (l*l)/(dx*dx) * dl  ! sum fractional area of each plume
    enddo
    C = Atot/cn   !Normalize C according to the defined total fraction (Atot)

    ! Make updraft area (UPA) a function of the buoyancy flux
    if ((landsea-1.5) .lt. zero) then  !land
       acfac = p5*tanh((fltv2 - 0.02_kind_phys)/0.05_kind_phys) + p5
    else
       acfac = p5*tanh((fltv2 - 0.012_kind_phys)/0.03_kind_phys) + p5
    endif
      
    !For hurricane tuning, add a windspeed-dependent adjustment to acfac that tapers off
    !the mass-flux scheme linearly above sfc wind speeds of 13 m/s.
    ac_wsp = one - min((max(wspd_pbl - 13.0_kind_phys, zero))/12._kind_phys, one)
    acfac  = min(acfac, ac_wsp)

    ! Find the portion of the total fraction (Atot) of each plume size:
    An2 = zero
    do ip=1,nup
       ! diameter of plume
       l  = minwidth + dl*real(ip-1,kind=kind_phys)
       N  = C*l**d                           ! number density of plume n
       UPA(1,ip) = N*l*l/(dx*dx) * dl        ! fractional area of plume n

       UPA(1,ip) = UPA(1,ip)*acfac
       An2 = An2 + UPA(1,ip)                 ! total fractional area of all plumes
       !if (upa(i,i)<zero)print*,"Neg plume area: plume size=",l," dl=",dl,&
       !    " area=",UPA(1,I),"; total=",An2," acfac=",acfac," N=",N
    end do

    if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
       print*,"===limiting factors on fractional area coverage:"
       print*,"total araa frac=",An2," atot=",atot
       print*,"upa=",UPA(1,1:nup)
       print*,"acfac fluxes=",acfac," fltv=",fltv
       print*,"acfac wspd=",ac_wsp," wspd_pbl=",wspd_pbl
    endif
    
    ! set initial conditions for updrafts
    z0    =50._kind_phys
    pwmin =0.1_kind_phys    ! was 0.5
    pwmax =0.4_kind_phys    ! was 3.0

    wstar =max(1.E-2_kind_phys,(gtr*fltv2*pblh)**(p333))
    qstar =max(flq,1e-5_kind_phys)/wstar
    thstar=flt/wstar

    csigma=1.34_kind_phys

    if (env_subs) then
       exc_fac = zero
    else
       if ((landsea-1.5).GE.zero) then
          !water: increase factor to compensate for decreased pwmin/pwmax
          !0.58*25 = 14.5
          !0.58*20 = 11.6
          !0.58*16 = 9.28
          !0.58*10 = 5.8
          !0.58*5  = 2.9
          exc_fac = 9.28_kind_phys
       else
          !land: no need to increase factor - already sufficiently large superadiabatic layers
          exc_fac = 0.58_kind_phys
       endif
    endif
    !decrease excess for large wind speeds
    exc_fac = exc_fac * ac_wsp

    !Note: sigmaW is typically about 0.5*wstar
    sigmaW =csigma*wstar*(z0/pblh)**(p333)*(one - 0.8_kind_phys*z0/pblh)
    sigmaQT=csigma*qstar*(z0/pblh)**(p333)
    sigmaTH=csigma*thstar*(z0/pblh)**(p333)

    !Note: Given the pwmin & pwmax set above, these max/mins are
    !      rarely exceeded. 
    wmin=MIN(sigmaW*pwmin,0.1_kind_phys)
    wmax=MIN(sigmaW*pwmax,0.5_kind_phys)

    if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
       print*,"===excess components:"
       print*,"wstar=",wstar,"qstar=",qstar,"thstar=",thstar
       print*,"csigma=",csigma,"exc_fac=",exc_fac,"wmin=",wmin
       print*,"sigmaw=",sigmaw,"sigmaqt=",sigmaqt,"sigmath=",sigmath
    endif

    !SPECIFY UPDRAFT PROPERTIES AT MODEL INTERFACE BETWEEN K = 1 & 2
    do ip=1,NUP
       wlv=wmin+(wmax-wmin)/real(NUP2,kind=kind_phys)*real(ip-1,kind=kind_phys)
       UPA(2,ip)    =UPA(1,ip)
       !SURFACE UPDRAFT VERTICAL VELOCITY
       UPW(2,ip)    =wmin + real(ip,kind=kind_phys)/real(NUP)*(wmax-wmin)
       UPU(2,ip)    =(u1(kts)    *dz1(kts+1)+u1(kts+1)    *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPV(2,ip)    =(v1(kts)    *dz1(kts+1)+v1(kts+1)    *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQC(2,ip)   =zero
       !UPQC(1,ip)  =(qc1(kts)*dz1(kts+1)+qc1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))

       exc_heat     =exc_fac*UPW(2,ip)*sigmaTH/sigmaW
       UPTHV(2,ip)  =(thv1(kts)  *dz1(kts+1)+thv1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1)) &
           &        + exc_heat
       UPTHL(2,ip)  =(thl1(kts)  *dz1(kts+1)+thl1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1)) &
           &        + exc_heat
       !calculate exc_moist by use of surface fluxes
       exc_moist    =exc_fac*UPW(2,ip)*sigmaQT/sigmaW
       UPQT(2,ip)   =(qt1(kts)   *dz1(kts+1)+qt1(kts+1)   *dz1(kts))/(dz1(kts)+dz1(kts+1))&
            &       + exc_moist
       UPQKE(2,ip)  =(qke1(kts)  *dz1(kts+1)+qke1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQSQ(2,ip)  =(qsq1(kts)  *dz1(kts+1)+qsq1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNC(2,ip)  =(qnc1(kts)  *dz1(kts+1)+qnc1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNI(2,ip)  =(qni1(kts)  *dz1(kts+1)+qni1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNWFA(2,ip)=(qnwfa1(kts)*dz1(kts+1)+qnwfa1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNIFA(2,ip)=(qnifa1(kts)*dz1(kts+1)+qnifa1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNBCA(2,ip)=(qnbca1(kts)*dz1(kts+1)+qnbca1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))
    enddo

    if ( mix_chem ) then
      do ip=1,NUP
        do ic = 1,nchem
          UPCHEM(2,ip,ic)=(chem1(kts,ic)*dz1(kts+1)+chem1(kts+1,ic)*dz1(kts))/(dz1(kts)+dz1(kts+1))
        enddo
      enddo
    endif
    if ( scalar_opt > 0 ) then
      do ip=1,NUP
        do ic = 1,nscalars
          upscalars(2,ip,ic)=(scalars(kts,ic)*dz1(kts+1)+scalars(kts+1,ic)*dz1(kts))/(dz1(kts)+dz1(kts+1))
        enddo
      enddo
   endif
   
    !Initialize environmental variables which can be modified by detrainment
    envm_thl(kts:kte)=thl1(kts:kte)
    envm_sqv(kts:kte)=qv1(kts:kte)
    envm_sqc(kts:kte)=qc1(kts:kte)
    envm_u(kts:kte)  =u1(kts:kte)
    envm_v(kts:kte)  =v1(kts:kte)
    !Interpolate some mass-layer variables to interface variables
    rhoi(kts)   = rho1(kts)
    dzi(kts)    = p5*dz1(kts)
    exneri(kts) = ex1(kts)
    do k=kts+1,kte
       rhoi(k)  = (rho1(k-1)*dz1(k)+rho1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))
       dzi(k)   = (dz1(k-1)*dz1(k)+dz1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))
       exneri(k)= (ex1(k-1)*dz1(k)+ex1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))
    enddo
   
    !dxsa is scale-adaptive factor governing the pressure-gradient term of the momentum transport
    dxsa = one - MIN(MAX((12000.0_kind_phys-dx)/(12000.0_kind_phys-3000.0_kind_phys), zero), one)

    ktop_plume = 0
    ! do integration  updraft
    do ip=1,NUP
       QCn = zero
       overshoot = 0 !int
       l  = minwidth + dl*real(ip-1,kind=kind_phys)    ! diameter of plume
       do k=kts+1,kte-2
          zagl = zw1(k) + p5*dz1(k)
          !Entrainment from Tian and Kuang (2016)
          !ENT(k,ip) = 0.35/(MIN(MAX(UPW(K-1,ip),0.75),1.9)*l)
          wmin   = p3 + l*0.0005_kind_phys
          !ENT(k,ip) = 0.33_kind_phys/(MIN(MAX(UPW(k,ip),wmin),one)*l)
          !ENT(k,ip) = (0.20*sqrt(qkebl))/(MIN(MAX(UPW(K-1,ip),wmin),one)*l)
          ! 0.34 ~ 1.62
          ! 0.33 ~ 1.57
          ! 0.30 ~ 1.43
          ! 0.28 ~ 1.33
          ! 0.26 ~ 1.24
          entfac = 0.21_kind_phys * min(1.57_kind_phys, max(1.30_kind_phys, sqrt(qkebl)))
          !entfac = 0.33_kind_phys
          !make entfac tend to original value (0.33) above the pblh:
          wt2    = min(one, max(zero, zagl - pblh)/500._kind_phys) !0 in pbl, 1 aloft
          entfac = entfac*(one-wt2) + wt2*0.33_kind_phys
          ENT(k,ip) = entfac/(MIN(MAX(UPW(K-1,ip),wmin),one)*l)
          
          !Entrainment from Negggers (2015, JAMES)
          !ENT(k,ip) = 0.02*l**-0.35 - 0.0009
          !ENT(k,ip) = 0.04*l**-0.50 - 0.0009   !more plume diversity
          !ENT(k,ip) = 0.04*l**-0.495 - 0.0009  !"neg1+"

          !Minimum background entrainment 
          ENT(k,ip) = max(ENT(k,ip),0.0003_kind_phys)
          !ENT(k,ip) = max(ENT(k,ip),0.05/zw1(k))  !not needed for Tian and Kuang

          !increase entrainment for plumes extending very high.
          IF(zw1(k) >= MIN(pblh+1500., 4000.))THEN
            ENT(k,ip)=ENT(k,ip) + (zw1(k)-MIN(pblh+1500.,4000.))*5.0E-6
          ENDIF

          !SPP
          ENT(k,ip) = ENT(k,ip) * (one - pattern_spp_pbl1(k))

          ENT(k,ip) = min(ENT(k,ip), 0.9_kind_phys/dz1(k))

          ! Define environment U & V at the model interface levels
          Uk     =(u1(k)*dz1(k+1)+u1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
          Ukm1   =(u1(k-1)*dz1(k)+u1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))
          Vk     =(v1(k)*dz1(k+1)+v1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
          Vkm1   =(v1(k-1)*dz1(k)+v1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))

          ! Linear entrainment:
          EntExp =ENT(k,ip)*dz1(k)
          EntExm =EntExp*p333    !reduce entrainment for momentum
          QTn    =UPQT(k,IP)   *(one-EntExp) + qt1(k)*EntExp
          THLn   =UPTHL(k,IP)  *(one-EntExp) + thl1(k)*EntExp
          Un     =UPU(k,IP)    *(one-EntExm) + u1(k)*EntExm + dxsa*pgfac*(Uk - Ukm1)
          Vn     =UPV(k,IP)    *(one-EntExm) + v1(k)*EntExm + dxsa*pgfac*(Vk - Vkm1)
          QKEn   =UPQKE(k,IP)  *(one-EntExp) + qke1(k)*EntExp
          QSQn   =UPQSQ(k,IP)  *(one-EntExp) + qsq1(k)*EntExp
          QNCn   =UPQNC(k,IP)  *(one-EntExp) + qnc1(k)*EntExp
          QNIn   =UPQNI(k,IP)  *(one-EntExp) + qni1(k)*EntExp
          QNWFAn =UPQNWFA(k,IP)*(one-EntExp) + qnwfa1(k)*EntExp
          QNIFAn =UPQNIFA(k,IP)*(one-EntExp) + qnifa1(k)*EntExp
          QNBCAn =UPQNBCA(k,IP)*(one-EntExp) + qnbca1(k)*EntExp

          !capture the updated qc, qt & thl modified by entranment alone,
          !since they will be modified later if condensation occurs.
          qc_ent =QCn
          qt_ent =QTn
          thl_ent=THLn

          ! Exponential Entrainment:
          !EntExp= exp(-ENT(K,IP)*(zw1(k)-zw1(k-1)))
          !QTn =qt1(K) *(1-EntExp)+UPQT(K-1,IP)*EntExp
          !THLn=thl1(K)*(1-EntExp)+UPTHL(K-1,IP)*EntExp
          !Un  =u1(K)  *(1-EntExp)+UPU(K-1,IP)*EntExp
          !Vn  =v1(K)  *(1-EntExp)+UPV(K-1,IP)*EntExp
          !QKEn=qke1(k)*(1-EntExp)+UPQKE(K-1,ip)*EntExp

          if ( mix_chem ) then
            do ic = 1,nchem
              ! Exponential Entrainment:
              !chemn(ic) = chem(k,ic)*(1-EntExp)+UPCHEM(K-1,ip,ic)*EntExp
              ! Linear entrainment:
              chemn(ic)=UPCHEM(k,ip,ic)*(one-EntExp) + chem1(k,ic)*EntExp
            enddo
          endif
          if ( scalar_opt > 0 ) then
            do ic = 1,nscalars
              scalarsn(ic)=upscalars(k,ip,ic)*(one-EntExp) + scalars(k,ic)*EntExp
            enddo
          endif
          
          ! Define pressure at model interface we just integrated to (k+1).
          Pk    =(pres1(k)*dz1(k+1)+pres1(k+1)*dz1(k))/(dz1(k)+dz1(k+1))
          ! Compute plume properties thvn and qcn
          call condensation_edmf(QTn,THLn,Pk,zw1(k+1),THVn,QCn)

          ! Define environment THV at the model interface levels k & k+1
          THVk  =(thv1(k)*dz1(k+1)+thv1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
          THVkm1=(thv1(k-1)*dz1(k)+thv1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))

!          B=g*(0.5*(THVn+UPTHV(k-1,ip))/thv1(k-1) - one)
          B=grav*(THVn/THVk - one)
          IF(B>0.)THEN
            BCOEFF = 0.15        !w typically stays < 2.5, so doesnt hit the limits nearly as much
          ELSE
            BCOEFF = 0.2 !0.33
          ENDIF

          ! Original StEM with exponential entrainment
          !EntW=exp(-2.*(Wb+Wc*ENT(K,ip))*(zw1(k)-zw1(k-1)))
          !Wn2=UPW(K-1,ip)**2*EntW + (1.-EntW)*0.5*Wa*B/(Wb+Wc*ENT(K,ip))
          ! Original StEM with linear entrainment
          !Wn2=UPW(K-1,ip)**2*(1.-EntExp) + EntExp*0.5*Wa*B/(Wb+Wc*ENT(K,ip))
          !Wn2=MAX(Wn2,zero)
          !WA: TEMF form
          !IF (B>zero .AND. UPW(K-1,ip) < 0.2 ) THEN
          IF (UPW(K-1,ip) < 0.2 ) THEN
             Wn = UPW(k,ip) + (-2. * ENT(k,ip) * UPW(k,ip) + BCOEFF*B / MAX(UPW(k,ip),0.2)) * MIN(dz1(k), 250.)
          ELSE
             Wn = UPW(k,ip) + (-2. * ENT(k,ip) * UPW(k,ip) + BCOEFF*B / UPW(k,ip)) * MIN(dz1(k), 250.)
          ENDIF
          !Do not allow a parcel to accelerate more than 1.25 m/s over 200 m.
          !Add max increase of 2.0 m/s for coarse vertical resolution.
          IF(Wn > UPW(k,ip) + MIN(1.25*dz1(k)/200., two) ) THEN
             Wn = UPW(k,ip) + MIN(1.25*dz1(k)/200., two)
          ENDIF
          !Add symmetrical max decrease in w
          IF(Wn < UPW(k,ip) - MIN(1.25*dz1(k)/200., two) ) THEN
             Wn = UPW(k,ip) - MIN(1.25*dz1(k)/200., two)
          ENDIF
          Wn = MIN(MAX(Wn, zero), 3.0_kind_phys)

          !Check to make sure that the plume made it up at least one level.
          !if it failed, then set nup2=0 and exit the mass-flux portion.
          if (k==kts+1 .and. wn < 1e-8) then
             nup2=0
             exit
          endif

          if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
            if (wn .ge. 3.0) then
              ! surface values
              print *," **** SUSPICIOUSLY LARGE W:"
              print *,' QCn:',QCn,' ENT=',ENT(k,ip),' Nup2=',Nup2
              print *,'pblh:',pblh,' Wn:',Wn,' UPW(k)=',UPW(k,ip)
              print *,'K=',k,' B=',B,' dz=',dz1(k)
            ENDIF
          ENDIF

          !Allow strongly forced plumes to overshoot if KE is sufficient
          IF (Wn <= zero .AND. overshoot == 0) THEN
             overshoot = 1
             IF ( THVk-THVkm1 .GT. zero ) THEN
                bvf = SQRT( gtr*(THVk-THVkm1)/dz1(k) )
                !vertical Froude number
                Frz = UPW(k,ip)/(bvf*dz1(k))
                !IF ( Frz >= 0.5 ) Wn =  MIN(Frz,one)*UPW(K-1,I)
                dzp = dz1(k)*MAX(MIN(Frz,one),zero) ! portion of highest layer the plume penetrates
             ENDIF
          ELSE
             dzp = dz1(k)
          ENDIF

          !minimize the plume penetratration in stratocu-topped PBL
          !IF (fltv2 < 0.06) THEN
          !   IF(ZW(k+1) >= pblh-200. .AND. qc1(k) > 1e-5 .AND. ip > 4) Wn=zero
          !ENDIF

          !Modify environment variables (representative of the model layer - envm*)
          !following the updraft dynamical detrainment of Asai and Kasahara (1967, JAS).
          !Reminder: w is limited to be non-negative (above)
          aratio   = MIN(UPA(k,IP)/(one-UPA(k,IP)), p5) !limit should never get hit
          detturb  = 0.00008
          oow      = -0.060/MAX(one,(0.5*(Wn+UPW(k,IP))))   !coef for dynamical detrainment rate
          detrate  = MIN(MAX(oow*(Wn-UPW(k,IP))/dz1(k), detturb), .0002) ! dynamical detrainment rate (m^-1)
          detrateUV= MIN(MAX(oow*(Wn-UPW(k,IP))/dz1(k), detturb), .0001) ! dynamical detrainment rate (m^-1) 
          envm_thl(k)=envm_thl(k) + (p5*(thl_ent + UPTHL(k,IP)) - thl1(k))*detrate*aratio*MIN(dzp,dzpmax)
          qv_ent = p5*(MAX(qt_ent-qc_ent,zero) + MAX(UPQT(k,IP)-UPQC(k,IP),zero))
          envm_sqv(k)=envm_sqv(k) + (qv_ent-qv1(K))*detrate*aratio*MIN(dzp,dzpmax)
          IF (UPQC(k,IP) > 1E-8) THEN
             IF (qc1(k) > 1E-6) THEN
                qc_grid = qc1(k)
             ELSE
                qc_grid = qc_bl1(K)
             ENDIF
             envm_sqc(k)=envm_sqc(k) + MAX(UPA(k,IP)*p5*(QCn + UPQC(k,IP)) - qc_grid, zero)*detrate*aratio*MIN(dzp,dzpmax)
          ENDIF
          envm_u(k)  =envm_u(k)   + (p5*(Un + UPU(k,IP)) - u1(K))*detrateUV*aratio*MIN(dzp,dzpmax)
          envm_v(k)  =envm_v(k)   + (p5*(Vn + UPV(k,IP)) - v1(K))*detrateUV*aratio*MIN(dzp,dzpmax)

          IF (Wn > 0.) THEN
             !Update plume variables at the k+1 index we just integrated to.
             UPW(k+1,IP)=Wn  !sqrt(Wn2)
             UPTHV(k+1,IP)=THVn
             UPTHL(k+1,IP)=THLn
             UPQT(k+1,IP)=QTn
             UPQC(k+1,IP)=QCn
             UPU(k+1,IP)=Un
             UPV(k+1,IP)=Vn
             UPQKE(k+1,IP)=QKEn
             UPQSQ(k+1,IP)=QSQn
             UPQNC(k+1,IP)=QNCn
             UPQNI(k+1,IP)=QNIn
             UPQNWFA(k+1,IP)=QNWFAn
             UPQNIFA(K+1,IP)=QNIFAn
             UPQNBCA(k+1,IP)=QNBCAn
             UPA(k+1,IP)=UPA(K,IP)
             if ( mix_chem ) then
               do ic = 1,nchem
                 UPCHEM(k+1,ip,ic) = chemn(ic)
               enddo
             endif
             if ( scalar_opt > 0 ) then
               do ic = 1,nscalars
                 upscalars(k+1,ip,ic) = scalarsn(ic)
               enddo
             endif
             ktop = MAX(ktop,k)
          ELSE
             if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
                print*,"plume #:",ip," ktop=",ktop
                print*,"area=",UPA(1:ktop,ip)
                print*,"w=",UPW(1:ktop,ip)
             endif
             exit  !exit k-loop
          END IF
       ENDDO !end k-loop

       IF (debug_mf == 1 .and. i==idbg .and. j==jdbg) THEN
          IF (MAXVAL(UPW(:,IP)) > 10.0 .OR. MINVAL(UPA(:,IP)) < zero .OR. &
              MAXVAL(UPA(:,IP)) > Atot .OR. NUP2 > 10) THEN
             ! surface values
             print *,'flq:',flq,' fltv:',fltv2,' Nup2=',Nup2
             print *,'pblh:',pblh,' wstar:',wstar,' ktop=',ktop
             print *,'sigmaW=',sigmaW,' sigmaTH=',sigmaTH,' sigmaQT=',sigmaQT
             ! means
             print *,'u:',u1
             print *,'v:',v1
             print *,'thl:',thl1
             print *,'UPA:',UPA(:,IP)
             print *,'UPW:',UPW(:,IP)
             print *,'UPTHL:',UPTHL(:,IP)
             print *,'UPQT:',UPQT(:,IP)
             print *,'ENT:',ENT(:,ip)
          ENDIF
       ENDIF
       ktop_plume(ip)=k !index where each individual plume stopped rising.
       excess_h =max(excess_h,exc_heat)
       excess_q =max(excess_q,exc_moist)
    ENDDO !end plume # loop
 ELSE
    !At least one of the conditions was not met for activating the MF scheme.
    NUP2=0.
 END IF !end criteria check for mass-flux scheme

 ktop=MIN(ktop,KTE-1)
 IF (ktop == 0) THEN
    ztop = zero
 ELSE
    ztop = zw1(ktop)
 ENDIF

 IF (nup2 > 0) THEN
    !Calculate the fluxes for each variable
    !All s_aw* variable are == 0 at k=1
    DO ip=1,NUP
       DO k=kts,kte-1
          s_aw1(k)   = s_aw1(k)    + rhoi(k)*UPA(K,ip)*UPW(K,ip)*Psig_w
          s_awthl1(k)= s_awthl1(k) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPTHL(K,ip)*Psig_w
          s_awqt1(k) = s_awqt1(k)  + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQT(K,ip)*Psig_w
          !to conform to grid mean properties, move qc to qv in grid mean
          !saturated layers, so total water fluxes are preserved but 
          !negative qc fluxes in unsaturated layers is reduced.
!         if (qc1(k) > 1e-12 .or. qc1(k+1) > 1e-12) then
             qc_plume = UPQC(K,ip)
!         else
!            qc_plume = zero
!         endif
          s_awqc1(k) = s_awqc1(k)  + rhoi(k)*UPA(K,ip)*UPW(K,ip)*qc_plume*Psig_w
          s_awqv1(k) = s_awqt1(k)  - s_awqc1(k)
       ENDDO
    ENDDO
    !momentum
    if (momentum_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awu1(k) = s_awu1(k)   + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPU(K,ip)*Psig_w
             s_awv1(k) = s_awv1(k)   + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPV(K,ip)*Psig_w
          enddo
       enddo
    endif
    !tke
    if (tke_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awqke1(k)= s_awqke1(k) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQKE(K,ip)*Psig_w
          enddo
       enddo
       !qsq
       if (bl_mynn_closure > 2.5) then
          do ip=1,nup
             do k=kts,kte-1
                s_awqsq1(k)= s_awqsq1(k) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQSQ(K,ip)*Psig_w
             enddo
          enddo
       endif
    endif
    !chem
    if ( mix_chem ) then
       do k=kts,kte-1
          do ip=1,nup
             do ic = 1,nchem
                s_awchem1(k,ic) = s_awchem1(k,ic) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPCHEM(K,ip,ic)*Psig_w
             enddo
          enddo
       enddo
    endif
    !scalars
    if ( scalar_opt > 0 ) then
       do k=kts,kte-1
          do ip=1,nup
             do ic = 1,nscalars
                s_awscalars1(k,ic) = s_awscalars1(k,ic) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*upscalars(k,ip,ic)*Psig_w
             enddo
          enddo
       enddo
    endif
    if (numcon_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awqnc1(k)  = s_awqnc1(k)   + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQNC(K,ip)*Psig_w
             s_awqni1(k)  = s_awqni1(k)   + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQNI(K,ip)*Psig_w
          enddo
       enddo
    endif
    if (aerosol_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awqnwfa1(k)= s_awqnwfa1(k) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQNWFA(K,ip)*Psig_w
             s_awqnifa1(k)= s_awqnifa1(k) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQNIFA(K,ip)*Psig_w
             s_awqnbca1(k)= s_awqnbca1(k) + rhoi(k)*UPA(K,ip)*UPW(K,ip)*UPQNBCA(K,ip)*Psig_w
          enddo
       enddo
    endif

   !Flux limiter: Check ratio of heat flux at top of first model layer
   !and at the surface. Make sure estimated flux out of the top of the
   !layer is < fluxportion*surface_heat_flux
   IF (abs(s_aw1(kts+1)) > 1e-8_kind_phys) THEN
      flx1 = max(s_aw1(kts+1)*(thv1(kts)-thv1(kts+1))/dzi(kts+1),1.0e-6_kind_phys)
   ELSE
      flx1 = zero
      !print*,"ERROR: s_aw(kts+1) == 0, NUP=",NUP," NUP2=",NUP2,&
      !       " superadiabatic=",superadiabatic," KTOP=",KTOP
   ENDIF
   adjustment=one
   flt2=max(fltv,zero)
   !Print*,"Flux limiter in MYNN-EDMF, adjustment=",fluxportion*flt/dz(kts)/flx1
   !Print*,"flt/dz=",flt/dz(kts)," flx1=",flx1," s_aw(kts+1)=",s_aw(kts+1)
   IF (flx1 > fluxportion*flt2/dz1(kts) .AND. flx1>zero) THEN
      adjustment = max(0.01, fluxportion*flt2/dz1(kts)/flx1)
      s_aw1      = s_aw1*adjustment
      s_awthl1   = s_awthl1*adjustment
      s_awqt1    = s_awqt1*adjustment
      s_awqc1    = s_awqc1*adjustment
      s_awqv1    = s_awqv1*adjustment
      if (numcon_opt > 0) then
         s_awqnc1   = s_awqnc1*adjustment
         s_awqni1   = s_awqni1*adjustment
      endif
      if (aerosol_opt > 0) then
         s_awqnwfa1 = s_awqnwfa1*adjustment
         s_awqnifa1 = s_awqnifa1*adjustment
         s_awqnbca1 = s_awqnbca1*adjustment
      endif
      IF (momentum_opt > 0) THEN
         s_awu1  = s_awu1*adjustment
         s_awv1  = s_awv1*adjustment
      ENDIF
      IF (tke_opt > 0) THEN
         s_awqke1= s_awqke1*adjustment
         if (bl_mynn_closure>2.5) s_awqsq1= s_awqsq1*adjustment
      ENDIF
      IF ( mix_chem ) THEN
         s_awchem1 = s_awchem1*adjustment
      ENDIF
      IF ( scalar_opt > 0 ) THEN
         s_awscalars1 = s_awscalars1*adjustment
      ENDIF
      UPA = UPA*adjustment
   ENDIF
   if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
      print*,"Flux adjustment:"
      print*,"adjustment=",adjustment," fluxportion=",fluxportion
      print*,"flt2=",flt2," flx1=",flx1
   endif
      
   !Calculate mean updraft properties for output:
   !all edmf_* variables are interpolated plume quantities to mass levels
   k=kts !at first mass level, can not interpolate--use first interface values
   do ip=1,nup
      edmf_a1(k)  =edmf_a1(k)  +upa(k+1,ip)
      edmf_w1(k)  =edmf_w1(k)  +upa(k+1,ip)*p5*UPW(k+1,ip)
      edmf_qt1(k) =edmf_qt1(k) +upa(k+1,ip)*UPQT(k+1,ip)
      edmf_thl1(k)=edmf_thl1(k)+upa(k+1,ip)*UPTHL(k+1,ip)
      edmf_ent1(k)=edmf_ent1(k)+upa(k+1,ip)*ENT(k+1,ip)
      edmf_qc1(k) =edmf_qc1(k) +upa(k+1,ip)*UPQC(k+1,ip)
   enddo
   do ip=1,nup
      do k=kts+1,ktop_plume(ip)-1  !within the plume, we can interpolate:
         upak        =(upa(k+1,ip)*dzi(k) + upa(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
         edmf_a1(k)  =edmf_a1(k)  +upak
         edmf_w1(k)  =edmf_w1(k)  +upak*(upw(k+1,ip)*dzi(k)  + upw(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
         edmf_qt1(k) =edmf_qt1(k) +upak*(upqt(k+1,ip)*dzi(k) + upqt(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
         edmf_thl1(k)=edmf_thl1(k)+upak*(upthl(k+1,ip)*dzi(k)+ upthl(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
         edmf_ent1(k)=edmf_ent1(k)+upak*(ent(k+1,ip)*dzi(k)  + ent(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
         if (abs(upqc(k,ip)) < 1e-10_kind_phys .and. upqc(k+1,ip) > zero) then
            !mass at cloud base is kept equal to the interface value
            edmf_qc1(k) =edmf_qc1(k) +upak*upqc(k+1,ip)
         else
            edmf_qc1(k) =edmf_qc1(k) +upak*(upqc(k+1,ip)*dzi(k) + upqc(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
         endif
      enddo
   enddo
   !Now use a single level at the top of the plume so there is no averaging with zeros above:
   do ip=1,nup
      k=ktop_plume(ip)
      edmf_a1(k)  =edmf_a1(k)  +upa(k,ip)
      edmf_w1(k)  =edmf_w1(k)  +upa(k,ip)*upw(k,ip)
      edmf_qt1(k) =edmf_qt1(k) +upa(k,ip)*upqt(k,ip)
      edmf_thl1(k)=edmf_thl1(k)+upa(k,ip)*upthl(k,ip)
      edmf_ent1(k)=edmf_ent1(k)+upa(k,ip)*ent(k,ip)
      edmf_qc1(k) =edmf_qc1(k) +upa(k,ip)*upqc(k,ip)
   enddo
   do k=kts,kte-1
      !Note that only edmf_a1 is multiplied by Psig_w. This helps with the scale/regime-awareness.
      if (edmf_a1(k)>0.) then
         edmf_w1(k)  =edmf_w1(k)/edmf_a1(k)
         edmf_qt1(k) =edmf_qt1(k)/edmf_a1(k)
         edmf_thl1(k)=edmf_thl1(k)/edmf_a1(k)
         edmf_ent1(k)=edmf_ent1(k)/edmf_a1(k)
         edmf_qc1(k) =edmf_qc1(k)/edmf_a1(k)
         edmf_a1(k)  =edmf_a1(k)*Psig_w
      endif
   enddo ! end k                                                                                                                                                                                     
   do k=kts,kte-1
      !FIND MAXIMUM MASS-FLUX IN THE COLUMN:
      if(edmf_a1(k)*edmf_w1(k) > maxmf) maxmf = edmf_a1(k)*edmf_w1(k)
      !tke production by udrafts
      !instead of dTKE/dt = 1/2 w^3, multiply by 2 for QKE.
      tkeprod_up(k)=(abs(edmf_w1(k))**3)*edmf_a1(k)/(b1*max(p5*(el1(k)+el1(k+1)),0.1)) 
   enddo

   !smoke/chem
   if ( mix_chem ) then
      do k=kts,kte-1
        do ip=1,nup
          do ic = 1,nchem
            upak 	    =(upa(k+1,ip)*dzi(k) + upa(k,ip)*dzi(k+1))/(dzi(k+1)+dzi(k))
            edmf_chem(k,ic) = edmf_chem(k,ic) + upak*(upchem(k+1,ip,ic)*dzi(k) + upchem(k,ip,ic)*dzi(k+1))/(dzi(k+1)+dzi(k))
          enddo
        enddo
      enddo
      do k=kts,kte-1
         if (edmf_a1(k)>0.) then
            do ic = 1,nchem
               edmf_chem(k,ic) = edmf_chem(k,ic)/edmf_a1(k)
            enddo
         endif
      enddo
   endif
   
   !Calculate the effects environmental subsidence.
   !All envi_*variables are valid at the interfaces, like the edmf_* variables
   IF (env_subs) THEN
      DO k=kts+1,kte-1
         !First, smooth the profiles of w & a, since sharp vertical gradients
         !in plume variables are not likely extended to env variables
         !Note1: w is treated as negative further below
         !Note2: both w & a will be transformed into env variables further below
         envi_w(k) = p333*(edmf_w1(k-1)+edmf_w1(k)+edmf_w1(k+1))
         envi_a(k) = p333*(edmf_a1(k-1)+edmf_a1(k)+edmf_a1(k+1))*adjustment
      ENDDO
      !define env variables at k=1 (top of first model layer)
      envi_w(kts) = edmf_w1(kts)
      envi_a(kts) = edmf_a1(kts)
      !define env variables at k=kte
      envi_w(kte) = zero
      envi_a(kte) = edmf_a1(kte)
      !define env variables at k=kte+1
      envi_w(kte+1) = zero
      envi_a(kte+1) = edmf_a1(kte)
      !Add limiter for very long time steps (i.e. dt > 300 s)
      !Note that this is not a robust check - only for violations in
      !   the first model level.
      IF (envi_w(kts) > 0.9*dz1(kts)/dt) THEN
         sublim = 0.9*dz1(kts)/dt/envi_w(kts)
      ELSE
         sublim = one
      ENDIF
      !Transform w & a into env variables
      DO k=kts,kte
         temp=envi_a(k)
         envi_a(k)=one-temp
         envi_w(k)=csub*sublim*envi_w(k)*temp/(1.-temp)
      ENDDO
      !calculate tendencies from subsidence and detrainment valid at the middle of
      !each model layer. The lowest model layer uses an assumes w=0 at the surface.
      sub_thl1(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*thl1(kts+1)-rho1(kts)*thl1(kts))/dzi(kts)/rhoi(kts)
      sub_sqv1(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*qv1(kts+1)-rho1(kts)*qv1(kts))/dzi(kts)/rhoi(kts)
      DO k=kts+1,kte-1
         dzi(k)     = 0.5*(dz1(k)+dz1(k+1))
         sub_thl1(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho1(k+1)*thl1(k+1)-rho1(k)*thl1(k))/dzi(k)/rhoi(k)
         sub_sqv1(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho1(k+1)*qv1(k+1)-rho1(k)*qv1(k))/dzi(k)/rhoi(k)
      ENDDO

      DO k=kts,KTE-1
         det_thl1(k)=Cdet*(envm_thl(k)-thl1(k))*envi_a(k)*Psig_w
         det_sqv1(k)=Cdet*(envm_sqv(k)-qv1(k))*envi_a(k)*Psig_w
         det_sqc1(k)=Cdet*(envm_sqc(k)-qc1(k))*envi_a(k)*Psig_w
      ENDDO

      IF (momentum_opt > 0) THEN
         sub_u1(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*u1(kts+1)-rho1(kts)*u1(kts))/dzi(kts)/rhoi(kts)
         sub_v1(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*v1(kts+1)-rho1(kts)*v1(kts))/dzi(kts)/rhoi(kts)
         DO k=kts+1,kte-1
            sub_u1(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                       (rho1(k+1)*u1(k+1)-rho1(k)*u1(k))/dzi(k)/rhoi(k)
            sub_v1(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                       (rho1(k+1)*v1(k+1)-rho1(k)*v1(k))/dzi(k)/rhoi(k)
         ENDDO

         DO k=kts,KTE-1
           det_u1(k) = Cdet*(envm_u(k)-u1(k))*envi_a(k)*Psig_w
           det_v1(k) = Cdet*(envm_v(k)-v1(k))*envi_a(k)*Psig_w
         ENDDO
       ENDIF
   ENDIF !end subsidence/env detranment

   !First, compute plume theta & total water at interfaces (qti)
   !These values do not need to be defined at k=kte (unused level).
   DO k=kts,kte-1
      edmf_th1(k)= edmf_thl1(k) + xlvcp/ex1(k)*edmf_qc1(K)
      qti(k+1)   =(qt1(k)*dz1(k+1)+qt1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
   ENDDO
   qti(kts)=qt1(kts)
   
!update CLDFRA_bl1, qc_bl1. They have already been defined in
!     mym_condensation. Here, a shallow-cu component is added, but no cumulus
!     clouds can be added at k=1 (start loop at k=2).
   cb_check = 0
   do k=kts+1,kte-2
      if (k > KTOP) exit
      if (edmf_qc1(k) > zero ) then !.and. (cldfra_bl1(k) < cf_thresh))THEN
         if (cb_check == 0) then
            mf_at_cb = edmf_a1(k)*edmf_w1(k)
            cb_check = 1
          endif
          !plume properties within mass layers.
          Aup = edmf_a1(k)
          THp = edmf_th1(k)
          QTp = edmf_qt1(k)
          QCp = edmf_qc1(k)
          dqwdz = (qti(k)-qti(k+1))/dz1(k) !pos = decreasing upward
          !convert TH to T
          !t = THp*exner(k)
          !SATURATED VAPOR PRESSURE
          esat = esat_blend(tk1(k))
          !SATURATED SPECIFIC HUMIDITY
          qsl=ep_2*esat/max(1.e-7_kind_phys,(pres1(k)-ep_3*esat)) 

          !COMPUTE CLDFRA & QC_BL FROM MASS-FLUX SCHEME and recompute vt & vq
          xl = xl_blend(tk1(k))               ! obtain blended heat capacity
          qsat_tk = qsat_blend(tk1(k),pres1(k))! get saturation water vapor mixing ratio
                                              !   at t and p
          rsl = xl*qsat_tk / (r_v*tk1(k)**2)  ! slope of C-C curve at t (abs temp)
                                              ! CB02, Eqn. 4
          cpm = cp + qt1(k)*cpv               ! CB02, sec. 2, para. 1
          a   = one/(one + xl*rsl/cpm)        ! CB02 variable "a"
          b9  = a*rsl                         ! CB02 variable "b" 

          q2p = xlvcp/ex1(k)
          pt  = thl1(k) +q2p*QCp*Aup ! potential temp (env + plume)
          bb  = b9*tk1(k)/pt ! bb is "b9" in BCMT95.  Their "b9" differs from
                             ! "b9" in CB02 by a factor
                             ! of T/theta.  Strictly, b9 above is formulated in
                             ! terms of sat. mixing ratio, but bb in BCMT95 is
                             ! cast in terms of sat. specific humidity.  The
                             ! conversion is neglected here.
          qww   = one+p608*qt1(k)
          alpha = p608*pt
          beta  = pt*xl/(tk1(k)*cp) - 1.61_kind_phys*pt
          !Buoyancy flux terms have been moved to the end of this section...

          !Now calculate convective component of the cloud fraction:
          if (a > zero) then
             f = MIN(one/a, four)      ! f is vertical profile scaling function (CB2005)
          else
             f = one
          endif

          !---CB form:
          !sigq = 3.5E-3 * Aup * edmf_w1(k) * f  ! convective component of sigma (CB2005)
          !sigq = SQRT(sigq**2 + sgm1(k)**2)     ! combined conv + stratus components

          !---Per S.DeRoode 2009?
          !sigq = nine * Aup * (QTp - qt1(k))

          !---Extended CB form, tauc is timescale similar to the eddy turnover timescale.
          tauc = 1800._kind_phys
          sigq = five * Aup * edmf_w1(k) * max(4e-3_kind_phys, tauc*max(zero, dqwdz))

          !constrain sigq wrt saturation:
          !sigq = max(sigq, qsat_tk*0.03_kind_phys)
          !constrain sigq wrt moisture excess in the updraft (deRoode)
          sigq = max(sigq, three * Aup * (QTp - qt1(k)))
          !sigq = SQRT(sigq**2 + sgm1(k)**2)     ! combined conv + stratus components
          !sigq = max(sigq, sgm1(k))             ! use max of conv + stratus components

          !qmq = a * (qt1(k) - qsat_tk)          ! saturation deficit/excess;
          qmq = qt1(k) - qsat_tk                ! saturation deficit/excess;
          if (qmq > zero) sigq = min(qsat_tk*0.015_kind_phys, sigq)
          Q1  = qmq/sigq                        !   the numerator of Q1

          if ((landsea-1.5).GE.zero) then   ! WATER
             !modified form from LES
             !cf_mf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.2)),0.01),0.6)
             !Original CB
             !cf_mf = min(max(p5 + 0.36_kind_phys * atan(1.55*Q1),0.01_kind_phys),0.8_kind_phys)
             !increase cf in cases of large mf at cloudbase
             cffac = min(max(zero, mf_at_cb - 0.05_kind_phys)/0.05_kind_phys, one)
             cf_mf = min(max(p5 + 0.36_kind_phys * atan(1.8_kind_phys*(Q1+p2+cffac)),0.01_kind_phys), one)
             cf_mf = max(cf_mf, 1.2_kind_phys * Aup)
             !cf_mf = max(cf_mf, 1.8_kind_phys * Aup)
             !cf_mf = min(cf_mf, 5.0 * Aup)
          else                              ! LAND
             !LES form
             !cf_mf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.4)),0.01),0.6)
             !Original CB
             !cf_mf = min(max(p5 + 0.36_kind_phys * atan(1.55*Q1),0.01_kind_phys),0.8_kind_phys
             cf_mf = min(max(p5 + 0.36_kind_phys * atan(1.8_kind_phys*(Q1+p2)),0.01_kind_phys), one)
             cf_mf = max(cf_mf, 1.8_kind_phys * Aup)
             !cf_mf = min(cf_mf, 5.0 * Aup)
          endif

          !if ( debug_mf == 1 .and. i==idbg .and. j==jdbg) then
          !   print*,"In MYNN-MF, macrophysics component"
          !   print*," env qt=",qt1(k)," qsat=",qsat_tk
          !   print*," k=",k," satdef=",QTp - qsat_tk," sgm=",sgm1(k)
          !   print*," sigq=",sigq," qmq=",qmq," tk=",tk1(k)
          !   print*," cf_mf=",cf_mf," cldfra_bl=",cldfra_bl1(k)," edmf_a1=",edmf_a1(k)
          !endif

          !Update cloud fractions and specific humidities in grid cells
          !where the mass-flux scheme is active. The specific humidities
          !are converted to grid means (not in-cloud quantities).
          qcfac = 1.20_kind_phys + 1.0_kind_phys*min(max(zero, mf_at_cb - 0.04_kind_phys)/0.05_kind_phys, one)
          if ((landsea-1.5).GE.zero) then  ! water
             !if (QCp * Aup > 5e-5) then
             !   qc_mf     = 1.86_kind_phys * (QCp * Aup) - 2.2e-5_kind_phys
             !else
             !   qc_mf     = 1.20_kind_phys * (QCp * Aup)
             !endif
             qc_mf     = qcfac * (QCp * Aup) 
          else                             ! land
             if (QCp * Aup > 5e-5) then
                qc_mf     = 1.86_kind_phys * (QCp * Aup) - 2.2e-5_kind_phys
             else
                qc_mf     = 1.20_kind_phys * (QCp * Aup)
             endif
          endif
            !In the condition of very large cloud fractions, where the instantaneous
            !condensed water in the plume is not a reasonable estimate of the
            !mixing ratio in a larger cloud, we must rely on a background estimate based
            !off of environmental qsat.
            wt2          = one - min(one, max(zero, (cf_mf - thirty))/(hundred-thirty))  !=1 for low cf_mf, =0 for cf_mf = 1.0
            qc_mf_min    = wt2*qsat_tk*0.01_kind_phys*cf_mf + (one-wt2)*qsat_tk*0.025_kind_phys*cf_mf
            qc_mf        = max(qc_mf, qc_mf_min)
            !Then blend with the stratus component:
            cf_strat     = cldfra_bl1(k)
            qc_strat     = qc_bl1(k)
            cldfra_bl1(k)= max(cf_mf, cf_strat)
            pct_mf       = cf_mf/cldfra_bl1(k)
            qc_bl1(k)    = qc_mf*pct_mf + (one-pct_mf)*qc_strat

            !Now recalculate the terms for the buoyancy flux for mass-flux clouds:
            !See mym_condensation for details on these formulations.
            !Use Bechtold and Siebesma (1998) piecewise estimation of Fng with 
            !limits ,since they really should be recalculated after all the other changes...:
            !Only overwrite vt & vq in non-stratus condition
            !if ((landsea-1.5).GE.zero) then      ! WATER
               Q1=max(Q1,-2.25)
            !else
            !   Q1=max(Q1,-2.0)
            !endif

            if (Q1 .ge. one) then
               Fng = one
            elseif (Q1 .ge. -1.7 .and. Q1 .lt. one) then
               Fng = EXP(-p4*(Q1-one))
            elseif (Q1 .ge. -2.5 .and. Q1 .lt. -1.7) then
               Fng = three + EXP(-3.8*(Q1+1.7))
            else
               Fng = min(23.9 + EXP(-1.6*(Q1+2.5)), 60.)
            endif

            !link the buoyancy flux function to active clouds only (c*Aup):
            vt1(k) = qww   - (1.5*Aup)*beta*bb*Fng - one
            vq1(k) = alpha + (1.5*Aup)*beta*a*Fng  - tv0
         endif !check for (qc in plume) .and. (cldfra_bl < threshold)
      enddo !k-loop

ENDIF  !end nup2 > 0

!modify output (negative: dry plume, positive: moist plume)
if (ktop > 0) then
   maxqc = maxval(edmf_qc1(1:ktop)) 
   if ( maxqc < 1.E-8) maxmf = -1.*maxmf
endif

!
! debugging
!
if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then

! means
!   print *,'u:',u1
!   print *,'v:',v1
!   print *,'thl:',thl1
!   print *,'thv:',thv1
!   print *,'qt:',qt1
!   print *,'p:',pres1
 
! updrafts
! DO ip=1,NUP2
!   print *,'up:A',ip
!   print *,UPA(:,ip)
!   print *,'up:W',ip
!   print*,UPW(:,ip)
!   print *,'up:thv',ip
!   print *,UPTHV(:,ip)
!   print *,'up:thl',ip 
!   print *,UPTHL(:,ip)
!   print *,'up:qt',ip
!   print *,UPQT(:,ip)
!   print *,'up:tQC',ip
!   print *,UPQC(:,ip)
!   print *,'up:ent',ip
!   print *,ENT(:,ip)
! ENDDO

   print*,"mean updraft properties at end of mf component:"
   print*,' edmf_a1',edmf_a1(1:max(1,ktop))
   print*,' edmf_w1',edmf_w1(1:max(1,ktop))
   print*,' edmf_qt1:',edmf_qt1(1:max(1,ktop))
   print*,' edmf_thl1:',edmf_thl1(1:max(1,ktop))
 
ENDIF !END Debugging


#ifdef HARDCODE_VERTICAL
# undef kts
# undef kte
#endif

END SUBROUTINE DMP_MF
!=================================================================
!>\ingroup gsd_mynn_edmf
!! This subroutine 
subroutine condensation_edmf(QT,THL,P,zagl,THV,QC)
!
! zero or one condensation for edmf: calculates THV and QC
!
  use module_bl_mynnedmf_common, only: p1000mb,xlvcp,rcp,rvovrd,&
       zero,one,p5,kind_phys
  
real(kind_phys),intent(in)   :: QT,THL,P,zagl
real(kind_phys),intent(inout):: THV
real(kind_phys),intent(inout):: QC

integer :: niter,ni
real(kind_phys):: diff,exn,t,th,qs,qcold

! constants used from module_model_constants.F
! p1000mb
! rcp ... Rd/cp
! xlv ... latent heat for water (2.5e6)
! cp
! rvord .. r_v/r_d  (1.6) 

! number of iterations
  niter=50
! minimum difference (usually converges in < 8 iterations with diff = 2e-5)
  diff=1.e-6_kind_phys

  EXN=(P/p1000mb)**rcp
  QC=zero  !better first guess QC is incoming from lower level, do not set to zero
  do ni=1,NITER
     T=EXN*THL + xlvcp/EXN*QC
     QS=qsat_blend(T,P)
     QCOLD=QC
     QC=p5*QC + p5*MAX((QT-QS),zero)
     if (abs(QC-QCOLD)<Diff) exit
  enddo

  T=EXN*THL + xlvcp/EXN*QC
  QS=qsat_blend(T,P)
  QC=max(QT-QS,zero)

  !Do not allow saturation below 100 m
  if(zagl < 100.)QC=zero

  !THV=(THL+xlv/cp*QC).*(1+(1-rvovrd)*(QT-QC)-QC);
  THV=(THL+xlvcp/EXN*QC)*(one+QT*(rvovrd-one)-rvovrd*QC)

!  IF (QC > zero) THEN
!    PRINT*,"EDMF SAT, p:",p," iterations:",ni
!    PRINT*," T=",T," THL=",THL," THV=",THV
!    PRINT*," QS=",QS," QT=",QT," QC=",QC,"ratio=",qc/qs
!  ENDIF

  !THIS BASICALLY GIVES THE SAME RESULT AS THE PREVIOUS LINE
  !TH = THL + xlvcp/EXN*QC
  !THV= TH*(1. + p608*QT - QC)

end subroutine condensation_edmf

!===============================================================
subroutine condensation_ddmf(qt,thl,p,zagl,thv,qc,qi,debug_dd)
!
! zero or one condensation for edmf: calculates thv, qc, and qi
!
  use module_bl_mynnedmf_common, only: p1000mb,xlvcp,xlscp,rcp, &
       rvovrd,p608,zero,one,p5,kind_phys
  
real(kind_phys),intent(in)   :: qt,thl,p,zagl
real(kind_phys),intent(inout):: thv,qc,qi

integer :: niter,ni,debug_dd
real(kind_phys):: diff,exn,t,th,qs,qx,qxold,frac_ice,frac_liq,thvin

! number of iterations
  niter=50
! minimum difference (usually converges in < 8 iterations with diff = 2e-5)
  diff=1.e-6

  if ((qi+qc) .gt. zero) then
     frac_ice = qi/(qi+qc)
     frac_liq = one - frac_ice
  else
     frac_ice = zero
     frac_liq = zero
  endif

  if (debug_dd == 1) then
     thvin = thv
     print*,"------- in consensation_ddmf-------------"
     print*,"input qc=",qc," qi=",qi," frac_liq=",frac_liq
  endif

  exn=(p/p1000mb)**rcp
  do ni=1,niter
     t=exn*thl + xlvcp/exn*qc + xlscp/exn*qi
     qs=qsat_blend(t,p)
     qxold=qc+qi
     qx=qc+qi
     qx=p5*qx + p5*max((qt-qs),zero)
     qc=frac_liq*qx
     qi=frac_ice*qx
     if (abs(qx-qxold)<diff) exit
  enddo

  t=exn*thl + xlvcp/exn*qc + xlscp/exn*qi
  qs=qsat_blend(t,p)
  qx=max(qt-qs,zero)
  qc=frac_liq*qx
  qi=frac_ice*qx

  !thv=(thl+xlv/cp*qc).*(1+(1-rvovrd)*(qt-qc)-qc)
!was this: thv=(thl+xlvcp*qc)*(1.+qt*(rvovrd-1.)-rvovrd*qc)
  !thv=(thl + xlvcp*qc + xlscp*qi)*(1. + qt*(rvovrd-1.)-rvovrd*qc)
  thv=(thl + xlvcp/exn*qc + xlscp/exn*qi)*(one + p608*(qt-qx) - qx)

!  if (qc > zero) then
!    print*,"edmf sat, p:",p," iterations:",ni
!    print*," t=",t," thl=",thl," thv=",thv
!    print*," qs=",qs," qt=",qt," qc=",qc,"ratio=",qc/qs
!  endif

  !this basically give the same result as the previous line
  !th = thl + xlv/cp/exn*qc
  !thv= th*(1. + p608*qt)

  !print *,'t,p,qt,qs,qc'
  !print *,t,p,qt,qs,qc

  if (debug_dd == 1) then
     print*,"output qc=",qc," qi=",qi," frac_liq=",frac_liq
     print*,"input thv=",thvin," out thv=",thv," diff=",thv-thvin
     print*,"------- exiting consensation_ddmf----------"
  endif

end subroutine condensation_ddmf

!===============================================================

subroutine condensation_edmf_r(QT,THL,P,zagl,THV,QC)
!                                                                                                
! zero or one condensation for edmf: calculates THL and QC                                       
! similar to condensation_edmf but with different inputs                                         
!
  use module_bl_mynnedmf_common, only: p1000mb,rcp,rvovrd,xlvcp,&
       one,zero,kind_phys
  
real(kind_phys),intent(in)   :: QT,THV,P,zagl
real(kind_phys),intent(inout):: THL, QC

integer :: niter,ni
real(kind_phys):: diff,exn,t,th,qs,qcold

! number of iterations                                                                           
  niter=50
! minimum difference                                                                             
  diff=2.e-5_kind_phys

  EXN=(P/p1000mb)**rcp
  ! assume first that th = thv                                                                   
  T = THV*EXN
  !QS = qsat_blend(T,P)                                                                          
  !QC = QS - QT                                                                                  

  QC=zero

  do ni=1,NITER
     QCOLD = QC
     T = EXN*THV/(one+QT*(rvovrd-one)-rvovrd*QC)
     QS=qsat_blend(T,P)
     QC= MAX((QT-QS), zero)
     if (abs(QC-QCOLD)<Diff) exit
  enddo
  THL = (T - xlvcp*QC)/EXN

end subroutine condensation_edmf_r

!===============================================================
! ddmp_mf is a downdraft mass flux scheme - analogous to the updrafts but
! flipped upsidedown and revised to be more specific for turbulence
! driven by radiative cooling at cloud tops. this scheme has been primarily
! designed for stratocumulus cloud conditions. This is based off of a
! combination of the original scheme [Wu et al. (2020, MWR)] but several
! modifications were made to better fit it in the MYNN-EDMF and to
! help the coupling to the Thompson microphysics scheme. Some of the
! primary design changes include:
!                                                                                                                                                                      
! 1) the spectral plume design, similar to the MYNNs updraft scheme
! 2) generalized to any cloud, not just PBL-topped clouds
! 3) always initialized at the cloud top
! 4) can penetrate all the way to the surface
! 5) generalized for mixed-phase clouds
! 6) mixes number concentration for double moment schemes
! 7) mixes aerosols
! 8) revised the conection to the solver

subroutine ddmp_mf(kts,kte,dt,dx,zw,dz,p,            &
              &u,v,th,thl,thv,tk,qt,qv,qc,qi,        &
              &qnc,qni,qnwfa,qnifa,                  &
              &qke,rho,exner,                        &
              &qc_bl1,qi_bl1,cldfra_bl1,             &
              &ust,flt,flq,fltv,pblh,kpbl,           &
              &edmf_a_dd,edmf_w_dd, edmf_qt_dd,      &
              &edmf_thl_dd,edmf_ent_dd,edmf_qc_dd,   &
              &sd_aw,sd_awthl,sd_awqt,               &
              &sd_awqv,sd_awqc,sd_awqi,              &
              &sd_awqnc,sd_awqni,                    &
              &sd_awqnwfa,sd_awqnifa,                &
              &sd_awu,sd_awv,sd_awqke,               &
              &tkeprod_dn,el,                        &
              &rthraten,psig,                        &
              &maxmf_dd,maxwidth_dd                  )

  use module_bl_mynnedmf_common, only: cp,grav,p608, &
       xlvcp,b1,one,two,zero,eight,p1,p2,p333,p5,    &
       p666,kind_phys

        integer, intent(in) :: kts,kte,kpbl
        real(kind_phys), dimension(kts:kte), intent(in) ::            &
            u,v,th,thl,tk,qt,qv,qc,qi,thv,p,qke,rho,exner,            &
            qnc,qni,qnwfa,qnifa,dz,qc_bl1,qi_bl1,cldfra_bl1,el
        real(kind_phys), dimension(kts:kte), intent(in) :: rthraten
        ! zw .. heights of the downdraft levels (edges of boxes)
        real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
        real(kind_phys), intent(in)  :: flt,flq,fltv
        real(kind_phys), intent(in)  :: dt,dx,ust,pblh,psig
   !outputs - downdraft properties
        real(kind_phys), dimension(kts:kte), intent(inout) ::         &
        edmf_a_dd,   edmf_w_dd,   edmf_qt_dd,  edmf_thl_dd,           &
        edmf_ent_dd, edmf_qc_dd,  tkeprod_dn
   !outputs - variables needed for solver (sd_aw - sum ai*wi, sd_awphi - sum ai*wi*phii)
        real(kind_phys), dimension(kts:kte+1) ::                      &
            sd_aw, sd_awthl, sd_awu, sd_awv,                          &
            sd_awqt, sd_awqc, sd_awqv, sd_awqi, sd_awqnc, sd_awqni,   &
            sd_awqnwfa, sd_awqnifa, sd_awqke, sd_aw2
   !outputs - diagnostics
        real(kind_phys), intent(inout) :: maxmf_dd,maxwidth_dd

   !downdraft properties
        integer, parameter::                &
            & ndd            = 5,           &  !number of downdrafts
            & kmin           = 3               !lowest k-level where downdrafts can start
        real(kind_phys),parameter ::        &
            & minddd         = 100.,        &  !min downdraft diameter (m)
            & maxddd         = 500.,        &  !max downdraft diameter (m)
            & zmin           = 50.,         &  !lowest height where downdrafts can start
            & dz200          = 200.,        &  !depth over which other parameters are normalized to
            & cooling_thresh = -0.00010        !activation threshold of downdrafts
                              !-0.000116 is ~ -10 C cooling at cloud top per day
                              !-0.000093 is ~  -8 C cooling at cloud top per day
                              !-0.000069 is ~  -6 C cooling at cloud top per day
        real(kind_phys)::                   &
            & maxdd2,                       &  !variable max downdraft diameter (m)
            &    ddd,                       &  !downdraft diameter
            &     dl,                       &  !diameter increment
            &    add,                       &  !total area of downdrafts
            & minexc,                       &  !minimum init downdraft temp pert
            & maxexc                           !maximum init downdraft temp pert
  ! k-index of downdraft starting height
        integer,         dimension(1:ndd) :: dd_initk
  ! downdraft column properties
        real(kind_phys), dimension(kts:kte+1,1:ndd) ::                &
            downw,downthl,downqt,downqc,downa,downu,downv,downthv,    &
            downqi,downqnc,downqni,downqnwfa,downqnifa
  ! entrainment variables
        real(kind_phys), dimension(kts:kte+1,1:ndd) :: ent
  ! internal variables
        real(kind_phys), dimension(kts:kte):: massflux
        integer :: k,i,ki,qltop,qlbase
        real(kind_phys):: qstar,thstar,sigmaw,sigmaqt,                &
            sigmath,z0,pwmin,pwmax,wmin,wmax,went,mindownw
        real(kind_phys):: qtn,thln,thvn,qcn,un,vn,qken,wn2,wn,        &
            qin,qncn,qnin,qnwfan,qnifan,                              &
            thvk,pk,entexp,entw,beta_dm,entexm,rho_int
        real(kind_phys):: jump_thv,jump_qt,jump_thetal,               &
            refthl,refthv,refthlv,refqt,refqc,refqi,refqnc,refqni,    &
            refqnwfa,refqnifa,refu,refv,refqke,refqt2,reftk,refp,     &
            qx_k,qx_km1,refthvm1,cftop,crate,ac_wsp,wspd_pbl,         &
            maxcooling
        real(kind_phys):: dthvx,tmp1,rcldb,ent_eff

  ! dd specific internal variables
        real(kind_phys):: radflux, f0, wstar_rad, dz_ent
        logical :: cloudflg
        logical :: singlelayer             !check for single or multi-layer clouds
        real(kind_phys):: sigq,xl,rsl,cpm,a,diffqt,dp,                &
            fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid,         &
            def_th,def_qt,frac_liq,frac_ice,qs1,qs2

  ! w parameters
        real(kind_phys),parameter ::                                  &
            &wa=1., wb=1.5, z00=100.
        real(kind_phys) :: buoy,bcoeff,ecoeff

  ! additional printouts for debugging
        integer, parameter :: debug_dd=0   !0:no debugging, 1:detailed output, 2:output for strange values only

  ! =================================================================
  ! initialize downdraft properties
   downw      =zero
   downthl    =zero
   downthv    =zero
   downqt     =zero
   downa      =zero
   downu      =zero
   downv      =zero
   downqc     =zero
   downqi     =zero
   downqnc    =zero
   downqni    =zero
   downqnwfa  =zero
   downqnifa  =zero
   ent        =zero
   dd_initk   =0

   edmf_a_dd  =zero
   edmf_w_dd  =zero
   edmf_qt_dd =zero
   edmf_thl_dd=zero
   edmf_ent_dd=zero
   edmf_qc_dd =zero

   sd_aw      =zero
   sd_awthl   =zero
   sd_awqt    =zero
   sd_awqv    =zero
   sd_awqc    =zero
   sd_awqi    =zero
   sd_awqnc   =zero
   sd_awqni   =zero
   sd_awqnwfa =zero
   sd_awqnifa =zero
   sd_awu     =zero
   sd_awv     =zero
   sd_awqke   =zero

   ! first, check for cloud tops. if more than one, find the one with maximum cooling at the cloud top.
   cloudflg   =.false.
   singlelayer=.false.
   qltop      =1     !initialize
   qlbase     =1     !initialize
   crate      =zero
   do k = max(kmin,kpbl-5),min(kpbl+20,kte-1)
      maxcooling = minval(rthraten(k-1:k+1))
      !---------------------------------criteria for downdraft activation:
      if ((cldfra_bl1(k).gt.p5 .and. cldfra_bl1(k+1).lt.p5) .and. &  !1) stratocu cloud top exists
          (maxcooling  .lt. cooling_thresh)) then                    !2) significant radiative cooling
         ! found sc cloud with significant radiative cooling
         cloudflg =.true.
         if (maxcooling < crate) then
            qltop    =k               !index for sc cloud top
            cftop    =cldfra_bl1(k)
            crate    =maxcooling      !maximum cooling rate at cloud top
            if (cldfra_bl1(k-1).lt.p2)singlelayer=.true.
         endif
      endif
   enddo

   ! determine maximum downdraft width and increments
   maxdd2  = min(maxddd, 1.2_kind_phys*dx)                                       !1) grid-scale dependence
   maxdd2  = min(maxdd2, -100._kind_phys*crate*3600._kind_phys + 200._kind_phys) !2) function of cloud top cooling rate
   maxdd2  = min(maxdd2, 80._kind_phys + p666*zw(qltop+1))                       !3) limited for low cloud-top heights
   maxwidth_dd = maxdd2
   dl      = max(maxdd2-minddd, zero)/real(ndd-1)
   
   !found sc cloud with conditions for downdrafts; compute downdrafts
   if (cloudflg .and. ( dl .gt. zero)) then
      do k = qltop, kts, -1
         qx_k =max(qt(k), qc_bl1(k) +qi_bl1(k))     !total cloud water and ice at k
         if (qx_k .gt. 1e-6_kind_phys) then
            qlbase = k ! index for sc cloud base
         endif
      enddo

      do i=1,ndd
         ! downdraft starts at the first layer interface below the cloud top
         dd_initk(i) = qltop
      enddo

      ! loop radflux
      f0 = zero
      do k = max(kmin,qltop-2), qltop+1
         radflux = rthraten(k) * exner(k)     ! converts theta/s to temperature/s
         dp      = p5 * (( p(k) - p(k+1) ) + ( p(k-1) - p(k) ))
         radflux = radflux * cp / grav * dp   ! converts k/s to w/m^2
         if ( radflux < zero ) f0 = abs(radflux) + f0
      enddo
      f0 = min(max(f0, one), 200._kind_phys)  !total radiative cooling (w/m2)

      !allow the total fractional area of the downdrafts (add) to be proportional
      !to the radiative forcing:
      !for  50 w/m2, add = 0.10
      !for 100 w/m2, add = 0.20
      !for 150 w/m2, add = 0.30
      add      = min( f0*0.003, 0.3_kind_phys)
      !taper off area for high wind speeds
      wspd_pbl = sqrt(max(u(kts)**2 + v(kts)**2, 0.01_kind_phys))
      ac_wsp   = one - min(max(zero, wspd_pbl - 10._kind_phys)/15._kind_phys, one)
      add      = add * min(ac_wsp, psig)

      !find inversion strength across cloud top entrainment zone--normalized to 200 m vertical grid spacing
      dz_ent      = p5 * (dz(qltop+1) + dz(qltop))
      jump_thv    = (thv(qltop+1) - thv(qltop)) / dz_ent * dz200
      jump_qt     = (qt(qltop+1)  - qt(qltop))	/ dz_ent * dz200
      jump_thetal = (thl(qltop+1) - thl(qltop))	/ dz_ent * dz200

      ki = qltop  !index of cloud top
      if (singlelayer) then !initialize dd with cloud properties (not using info from above)
         refthl  = thl(ki)
         refthlv = refthl*(1.+p608*qv(ki))
         refthv  = thv(ki)
         refthvm1= thv(ki-1)
         reftk   = tk(ki)
         refqt   = qt(ki)
         refqc   = qc(ki)
         refqi   = qi(ki)
         refqnc  = qnc(ki)
         refqni  = qni(ki)
         refqnwfa= qnwfa(ki)
         refqnifa= qnifa(ki)
         refu    = u(ki)
         refv    = v(ki)
         refqke  = qke(ki)
         refp    = p(ki)
      else                  !initialize dd with avg of in-cloud properties (at & below) cloud top 
         refthl  = (thl(ki-1)*dz(ki)   + thl(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refthlv = refthl*(1.+p608*(qv(ki-1)*dz(ki) + qv(ki)*dz(ki-1))/(dz(ki)+dz(ki-1)))
         refthv  = (thv(ki-1)*dz(ki)   + thv(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refthvm1= (thv(ki-2)*dz(ki-1) + thv(ki-1)*dz(ki-2))/(dz(ki-1)+dz(ki-2)) 
         reftk   = (tk(ki-1)*dz(ki)    + tk(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
         refqt   = (qt(ki-1)*dz(ki)    + qt(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
         refqc   = (qc(ki-1)*dz(ki)    + qc(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
         refqi   = (qi(ki-1)*dz(ki)    + qi(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
         refqnc  = (qnc(ki-1)*dz(ki)   + qnc(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refqni  = (qni(ki-1)*dz(ki)   + qni(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refqnwfa= (qnwfa(ki-1)*dz(ki) + qnwfa(ki)*dz(ki-1))/(dz(ki)+dz(ki-1))
         refqnifa= (qnifa(ki-1)*dz(ki) + qnifa(ki)*dz(ki-1))/(dz(ki)+dz(ki-1))
         refu    = (u(ki-1)*dz(ki)     + u(ki)*dz(ki-1))    /(dz(ki)+dz(ki-1))
         refv    = (v(ki-1)*dz(ki)     + v(ki)*dz(ki-1))    /(dz(ki)+dz(ki-1))
         refqke  = (qke(ki-1)*dz(ki)   + qke(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refp    = (p(ki-1)*dz(ki)     + p(ki)*dz(ki-1))    /(dz(ki)+dz(ki-1))
      endif

      ! w* from radiative forcing (m/s)
      !wstar_rad =    ( grav * dz_ent * f0 / (refthl * rho(ki) * cp) ) **p333
      wstar_rad = 1.25_kind_phys * ( grav*dz200 * f0 / (refthl * rho(ki) * cp) ) **p333
      wstar_rad = min(max(wstar_rad, p1), 1.5_kind_phys)      ! (m/s)

      !entrainment efficiency
      dthvx     = (thl(ki+2) + th(ki+2)*p608*qt(ki+2)) - &
                  (thl(ki)   + th(ki)  *p608*qt(ki))
      dthvx     = dthvx * dz200/dz_ent  !normalized gradient
      dthvx     = max(dthvx, p1)
      rcldb     = max(refqc+refqi, qc_bl1(ki)+qi_bl1(ki))  !any type of cloud water (sgs or resolved) 
      tmp1      = xlvcp * rcldb/(exner(ki)*dthvx)
      !entrainment efficiency: originally from Nichols and Turton (1986), where a2 = 60, but lowered
      !here to 8, as in Grenier and Bretherton (2001).
      ent_eff   = max(zero, min(one, p2 + p2*eight*tmp1))

      went      = p333 * ent_eff * wstar_rad
      qstar     = went * jump_qt / wstar_rad
      thstar    = -f0/rho(ki)/cp/wstar_rad - went*jump_thv/wstar_rad

      sigmaw    = 0.2 * wstar_rad   ! 0.8*wstar_rad ! tuning parameter ! 0.5 was good
      sigmaqt   = 50. * abs(qstar)  ! 50 was good
      sigmath   = 1.5 * abs(thstar) ! 0.5 was good

      pwmin     = -1. ! drawing from the negative tail -3sigma to -1sigma
      pwmax     = -3.
      wmin      = max(min(sigmaw*pwmin, -0.1), -0.2)
      wmax      = max(min(sigmaw*pwmax, -0.2), -0.5)

      !initialize downdraft temperature deficit def_th
      minexc = min(-0.06     , 0.5*(refthvm1 - refthv)) !increase prob to go down 1 level
      maxexc = min(minexc-0.1, -0.6          )
      maxexc = max(maxexc    , thv(1) - refthlv)        !init parcel temp >= thv_sfc
      def_th = max(min(0.05*(-0.3)*sigmath/sigmaw, minexc), maxexc)

      !initialize downdraft moisture deficit def_qt (consistent with def_th)
      qs1    = qsat_blend(reftk       ,refp)
      qs2    = qsat_blend(reftk+def_th,refp)
      refqt2 = refqt*qs2/qs1
      !def_qt = min(refqt2 - refqt, zero) !small negative moisture perturbation
      def_qt = min(qs2 - qs1,  -0.33*(refqc+refqi)) !small negative moisture perturbation
      !def_qt = 0.5*(-0.3)*sigmaqt/sigmaw

      if (refqc+refqi > zero) then
         frac_liq = refqc/(refqc+refqi)
         frac_ice = one - frac_liq
      else
         frac_liq = zero
         frac_ice = zero
      endif

      if (debug_dd .eq. 1) then
         print*,"-------------------------------------"
         print*,"found conditions for downdraft mixing"
         print*,"qltop=",qltop," qlbase=",qlbase
         print*,"q*=",qstar," theta*=",thstar," went=",went
         print*,"f0=",f0," jump_thv=",jump_thv,"w*rad=",wstar_rad
         print*,"u*=",ust," jump_qt=",jump_qt," fltv=",fltv
         print*,"grav=",grav," pblh=",pblh," thv(1)=",thv(1)
         print*,"p(1)=",p(1)," jump_thetal=",jump_thetal
         print*,"sigmaw=",sigmaw," wmin=",wmin," wmax=",wmax
      endif

      !intialize downdrafts
      !as of now, all downdrafts start at the same height
      do i=1,ndd
         ki = dd_initk(i)

         !initialize a negative (downward) vertical velocity
         downw(ki,i)=wmin-(abs(wmax)-abs(wmin))/real(ndd-1)*(i-1)

         !perhaps multiply downa by cloud fraction, so it's impact will diminish if
         !clouds are mixed away over the course of the longer radiation time step(?)
         downa(ki,i)     = add/real(ndd)
         downu(ki,i)     = refu
         downv(ki,i)     = refv
         downqc(ki,i)    = max(refqc + def_qt*frac_liq, zero)
         downqi(ki,i)    = max(refqi + def_qt*frac_ice, zero)
         downqnc(ki,i)   = refqnc
         downqni(ki,i)   = refqni
         downqnwfa(ki,i) = refqnwfa
         downqnifa(ki,i) = refqnifa
         downqt(ki,i)    = refqt  + def_qt
         downthv(ki,i)   = refthv + def_th*real(i)/real(ndd)
         downthl(ki,i)   = refthl + def_th*real(i)/real(ndd)

         !input :: qt,thv,p,zagl,  output :: thl, qc
!         pk  =(p(ki-1)*dz(ki)+p(ki)*dz(ki-1))/(dz(ki)+dz(ki-1))
!         call condensation_edmf_r(downqt(ki,i),   &
!              &        downthl(ki,i),pk,zw(ki),   &
!              &     downthv(ki,i),downqc(ki,i)    )
      enddo

      if (debug_dd .eq. 1) then
         print*,"=====initialized downdraft properties===="
         print*,"downw=",downw(ki,:)
         print*,"downa=",downa(ki,:)
         print*,"downthv=",downthv(ki,:)
         print*,"def_qt=",def_qt," def_th=",def_th," tot area=",add
         print*,"begin integration of downdrafts (smallest to largest):"
      endif

      do i=1,ndd
         ddd  = minddd + dl*real(i-1) !downdraft diameter
         !initialize input into condensation routine
         thvn = downthv(ki,i)
         qcn  = downqc(ki,i)
         qin  = downqi(ki,i)
         do k=dd_initk(i)-1,kts+1,-1
            !entrainment from tian and kuang (2016), with constraints
            wmin = p1 + ddd*0.0005
            !ent(k,i) = 0.33/(min(max(abs(downw(k+1,i)),wmin),0.9)*ddd)
            !ent(k,i) = 0.53/(min(max(abs(downw(k+1,i)),wmin),0.6)*ddd)
            !was       0.26  before reducing downdraft sizes
            ent(k,i) = 0.18/(min(max(abs(downw(k+1,i)),wmin),0.6)*ddd)
            
            !minimum background entrainment and 1/z enhancement near surface
            ent(k,i) = max(ent(k,i),0.0003)
            ent(k,i) = max(ent(k,i),0.06/zw(k))
            ent(k,i) = min(ent(k,i),0.9/dz(k)) ! liberal check for numerical stability

            !starting at the first interface level below cloud top
            !entexp=exp(-ent(k,i)*dz(k))
            !entexp_m=exp(-ent(k,i)/3.*dz(k))
            entexp  =ent(k,i)*dz(k)        !for all scalars
            entexm  =ent(k,i)*dz(k)*p5     !for momentum

            qtn     =downqt(k+1,i)   *(one-entexp) + qt(k)*entexp
            thln    =downthl(k+1,i)  *(one-entexp) + thl(k)*entexp
            qncn    =downqnc(k+1,i)  *(one-entexp) + qnc(k)*entexp
            qnin    =downqni(k+1,i)  *(one-entexp) + qni(k)*entexp
            qnwfan  =downqnwfa(k+1,i)*(one-entexp) + qnwfa(k)*entexp
            qnifan  =downqnifa(k+1,i)*(one-entexp) + qnifa(k)*entexp
            un      =downu(k+1,i)    *(one-entexm) + u(k)*entexm
            vn      =downv(k+1,i)    *(one-entexm) + v(k)*entexm
            !qken=downqke(k-1,i)*(1.-entexp) + qke(k)*entexp
!            qtn =downqt(k+1,i) +(qt(k) -downqt(k+1,i)) *(1.-entexp)
!            thln=downthl(k+1,i)+(thl(k)-downthl(k+1,i))*(1.-entexp)
!            un  =downu(k+1,i)  +(u(k)  -downu(k+1,i))*(1.-entexp_m)
!            vn  =downv(k+1,i)  +(v(k)  -downv(k+1,i))*(1.-entexp_m)

            ! given new p & z, solve for thvn, qcn, and qin
            pk  =(p(k-1)*dz(k)+p(k)*dz(k-1))/(dz(k)+dz(k-1))
            call condensation_ddmf(qtn,thln,pk,zw(k),thvn,qcn,qin,debug_dd)

            !calculate buoyancy
            thvk =(thv(k-1)*dz(k)+thv(k)*dz(k-1))/(dz(k)+dz(k-1))
            buoy =grav*(thvn/thvk - one)
            !downdraft decelleration to limit large tendencies at the surface
            if (buoy > zero) then
               bcoeff =  0.2  !0.15 = same as in updrafts
               ecoeff = -2.0
            else
               !for numerical stability reasons, apply brakes near the surface.
               bcoeff = 0.2*(one - exp(-zw(k)/(0.15_kind_phys*pblh)))
               ecoeff = -2.0_kind_phys - exp(-zw(k)/(0.15_kind_phys*pblh))
               buoy   = buoy*(one - exp(-zw(k)/(0.15_kind_phys*pblh)))
            endif
            if (k==kts+2) buoy=max(buoy,0.01)
            if (k==kts+1) buoy=max(buoy,0.05)

!            beta_dm = 2*wb*ent(k,i) + 0.5/(zw(k)-dz(k)) * &
!                 &    max(1. - exp((zw(k) -dz(k))/z00 - 1. ) , 0.)
!            entw=exp(-beta_dm * dz(k))
!            entw=entexp
!            if (beta_dm >0) then
!               wn2=downw(k+1,i)**2*entw - wa*buoy/beta_dm * (1. - entw)
!            else
!               wn2=downw(k+1,i)**2      - 2.*wa*buoy*dz(k)
!            end if

            mindownw = min(downw(k+1,i),-0.2_kind_phys)
            wn = downw(k+1,i) + (ecoeff*ent(k,i)*downw(k+1,i)        &
                              - bcoeff*buoy/mindownw)*min(dz(k), 250.)

            !do not allow a parcel to accelerate more than 1.25 m/s over 200 m.
            !add max acceleration of -2.0 m/s for coarse vertical resolution.
            if (wn < downw(k+1,i) - min(1.25*dz(k)/200., 2.0))then
                wn = downw(k+1,i) - min(1.25*dz(k)/200., 2.0)
            endif
            !add symmetrical max decrease in velocity (less negative)
            if (wn > downw(k+1,i) + min(1.25*dz(k)/200., 2.0))then
                wn = downw(k+1,i) + min(1.25*dz(k)/200., 2.0)
            endif
            wn = max(min(wn,zero), -2.0_kind_phys)

            if (debug_dd .eq.1) then
               print *, "downdraft #:", i," diameter:",ddd,"==========="
               print *, "  k       =",      k,      " z    =", zw(k)
               print *, "  ent     =",ent(k,i),     " bouy =", buoy
               print *, "  downthv =",   thvn,      " thvk =", thvk
               print *, "  downthl =",   thln,      " thl  =", thl(k)
               print *, "  downqt  =",   qtn ,      " qt   =", qt(k)
               print *, "  downqc  =",   qcn ,      " qc   =", qc(k)
               print *, "  downw+1 =",downw(k+1,i), " wn2  =", wn
            endif

            if (wn .lt. 0.) then !terminate when velocity is too small
               downw(k,i)     = wn !-sqrt(wn2)
               downthv(k,i)   = thvn
               downthl(k,i)   = thln
               downqt(k,i)    = qtn
               downqc(k,i)    = qcn
               downqi(k,i)    = qin
               downqnc(k,i)   = qncn
               downqni(k,i)   = qnin
               downqnwfa(k,i) = qnwfan
               downqnifa(k,i) = qnifan
               downu(k,i)     = un
               downv(k,i)     = vn
               downa(k,i)     = downa(k+1,i)
            else
               !downdrafts must go at least 2 levels
!               if (dd_initk(i) - k .lt. 2) then
!                  downw(:,i)  = zero
!                  downthv(:,i)= zero
!                  downthl(:,i)= zero
!                  downqt(:,i) = zero
!                  downqc(:,i) = zero
!                  downu(:,i)  = zero
!                  downv(:,i)  = zero
!               endif
               exit
            endif
         enddo
      enddo

      downw(1,:) = zero !make sure downdrafts do not penetrate the surface
      downa(1,:) = zero

      !
      ! combine all downdrafts into one averaged downdraft
      !
      do k=qltop,kts,-1
         do i=1,ndd
            edmf_a_dd(k)  =edmf_a_dd(k)  +downa(k,i)
            edmf_w_dd(k)  =edmf_w_dd(k)  +downa(k,i)*downw(k,i)
            edmf_qt_dd(k) =edmf_qt_dd(k) +downa(k,i)*downqt(k,i)
            edmf_thl_dd(k)=edmf_thl_dd(k)+downa(k,i)*downthl(k,i)
            edmf_ent_dd(k)=edmf_ent_dd(k)+downa(k,i)*ent(k,i)
            edmf_qc_dd(k) =edmf_qc_dd(k) +downa(k,i)*downqc(k,i)
         enddo

         if (edmf_a_dd(k) > 0.) then
            edmf_w_dd(k)  =edmf_w_dd(k)  /edmf_a_dd(k)
            edmf_qt_dd(k) =edmf_qt_dd(k) /edmf_a_dd(k)
            edmf_thl_dd(k)=edmf_thl_dd(k)/edmf_a_dd(k)
            edmf_ent_dd(k)=edmf_ent_dd(k)/edmf_a_dd(k)
            edmf_qc_dd(k) =edmf_qc_dd(k) /edmf_a_dd(k)
         endif
         !instead of dTKE/dt = 1/2 w^3, multiply by 1 for QKE.
         tkeprod_dn(k)=(abs(edmf_w_dd(k))**3)*edmf_a_dd(k)/(b1*max(el(k),p5))
         tkeprod_dn(k)=min(0.0004_kind_phys, tkeprod_dn(k))
      enddo
      ! add tke source for entrainment at layer above cloud. use same area
      ! above cloud as used in the initialized downdraft area.
      tkeprod_dn(qltop+1)=went**3*edmf_a_dd(qltop)/(b1*max(el(qltop+1),p5))
      massflux = edmf_a_dd * edmf_w_dd
      maxmf_dd = abs(minval(massflux))

      if (debug_dd .eq. 2) then
         do k=kts,kte
            if (tkeprod_dn(k) > 0.0004_kind_phys .or. edmf_w_dd(k) < -two) then
               print*,"-------------------------------------"
               print*,"found strange behavior in downdrafts at",k
               print*,"qltop=",qltop," qlbase=",qlbase," a(top)=",edmf_a_dd(qltop)
               print*,"q*=",qstar," theta*=",thstar," went=",went
               print*,"f0=",f0," jump_thv=",jump_thv,"w*rad=",wstar_rad
               print*,"u*=",ust," jump_qt=",jump_qt," fltv=",fltv
               print*,"grav=",grav," pblh=",pblh," thv(1)=",thv(1)
               print*,"p(1)=",p(1)," jump_thl=",jump_thetal," ent_eff=",ent_eff
               print*,"sigmaw=",sigmaw," wmin=",wmin," wmax=",wmax
               print*,"tke prod=",tkeprod_dn(k)," w=",edmf_w_dd(k)," a(k)=",edmf_a_dd(k)
            endif
         enddo
      endif
      ! compute variables needed for solver
      !
      do k=kts,qltop
         rho_int = (rho(k)*dz(k+1)+rho(k+1)*dz(k))/(dz(k+1)+dz(k))
         do i=1,ndd
            sd_aw(k)   =sd_aw(k)   +rho_int*downa(k,i)*abs(downw(k,i))
            sd_awthl(k)=sd_awthl(k)+rho_int*downa(k,i)*downw(k,i)*downthl(k,i)
            sd_awqt(k) =sd_awqt(k) +rho_int*downa(k,i)*downw(k,i)*downqt(k,i)
            sd_awqc(k) =sd_awqc(k) +rho_int*downa(k,i)*downw(k,i)*downqc(k,i)
            sd_awqi(k) =sd_awqi(k) +rho_int*downa(k,i)*downw(k,i)*downqi(k,i)
            sd_awqnc(k)=sd_awqnc(k)+rho_int*downa(k,i)*downw(k,i)*downqnc(k,i)
            sd_awqni(k)=sd_awqni(k)+rho_int*downa(k,i)*downw(k,i)*downqni(k,i)
            sd_awqnwfa(k)=sd_awqnwfa(k)+rho_int*downa(k,i)*downw(k,i)*downqnwfa(k,i)
            sd_awqnifa(k)=sd_awqnifa(k)+rho_int*downa(k,i)*downw(k,i)*downqnifa(k,i)
            sd_awu(k)  =sd_awu(k)  +rho_int*downa(k,i)*downw(k,i)*downu(k,i)
            sd_awv(k)  =sd_awv(k)  +rho_int*downa(k,i)*downw(k,i)*downv(k,i)
         enddo
         sd_awqv(k) = sd_awqt(k)  - sd_awqc(k) - sd_awqi(k)
      enddo

   endif ! end cloud flag

end subroutine ddmp_mf
!=============================================================== 

SUBROUTINE SCALE_AWARE(dx,pblh,Psig_bl,Psig_shcu)

    !---------------------------------------------------------------
    !             NOTES ON SCALE-AWARE FORMULATION
    !
    !     scale-aware factor (Psig) here, taken from Honnert et al. (2011,
    !     JAS) and/or from Hyeyum Hailey Shin and Song-You Hong (2013, JAS)
    !
    ! Psig_bl tapers local mixing
    ! Psig_shcu tapers nonlocal mixing

  use module_bl_mynnedmf_common, only: zero,one,ten,kind_phys
  
    real(kind_phys), intent(in)  :: dx,pblh
    real(kind_phys), intent(out) :: Psig_bl,Psig_shcu
    real(kind_phys)              :: dxdh

    Psig_bl=one
    Psig_shcu=one
    dxdh=MAX(2.5_kind_phys*dx,ten)/MIN(PBLH,3000._kind_phys)
    ! Honnert et al. 2011, TKE in PBL  *** original form used until 201605
    !Psig_bl= ((dxdh**2) + 0.07*(dxdh**0.667))/((dxdh**2) + &
    !         (3./21.)*(dxdh**0.67) + (3./42.))
    ! Honnert et al. 2011, TKE in entrainment layer
    !Psig_bl= ((dxdh**2) + (4./21.)*(dxdh**0.667))/((dxdh**2) + &
     !        (3./20.)*(dxdh**0.67) + (7./21.))
    ! New form to preseve parameterized mixing - only down 5% at dx = 750 m
    Psig_bl= ((dxdh**2) + 0.106_kind_phys*(dxdh**0.667_kind_phys))/((dxdh**2) + &
         0.066_kind_phys*(dxdh**0.667_kind_phys) + 0.071_kind_phys)

    !assume a 500 m cloud depth for shallow-cu clods
    dxdh=MAX(2.5_kind_phys*dx,ten)/MIN(PBLH+500._kind_phys,3500._kind_phys)
    ! Honnert et al. 2011, TKE in entrainment layer *** original form used until 201605
    !Psig_shcu= ((dxdh**2) + (4./21.)*(dxdh**0.667))/((dxdh**2) + &
    !         (3./20.)*(dxdh**0.67) + (7./21.))

    ! Honnert et al. 2011, TKE in cumulus
    !Psig(i)= ((dxdh**2) + 1.67*(dxdh**1.4))/((dxdh**2) +1.66*(dxdh**1.4) +
    !0.2)

    ! Honnert et al. 2011, w'q' in PBL
    !Psig(i)= 0.5 + 0.5*((dxdh**2) + 0.03*(dxdh**1.4) -
    !(4./13.))/((dxdh**2) + 0.03*(dxdh**1.4) + (4./13.))
    ! Honnert et al. 2011, w'q' in cumulus
    !Psig(i)= ((dxdh**2) - 0.07*(dxdh**1.4))/((dxdh**2) -0.07*(dxdh**1.4) +
    !0.02)

    ! Honnert et al. 2011, q'q' in PBL
    !Psig(i)= 0.5 + 0.5*((dxdh**2) + 0.25*(dxdh**0.667) -0.73)/((dxdh**2)
    !-0.03*(dxdh**0.667) + 0.73)
    ! Honnert et al. 2011, q'q' in cumulus
    !Psig(i)= ((dxdh**2) - 0.34*(dxdh**1.4))/((dxdh**2) - 0.35*(dxdh**1.4)
    !+ 0.37)

    ! Hyeyum Hailey Shin and Song-You Hong 2013, TKE in PBL (same as Honnert's above)
    !Psig_shcu= ((dxdh**2) + 0.070*(dxdh**0.667))/((dxdh**2)
    !+0.142*(dxdh**0.667) + 0.071)
    ! Hyeyum Hailey Shin and Song-You Hong 2013, TKE in entrainment zone  *** switch to this form 201605
    Psig_shcu= ((dxdh**2) + 0.145_kind_phys*(dxdh**0.667_kind_phys))/((dxdh**2) + &
         0.172_kind_phys*(dxdh**0.667_kind_phys) + 0.170_kind_phys)

    ! Hyeyum Hailey Shin and Song-You Hong 2013, w'theta' in PBL
    !Psig(i)= 0.5 + 0.5*((dxdh**2) -0.098)/((dxdh**2) + 0.106) 
    ! Hyeyum Hailey Shin and Song-You Hong 2013, w'theta' in entrainment zone
    !Psig(i)= 0.5 + 0.5*((dxdh**2) - 0.112*(dxdh**0.25) -0.071)/((dxdh**2)
    !+ 0.054*(dxdh**0.25) + 0.10)

    !print*,"in scale_aware; dx, dxdh, Psig(i)=",dx,dxdh,Psig(i)
    !If(Psig_bl(i) < 0.0 .OR. Psig(i) > 1.)print*,"dx, dxdh, Psig(i)=",dx,dxdh,Psig_bl(i) 
    If(Psig_bl > one) Psig_bl=one
    If(Psig_bl < zero) Psig_bl=zero

    If(Psig_shcu > one) Psig_shcu=one
    If(Psig_shcu < zero) Psig_shcu=zero

  END SUBROUTINE SCALE_AWARE

! =====================================================================
!>\ingroup gsd_mynn_edmf
!! \author JAYMES- added 22 Apr 2015
!! This function calculates saturation vapor pressure.  Separate ice and liquid functions
!! are used (identical to those in module_mp_thompson.F, v3.6). Then, the
!! final returned value is a temperature-dependant "blend". Because the final
!! value is "phase-aware", this formulation may be preferred for use throughout
!! the module (replacing "svp").
  FUNCTION esat_blend(t) 

      use module_bl_mynnedmf_common, only: t0c,tice,one,kind_phys
      
      real(kind_phys), intent(in):: t
      real(kind_phys):: esat_blend,XC,ESL,ESI,chi
      !liquid
      real(kind_phys), parameter:: J0= .611583699E03
      real(kind_phys), parameter:: J1= .444606896E02
      real(kind_phys), parameter:: J2= .143177157E01
      real(kind_phys), parameter:: J3= .264224321E-1
      real(kind_phys), parameter:: J4= .299291081E-3
      real(kind_phys), parameter:: J5= .203154182E-5
      real(kind_phys), parameter:: J6= .702620698E-8
      real(kind_phys), parameter:: J7= .379534310E-11
      real(kind_phys), parameter:: J8=-.321582393E-13
      !ice
      real(kind_phys), parameter:: K0= .609868993E03
      real(kind_phys), parameter:: K1= .499320233E02
      real(kind_phys), parameter:: K2= .184672631E01
      real(kind_phys), parameter:: K3= .402737184E-1
      real(kind_phys), parameter:: K4= .565392987E-3
      real(kind_phys), parameter:: K5= .521693933E-5
      real(kind_phys), parameter:: K6= .307839583E-7
      real(kind_phys), parameter:: K7= .105785160E-9
      real(kind_phys), parameter:: K8= .161444444E-12

      XC=MAX(-80._kind_phys, t - t0c) !note t0c = 273.15, tice is set in module mynn_common to 240

! For 240 < t < 268.16 K, the vapor pressures are "blended" as a function of temperature, 
! using the approach similar to Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting 
! values are returned from the function.
      IF (t .GE. (t0c-6._kind_phys)) THEN
          esat_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
      ELSE IF (t .LE. tice) THEN
          esat_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
      ELSE
          ESL = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          chi = ((t0c-6._kind_phys) - t)/((t0c-6._kind_phys) - tice)
          esat_blend = (one-chi)*ESL  + chi*ESI
      END IF

  END FUNCTION esat_blend

! ====================================================================

!>\ingroup gsd_mynn_edmf
!! This function extends function "esat" and returns a "blended"
!! saturation mixing ratio. Tice currently set to 240 K, t0c = 273.15 K.
!!\author JAYMES
  FUNCTION qsat_blend(t, P)

      use module_bl_mynnedmf_common, only: t0c,tice,one,kind_phys

      real(kind_phys), intent(in):: t, P
      real(kind_phys):: qsat_blend,XC,ESL,ESI,RSLF,RSIF,chi
      !liquid
      real(kind_phys), parameter:: J0= .611583699E03
      real(kind_phys), parameter:: J1= .444606896E02
      real(kind_phys), parameter:: J2= .143177157E01
      real(kind_phys), parameter:: J3= .264224321E-1
      real(kind_phys), parameter:: J4= .299291081E-3
      real(kind_phys), parameter:: J5= .203154182E-5
      real(kind_phys), parameter:: J6= .702620698E-8
      real(kind_phys), parameter:: J7= .379534310E-11
      real(kind_phys), parameter:: J8=-.321582393E-13
      !ice
      real(kind_phys), parameter:: K0= .609868993E03
      real(kind_phys), parameter:: K1= .499320233E02
      real(kind_phys), parameter:: K2= .184672631E01
      real(kind_phys), parameter:: K3= .402737184E-1
      real(kind_phys), parameter:: K4= .565392987E-3
      real(kind_phys), parameter:: K5= .521693933E-5
      real(kind_phys), parameter:: K6= .307839583E-7
      real(kind_phys), parameter:: K7= .105785160E-9
      real(kind_phys), parameter:: K8= .161444444E-12

      XC=MAX(-80._kind_phys,t - t0c)

      IF (t .GE. (t0c-6._kind_phys)) THEN
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESL  = min(ESL, P*0.15_kind_phys) ! Even with P=1050mb and T=55C, the sat. vap. pres only contributes to ~15% of total pres.
          qsat_blend = 0.622_kind_phys*ESL/max(P-ESL, 1e-5_kind_phys) 
      ELSE IF (t .LE. tice) THEN
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          ESI  = min(ESI, P*0.15_kind_phys)
          qsat_blend = 0.622*ESI/max(P-ESI, 1e-5_kind_phys)
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESL  = min(ESL, P*0.15_kind_phys)
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          ESI  = min(ESI, P*0.15_kind_phys)
          RSLF = 0.622_kind_phys*ESL/max(P-ESL, 1e-5_kind_phys)
          RSIF = 0.622_kind_phys*ESI/max(P-ESI, 1e-5_kind_phys)
!          chi  = (268.16-t)/(268.16-240.)
          chi  = ((t0c-6._kind_phys) - t)/((t0c-6._kind_phys) - tice) 
          qsat_blend = (one-chi)*RSLF + chi*RSIF
      END IF

  END FUNCTION qsat_blend

! ===================================================================

!>\ingroup gsd_mynn_edmf
!! This function interpolates the latent heats of vaporization and sublimation into
!! a single, temperature-dependent, "blended" value, following 
!! Chaboureau and Bechtold (2002) \cite Chaboureau_2002, Appendix.
!!\author JAYMES
  FUNCTION xl_blend(t)

    use module_bl_mynnedmf_common, only: xlv,xls,cpv,cliq,cice,t0c,tice,&
         one,kind_phys

      real(kind_phys), intent(in):: t
      real(kind_phys):: xl_blend,xlvt,xlst,chi
      !note: t0c = 273.15, tice is set in mynn_common

      IF (t .GE. t0c) THEN
          xl_blend = xlv + (cpv-cliq)*(t-t0c)  !vaporization/condensation
      ELSE IF (t .LE. tice) THEN
          xl_blend = xls + (cpv-cice)*(t-t0c)  !sublimation/deposition
      ELSE
          xlvt = xlv + (cpv-cliq)*(t-t0c)  !vaporization/condensation
          xlst = xls + (cpv-cice)*(t-t0c)  !sublimation/deposition
!          chi  = (273.16-t)/(273.16-240.)
          chi  = (t0c - t)/(t0c - tice)
          xl_blend = (one-chi)*xlvt + chi*xlst     !blended
      END IF

  END FUNCTION xl_blend

! ===================================================================

  FUNCTION phim(zet)
     ! New stability function parameters for momentum (Puhales, 2020, WRF 4.2.1)
     ! The forms in unstable conditions (z/L < 0) use Grachev et al. (2000), which are a blend of 
     ! the classical (Kansas) forms (i.e., Paulson 1970, Dyer and Hicks 1970), valid for weakly 
     ! unstable conditions (-1 < z/L < 0). The stability functions for stable conditions use an
     ! updated form taken from Cheng and Brutsaert (2005), which extends the validity into very
     ! stable conditions [z/L ~ O(10)].
      use module_bl_mynnedmf_common, only: cphm_unst,zero,one,two,p333,p5,p666,kind_phys

      real(kind_phys), intent(in):: zet
      real(kind_phys):: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      real(kind_phys), parameter :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      real(kind_phys), parameter :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      real(kind_phys), parameter :: am_unst=10., ah_unst=34.
      real(kind_phys):: phi_m,phim

      if ( zet >= zero ) then
         dummy_0 = one+zet**bm_st
         dummy_1 = zet+dummy_0**(rbm_st)
         dummy_11= one+dummy_0**(rbm_st-1)*zet**(bm_st-one)
         dummy_2 = (-am_st/dummy_1)*dummy_11
         phi_m   = one-zet*dummy_2
      else
         dummy_0 = (one-cphm_unst*zet)**0.25
         phi_m   = one/dummy_0
         dummy_psi = two*log(p5*(one+dummy_0))+log(p5*(one+dummy_0**2))-two*atan(dummy_0)+1.570796

         dummy_0 = (one-am_unst*zet)          ! parentesis arg
         dummy_1 = dummy_0**p333          ! y
         dummy_11=-p333*am_unst*dummy_0**(-p666) ! dy/dzet
         dummy_2 = p333*(dummy_1**2+dummy_1+one)      ! f
         dummy_22= p333*dummy_11*(two*dummy_1+one)    ! df/dzet
         dummy_3 = 0.57735*(two*dummy_1+one)  ! g
         dummy_33= 1.1547*dummy_11            ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44= (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(one+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = one/(one+dummy_0) ! denon
         dummy_11= two*zet         ! denon/dzet
         dummy_2 = ((one-phi_m)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22= -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phi_m = one-zet*(dummy_2+dummy_22)
      end if

      !phim = phi_m - zet
      phim = phi_m

  END FUNCTION phim
! ===================================================================

  FUNCTION phih(zet)
    ! New stability function parameters for heat (Puhales, 2020, WRF 4.2.1)
    ! The forms in unstable conditions (z/L < 0) use Grachev et al. (2000), which are a blend of
    ! the classical (Kansas) forms (i.e., Paulson 1970, Dyer and Hicks 1970), valid for weakly
    ! unstable conditions (-1 < z/L < 0). The stability functions for stable conditions use an
    ! updated form taken from Cheng and Brutsaert (2005), which extends the validity into very
    ! stable conditions [z/L ~ O(10)].
      use module_bl_mynnedmf_common, only: cphh_unst,zero,one,two,p333,p5,p666,kind_phys

      real(kind_phys), intent(in):: zet
      real(kind_phys):: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      real(kind_phys), parameter :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      real(kind_phys), parameter :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      real(kind_phys), parameter :: am_unst=10., ah_unst=34.
      real(kind_phys):: phh,phih

      if ( zet >= zero ) then
         dummy_0 = one+zet**bh_st
         dummy_1 = zet+dummy_0**(rbh_st)
         dummy_11= one+dummy_0**(rbh_st-one)*zet**(bh_st-one)
         dummy_2 = (-ah_st/dummy_1)*dummy_11
         phih    = one-zet*dummy_2
      else
         dummy_0 = (one-cphh_unst*zet)**p5
         phh     = one/dummy_0
         dummy_psi = two*log(p5*(one+dummy_0))

         dummy_0 = (one-ah_unst*zet)         ! parenthesis arg
         dummy_1 = dummy_0**p333         ! y
         dummy_11=-p333*ah_unst*dummy_0**(-p666) ! dy/dzet
         dummy_2 = p333*(dummy_1**2+dummy_1+one)      ! f
         dummy_22= p333*dummy_11*(two*dummy_1+one)    ! df/dzet
         dummy_3 = 0.57735*(two*dummy_1+one) ! g
         dummy_33= 1.1547*dummy_11           ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44= (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(one+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = one/(one+dummy_0)         ! denon
         dummy_11= two*zet                   ! ddenon/dzet
         dummy_2 = ((one-phh)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22= -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phih = one-zet*(dummy_2+dummy_22)
      end if

end function phih
! ==================================================================
 subroutine topdown_cloudrad(kts,kte,                         &
               &dz1,zw,fltv,u1,v1,xland,kpbl,pblh,            &
               &sqc,sqi,sqw,thl,th1,ex1,pres1,rho1,thv,       &
               &cldfra_bl1,qc_bl,qi_bl,rthraten,              &
               &tkeprod_dn,psig,                              &
               &maxtkeprod,cldtop_cooling,ent_eff)

   use module_bl_mynnedmf_common, only: cp,ep_2,xlv,r_d,p608, &
        xlvcp,grav,zero,one,two,three,eight,p1,p2,p25,p3,     &
        p333,p5,kind_phys
   
    !input
    integer,         intent(in) :: kte,kts
    real(kind_phys), dimension(kts:kte), intent(in) :: dz1,sqc,sqi,sqw,&
          thl,th1,ex1,pres1,rho1,thv,cldfra_bl1,qc_bl,qi_bl
    real(kind_phys), dimension(kts:kte), intent(in) :: rthraten
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in) :: pblh,fltv,u1,v1,psig
    real(kind_phys), intent(in) :: xland
    integer        , intent(in) :: kpbl
    !output - production of tke from cloud-top radiative cooling
    real(kind_phys), dimension(kts:kte), intent(inout) :: tkeprod_dn
    real(kind_phys), intent(inout) :: maxtkeprod
    real(kind_phys), intent(inout) :: cldtop_cooling
    real(kind_phys), intent(inout) :: ent_eff
    !local
    real(kind_phys),parameter:: dz500   = 500.  !scale height for pbl
    real(kind_phys),parameter:: dz200   = 200.  !scale height for entrainment layer
    real(kind_phys),parameter:: cooling_thresh = -0.000069  !activation threshold of turbulence
                                                !-0.000116 is ~ -10 C cooling at cloud top per day
                                                !-0.000093 is ~  -8 C cooling at cloud top per day
                                                !-0.000069 is ~  -6 C cooling at cloud top per day
    real(kind_phys), dimension(kts:kte) :: zfac,zfacent,tkeprodorig
    real(kind_phys) :: dthvx,tmp1,zfacent_max,zagl
    real(kind_phys) :: temps,templ,zl1,zb1,dz_ent,wstar_rad,maxcooling
    real(kind_phys) :: f0,radflux,went,rcldb,rvls,minrad,zminrad,dp
    real(kind_phys) :: wspd_pbl,ac_wsp
    real(kind_phys), parameter :: pfac =2.0, zfmin = 0.01, phifac=8.0
    integer :: k,kk,kminrad
    logical :: cloudflg
    logical,parameter:: debug=.false.

    cloudflg       = .false.
    minrad         = zero
    kminrad        = kpbl
    zminrad        = pblh
    !save the tke production from the downdraft scheme for comparison
    tkeprodorig(kts:kte)=tkeprod_dn(kts:kte)
    
    !check for stratocumulus- or fog-tops with cooling at top
    do k = 1,min(kpbl+10,kte-1)
       maxcooling = minval(rthraten(max(kts,k-1):k+1))
       !---------------------------------criteria for downdraft activa tion:
       if ((cldfra_bl1(k).gt.p5 .and. cldfra_bl1(k+1).lt.p5) .and. &  !1) stratocu cloud top exists
           (maxcooling  .lt. cooling_thresh)) then                    !2) significant radiative cooling
          ! found sc cloud with significant radiative cooling
          cloudflg =.true.
          if (maxcooling < minrad) then
             minrad    =maxcooling      !maximum cooling rate at cloud top
             kminrad  = k
             zminrad  = zw(k) + p5*dz1(k) !Best estimate of height of TKE source (top of downdrafts)
          endif
       endif
    enddo

    if (cloudflg) then
       zl1     = p25*dz1(kts)
       k       = kminrad

       templ   = thl(k)*ex1(k)
       !rvls is ws at full level
       rvls    = 100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep_2/pres1(k+1))
       temps   = templ + (sqw(k)-rvls)/(cp/xlv  +  ep_2*xlv*rvls/(r_d*templ**2))
       rvls    = 100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep_2/pres1(k+1))
       rcldb   = max(sqw(k)-rvls,zero)         !resolved cloud water
       rcldb   = max(rcldb, qc_bl(k)+qi_bl(k)) !any type of cloud water (sgs or resolved)

       !entrainment efficiency
       dthvx   = (thl(k+2) + th1(k+2)*p608*sqw(k+2)) - &
                 (thl(k)   + th1(k)  *p608*sqw(k))
       !reduce vertical resolution sensitivity; normalize to 200 m depth
       dz_ent  = (dz1(k  )*dz1(k+1)+dz1(k+1)*dz1(k  ))/(dz1(k+1)+dz1(k)) + &
                 (dz1(k+1)*dz1(k+2)+dz1(k+2)*dz1(k+1))/(dz1(k+1)+dz1(k+2)) 
       dthvx   = dthvx * dz200/dz_ent  !normalized gradient
       dthvx   = max(dthvx, p1)
       tmp1    = xlvcp * rcldb/(ex1(k)*dthvx)
       !entrainment efficiency: originally from Nichols and Turton (1986), where a2 = 60, but lowered
       !here to 8, as in Grenier and Bretherton (2001).
       ent_eff = min(two, p2 + p2*eight*tmp1)

       f0      = zero
       do kk = max(kts,kminrad-1),kminrad+1
          radflux = rthraten(kk)*ex1(kk)             !converts theta/s to temp/s
          dp      = rho1(k)*grav*dz1(k)
          radflux = radflux * cp / grav * dp         !converts temp/s to w/m^2 over layer
          if (radflux < zero ) then
             f0   = abs(radflux)+f0                  !sum of cooling (w/m^2) over relevant layers
          endif
       enddo

       !more strict limits over land to reduce stable-layer mixouts
       if ((xland-1.5).ge.zero) then     ! water
          f0     = min(f0, 150.0_kind_phys)
       else                              ! land - limited considerably
          f0     = min(p25*f0, 30._kind_phys)
       endif
       cldtop_cooling = abs(minrad) * 3600.   ! (K/hr) max cooling 
       
       !entrainment from pbl top thermals
       wstar_rad = 1.25_kind_phys * ( grav/thv(k) * dz200 * f0 / (rho1(k) * cp) ) **p333
       wstar_rad = min(max(wstar_rad, p1), three)
       dthvx     = max(thv(k+1)-thv(k), p1)
       went      = p333*ent_eff*wstar_rad     !entrainment velocity (m/s)

       !compute normalized analytical profiles (max=1). constrain the profile of tke production to not
       !always go to the surface if forcing isnt strong enough, i.e., depth ~ wstar * Timescale of descent
       zb1       = max(zl1, zminrad - wstar_rad*800._kind_phys)
       do k = kts,kminrad+3
          !for fog at k=1, make min height above ground (zagl) larger than first model deptk
          zagl       = zw(k) + p5*dz1(k)
          if (k==kts .and. kminrad==kts) zagl = 0.4_kind_phys * dz1(kts)
          !analytic vertical profile
          zfac(k)    = min(max((one - max(zero,zagl-zb1)/(zminrad-zl1)), zfmin), one)
          zfacent(k) = max((zminrad-zagl)/zminrad, zero)*(one-zfac(k))**3
       enddo
       zfacent_max = maxval(zfacent(kts:kminrad+3))
       zfacent     = zfacent/max(zfacent_max, 1e-5_kind_phys)  !normalize zfacent

       !calculate tke production
       do k = kts,kminrad+3
          tkeprod_dn(k) = max(tkeprod_dn(k), p3*ent_eff*wstar_rad**3/dz500*zfacent(k))
          tkeprod_dn(k) = max(tkeprod_dn(k), zero)
          tkeprod_dn(k) = min(0.0004_kind_phys, tkeprod_dn(k))
       enddo
       !make sure there is some TKE source at the max cooling level
       tkeprod_dn(kminrad)=max(tkeprod_dn(kminrad), p333*tkeprod_dn(max(kts,kminrad-1)))
       !make sure there is a small TKE source above the max cooling level
       tkeprod_dn(kminrad+1)=max(tkeprod_dn(kminrad+1), min(p5,ent_eff)*tkeprod_dn(max(kts,kminrad)))
       
       !taper off tke production in high-wind conditions and at high resolutions
       wspd_pbl   = sqrt(max(u1**2 + v1**2, 0.01_kind_phys))
       ac_wsp     = one - min(max(zero, wspd_pbl - 10._kind_phys)/15._kind_phys, one)
       tkeprod_dn = tkeprod_dn * min(ac_wsp, psig)
    endif !end cloud check

    !calculate the maximum tke production in the profile and convert to an hourly rate.
    maxtkeprod = maxval(tkeprod_dn)*3600.
    
    if (cloudflg .and. debug) then
       !debug printouts
       print *, '-----------------------------------------------------------------'
       if ((xland-1.5).ge.zero) then
          print *, ' over water, kminrad=',kminrad
       else
          print *, ' over land, kminrad=',kminrad
       endif
       print *, ' ctop cooling (K/hr)=',minrad*3600.,' rad cooling (W/m^2)=',f0
       print *, ' ent_eff=',ent_eff,' went=',went,' w*rad=',wstar_rad
       print *, '  k    cldfra   dT/dt(hr)  dTKE/dt_orig dTKE/dt(hr)  zfacent'
       print *, '-----------------------------------------------------------------'
       do k = kminrad+3,kts,-1
          print '(I4,5F11.5)', k, cldfra_bl1(k), rthraten(k)*3600., tkeprodorig(k)*3600., tkeprod_dn(k)*3600.,  zfacent(k) 
       enddo
    endif
    
 END SUBROUTINE topdown_cloudrad
! ==================================================================
! ===================================================================
! ===================================================================

END MODULE module_bl_mynnedmf
