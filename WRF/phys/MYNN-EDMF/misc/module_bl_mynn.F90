!>\file module_bl_mynn.F90
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
!***********************************************************************
! ==================================================================
! Notes on original implementation into WRF-ARW
! changes to original code:
! 1. code is 1D (in z)
! 2. option to advect TKE, but not the covariances and variances
! 3. Cranck-Nicholson replaced with the implicit scheme
! 4. removed terrain-dependent grid since input in WRF in actual
!    distances in z[m]
! 5. cosmetic changes to adhere to WRF standard (remove common blocks,
!            intent etc)
!-------------------------------------------------------------------
! Further modifications post-implementation
!
! 1. Addition of BouLac mixing length in the free atmosphere.
! 2. Changed the turbulent mixing length to be integrated from the
!    surface to the top of the BL + a transition layer depth.
! v3.4.1:    Option to use Kitamura/Canuto modification which removes 
!            the critical Richardson number and negative TKE (default).
!            Hybrid PBL height diagnostic, which blends a theta-v-based
!            definition in neutral/convective BL and a TKE-based definition
!            in stable conditions.
!            TKE budget output option
! v3.5.0:    TKE advection option (bl_mynn_tkeadvect)
! v3.5.1:    Fog deposition related changes.
! v3.6.0:    Removed fog deposition from the calculation of tendencies
!            Added mixing of qc, qi, qni
!            Added output for wstar, delta, TKE_PBL, & KPBL for correct 
!                   coupling to shcu schemes  
! v3.8.0:    Added subgrid scale cloud output for coupling to radiation
!            schemes (activated by setting icloud_bl =1 in phys namelist).
!            Added WRF_DEBUG prints (at level 3000)
!            Added Tripoli and Cotton (1981) correction.
!            Added namelist option bl_mynn_cloudmix to test effect of mixing
!                cloud species (default = 1: on). 
!            Added mass-flux option (bl_mynn_edmf, = 1 for DMP mass-flux, 0: off).
!                Related options: 
!                 bl_mynn_edmf_mom = 1 : activate momentum transport in MF scheme
!                 bl_mynn_edmf_tke = 1 : activate TKE transport in MF scheme
!            Added mixing length option (bl_mynn_mixlength, see notes below)
!            Added more sophisticated saturation checks, following Thompson scheme
!            Added new cloud PDF option (bl_mynn_cloudpdf = 2) from Chaboureau
!                and Bechtold (2002, JAS, with mods) 
!            Added capability to mix chemical species when env variable
!                WRF_CHEM = 1, thanks to Wayne Angevine.
!            Added scale-aware mixing length, following Junshi Ito's work
!                Ito et al. (2015, BLM).
! v3.9.0    Improvement to the mass-flux scheme (dynamic number of plumes,
!                better plume/cloud depth, significant speed up, better cloud
!                fraction). 
!            Added Stochastic Parameter Perturbation (SPP) implementation.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid clouds.
! v.4.0      Removed or added alternatives to WRF-specific functions/modules
!                for the sake of portability to other models.
!                the sake of portability to other models.
!            Further refinement of mass-flux scheme from SCM experiments with
!                Wayne Angevine: switch to linear entrainment and back to
!                Simpson and Wiggert-type w-equation.
!            Addition of TKE production due to radiation cooling at top of 
!                clouds (proto-version); not activated by default.
!            Some code rewrites to move if-thens out of loops in an attempt to
!                improve computational efficiency.
!            New tridiagonal solver, which is supposedly 14% faster and more
!                conservative. Impact seems very small.
!            Many miscellaneous tweaks to the mixing lengths and stratus
!                component of the subgrid-scale (SGS) clouds.
! v4.1       Big improvements in downward SW radiation due to revision of subgrid clouds
!                - better cloud fraction and subgrid scale mixing ratios.
!                - may experience a small cool bias during the daytime now that high 
!                  SW-down bias is greatly reduced...
!            Some tweaks to increase the turbulent mixing during the daytime for
!                bl_mynn_mixlength option 2 to alleviate cool bias (very small impact).
!            Improved ensemble spread from changes to SPP in MYNN
!                - now perturbing eddy diffusivity and eddy viscosity directly
!                - now perturbing background rh (in SGS cloud calc only)
!                - now perturbing entrainment rates in mass-flux scheme
!            Added IF checks (within IFDEFS) to protect mixchem code from being used
!                when HRRR smoke is used (no impact on regular non-wrf chem use)
!            Important bug fix for wrf chem when transporting chemical species in MF scheme
!            Removed 2nd mass-flux scheme (no only bl_mynn_edmf = 1, no option 2)
!            Removed unused stochastic code for mass-flux scheme
!            Changed mass-flux scheme to be integrated on interface levels instead of
!                mass levels - impact is small
!            Added option to mix 2nd moments in MYNN as opposed to the scalar_pblmix option.
!                - activated with bl_mynn_mixscalars = 1; this sets scalar_pblmix = 0
!                - added tridagonal solver used in scalar_pblmix option to duplicate tendencies
!                - this alone changes the interface call considerably from v4.0.
!            Slight revision to TKE production due to radiation cooling at top of clouds
!            Added the non-Guassian buoyancy flux function of Bechtold and Siebesma (1998, JAS).
!                - improves TKE in SGS clouds
!            Added heating due to dissipation of TKE (small impact, maybe + 0.1 C daytime PBL temp)
!            Misc changes made for FV3/MPAS compatibility
! v4.2       A series of small tweaks to help reduce a cold bias in the PBL:
!                - slight increase in diffusion in convective conditions
!                - relaxed criteria for mass-flux activation/strength
!                - added capability to cycle TKE for continuity in hourly updating HRRR
!                - added effects of compensational environmental subsidence in mass-flux scheme,
!                  which resulted in tweaks to detrainment rates.
!            Bug fix for diagnostic-decay of SGS clouds - noticed by Greg Thompson. This has
!                a very small, but primarily  positive, impact on SW-down biases.
!            Tweak to calculation of KPBL - urged by Laura Fowler - to make more intuitive.
!            Tweak to temperature range of blending for saturation check (water to ice). This
!                slightly reduces excessive SGS clouds in polar region. No impact warm clouds. 
!            Added namelist option bl_mynn_output (0 or 1) to suppress or activate the
!                allocation and output of 10 3D variables. Most people will want this
!                set to 0 (default) to save memory and disk space.
!            Added new array qi_bl as opposed to using qc_bl for both SGS qc and qi. This
!                gives us more control of the magnitudes which can be confounded by using
!                a single array. As a results, many subroutines needed to be modified,
!                especially mym_condensation.
!            Added the blending of the stratus component of the SGS clouds to the mass-flux
!                clouds to account for situations where stratus and cumulus may exist in the
!                grid cell.
!            Misc small-impact bugfixes:
!                1) dz was incorrectly indexed in mym_condensation
!                2) configurations with icloud_bl = 0 were using uninitialized arrays
! v4.5 / CCPP
!            This version includes many modifications that proved valuable in the global
!            framework and removes some key lingering bugs in the mixing of chemical species.
!            TKE Budget output fixed (Puhales, 2020-12)
!            New option for stability function: (Puhales, 2020-12)
!                bl_mynn_stfunc = 0 (original, Kansas-type function, Paulson, 1970 )
!                bl_mynn_stfunc = 1 (expanded range, same as used for Jimenez et al (MWR)
!                see the Technical Note for this implementation (small impact).
!            Improved conservation of momentum and higher-order moments.
!            Important bug fixes for mixing of chemical species.
!            Addition of pressure-gradient effects on updraft momentum transport.
!            Addition of bl_mynn_closure option = 2.5, 2.6, or 3.0
!            Addition of higher-order moments for sigma when using 
!                bl_mynn_cloudpdf = 2 (Chab-Becht).
!            Removed WRF_CHEM dependencies.
!            Many miscellaneous tweaks.
! v4.6 / CCPP
!            Some code optimization. Removed many conditions from loops. Redesigned the mass-
!                flux scheme to use 8 plumes instead of a variable n plumes. This results in
!                the removal of the output variable "nudprafts" and adds maxwidth and ztop_plume.
!            Revision option bl_mynn_cloudpdf = 2, which now ensures cloud fractions for all
!                optically relevant mixing ratios (tip from Greg Thompson). Also, added flexibility
!                for tuning near-surface cloud fractions to remove excess fog/low ceilings.
!            Now outputs all SGS cloud mixing ratios as grid-mean values, not in-cloud. This 
!                results in a change in the pre-radiation code to no longer multiply mixing ratios
!                by cloud fractions.
!            Bug fix for the momentum transport.
!            Lots of code cleanup: removal of test code, comments, changing text case, etc.
!            Many misc tuning/tweaks.
!
! Many of these changes are now documented in references listed above.
!====================================================================

MODULE module_bl_mynn

 use module_bl_mynn_common,only:                        &
       cp        , cpv       , cliq       , cice      , &
       p608      , ep_2      , ep_3       , gtr       , &
       grav      , g_inv     , karman     , p1000mb   , &
       rcp       , r_d       , r_v        , rk        , &
       rvovrd    , svp1      , svp2       , svp3      , &
       xlf       , xlv       , xls        , xlscp     , &
       xlvcp     , tv0       , tv1        , tref      , &
       zero      , half      , one        , two       , &
       onethird  , twothirds , tkmin      , t0c       , &
       tice      , kind_phys


 IMPLICIT NONE

!===================================================================
! From here on, these are MYNN-specific parameters:
! The parameters below depend on stability functions of module_sf_mynn.
 real(kind_phys), parameter :: cphm_st=5.0, cphm_unst=16.0, &
                               cphh_st=5.0, cphh_unst=16.0

! Closure constants
 real(kind_phys), parameter ::  &
      &pr  =  0.74,             &
      &g1  =  0.235,            &  ! NN2009 = 0.235
      &b1  = 24.0,              &
      &b2  = 15.0,              &  ! CKmod     NN2009
      &c2  =  0.729,            &  ! 0.729, & !0.75, &
      &c3  =  0.340,            &  ! 0.340, & !0.352, &
      &c4  =  0.0,              &
      &c5  =  0.2,              &
      &a1  = b1*( 1.0-3.0*g1 )/6.0, &
!      &c1  = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) ), &
      &c1  = g1 -1.0/( 3.0*a1*2.88449914061481660), &
      &a2  = a1*( g1-c1 )/( g1*pr ), &
      &g2  = b2/b1*( 1.0-c3 ) +2.0*a1/b1*( 3.0-2.0*c2 )

 real(kind_phys), parameter ::  &
      &cc2 =  1.0-c2,           &
      &cc3 =  1.0-c3,           &
      &e1c =  3.0*a2*b2*cc3,    &
      &e2c =  9.0*a1*a2*cc2,    &
      &e3c =  9.0*a2*a2*cc2*( 1.0-c5 ), &
      &e4c = 12.0*a1*a2*cc2,    &
      &e5c =  6.0*a1*a1

! Constants for min tke in elt integration (qmin), max z/L in els (zmax), 
! and factor for eddy viscosity for TKE (Kq = Sqfac*Km):
 real(kind_phys), parameter :: qmin=0.0, zmax=1.0, Sqfac=3.0
! Note that the following mixing-length constants are now specified in mym_length
!     &cns=3.5, alp1=0.23, alp2=0.3, alp3=3.0, alp4=10.0, alp5=0.2

 real(kind_phys), parameter :: qkemin=1.e-5
 real(kind_phys), parameter :: tliq = 269. !all hydrometeors are liquid when T > tliq

! Constants for cloud PDF (mym_condensation)
 real(kind_phys), parameter :: rr2=0.7071068, rrp=0.3989423

!>Use Canuto/Kitamura mod (remove Ric and negative TKE) (1:yes, 0:no)
!!For more info, see Canuto et al. (2008 JAS) and Kitamura (Journal of the 
!!Meteorological Society of Japan, Vol. 88, No. 5, pp. 857-864, 2010).
!!Note that this change required further modification of other parameters
!!above (c2, c3). If you want to remove this option, set c2 and c3 constants 
!!(above) back to NN2009 values (see commented out lines next to the
!!parameters above). This only removes the negative TKE problem
!!but does not necessarily improve performance - neutral impact.
 real(kind_phys), parameter :: CKmod=1.

!For calculating the b-f freq and Ri, either use the buoyancy flux functions
!(true) or use the direct calc of thlv with resolved and sgs clouds (false).  
 logical, parameter :: use_buoy=.false.

!>Use Ito et al. (2015, BLM) scale-aware (0: no, 1: yes). Note that this also has impacts
!!on the cloud PDF and mass-flux scheme, using LES-derived similarity function.
 real(kind_phys), parameter :: scaleaware=1.

!>Of the following the options, use one OR the other, not both.
!>Adding top-down diffusion driven by cloud-top radiative cooling
 integer, parameter :: bl_mynn_topdown = 0
!>Option to activate downdrafts, from Elynn Wu (0: deactive, 1: active)
 integer, parameter :: bl_mynn_edmf_dd = 0

!>Option to activate heating due to dissipation of TKE (1: active, 0: off)
 integer, parameter :: dheat_opt = 1

!Option to activate environmental subsidence in mass-flux scheme
 logical, parameter :: env_subs = .false.

!Option to switch flux-profile relationship for surface (from Puhales et al. 2020)
!0: use original Dyer-Hicks, 1: use Cheng-Brustaert and Blended COARE
 integer, parameter :: bl_mynn_stfunc = 1

!option to print out more stuff for debugging purposes
 logical, parameter :: debug_code = .false.
 integer, parameter :: idbg = 452 !specific i-point to write out
 integer, parameter :: jdbg = 272 !specific i-point to write out

!Used in WRF-ARW module_physics_init.F
 integer :: mynn_level


CONTAINS

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
             qnbca1             , ozone1            , p1                , &
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
             nchem              , kdvel             , ndvel             , &
             chem               , vdep              , frp               , &
             emis_ant_no        , mix_chem,enh_mix  , rrfs_sd           , &
             smoke_dbg          ,                                         &
             !higher-order moments
             tsq1               , qsq1              , cov1              , &
             !tendencies
             du1                , dv1               , dth1              , &
             dqv1               , dqc1              , dqi1              , &
             dqnc1              , dqni1             , dqs1              , &
             dqnwfa1            , dqnifa1           , dqnbca1           , &
             dozone1            , rthraten1         ,                     &
             !2d output
             pblh               , kpbl              ,                     &
             maxwidth           , maxmf             , ztop_plume        , &
             !tke budget arrays
             dqke1              , qwt1              , qshear1           , &
             qbuoy1             , qdiss1            ,                     &
             !subgrid clouds
             qc_bl1             , qi_bl1            , cldfra_bl1        , &
             !namelist configurations option
             bl_mynn_tkeadvect  , tke_budget        , bl_mynn_cloudpdf  , &
             bl_mynn_mixlength  , icloud_bl         , closure           , &
             bl_mynn_edmf       , bl_mynn_edmf_mom  , bl_mynn_edmf_tke  , &
             bl_mynn_mixscalars , bl_mynn_output    , bl_mynn_cloudmix  , &
             bl_mynn_mixqt      ,                                         &
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

 integer, intent(in) :: initflag
 logical, intent(in) :: restart,cycling
 integer, intent(in) :: tke_budget
 integer, intent(in) :: bl_mynn_cloudpdf
 integer, intent(in) :: bl_mynn_mixlength
 integer, intent(in) :: bl_mynn_edmf
 logical, intent(in) :: bl_mynn_tkeadvect
 integer, intent(in) :: bl_mynn_edmf_mom
 integer, intent(in) :: bl_mynn_edmf_tke
 integer, intent(in) :: bl_mynn_mixscalars
 integer, intent(in) :: bl_mynn_output
 integer, intent(in) :: bl_mynn_cloudmix
 integer, intent(in) :: bl_mynn_mixqt
 integer, intent(in) :: icloud_bl
 real(kind_phys), intent(in) :: closure

 logical, intent(in) :: FLAG_QI,FLAG_QNI,FLAG_QC,FLAG_QNC,&
                        FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA, &
                        FLAG_OZONE,FLAG_QS

 logical, intent(in) :: mix_chem,enh_mix,rrfs_sd,smoke_dbg

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
 real(kind_phys), intent(inout) :: pblh
 real(kind_phys), intent(inout) :: maxmf,maxwidth,ztop_plume
 integer,         intent(in)    :: i,j
 integer,         intent(inout) :: kpbl
!local
 real(kind_phys) :: psig_bl,psig_shcu,rmol
 integer :: imd,jmd
 integer :: k,kproblem

!column variables (all end with a "1")
 real(kind_phys), dimension(kts:kte), intent(in)      ::            &
       dz1,u1,v1,w1,th1,p1,ex1,rho1,tk1,rthraten1
 real(kind_phys), dimension(kts:kte), intent(inout)   ::            &
       sqv1,sqc1,sqi1,sqs1,qni1,qnc1,qnwfa1,qnifa1,qnbca1,ozone1,   &
       qke1,tsq1,qsq1,cov1,qke_adv1,                                &
       sh1,sm1,el1,                                                 & !interface, but kte+1 not included
       du1,dv1,dth1,dqv1,dqc1,dqi1,dqs1,                            &
       dqni1,dqnc1,dqnwfa1,dqnifa1,dqnbca1,dozone1,                 &
       qc_bl1,qi_bl1,cldfra_bl1,edmf_a1,edmf_w1,                    &
       edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1,                       &
       sub_thl1,sub_sqv1,det_thl1,det_sqv1
 real(kind_phys), dimension(kts:kte), intent(inout)   ::            &
       qwt1,qshear1,qbuoy1,qdiss1,dqke1
 real(kind_phys), dimension(kts:kte), intent(out)     ::            & !interface
       kh1,km1
!local
 real(kind_phys), dimension(kts:kte)                  ::            &
       qc_bl1_old,qi_bl1_old,cldfra_bl1_old,dummy1,dummy2,          &
       diss_heat1,                                                  &
       thl1,thv1,thlv1,qv1,qc1,qi1,qs1,sqw1,                        &
       dfm1, dfh1, dfq1, tcd1, qcd1,                                &
       pdk1, pdt1, pdq1, pdc1,                                      &
       vt1, vq1, sgm1, kzero1

! smoke/chemical arrays
 integer, intent(in) ::   nchem, kdvel, ndvel
 real(kind_phys), dimension(kts:kte,nchem), intent(inout) :: chem
 real(kind_phys), dimension(ndvel), intent(in)    :: vdep
 real(kind_phys),                   intent(in)    :: frp,emis_ant_no
 real(kind_phys), dimension(kts:kte+1,nchem)      :: s_awchem1
 integer :: ic

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
       s_awqv1,s_awqc1,s_awu1,s_awv1,s_awqke1,                      &
       s_awqnc1,s_awqni1,s_awqnwfa1,s_awqnifa1,                     &
       s_awqnbca1
 real(kind_phys), dimension(kts:kte+1) ::                           & !interface
       sd_aw1,sd_awthl1,sd_awqt1,                                   &
       sd_awqv1,sd_awqc1,sd_awu1,sd_awv1,sd_awqke1,                 &
       sd_awqi1,sd_awqnc1,sd_awqni1,                                &
       sd_awqnwfa1,sd_awqnifa1
 real(kind_phys), dimension(kts:kte+1) :: zw1              !interface
 real(kind_phys) :: cpm,sqcg,flt,fltv,flq,flqv,flqc,                &
       pmz,phh,exnerg,zet,phi_m,                                    &
       afk,abk,ts_decay, qc_bl2, qi_bl2,                            &
       th_sfc,wsp
 integer:: ktop_plume

!top-down diffusion
 real(kind_phys) :: maxKHtopdown
 real(kind_phys), dimension(kts:kte) :: KHtopdown
!mass flux tke production
 real(kind_phys), dimension(kts:kte) :: TKEprod_dn,TKEprod_up

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
              (ps-p1(kts)) < zero) then
             kproblem = k
             problem = .true.
             print*,"Incoming problem at: i=",i," j=",j," k=",k
             print*," QFX=",qfx," HFX=",hfx
             print*," wsp=",wsp," T=",tk1(k)
             print*," qv=",sqv1(k)," qc=",sqc1(k)
             print*," u*=",ust," wspd=",wspd
             print*," xland=",xland," ts=",ts
             print*," z/L=",half*dz1(1)*rmol," ps=",ps,"delp1=",ps-p1(kts)
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
          print*,"====p:",p1(max(kproblem-3,1):min(kproblem+3,kte))
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
    maxKHtopdown =zero
    kzero1       =zero

    ! DH* CHECK HOW MUCH OF THIS INIT IF-BLOCK IS ACTUALLY NEEDED FOR RESTARTS
!> - Within the MYNN-EDMF, there is a dependecy check for the first time step,
!! If true, a three-dimensional initialization loop is entered. Within this loop,
!! several arrays are initialized and k-oriented (vertical) subroutines are called 
!! at every i and j point, corresponding to the x- and y- directions, respectively.  
    IF (initflag > 0 .and. .not.restart) THEN

       !Test to see if we want to initialize qke
       IF ( (restart .or. cycling)) THEN
          IF (MAXVAL(qke1(:)) < 0.0002) THEN
             INITIALIZE_QKE = .TRUE.
             !print*,"QKE is too small, must initialize"
          ELSE
             INITIALIZE_QKE = .FALSE.
             !print*,"Using background QKE, will not initialize"
          ENDIF
       ELSE ! not cycling or restarting:
          INITIALIZE_QKE = .TRUE.
          !print*,"not restart nor cycling, must initialize QKE"
       ENDIF
 
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
          thv1(k)=th1(k)*(one+p608*sqv1(k))
          !keep snow out for now - increases ceiling bias
          sqw1(k)=sqv1(k)+sqc1(k)+sqi1(k)!+sqs1(k)
          thl1(k)=th1(k) - xlvcp/ex1(k)*sqc1(k) &
               &         - xlscp/ex1(k)*(sqi1(k))!+sqs1(k))
          !Use form from Tripoli and Cotton (1981) with their
          !suggested min temperature to improve accuracy.
          !thl1(k)=th1(k)*(1.- xlvcp/MAX(tk1(k),TKmin)*sqc1(k) &
          !    &             - xlscp/MAX(tk1(k),TKmin)*sqi1(k))
          thlv1(k)=thl1(k)*(one+p608*sqv1(k))
          zw1(k+1)=zw1(k)+dz1(k)
       enddo

       if (INITIALIZE_QKE) then
          !Initialize tke for initial PBLH calc only - using
          !simple PBLH form of Koracin and Berkowicz (1988, BLM)
          !to linearly taper off tke towards top of PBL.
          do k=kts,kte
             qke1(k)=5.*ust * MAX((ust*700. - zw1(k))/(MAX(ust,0.01)*700.), 0.01)
          enddo
       endif


!>  - Call get_pblh() to calculate hybrid (\f$\theta_{v}-TKE\f$) PBL height.
       CALL GET_PBLH(KTS,KTE,PBLH,thv1,qke1,zw1,dz1,xland,kpbl)
             
!>  - Call scale_aware() to calculate similarity functions for scale-adaptive control
!! (\f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$).
       IF (scaleaware > zero) THEN
          CALL SCALE_AWARE(dx,PBLH,Psig_bl,Psig_shcu)
       ELSE
          Psig_bl   = one
          Psig_shcu = one
       ENDIF

       ! DH* CHECK IF WE CAN DO WITHOUT CALLING THIS ROUTINE FOR RESTARTS
!>  - Call mym_initialize() to initializes the mixing length, TKE, \f$\theta^{'2}\f$,
!! \f$q^{'2}\f$, and \f$\theta^{'}q^{'}\f$. These variables are calculated after 
!! obtaining prerequisite variables by calling the following subroutines from 
!! within mym_initialize(): mym_level2() and mym_length().
       CALL mym_initialize (                & 
            &kts,kte,xland,                 &
            &dz1, dx, zw1,                  &
            &u1, v1, thl1, sqv1,            &
            &PBLH, th1, thv1, thlv1,        &
            &sh1, sm1,                      &
            &ust, rmol,                     &
            &el1, qke1, tsq1, qsq1, cov1,   &
            &Psig_bl, cldfra_bl1,           &
            &bl_mynn_mixlength,             &
            &edmf_w1,edmf_a1,               &
            &edmf_w_dd1,edmf_a_dd1,         &
            &INITIALIZE_QKE,                &
            &spp_pbl,pattern_spp_pbl1       )

       IF (.not.restart) THEN
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
    if (icloud_bl > 0) then
       cldfra_bl1_old(kts:kte)=cldfra_bl1(kts:kte)
       qc_bl1_old(kts:kte)    =qc_bl1(kts:kte)
       qi_bl1_old(kts:kte)    =qi_bl1(kts:kte)
    else
       cldfra_bl1    =zero
       qc_bl1        =zero
       qi_bl1        =zero
       cldfra_bl1_old=zero
       qc_bl1_old    =zero
       qi_bl1_old    =zero
    endif
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
    s_aw1      =zero
    s_awthl1   =zero
    s_awqt1    =zero
    s_awqv1    =zero
    s_awqc1    =zero
    s_awu1     =zero
    s_awv1     =zero
    s_awqke1   =zero
    s_awqnc1   =zero
    s_awqni1   =zero
    s_awqnwfa1 =zero
    s_awqnifa1 =zero
    s_awqnbca1 =zero
    s_awchem1(kts:kte+1,1:nchem) = zero
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
       thv1(k)=th1(k)*(one+p608*sqv1(k))
    enddo ! end k

!>  - Call get_pblh() to calculate the hybrid \f$\theta_{v}-TKE\f$
!! PBL height diagnostic.
    CALL GET_PBLH(KTS,KTE,PBLH,thv1,qke1,zw1,dz1,xland,KPBL)

!>  - Call scale_aware() to calculate the similarity functions,
!! \f$P_{\sigma-PBL}\f$ and \f$P_{\sigma-shcu}\f$, to control 
!! the scale-adaptive behaviour for the local and nonlocal 
!! components, respectively.
    if (scaleaware > 0.) then
       call SCALE_AWARE(dx,PBLH,Psig_bl,Psig_shcu)
    else
       Psig_bl=one
       Psig_shcu=one
    endif

    sqcg= zero   !ill-defined variable; qcg has been removed
    cpm=cp*(one + 0.84*max(sqv1(kts),1e-8))
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
 if (i==idbg .and. j==jdbg) then
    print*,"===incoming forcing in mynn driver:"
    print*,"fltv=",fltv,"flt=",flt,"hfx=",hfx
    print*,"flqv=",flqv,"flq=",flq,"qfx=",qfx
    print*,"cpm=",cpm,"rho=",rho1(kts),"ts=",ts
    print*,"cp=",cp,"qv=",sqv1(kts),"th1=",th1(kts)
 endif
    ! Update 1/L using updated sfc heat flux and friction velocity
    rmol = -karman*gtr*fltv/max(ust**3,1.0e-6_kind_phys)
    zet = half*dz1(kts)*rmol
    zet = max(zet, -20._kind_phys)
    zet = min(zet,  20._kind_phys)
    !if(i.eq.idbg)print*,"updated z/L=",zet
    if (bl_mynn_stfunc == 0) then
       !Original Kansas-type stability functions
       if ( zet >= zero ) then
          pmz = one + (cphm_st-one) * zet
          phh = one +  cphh_st      * zet
       else
          pmz = one/    (one-cphm_unst*zet)**0.25 - zet
          phh = one/sqrt(one-cphh_unst*zet)
       end if
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
         &p1,ex1,tsq1,qsq1,cov1,                      &
         &sh1,el1,bl_mynn_cloudpdf,                   &
         &qc_bl1,qi_bl1,cldfra_bl1,                   &
         &pblh,hfx,                                   &
         &vt1, vq1, th1, sgm1,                        &
         &spp_pbl, pattern_spp_pbl1                   )

!>  - Add TKE source driven by cloud top cooling
!!  Calculate the buoyancy production of TKE from cloud-top cooling when
!! \p bl_mynn_topdown =1.
    if (bl_mynn_topdown.eq.1) then
       call topdown_cloudrad(kts,kte,dz1,zw1,fltv,     &
            &xland,kpbl,PBLH,                          &
            &sqc1,sqi1,sqw1,thl1,th1,ex1,p1,rho1,thv1, &
            &cldfra_bl1,rthraten1,                     &
            &maxKHtopdown,KHtopdown,TKEprod_dn         )
    else
       maxKHtopdown = zero
       KHtopdown    = zero
       TKEprod_dn   = zero
    endif

    if (bl_mynn_edmf > 0) then
       !PRINT*,"Calling DMP Mass-Flux"
       call DMP_mf(i,j,                               &
            &kts,kte,delt,zw1,dz1,p1,rho1,            &
            &bl_mynn_edmf_mom,                        &
            &bl_mynn_edmf_tke,                        &
            &bl_mynn_mixscalars,                      &
            &u1,v1,w1,th1,thl1,thv1,tk1,              &
            &sqw1,sqv1,sqc1,qke1,                     &
            &qnc1,qni1,qnwfa1,qnifa1,qnbca1,          &
            &ex1,vt1,vq1,sgm1,                        &
            &ust,flt,fltv,flq,flqv,                   &
            &pblh,kpbl,dx,                            &
            &xland,th_sfc,                            &
            ! now outputs - tendencies
            ! &,dth1mf,dqv1mf,dqc1mf,du1mf,dv1mf         &
            ! outputs - updraft properties
            &edmf_a1,edmf_w1,edmf_qt1,                &
            &edmf_thl1,edmf_ent1,edmf_qc1,            &
            ! for the solver
            &s_aw1,s_awthl1,s_awqt1,                  &
            &s_awqv1,s_awqc1,                         &
            &s_awu1,s_awv1,s_awqke1,                  &
            &s_awqnc1,s_awqni1,                       &
            &s_awqnwfa1,s_awqnifa1,s_awqnbca1,        &
            &sub_thl1,sub_sqv1,                       &
            &sub_u1,sub_v1,                           &
            &det_thl1,det_sqv1,det_sqc1,              &
            &det_u1,det_v1,                           &
            ! chem/smoke mixing
            &nchem,chem,s_awchem1,                    &
            &mix_chem,                                &
            &qc_bl1,cldfra_bl1,                       &
            &qc_bl1_old,cldfra_bl1_old,               &
            &FLAG_QC,FLAG_QI,                         &
            &FLAG_QNC,FLAG_QNI,                       &
            &FLAG_QNWFA,FLAG_QNIFA,FLAG_QNBCA,        &
            &Psig_shcu,                               &
            &maxwidth,ktop_plume,                     &
            &maxmf,ztop_plume,                        &
            &spp_pbl,pattern_spp_pbl1,                &
            &TKEprod_up,el1                           )
    else
       TKEprod_up = zero
    endif

    if (bl_mynn_edmf_dd == 1) then
       call ddmp_mf(kts,kte,delt,dx,zw1,dz1,p1,       &
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
            &rthraten1                                )
    else
       tkeprod_dn = zero
    endif

    !Capability to substep the eddy-diffusivity portion
    !do nsub = 1,2
    delt2 = delt !*0.5    !only works if topdown=0

    !thlv with updated subgrid clouds
    do k=kts,kte
       thlv1(k)=(th1(k) - xlvcp/ex1(k)*max(qc_bl1(k),sqc1(k))          &
                &       - xlscp/ex1(k)*max(qi_bl1(k),sqi1(k)+sqs1(k))) &
                &       * (one+p608*sqv1(k))
    enddo

    call mym_turbulence(                                 &
            &kts,kte,xland,closure,                      &
            &dz1, dx, zw1,                               &
            &u1, v1, thl1, thv1, thlv1,                  &
            &sqc1, sqw1,                                 &
            &qke1, tsq1, qsq1, cov1,                     &
            &vt1, vq1,                                   &
            &rmol, flt, fltv, flq,                       &
            &pblh, th1,                                  &
            &sh1,sm1,el1,                                &
            &Dfm1,Dfh1,Dfq1,                             &
            &Tcd1,Qcd1,Pdk1,                             &
            &Pdt1,Pdq1,Pdc1,                             &
            &qWT1,qSHEAR1,qBUOY1,qDISS1,                 &
            &tke_budget,                                 &
            &Psig_bl,Psig_shcu,                          &
            &cldfra_bl1,bl_mynn_mixlength,               &
            &edmf_w1,edmf_a1,                            &
            &edmf_w_dd1,edmf_a_dd1,                      &
            &TKEprod_dn,TKEprod_up,                      &
            &spp_pbl,pattern_spp_pbl1                    )

!>  - Call mym_predict() to solve TKE and 
!! \f$\theta^{'2}, q^{'2}, and \theta^{'}q^{'}\f$
!! for the following time step.
    call mym_predict(kts,kte,closure,                    &
            &delt2, dz1,                                 &
            &ust, flt, flq, pmz, phh,                    &
            &el1, dfq1, rho1, pdk1, pdt1, pdq1, pdc1,    &
            &qke1, tsq1, qsq1, cov1,                     &
            &s_aw1, s_awqke1, bl_mynn_edmf_tke,          &
            &qWT1, qDISS1, tke_budget                    )

    if (dheat_opt > 0) then
       do k=kts,kte-1
          ! Set max dissipative heating rate to 7.2 K per hour
          diss_heat1(k) = MIN(MAX(1.0*(qke1(k)**1.5)/(b1*MAX(half*(el1(k)+el1(k+1)),one))/cp, 0.0),0.002)
          ! Limit heating above 100 mb:
          diss_heat1(k) = diss_heat1(k) * exp(-10000./MAX(p1(k),one)) 
       enddo
       diss_heat1(kte) = 0.
    else
       diss_heat1 = 0.
    endif

!>  - Call mynn_tendencies() to solve for tendencies of 
!! \f$U, V, \theta, q_{v}, q_{c}, and q_{i}\f$.
    call mynn_tendencies(kts,kte,i,                      &
            &delt, dz1, rho1,                            &
            &u1, v1, th1, tk1, qv1,                      &
            &qc1, qi1, kzero1, qnc1, qni1,               & !kzero replaces qs1 - not mixing snow
            &ps, p1, ex1, thl1,                          &
            &sqv1, sqc1, sqi1, kzero1, sqw1,             & !kzero replaces sqs - not mixing snow
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
            &bl_mynn_mixscalars                          )

    if ( mix_chem ) then
       if ( rrfs_sd ) then 
          call mynn_mix_chem(kts,kte,i,                  &
               &delt, dz1, pblh,                         &
               &nchem, kdvel, ndvel,                     &
               &chem, vdep,                              &
               &rho1, flt,                               &
               &tcd1, qcd1,                              &
               &dfh1,                                    &
               &s_aw1,s_awchem1,                         &
               &emis_ant_no,                             &
               &frp, rrfs_sd,                            &
               &enh_mix, smoke_dbg                       )
       else
          call mynn_mix_chem(kts,kte,i,                  &
               &delt, dz1, pblh,                         &
               &nchem, kdvel, ndvel,                     &
               &chem, vdep,                              &
               &rho1, flt,                               &
               &tcd1, qcd1,                              &
               &dfh1,                                    &
               &s_aw1,s_awchem1,                         &
               &zero,                                    &
               &zero, rrfs_sd,                           &
               &enh_mix, smoke_dbg                       )
       endif
       do ic = 1,nchem
          do k = kts,kte
             chem(k,ic) = max(1.e-12, chem(k,ic))
          enddo
       enddo
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
          dummy1(k) = half*(qSHEAR1(k)+qSHEAR1(k+1)) !!! unstaggering in z
          dummy2(k) = half*(qBUOY1(k)+qBUOY1(k+1)) !!! unstaggering in z
          dqke1(k)  = half*(qke1(k)-dqke1(k))/delt
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
          IF (icloud_bl > 0) then
             IF (cldfra_bl1(k) < zero .OR. cldfra_bl1(k)> one)THEN
                PRINT*,"SUSPICIOUS VALUES: CLDFRA_BL=",cldfra_bl1(k),&
                                             " qc_bl=",qc_bl1(k)
             endif
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
       &            kts,kte,xland,                            &
       &            dz, dx, zw,                               &
       &            u, v, thl, qw,                            &
!       &            ust, rmo, pmz, phh, flt, flq,             &
       &            zi, theta, thv, thlv, sh, sm,             &
       &            ust, rmol, el,                            &
       &            Qke, Tsq, Qsq, Cov, Psig_bl, cldfra_bl1,  &
       &            bl_mynn_mixlength,                        &
       &            edmf_w1,edmf_a1,                          &
       &            edmf_w_dd1,edmf_a_dd1,                    &
       &            INITIALIZE_QKE,                           &
       &            spp_pbl,pattern_spp_pbl1                  )
!
!-------------------------------------------------------------------

    integer, intent(in)           :: kts,kte
    integer, intent(in)           :: bl_mynn_mixlength
    logical, intent(in)           :: INITIALIZE_QKE
!    real(kind_phys), intent(in)   :: ust, rmol, pmz, phh, flt, flq
    real(kind_phys), intent(in)   :: rmol, Psig_bl, xland
    real(kind_phys), intent(in)   :: dx, ust, zi
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), dimension(kts:kte),   intent(in) :: u,v,thl,&
         &thlv,qw,cldfra_bl1,edmf_w1,edmf_a1,edmf_w_dd1,edmf_a_dd1
    real(kind_phys), dimension(kts:kte),   intent(inout) :: tsq,qsq,cov
    real(kind_phys), dimension(kts:kte),   intent(inout) :: el,qke
    real(kind_phys), dimension(kts:kte) ::                       &
         &ql,pdk,pdt,pdq,pdc,dtl,dqw,dtv,                        &
         &gm,gh,sm,sh,qkw,vt,vq
    integer :: k,l,lmax
    real(kind_phys):: phm,vkz,elq,elv,b1l,b2l,pmz=1.,phh=1.,     &
         &flt=0.,fltv=0.,flq=0.,tmpq
    real(kind_phys), dimension(kts:kte) :: theta,thv
    real(kind_phys), dimension(kts:kte) :: pattern_spp_pbl1
    integer ::spp_pbl

!> - At first ql, vt and vq are set to zero.
    DO k = kts,kte
       ql(k) = zero
       vt(k) = zero
       vq(k) = zero
    END DO
!
!> - Call mym_level2() to calculate the stability functions at level 2.
    CALL mym_level2 ( kts,kte,                      &
         &            dz,                           &
         &            u, v, thl, thv, thlv, qw,     &
         &            ql, vt, vq,                   &
         &            dtl, dqw, dtv, gm, gh, sm, sh )
!
!   **  Preliminary setting  **

    el(kts) = zero
    IF (INITIALIZE_QKE) THEN
       !qke(kts) = ust**2 * ( b1*pmz )**(2.0/3.0)
       qke(kts) = 1.5 * ust**2 * ( b1*pmz )**(2.0/3.0)
       DO k = kts+1,kte
          !qke(k) = 0.0
          !linearly taper off towards top of pbl
          qke(k)=qke(kts)*MAX((ust*700. - zw(k))/(MAX(ust,0.01)*700.), 0.01)
       ENDDO
    ENDIF
!
    phm      = phh*b2 / ( b1*pmz )**(1.0/3.0)
    tsq(kts) = phm*( flt/ust )**2
    qsq(kts) = phm*( flq/ust )**2
    cov(kts) = phm*( flt/ust )*( flq/ust )
!
    DO k = kts+1,kte
       vkz = karman*zw(k)
       el(k) = vkz/( one + vkz/100.0 )
!       qke(k) = 0.0
!
       tsq(k) = 0.0
       qsq(k) = 0.0
       cov(k) = 0.0
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
            &            rmol, flt, fltv, flq,    &
            &            vt, vq,                  &
            &            u, v, qke,               &
            &            dtv,                     &
            &            el,                      &
            &            zi,theta,                &
            &            qkw,Psig_bl,cldfra_bl1,  &
            &            bl_mynn_mixlength,       &
            &            edmf_w1,edmf_a1,         &
            &            edmf_w_dd1,edmf_a_dd1    )
!
       DO k = kts+1,kte
          elq = el(k)*qkw(k)
          pdk(k) = elq*( sm(k)*gm(k) + &
               &         sh(k)*gh(k) )
          pdt(k) = elq*  sh(k)*dtl(k)**2
          pdq(k) = elq*  sh(k)*dqw(k)**2
          pdc(k) = elq*  sh(k)*dtl(k)*dqw(k)
       END DO
!
!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
       vkz = karman*half*dz(kts)
       elv = half*( el(kts+1)+el(kts) ) /  vkz
       IF (INITIALIZE_QKE)THEN 
          !qke(kts) = ust**2 * ( b1*pmz*elv    )**(2.0/3.0)
          qke(kts) = 1.0 * MAX(ust,0.02)**2 * ( b1*pmz*elv    )**(2.0/3.0) 
       ENDIF

       phm      = phh*b2 / ( b1*pmz/elv**2 )**(1.0/3.0)
       tsq(kts) = phm*( flt/ust )**2
       qsq(kts) = phm*( flq/ust )**2
       cov(kts) = phm*( flt/ust )*( flq/ust )

       DO k = kts+1,kte-1
          b1l = b1*0.25*( el(k+1)+el(k) )
          !tmpq=MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin)
          !add MIN to limit unreasonable QKE
          tmpq=MIN(MAX(b1l*( pdk(k+1)+pdk(k) ),qkemin),125.)
          !print*,'tmpq=',tmpq,pdk(k+1),pdk(k)
          IF (INITIALIZE_QKE)THEN
             qke(k) = tmpq**twothirds
          ENDIF

          IF ( qke(k) .LE. zero ) THEN
             b2l = 0.0
          ELSE
             b2l = b2*( b1l/b1 ) / SQRT( qke(k) )
          END IF

          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO

    END DO

!!    qke(kts)=qke(kts+1)
!!    tsq(kts)=tsq(kts+1)
!!    qsq(kts)=qsq(kts+1)
!!    cov(kts)=cov(kts+1)

    IF (INITIALIZE_QKE)THEN
       qke(kts)=0.5*(qke(kts)+qke(kts+1))
       qke(kte)=qke(kte-1)
    ENDIF
    tsq(kte)=tsq(kte-1)
    qsq(kte)=qsq(kte-1)
    cov(kte)=cov(kte-1)

!
!    RETURN

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
       &            dz,                           &
       &            u, v, thl, thv, thlv,         &
       &            qw, ql, vt, vq,               &
       &            dtl, dqw, dtv, gm, gh, sm, sh )
!
!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), dimension(kts:kte), intent(in)  :: dz
    real(kind_phys), dimension(kts:kte), intent(in)  :: u,v, &
         &thl,qw,ql,vt,vq,thv,thlv
    real(kind_phys), dimension(kts:kte), intent(out) ::      &
         &dtl,dqw,dtv,gm,gh,sm,sh

    integer :: k

    real(kind_phys):: rfc,f1,f2,rf1,rf2,smc,shc,             &
         &ri1,ri2,ri3,ri4,duz,dtz,dqz,vtt,vqq,dtq,dzk,       &
         &afk,abk,ri,rf

    real(kind_phys):: a2fac

!    ev  = 2.5e6
!    tv0 = 0.61*tref
!    tv1 = 1.61*tref
!    gtr = 9.81/tref
!
!intialize output
    dtl(kts)=0.0
    dqw(kts)=0.0
    dtv(kts)=0.0
    gm(kts)=0.0
    gh(kts)=0.0
    sm(kts)=0.0
    sh(kts)=0.0

    do k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = one -afk
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       dtz = ( thl(k)-thl(k-1) )/( dzk )
       dqz = ( qw(k)-qw(k-1) )/( dzk )
!
       vtt =  one +vt(k)*abk +vt(k-1)*afk  ! Beta-theta in NN09, Eq. 39
       vqq =  tv0 +vq(k)*abk +vq(k-1)*afk  ! Beta-q
       if (use_buoy) then
          !use the buoyancy flux functions
          dtq =  vtt*dtz +vqq*dqz
       else
          !alternatively, use theta-l-v with the SGS clouds
          dtq = ( thlv(k)-thlv(k-1) )/( dzk )
       endif
       dtl(k) =  dtz
       dqw(k) =  dqz
       dtv(k) =  dtq
!
       gm(k)  =  duz
       gh(k)  = -dtq*gtr
!
!   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )

       !a2fac is needed for the Canuto/Kitamura mod
       if (CKmod .eq. 1) then
          a2fac = 1./(1. + MAX(ri,0.0))
       else
          a2fac = 1.
       endif

       rfc = g1/( g1+g2 )
       f1  = b1*( g1-c1 ) +3.0*a2*a2fac *( one    -c2 )*( one-c5 ) &
    &                     +2.0*a1*( 3.0-2.0*c2 )
       f2  = b1*( g1+g2 ) -3.0*a1*( one    -c2 )
       rf1 = b1*( g1-c1 )/f1
       rf2 = b1*  g1     /f2
       smc = a1 /(a2*a2fac)*  f1/f2
       shc = 3.0*(a2*a2fac)*( g1+g2 )

       ri1 = 0.5/smc
       ri2 = rf1*smc
       ri3 = 4.0*rf2*smc -2.0*ri2
       ri4 = ri2**2

!   **  Flux Richardson number  **
       rf = MIN( ri1*( ri + ri2-SQRT(ri**2 - ri3*ri + ri4) ), rfc )
!
       sh(k) = shc*( rfc-rf )/( one-rf )
       sm(k) = smc*( rf1-rf )/( rf2-rf ) * sh(k)
    enddo
!
!    RETURN

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
!       elt(nx,ny)      : Length scale depending on the PBL depth    (m)
!       vsc(nx,ny)      : Velocity scale q_c                       (m/s)
!                         at first, used for computing elt
!
!     NOTE: the mixing lengths are meant to be calculated at the full-
!           sigmal levels (or interfaces beween the model layers).
!
!>\ingroup gsd_mynn_edmf
!! This subroutine calculates the mixing lengths.
  SUBROUTINE  mym_length (                     & 
    &            kts,kte,xland,                &
    &            dz, dx, zw,                   &
    &            rmol, flt, fltv, flq,         &
    &            vt, vq,                       &
    &            u1, v1, qke,                  &
    &            dtv,                          &
    &            el,                           &
    &            pblh, theta, qkw,             &
    &            Psig_bl, cldfra_bl1,          &
    &            bl_mynn_mixlength,            &
    &            edmf_w1,edmf_a1,              &
    &            edmf_w_dd1,edmf_a_dd1         )
    
!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in)   :: bl_mynn_mixlength
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in) :: rmol,flt,fltv,flq,Psig_bl,xland
    real(kind_phys), intent(in) :: dx,pblh
    real(kind_phys), dimension(kts:kte), intent(in)   :: u1,v1,  &
         &qke,vt,vq,cldfra_bl1,edmf_w1,edmf_a1,edmf_w_dd1,edmf_a_dd1
    real(kind_phys), dimension(kts:kte), intent(out)  :: qkw, el
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
    real(kind_phys), parameter :: minpblh     = 300.  !< min mixed-layer height
    real(kind_phys), parameter :: maxdz       = 750.  !< max (half) transition layer depth
                                     !! =0.3*2500 m PBLH, so the transition
                                     !! layer stops growing for PBLHs > 2.5 km.
    real(kind_phys), parameter :: mindz       = 300.  !< 300  !min (half) transition layer depth

    !SURFACE LAYER LENGTH SCALE MODS TO REDUCE IMPACT IN UPPER BOUNDARY LAYER
    real(kind_phys), parameter :: ZSLH        = 100.  !< Max height correlated to surface conditions (m)
    real(kind_phys), parameter :: CSL         = 2.    !< CSL = constant of proportionality to L O(1)
    real(kind_phys), parameter :: qkw_elb_min = 0.18

    integer :: i,j,k
    real(kind_phys):: afk,abk,zwk,zwk1,dzk,qdz,vflx,bv,tau_cloud,      &
           & wstar,elb,els,elf,el_stab,el_mf,el_stab_mf,elb_mf,elt_max,&
           & PBLH_PLUS_ENT,Uonset,Ugrid,wt_u1,wt_u2,el_les,qkw_mf
    real(kind_phys), parameter :: ctau = 1000. !constant for tau_cloud

!    tv0 = 0.61*tref
!    gtr = 9.81/tref

    SELECT CASE(bl_mynn_mixlength)

      CASE (0) ! ORIGINAL MYNN MIXING LENGTH + BouLac

        cns  = 2.7
        alp1 = 0.23
        alp2 = one
        alp3 = 5.0
        alp4 = 100.
        alp5 = 0.3

        ! Impose limits on the height integration for elt and the transition layer depth
        pblh2= min(10000.,zw(kte-2))  !originally integrated to model top, not just 10 km.
        h1   = max(0.3*pblh2,mindz)
        h1   = min(h1,maxdz)         ! 1/2 transition layer depth
        h2   = h1/2.0                ! 1/4 transition layer depth

        qkw(kts) = SQRT(MAX(qke(kts), qkemin))
        DO k = kts+1,kte
           afk = dz(k)/( dz(k)+dz(k-1) )
           abk = one -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk, qkemin))
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5        

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k   = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. pblh2+h1)
           dzk = half*( dz(k)+dz(k-1) )
           qdz = MAX( qkw(k)-qmin, 0.03 )*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt  = alp1*elt/vsc
        vflx = ( vt(kts)+one )*flt +( vq(kts)+tv0 )*flq
        vsc  = ( gtr*elt*MAX( vflx, zero ) )**onethird

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
              elb = 1.0e10
              elf = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmol .GT. 0.0 ) THEN
              els  = karman*zwk/(one+cns*MIN( zwk*rmol, zmax ))
           ELSE
              els  =  karman*zwk*( one - alp4* zwk*rmol )**0.2
           END IF

           !   ** HARMONC AVERGING OF MIXING LENGTH SCALES:
           !       el(k) =    MIN(elb/( elb/elt+elb/els+one ),elf)
           !       el(k) =    elb/( elb/elt+elb/els+one )

           wt=half*TANH((zwk - (pblh2+h1))/h2) + half

           el(k) = MIN(elb/( elb/elt+elb/els+one ),elf)

        END DO

      CASE (1) !NONLOCAL (using BouLac) FORM OF MIXING LENGTH

        !wt_u* are for hurricane tuning, meant to reduce diffusion in hurricanes
        ugrid = sqrt(u1(kts)**2 + v1(kts)**2)
        uonset= 20.
        wt_u1 = one - 0.2*min(1.0, max(zero, ugrid - uonset)/50.0) !reduce to 0.8
        wt_u2 = one - 0.4*min(1.0, max(zero, ugrid - uonset)/50.0) !reduce to 0.6
        cns   = 3.5
        alp1  = 0.23
        alp2  = 0.3
        alp3  = 2.5 * wt_u2 !taper off bouyancy enhancement in shear-driven pbls
        alp4  = 5.0
        alp5  = 0.3
        alp6  = 50.

        ! Impose limits on the height integration for elt and the transition layer depth
        pblh2= max(pblh,300.) !minpblh)
        h1   = max(0.3*pblh2,300.)
        h1   = min(h1,600.)          ! 1/2 transition layer depth
        h2   = h1/2.0                ! 1/4 transition layer depth

        qtke(kts)  =max(half*qke(kts), half*qkemin) !tke at full sigma levels
        thetaw(kts)=theta(kts)            !theta at full-sigma levels
        qkw(kts)   =sqrt(max(qke(kts), qkemin))

        DO k = kts+1,kte
           afk      = dz(k)/( dz(k)+dz(k-1) )
           abk      = one -afk
           qkw(k)   = sqrt(max(qke(k)*abk+qke(k-1)*afk, qkemin))
           qtke(k)  = max(half*(qkw(k)**2), 0.005) ! q -> TKE
           thetaw(k)= theta(k)*abk + theta(k-1)*afk
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        k   = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. pblh2+h1)
           dzk = half*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.01 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        if ((xland-1.5).GE.zero) then !hurricane tuning, over water only
           elt_max=350.+100.*min(one, max(zero, ugrid - 50.0)/25.0)
        else
           elt_max=400.
        endif
        elt = MIN( MAX( alp1*elt/vsc, 8.), elt_max)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+one )*flt + ( vq(kts)+tv0 )*flq
        vflx= fltv
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**onethird

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)         !full-sigma levels

        ! COMPUTE BouLac mixing length
        CALL boulac_length(kts,kte,zw,dz,qtke,thetaw,elBLmin,elBLavg)

        DO k = kts+1,kte
           zwk    = zw(k)          !full-sigma levels
           qkw_mf = max(edmf_a1(k-1)*edmf_w1(k-1),                &
                  & abs(edmf_a_dd1(k-1)*edmf_w_dd1(k-1)))

           !   **  Length scale limited by the buoyancy effect  **
           IF ( dtv(k) .GT. 0.0 ) THEN
              bv     = max( sqrt( gtr*dtv(k) ), 0.001)
              elb_mf = alp6*qkw_mf/bv
              elb_mf = elb_mf / (1. + (elb_mf/100.))
              elb    = alp2*(max(qkw(k),qkw_elb_min))/bv    &
                     &  *( one + alp3*SQRT( vsc/(bv*elt) ) )
              elb    = max(elb, elb_mf)
              elb    = MIN(elb, zwk)
              elf    = one * max(qkw(k), qkw_elb_min)/bv
              elf    = max(elf, elb_mf)
              elBLavg(k) = MAX(elBLavg(k), elb_mf)
           ELSE
              elb    = 1.0e10
              elf    = elb
           ENDIF

           !   **  Length scale in the surface layer  **
           IF ( rmol .GT. 0.0 ) THEN
              els  = karman*zwk/(one+cns*MIN( zwk*rmol, zmax ))
           ELSE
              els  = karman*zwk*( one - alp4* zwk*rmol)**0.2
           END IF

           !   ** NOW BLEND THE MIXING LENGTH SCALES:
           wt=half*TANH((zwk - (pblh2+h1))/h2) + half

           !add blending to use BouLac mixing length in free atmos;
           !defined relative to the PBLH (pblh) + transition layer (h1)
           !el(k) = MIN(elb/( elb/elt+elb/els+one ),elf)
           !try squared-blending - but take out elb (makes it underdiffusive)
           !el(k) = SQRT( els**2/(1. + (els**2/elt**2) +(els**2/elb**2)))
           el(k) = sqrt( els**2/(1. + (els**2/elt**2)))
           el(k) = min(el(k), elb)
           el(k) = min(el(k), elf)  !elf can be smaller than elb in upper pbl
           if ((xland-1.5).GE.zero) then !hurricane tuning, over water only
              el(k)=el(k)*wt_u1
           endif
           el(k) = el(k)*(1.-wt) + alp5*elBLavg(k)*wt

           !if (el(k) > 1000.) then
           !   print*,"big ML at k=",k," el=",el(k),"elb=",elb," elBL=",elBLavg(k), &
           !   " elt=",elt," els=",els," N=",bv," qtke=",qtke(k), &
           !   " pblh=",pblh2," zwk=",zwk," wt=",wt
           !endif

           ! include scale-awareness, except for original MYNN
           el(k) = el(k)*Psig_bl

        END DO

     CASE (2) !Local (mostly) mixing length formulation

        Uonset = 3.5 + dz(kts)*0.1
        Ugrid  = sqrt(u1(kts)**2 + v1(kts)**2)
        cns  = 3.5 !JOE-test  * (one - MIN(MAX(Ugrid - Uonset, 0.0)/10.0, 1.0))
        alp1 = 0.22
        alp2 = 0.30
        alp3 = 2.5
        alp4 = 5.0
        alp5 = alp2 !like alp2, but for free atmosphere
        alp6 = 50.0 !used for MF mixing length

        ! Impose limits on the height integration for elt and the transition layer depth
        !pblh2=MAX(pblh,minpblh)
        pblh2=MAX(pblh,    300.)
        !h1=MAX(0.3*pblh2,mindz)
        !h1=MIN(h1,maxdz)         ! 1/2 transition layer depth
        h1=MAX(0.3*pblh2,300.)
        h1=MIN(h1,600.)
        h2=h1*half                ! 1/4 transition layer depth

        qtke(kts)=MAX(half*qke(kts), half*qkemin) !tke at full sigma levels
        qkw(kts) = SQRT(MAX(qke(kts), qkemin))

        DO k = kts+1,kte
           afk    = dz(k)/( dz(k)+dz(k-1) )
           abk    = one -afk
           qkw(k) = SQRT(MAX(qke(k)*abk+qke(k-1)*afk, qkemin))
           qtke(k)= half*qkw(k)**2  ! qkw -> TKE
        END DO

        elt = 1.0e-5
        vsc = 1.0e-5

        !   **  Strictly, zwk*h(i,j) -> ( zwk*h(i,j)+z0 )  **
        PBLH_PLUS_ENT = MAX(pblh+h1, 100.)
        k = kts+1
        zwk = zw(k)
        DO WHILE (zwk .LE. PBLH_PLUS_ENT)
           dzk = half*( dz(k)+dz(k-1) )
           qdz = min(max( qkw(k)-qmin, 0.03 ), 30.0)*dzk
           elt = elt +qdz*zwk
           vsc = vsc +qdz
           k   = k+1
           zwk = zw(k)
        END DO

        elt = MIN( MAX(alp1*elt/vsc, 10.), 400.)
        !avoid use of buoyancy flux functions which are ill-defined at the surface
        !vflx = ( vt(kts)+one )*flt +( vq(kts)+tv0 )*flq
        vflx = fltv
        vsc = ( gtr*elt*MAX( vflx, 0.0 ) )**onethird

        !   **  Strictly, el(i,j,1) is not zero.  **
        el(kts) = 0.0
        zwk1    = zw(kts+1)

        DO k = kts+1,kte
           zwk = zw(k)              !full-sigma levels
           dzk = 0.5*( dz(k)+dz(k-1) )
           cldavg = 0.5*(cldfra_bl1(k-1)+cldfra_bl1(k))
           qkw_mf = max(edmf_a1(k-1)*edmf_w1(k-1),               &
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

              !tau_cloud = MIN(MAX(0.5*pblh/((gtr*pblh*MAX(vflx,1.0e-4))**onethird),30.),150.)
              wstar = 1.25*(gtr*pblh*MAX(vflx,1.0e-4))**onethird
              tau_cloud = MIN(MAX(ctau * wstar/grav, 30.), 150.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt=half*TANH((zwk - (pblh2+h1))/h2) + half
              tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              elf = MIN(MAX(tau_cloud*SQRT(MIN(qtke(k),40.)), &
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
              !tau_cloud = MIN(MAX(0.5*pblh/((gtr*pblh*MAX(vflx,1.0e-4))**onethird),50.),150.)
              wstar     = 1.25*(gtr*pblh*MAX(vflx,1.0e-4))**onethird
              tau_cloud = MIN(MAX(ctau * wstar/grav, 50.), 200.)
              !minimize influence of surface heat flux on tau far away from the PBLH.
              wt        = half*TANH((zwk - (pblh2+h1))/h2) + half
              !tau_cloud = tau_cloud*(1.-wt) + 50.*wt
              tau_cloud = tau_cloud*(1.-wt) + MAX(100.,dzk*0.25)*wt

              elb       = MIN(tau_cloud*SQRT(MIN(qtke(k),40.)), zwk)
              !elf = elb
              elf       = elb !/(1. + (elb/800.))  !bound free-atmos mixing length to < 800 m.
              elb_mf    = elb
         END IF
         elf    = elf/(1. + (elf/800.))  !bound free-atmos mixing length to < 800 m.
         elb_mf = MAX(elb_mf, 0.01) !to avoid divide-by-zero below

         !   **  Length scale in the surface layer  **
         IF ( rmol .GT. 0.0 ) THEN
            els  = karman*zwk/(one+cns*MIN( zwk*rmol, zmax ))
         ELSE
            els  =  karman*zwk*( one - alp4* zwk*rmol)**0.2
         END IF

         !   ** NOW BLEND THE MIXING LENGTH SCALES:
         wt=half*TANH((zwk - (pblh2+h1))/h2) + half

         !try squared-blending
         el(k) = SQRT( els**2/(1. + (els**2/elt**2) +(els**2/elb_mf**2)))
         el(k) = el(k)*(1.-wt) + elf*wt

         ! include scale-awareness. For now, use simple asymptotic kz -> 12 m (should be ~dz).
         el_les= MIN(els/(1. + (els/12.)), elb_mf)
         el(k) = el(k)*Psig_bl + (1.-Psig_bl)*el_les

       END DO

    END SELECT


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
     zup=0.
     dlu=zw(kte+1)-zw(k)-dz(k)*0.5
     zzz=0.
     zup_inf=0.
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
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt*0.5 ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer k to izz+1
              !print*,"  PE=",zup," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zup .and. qtke(k).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(k)) + &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(k)))**2 + &
                      &       2.*bbb*beta*(qtke(k)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k))then
                       tl=(qtke(k)-zup_inf)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=0.
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
     zdo=0.
     zdo_sup=0.
     dld=zw(k)
     zzz=0.

     !print*,"FINDING Ddown, k=",k," zwk=",zw(k)
     if (k .gt. kts) then  !cant integrate downwards from lowest level

        found = 0
        izz=k
        DO WHILE (found .EQ. 0)

           if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(k)*dzt
              !print*,"  ",k,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt*0.5
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(k)," z=",zw(izz)
              if (qtke(k).lt.zdo .and. qtke(k).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(k))+ &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(k)))**2 + &
                      &       2.*bbb*beta*(qtke(k)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(k)) then
                       tl=(qtke(k)-zdo_sup)/(beta*(theta(izz)-theta(k)))
                    else
                       tl=0.
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
     dlu=MAX(0.1,MIN(dlu,1000.))
     dld=MAX(0.1,MIN(dld,1000.))
     lb2 = sqrt(dlu*dld)    !average - biased towards smallest
     !lb2 = 0.5*(dlu+dld)   !average

     if (k .eq. kte) then
        lb1 = 0.
        lb2 = 0.
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
        zup=0.
        dlu(iz)=zw(kte+1)-zw(iz)-dz(iz)*0.5
        zzz=0.
        zup_inf=0.
        beta=gtr           !Buoyancy coefficient (g/tref)

        !print*,"FINDING Dup, k=",iz," zw=",zw(iz)

        if (iz .lt. kte) then      !cant integrate upwards from highest level

          found = 0
          izz=iz
          DO WHILE (found .EQ. 0)

            if (izz .lt. kte) then
              dzt=dz(izz)                    ! layer depth above
              zup=zup-beta*theta(iz)*dzt     ! initial PE the parcel has at iz
              !print*,"  ",iz,izz,theta(izz),dz(izz)
              zup=zup+beta*(theta(izz+1)+theta(izz))*dzt*0.5 ! PE gained by lifting a parcel to izz+1
              zzz=zzz+dzt                   ! depth of layer iz to izz+1
              !print*,"  PE=",zup," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zup .and. qtke(iz).ge.zup_inf) then
                 bbb=(theta(izz+1)-theta(izz))/dzt
                 if (bbb .ne. 0.) then
                    !fractional distance up into the layer where TKE becomes < PE
                    tl=(-beta*(theta(izz)-theta(iz)) + &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(iz)))**2 + &
                      &       2.*bbb*beta*(qtke(iz)-zup_inf))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz))then
                       tl=(qtke(iz)-zup_inf)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=0.
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
        zdo=0.
        zdo_sup=0.
        dld(iz)=zw(iz)
        zzz=0.

        !print*,"FINDING Ddown, k=",iz," zwk=",zw(iz)
        if (iz .gt. kts) then  !cant integrate downwards from lowest level

          found = 0
          izz=iz       
          DO WHILE (found .EQ. 0) 

            if (izz .gt. kts) then
              dzt=dz(izz-1)
              zdo=zdo+beta*theta(iz)*dzt
              !print*,"  ",iz,izz,theta(izz),dz(izz-1)
              zdo=zdo-beta*(theta(izz-1)+theta(izz))*dzt*0.5
              zzz=zzz+dzt
              !print*,"  PE=",zdo," TKE=",qtke(iz)," z=",zw(izz)
              if (qtke(iz).lt.zdo .and. qtke(iz).ge.zdo_sup) then
                 bbb=(theta(izz)-theta(izz-1))/dzt
                 if (bbb .ne. 0.) then
                    tl=(beta*(theta(izz)-theta(iz))+ &
                      & sqrt( max(0.,(beta*(theta(izz)-theta(iz)))**2 + &
                      &       2.*bbb*beta*(qtke(iz)-zdo_sup))))/bbb/beta
                 else
                    if (theta(izz) .ne. theta(iz)) then
                       tl=(qtke(iz)-zdo_sup)/(beta*(theta(izz)-theta(iz)))
                    else
                       tl=0.
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
        dlu(iz)=MAX(0.1, dlu(iz)/(1. + (dlu(iz)/Lmax)) )
        dld(iz)=MAX(0.1, dld(iz)/(1. + (dld(iz)/Lmax)) )

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
    &            xland,closure,                               &
    &            dz, dx, zw,                                  &
    &            u, v, thl, thv, thlv, ql, qw,                &
    &            qke, tsq, qsq, cov,                          &
    &            vt, vq,                                      &
    &            rmol, flt, fltv, flq,                        &
    &            pblh,theta,                                  &
    &            sh, sm,                                      &
    &            El,                                          &
    &            Dfm, Dfh, Dfq, Tcd, Qcd, Pdk, Pdt, Pdq, Pdc, &
    &		 qWT1,qSHEAR1,qBUOY1,qDISS1,                  &
    &            tke_budget,                                  &
    &            Psig_bl,Psig_shcu,cldfra_bl1,                &
    &            bl_mynn_mixlength,                           &
    &            edmf_w1,edmf_a1,                             &
    &            edmf_w_dd1,edmf_a_dd1,                       &
    &            TKEprod_dn,TKEprod_up,                       &
    &            spp_pbl,pattern_spp_pbl1                     )

!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in)               :: bl_mynn_mixlength,tke_budget
    real(kind_phys), intent(in)       :: closure
    real(kind_phys), dimension(kts:kte),   intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in)       :: rmol,flt,fltv,flq,                &
         &Psig_bl,Psig_shcu,xland,dx,pblh
    real(kind_phys), dimension(kts:kte), intent(in) :: u,v,thl,thv,        &
         &thlv,qw,ql,vt,vq,qke,tsq,qsq,cov,cldfra_bl1,edmf_w1,edmf_a1,     &
         &edmf_w_dd1,edmf_a_dd1,TKEprod_dn,TKEprod_up

    real(kind_phys), dimension(kts:kte), intent(out) :: dfm,dfh,dfq,       &
         &pdk,pdt,pdq,pdc,tcd,qcd,el

    real(kind_phys), dimension(kts:kte), intent(inout) ::                  &
         qWT1,qSHEAR1,qBUOY1,qDISS1
    real(kind_phys):: q3sq_old,dlsq1,qWTP_old,qWTP_new
    real(kind_phys):: dudz,dvdz,dTdz,upwp,vpwp,Tpwp

    real(kind_phys), dimension(kts:kte) :: qkw,dtl,dqw,dtv,gm,gh,sm,sh

    integer :: k
!    real(kind_phys):: cc2,cc3,e1c,e2c,e3c,e4c,e5c
    real(kind_phys):: e6c,dzk,afk,abk,vtt,vqq,                             &
         &cw25,clow,cupp,gamt,gamq,smd,gamv,elq,elh

    real(kind_phys):: cldavg
    real(kind_phys), dimension(kts:kte), intent(in) :: theta

    real(kind_phys)::  a2fac, duz, ri !JOE-Canuto/Kitamura mod

    real:: auh,aum,adh,adm,aeh,aem,Req,Rsl,Rsl2,                           &
           gmelq,sm20,sh20,sm25max,sh25max,sm25min,sh25min,                &
           sm_pbl,sh_pbl,pblh2,wt,slht,wtpr,mfmax

    DOUBLE PRECISION  q2sq, t2sq, r2sq, c2sq, elsq, gmel, ghel
    DOUBLE PRECISION  q3sq, t3sq, r3sq, c3sq, dlsq, qdiv
    DOUBLE PRECISION  e1, e2, e3, e4, enum, eden, wden

!   Stochastic
    integer,         intent(in)                   :: spp_pbl
    real(kind_phys), dimension(kts:kte)           :: pattern_spp_pbl1
    real(kind_phys):: Prnum, shb, Prlim
    real(kind_phys), parameter :: Prlimit = 5.0

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
!

    CALL mym_level2 (kts,kte,                   &
    &            dz,                            &
    &            u, v, thl, thv, thlv,          &
    &            qw, ql, vt, vq,                &
    &            dtl, dqw, dtv, gm, gh, sm, sh  )
!
    CALL mym_length (                           &
    &            kts,kte,xland,                 &
    &            dz, dx, zw,                    &
    &            rmol, flt, fltv, flq,          &
    &            vt, vq,                        &
    &            u, v, qke,                     &
    &            dtv,                           &
    &            el,                            &
    &            pblh,theta,                    &
    &            qkw,Psig_bl,cldfra_bl1,        &
    &            bl_mynn_mixlength,             &
    &            edmf_w1,edmf_a1,               &
    &            edmf_w_dd1,edmf_a_dd1          )
!

    DO k = kts+1,kte
       dzk = 0.5  *( dz(k)+dz(k-1) )
       afk = dz(k)/( dz(k)+dz(k-1) )
       abk = one -afk
       elsq = el (k)**2
       q3sq = qkw(k)**2
       q2sq = b1*elsq*( sm(k)*gm(k)+sh(k)*gh(k) )

       sh20 = MAX(sh(k), 1e-5)
       sm20 = MAX(sm(k), 1e-5)
       sh(k)= MAX(sh(k), 1e-5)

       !Canuto/Kitamura mod
       duz = ( u(k)-u(k-1) )**2 +( v(k)-v(k-1) )**2
       duz =   duz                    /dzk**2
       !   **  Gradient Richardson number  **
       ri = -gh(k)/MAX( duz, 1.0e-10 )
       if (CKmod .eq. 1) then
          a2fac = one/(one + max(ri,zero))
       else
          a2fac = one
       endif
       !end Canuto/Kitamura mod

       !level 2.0 Prandtl number
       !Prnum = MIN(sm20/sh20, 4.0)
       !The form of Zilitinkevich et al. (2006) but modified
       !half-way towards Esau and Grachev (2007, Wind Eng)
       !Prnum = MIN(0.76 + 3.0*MAX(ri,0.0), Prlimit)
       Prnum = MIN(0.76_kind_phys + 4.0_kind_phys*MAX(ri,zero), Prlimit)
       !Prnum = MIN(0.76 + 5.0*MAX(ri,0.0), Prlimit)

       !TZC - Kondo Correction
       if (ri >= one) then
          ! Kh/Km = 1/(7*Ri)
          Prlim = 7._kind_phys*ri
       elseif (ri >= 0.01 .and. ri <= one) then
          ! Kh/Km(i,k) = 1/(6.873*Ri + 1/(6.873*Ri))
          Prlim = (6.873_kind_phys*ri + one/(6.873_kind_phys*ri))
       else
          ! no Pr limit required?
          Prlim = Prlimit
       end if
!     
!  Modified: Dec/22/2005, from here, (dlsq -> elsq)
       gmel = gm (k)*elsq
       ghel = gh (k)*elsq
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
          !sm(k) = Prnum*sh(k)
          !sm(k) = sm(k) * qdiv

          !Use level 2.0 functions as in original MYNN
          sh(k) = sh(k) * qdiv
          sm(k) = sm(k) * qdiv
        !  !sm_pbl = sm(k) * qdiv
        !
        !  !Or, use the simple Pr relationship
        !  sm(k) = Prnum*sh(k)
        !
        !  !or blend them:
        !  pblh2   = MAX(pblh, 300.)
        !  wt    =.5*TANH((zw(k) - pblh2)/200.) + .5
        !  sm(k) = sm_pbl*(1.-wt) + sm(k)*wt

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
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden  - retro 5
          !sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
          !sm(k) = Prnum*sh(k)
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
          sm(k) = q3sq*a1*( e3-3.0*c1*e4       )/eden
        !  sm_pbl = q3sq*a1*( e3-3.0*c1*e4       )/eden
          !!JOE-Canuto/Kitamura mod
          !!sh(k) = q3sq*a2*( e2+3.0*c1*e5c*gmel )/eden
          sh(k) = q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel )/eden
        !  sm(k) = Prnum*sh(k)

        !  !or blend them:
        !  pblh2   = MAX(pblh, 300.)
        !  wt    = .5*TANH((zw(k) - pblh2)/200.) + .5
        !  sm(k) = sm_pbl*(1.-wt) + sm(k)*wt
       END IF !end Helfand & Labraga check

       !Impose broad limits on Sh and Sm:
       gmelq    = MAX(gmel/q3sq, 1e-8_kind_phys)
       sm25max  = 4._kind_phys  !MIN(sm20*3.0, SQRT(.1936/gmelq))
       sh25max  = 4._kind_phys  !MIN(sh20*3.0, 0.76*b2)
       sm25min  = zero !MAX(sm20*0.1, 1e-6)
       sh25min  = zero !MAX(sh20*0.1, 1e-6)

       !JOE: Level 2.5 debug prints
       ! HL88 , lev2.5 criteria from eqs. 3.17, 3.19, & 3.20
       IF ( debug_code ) THEN
         IF ((sh(k)<sh25min .OR. sm(k)<sm25min .OR. &
              sh(k)>sh25max .OR. sm(k)>sm25max) ) THEN
           print*,"In mym_turbulence 2.5: k=",k
           print*," sm=",sm(k)," sh=",sh(k)
           print*," ri=",ri," Pr=",sm(k)/MAX(sh(k),1e-8)
           print*," gm=",gm(k)," gh=",gh(k)
           print*," q2sq=",q2sq," q3sq=",q3sq, q3sq/q2sq
           print*," qke=",qke(k)," el=",el(k)
           print*," PBLH=",pblh," u=",u(k)," v=",v(k)
           print*," SMnum=",q3sq*a1*( e3-3.0*c1*e4)," SMdenom=",eden
           print*," SHnum=",q3sq*(a2*a2fac)*( e2+3.0*c1*e5c*gmel ),&
                  " SHdenom=",eden
         ENDIF
       ENDIF

       !Enforce constraints for level 2.5 functions
       IF ( sh(k) > sh25max ) sh(k) = sh25max
       IF ( sh(k) < sh25min ) sh(k) = sh25min
       !IF ( sm(k) > sm25max ) sm(k) = sm25max
       !IF ( sm(k) < sm25min ) sm(k) = sm25min

       !shb   = max(sh(k), 0.002)
       shb   = max(sh(k), 0.02)
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
          cw25 = e1*( e2 + 3.0*c1*e5c*gmel*qdiv**2 )/( 3.0*eden )
!
!     **  Limitation on q, instead of L/q  **
          dlsq =  elsq
          IF ( q3sq/dlsq .LT. -gh(k) ) q3sq = -dlsq*gh(k)
!
!     **  Limitation on c3sq (0.12 =< cw =< 0.76) **
          ! Use Janjic's (2001; p 13-17) methodology (eqs 4.11-414 and 5.7-5.10)
          ! to calculate an exact limit for c3sq:
          auh = 27.*a1*((a2*a2fac)**2)*b2*(gtr)**2
          aum = 54.*(a1**2)*(a2*a2fac)*b2*c1*(gtr)
          adh = 9.*a1*((a2*a2fac)**2)*(12.*a1 + 3.*b2)*(gtr)**2
          adm = 18.*(a1**2)*(a2*a2fac)*(b2 - 3.*(a2*a2fac))*(gtr)

          aeh = (9.*a1*((a2*a2fac)**2)*b1 +9.*a1*((a2*a2fac)**2)* &
                (12.*a1 + 3.*b2))*(gtr)
          aem = 3.*a1*(a2*a2fac)*b1*(3.*(a2*a2fac) + 3.*b2*c1 + &
                (18.*a1*c1 - b2)) + &
                (18.)*(a1**2)*(a2*a2fac)*(b2 - 3.*(a2*a2fac))

          Req = -aeh/aem
          Rsl = (auh + aum*Req)/(3.*adh + 3.*adm*Req)
          !For now, use default values, since tests showed little/no sensitivity
          Rsl = .12_kind_phys   !lower limit
          Rsl2= one - 2.*Rsl    !upper limit
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
             !JOE: test dynamic limits
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
       cldavg = half*(cldfra_bl1(k-1) + cldfra_bl1(k))
       mfmax  = max(edmf_a1(k)*edmf_w1(k),abs(edmf_a_dd1(k)*edmf_w_dd1(k)))
       ! impose minimum for mass-flux columns
       sm(k) = max(sm(k), 0.04*min(10.*mfmax,one) )
       sh(k) = max(sh(k), 0.04*min(10.*mfmax,one) )
       ! impose minimum for clouds
       sm(k) = max(sm(k), 0.04*min(cldavg,one) )
       sh(k) = max(sh(k), 0.04*min(cldavg,one) )
!
       elq = el(k)*qkw(k)
       elh = elq*qdiv

       ! Production of TKE (pdk), T-variance (pdt),
       ! q-variance (pdq), and covariance (pdc)
       pdk(k) = elq*( sm(k)*gm(k)                &
            &        +sh(k)*gh(k)+gamv ) +       &
            &    0.5*TKEprod_dn(k)       +       & ! xmchen
            &    0.5*TKEprod_up(k)
       pdt(k) = elh*( sh(k)*dtl(k)+gamt )*dtl(k)
       pdq(k) = elh*( sh(k)*dqw(k)+gamq )*dqw(k)
       pdc(k) = elh*( sh(k)*dtl(k)+gamt )        &
            &   *dqw(k)*0.5                      &
            & + elh*( sh(k)*dqw(k)+gamq )*dtl(k)*0.5

       ! Contergradient terms
       tcd(k) = elq*gamt
       qcd(k) = elq*gamq

       ! Eddy Diffusivity/Viscosity divided by dz
       dfm(k) = elq*sm(k) / dzk
       dfh(k) = elq*sh(k) / dzk
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
       qBUOY1(k) = elq*(sh(k)*gh(k)+gamv)   +         &
       &           0.5*TKEprod_dn(k)        +         & ! xmchen
       &           0.5*TKEprod_up(k) 

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
       &            delt,                                               &
       &            dz,                                                 &
       &            ust, flt, flq, pmz, phh,                            &
       &            el,  dfq, rho,                                      &
       &            pdk, pdt, pdq, pdc,                                 &
       &            qke, tsq, qsq, cov,                                 &
       &            s_aw1,s_awqke1,bl_mynn_edmf_tke,                    &
       &            qWT1, qDISS1,tke_budget                             )

!-------------------------------------------------------------------
    integer, intent(in) :: kts,kte    

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), intent(in)    :: closure
    integer, intent(in) :: bl_mynn_edmf_tke,tke_budget
    real(kind_phys), dimension(kts:kte), intent(in) :: dz, dfq, el, rho
    real(kind_phys), dimension(kts:kte), intent(inout) :: pdk, pdt, pdq, pdc
    real(kind_phys), intent(in)    :: flt, flq, pmz, phh
    real(kind_phys), intent(in)    :: ust, delt
    real(kind_phys), dimension(kts:kte), intent(inout) :: qke,tsq, qsq, cov
! WA 8/3/15
    real(kind_phys), dimension(kts:kte+1), intent(inout) :: s_awqke1,s_aw1
    
    !!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB 
    real(kind_phys), dimension(kts:kte), intent(out) :: qWT1, qDISS1
    real(kind_phys), dimension(kts:kte) :: tke_up,dzinv  
    !! >> EOB
    
    integer :: k
    real(kind_phys), dimension(kts:kte) :: qkw, bp, rp, df3q
    real(kind_phys):: vkz,pdk1,phm,pdt1,pdq1,pdc1,b1l,b2l,onoff
    real(kind_phys), dimension(kts:kte) :: dtz
    real(kind_phys), dimension(kts:kte) :: a,b,c,d,x

    real(kind_phys), dimension(kts:kte) :: rhoinv
    real(kind_phys), dimension(kts:kte+1) :: rhoz,kqdz,kmdz

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    IF (bl_mynn_edmf_tke == 0) THEN
       onoff=zero
    ELSE
       onoff=one
    ENDIF

!   **  Strictly, vkz*h(i,j) -> karman*( 0.5*dz(1)*h(i,j)+z0 )  **
    vkz = karman*0.5*dz(kts)
!
!   **  dfq for the TKE is 3.0*dfm.  **
!
    DO k = kts,kte
!!       qke(k) = MAX(qke(k), zero)
       qkw(k) = SQRT( MAX( qke(k), zero ) )
       df3q(k)=Sqfac*dfq(k)
       dtz(k)=delt/dz(k)
    END DO
!
!JOE-add conservation + stability criteria
    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhoz(kts)  =rho(kts)
    rhoinv(kts)=1./rho(kts)
    kqdz(kts)  =rhoz(kts)*df3q(kts)
    kmdz(kts)  =rhoz(kts)*dfq(kts)
    DO k=kts+1,kte
       rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4)
       rhoinv(k)=1./MAX(rho(k),1E-4)
       kqdz(k)  = rhoz(k)*df3q(k) ! for TKE
       kmdz(k)  = rhoz(k)*dfq(k)  ! for T'2, q'2, and T'q'
    ENDDO
    rhoz(kte+1)=rhoz(kte)
    kqdz(kte+1)=rhoz(kte+1)*df3q(kte)
    kmdz(kte+1)=rhoz(kte+1)*dfq(kte)

    !stability criteria for mf
    DO k=kts+1,kte-1
       kqdz(k) = MAX(kqdz(k),  0.5* s_aw1(k))
       kqdz(k) = MAX(kqdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)))
       kmdz(k) = MAX(kmdz(k),  0.5* s_aw1(k))
       kmdz(k) = MAX(kmdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)))
!       kqdz(k) = max(kqdz(k),  0.5*(s_aw1(k)+sd_aw1(k)))
!       kqdz(k) = max(kqdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)) -0.5*(sd_aw1(k)-sd_aw1(k+1)) )
!       kmdz(k) = max(kmdz(k),  0.5*(s_aw1(k)+sd_aw1(k)))
!       kmdz(k) = max(kmdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)) -0.5*(sd_aw1(k)-sd_aw1(k+1)) )
    ENDDO
    !end conservation mods

    pdk1 = 2.0*ust**3*pmz/( vkz )
    phm  = 2.0/ust   *phh/( vkz )
    pdt1 = phm*flt**2
    pdq1 = phm*flq**2
    pdc1 = phm*flt*flq
!
!   **  pdk(1)+pdk(2) corresponds to pdk1.  **
    pdk(kts) = pdk1 - pdk(kts+1)

!!    pdt(kts) = pdt1 -pdt(kts+1)
!!    pdq(kts) = pdq1 -pdq(kts+1)
!!    pdc(kts) = pdc1 -pdc(kts+1)
    pdt(kts) = pdt(kts+1)
    pdq(kts) = pdq(kts+1)
    pdc(kts) = pdc(kts+1)
!
!   **  Prediction of twice the turbulent kinetic energy  **
!!    DO k = kts+1,kte-1
    DO k = kts,kte-1
       b1l = b1*0.5*( el(k+1)+el(k) )
       bp(k) = 2.*qkw(k) / b1l
       rp(k) = pdk(k+1) + pdk(k)
    END DO

!!    a(1)=0.
!!    b(1)=1.
!!    c(1)=-1.
!!    d(1)=0.

! Since df3q(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*df3q(k+1)+bp(k)*delt.
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

!!    DO k=kts+1,kte-1
!!       a(k-kts+1)=-dtz(k)*df3q(k)
!!       b(k-kts+1)=1.+dtz(k)*(df3q(k)+df3q(k+1))
!!       c(k-kts+1)=-dtz(k)*df3q(k+1)
!!       d(k-kts+1)=rp(k)*delt + qke(k) - qke(k)*bp(k)*delt
!!    ENDDO

!! "no flux at top"
!    a(kte)=-1. !0.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.
!! "prescribed value"
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qke(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)

    DO k=kts,kte
!       qke(k)=max(d(k-kts+1), qkemin)
       qke(k)=max(x(k), qkemin)
       qke(k)=min(qke(k), 150.)
    ENDDO
      
   
!!  TKE budget  (Puhales, 2020, WRF 4.2.1)  << EOB 
    IF (tke_budget .eq. 1) THEN
       !! TKE Vertical transport << EOBvt
        tke_up=0.5*qke
        dzinv=1./dz
        k=kts
        qWT1(k)=dzinv(k)*(                                           &
            &  (kqdz(k+1)*(tke_up(k+1)-tke_up(k))-kqdz(k)*tke_up(k)) &
            &  + 0.5*rhoinv(k)*(s_aw1(k+1)*tke_up(k+1)               &
            &  +      (s_aw1(k+1)-s_aw1(k))*tke_up(k)                &
            &  +      (s_awqke1(k)-s_awqke1(k+1)))*onoff) !unstaggered
        DO k=kts+1,kte-1
            qWT1(k)=dzinv(k)*(                                       &
            & (kqdz(k+1)*(tke_up(k+1)-tke_up(k))-kqdz(k)*(tke_up(k)-tke_up(k-1))) &
            &  + 0.5*rhoinv(k)*(s_aw1(k+1)*tke_up(k+1)               &
            &  +      (s_aw1(k+1)-s_aw1(k))*tke_up(k)                &
            &  -                  s_aw1(k)*tke_up(k-1)               &
            &  +      (s_awqke1(k)-s_awqke1(k+1)))*onoff) !unstaggered
        ENDDO
        k=kte
        qWT1(k)=dzinv(k)*(-kqdz(k)*(tke_up(k)-tke_up(k-1))           &
            &  + 0.5*rhoinv(k)*(-s_aw1(k)*tke_up(k)-s_aw1(k)*tke_up(k-1)+s_awqke1(k))*onoff) !unstaggered
        !!  >> EOBvt
        qDISS1=bp*tke_up !! TKE dissipation rate !unstaggered
    END IF
!! >> EOB 
   
    IF ( closure > 2.5 ) THEN

       !   **  Prediction of the moisture variance  **
       DO k = kts,kte-1
          b2l   = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
          rp(k) = pdq(k+1) + pdq(k)
       END DO

       !zero gradient for qsq at bottom and top
       !a(1)=0.
       !b(1)=1.
       !c(1)=-1.
       !d(1)=0.

       ! Since dfq(kts)=0.0, a(1)=0.0 and b(1)=1.+dtz(k)*dfq(k+1)+bp(k)*delt.
       DO k=kts,kte-1
          a(k)=   - dtz(k)*kmdz(k)*rhoinv(k)
          b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + qsq(k)
       ENDDO

       a(kte)=-1. !0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=0.

!       CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
       
       DO k=kts,kte
          !qsq(k)=d(k-kts+1)
          qsq(k)=MAX(x(k),1e-17)
       ENDDO
    ELSE
       !level 2.5 - use level 2 diagnostic
       DO k = kts,kte-1
          IF ( qkw(k) .LE. zero ) THEN
             b2l = zero
          ELSE
             b2l = b2*0.25*( el(k+1)+el(k) )/qkw(k)
          END IF
          qsq(k) = b2l*( pdq(k+1)+pdq(k) )
       END DO
       qsq(kte)=qsq(kte-1)
    END IF
!!!!!!!!!!!!!!!!!!!!!!end level 2.6   

    IF ( closure .GE. 3.0 ) THEN
!
!   **  dfq for the scalar variance is 1.0*dfm.  **
!
!   **  Prediction of the temperature variance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
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
          b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
          c(k)=   - dtz(k)*kmdz(k+1)*rhoinv(k)
          d(k)=rp(k)*delt + tsq(k)
       ENDDO

!!       DO k=kts+1,kte-1
!!          a(k-kts+1)=-dtz(k)*dfq(k)
!!          b(k-kts+1)=1.+dtz(k)*(dfq(k)+dfq(k+1))
!!          c(k-kts+1)=-dtz(k)*dfq(k+1)
!!          d(k-kts+1)=rp(k)*delt + tsq(k) - tsq(k)*bp(k)*delt
!!       ENDDO

       a(kte)=-1. !0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=0.
       
!       CALL tridiag(kte,a,b,c,d)
       CALL tridiag2(kte,a,b,c,d,x)

       DO k=kts,kte
!          tsq(k)=d(k-kts+1)
           tsq(k)=x(k)
       ENDDO

!   **  Prediction of the temperature-moisture covariance  **
!!       DO k = kts+1,kte-1
       DO k = kts,kte-1
          b2l = b2*0.5*( el(k+1)+el(k) )
          bp(k) = 2.*qkw(k) / b2l
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
          b(k)=1. + dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k) + bp(k)*delt
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
       b(kte)=1.
       c(kte)=0.
       d(kte)=0.

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
!
          tsq(k) = b2l*( pdt(k+1)+pdt(k) )
          cov(k) = b2l*( pdc(k+1)+pdc(k) )
       END DO
       
       tsq(kte)=tsq(kte-1)
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
    &            spp_pbl,pattern_spp_pbl1  )

!-------------------------------------------------------------------

    integer, intent(in)   :: kts,kte, bl_mynn_cloudpdf

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), intent(in)      :: HFX1,xland
    real(kind_phys), intent(in)      :: dx,pblh1
    real(kind_phys), dimension(kts:kte), intent(in) :: dz
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), dimension(kts:kte), intent(in) :: p,exner,thl,qw,   &
         &qv,qc,qi,qs,tsq,qsq,cov,th

    real(kind_phys), dimension(kts:kte), intent(inout) :: vt,vq,sgm

    real(kind_phys), dimension(kts:kte) :: alp,a,bet,b,ql,q1,RH
    real(kind_phys), dimension(kts:kte), intent(out) :: qc_bl1,qi_bl1, &
         &cldfra_bl1
    DOUBLE PRECISION :: t3sq, r3sq, c3sq

    real(kind_phys):: qsl,esat,qsat,dqsl,cld0,q1k,qlk,eq1,qll,           &
         &q2p,pt,rac,qt,t,xl,rsl,cpm,Fng,qww,alpha,beta,bb,              &
         &ls,wt,wt2,qpct,cld_factor,fac_damp,liq_frac,ql_ice,ql_water,   &
         &qmq,qsat_tk,q1_rh,rh_hack,zsl,maxqc,cldfra_rh,cldfra_qsq,      &
         &cldfra_rh0,cldfra_rh1,cldfra_qsq0,cldfra_qsq1
    real(kind_phys), parameter :: qpct_sfc=0.015
    real(kind_phys), parameter :: qpct_pbl=0.025
    real(kind_phys), parameter :: qpct_trp=0.030
    real(kind_phys), parameter :: rhcrit  =0.83 !for cloudpdf = 2
    real(kind_phys), parameter :: rhmax   =1.10 !for cloudpdf = 2
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
           sgm(k) = SQRT( MAX( (alp(k)**2 * MAX(el(k)**2,0.1) * &
                             b2 * MAX(Sh(k),0.03))/4. * &
                      (dqw/dzk - bet(k)*(dth/dzk ))**2 , 1.0e-10) )
           qmq   = qw(k) -qsl
           q1(k) = qmq / sgm(k)
           cldfra_bl1(K) = 0.5*( one+erf( q1(k)*rr2 ) )

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
        pblh2=MAX(10._kind_phys,pblh1)
        DO k = kts,kte-1
           zagl   = zw(k) + half*dz(k)

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
           qw_pert= qw(k) + qw(k)*0.5*pattern_spp_pbl1(k)*real(spp_pbl)

           !This form of qmq (the numerator of Q1) no longer uses the a(k) factor
           qmq    = qw_pert - qsat_tk          ! saturation deficit/excess;

           !Use the form of Eq. (6) in Chaboureau and Bechtold (2002)
           !except neglect all but the first term for sig_r
           r3sq   = max( qsq(k), zero )
           !Calculate sigma using higher-order moments:
           sgm(k) = SQRT( r3sq )
           !Set constraints on sigma relative to total water
           sgm(k) = min( sgm(k), qw(k)*onethird )
           
           !introduce vertical grid spacing dependence on min sgm
           wt     = min(one, max(zero, dz(k)-100.)/500.) !=0 for dz < 100 m, =1 for dz > 600 m
           sgm(k) = sgm(k) + sgm(k)*0.2*wt !inflate sgm for coarse dz

           !allow minimum sgm to vary with z.
           wt     = min(one, max(zero, (zagl - (pblh2+10.)))/300.) !0 in pbl, 1 aloft
           qpct   = qpct_pbl*(one-wt) + qpct_trp*wt
           zsl    = min(150., max(50., 0.1*pblh2))        !crude ekman layer
           wt     = min(one, max(zero, zagl - zsl)/200.)  !0 near sfc, 1 above 
           qpct   = qpct_sfc*(one-wt) + qpct*wt
           sgm(k) = max( sgm(k), qw(k)*qpct )

           !in saturated conditions, apply lower limit on sgm
           if (qmq .ge. zero) sgm(k) = max(0.02*qw(k), sgm(k))
           
           q1(k)  = qmq  / sgm(k)  ! Q1, the normalized saturation

           !Add condition for falling/settling into low-RH layers, so at least
           !some cloud fraction is applied for all qc, qs, and qi.
           rh_hack= rh(k)
           wt2    = min(one, max(zero, zagl - pblh2)/300.) !0 in pbl, 1 aloft
           !ensure adequate RH & q1 when qi is at least 1e-9 (above the PBLH)
           if ((qi(k)+qs(k))>1.e-9 .and. (zagl .gt. pblh2)) then
              rh_hack =min(rhmax, rhcrit + wt2*0.045*(9.0 + log10(qi(k)+qs(k))))
              rh(k)   =max(rh(k), rh_hack)
              !add rh-based q1
              q1_rh   =-3. + 3.*(rh(k)-rhcrit)/(one-rhcrit)
              q1(k)   =max(q1_rh, q1(k) )
           endif
           !ensure adequate rh & q1 when qc is at least 1e-6 (above the PBLH)
           if (qc(k)>1.e-6 .and. (zagl .gt. pblh2)) then
              rh_hack =min(rhmax, rhcrit + wt2*0.08*(6.0 + log10(qc(k))))
              rh(k)   =max(rh(k), rh_hack)
              !add rh-based q1
              q1_rh   =-3. + 3.*(rh(k)-rhcrit)/(one-rhcrit)
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
           !cldfra_bl1(k) = max(zero, min(one, half+0.36*atan(1.65*q1k)))

           !For clouds within the pbl, force higher saturation to make clouds
           wt2           = min(one, max(zero, (zagl - (pblh1-100.))/200.)) !0 in pbl, 1 aloft

           cldfra_qsq0   = max(zero, min(one, half+0.36*atan(3.1*(q1k))))
           cldfra_qsq1   = max(zero, min(one, half+0.36*atan(2.1*(q1k+0.2))))
           cldfra_qsq    = cldfra_qsq0*(one-wt2) + cldfra_qsq1*wt2

           !For ceiling detection, apply minimum rh-based cloud fraction
           cldfra_rh0    = min(one, max(zero, 0.56*tanh((rh(k)-0.976)/0.030)+half))
           cldfra_rh1    = min(one, max(zero, 0.55*tanh((rh(k)-0.955)/0.055)+half))
           cldfra_rh     = cldfra_rh0*(one-wt2) + cldfra_rh1*wt2

           cldfra_bl1(k) = max(cldfra_qsq, cldfra_rh)
           
           ! Specify hydrometeors
           ! The cloud water formulations are taken from CB02, Eq. 8.
           maxqc = max(qw(k) - qsat_tk, zero)
           if (q1k < zero) then        !unsaturated
              !orig: ql_water = sgm(k)*exp(1.2*q1k-one)
              !orig: ql_ice   = sgm(k)*exp(1.2*q1k-one)
              ql_water = min(sgm(k),0.03*qw(k))*cldfra_bl1(k)
              ql_ice   = min(sgm(k),0.03*qw(k))*cldfra_bl1(k)
           elseif (q1k > 2.) then !supersaturated
              ql_water = min(sgm(k)*q1k, maxqc)
              ql_ice   =     sgm(k)*q1k
           else                      !slightly saturated (0 > q1 < 2)
              !orig: ql_water = min(sgm(k)*(exp(-1.) + 0.66*q1k + 0.086*q1k**2), maxqc)
              !orig: ql_ice   =     sgm(k)*(exp(-1.) + 0.66*q1k + 0.086*q1k**2)
              ql_water = min(sgm(k),0.02*qw(k))*(exp(-1.) + 0.66*q1k + 0.086*q1k**2)
              ql_ice   = min(sgm(k),0.02*qw(k))*(exp(-1.) + 0.66*q1k + 0.086*q1k**2)
           endif

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
              Fng = 3.0 + exp(-3.8*(q1k+1.7))
           else
              Fng = min(23.9 + exp(-1.6*(q1k+2.5)), 60._kind_phys)
           endif

           cfmax = min(cldfra_bl1(k), 0.6_kind_phys)
           !Further limit the cf going into vt & vq near the surface
           zsl   = min(max(25., 0.1*pblh2), 100.)
           wt    = min(zagl/zsl, one) !=0 at z=0 m, =1 above ekman layer
           cfmax = cfmax*wt

           bb = b(k)*t/th(k) ! bb is "b" in BCMT95.  Their "b" differs from
                             ! "b" in CB02 (i.e., b(k) above) by a factor
                             ! of T/theta.  Strictly, b(k) above is formulated in
                             ! terms of sat. mixing ratio, but bb in BCMT95 is
                             ! cast in terms of sat. specific humidity.  The
                             ! conversion is neglected here.
           qww   = one + 0.61*qw(k)
           alpha = 0.61*th(k)
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
       &delt,dz,rho,                          &
       &u,v,th,tk,qv,qc,qi,qs,qnc,qni,        &
       &psfc,p,exner,                         &
       &thl,sqv,sqc,sqi,sqs,sqw,              &
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
       &bl_mynn_mixscalars                    )

!-------------------------------------------------------------------
    integer, intent(in) :: kts,kte,i

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    integer, intent(in) :: bl_mynn_cloudmix,bl_mynn_mixqt,                &
                           bl_mynn_edmf,bl_mynn_edmf_mom,                 &
                           bl_mynn_mixscalars
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
         &qs,qni,qnc,rho,p,exner,dfq,dz,tsq,qsq,cov,tcd,qcd,              &
         &cldfra_bl1,diss_heat
    real(kind_phys), dimension(kts:kte), intent(inout) :: thl,sqw,sqv,sqc,&
         &sqi,sqs,qnwfa,qnifa,qnbca,ozone,dfm,dfh
    real(kind_phys), dimension(kts:kte), intent(inout) :: du,dv,dth,dqv,  &
         &dqc,dqi,dqs,dqni,dqnc,dqnwfa,dqnifa,dqnbca,dozone
    real(kind_phys), intent(in) :: flt,flq,flqv,flqc,uoce,voce
    real(kind_phys), intent(in) :: ust,delt,psfc,wspd
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
    real(kind_phys), dimension(kts:kte) :: a,b,c,d,x
    real(kind_phys), dimension(kts:kte+1) :: rhoz,                        & !rho on model interface
          &khdz,kmdz
    real(kind_phys):: rhs,gfluxm,gfluxp,dztop,maxdfh,mindfh,maxcf,maxKh,zw
    real(kind_phys):: t,esat,qsl,onoff,kh,km,dzk,rhosfc
    real(kind_phys):: ustdrag,ustdiff,qvflux
    real(kind_phys):: th_new,portion_qc,portion_qi,condensate,qsat
    integer :: k,kk

    !Activate nonlocal mixing from the mass-flux scheme for
    !number concentrations and aerosols (0.0 = no; 1.0 = yes)
    real(kind_phys), parameter :: nonloc = 1.0
    real(kind_phys), parameter :: nc_min = 100.0
    real(kind_phys), parameter :: ni_min = 1e-6
      
    dztop=.5*(dz(kte)+dz(kte-1))

    ! REGULATE THE MOMENTUM MIXING FROM THE MASS-FLUX SCHEME (on or off)
    ! Note that s_awu and s_awv already come in as 0.0 if bl_mynn_edmf_mom == 0, so
    ! we only need to zero-out the MF term
    IF (bl_mynn_edmf_mom == 0) THEN
       onoff=zero
    ELSE
       onoff=one
    ENDIF

    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhosfc     = psfc/(R_d*(tk(kts)+p608*qv(kts)))
    dtz(kts)   = delt/dz(kts)
    rhoz(kts)  = rho(kts)
    rhoinv(kts)= 1./rho(kts)
    khdz(kts)  = rhoz(kts)*dfh(kts)
    kmdz(kts)  = rhoz(kts)*dfm(kts)
    DO k=kts+1,kte
       dtz(k)   = delt/dz(k)
       rhoz(k)  = (rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4)
       rhoinv(k)= 1./MAX(rho(k),1E-4)
       dzk      = 0.5  *( dz(k)+dz(k-1) )
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
    if ( delp(kts) < 0.5*delp(kts+1) )delp(kts)=0.5*delp(kts+1)

    !stability criteria for mf
    DO k=kts+1,kte-1
       khdz(k) = max(khdz(k),  0.5*(s_aw1(k) +sd_aw1(k)))
       khdz(k) = max(khdz(k), -0.5*(s_aw1(k) -s_aw1(k+1))  &
                              -0.5*(sd_aw1(k)-sd_aw1(k+1)) )
       kmdz(k) = max(kmdz(k),  0.5*(s_aw1(k) +sd_aw1(k)))
       kmdz(k) = max(kmdz(k), -0.5*(s_aw1(k) -s_aw1(k+1))  &
                              -0.5*(sd_aw1(k)-sd_aw1(k+1)) )
    ENDDO

    ustdrag = MIN(ust*ust,0.99)/wspd  ! limit at ~ 20 m/s
    ustdiff = MIN(ust*ust,0.01)/wspd  ! limit at ~ 2 m/s
    dth(kts:kte) = zero  ! must initialize for moisture_check routine

!!============================================
!! u
!!============================================

    k=kts

!rho-weighted (drag in b-vector):
    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(kmdz(k+1)+rhosfc*ust**2/wspd)*rhoinv(k)       &
           & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff              &
           & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k) &
           & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff              &
           & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    d(k)=u(k)  + dtz(k)*uoce*ust**2/wspd                          &
           & - dtz(k)*rhoinv(k)*s_awu1(k+1)*onoff                 &
           & + dtz(k)*rhoinv(k)*sd_awu1(k+1)*onoff                &
           & + sub_u(k)*delt + det_u(k)*delt

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)                           &
           &  + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)*onoff               &
           &  + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)*onoff
       b(k)=1.+ dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k)              &
           &  + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*onoff  &
           &  + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))*onoff
       c(k)=  - dtz(k)*kmdz(k+1)*rhoinv(k)                        &
           &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff             &
           &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
       d(k)=u(k) + dtz(k)*rhoinv(k)*(s_awu1(k)-s_awu1(k+1))*onoff &
           &  - dtz(k)*rhoinv(k)*(sd_awu1(k)-sd_awu1(k+1))*onoff  &
           &  + sub_u(k)*delt + det_u(k)*delt
    enddo

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
    a(kte)=0
    b(kte)=1.
    c(kte)=0.
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

    k=kts

!rho-weighted (drag in b-vector):
    a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(kmdz(k+1) + rhosfc*ust**2/wspd)*rhoinv(k)    &
        &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff               &
        &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)                           &
        &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff               &
        &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
    d(k)=v(k)  + dtz(k)*voce*ust**2/wspd                         &
        &  - dtz(k)*rhoinv(k)*s_awv1(k+1)*onoff                  &
        &  + dtz(k)*rhoinv(k)*sd_awv1(k+1)*onoff                 &
        &  + sub_v(k)*delt + det_v(k)*delt

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*kmdz(k)*rhoinv(k)                          &
         & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)*onoff                 &
         & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)*onoff
       b(k)=1.+dtz(k)*(kmdz(k)+kmdz(k+1))*rhoinv(k)              &
         & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*onoff    &
         & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))*onoff
       c(k)=  -dtz(k)*kmdz(k+1)*rhoinv(k)                        &
         & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*onoff               &
         & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)*onoff
       d(k)=v(k) + dtz(k)*rhoinv(k)*(s_awv1(k)-s_awv1(k+1))*onoff&
         & - dtz(k)*rhoinv(k)*(sd_awv1(k)-sd_awv1(k+1))*onoff    &
         & + sub_v(k)*delt + det_v(k)*delt
    enddo

!! no flux at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=0.

!! specified gradient at the top
!    a(kte)=-1.
!    b(kte)=1.
!    c(kte)=0.
!    d(kte)=gradv_top*dztop

!! prescribed value
    a(kte)=0
    b(kte)=1.
    c(kte)=0.
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
    k=kts

!rho-weighted: rhosfc*x*rhoinv(k)
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)            &
       &   - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
       &   - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                      &
       &   - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
       &   - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=thl(k) + dtz(k)*rhosfc*flt*rhoinv(k) + tcd(k)*delt &
       &   - dtz(k)*rhoinv(k)*s_awthl1(k+1)                 &
       &   + dtz(k)*rhoinv(k)*sd_awthl1(k+1)                &
       & + diss_heat(k)*delt + sub_thl(k)*delt + det_thl(k)*delt

    do k=kts+1,kte-1
       a(k)= -dtz(k)*khdz(k)*rhoinv(k)                      &
       &    + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                 &
       &    + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)         &
       &  + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))      &
       &  + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)= -dtz(k)*khdz(k+1)*rhoinv(k)                    &
       &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                 &
       &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=thl(k) + tcd(k)*delt                            &
       & + dtz(k)*rhoinv(k)*(s_awthl1(k)-s_awthl1(k+1))     &
       & - dtz(k)*rhoinv(k)*(sd_awthl1(k)-sd_awthl1(k+1))   &
       & +     diss_heat(k)*delt                            &
       & +     sub_thl(k)*delt + det_thl(k)*delt
    enddo

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
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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

    k=kts

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)         &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=sqw(k)  + dtz(k)*rhosfc*flq*rhoinv(k) + qcd(k)*delt &
       &  - dtz(k)*rhoinv(k)*s_awqt1(k+1)                &
       &  + dtz(k)*rhoinv(k)*sd_awqt1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                  &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                 &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)      &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))    &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=sqw(k) + qcd(k)*delt                         &
       & + dtz(k)*rhoinv(k)*(s_awqt1(k)-s_awqt1(k+1))    &
       & - dtz(k)*rhoinv(k)*(sd_awqt1(k)-sd_awqt1(k+1))
    enddo

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
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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

    k=kts

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                  &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=sqc(k)  + dtz(k)*rhosfc*flqc*rhoinv(k) + qcd(k)*delt &
    &  - dtz(k)*rhoinv(k)*s_awqc1(k+1)                  &
    &  + dtz(k)*rhoinv(k)*sd_awqc1(k+1)                 &
    &  + det_sqc(k)*delt

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                 &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))   &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)               &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)              &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=sqc(k) + qcd(k)*delt                        &
       & + dtz(k)*rhoinv(k)*(s_awqc1(k)-s_awqc1(k+1))   &
       & - dtz(k)*rhoinv(k)*(sd_awqc1(k)-sd_awqc1(k+1)) &
       & + det_sqc(k)*delt
    enddo

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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

!rho-weighted:
    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
    & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                 &
    & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                  &
    & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                 &
    & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=sqv(k)  + dtz(k)*rhosfc*qvflux*rhoinv(k) + qcd(k)*delt &
    & - dtz(k)*rhoinv(k)*s_awqv1(k+1)                   &
    & + dtz(k)*rhoinv(k)*sd_awqv1(k+1)                  &
    & + sub_sqv(k)*delt + det_sqv(k)*delt

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                 &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))   &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)               &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)              &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=sqv(k) + qcd(k)*delt                        &
       & + dtz(k)*rhoinv(k)*(s_awqv1(k)-s_awqv1(k+1))   &
       & - dtz(k)*rhoinv(k)*(sd_awqv1(k)-sd_awqv1(k+1)) &
       & + sub_sqv(k)*delt + det_sqv(k)*delt
    enddo

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
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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

    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)      &
!    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                &
!    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)               &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=sqi(k)                                       &
!    &  - dtz(k)*rhoinv(k)*s_awqi1(k+1)                 &
    &  + dtz(k)*rhoinv(k)*sd_awqi1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)               &
!       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)               &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)   &
!       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))   &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)             &
!       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)             &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=sqi(k)                                    &
!       & + dtz(k)*rhoinv(k)*(s_awqi1(k)-s_awqi1(k+1))   &
       & - dtz(k)*rhoinv(k)*(sd_awqi1(k)-sd_awqi1(k+1))
    enddo

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
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
    d(k)=sqs(k)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=sqs(k)
    ENDDO

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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
      bl_mynn_mixscalars > 0) THEN

    DO k=kts,kte
       qni2(k)=max(qni(k), zero)
       !enforce minimum number concentration
       if (sqi(k) .gt. 1e-12)qni2(k)=max(qni2(k), ni_min)
    ENDDO
      
    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                  &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=qni2(k)                                        &
    &  - dtz(k)*rhoinv(k)*s_awqni1(k+1)                 &
    &  + dtz(k)*rhoinv(k)*sd_awqni1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                 &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))   &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)               &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)              &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=qni2(k)                                     &
       & + dtz(k)*rhoinv(k)*(s_awqni1(k)-s_awqni1(k+1)) &
       & - dtz(k)*rhoinv(k)*(sd_awqni1(k)-sd_awqni1(k+1))
    enddo

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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
      bl_mynn_mixscalars > 0) THEN

    DO k=kts,kte
       qnc2(k)=max(qnc(k),zero)
       !enforce minimum number concentration
       if (sqc(k) .gt. 1e-12)qnc2(k)=max(qnc2(k), nc_min)
    ENDDO
      
    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)        &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                  &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=qnc2(k)                                        &
    &  - dtz(k)*rhoinv(k)*s_awqnc1(k+1)                 &
    &  + dtz(k)*rhoinv(k)*sd_awqnc1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                 &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)     &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))   &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)               &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)              &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=qnc2(k)                                     &
       & + dtz(k)*rhoinv(k)*(s_awqnc1(k)-s_awqnc1(k+1)) &
       & - dtz(k)*rhoinv(k)*(sd_awqnc1(k)-sd_awqnc1(k+1))
    enddo

!! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
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
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNWFA .AND. &
      bl_mynn_mixscalars > 0) THEN

    do k=kts,kte
       qnwfa2(k)=max(qnwfa(k),zero)
    enddo
      
    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)            &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                    &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                      &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                    &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=qnwfa2(k)                                          &
    &  - dtz(k)*rhoinv(k)*s_awqnwfa1(k+1)                   &
    &  + dtz(k)*rhoinv(k)*sd_awqnwfa1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                     &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                    &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)         &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))       &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=qnwfa2(k)                                       &
       & + dtz(k)*rhoinv(k)*(s_awqnwfa1(k)-s_awqnwfa1(k+1)) &
       & - dtz(k)*rhoinv(k)*(sd_awqnwfa1(k)-sd_awqnwfa1(k+1))
    enddo

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnwfa2(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnwfa2(k)=d(k)
       qnwfa2(k)=max(x(k),zero)
    ENDDO

ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnwfa2=qnwfa
ENDIF

!============================================
! Ice-friendly aerosols ( qnifa ).
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNIFA .AND. &
      bl_mynn_mixscalars > 0) THEN

    do k=kts,kte
       qnifa2(k)=max(qnifa(k),zero)
    enddo

    k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)            &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                    &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                      &
    &  - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                    &
    &  - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
    d(k)=qnifa2(k)                                          &
    &  - dtz(k)*rhoinv(k)*s_awqnifa1(k+1)                   &
    &  + dtz(k)*rhoinv(k)*sd_awqnifa1(k+1)

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                     &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)                    &
       & + 0.5*dtz(k)*rhoinv(k)*sd_aw1(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)         &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))       &
       & + 0.5*dtz(k)*rhoinv(k)*(sd_aw1(k)-sd_aw1(k+1))
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                   &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)                  &
       & - 0.5*dtz(k)*rhoinv(k)*sd_aw1(k+1)
       d(k)=qnifa2(k)                                       &
       & + dtz(k)*rhoinv(k)*(s_awqnifa1(k)-s_awqnifa1(k+1)) &
       & - dtz(k)*rhoinv(k)*(sd_awqnifa1(k)-sd_awqnifa1(k+1))
    enddo

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnifa2(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnifa2(k)=d(k-kts+1)
       qnifa2(k)=max(x(k),zero)
    ENDDO

ELSE
    !If not mixing aerosols, set "updated" array equal to original array
    qnifa2=qnifa
ENDIF

!============================================
! Black-carbon aerosols ( qnbca ).           
!============================================
IF (bl_mynn_cloudmix > 0 .AND. FLAG_QNBCA .AND. &
      bl_mynn_mixscalars > 0) THEN

   k=kts

    a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
    b(k)=1.+dtz(k)*(khdz(k) + khdz(k+1))*rhoinv(k)           &
    & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*nonloc
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                       &
    & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*nonloc
    d(k)=qnbca(k)  - dtz(k)*rhoinv(k)*s_awqnbca1(k+1)*nonloc

    do k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)                      &
       & + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)*nonloc
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)          &
       & + 0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))*nonloc
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)                    &
       & - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)*nonloc
       d(k)=qnbca(k) + dtz(k)*rhoinv(k)*(s_awqnbca1(k)-s_awqnbca1(k+1))*nonloc
    enddo

! prescribed value
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=qnbca(kte)

!    CALL tridiag(kte,a,b,c,d)
   CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !qnbca2(k)=d(k-kts+1)
       qnbca2(k)=x(k)
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
    b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k)
    c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
    d(k)=ozone(k)

    DO k=kts+1,kte-1
       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)
       d(k)=ozone(k)
    ENDDO

! prescribed value                                                                                                           
    a(kte)=0.
    b(kte)=1.
    c(kte)=0.
    d(kte)=ozone(kte)

!    CALL tridiag(kte,a,b,c,d)
    CALL tridiag2(kte,a,b,c,d,x)
!    CALL tridiag3(kte,a,b,c,d,x)

    DO k=kts,kte
       !ozone2(k)=d(k-kts+1)
       dozone(k)=(x(k)-ozone(k))/delt
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
      IF (FLAG_QNC .AND. bl_mynn_mixscalars > 0) THEN
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
      IF (FLAG_QNI .AND. bl_mynn_mixscalars > 0) THEN
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
         Dth(k)=(thl(k) + xlvcp/exner(k)*sqc2(k)          &
           &            + xlscp/exner(k)*(sqi2(k))        & !+sqs2(k)) &
           &            - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy:
         !Dth(k)=(thl(k)*(1.+ xlvcp/MAX(tk(k),TKmin)*sqc(k)  &
         !  &               + xlscp/MAX(tk(k),TKmin)*sqi(k)) &
         !  &               - th(k))/delt
      ENDDO
    ELSE
      DO k=kts,kte
         Dth(k)=(thl(k)+xlvcp/exner(k)*sqc2(k) - th(k))/delt
         !Use form from Tripoli and Cotton (1981) with their
         !suggested min temperature to improve accuracy.
         !Dth(k)=(thl(k)*(1.+ xlvcp/MAX(tk(k),TKmin)*sqc(k))  &
         !&               - th(k))/delt
      ENDDO
    ENDIF

    !===================
    ! AEROSOL TENDENCIES
    !===================
    IF (FLAG_QNWFA .AND. FLAG_QNIFA .AND. &
        bl_mynn_mixscalars > 0) THEN
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
    IF (FLAG_QNBCA .AND. bl_mynn_mixscalars > 0) THEN
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
           if( qv(k) .gt. 2.0*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20_kind_phys,sum)
        if( aa .lt. half ) then
            do k = 1, kte
               if( qv(k) .gt. 2.0*qvmin ) then
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
       nchem, kdvel, ndvel,               &
       chem1, vd1,                        &
       rho,                               &
       flt, tcd, qcd,                     &
       dfh,                               &
       s_aw1, s_awchem1,                  &
       emis_ant_no, frp, rrfs_sd,         &
       enh_mix, smoke_dbg                 )

!-------------------------------------------------------------------
    integer, intent(in) :: kts,kte,i
    real(kind_phys), dimension(kts:kte), intent(in) :: dfh,dz,tcd,qcd
    real(kind_phys), dimension(kts:kte), intent(in) :: rho
    real(kind_phys), intent(in)    :: flt
    real(kind_phys), intent(in)    :: delt,pblh
    integer, intent(in) :: nchem, kdvel, ndvel
    real(kind_phys), dimension( kts:kte+1), intent(in) :: s_aw1
    real(kind_phys), dimension( kts:kte, nchem ), intent(inout) :: chem1
    real(kind_phys), dimension( kts:kte+1,nchem), intent(in) :: s_awchem1
    real(kind_phys), dimension( ndvel ), intent(in) :: vd1
    real(kind_phys), intent(in) :: emis_ant_no,frp
    logical, intent(in) :: rrfs_sd,enh_mix,smoke_dbg
!local vars

    real(kind_phys), dimension(kts:kte)     :: dtz
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
    real(kind_phys), parameter :: NO_threshold    = 10.0     ! For anthropogenic sources
    real(kind_phys), parameter :: frp_threshold   = 10.0     ! RAR 02/11/22: I increased the frp threshold to enhance mixing over big fires
    real(kind_phys), parameter :: pblh_threshold  = 100.0

    dztop=.5*(dz(kte)+dz(kte-1))

    DO k=kts,kte
       dtz(k)=delt/dz(k)
    ENDDO

    !Prepare "constants" for diffusion equation.
    !khdz = rho*Kh/dz = rho*dfh
    rhoz(kts)  =rho(kts)
    rhoinv(kts)=1./rho(kts)
    khdz(kts)  =rhoz(kts)*dfh(kts)

    DO k=kts+1,kte
       rhoz(k)  =(rho(k)*dz(k-1) + rho(k-1)*dz(k))/(dz(k-1)+dz(k))
       rhoz(k)  =  MAX(rhoz(k),1E-4)
       rhoinv(k)=1./MAX(rho(k),1E-4)
       dzk      = 0.5  *( dz(k)+dz(k-1) )
       khdz(k)  = rhoz(k)*dfh(k)
    ENDDO
    rhoz(kte+1)=rhoz(kte)
    khdz(kte+1)=rhoz(kte+1)*dfh(kte)

    !stability criteria for mf
    DO k=kts+1,kte-1
       khdz(k) = MAX(khdz(k),  0.5*s_aw1(k))
       khdz(k) = MAX(khdz(k), -0.5*(s_aw1(k)-s_aw1(k+1)))
    ENDDO

    !Enhanced mixing over fires
    IF ( rrfs_sd .and. enh_mix ) THEN
       DO k=kts+1,kte-1
          khdz_old  = khdz(k)
          khdz_back = pblh * 0.15 / dz(k)
          !Modify based on anthropogenic emissions of NO and FRP
          IF ( pblh < pblh_threshold ) THEN
             IF ( emis_ant_no > NO_threshold ) THEN
                khdz(k) = MAX(1.1*khdz(k),sqrt((emis_ant_no / NO_threshold)) / dz(k) * rhoz(k)) ! JLS 12/21/21
!                khdz(k) = MAX(khdz(k),khdz_back)
             ENDIF
             IF ( frp > frp_threshold ) THEN
                kmaxfire = ceiling(log(frp))
                khdz(k) = MAX(1.1*khdz(k), (1. - k/(kmaxfire*2.)) * ((log(frp))**2.- 2.*log(frp)) / dz(k)*rhoz(k)) ! JLS 12/21/21
!                khdz(k) = MAX(khdz(k),khdz_back)
             ENDIF
          ENDIF
       ENDDO
    ENDIF

  !============================================
  ! Patterned after mixing of water vapor in mynn_tendencies.
  !============================================

    DO ic = 1,nchem
       k=kts

       a(k)=  -dtz(k)*khdz(k)*rhoinv(k)
       b(k)=1.+dtz(k)*(khdz(k+1)+khdz(k))*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)
       c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k)           - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)
       d(k)=chem1(k,ic) & !dtz(k)*flt  !neglecting surface sources 
            & - dtz(k)*vd1(ic)*chem1(k,ic) &
            & - dtz(k)*rhoinv(k)*s_awchem1(k+1,ic)

       DO k=kts+1,kte-1
          a(k)=  -dtz(k)*khdz(k)*rhoinv(k)     + 0.5*dtz(k)*rhoinv(k)*s_aw1(k)
          b(k)=1.+dtz(k)*(khdz(k)+khdz(k+1))*rhoinv(k) + &
             &    0.5*dtz(k)*rhoinv(k)*(s_aw1(k)-s_aw1(k+1))
          c(k)=  -dtz(k)*khdz(k+1)*rhoinv(k) - 0.5*dtz(k)*rhoinv(k)*s_aw1(k+1)
          d(k)=chem1(k,ic) + dtz(k)*rhoinv(k)*(s_awchem1(k,ic)-s_awchem1(k+1,ic))
       ENDDO

      ! prescribed value at top
       a(kte)=0.
       b(kte)=1.
       c(kte)=0.
       d(kte)=chem1(kte,ic)

       CALL tridiag3(kte,a,b,c,d,x)

       DO k=kts,kte
          chem1(k,ic)=x(k)
       ENDDO
    ENDDO

  END SUBROUTINE mynn_mix_chem

! ==================================================================
!>\ingroup gsd_mynn_edmf
  SUBROUTINE retrieve_exchange_coeffs(kts,kte,dfm,dfh,dz,km1,kh1)

!-------------------------------------------------------------------

    integer , intent(in) :: kts,kte

    real(kind_phys), dimension(kts:kte), intent(in)  :: dz,dfm,dfh
    real(kind_phys), dimension(kts:kte), intent(out) :: km1, kh1


    integer :: k
    real(kind_phys):: dzk

    km1(kts)=0.
    kh1(kts)=0.

    do k=kts+1,kte
       dzk   = 0.5 *( dz(k)+dz(k-1) )
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

    integer, intent(in):: n
    real(kind_phys), dimension(n), intent(in) :: a,b
    real(kind_phys), dimension(n), intent(inout) :: c,d
    
    integer :: i
    real(kind_phys):: p
    real(kind_phys), dimension(n) :: q
    
    c(n)=0.
    q(1)=-c(1)/b(1)
    d(1)=d(1)/b(1)
    
    DO i=2,n
       p=1./(b(i)+a(i)*q(i-1))
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
      implicit none
!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!      b - the main diagonal
!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!      d - right part
!      x - the answer
!      kte - number of unknowns (levels)

        integer,intent(in) :: kte
        real(kind_phys), dimension(kte), intent(in) :: a,b,c,d
        real(kind_phys), dimension(kte), intent(out):: x
        real(kind_phys), dimension(kte)  :: cp,dp
        real(kind_phys):: m
        integer :: k

        ! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
        ! solve for vectors c-prime and d-prime
        do k = 2,kte
           m = b(k)-cp(k-1)*a(k)
           cp(k) = c(k)/m
           dp(k) = (d(k)-dp(k-1)*a(k))/m
        enddo
        ! initialize x
        x(kte) = dp(kte)
        ! solve for x from the vectors c-prime and d-prime
        do k = kte-1, 1, -1
           x(k) = dp(k)-cp(k)*x(k+1)
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

       implicit none
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
!! NOTES ON THE PBLH FORMULATION: The 1.5-theta-increase method defines
!!PBL heights as the level at.
!!which the potential temperature first exceeds the minimum potential.
!!temperature within the boundary layer by 1.5 K. When applied to.
!!observed temperatures, this method has been shown to produce PBL-
!!height estimates that are unbiased relative to profiler-based.
!!estimates (Nielsen-Gammon et al. 2008 \cite Nielsen_Gammon_2008). 
!! However, their study did not
!!include LLJs. Banta and Pichugina (2008) \cite Pichugina_2008  show that a TKE-based.
!!threshold is a good estimate of the PBL height in LLJs. Therefore,
!!a hybrid definition is implemented that uses both methods, weighting
!!the TKE-method more during stable conditions (PBLH < 400 m).
!!A variable tke threshold (TKEeps) is used since no hard-wired
!!value could be found to work best in all conditions.
!>\section gen_get_pblh  GSD get_pblh General Algorithm
!> @{
  SUBROUTINE GET_PBLH(KTS,KTE,pblh,thv1,qke1,zw1,dz1,landsea,kpbl)

    !---------------------------------------------------------------
    !             NOTES ON THE PBLH FORMULATION
    !
    !The 1.5-theta-increase method defines PBL heights as the level at 
    !which the potential temperature first exceeds the minimum potential 
    !temperature within the boundary layer by 1.5 K. When applied to 
    !observed temperatures, this method has been shown to produce PBL-
    !height estimates that are unbiased relative to profiler-based 
    !estimates (Nielsen-Gammon et al. 2008). However, their study did not
    !include LLJs. Banta and Pichugina (2008) show that a TKE-based 
    !threshold is a good estimate of the PBL height in LLJs. Therefore,
    !a hybrid definition is implemented that uses both methods, weighting
    !the TKE-method more during stable conditions (PBLH < 400 m).
    !A variable tke threshold (TKEeps) is used since no hard-wired
    !value could be found to work best in all conditions.
    !---------------------------------------------------------------

    integer,intent(in) :: KTS,KTE

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

    real(kind_phys), intent(out) :: pblh
    real(kind_phys), intent(in) :: landsea
    real(kind_phys), dimension(kts:kte), intent(in) :: thv1, qke1, dz1
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw1
    !LOCAL VARS
    real(kind_phys)::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv
    real(kind_phys):: delt_thv   !delta theta-v; dependent on land/sea point
    real(kind_phys), parameter :: sbl_lim  = 200. !upper limit of stable BL height (m).
    real(kind_phys), parameter :: sbl_damp = 400. !transition length for blending (m).
    integer :: I,J,K,kthv,ktke,kpbl

    !Initialize kpbl
    kpbl = 2

    !> - FIND MIN THETAV IN THE LOWEST 200 M AGL
    k = kts+1
    kthv = 1
    minthv = 9.E9
    DO WHILE (zw1(k) .LE. 200.)
    !DO k=kts+1,kte-1
       IF (minthv > thv1(k)) then
           minthv = thv1(k)
           kthv = k
       ENDIF
       k = k+1
       !IF (zw1(k) .GT. sbl_lim) exit
    ENDDO

    !> - FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
    pblh=0.
    k = kthv+1
    IF((landsea-1.5).GE.zero)THEN
        ! WATER
        delt_thv = one
    ELSE
        ! LAND
        delt_thv = 1.25
    ENDIF

    pblh=0.
    k = kthv+1
!    DO WHILE (pblh .EQ. 0.) 
    DO k=kts+1,kte-1
       IF (thv1(k) .GE. (minthv + delt_thv))THEN
          pblh = zw1(k) - dz1(k-1)* &
             & MIN((thv1(k)-(minthv + delt_thv))/ &
             & MAX(thv1(k)-thv1(k-1),1E-6),one)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) pblh = zw1(kts+1) !EXIT SAFEGUARD
       IF (pblh .NE. zero) exit
    ENDDO
    !print*,"IN GET_PBLH:",thsfc,pblh

    !> - FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
    !! THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
    !!THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE 
    !!WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
    ktke   = 1
    maxqke = MAX(qke1(kts),0.)
    !Use 5% of tke max (Kosovic and Curry, 2000; JAS)
    !TKEeps = maxtke/20. = maxqke/40.
    TKEeps = maxqke/40.
    TKEeps = MAX(TKEeps,0.01) !0.025) 
    PBLH_TKE=0.

    k = ktke+1
!    DO WHILE (PBLH_TKE .EQ. 0.) 
    DO k=kts+1,kte-1
       !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
       qtke  =MAX(0.5*qke1(k)  ,0.)      ! maximum TKE
       qtkem1=MAX(0.5*qke1(k-1),0.)
       IF (qtke .LE. TKEeps) THEN
           PBLH_TKE = zw1(k) - dz1(k-1)* &
             & MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1E-6), one)
           !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
           PBLH_TKE = MAX(PBLH_TKE,zw1(kts+1))
           !print *,"PBLH_TKE:",i,PBLH_TKE, Qke1(k)/2., zw1(kts+1)
       ENDIF
       !k = k+1
       IF (k .EQ. kte-1) PBLH_TKE = zw1(kts+1) !EXIT SAFEGUARD
       IF (PBLH_TKE .NE. 0.) exit
    ENDDO

    !> - The TKE-based PBLH can (rarely) become very large 
    !! in grid points with deep convection (> 8 km!),
    !! so a limit is imposed to not let PBLH_TKE exceed the
    !! theta_v-based PBL height +/- 350 m.
    !! This has no impact on 99% of the domain.
    PBLH_TKE = MIN(PBLH_TKE,pblh+350.)
    PBLH_TKE = MAX(PBLH_TKE,MAX(pblh-350.,10.))

    wt=.5*TANH((pblh - sbl_lim)/sbl_damp) + .5
    IF (maxqke <= TKEeps) THEN
       !Cold pool situation - default to theta_v-based def
    ELSE
       !BLEND THE TWO PBLH TYPES HERE: 
       pblh=PBLH_TKE*(1.-wt) + pblh*wt
    ENDIF

    !Compute kpbl
    DO k=kts+1,kte-1
       IF ( zw1(k) >= pblh) THEN
          kpbl = k-1
          exit
       ENDIF
    ENDDO

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
!!  -# Chaboureau-Bechtold cloud fraction & coupling to radiation (when icloud_bl > 0)
!!  -# some extra limits for numerical stability
!!
!! This scheme remains under development, so consider it experimental code. 
!!
  SUBROUTINE DMP_mf(i,j,                           &
                 & kts,kte,dt,zw1,dz1,p1,rho1,     &
                 & momentum_opt,                   &
                 & tke_opt,                        &
                 & scalar_opt,                     &
                 & u1,v1,w1,th1,thl1,thv1,tk1,     &
                 & qt1,qv1,qc1,qke1,               &
                 & qnc1,qni1,qnwfa1,qnifa1,qnbca1, &
                 & ex1,vt1,vq1,sgm1,               &
                 & ust,flt,fltv,flq,flqv,          &
                 & pblh,kpbl,dx,landsea,ts,        &
            ! outputs - updraft properties   
                 & edmf_a1,edmf_w1,                &
                 & edmf_qt1,edmf_thl1,             &
                 & edmf_ent1,edmf_qc1,             &
            ! outputs - variables needed for solver 
                 & s_aw1,s_awthl1,s_awqt1,         &
                 & s_awqv1,s_awqc1,                &
                 & s_awu1,s_awv1,s_awqke1,         &
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
            ! inputs for stochastic perturbations
                 & spp_pbl,pattern_spp_pbl1,       &
                 & tkeprod_up,el1                  )

#ifdef HARDCODE_VERTICAL
# define kts 1
# define kte HARDCODE_VERTICAL
#endif

! configuration inputs:
 integer, intent(in) ::                                            &
      kts,                kte,                  kpbl,              &
      momentum_opt,       tke_opt,              scalar_opt,        &
      spp_pbl,            i,                    j
 real(kind_phys), dimension(kts:kte), intent(in)  :: pattern_spp_pbl1
! state variables
 real(kind_phys), dimension(kts:kte), intent(in)  ::               &
      &u1,v1,w1,th1,thl1,tk1,qt1,qv1,qc1,                          &
      &ex1,dz1,thv1,p1,rho1,qke1,qnc1,qni1,                        &
      &qnwfa1,qnifa1,qnbca1,el1
 real(kind_phys),dimension(kts:kte+1), intent(in) ::               &
      &zw1
 real(kind_phys), intent(in)                      ::               &
      &flt,fltv,flq,flqv,Psig_shcu,                                &
      &landsea,ts,dx,dt,ust,pblh
 logical, optional :: F_QC,F_QI,F_QNC,F_QNI,F_QNWFA,F_QNIFA,F_QNBCA

 ! outputs - updraft properties
 real(kind_phys),dimension(kts:kte), intent(inout) ::              &
      &edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1
 ! add one local edmf variable:
  real(kind_phys),dimension(kts:kte) :: edmf_th1
 ! output
 integer,         intent(inout) :: ktop
 real(kind_phys), intent(inout) :: maxmf,ztop,maxwidth
 ! outputs - variables needed for solver: sum ai*rho*wis_awphi
 real(kind_phys),dimension(kts:kte+1), intent(inout) ::            &
      &s_aw1,s_awthl1,s_awqt1,s_awqv1,s_awqc1,s_awqnc1,s_awqni1,   &
      &s_awqnwfa1,s_awqnifa1,s_awqnbca1,s_awu1,s_awv1,             &
      &s_awqke1
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
      &UPA,UPU,UPV,UPTHV,UPQKE,UPQNC,                              &
      &UPQNI,UPQNWFA,UPQNIFA,UPQNBCA
 ! entrainment variables
 real(kind_phys),dimension(kts:kte,1:NUP) :: ENT
 ! internal variables
 integer :: k,ip,k50
 real(kind_phys):: fltv2,wstar,qstar,thstar,sigmaW,sigmaQT,        &
      &sigmaTH,z0,pwmin,pwmax,wmin,wmax,wlv,Psig_w,maxw,maxqc,wpbl
 real(kind_phys):: B,QTn,THLn,THVn,QCn,Un,Vn,QKEn,QNCn,QNIn,       &
      &  QNWFAn,QNIFAn,QNBCAn,                                     &
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
 real(kind_phys),dimension(kts:kte,   nchem) :: chem1
 real(kind_phys),dimension(kts:kte+1, nchem) :: s_awchem1
 real(kind_phys),dimension(nchem)            :: chemn
 real(kind_phys),dimension(kts:kte+1,1:NUP, nchem) :: UPCHEM
 integer :: ic
 real(kind_phys),dimension(kts:kte,   nchem) :: edmf_chem
 logical, intent(in) :: mix_chem

 logical :: superadiabatic

 ! Varaibles for mass flux cloud fraction
 real(kind_phys),dimension(kts:kte), intent(inout) :: vt1, vq1, sgm1
 real(kind_phys):: sigq,xl,rsl,cpm,a,qmq,mf_cf,Aup,Q1,diffqt,qsat_tk,&
         Fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid,           &
         Ac_mf,Ac_strat,qc_mf
 real(kind_phys), parameter :: cf_thresh = 0.5 ! only overwrite stratus CF less than this value

 ! Variables for plume interpolation/saturation check
 real(kind_phys),dimension(kts:kte) :: exneri,dzi,rhoz
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
 !Subsidence
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

 !set for debugging at specific point
 !integer, parameter::idbg = 452, jdbg = 272
      
! Inititialize 2d ouput
 ktop      =0    !integer
 ztop      =zero
 maxmf     =zero
 maxwidth  =zero 
 ! Initialize individual updraft properties
 UPW       =zero
 UPTHL     =zero
 UPTHV     =zero
 UPQT      =zero
 UPA       =zero
 UPU       =zero
 UPV       =zero
 UPQC      =zero
 UPQV      =zero
 UPQKE     =zero
 UPQNC     =zero
 UPQNI     =zero
 UPQNWFA   =zero
 UPQNIFA   =zero
 UPQNBCA   =zero
 if ( mix_chem ) then
    UPCHEM(kts:kte+1,1:NUP,1:nchem)=zero
 endif
 ENT       =0.001
 ! Initialize mean updraft properties
 edmf_a1   =zero
 edmf_w1   =zero
 edmf_qt1  =zero
 edmf_thl1 =zero
 edmf_ent1 =zero
 edmf_qc1  =zero
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
 s_awqnc1  =zero
 s_awqni1  =zero
 s_awqnwfa1=zero
 s_awqnifa1=zero
 s_awqnbca1=zero
 if ( mix_chem ) then
    s_awchem1(kts:kte+1,1:nchem) = zero
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
    print*,"fltv=",fltv," Psig_shcu=",Psig_shcu," wspd=",sqrt(max(u1(kts)**2 + v1(kts)**2, 0.01_kind_phys))
 endif
      
 ! Taper off MF scheme when significant resolved-scale motions
 ! are present This function needs to be asymetric...
 maxw        = zero
 cloud_base  = 9000.0
 do k=1,kte-1
    zagl = zw1(k) + half*dz1(k)

    if (zagl > (pblh + 500.)) exit

    wpbl = w1(k)
    if (w1(k) < zero)wpbl = 2.*w1(k)
    maxw = max(maxw,abs(wpbl))

    !Find highest k-level below 50m AGL
    if (zagl <= 50.)k50=k

    !Search for cloud base
    qc_sgs = max(qc1(k), qc_bl1(k))
    if (qc_sgs> 1E-5 .and. (cldfra_bl1(k) .ge. 0.5) .and. cloud_base == 9000.0) then
       cloud_base = zw1(k) !height at interface below cloud mass level
    endif
 enddo

 !do nothing for small w (< 1 m/s), but linearly taper off for w > 1.0 m/s
 maxw = max(zero, maxw - one)
 Psig_w = max(zero, one - maxw)
 Psig_w = min(Psig_w, Psig_shcu)

 !Completely shut off MF scheme for strong resolved-scale vertical velocities.
 fltv2 = fltv
 if (Psig_w == zero .and. fltv > zero) fltv2 = -1.*fltv

 if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
    print*,"===criteria for small w in pbl:"
    print*,"maxw=",maxw," Psig_w=",Psig_w," k50=",k50
 endif
      
 ! If surface buoyancy is positive we do integration, otherwise no.
 ! Also, ensure that it is at least slightly superadiabatic up through 50 m
 superadiabatic = .false.
 if ((landsea-1.5).ge.zero) then
    hux = -0.001 ! WATER   ! dT/dz must be < - 0.1 K per 100 m.
 else
    hux = -0.003 ! LAND    ! dT/dz must be < - 0.3 K per 100 m.
 endif
 tvs = ts*(one+p608*qv1(kts))
 do k=1,max(1,k50) 
    if (k == 1) then
       dthvdz = (thv1(k)-tvs)/(half*dz1(k))
       if (dthvdz < hux) then
          superadiabatic = .true.
       else
          superadiabatic = .false.
          exit
       endif
    else
       hux = -0.0005  !allow for smaller superadiabatic layers above the surface
       dthvdz = (thv1(k)-thv1(k-1))/(0.5*(dz1(k)+dz1(k-1)))
       if (dthvdz < hux) then
          superadiabatic = .true.
       else
          superadiabatic = .false.
          exit
       endif
    endif
 enddo

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
    maxwidth_cld = min(lmax, max(0.5_kind_phys*cloud_base, 400._kind_phys))
 else                               !water
    maxwidth_cld = min(lmax, max(0.9_kind_phys*cloud_base, 400._kind_phys))
 endif
 ! Criteria (4)
 wspd_pbl=sqrt(max(u1(kts)**2 + v1(kts)**2, 0.01_kind_phys))
 !Note: area fraction (acfac) is modified below
 ! Criteria (5) - function of fltv
 if ((landsea-1.5) .lt. zero) then  !land
    maxwidth_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.040)/0.04) + .5),1000._kind_phys), zero)
 else                             !water
    maxwidth_flx = MAX(MIN(1000.*(0.6*tanh((fltv - 0.007)/0.02) + .5),1000._kind_phys), zero)
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

 !allow min plume size to increase in large flux conditions (eddy diffusivity should be
 !large enough to handle the representation of small plumes).
 !if (maxwidth .ge. (lmax - one) .and. fltv .gt. 0.2)minwidth = lmin + dlmin*min((fltv-0.2)/0.3, one) 

 if (maxwidth .le. minwidth) then ! deactivate MF component
    nup2 = 0
    maxwidth = zero
 endif

 !Begin plume processing if passes criteria
 if ( fltv2 > 0.002 .AND. (maxwidth > minwidth) .AND. superadiabatic) then

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
    acfac = half*tanh((fltv2 - 0.02)/0.05) + half

    !add a windspeed-dependent adjustment to acfac that tapers off
    !the mass-flux scheme linearly above sfc wind speeds of 13 m/s.
    !Note: this effect may be better represented by an increase in
    !entrainment rate for high wind consitions (more ambient turbulence).
    if (wspd_pbl .le. 10.) then
       ac_wsp = one
    else
       ac_wsp = one - min((max(wspd_pbl - 13.0, zero))/10., one)
    endif
    !acfac  = acfac * ac_wsp
    acfac  = min(acfac, ac_wsp)

    ! Find the portion of the total fraction (Atot) of each plume size:
    An2 = 0.
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
    z0=50.
    pwmin=0.1       ! was 0.5
    pwmax=0.4       ! was 3.0

    wstar=max(1.E-2,(gtr*fltv2*pblh)**(onethird))
    qstar=max(flq,1.0E-5)/wstar
    thstar=flt/wstar

    if ((landsea-1.5) .ge. zero) then
       csigma = 1.34   ! WATER
    else
       csigma = 1.34   ! LAND
    endif

    if (env_subs) then
       exc_fac = zero
    else
       if ((landsea-1.5).GE.zero) then
          !water: increase factor to compensate for decreased pwmin/pwmax
          exc_fac = 0.58*4.0
       else
          !land: no need to increase factor - already sufficiently large superadiabatic layers
          exc_fac = 0.58
       endif
    endif
    !decrease excess for large wind speeds
    exc_fac = exc_fac * ac_wsp

    !Note: sigmaW is typically about 0.5*wstar
    sigmaW =csigma*wstar*(z0/pblh)**(onethird)*(one - 0.8*z0/pblh)
    sigmaQT=csigma*qstar*(z0/pblh)**(onethird)
    sigmaTH=csigma*thstar*(z0/pblh)**(onethird)

    !Note: Given the pwmin & pwmax set above, these max/mins are
    !      rarely exceeded. 
    wmin=MIN(sigmaW*pwmin,0.1)
    wmax=MIN(sigmaW*pwmax,0.5)

    if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
       print*,"===excess components:"
       print*,"wstar=",wstar,"qstar=",qstar,"thstar=",thstar
       print*,"csigma=",csigma,"exc_fac=",exc_fac,"wmin=",wmin
       print*,"sigmaw=",sigmaw,"sigmaqt=",sigmaqt,"sigmath=",sigmath
    endif

    !SPECIFY SURFACE UPDRAFT PROPERTIES AT MODEL INTERFACE BETWEEN K = 1 & 2
    do ip=1,NUP
       wlv=wmin+(wmax-wmin)/real(NUP2,kind=kind_phys)*real(ip-1,kind=kind_phys)

       !SURFACE UPDRAFT VERTICAL VELOCITY
       UPW(1,ip)    =wmin + real(ip,kind=kind_phys)/real(NUP)*(wmax-wmin)
       UPU(1,ip)    =(u1(kts)    *dz1(kts+1)+u1(kts+1)    *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPV(1,ip)    =(v1(kts)    *dz1(kts+1)+v1(kts+1)    *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQC(1,ip)   =zero
       !UPQC(1,ip)  =(qc1(kts)*dz1(kts+1)+qc1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))

       exc_heat     =exc_fac*UPW(1,ip)*sigmaTH/sigmaW
       UPTHV(1,ip)  =(thv1(kts)  *dz1(kts+1)+thv1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1)) &
           &        + exc_heat
       UPTHL(1,ip)  =(thl1(kts)  *dz1(kts+1)+thl1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1)) &
           &        + exc_heat
       !calculate exc_moist by use of surface fluxes
       exc_moist    =exc_fac*UPW(1,ip)*sigmaQT/sigmaW
       UPQT(1,ip)   =(qt1(kts)   *dz1(kts+1)+qt1(kts+1)   *dz1(kts))/(dz1(kts)+dz1(kts+1))&
            &       + exc_moist
       UPQKE(1,ip)  =(qke1(kts)  *dz1(kts+1)+qke1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNC(1,ip)  =(qnc1(kts)  *dz1(kts+1)+qnc1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNI(1,ip)  =(qni1(kts)  *dz1(kts+1)+qni1(kts+1)  *dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNWFA(1,ip)=(qnwfa1(kts)*dz1(kts+1)+qnwfa1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNIFA(1,ip)=(qnifa1(kts)*dz1(kts+1)+qnifa1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))
       UPQNBCA(1,ip)=(qnbca1(kts)*dz1(kts+1)+qnbca1(kts+1)*dz1(kts))/(dz1(kts)+dz1(kts+1))
    enddo

    if ( mix_chem ) then
      do ip=1,NUP
        do ic = 1,nchem
          UPCHEM(1,ip,ic)=(chem1(kts,ic)*dz1(kts+1)+chem1(kts+1,ic)*dz1(kts))/(dz1(kts)+dz1(kts+1))
        enddo
      enddo
    endif

    !Initialize environmental variables which can be modified by detrainment
    envm_thl(kts:kte)=thl1(kts:kte)
    envm_sqv(kts:kte)=qv1(kts:kte)
    envm_sqc(kts:kte)=qc1(kts:kte)
    envm_u(kts:kte)  =u1(kts:kte)
    envm_v(kts:kte)  =v1(kts:kte)
    do k=kts,kte-1
       rhoz(k)  = (rho1(k)*dz1(k+1)+rho1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
    enddo
    rhoz(kte) = rho1(kte)

    !dxsa is scale-adaptive factor governing the pressure-gradient term of the momentum transport
    dxsa = one - MIN(MAX((12000.0-dx)/(12000.0-3000.0), zero), one)

    ! do integration  updraft
    do ip=1,NUP
       QCn = zero
       overshoot = 0 !int
       l  = minwidth + dl*real(ip-1,kind=kind_phys)    ! diameter of plume
       do k=kts+1,kte-1
          !Entrainment from Tian and Kuang (2016)
          !ENT(k,ip) = 0.35/(MIN(MAX(UPW(K-1,ip),0.75),1.9)*l)
          wmin = 0.3_kind_phys + l*0.0005_kind_phys !* MAX(pblh-zw1(k+1), zero)/pblh
          ENT(k,ip) = 0.33_kind_phys/(MIN(MAX(UPW(K-1,ip),wmin),0.9_kind_phys)*l)

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

          ENT(k,ip) = min(ENT(k,ip),0.9/(zw1(k+1)-zw1(k)))

          ! Define environment U & V at the model interface levels
          Uk     =(u1(k)*dz1(k+1)+u1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
          Ukm1   =(u1(k-1)*dz1(k)+u1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))
          Vk     =(v1(k)*dz1(k+1)+v1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
          Vkm1   =(v1(k-1)*dz1(k)+v1(k)*dz1(k-1))/(dz1(k-1)+dz1(k))

          ! Linear entrainment:
          EntExp =ENT(K,ip)*(zw1(k+1)-zw1(k))
          EntExm =EntExp*0.3333    !reduce entrainment for momentum
          QTn    =UPQT(k-1,IP)   *(1.-EntExp) + qt1(k)*EntExp
          THLn   =UPTHL(k-1,IP)  *(1.-EntExp) + thl1(k)*EntExp
          Un     =UPU(k-1,IP)    *(1.-EntExm) + u1(k)*EntExm + dxsa*pgfac*(Uk - Ukm1)
          Vn     =UPV(k-1,IP)    *(1.-EntExm) + v1(k)*EntExm + dxsa*pgfac*(Vk - Vkm1)
          QKEn   =UPQKE(k-1,IP)  *(1.-EntExp) + qke1(k)*EntExp
          QNCn   =UPQNC(k-1,IP)  *(1.-EntExp) + qnc1(k)*EntExp
          QNIn   =UPQNI(k-1,IP)  *(1.-EntExp) + qni1(k)*EntExp
          QNWFAn =UPQNWFA(k-1,IP)*(1.-EntExp) + qnwfa1(k)*EntExp
          QNIFAn =UPQNIFA(k-1,IP)*(1.-EntExp) + qnifa1(k)*EntExp
          QNBCAn =UPQNBCA(k-1,IP)*(1.-EntExp) + qnbca1(k)*EntExp

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
              chemn(ic)=UPCHEM(k-1,ip,ic)*(1.-EntExp) + chem1(k,ic)*EntExp
            enddo
          endif

          ! Define pressure at model interface
          Pk    =(p1(k)*dz1(k+1)+p1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
          ! Compute plume properties thvn and qcn
          call condensation_edmf(QTn,THLn,Pk,zw1(k+1),THVn,QCn)

          ! Define environment THV at the model interface levels
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
             Wn = UPW(K-1,ip) + (-2. * ENT(K,ip) * UPW(K-1,ip) + BCOEFF*B / MAX(UPW(K-1,ip),0.2)) * MIN(zw1(k)-zw1(k-1), 250.)
          ELSE
             Wn = UPW(K-1,ip) + (-2. * ENT(K,ip) * UPW(K-1,ip) + BCOEFF*B / UPW(K-1,ip)) * MIN(zw1(k)-zw1(k-1), 250.)
          ENDIF
          !Do not allow a parcel to accelerate more than 1.25 m/s over 200 m.
          !Add max increase of 2.0 m/s for coarse vertical resolution.
          IF(Wn > UPW(K-1,ip) + MIN(1.25*(zw1(k)-zw1(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,ip) + MIN(1.25*(zw1(k)-zw1(k-1))/200., 2.0)
          ENDIF
          !Add symmetrical max decrease in w
          IF(Wn < UPW(K-1,ip) - MIN(1.25*(zw1(k)-zw1(k-1))/200., 2.0) ) THEN
             Wn = UPW(K-1,ip) - MIN(1.25*(zw1(k)-zw1(k-1))/200., 2.0)
          ENDIF
          Wn = MIN(MAX(Wn, zero), 3.0_kind_phys)

          !Check to make sure that the plume made it up at least one level.
          !if it failed, then set nup2=0 and exit the mass-flux portion.
          IF (k==kts+1 .AND. Wn == zero) THEN
             NUP2=0
             exit
          ENDIF

          IF (debug_mf == 1 .and. i==idbg .and. j==jdbg) THEN
            IF (Wn .GE. 3.0) THEN
              ! surface values
              print *," **** SUSPICIOUSLY LARGE W:"
              print *,' QCn:',QCn,' ENT=',ENT(k,ip),' Nup2=',Nup2
              print *,'pblh:',pblh,' Wn:',Wn,' UPW(k-1)=',UPW(K-1,ip)
              print *,'K=',k,' B=',B,' dz=',zw1(k)-zw1(k-1)
            ENDIF
          ENDIF

          !Allow strongly forced plumes to overshoot if KE is sufficient
          IF (Wn <= zero .AND. overshoot == 0) THEN
             overshoot = 1
             IF ( THVk-THVkm1 .GT. zero ) THEN
                bvf = SQRT( gtr*(THVk-THVkm1)/dz1(k) )
                !vertical Froude number
                Frz = UPW(K-1,ip)/(bvf*dz1(k))
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
          aratio   = MIN(UPA(K-1,IP)/(1.-UPA(K-1,IP)), 0.5) !limit should never get hit
          detturb  = 0.00008
          oow      = -0.060/MAX(one,(0.5*(Wn+UPW(K-1,IP))))   !coef for dynamical detrainment rate
          detrate  = MIN(MAX(oow*(Wn-UPW(K-1,IP))/dz1(k), detturb), .0002) ! dynamical detrainment rate (m^-1)
          detrateUV= MIN(MAX(oow*(Wn-UPW(K-1,IP))/dz1(k), detturb), .0001) ! dynamical detrainment rate (m^-1) 
          envm_thl(k)=envm_thl(k) + (0.5*(thl_ent + UPTHL(K-1,IP)) - thl1(k))*detrate*aratio*MIN(dzp,dzpmax)
          qv_ent = 0.5*(MAX(qt_ent-qc_ent,0.) + MAX(UPQT(K-1,IP)-UPQC(K-1,IP),0.))
          envm_sqv(k)=envm_sqv(k) + (qv_ent-qv1(K))*detrate*aratio*MIN(dzp,dzpmax)
          IF (UPQC(K-1,IP) > 1E-8) THEN
             IF (qc1(k) > 1E-6) THEN
                qc_grid = qc1(k)
             ELSE
                qc_grid = qc_bl1(K)
             ENDIF
             envm_sqc(k)=envm_sqc(k) + MAX(UPA(K-1,IP)*0.5*(QCn + UPQC(K-1,IP)) - qc_grid, zero)*detrate*aratio*MIN(dzp,dzpmax)
          ENDIF
          envm_u(k)  =envm_u(k)   + (0.5*(Un + UPU(K-1,IP)) - u1(K))*detrateUV*aratio*MIN(dzp,dzpmax)
          envm_v(k)  =envm_v(k)   + (0.5*(Vn + UPV(K-1,IP)) - v1(K))*detrateUV*aratio*MIN(dzp,dzpmax)

          IF (Wn > 0.) THEN
             !Update plume variables at current k index
             UPW(K,IP)=Wn  !sqrt(Wn2)
             UPTHV(K,IP)=THVn
             UPTHL(K,IP)=THLn
             UPQT(K,IP)=QTn
             UPQC(K,IP)=QCn
             UPU(K,IP)=Un
             UPV(K,IP)=Vn
             UPQKE(K,IP)=QKEn
             UPQNC(K,IP)=QNCn
             UPQNI(K,IP)=QNIn
             UPQNWFA(K,IP)=QNWFAn
             UPQNIFA(K,IP)=QNIFAn
             UPQNBCA(K,IP)=QNBCAn
             UPA(K,IP)=UPA(K-1,IP)
             IF ( mix_chem ) THEN
               do ic = 1,nchem
                 UPCHEM(k,ip,ic) = chemn(ic)
               enddo
             ENDIF
             ktop = MAX(ktop,k)
          ELSE
             if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
                print*,"plume #:",ip," ktop=",ktop
                print*,"area=",UPA(1:ktop,ip)
                print*,"w=",UPW(1:ktop,ip)
             endif
             exit  !exit k-loop
          END IF
       ENDDO

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
          s_aw1(k+1)   = s_aw1(k+1)    + rhoz(k)*UPA(K,ip)*UPW(K,ip)*Psig_w
          s_awthl1(k+1)= s_awthl1(k+1) + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPTHL(K,ip)*Psig_w
          s_awqt1(k+1) = s_awqt1(k+1)  + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQT(K,ip)*Psig_w
          !to conform to grid mean properties, move qc to qv in grid mean
          !saturated layers, so total water fluxes are preserved but 
          !negative qc fluxes in unsaturated layers is reduced.
!         if (qc1(k) > 1e-12 .or. qc1(k+1) > 1e-12) then
             qc_plume = UPQC(K,ip)
!         else
!            qc_plume = zero
!         endif
          s_awqc1(k+1) = s_awqc1(k+1)  + rhoz(k)*UPA(K,ip)*UPW(K,ip)*qc_plume*Psig_w
          s_awqv1(k+1) = s_awqt1(k+1)  - s_awqc1(k+1)
       ENDDO
    ENDDO
    !momentum
    if (momentum_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awu1(k+1) = s_awu1(k+1)   + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPU(K,ip)*Psig_w
             s_awv1(k+1) = s_awv1(k+1)   + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPV(K,ip)*Psig_w
          enddo
       enddo
    endif
    !tke
    if (tke_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awqke1(k+1)= s_awqke1(k+1) + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQKE(K,ip)*Psig_w
          enddo
       enddo
    endif
    !chem
    if ( mix_chem ) then
       do k=kts,kte-1
          do ip=1,nup
             do ic = 1,nchem
                s_awchem1(k+1,ic) = s_awchem1(k+1,ic) + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPCHEM(K,ip,ic)*Psig_w
             enddo
          enddo
       enddo
    endif

    if (scalar_opt > 0) then
       do ip=1,nup
          do k=kts,kte-1
             s_awqnc1(k+1)  = s_awqnc1(K+1)   + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQNC(K,ip)*Psig_w
             s_awqni1(k+1)  = s_awqni1(K+1)   + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQNI(K,ip)*Psig_w
             s_awqnwfa1(k+1)= s_awqnwfa1(K+1) + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQNWFA(K,ip)*Psig_w
             s_awqnifa1(k+1)= s_awqnifa1(K+1) + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQNIFA(K,ip)*Psig_w
             s_awqnbca1(k+1)= s_awqnbca1(K+1) + rhoz(k)*UPA(K,ip)*UPW(K,ip)*UPQNBCA(K,ip)*Psig_w
          enddo
       enddo
    endif

   !Flux limiter: Check ratio of heat flux at top of first model layer
   !and at the surface. Make sure estimated flux out of the top of the
   !layer is < fluxportion*surface_heat_flux
   IF (s_aw1(kts+1) /= 0.) THEN
      dzi(kts) = half*(dz1(kts)+dz1(kts+1)) !dz centered at model interface
      flx1 = max(s_aw1(kts+1)*(thv1(kts)-thv1(kts+1))/dzi(kts),1.0e-6_kind_phys)
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
      s_awqnc1   = s_awqnc1*adjustment
      s_awqni1   = s_awqni1*adjustment
      s_awqnwfa1 = s_awqnwfa1*adjustment
      s_awqnifa1 = s_awqnifa1*adjustment
      s_awqnbca1 = s_awqnbca1*adjustment
      IF (momentum_opt > 0) THEN
         s_awu1  = s_awu1*adjustment
         s_awv1  = s_awv1*adjustment
      ENDIF
      IF (tke_opt > 0) THEN
         s_awqke1= s_awqke1*adjustment
      ENDIF
      IF ( mix_chem ) THEN
         s_awchem1 = s_awchem1*adjustment
      ENDIF
      UPA = UPA*adjustment
   ENDIF
   if (debug_mf == 1 .and. i==idbg .and. j==jdbg) then
      print*,"Flux adjustment:"
      print*,"adjustment=",adjustment," fluxportion=",fluxportion
      print*,"flt2=",flt2," flx1=",flx1
   endif
      
   !Calculate mean updraft properties for output:
   !all edmf_* variables at k=1 correspond to the interface at top of first model layer
   do ip=1,nup
      do k=kts,kte-1
         edmf_a1(k)  =edmf_a1(k)  +UPA(k,ip)
         edmf_w1(k)  =edmf_w1(k)  +UPA(k,ip)*UPW(k,ip)
         edmf_qt1(k) =edmf_qt1(k) +UPA(k,ip)*UPQT(k,ip)
         edmf_thl1(k)=edmf_thl1(k)+UPA(k,ip)*UPTHL(k,ip)
         edmf_ent1(k)=edmf_ent1(k)+UPA(k,ip)*ENT(k,ip)
         edmf_qc1(k) =edmf_qc1(k) +UPA(k,ip)*UPQC(k,ip)
      enddo
   enddo
   do k=kts,kte-1
      !Note that only edmf_a1 is multiplied by Psig_w. This takes care of the
      !scale-awareness of the subsidence below:
      if (edmf_a1(k)>0.) then
         edmf_w1(k)  =edmf_w1(k)/edmf_a1(k)
         edmf_qt1(k) =edmf_qt1(k)/edmf_a1(k)
         edmf_thl1(k)=edmf_thl1(k)/edmf_a1(k)
         edmf_ent1(k)=edmf_ent1(k)/edmf_a1(k)
         edmf_qc1(k) =edmf_qc1(k)/edmf_a1(k)
         edmf_a1(k)  =edmf_a1(k)*Psig_w
         !FIND MAXIMUM MASS-FLUX IN THE COLUMN:
         if(edmf_a1(k)*edmf_w1(k) > maxmf) maxmf = edmf_a1(k)*edmf_w1(k)
      endif
      !instead of dTKE/dt = 1/2 w^3, multiply by 2 for QKE.
      tkeprod_up(k)=(abs(edmf_w1(k))**3)*edmf_a1(k)/(b1*max(el1(k),0.1)) 
   enddo ! end k

   !smoke/chem
   if ( mix_chem ) then
      do k=kts,kte-1
        do ip=1,nup
          do ic = 1,nchem
            edmf_chem(k,ic) = edmf_chem(k,ic) + UPA(k,ip)*UPCHEM(k,ip,ic)
          enddo
        enddo
      enddo
      do k=kts,kte-1
        if (edmf_a1(k)>0.) then
          do ic = 1,nchem
            edmf_chem(k,ic) = edmf_chem(k,ic)/edmf_a1(k)
          enddo
        endif
      enddo ! end k
   endif

   !Calculate the effects environmental subsidence.
   !All envi_*variables are valid at the interfaces, like the edmf_* variables
   IF (env_subs) THEN
      DO k=kts+1,kte-1
         !First, smooth the profiles of w & a, since sharp vertical gradients
         !in plume variables are not likely extended to env variables
         !Note1: w is treated as negative further below
         !Note2: both w & a will be transformed into env variables further below
         envi_w(k) = onethird*(edmf_w1(k-1)+edmf_w1(k)+edmf_w1(k+1))
         envi_a(k) = onethird*(edmf_a1(k-1)+edmf_a1(k)+edmf_a1(k+1))*adjustment
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
      dzi(kts)     = 0.5*(dz1(kts)+dz1(kts+1))
      sub_thl1(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*thl1(kts+1)-rho1(kts)*thl1(kts))/dzi(kts)/rhoz(kts)
      sub_sqv1(kts)= 0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*qv1(kts+1)-rho1(kts)*qv1(kts))/dzi(kts)/rhoz(kts)
      DO k=kts+1,kte-1
         dzi(k)     = 0.5*(dz1(k)+dz1(k+1))
         sub_thl1(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho1(k+1)*thl1(k+1)-rho1(k)*thl1(k))/dzi(k)/rhoz(k)
         sub_sqv1(k)= 0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                      (rho1(k+1)*qv1(k+1)-rho1(k)*qv1(k))/dzi(k)/rhoz(k)
      ENDDO

      DO k=kts,KTE-1
         det_thl1(k)=Cdet*(envm_thl(k)-thl1(k))*envi_a(k)*Psig_w
         det_sqv1(k)=Cdet*(envm_sqv(k)-qv1(k))*envi_a(k)*Psig_w
         det_sqc1(k)=Cdet*(envm_sqc(k)-qc1(k))*envi_a(k)*Psig_w
      ENDDO

      IF (momentum_opt > 0) THEN
         sub_u1(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*u1(kts+1)-rho1(kts)*u1(kts))/dzi(kts)/rhoz(kts)
         sub_v1(kts)=0.5*envi_w(kts)*envi_a(kts)*                               &
                     (rho1(kts+1)*v1(kts+1)-rho1(kts)*v1(kts))/dzi(kts)/rhoz(kts)
         DO k=kts+1,kte-1
            sub_u1(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                       (rho1(k+1)*u1(k+1)-rho1(k)*u1(k))/dzi(k)/rhoz(k)
            sub_v1(k)=0.5*(envi_w(k)+envi_w(k-1))*0.5*(envi_a(k)+envi_a(k-1)) * &
                       (rho1(k+1)*v1(k+1)-rho1(k)*v1(k))/dzi(k)/rhoz(k)
         ENDDO

         DO k=kts,KTE-1
           det_u1(k) = Cdet*(envm_u(k)-u1(k))*envi_a(k)*Psig_w
           det_v1(k) = Cdet*(envm_v(k)-v1(k))*envi_a(k)*Psig_w
         ENDDO
       ENDIF
   ENDIF !end subsidence/env detranment

   !First, compute exner, plume theta, and dz centered at interface
   !Here, k=1 is the top of the first model layer. These values do not 
   !need to be defined at k=kte (unused level).
   DO k=kts,kte-1
      exneri(k)  = (ex1(k)*dz1(k+1)+ex1(k+1)*dz1(k))/(dz1(k+1)+dz1(k))
      edmf_th1(k)= edmf_thl1(k) + xlvcp/ex1(k)*edmf_qc1(K)
      dzi(k)     = 0.5*(dz1(k)+dz1(k+1))
   ENDDO

!JOE: ADD CLDFRA_bl1, qc_bl1. Note that they have already been defined in
!     mym_condensation. Here, a shallow-cu component is added, but no cumulus
!     clouds can be added at k=1 (start loop at k=2).
   do k=kts+1,kte-2
      if (k > KTOP) exit
         if(0.5*(edmf_qc1(k)+edmf_qc1(k-1))>zero .and. (cldfra_bl1(k) < cf_thresh))THEN
            !interpolate plume quantities to mass levels
            Aup = (edmf_a1(k)*dzi(k-1)+edmf_a1(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            THp = (edmf_th1(k)*dzi(k-1)+edmf_th1(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            QTp = (edmf_qt1(k)*dzi(k-1)+edmf_qt1(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            !convert TH to T
!            t = THp*exner(k)
            !SATURATED VAPOR PRESSURE
            esat = esat_blend(tk1(k))
            !SATURATED SPECIFIC HUMIDITY
            qsl=ep_2*esat/max(1.e-7,(p1(k)-ep_3*esat)) 

            !condensed liquid in the plume on mass levels
            if (edmf_qc1(k)>zero .and. edmf_qc1(k-1)>zero) then
              QCp = (edmf_qc1(k)*dzi(k-1)+edmf_qc1(k-1)*dzi(k))/(dzi(k-1)+dzi(k))
            else
              QCp = max(edmf_qc1(k),edmf_qc1(k-1))
            endif

            !COMPUTE CLDFRA & QC_BL FROM MASS-FLUX SCHEME and recompute vt & vq
            xl = xl_blend(tk1(k))               ! obtain blended heat capacity
            qsat_tk = qsat_blend(tk1(k),p1(k))  ! get saturation water vapor mixing ratio
                                                !   at t and p
            rsl = xl*qsat_tk / (r_v*tk1(k)**2)  ! slope of C-C curve at t (abs temp)
                                                ! CB02, Eqn. 4
            cpm = cp + qt1(k)*cpv               ! CB02, sec. 2, para. 1
            a   = 1./(1. + xl*rsl/cpm)          ! CB02 variable "a"
            b9  = a*rsl                         ! CB02 variable "b" 

            q2p  = xlvcp/ex1(k)
            pt = thl1(k) +q2p*QCp*Aup ! potential temp (env + plume)
            bb = b9*tk1(k)/pt ! bb is "b9" in BCMT95.  Their "b9" differs from
                           ! "b9" in CB02 by a factor
                           ! of T/theta.  Strictly, b9 above is formulated in
                           ! terms of sat. mixing ratio, but bb in BCMT95 is
                           ! cast in terms of sat. specific humidity.  The
                           ! conversion is neglected here.
            qww   = 1.+0.61*qt1(k)
            alpha = 0.61*pt
            beta  = pt*xl/(tk1(k)*cp) - 1.61*pt
            !Buoyancy flux terms have been moved to the end of this section...

            !Now calculate convective component of the cloud fraction:
            if (a > zero) then
               f = MIN(1.0/a, 4.0)              ! f is vertical profile scaling function (CB2005)
            else
               f = one
            endif

            !CB form:
            !sigq = 3.5E-3 * Aup * 0.5*(edmf_w1(k)+edmf_w1(k-1)) * f  ! convective component of sigma (CB2005)
            !sigq = SQRT(sigq**2 + sgm1(k)**2)    ! combined conv + stratus components
            !Per S.DeRoode 2009?
            !sigq = 5. * Aup * (QTp - qt1(k))
            sigq = 10. * Aup * (QTp - qt1(k)) 
            !constrain sigq wrt saturation:
            sigq = max(sigq, qsat_tk*0.02 )
            sigq = min(sigq, qsat_tk*0.25 )

            qmq = a * (qt1(k) - qsat_tk)          ! saturation deficit/excess;
            Q1  = qmq/sigq                        !   the numerator of Q1

            if ((landsea-1.5).GE.zero) then   ! WATER
               !modified form from LES
               !mf_cf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.2)),0.01),0.6)
               !Original CB
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.8)
               mf_cf = max(mf_cf, 1.2 * Aup)
               !mf_cf = min(mf_cf, 5.0 * Aup)
            else                              ! LAND
               !LES form
               !mf_cf = min(max(0.5 + 0.36 * atan(1.20*(Q1+0.4)),0.01),0.6)
               !Original CB
               mf_cf = min(max(0.5 + 0.36 * atan(1.55*Q1),0.01),0.8)
               mf_cf = max(mf_cf, 1.8 * Aup)
               !mf_cf = min(mf_cf, 5.0 * Aup)
            endif

            !if ( debug_mf == 1 .and. i==idbg .and. j==jdbg) then
            !   print*,"In MYNN-MF, macrophysics component"
            !   print*," env qt=",qt1(k)," qsat=",qsat_tk
            !   print*," k=",k," satdef=",QTp - qsat_tk," sgm=",sgm1(k)
            !   print*," sigq=",sigq," qmq=",qmq," tk=",tk1(k)
            !   print*," mf_cf=",mf_cf," cldfra_bl=",cldfra_bl1(k)," edmf_a1=",edmf_a1(k)
            !endif

            ! Update cloud fractions and specific humidities in grid cells
            ! where the mass-flux scheme is active. The specific humidities
            ! are converted to grid means (not in-cloud quantities).
            if ((landsea-1.5).GE.zero) then  ! water
               if (QCp * Aup > 5e-5) then
                  qc_bl1(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1(k) = 1.18 * (QCp * Aup)
               endif
               cldfra_bl1(k)= mf_cf
               Ac_mf        = mf_cf
            else                             ! land
               if (QCp * Aup > 5e-5) then
                  qc_bl1(k) = 1.86 * (QCp * Aup) - 2.2e-5
               else
                  qc_bl1(k) = 1.18 * (QCp * Aup)
               endif
               cldfra_bl1(k)= mf_cf
               Ac_mf        = mf_cf
            endif

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
               Fng = EXP(-0.4*(Q1-one))
            elseif (Q1 .ge. -2.5 .and. Q1 .lt. -1.7) then
               Fng = 3.0 + EXP(-3.8*(Q1+1.7))
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
!   print *,'p:',p1
 
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
  diff=1.e-6

  EXN=(P/p1000mb)**rcp
  !QC=0.  !better first guess QC is incoming from lower level, do not set to zero
  do ni=1,NITER
     T=EXN*THL + xlvcp*QC
     QS=qsat_blend(T,P)
     QCOLD=QC
     QC=0.5*QC + 0.5*MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo

  T=EXN*THL + xlvcp*QC
  QS=qsat_blend(T,P)
  QC=max(QT-QS,0.)

  !Do not allow saturation below 100 m
  if(zagl < 100.)QC=0.

  !THV=(THL+xlv/cp*QC).*(1+(1-rvovrd)*(QT-QC)-QC);
  THV=(THL+xlvcp*QC)*(1.+QT*(rvovrd-1.)-rvovrd*QC)

!  IF (QC > zero) THEN
!    PRINT*,"EDMF SAT, p:",p," iterations:",ni
!    PRINT*," T=",T," THL=",THL," THV=",THV
!    PRINT*," QS=",QS," QT=",QT," QC=",QC,"ratio=",qc/qs
!  ENDIF

  !THIS BASICALLY GIVE THE SAME RESULT AS THE PREVIOUS LINE
  !TH = THL + xlv/cp/EXN*QC
  !THV= TH*(1. + p608*QT)

  !print *,'t,p,qt,qs,qc'
  !print *,t,p,qt,qs,qc 


end subroutine condensation_edmf

!===============================================================
subroutine condensation_ddmf(qt,thl,p,zagl,thv,qc,qi,debug_dd)
!
! zero or one condensation for edmf: calculates thv, qc, and qi
!
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
     t=exn*thl + xlvcp*qc + xlscp*qi
     qs=qsat_blend(t,p)
     qxold=qc+qi
     qx=qc+qi
     qx=0.5*qx + 0.5*max((qt-qs),0.)
     qc=frac_liq*qx
     qi=frac_ice*qx
     if (abs(qx-qxold)<diff) exit
  enddo

  t=exn*thl + xlvcp*qc + xlscp*qi
  qs=qsat_blend(t,p)
  qx=max(qt-qs,0.)
  qc=frac_liq*qx
  qi=frac_ice*qx

  !thv=(thl+xlv/cp*qc).*(1+(1-rvovrd)*(qt-qc)-qc)
!was this: thv=(thl+xlvcp*qc)*(1.+qt*(rvovrd-1.)-rvovrd*qc)
  !thv=(thl + xlvcp*qc + xlscp*qi)*(1. + qt*(rvovrd-1.)-rvovrd*qc)
  thv=(thl + xlvcp*qc + xlscp*qi)*(1. + p608*(qt-qx))

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
real(kind_phys),intent(in)   :: QT,THV,P,zagl
real(kind_phys),intent(inout):: THL, QC

integer :: niter,ni
real(kind_phys):: diff,exn,t,th,qs,qcold

! number of iterations                                                                           
  niter=50
! minimum difference                                                                             
  diff=2.e-5

  EXN=(P/p1000mb)**rcp
  ! assume first that th = thv                                                                   
  T = THV*EXN
  !QS = qsat_blend(T,P)                                                                          
  !QC = QS - QT                                                                                  

  QC=0.

  do ni=1,NITER
     QCOLD = QC
     T = EXN*THV/(1.+QT*(rvovrd-1.)-rvovrd*QC)
     QS=qsat_blend(T,P)
     QC= MAX((QT-QS),0.)
     if (abs(QC-QCOLD)<Diff) exit
  enddo
  THL = (T - xlv/cp*QC)/EXN

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
              &rthraten                              )

        integer, intent(in) :: kts,kte,kpbl
        real(kind_phys), dimension(kts:kte), intent(in) ::            &
            u,v,th,thl,tk,qt,qv,qc,qi,thv,p,qke,rho,exner,            &
            qnc,qni,qnwfa,qnifa,dz,qc_bl1,qi_bl1,cldfra_bl1,el
        real(kind_phys), dimension(kts:kte), intent(in) :: rthraten
        ! zw .. heights of the downdraft levels (edges of boxes)
        real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
        real(kind_phys), intent(in)  :: flt,flq,fltv
        real(kind_phys), intent(in)  :: dt,dx,ust,pblh
  ! outputs - downdraft properties
        real(kind_phys), dimension(kts:kte), intent(inout) ::         &
        edmf_a_dd,   edmf_w_dd,   edmf_qt_dd,  edmf_thl_dd,           &
        edmf_ent_dd, edmf_qc_dd,  tkeprod_dn

  ! outputs - variables needed for solver (sd_aw - sum ai*wi, sd_awphi - sum ai*wi*phii)
        real(kind_phys), dimension(kts:kte+1) ::                      &
            sd_aw, sd_awthl, sd_awu, sd_awv,                          &
            sd_awqt, sd_awqc, sd_awqv, sd_awqi, sd_awqnc, sd_awqni,   &
            sd_awqnwfa, sd_awqnifa, sd_awqke, sd_aw2

   !downdraft properties
        integer, parameter::                                          &
            & ndd    = 5,               &  !number of downdrafts
            & kmin   = 3                   !lowest k-level where downdrafts can start
        real(kind_phys),parameter ::                                  &
            & minddd = 400.,            &  !min downdraft diameter (m)
            & maxddd = 1000,            &  !max downdraft diameter (m)
            & zmin   = 50.,             &  !lowest height where downdrafts can start
            & dz200  = 200.                !depth over which other parameters are normalized to
        real(kind_phys)::                                             &
            & maxdd2,                   &  !variable max downdraft diameter (m)
            &    ddd,                   &  !downdraft diameter
            &     dl,                   &  !diameter increment
            &    add,                   &  !total area of downdrafts
            & minexc,                   &  !minimum init downdraft temp pert
            & maxexc                       !maximum init downdraft temp pert
  ! k-index of downdraft starting height
        integer,         dimension(1:ndd) :: dd_initk
  ! downdraft column properties
        real(kind_phys), dimension(kts:kte+1,1:ndd) ::                &
            downw,downthl,downqt,downqc,downa,downu,downv,downthv,    &
            downqi,downqnc,downqni,downqnwfa,downqnifa
  ! entrainment variables
        real(kind_phys), dimension(kts:kte+1,1:ndd) :: ent
  ! internal variables
        integer :: k,i,ki,qltop,qlbase
        real(kind_phys):: qstar,thstar,sigmaw,sigmaqt,                &
            sigmath,z0,pwmin,pwmax,wmin,wmax,went,mindownw
        real(kind_phys):: qtn,thln,thvn,qcn,un,vn,qken,wn2,wn,        &
            qin,qncn,qnin,qnwfan,qnifan,                              &
            thvk,pk,entexp,entw,beta_dm,entexm,rho_int
        real(kind_phys):: jump_thv,jump_qt,jump_thetal,               &
            refthl,refthv,refthlv,refqt,refqc,refqi,refqnc,refqni,    &
            refqnwfa,refqnifa,refu,refv,refqke,refqt2,reftk,refp,     &
            qx_k,qx_km1,refthvm1,cftop,ac_wsp,wspd_pbl

  ! dd specific internal variables
        real(kind_phys):: radflux, f0, wstar_rad, dz_ent
        logical :: cloudflg
        logical :: singlelayer             !check for single or multi-layer clouds
        real(kind_phys):: sigq,xl,rsl,cpm,a,mf_cf,diffqt,             &
            fng,qww,alpha,beta,bb,f,pt,t,q2p,b9,satvp,rhgrid,         &
            def_th,def_qt,frac_liq,frac_ice,qs1,qs2

  ! w parameters
        real(kind_phys),parameter ::                                  &
            &wa=1., wb=1.5, z00=100.
        real(kind_phys) :: buoy,bcoeff,ecoeff

  ! additional printouts for debugging
        integer, parameter :: debug_dd=0

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
   tkeprod_dn =zero
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

  ! determine maximum downdraft width and increments
   maxdd2     =min(maxddd, 1.2*dx)
   dl         =max(maxdd2-minddd, zero)/real(ndd-1)

  ! first, check for stratocumulus-topped pbl with cooling at the cloud top.
   cloudflg   =.false.
   singlelayer=.false.
   qltop      =1     !initialize
   qlbase     =1     !initialize
   do k = max(kmin+1,kpbl-2),kpbl+10
      qx_k    =max(qt(k)  , qc_bl1(k)  +qi_bl1(k))   !total clouds at k
      qx_km1  =max(qt(k-1), qc_bl1(k-1)+qi_bl1(k-1)) !total clouds at k-1
!      if (qx_k  .lt. 1.e-6 .and. cldfra_bl(k-1).lt.0.1 .and. &
!          qx_km1.gt. 1.e-6 .and. cldfra_bl(k)  .gt.0.5 .and. &
      !---------------------------------criteria for downdraft activation
      if ((cldfra_bl1(k-1).gt.0.5       ) .and. & !1) stratocumulus exists
!          (cldfra_bl1(k)  .lt.0.1       ) .and. & !2) clear above
          (rthraten(k-1) .lt. -0.000025)  .and. & !3) radiative cooling
          (           dl .gt.zero       )) then   !4) widest downdraft > thinnest
         cloudflg =.true. ! found sc cloud
         qltop    =k-1    ! index for sc cloud top
         cftop    =cldfra_bl1(k-1)
         if (cldfra_bl1(k-2).lt.0.1)singlelayer=.true.
      endif
   enddo

   !found sc cloud, trigger downdrafts
   if (cloudflg) then
      do k = qltop, kts, -1
         qx_k    =max(qt(k)  , qc_bl1(k) +qi_bl1(k))   !total clouds at k
         if (qx_k .gt. 1e-6) then
            qlbase = k ! index for sc cloud base
         endif
      enddo

      do i=1,ndd
         ! downdraft starts at the first layer interface below the cloud top
         dd_initk(i) = qltop
      enddo

      ! loop radflux
      f0 = zero
      do k = kmin, qltop
         radflux = rthraten(k) * exner(k) ! converts theta/s to temperature/s
         radflux = radflux * cp / grav * ( p(k) - p(k+1) ) ! converts k/s to w/m^2
         if ( radflux < zero ) f0 = abs(radflux) + f0
      enddo
      f0 = min(max(f0, 1.0_kind_phys), 200._kind_phys)  !total radiative cooling (w/m2)

      !allow the total fractional area of the downdrafts (add) to be proportional
      !to the radiative forcing:
      !for  50 w/m2, add = 0.20
      !for 100 w/m2, add = 0.30
      !for 150 w/m2, add = 0.40
      add = min( 0.1 + f0*0.002, 0.3_kind_phys)
      !taper off area for high wind speeds
      wspd_pbl=SQRT(MAX(u(kts)**2 + v(kts)**2, 0.01_kind_phys))
      if (wspd_pbl .le. 10.) then
         ac_wsp = one
      else
         ac_wsp = one - min((wspd_pbl - 10.0)/15., 1.0_kind_phys)
      endif
      add  = add * ac_wsp

      !find inversion strength across cloud top entrainment zone--normalized to 200 m vertical grid spacing
      dz_ent      = 0.5 * (dz(qltop+1) + dz(qltop))
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
         refthl  = (thl(ki-1)*dz(ki) + thl(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refthlv = refthl*(1.+p608*(qv(ki-1)*dz(ki) + qv(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1)))
         refthv  = (thv(ki-1)*dz(ki) + thv(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refthvm1= (thv(ki-2)*dz(ki-1) + thv(ki-1)*dz(ki-2)) /(dz(ki-1)+dz(ki-2)) 
         reftk   = (tk(ki-1)*dz(ki)  + tk(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refqt   = (qt(ki-1)*dz(ki)  + qt(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refqc   = (qc(ki-1)*dz(ki)  + qc(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refqi   = (qi(ki-1)*dz(ki)  + qi(ki)*dz(ki-1))  /(dz(ki)+dz(ki-1))
         refqnc  = (qnc(ki-1)*dz(ki) + qnc(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refqni  = (qni(ki-1)*dz(ki) + qni(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refqnwfa= (qnwfa(ki-1)*dz(ki) + qnwfa(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refqnifa= (qnifa(ki-1)*dz(ki) + qnifa(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refu    = (u(ki-1)*dz(ki)   + u(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
         refv    = (v(ki-1)*dz(ki)   + v(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
         refqke  = (qke(ki-1)*dz(ki) + qke(ki)*dz(ki-1)) /(dz(ki)+dz(ki-1))
         refp    = (p(ki-1)*dz(ki)   + p(ki)*dz(ki-1))   /(dz(ki)+dz(ki-1))
      endif

      ! w* from radiative forcing
      !wstar_rad = ( grav * dz_ent * f0 / (refthl * rho(ki) * cp) ) **onethird
      wstar_rad = 10.* ( grav*dz200 * f0 / (refthl * rho(ki) * cp) ) **onethird
      wstar_rad = min(max(wstar_rad, 0.1), 3.0)
      ! note: since dz_ent cancels, went is not a function of dz_ent
      !went      = thv(1) / ( grav * jump_thv * dz_ent ) * &
      went      = thv(1) / ( grav * jump_thv * dz200 ) * &
                  (0.5 * wstar_rad**3 )
      qstar     = abs(went*jump_qt/wstar_rad)
      thstar    = f0/rho(ki)/cp/wstar_rad - went*jump_thv/wstar_rad

      sigmaw    = 0.3 * wstar_rad ! 0.8*wstar_rad ! tuning parameter ! 0.5 was good
      sigmaqt   = 40. * qstar  ! 50 was good
      sigmath   = 1.0 * thstar ! 0.5 was good

      pwmin     = -1. ! drawing from the negative tail -3sigma to -1sigma
      pwmax     = -3.
      wmin      = max(min(sigmaw*pwmin, -0.1), -0.2)
      wmax      = max(min(sigmaw*pwmax, -0.2), -0.5)

      !initialize downdraft temperature deficit def_th
      minexc = min(-0.05     , 0.5*(refthvm1 - refthv)) !increase prob to go down 1 level
      maxexc = min(minexc-0.1, -0.5          )
      maxexc = max(maxexc    , thv(1) - refthlv)  !init parcel temp >= thv_sfc
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
         print*,"found conditions for downdraft mixing"
         print*,"qltop=",qltop," qlbase=",qlbase
         print*,"qstar=",qstar," thstar=",thstar," went=",went
         print*,"f0=",f0," jump_thv=",jump_thv
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

      if (debug_dd .eq.1) then
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
            wmin = 0.1 + ddd*0.0005
            !ent(k,i) = 0.33/(min(max(abs(downw(k+1,i)),wmin),0.9)*ddd)
            ent(k,i) = 0.53/(min(max(abs(downw(k+1,i)),wmin),0.6)*ddd)

            !minimum background entrainment and 1/z enhancement near surface
            ent(k,i) = max(ent(k,i),0.0003)
            ent(k,i) = max(ent(k,i),0.06/zw(k))
            ent(k,i) = min(ent(k,i),0.9/dz(k)) ! liberal check for numerical stability

            !starting at the first interface level below cloud top
            !entexp=exp(-ent(k,i)*dz(k))
            !entexp_m=exp(-ent(k,i)/3.*dz(k))
            entexp  =ent(k,i)*dz(k)        !for all scalars
            entexm  =ent(k,i)*dz(k)*0.5    !test for momentum

            qtn     =downqt(k+1,i) *(1.-entexp) + qt(k)*entexp
            thln    =downthl(k+1,i)*(1.-entexp) + thl(k)*entexp
            qncn    =downqnc(k+1,i) *(1.-entexp) + qnc(k)*entexp
            qnin    =downqni(k+1,i) *(1.-entexp) + qni(k)*entexp
            qnwfan  =downqnwfa(k+1,i) *(1.-entexp) + qnwfa(k)*entexp
            qnifan  =downqnifa(k+1,i) *(1.-entexp) + qnifa(k)*entexp
            un      =downu(k+1,i)  *(1.-entexm) + u(k)*entexm
            vn      =downv(k+1,i)  *(1.-entexm) + v(k)*entexm
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
               bcoeff = 0.2*min(zw(k)/pblh, one)
               ecoeff = -2.0 - max(0.5*pblh-zw(k), zero)/(0.5*pblh)
               buoy   = buoy*min(zw(k)/pblh, one)
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

            mindownw = min(downw(k+1,i),-0.2)
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
            wn = max(min(wn,zero), -3.0)

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
               downw(k,i)  = wn !-sqrt(wn2)
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
   endif ! end cloud flag

   downw(1,:) = 0. !make sure downdrafts do not penetrate the surface
   downa(1,:) = 0.

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
      !instead of dTKE/dt = 1/2 w^3, multiply by 2 for QKE.
      tkeprod_dn(k)=(abs(edmf_w_dd(k))**3)*edmf_a_dd(k)/(b1*max(el(k),0.2))
   enddo
   ! add tke source for entrainment at layer above cloud. use same area
   ! above cloud as used in the initialized downdraft area.
   tkeprod_dn(qltop+1)=went*edmf_a_dd(qltop)/(b1*max(el(qltop+1),0.1)) 

   !
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

end subroutine ddmp_mf
!=============================================================== 

SUBROUTINE SCALE_AWARE(dx,pblh,Psig_bl,Psig_shcu)

    !---------------------------------------------------------------
    !             NOTES ON SCALE-AWARE FORMULATION
    !
    !JOE: add scale-aware factor (Psig) here, taken from Honnert et al. (2011,
    !     JAS) and/or from Hyeyum Hailey Shin and Song-You Hong (2013, JAS)
    !
    ! Psig_bl tapers local mixing
    ! Psig_shcu tapers nonlocal mixing

    real(kind_phys), intent(in)  :: dx,pblh
    real(kind_phys), intent(out) :: Psig_bl,Psig_shcu
    real(kind_phys)              :: dxdh

    Psig_bl=one
    Psig_shcu=one
    dxdh=MAX(2.5*dx,10.)/MIN(PBLH,3000.)
    ! Honnert et al. 2011, TKE in PBL  *** original form used until 201605
    !Psig_bl= ((dxdh**2) + 0.07*(dxdh**0.667))/((dxdh**2) + &
    !         (3./21.)*(dxdh**0.67) + (3./42.))
    ! Honnert et al. 2011, TKE in entrainment layer
    !Psig_bl= ((dxdh**2) + (4./21.)*(dxdh**0.667))/((dxdh**2) + &
     !        (3./20.)*(dxdh**0.67) + (7./21.))
    ! New form to preseve parameterized mixing - only down 5% at dx = 750 m
     Psig_bl= ((dxdh**2) + 0.106*(dxdh**0.667))/((dxdh**2) +0.066*(dxdh**0.667) + 0.071)

    !assume a 500 m cloud depth for shallow-cu clods
    dxdh=MAX(2.5*dx,10.)/MIN(PBLH+500.,3500.)
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
    Psig_shcu= ((dxdh**2) + 0.145*(dxdh**0.667))/((dxdh**2) +0.172*(dxdh**0.667) + 0.170)

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

      IMPLICIT NONE
      
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

      XC=MAX(-80.,t - t0c) !note t0c = 273.15, tice is set in module mynn_common to 240

! For 240 < t < 268.16 K, the vapor pressures are "blended" as a function of temperature, 
! using the approach similar to Chaboureau and Bechtold (2002), JAS, p. 2363.  The resulting 
! values are returned from the function.
      IF (t .GE. (t0c-6.)) THEN
          esat_blend = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8))))))) 
      ELSE IF (t .LE. tice) THEN
          esat_blend = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
      ELSE
          ESL = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESI = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          chi = ((t0c-6.) - t)/((t0c-6.) - tice)
          esat_blend = (1.-chi)*ESL  + chi*ESI
      END IF

  END FUNCTION esat_blend

! ====================================================================

!>\ingroup gsd_mynn_edmf
!! This function extends function "esat" and returns a "blended"
!! saturation mixing ratio. Tice currently set to 240 K, t0c = 273.15 K.
!!\author JAYMES
  FUNCTION qsat_blend(t, P)

      IMPLICIT NONE

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

      XC=MAX(-80.,t - t0c)

      IF (t .GE. (t0c-6.)) THEN
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESL  = min(ESL, P*0.15) ! Even with P=1050mb and T=55C, the sat. vap. pres only contributes to ~15% of total pres.
          qsat_blend = 0.622*ESL/max(P-ESL, 1e-5) 
      ELSE IF (t .LE. tice) THEN
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          ESI  = min(ESI, P*0.15)
          qsat_blend = 0.622*ESI/max(P-ESI, 1e-5)
      ELSE
          ESL  = J0+XC*(J1+XC*(J2+XC*(J3+XC*(J4+XC*(J5+XC*(J6+XC*(J7+XC*J8)))))))
          ESL  = min(ESL, P*0.15)
          ESI  = K0+XC*(K1+XC*(K2+XC*(K3+XC*(K4+XC*(K5+XC*(K6+XC*(K7+XC*K8)))))))
          ESI  = min(ESI, P*0.15)
          RSLF = 0.622*ESL/max(P-ESL, 1e-5)
          RSIF = 0.622*ESI/max(P-ESI, 1e-5)
!          chi  = (268.16-t)/(268.16-240.)
          chi  = ((t0c-6.) - t)/((t0c-6.) - tice) 
         qsat_blend = (1.-chi)*RSLF + chi*RSIF
      END IF

  END FUNCTION qsat_blend

! ===================================================================

!>\ingroup gsd_mynn_edmf
!! This function interpolates the latent heats of vaporization and sublimation into
!! a single, temperature-dependent, "blended" value, following 
!! Chaboureau and Bechtold (2002) \cite Chaboureau_2002, Appendix.
!!\author JAYMES
  FUNCTION xl_blend(t)

      IMPLICIT NONE

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
          xl_blend = (1.-chi)*xlvt + chi*xlst     !blended
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
      IMPLICIT NONE

      real(kind_phys), intent(in):: zet
      real(kind_phys):: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      real(kind_phys), parameter :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      real(kind_phys), parameter :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      real(kind_phys), parameter :: am_unst=10., ah_unst=34.
      real(kind_phys):: phi_m,phim

      if ( zet >= zero ) then
         dummy_0=1+zet**bm_st
         dummy_1=zet+dummy_0**(rbm_st)
         dummy_11=1+dummy_0**(rbm_st-1)*zet**(bm_st-1)
         dummy_2=(-am_st/dummy_1)*dummy_11
         phi_m = 1-zet*dummy_2
      else
         dummy_0 = (one-cphm_unst*zet)**0.25
         phi_m = 1./dummy_0
         dummy_psi = 2.*log(0.5*(1.+dummy_0))+log(0.5*(1.+dummy_0**2))-2.*atan(dummy_0)+1.570796

         dummy_0=(1.-am_unst*zet)          ! parentesis arg
         dummy_1=dummy_0**0.333333         ! y
         dummy_11=-0.33333*am_unst*dummy_0**(-0.6666667) ! dy/dzet
         dummy_2 = 0.33333*(dummy_1**2.+dummy_1+1.)    ! f
         dummy_22 = 0.3333*dummy_11*(2.*dummy_1+1.)    ! df/dzet
         dummy_3 = 0.57735*(2.*dummy_1+1.) ! g
         dummy_33 = 1.1547*dummy_11        ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44 = (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(1.+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = 1./(1.+dummy_0) ! denon
         dummy_11 = 2.*zet         ! denon/dzet
         dummy_2 = ((1-phi_m)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22 = -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phi_m = 1.-zet*(dummy_2+dummy_22)
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
      IMPLICIT NONE

      real(kind_phys), intent(in):: zet
      real(kind_phys):: dummy_0,dummy_1,dummy_11,dummy_2,dummy_22,dummy_3,dummy_33,dummy_4,dummy_44,dummy_psi
      real(kind_phys), parameter :: am_st=6.1, bm_st=2.5, rbm_st=1./bm_st
      real(kind_phys), parameter :: ah_st=5.3, bh_st=1.1, rbh_st=1./bh_st
      real(kind_phys), parameter :: am_unst=10., ah_unst=34.
      real(kind_phys):: phh,phih

      if ( zet >= zero ) then
         dummy_0=1+zet**bh_st
         dummy_1=zet+dummy_0**(rbh_st)
         dummy_11=1+dummy_0**(rbh_st-1)*zet**(bh_st-1)
         dummy_2=(-ah_st/dummy_1)*dummy_11
         phih = 1-zet*dummy_2
      else
         dummy_0 = (one-cphh_unst*zet)**0.5
         phh = one/dummy_0
         dummy_psi = 2.*log(0.5*(1.+dummy_0))

         dummy_0=(one-ah_unst*zet)          ! parentesis arg
         dummy_1=dummy_0**0.333333         ! y
         dummy_11=-0.33333*ah_unst*dummy_0**(-0.6666667) ! dy/dzet
         dummy_2 = 0.33333*(dummy_1**2.+dummy_1+one)    ! f
         dummy_22 = 0.3333*dummy_11*(2.*dummy_1+one)    ! df/dzet
         dummy_3 = 0.57735*(2.*dummy_1+one) ! g
         dummy_33 = 1.1547*dummy_11        ! dg/dzet
         dummy_4 = 1.5*log(dummy_2)-1.73205*atan(dummy_3)+1.813799364 !psic
         dummy_44 = (1.5/dummy_2)*dummy_22-1.73205*dummy_33/(1.+dummy_3**2)! dpsic/dzet

         dummy_0 = zet**2
         dummy_1 = one/(one+dummy_0)         ! denon
         dummy_11 = 2.*zet                 ! ddenon/dzet
         dummy_2 = ((1-phh)/zet+dummy_11*dummy_4+dummy_0*dummy_44)*dummy_1
         dummy_22 = -dummy_11*(dummy_psi+dummy_0*dummy_4)*dummy_1**2

         phih = 1.-zet*(dummy_2+dummy_22)
      end if

END FUNCTION phih
! ==================================================================
 SUBROUTINE topdown_cloudrad(kts,kte,                         &
               &dz1,zw,fltv,xland,kpbl,PBLH,                  &
               &sqc,sqi,sqw,thl,th1,ex1,p1,rho1,thv,          &
               &cldfra_bl1,rthraten,                          &
               &maxKHtopdown,KHtopdown,TKEprodTD              )

    !input
    integer,         intent(in) :: kte,kts
    real(kind_phys), dimension(kts:kte), intent(in) :: dz1,sqc,sqi,sqw,&
          thl,th1,ex1,p1,rho1,thv,cldfra_bl1
    real(kind_phys), dimension(kts:kte), intent(in) :: rthraten
    real(kind_phys), dimension(kts:kte+1), intent(in) :: zw
    real(kind_phys), intent(in) :: pblh,fltv
    real(kind_phys), intent(in) :: xland
    integer        , intent(in) :: kpbl
    !output
    real(kind_phys), intent(out) :: maxKHtopdown
    real(kind_phys), dimension(kts:kte), intent(out) :: KHtopdown,TKEprodTD
    !local
    real(kind_phys), dimension(kts:kte) :: zfac,wscalek2,zfacent
    real(kind_phys) :: bfx0,wm2,wm3,bfxpbl,dthvx,tmp1
    real(kind_phys) :: temps,templ,zl1,wstar3_2
    real(kind_phys) :: ent_eff,radsum,radflux,we,rcldb,rvls,minrad,zminrad
    real(kind_phys), parameter :: pfac =2.0, zfmin = 0.01, phifac=8.0
    integer :: k,kk,kminrad
    logical :: cloudflg

    cloudflg=.false.
    minrad=100.
    kminrad=kpbl
    zminrad=PBLH
    KHtopdown(kts:kte)=zero
    TKEprodTD(kts:kte)=zero
    maxKHtopdown=zero

    !CHECK FOR STRATOCUMULUS-TOPPED BOUNDARY LAYERS
    DO kk = MAX(1,kpbl-2),kpbl+3
       if (sqc(kk).gt. 1.e-6 .OR. sqi(kk).gt. 1.e-6 .OR. &
           cldfra_bl1(kk).gt.0.5) then
          cloudflg=.true.
       endif
       if (rthraten(kk) < minrad)then
          minrad=rthraten(kk)
          kminrad=kk
          zminrad=zw(kk) + 0.5*dz1(kk)
       endif
    ENDDO

    IF (MAX(kminrad,kpbl) < 2)cloudflg = .false.
    IF (cloudflg) THEN
       zl1 = dz1(kts)
       k = MAX(kpbl-1, kminrad-1)
       !Best estimate of height of TKE source (top of downdrafts):
       !zminrad = 0.5*pblh(i) + 0.5*zminrad

       templ=thl(k)*ex1(k)
       !rvls is ws at full level
       rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep_2/p1(k+1))
       temps=templ + (sqw(k)-rvls)/(cp/xlv  +  ep_2*xlv*rvls/(r_d*templ**2))
       rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep_2/p1(k+1))
       rcldb=max(sqw(k)-rvls,0.)

       !entrainment efficiency
       dthvx     = (thl(k+2) + th1(k+2)*p608*sqw(k+2)) &
                 - (thl(k)   + th1(k)  *p608*sqw(k))
       dthvx     = max(dthvx,0.1)
       tmp1      = xlvcp * rcldb/(ex1(k)*dthvx)
       !Originally from Nichols and Turton (1986), where a2 = 60, but lowered
       !here to 8, as in Grenier and Bretherton (2001).
       ent_eff   = 0.2 + 0.2*8.*tmp1

       radsum=0.
       DO kk = MAX(1,kpbl-3),kpbl+3
          radflux=rthraten(kk)*ex1(kk)         !converts theta/s to temp/s
          radflux=radflux*cp/grav*(p1(kk)-p1(kk+1)) ! converts temp/s to W/m^2
          if (radflux < zero ) radsum=abs(radflux)+radsum
       ENDDO

       !More strict limits over land to reduce stable-layer mixouts
       if ((xland-1.5).GE.0)THEN      ! WATER
          radsum=MIN(radsum,90.0)
          bfx0 = max(radsum/rho1(k)/cp, zero)
       else                           ! LAND
          radsum=MIN(0.25*radsum,30.0)!practically turn off over land
          bfx0 = max(radsum/rho1(k)/cp - max(fltv,zero), zero)
       endif

       !entrainment from PBL top thermals
       wm3    = grav/thv(k)*bfx0*MIN(pblh,1500.) ! this is wstar3(i)
       wm2    = wm2 + wm3**twothirds
       bfxpbl = - ent_eff * bfx0
       dthvx  = max(thv(k+1)-thv(k),0.1)
       we     = max(bfxpbl/dthvx,-sqrt(wm3**twothirds))

       DO kk = kts,kpbl+3
          !Analytic vertical profile
          zfac(kk) = min(max((1.-(zw(kk+1)-zl1)/(zminrad-zl1)),zfmin),1.)
          zfacent(kk) = 10.*MAX((zminrad-zw(kk+1))/zminrad,zero)*(1.-zfac(kk))**3

          !Calculate an eddy diffusivity profile (not used at the moment)
          wscalek2(kk) = (phifac*karman*wm3*(zfac(kk)))**onethird
          !Modify shape of Kh to be similar to Lock et al (2000): use pfac = 3.0
          KHtopdown(kk) = wscalek2(kk)*karman*(zminrad-zw(kk+1))*(1.-zfac(kk))**3 !pfac
          KHtopdown(kk) = MAX(KHtopdown(kk),zero)

          !Calculate TKE production = 2(g/TH)(w'TH'), where w'TH' = A(TH/g)wstar^3/PBLH,
          !A = ent_eff, and wstar is associated with the radiative cooling at top of PBL.
          !An analytic profile controls the magnitude of this TKE prod in the vertical.
          TKEprodTD(kk)=2.*ent_eff*wm3/MAX(pblh,100.)*zfacent(kk)
          TKEprodTD(kk)= MAX(TKEprodTD(kk),zero)
       ENDDO
    ENDIF !end cloud check
    maxKHtopdown=MAXVAL(KHtopdown(:))

 END SUBROUTINE topdown_cloudrad
! ==================================================================
! ===================================================================
! ===================================================================

END MODULE module_bl_mynn
