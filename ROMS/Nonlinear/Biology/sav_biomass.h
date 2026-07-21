      SUBROUTINE sav_biomass_sub(ng, Istr, Iend, LBi, UBi,              &
     &                        pmonth, wtemp,                            &
     &                        PARz, DINwcr_loc,                         &
#ifdef WET_DRY
     &                        rmask_wet_loc,                            &
#endif 
#if defined VEGETATION && defined VEG_BIOMASS
     &                        pdens_loc, phght_loc, pdiam_loc,          & 
#endif 
     &                        DINsed_loc, DINwcr_sav_loc, DOwcr_loc,    &
     &                        CO2wcr_loc, LDeCwcr_loc, LDeNwcr_loc,     &
     &                        agb_loc, bgb_loc, epb_loc,                &
     &                        pp_loc, agm_loc, agar_loc,                &
     &                        agbr_loc, sears_loc, agbg_loc,            &
     &                        bgag_loc, bgr_loc, bgm_loc)
!                                                                      ! 
!**********************************************************************!
!****************************************** John C. Warner ************!
!****************************************** Neil K. Ganju *************!
!****************************************** Jeremy Testa **************!
!****************************************** Tarandeep S. Kalra ********!
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     t          Water temperature (Celsius).                          !
!     PArout     PAR at depth (mu E m-2 s-1).                          !
!     DINwcr_loc Dissolved Inorganic N in water col. (mu M)            !
!     DINsed_loc Dissolved Inorganic N in sediment col. (mu M)         !
!#ifdef WET_DRY                                                        !
!     rmask_wet_loc wetting-drying mask                                !
!#endif                                                                !
!#ifdef VEG_BIOMASS                                                    !
!     pdens_loc  Vector of SAV density (stems m-2)                     !
!     phght_loc  Vector of SAV height (m)                              !
!     pdiam_loc  Vector of SAV diameter(width) (m)                     !
!#endif 
!     DINwcr_sav_loc Change in dissolved Inorganic N in water col.     !
!                    due to SAV model (mu M)                           !
!     DOwcr_loc   O2 interaction with bed                              !
!     CO2wcr_loc  CO2 interactions with bed                            !
!     LdeCwcr_loc Labile detrital carbon in fennel                     !
!     LdeNwcr_loc Labile detrital carbon in fennel                     !
!     agb_loc     Vector of above ground biomass  (mmol N m-2)         !
!     bgb_loc     Vector of below ground biomass  (mmol N m-2)         !
!     epb_loc     Vector of epiphyte biomass      (mmol N m-2)         !
!                                                                      !
!  Estuarine SAV Model developed by                                    !
!  Jeremy Testa, May 2015, Chesapeake Biological Laboratory            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!      Kalra, T. S., Ganju, N. K., and Testa, J. M.: Development of a  !
!      Submerged Aquatic Vegetation Growth Model in a Coupled          !
!      Wave-Current-Sediment-Transport Modeling System (COAWST v3.5),  !
!      Geosci. Model Dev. Discuss. (in review 2019).                   !
!                                                                      !
!      Madden, C. J., Kemp, W. M. June 1996: Ecosystem Model of an     !
!      Estuarine Submersed Plant Community: Calibration and Simulation !
!      of Eutrophication Responses: Estuarine Research Foundation      !
!      Vol. 19, No. 2B, p. 457-474                                     !
!                                                                      !
!      Cerco. C.F. and K. Moore. 2001. System-wide submerged aquatic   !
!      vegetation model for Chesapeake Bay. Estuaries 24: 522-534      !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
      USE mod_biology
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi
      integer, intent(in) :: Istr, Iend
      real(r8), intent(in) :: pmonth
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: wtemp(LBi:)
      real(r8), intent(in) :: PARz(LBi:)
      real(r8), intent(inout) :: DINwcr_loc(LBi:)
      real(r8), intent(inout) :: DINsed_loc(LBi:)
#  ifdef WET_DRY
      real(r8), intent(in)    :: rmask_wet_loc(LBi:)
#  endif 
#  if defined VEGETATION && defined VEG_BIOMASS
      real(r8), intent(in)    :: pdiam_loc(LBi:)
      real(r8), intent(inout) :: pdens_loc(LBi:)
      real(r8), intent(inout) :: phght_loc(LBi:)
#  endif 
      real(r8), intent(inout) :: DINwcr_sav_loc(LBi:)
      real(r8), intent(inout) :: DOwcr_loc(LBi:)
      real(r8), intent(inout) :: CO2wcr_loc(LBi:)
      real(r8), intent(inout) :: LdeCwcr_loc(LBi:)
      real(r8), intent(inout) :: LdeNwcr_loc(LBi:)
      real(r8), intent(inout) :: agb_loc(LBi:)
      real(r8), intent(inout) :: bgb_loc(LBi:)
      real(r8), intent(inout) :: epb_loc(LBi:)
      real(r8), intent(inout) :: pp_loc(LBi:)
      real(r8), intent(inout) :: agm_loc(LBi:)
      real(r8), intent(inout) :: agar_loc(LBi:)
      real(r8), intent(inout) :: agbr_loc(LBi:)
      real(r8), intent(inout) :: sears_loc(LBi:)
      real(r8), intent(inout) :: agbg_loc(LBi:)
      real(r8), intent(inout) :: bgag_loc(LBi:)
      real(r8), intent(inout) :: bgr_loc(LBi:)
      real(r8), intent(inout) :: bgm_loc(LBi:)
#  else
      real(r8), intent(in) :: wtemp(LBi:UBi)
      real(r8), intent(in) :: PARz(LBi:UBi)
      real(r8), intent(inout) :: DINwcr_loc(LBi:UBi)
      real(r8), intent(inout) :: DINsed_loc(LBi:UBi)
#  ifdef WET_DRY
      real(r8), intent(in)    :: rmask_wet_loc(LBi:UBi)
#  endif 
#  if defined VEGETATION && defined VEG_BIOMASS
      real(r8), intent(in)    :: pdiam_loc(LBi:UBi)
      real(r8), intent(inout) :: pdens_loc(LBi:UBi)
      real(r8), intent(inout) :: phght_loc(LBi:UBi)
#  endif 
      real(r8), intent(inout) :: DINwcr_sav_loc(LBi:UBi)
      real(r8), intent(inout) :: DOwcr_loc(LBi:UBi)
      real(r8), intent(inout) :: CO2wcr_loc(LBi:UBi)
      real(r8), intent(inout) :: LdeCwcr_loc(LBi:UBi)
      real(r8), intent(inout) :: LdeNwcr_loc(LBi:UBi)
      real(r8), intent(inout) :: agb_loc(LBi:UBi)
      real(r8), intent(inout) :: bgb_loc(LBi:UBi)
      real(r8), intent(inout) :: epb_loc(LBi:UBi)
      real(r8), intent(inout) :: pp_loc(LBi:UBi)
      real(r8), intent(inout) :: agm_loc(LBi:UBi)
      real(r8), intent(inout) :: agar_loc(LBi:UBi)
      real(r8), intent(inout) :: agbr_loc(LBi:UBi)
      real(r8), intent(inout) :: sears_loc(LBi:UBi)
      real(r8), intent(inout) :: agbg_loc(LBi:UBi)
      real(r8), intent(inout) :: bgag_loc(LBi:UBi)
      real(r8), intent(inout) :: bgr_loc(LBi:UBi)
      real(r8), intent(inout) :: bgm_loc(LBi:UBi)
#  endif
!
!-------------------------------------------------------
!  Local variable declarations.
!-------------------------------------------------------
!
      integer :: i, thresh, thresh2
      integer, parameter  :: one=1
      real(r8) :: ua, day_year
      real(r8) :: llim, knt, nlim, sed_frc
      real(r8) :: lmba
      real(r8) :: ppbm
!     real(r8) :: pp, agm, agar
!     real(r8) :: agbr, sears, agbg, bgag
!     real(r8) :: bgr, bgm
      real(r8) :: cff,cff1
      real(r8) :: temp
      real(r8) :: dtdays
!
      real(r8) :: pdens, phght
      real(r8), parameter :: pdiam_const=0.003_r8
      real(r8), parameter :: cm_m=0.01_r8
#  ifdef SAV_BIOMASS_DAILY_UPDATE 
      real(r8) :: time_loc
      real(r8) :: freq1, freq2
#  endif 
!
      real(r8), parameter :: Inival=0.0_r8
      real(r8), parameter :: eps=1.0e-10_r8
      real(r8), parameter :: gr2mmol=1000.0_r8/14.007_r8
      real(r8), parameter :: C2N_ratio=30.0_r8
      real(r8), parameter :: gr2mmolC=1000.0_r8/12.011_r8
      real(r8), parameter :: molNmolC=(14.007_r8/12.011_r8)*30.0_r8
      real(r8), parameter :: pqrq=1.0_r8
      real(r8), parameter :: molNmgC_epb=6.625_r8
      real(r8), parameter :: mx_frc=0.2_r8
      real(r8), parameter :: ks_frc=15.0_r8
!
!-------------------------------------------------------
!  Epiphyte Model
!-------------------------------------------------------
!
     real(r8) :: ua_epb, nlim_epb, llim_epb
     real(r8) :: lmba_epb,  pp_epb
     real(r8) :: bresp_epb, aresp_epb
     real(r8) :: grz_epb,   mort_epb
     real(r8) :: cffepb,   cff1epb 
     real(r8) :: frc_slgh, epb_slgh
     real(r8) :: PARepb,   blade_area
!
!-------------------------------------------------------
! Initialize some arrays at each time 
!-------------------------------------------------------
!
      DO i=Istr, Iend
        DINwcr_sav_loc(i)=Inival
        DOwcr_loc(i)=Inival
        CO2wcr_loc(i)=Inival
        LDeNwcr_loc(i)=Inival
        LDeCwcr_loc(i)=Inival
      END DO
!
#  ifdef SAV_BIOMASS_DAILY_UPDATE 
!
!-------------------------------------------------------
! To only update the vegetation growth every 24 hours
!-------------------------------------------------------
!
      time_loc=dt(ng)*iic(ng)
      freq1=MOD(time_loc,day2sec)
      freq2=1.0_r8-CEILING(freq1/day2sec)
#  endif
!
!-----------------------------------------------------------------------
!     pmonth contains the month information calculated in estuarybgc.h
!     compute more constants 
!-----------------------------------------------------------------------
!
      day_year = (pmonth - 52.0_r8)*365.0_r8
      dtdays=dt(ng)*sec2day
      knt=knwc(ng)/knsed(ng)
!       
!     Use temp to convert from grams carbon to mmol nitrogen if needed
!     Taran commented these 
!     temp=gr2mmol/C2N_ratio
!
!---------------------------------------------------------------------
!     Enter the loop to get SAV Biomass computations 
!---------------------------------------------------------------------
! 
      DO i=Istr, Iend
        IF (DINsed_loc(i).lt.1.0_r8) THEN
          DINsed_loc(i)=1.0_r8
        END IF
!
!---------------------------------------------------------------------
!     Initialize local variables and arrays
!---------------------------------------------------------------------
!
!--------------------------------------------------------
! Taran commented these, may be laterdelete
        lmba  = 1.0_r8
        lmba_epb   = 1.0_r8
        blade_area = 1.0_r8
        frc_slgh   = 1.0_r8
        epb_slgh   = 1.0_r8
!
!---------------------------------------------------------------------
!     AGB, BGB are multiplied by wetting-drying masks   
!---------------------------------------------------------------------
!
#  ifdef WET_DRY
        agb_loc(i)=agb_loc(i)*rmask_wet_loc(i)
        bgb_loc(i)=bgb_loc(i)*rmask_wet_loc(i)
#  endif
!
#  ifdef SAV_BIOMASS_DAILY_UPDATE
!
!---------------------------------------------------------------------
!       Compute SAV density, height, leaf area as below
!       Only update the SAV density and height every 24 hours
!---------------------------------------------------------------------
!
        IF(freq2.eq.1.0_r8) THEN
          pdens=4.4554_r8*agb_loc(i)
          phght=cm_m*45.0_r8*(agb_loc(i)/(120.0_r8+agb_loc(i)))
        ENDIF 
#  else
!---------------------------------------------------------------------
!       Update SAV density, height at model time step
!---------------------------------------------------------------------
!
        pdens=4.4554_r8*agb_loc(i)
        phght=cm_m*45.0_r8*(agb_loc(i)/(120.0_r8+agb_loc(i)))
#  endif
!
#  if defined VEGETATION && defined VEG_BIOMASS
!
!---------------------------------------------------------------------
!       If VEG feedback is on, use the user input diameter/width
!---------------------------------------------------------------------
!
        blade_area=0.02_r8*pdens*phght*pdiam_loc(i)
#  else
!
!---------------------------------------------------------------------
!       If VEG feedback is off, use a constant SAV width 
!---------------------------------------------------------------------
!
        blade_area=0.02_r8*pdens*phght*pdiam_const
#  endif 
!
!---------------------------------------------------------------------
!       Epiphyte Model
!       Add epiphyte effects on light available to leaves
!       Compute epiphyte density in g C/m2 SAV leaves
!       epb_lf = epb_loc(i)/(blade_area*10000.0_r8)*molNmgC_epb
!---------------------------------------------------------------------
!
        PARepb = PARz(i)*EXP(-0.42_r8*(2.5_r8*epb_loc(i)))
!
        llim = PARepb/(PARepb+ki(ng))
        nlim=(DINwcr_loc(i)+(knt*DINsed_loc(i)))/                        &
     &       (knwc(ng)+(DINwcr_loc(i)+(knt*DINsed_loc(i))))
!
!-----------------------------------------------------------------------
!    SAV Growth Rate, including self-shading effects
!-----------------------------------------------------------------------
!
        lmba=1.0_r8-(agb_loc(i)/lmbmax(ng))**2
!
        IF(GMODopt(ng).eq.1)THEN
          ua=lmba*scl(ng)*(EXP(thta(ng)*(wtemp(i)-Topt(ng)))-1.0_r8)
        ELSE
          ua=lmba*scl2(ng)*thta2(ng)**(wtemp(i)-Topt(ng))
        ENDIF
!
!-----------------------------------------------------------------------
!       Compute seasonal thresholds for N translocation
!-----------------------------------------------------------------------
! 
        IF ( wtemp(i).gt.Tcrit(ng) ) THEN
          thresh=1.0_r8
        ELSE
          thresh=0.0_r8
        ENDIF
!
        IF ( day_year.gt.273.0_r8 ) THEN
          thresh2=1.0_r8
        ELSE
          thresh2=0.0_r8
        ENDIF
!
!---------------------------------------------------------------------
!    Epiphyte Model, Jeremy Testa and Amanda Moore, Apr-June 2017
!    Model Includes Growth and mortality terms (grazing and sloughing)
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!       Compute light and nutrient limitation terms
!---------------------------------------------------------------------
!
        llim_epb=PARz(i)/(PARz(i)+kl_epb(ng))
        nlim_epb=DINwcr_loc(i)/(kn_epb(ng)+DINwcr_loc(i))
!
!---------------------------------------------------------------------
!       Maximum Growth
!---------------------------------------------------------------------
!
        lmba_epb = 1.0_r8-( (epb_loc(i)/lmbmax_epb(ng))**2 )
        IF(GMODopt(ng).eq.1)THEN
          ua_epb=lmba_epb*scl2_epb(ng)*( EXP(arc_epb(ng)*               &
     &	                          (wtemp(i)-Topt_epb(ng)))-1.0_r8 )
        ELSE
          ua_epb=lmba_epb*scl2_epb(ng)*thta2_epb(ng)**( wtemp(i)-       &
     &                                                 Topt_epb(ng) )
        ENDIF
!
!---------------------------------------------------------------------
!       Realized Growth
!---------------------------------------------------------------------
!
        pp_epb=epb_loc(i)*(ua_epb*MIN(nlim_epb,llim_epb))
!
!---------------------------------------------------------------------
!       Active Respiration
!---------------------------------------------------------------------
!
        aresp_epb=pp_epb*arsc_epb(ng)*EXP( arc_epb(ng)*wtemp(i) )
!
!---------------------------------------------------------------------
!       Basal Respiration
!---------------------------------------------------------------------
!
        bresp_epb=epb_loc(i)*bsrc_epb(ng)*EXP( brc_epb(ng)*wtemp(i) )
!
!---------------------------------------------------------------------
!       Mortality rate if no sloughing
!---------------------------------------------------------------------
!
        mort_epb=epb_loc(i)*kmort_epb(ng)
!
!---------------------------------------------------------------------
!       Grazing
!---------------------------------------------------------------------
!
        grz_epb=grzmx_epb(ng)*( 1.0_r8-EXP(-grzk_epb(ng)*epb_loc(i)) )
!
!-----------------------------------------------------------------------
!  Update Epiphyte Biomass  (mmol N m-2)
!-----------------------------------------------------------------------
!
        cff1epb=(pp_epb-aresp_epb-bresp_epb-mort_epb-grz_epb)
!
#  ifdef WET_DRY
        cff1epb=cff1epb*rmask_wet_loc(i) 
#  endif 
        epb_loc(i)=epb_loc(i)+cff1epb*dtdays
!
        IF( agb_loc(i).lt.0.1_r8 ) THEN
          epb_loc(i)=0.0001_r8
        ELSE
          epb_loc(i)=epb_loc(i)
        ENDIF
!
!-----------------------------------------------------------------------
!  SAV Primary production rate
!-----------------------------------------------------------------------
!
        pp_loc(i)=agb_loc(i)*(ua*MIN(llim,nlim))
!
!-----------------------------------------------------------------------
!  Above ground mortality and above ground active respiration
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!       Original quadratic sloughing term
!-----------------------------------------------------------------------
!
!       agm_loc(i)=(kmag(ng)*agb_loc(i))**2.0_r8
!
!-----------------------------------------------------------------------
!       Linear sloughing term
!-----------------------------------------------------------------------
!
        agm_loc(i)=kmag(ng)*agb_loc(i)
!
!-----------------------------------------------------------------------
!       Active respiration
!-----------------------------------------------------------------------
!
        agar_loc(i)=pp_loc(i)*arsc(ng)*EXP(arc(ng)*wtemp(i))
!
!-----------------------------------------------------------------------
!  Above ground basal respiration and seasonal root storage of carbon
!-----------------------------------------------------------------------
!
        agbr_loc(i)=agb_loc(i)*(bsrc(ng)*EXP(brc(ng)*wtemp(i)))
        sears_loc(i)=agb_loc(i)*RtSttl(ng)*thresh2
!
!-----------------------------------------------------------------------
!  Translocation of above ground biomass to below ground
!-----------------------------------------------------------------------
!
        agbg_loc(i)=pp_loc(i)*DOWNt(ng)
!
!-----------------------------------------------------------------------
!  Translocation of below ground biomass to above ground
!-----------------------------------------------------------------------
!
        bgag_loc(i)=bgb_loc(i)*trns(ng)*thresh
!
!-----------------------------------------------------------------------
!  Below ground biomass respiration and Below ground biomass mortality
!-----------------------------------------------------------------------
!
        bgr_loc(i)=bgb_loc(i)*bsrc(ng)*EXP(brc(ng)*wtemp(i))
        bgm_loc(i)=bgb_loc(i)*(0.01_r8*EXP(km(ng)*wtemp(i)))
!
!-----------------------------------------------------------------------
!       If using carbon units, multiply this term by temp
!-----------------------------------------------------------------------
!
        DINsed_loc(i)=DINsed_loc(i)+(bgr_loc(i)+bgm_loc(i))*dtdays
!
!-----------------------------------------------------------------------
!  Compute new AGB biomass (mmol N m-2)
!-----------------------------------------------------------------------
!
        cff=pp_loc(i)*dtdays
        cffepb=pp_epb*dtdays
!
!-----------------------------------------------------------------------
!  If pp>(available nutrients), no growth
!  assumes no new growth if growth exceeds available nutrients
!  April 2017 - remove this limitation, let SAV take up sediment nutrients
!  if no wq nutrients
!-----------------------------------------------------------------------
!
        IF(cff.lt.DINwcr_loc(i))THEN     !This was origonally 'gt'
          cff1=pp_loc(i)+bgag_loc(i)-agm_loc(i)-agar_loc(i)             &
     &                  -agbr_loc(i)-sears_loc(i)-agbg_loc(i)
        ELSE
          cff1=bgag_loc(i)-agm_loc(i)-agar_loc(i)-agbr_loc(i)           &
     &                                    -sears_loc(i)-agbg_loc(i)
        ENDIF
!
!-----------------------------------------------------------------------
!       Remove epiphytes with seagrass sloughing loss
!-----------------------------------------------------------------------
!       SAV fraction sloughed
!-----------------------------------------------------------------------
!
        frc_slgh=(agm_loc(i)*dtdays)/(agb_loc(i)+eps)  
!
!-----------------------------------------------------------------------
!       EPB loss proportional
!-----------------------------------------------------------------------
!
        epb_slgh=epb_loc(i)*frc_slgh              
!
!-----------------------------------------------------------------------
!       Update EPB Biomass
!-----------------------------------------------------------------------
!
        IF(epb_slgh.gt.0.0_r8)THEN 
           epb_loc(i)=epb_loc(i)-epb_slgh  
        ENDIF
!
!-----------------------------------------------------------------------
!       If SAV mortality exceeds biomass, epiphytes are kept > 0
!-----------------------------------------------------------------------
!
        epb_loc(i)=MAX(epb_loc(i),0.0001_r8)
!
!-----------------------------------------------------------------------
!  Update AGB Biomass  (mmol N m-2)
!-----------------------------------------------------------------------
!
!#  ifdef WET_DRY
!        cff1=cff1*rmask_wet_loc(i)
!#  endif 
        agb_loc(i)=agb_loc(i)+cff1*dtdays
!
#  ifdef SAV_BIOMASS_DAILY_UPDATE
!
!-----------------------------------------------------------------------
!  Compute stems per square meter (mean of 6 Chinco cores = 1100)
!  Formulation based upon  Krause-Jensen et al 2000 and lit synthesis
!  Based on the updated agb_loc
!-----------------------------------------------------------------------
!  Empirical computation of shoot length from above ground biomass (cm)
!  Update at 24 hour frequency
!-----------------------------------------------------------------------
!
        IF(freq2.eq.1.0_r8) THEN
!
!---------------------------------------------------------------------
!       Only update the SAV density and height every 24 hours
!---------------------------------------------------------------------
!
          pdens=4.4554_r8*agb_loc(i)
          phght=cm_m*45.0_r8*(agb_loc(i)/(120.0_r8+agb_loc(i)))
!
!---------------------------------------------------------------------
!       Save updated SAV density, height to feedback main routines 
!---------------------------------------------------------------------
!
          pdens_loc(i)=pdens
          phght_loc(i)=phght
        ENDIF
#  else
!
!---------------------------------------------------------------------
!  Update at model time step 
!---------------------------------------------------------------------
!
        pdens=4.4554_r8*agb_loc(i)
        phght=cm_m*45.0_r8*(agb_loc(i)/(120.0_r8+agb_loc(i)))
!
!---------------------------------------------------------------------
!       Save updated SAV density, height to feedback main routines 
!---------------------------------------------------------------------
!
        pdens_loc(i)=pdens
        phght_loc(i)=phght
#  endif
!
#  if defined VEGETATION && defined VEG_BIOMASS
!
!---------------------------------------------------------------------
!       Updated blade area  
!       If VEG feedback is on, use the user input SAV diameter/width 
!---------------------------------------------------------------------
!
        blade_area=0.02_r8*pdens*phght*pdiam_loc(i)
!
#  else
!
!---------------------------------------------------------------------
!       Updated blade area  
!       If VEG feedback is off, use a constant SAV width 
!---------------------------------------------------------------------
!
        blade_area=0.02_r8*pdens*phght*pdiam_const
#  endif 
!
!-----------------------------------------------------------------------
!  Updating Nitrogen in water column with N uptake by plant and N
!  released from SAV respiration and mortality
!  (temp--> converts gram Carbon units to mmol Nitrogen units)
!  If using carbon units, multiply this term by temp
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  Partition nutrient uptake between sediment and water-column
!  via logistic function.
!-----------------------------------------------------------------------
!
        cff1=1.0_r8/(1.0_r8+exp(-mx_frc*(DINwcr_loc(i)-ks_frc)))
        sed_frc=1.0_r8-cff1
!
!-----------------------------------------------------------------------
!       Add seagrass effects first
!-----------------------------------------------------------------------
!
        IF(cff.gt.DINwcr_loc(i))THEN
          DINwcr_sav_loc(i)=DINwcr_sav_loc(i)+(agar_loc(i)+             &
     &                                      agbr_loc(i))*dtdays*        &
     &                                      (1.0_r8-sed_frc)
        ELSE
          DINwcr_sav_loc(i)=DINwcr_sav_loc(i)+(agar_loc(i)+agbr_loc(i)- &
     &                                         pp_loc(i))*dtdays*       &
     &                                         (1.0_r8-sed_frc) 
        ENDIF
!
!----------------------------------------------------------------------
!       SAV sediment N uptake
!----------------------------------------------------------------------
!
        cff1=(agar_loc(i)+agbr_loc(i)-pp_loc(i))*dtdays*sed_frc
#  ifdef WET_DRY
        cff1=cff1*rmask_wet_loc(i)
#  endif
        DINsed_loc(i)=DINsed_loc(i)+cff1
!
!-----------------------------------------------------------------------
!       Next include epiphyte effects second
!-----------------------------------------------------------------------
!
        IF(cff.gt.DINwcr_loc(i))THEN
          DINwcr_sav_loc(i)=DINwcr_sav_loc(i)+(aresp_epb+               &
     &                                          bresp_epb)*dtdays
        ELSE
          DINwcr_sav_loc(i)=DINwcr_sav_loc(i)+(aresp_epb+bresp_epb-     &
     &                                          pp_epb)*dtdays
        ENDIF
#  ifdef WET_DRY
        DINwcr_sav_loc(i)=DINwcr_sav_loc(i)*rmask_wet_loc(i)
#  endif
!
!-----------------------------------------------------------------------
!  Compute updated BGB Biomass  (mmol N m-2)
!-----------------------------------------------------------------------
!
        cff1=sears_loc(i)+agbg_loc(i)-bgag_loc(i)-bgm_loc(i)-bgr_loc(i)
!#  ifdef WET_DRY
!        bgb_loc(i)=cff1*rmask_wet_loc(i)
!#  endif
        bgb_loc(i)=bgb_loc(i)+cff1*dtdays
!
!-----------------------------------------------------------------------
!  Compute Primary Production per unit biomass
!-----------------------------------------------------------------------
!
        ppbm=(pp_loc(i)/(agb_loc(i)+eps))
!
!-----------------------------------------------------------------------
!  O2 and CO2 interactions with bed
!  If using carbon units, multiply DOwcr_loc and CO2wcr_loc by gr2mmolC
!  (instead of molNmolC)
!-----------------------------------------------------------------------
!     Seagrass  Oxygen
!-----------------------------------------------------------------------
!
        DOwcr_loc(i)=DOwcr_loc(i)+(pp_loc(i)-agar_loc(i)                &
     &                           -agbr_loc(i))*molNmolC*pqrq*dtdays
!
!-----------------------------------------------------------------------
!     Epiphytes Oxygen
!-----------------------------------------------------------------------
!
        DOwcr_loc(i)=DOwcr_loc(i)+(pp_epb-aresp_epb-bresp_epb)*         &
     &                               molNmgC_epb*pqrq*dtdays
#  ifdef WET_DRY
        DOwcr_loc(i)=DOwcr_loc(i)*rmask_wet_loc(i)
#  endif
!
!-----------------------------------------------------------------------
!     Seagrass CO2
!-----------------------------------------------------------------------
!
        CO2wcr_loc(i)=CO2wcr_loc(i)+(agar_loc(i)+agbr_loc(i)-pp_loc(i)) &
     &                             *molNmolC*pqrq*dtdays
!
!-----------------------------------------------------------------------
!     Epiphyte CO2
!-----------------------------------------------------------------------
!
        CO2wcr_loc(i)=CO2wcr_loc(i)+(aresp_epb+bresp_epb-pp_epb)        &
     &                             *molNmgC_epb*pqrq*dtdays
#  ifdef WET_DRY
        CO2wcr_loc(i)=CO2wcr_loc(i)*rmask_wet_loc(i)
#  endif
!
!-----------------------------------------------------------------------
!  Labile detrital carbon and nitrogen interactions with fennel
!  If using carbon units, multiply LDeNwcr_loc by temp and
!  LDeCwcr_loc by gr2mmolC (instead of molNmolC)
!  Both seagrass and epiphytes are included
!-----------------------------------------------------------------------
!
        LDeNwcr_loc(i)=LDeNwcr_loc(i)+(agm_loc(i))*dtdays +             &
     &                  (mort_epb+grz_epb)*dtdays  
#  ifdef WET_DRY
        LDeNwcr_loc(i)=LDeNwcr_loc(i)*rmask_wet_loc(i)
#  endif
!
        LDeCwcr_loc(i)=LDeCwcr_loc(i)+(agm_loc(i))*molNmolC*dtdays +    &
     &                      (mort_epb+grz_epb)*molNmgC_epb*dtdays
#  ifdef WET_DRY
        LDeCwcr_loc(i)=LDeCwcr_loc(i)*rmask_wet_loc(i)
#  endif
!
      END DO
!
      RETURN
      END SUBROUTINE SAV_BIOMASS_SUB
