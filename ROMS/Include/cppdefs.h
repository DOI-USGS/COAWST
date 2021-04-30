/*
** Include file "cppdefs.h"
**
** svn $Id: cppdefs.h 1054 2021-03-06 19:47:12Z arango $
********************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**                                                                           **
** The following is short description of all available CPP options.          **
**                                                                           **
** OPTIONS associated with momentum equations:                               **
**                                                                           **
**   The default horizontal advection is 3rd-order upstream bias for         **
**   3D momentum and 4th-order centered for 2D momentum. The default         **
**   vertical advection is 4th-order centered for 3D momentum. If this       **
**   is the case, no flags for momentum advection need to be activated.      **
**                                                                           **
**   The 3rd-order upstream split advection (UV_U3ADV_SPLIT) can be used     **
**   to correct for the spurious mixing of the advection operator in         **
**   terrain-following coordinates. If this is the case, the advection       **
**   operator is split in advective and viscosity components and several     **
**   internal flags are activated in "globaldefs.h".  Notice that            **
**   horizontal and vertical advection of momentum is 4th-order centered     **
**   plus biharmonic viscosity to correct for spurious mixing. The total     **
**   time-dependent horizontal mixing coefficient are computed in            **
**   "hmixing.F".                                                            **
**                                                                           **
**   WARNING:  Use the splines vertical advection option (UV_SADVECTION)     **
**             only in idealized, high vertical resolution applications.     **
**                                                                           **
** UV_ADV                  to turn ON or OFF advection terms                 **
** UV_COR                  to turn ON or OFF Coriolis term                   **
** UV_U3ADV_SPLIT          if 3rd-order upstream split momentum advection    **
** UV_C2ADVECTION          to turn ON or OFF 2nd-order centered advection    **
** UV_C4ADVECTION          to turn ON or OFF 4th-order centered advection    **
** UV_SADVECTION           to turn ON or OFF splines vertical advection      **
** UV_VIS2                 to turn ON or OFF harmonic horizontal mixing      **
** UV_VIS4                 to turn ON or OFF biharmonic horizontal mixing    **
** UV_SMAGORINSKY          to turn ON or OFF Smagorinsky-like viscosity      **
** UV_DRAG_GRID            if spatially varying bottom friction parameters   **
** UV_LOGDRAG              to turn ON or OFF logarithmic bottom friction     **
** UV_LDRAG                to turn ON or OFF linear bottom friction          **
** UV_QDRAG                to turn ON or OFF quadratic bottom friction       **
** UV_WAVEDRAG             to turn ON or OFF extra linear bottom wave drag   **
** SPLINES_VVISC           if splines reconstruction of vertical viscosity   **
**                                                                           **
** OPTIONS associated with tracers equations:                                **
**                                                                           **
** OPTIONS associated with tracers equations:                                **
**                                                                           **
** TS_DIF2                 to turn ON or OFF harmonic horizontal mixing      **
** TS_DIF4                 to turn ON or OFF biharmonic horizontal mixing    **
** TS_SMAGORINSKY          to turn ON or OFF Smagorinsky-like diffusion      **
** TS_FIXED                if diagnostic run, no evolution of tracers        **
** T_PASSIVE               if inert passive tracers (dyes, etc)              **
** AGE_MEAN                if computing Mean Age of inert passive tracers    **
** NONLIN_EOS              if using nonlinear equation of state              **
** QCORRECTION             if net heat flux correction                       **
** SALINITY                if having salinity                                **
** SCORRECTION             if freshwater flux correction                     **
** SSSC_THRESHOLD          limit on freshwater flux correction               **
** SOLAR_SOURCE            if solar radiation source term                    **
** SPLINES_VDIFF           if splines reconstruction of vertical diffusion   **
** SRELAXATION             if salinity relaxation as a freshwater flux       **
** TRC_PSOURCE             if source of inert passive tracers (dyes, etc)    **
** ONE_TRACER_SOURCE       if one value per tracer for all sources           **
** TWO_D_TRACER_SOURCE     if one value per tracer per source                **
** WTYPE_GRID              to turn ON spatially varying Jerlov water type    **
**                                                                           **
** OPTION to suppress further surface cooling if the SST is at freezing      **
** point or below and the net surface heat flux is cooling:                  **
**                                                                           **
** LIMIT_STFLX_COOLING     to suppress SST cooling below freezing point      **
**                                                                           **
** OPTIONS for MPDATA 3D Advection: Hadvection(itrc,ng)%MPDATA and           **
**                                  Vadvection(itrc,ng)%MPDATA switches      **
**                                                                           **
** TS_MPDATA_LIMIT         to limit upwind corrector fluxes for stability    **
**                                                                           **
** Pressure gradient algorithm OPTIONS:                                      **
**                                                                           **
**   If no option is selected, the pressure gradient term is computed using  **
**   standard density Jacobian algorithm. Notice that there are two quartic  **
**   pressure Jacobian options. They differ on how the WENO reconciliation   **
**   step is done and in the monotonicity constraining algorithms.           **
**                                                                           **
** DJ_GRADPS               if splines density Jacobian (Shchepetkin, 2000)   **
** PJ_GRADP                if finite volume Pressure Jacobian (Lin,1997)     **
** PJ_GRADPQ2              if quartic 2 Pressure Jacobian (Shchepetkin,2000) **
** PJ_GRADPQ4              if quartic 4 Pressure Jacobian (Shchepetkin,2000) **
** WJ_GRADP                if weighted density Jacobian (Song,1998)          **
**                                                                           **
** ATM_PRESS               to impose atmospheric pressure onto sea surface   **
** PRESS_COMPENSATE        to compensate for boundary without ATM pressure   **
**                                                                           **
** OPTIONS for surface fluxes formulation using atmospheric boundary layer   **
** (Fairall et al, 1996):                                                    **
**                                                                           **
**   There are three ways to provide longwave radiation in the atmospheric   **
**   boundary layer: (1) Compute the net longwave radiation internally using **
**   the Berliand (1952) equation (LONGWAVE) as function of air temperature, **
**   sea surface temperature, relative humidity, and cloud fraction;         **
**   (2) provide (read) longwave downwelling radiation only  and then add    **
**   outgoing longwave radiation (LONGWAVE_OUT) as a function of the model   **
**   sea surface temperature; (3) provide net longwave radiation (default).  **
**                                                                           **
** BULK_FLUXES             if bulk fluxes computation                        **
** COOL_SKIN               if cool skin correction                           **
** LONGWAVE                if computing net longwave radiation               **
** LONGWAVE_OUT            if computing outgoing longwave radiation          **
** EMINUSP                 if computing E-P                                  **
** EMINUSP_SSH             if computing changes in SSH due to E-P            **
** RUNOFF                  if adding runoff as a second rain field           **
** RUNOFF_SSH              if adjusting zeta based on runoff field           **
** WIND_MINUS_CURRENT      if compute effective wind by removing current     **
**                                                                           **
** OPTIONS for wave roughness formulation in bulk fluxes:                    **
**                                                                           **
** COARE_TAYLOR_YELLAND      Taylor and Yelland (2001) relation              **
** COARE_OOST                Oost et al (2002) relation                      **
** DRENNAN                   Drennan (2003) relation                         **
** DEEPWATER_WAVES           Deep water waves approximation                  **
**                                                                           **
** OPTIONS for shortwave radiation:                                          **
**                                                                           **
**   The shortwave radiation can be computed using the global albedo         **
**   equation with a cloud correction. Alternatively, input shortwave        **
**   radiation data computed from averaged data (with snapshots greater      **
**   or equal than 24 hours) can be modulated by the local diurnal cycle     **
**   which is a function longitude, latitude and day-of-year.                **
**                                                                           **
** ALBEDO                  if albedo equation for shortwave radiation        **
** DIURNAL_SRFLUX          to impose shortwave radiation local diurnal cycle **
**                                                                           **
** Model configuration OPTIONS:                                              **
**                                                                           **
** SOLVE3D                 if solving 3D primitive equations                 **
** CURVGRID                if curvilinear coordinates grid                   **
** MASKING                 if land/sea masking                               **
** BODYFORCE               if applying stresses as bodyforces                **
** PROFILE                 if time profiling                                 **
** AVERAGES                if writing out NLM time-averaged data             **
** AVERAGES_DETIDE         if writing out NLM time-averaged detided fields   **
** AD_AVERAGES             if writing out ADM time-averaged data             **
** RP_AVERAGES             if writing out TLM time-averaged data             **
** TL_AVERAGES             if writing out ADM time-averaged data             **
** DIAGNOSTICS_BIO         if writing out biological diagnostics             **
** DIAGNOSTICS_UV          if writing out momentum diagnostics               **
** DIAGNOSTICS_TS          if writing out tracer diagnostics                 **
** ICESHELF                if including ice shelf cavities                   **
** SINGLE_PRECISION        if single precision arithmetic numerical kernel   **
** SPHERICAL               if analytical spherical grid                      **
** STATIONS                if writing out station data                       **
** STATIONS_CGRID          if extracting data at native C-grid               **
**                                                                           **
** OPTIONS for Lagrangian drifters:                                          **
**                                                                           **
** FLOATS                  to activate simulated Lagrangian drifters         **
** FLOAT_OYSTER            to activate oyster model behavior in floats       **
** FLOAT_STICKY            to reflect/stick floats that hit surface/bottom   **
** FLOAT_VWALK             if vertical random walk                           **
** VWALK_FORWARD           if forward time stepping vertical random walk     **
** DIAPA                   to simulate diapa                                 **
**                                                                           **
** OPTIONS for analytical fields configuration:                              **
**                                                                           **
**    Any of the analytical expressions are coded in "analytical.F".         **
**                                                                           **
** ANA_BIOLOGY             if analytical biology initial conditions          **
** ANA_BPFLUX              if analytical bottom passive tracers fluxes       **
** ANA_BSFLUX              if analytical bottom salinity flux                **
** ANA_BTFLUX              if analytical bottom temperature flux             **
** ANA_CLOUD               if analytical cloud fraction                      **
** ANA_DIAG                if customized diagnostics                         **
** ANA_DQDSST              if analytical surface heat flux sensitivity to SST**
** ANA_DRAG                if analytical spatially varying drag parameters   **
** ANA_FSOBC               if analytical free-surface boundary conditions    **
** ANA_GRID                if analytical model grid set-up                   **
** ANA_HUMIDITY            if analytical surface air humidity                **
** ANA_INITIAL             if analytical initial conditions                  **
** ANA_M2CLIMA             if analytical 2D momentum climatology             **
** ANA_M2OBC               if analytical 2D momentum boundary conditions     **
** ANA_M3CLIMA             if analytical 3D momentum climatology             **
** ANA_M3OBC               if analytical 3D momentum boundary conditions     **
** ANA_MASK                if analytical Land/Sea masking                    **
** ANA_NUDGCOEF            if analytical climatology nudging coefficients    **
** ANA_PAIR                if analytical surface air pressure                **
** ANA_PASSIVE             if analytical inert tracers initial conditions    **
** ANA_PERTURB             if analytical perturbation of initial conditions  **
** ANA_PSOURCE             if analytical point Sources/Sinks                 **
** ANA_PTOBC               if analytical passive tracers boundary conditions **
** ANA_RAIN                if analytical rain fall rate                      **
** ANA_SEDIMENT            if analytical sediment initial fields             **
** ANA_SMFLUX              if analytical surface momentum stress             **
** ANA_SNOW                for analytic snowfall rate                        **
** ANA_SPFLUX              if analytical surface passive tracers fluxes      **
** ANA_SPINNING            if analytical time-varying rotation force         **
** ANA_SPONGE              if analytical enhanced viscosity/diffusion sponge **
** ANA_SRFLUX              if analytical surface shortwave radiation flux    **
** ANA_SSFLUX              if analytical surface salinity flux               **
** ANA_SSH                 if analytical sea surface height                  **
** ANA_SSS                 if analytical sea surface salinity                **
** ANA_SST                 if analytical sea surface temperature, SST        **
** ANA_STFLUX              if analytical surface net heat flux               **
** ANA_TAIR                if analytical surface air temperature             **
** ANA_TCLIMA              if analytical tracers climatology                 **
** ANA_TOBC                if analytical tracers boundary conditions         **
** ANA_TRC_PSOURCE         if analytical point sources of inert tracers      **
** ANA_VMIX                if analytical vertical mixing coefficients        **
** ANA_WINDS               if analytical surface winds                       **
** ANA_WWAVE               if analytical wind induced waves                  **
**                                                                           **
** OPTIONS for horizontal mixing of momentum:                                **
**                                                                           **
** VISC_GRID               to scale viscosity coefficient by grid size       **
** MIX_S_UV                if mixing along constant S-surfaces               **
** MIX_GEO_UV              if mixing on geopotential (constant Z) surfaces   **
**                                                                           **
** OPTIONS for horizontal mixing of tracers:                                 **
**                                                                           **
** DIFF_GRID               to scale diffusion coefficients by grid size      **
** MIX_S_TS                if mixing along constant S-surfaces               **
** MIX_GEO_TS              if mixing on geopotential (constant Z) surfaces   **
** MIX_ISO_TS              if mixing on epineutral (constant RHO) surfaces   **
** TS_MIX_CLIMA            if diffusion of tracer perturbation (t-tclm)      **
** TS_MIX_MAX_SLOPE        if maximum slope in epineutral diffusion          **
** TS_MIX_MIN_STRAT        if minimum stratification in epineutral diffusion **
** TS_MIX_STABILITY        if weighting diffusion between two time levels    **
**                                                                           **
** OPTIONS for vertical turbulent mixing scheme of momentum and tracers      **
** (activate only one closure):                                              **
**                                                                           **
** BVF_MIXING              if Brunt-Vaisala frequency mixing                 **
** GLS_MIXING              if Generic Length-Scale mixing closure            **
** MY25_MIXING             if Mellor/Yamada Level-2.5 closure                **
** LMD_MIXING              if Large et al. (1994) interior closure           **
**                                                                           **
** LIMIT_VDIFF             to impose an upper limit on vertical diffusion    **
** LIMIT_VVISC             to impose an upper limit on vertical viscosity    **
**                                                                           **
** OPTIONS for the Generic Length-Scale closure (Warner et al., 2005):       **
**                                                                           **
**   The default horizontal advection is third-order upstream bias.  The     **
**   default vertical advection is 4th-order centered advection.             **
**                                                                           **
** CANUTO_A                if Canuto A-stability function formulation        **
** CANUTO_B                if Canuto B-stability function formulation        **
** CHARNOK                 if Charnok surface roughness from wind stress     **
** CRAIG_BANNER            if Craig and Banner wave breaking surface flux    **
** KANTHA_CLAYSON          if Kantha and Clayson stability function          **
** K_C2ADVECTION           if 2nd-order centered advection                   **
** K_C4ADVECTION           if 4th-order centered advection                   **
** N2S2_HORAVG             if horizontal smoothing of buoyancy/shear         **
** RI_SPLINES              if splines reconstruction for vertical sheer      **
** ZOS_HSIG                if surface roughness from wave amplitude          **
** TKE_WAVEDISS            if wave breaking surface flux from wave amplitude **
**                                                                           **
** OPTIONS for the Mellor/Yamada level 2.5 closure:                          **
**                                                                           **
**   The default horizontal advection is third-order upstream bias.  The     **
**   default vertical advection is 4th-order centered advection.             **
**                                                                           **
** N2S2_HORAVG             if horizontal smoothing of buoyancy/shear         **
** KANTHA_CLAYSON          if Kantha and Clayson stability function          **
** K_C2ADVECTION           if 2nd-order centered advection                   **
** K_C4ADVECTION           if 4th-order centered advection                   **
** RI_SPLINES              if splines reconstruction for vertical sheer      **
**                                                                           **
** OPTIONS for the Large et al. (1994) K-profile parameterization mixing:    **
** mixing:                                                                   **
**                                                                           **
** LMD_BKPP                if bottom boundary layer KPP mixing               **
** LMD_CONVEC              to add convective mixing due to shear instability **
** LMD_DDMIX               to add double-diffusive mixing                    **
** LMD_NONLOCAL            if nonlocal transport                             **
** LMD_RIMIX               to add diffusivity due to shear instability       **
** LMD_SHAPIRO             if Shapiro filtering boundary layer depth         **
** LMD_SKPP                if surface boundary layer KPP mixing              **
** M2TIDE_DIFF             to add simulated tidal diffusion                  **
** RI_SPLINES              if splines reconstruction for Richardson Number   **
**                                                                           **
** OPTIONS in the K-profile parameterization to activate smoothing of        **
** Richardson number, if RI_SPLINES is not activated:                        **
**                                                                           **
** RI_HORAVG               if horizontal Richardson number smoothing         **
** RI_VERAVG               if vertical   Richardson number smoothing         **
**                                                                           **
** OPTIONS for Meinte Blass bottom boundary layer closure:                   **
**                                                                           **
**   The Options MB_Z0BL and MB_Z0RIP should be activated concurrently.      **
**                                                                           **
** MB_BBL                  if Meinte Blaas BBL closure                       **
** MB_CALC_ZNOT            if computing bottom roughness internally          **
** MB_CALC_UB              if computing bottom orbital velocity internally   **
** MB_Z0BIO                if biogenic bedform roughness for ripples         **
** MB_Z0BL                 if bedload roughness for ripples                  **
** MB_Z0RIP                if bedform roughness for ripples                  **
**                                                                           **
** OPTIONS for Styles and Glenn (2000) bottom boundary layer closure:        **
**                                                                           **
** SG_BBL                  if Styles and Glenn (2000) BBL closure            **
** SG_CALC_ZNOT            if computing bottom roughness internally          **
** SG_CALC_UB              if computing bottom orbital velocity internally   **
** SG_LOGINT               if logarithmic interpolation of (Ur,Vr)           **
**                                                                           **
** OPTIONS for the Sherwood/Signell/Warner bottom boundary layer closure:    **
**                                                                           **
** SSW_BBL                 if Sherwood et al. BBL closure                    **
** SSW_CALC_ZNOT           if computing bottom roughness internally          **
** SSW_LOGINT              if logarithmic interpolation of (Ur,Vr)           **
** SSW_CALC_UB             if computing bottom orbital velocity internally   **
** SSW_FORM_DRAG_COR       to activate form drag coefficient                 **
** SSW_ZOBIO               if biogenic bedform roughness from ripples        **
** SSW_ZOBL                if bedload roughness for ripples                  **
** SSW_ZORIP               if bedform roughness from ripples                 **
**                                                                           **
** Lateral boundary conditions OPTIONS:                                      **
**                                                                           **
** IMPLICIT_NUDGING        if implicit nudging term in momentum radiation    **
** RADIATION_2D            if tangential phase speed in radiation conditions **
**                                                                           **
** OPTIONS for tidal forcing at open boundaries:                             **
**                                                                           **
**   The tidal data is processed in terms of tidal components, classified by **
**   period. The tidal forcing is computed for the full horizontal grid. If  **
**   requested, the tidal forcing is added to the processed open boundary    **
**   data.                                                                   **
**                                                                           **
**   Both tidal elevation and tidal currents are required to force the model **
**   properly. However, if only the tidal elevation is available, the tidal  **
**   currents at the open boundary can be estimated by reduced physics. Only **
**   the pressure gradient, Coriolis, and surface and bottom stresses terms  **
**   are considered at the open boundary. See "u2dbc_im.F" or "v2dbc_im.F"   **
**   for details. Notice that there is an additional option (FSOBC_REDUCED)  **
**   for the computation of the pressure gradient term in both Flather or    **
**   reduced physics conditions (*_M2FLATHER, *_M2REDUCED).                  **
**                                                                           **
** SSH_TIDES               if imposing tidal elevation                       **
** UV_TIDES                if imposing tidal currents                        **
** POT_TIDES               if imposing potential tides                       **
** RAMP_TIDES              if ramping (over one day) tidal forcing           **
** FSOBC_REDUCED           if SSH data and reduced physics conditions        **
** ADD_FSOBC               to add tidal elevation to processed OBC data      **
** ADD_M2OBC               to add tidal currents  to processed OBC data      **
**                                                                           **
** OPTIONS for reading and processing of climatological fields:              **
**                                                                           **
** OCLIMATOLOGY            if processing 3D vertical momentum climatology    **
** AKTCLIMATOLOGY          if processing 3D vertical salinity diffustion     **
**                                                                           **
** ROMS/TOMS driver OPTIONS:                                                 **
**                                                                           **
** ADM_DRIVER                 if generic adjoint model                       **
** AD_SENSITIVITY             if adjoint sensitivity                         **
** AFT_EIGENMODES             if adjoint finite time eingenmodes             **
** ARRAY_MODES                if 4D-Var representer matrix array modes       **
** CLIPPING                   if R4D-Var representer matrix clipping analysis**
** CORRELATION                if background-error correlation model          **
** ENSEMBLE                   if ensemble prediction                         **
** EVOLVED_LCZ                if 4D-Var evolved Hessian singular vectors     **
** FORCING_SV                 if forcing singular vectors driver             **
** FT_EIGENMODES              if finite time eingenmodes: normal modes       **
** HESSIAN_FSV                if Hessian forcing singular vectors            **
** HESSIAN_SO                 if Hessian stochastic optimals                 **
** HESSIAN_SV                 if Hessian singular vectors                    **
** INNER_PRODUCT              if tangent/adjoint inner product check         **
** I4DVAR                     if incremental 4D-Var data assimilation        **
** I4DVAR_ANA_SENSITIVITY     if I4D-Var observations sensitivity            **
** JEDI                       if using Joint Effort for DA Integration       **
** LCZ_FINAL                  it computing 4D-Var Hessian singular vectors   **
** OPT_OBSERVATIONS           if optimal observations                        **
** OPT_PERTURBATION           if optimal perturbations, singular vectors     **
** PICARD_TEST                if representer tangent linear model test       **
** PSEUDOSPECTRA              if pseudospectra of tangent linear resolvant   **
** RBL4DVAR                   if weak constraint RBL4D-Var data assimilation **
** RBL4DVAR_ANA_SENSITIVITY   if RBL4D-Var analysis observation sensitivity  **
** RBL4DVAR_FCT_SENSITIVITY   if RBL4D-Var forecast observation sensitivity  **
** RPM_DRIVER                 if generic representers model                  **
** R_SYMMETRY                 if representer matrix symmetry test            **
** R4DVAR                     if R4D-Var data assimilation                   **
** R4DVAR_ANA_SENSITIVITY     if R4D-Var analysis observation sensitivity    **
** SANITY_CHECK               if tangent/adjoint codes sanity check          **
** SO_SEMI                    if stochastic optimals driver, semi-norm       **
** SO_TRACE                   if stochastic optimals, randomized trace       **
** SPLIT_I4DVAR               if split I4D-Var data assimilation             **
** SPLIT_RBL4DVAR             if split RBL4D-Var data assimilation           **
** SPLIT_R4DVAR               if split R4D-Var data assimilation             **
** SPLIT_SP4DVAR              if split SP4D-Var data assimilation            **
** SP4DVAR                    if Saddle-Point 4D-Var data assimilation       **
** STOCHASTIC_OPT             if stochastic optimals                         **
** TLM_CHECK                  if tangent linear model linearization check    **
** TLM_DRIVER                 if generic tangent linear model driver         **
**                                                                           **
** OPTIONS associated with tangent linear, representer and adjoint models:   **
**                                                                           **
** AD_IMPULSE              to force adjoint model with intermittent impulses **
** ADJUST_BOUNDARY         if including boundary conditions in 4DVar state   **
** ADJUST_STFLUX           if including surface tracer flux in 4DVar state   **
** ADJUST_WSTRESS          if including wind-stress in 4DVar state           **
** ARRAY_MODES_SPLIT       to separate analysis due to IC, forcing, and OBC  **
** BALANCE_OPERATOR        if error covariance multivariate balance term     **
** BEOFS_ONLY              if computing EOFs of background error covariance  **
** BGQC                    if background quality control of observations     **
** BNORM                   if Background norm Hessian singular vectors       **
** CELERITY_WRITE          if writing radiation celerity in forward file     **
** CLIPPING_SPLIT          to separate analysis due to IC, forcing, and OBC  **
** DATALESS_LOOPS          if testing convergence of Picard iterations       **
** ENKF_RESTART            if writting restart fields for EnKF               **
** FORWARD_FLUXES          if using NLM trajectory surface fluxes            **
** FORWARD_MIXING          if processing forward vertical mixing coefficient **
** FORWARD_WRITE           if writing out forward solution, basic state      **
** FORWARD_READ            if reading in  forward solution, basic state      **
** FORWARD_RHS             if processing forward right-hand-side terms       **
** GEOPOTENTIAL_HCONV      if horizontal convolutions along geopotentials    **
** IMPACT_INNER            to write observations impacts for each inner loop **
** IMPLICIT_VCONV          if implicit vertical convolution algorithm        **
** IMPULSE                 if processing adjoint impulse forcing             **
** MINRES                  if Minimal Residual Method for 4DVar minimization **
** MULTIPLE_TLM            if multiple TLM history files in 4DVAR            **
** NLM_OUTER               if nonlinear model as basic state in outer loop   **
** OBS_IMPACT              if observation impact to 4DVAR data assimilation  **
** OBS_IMPACT_SPLIT        to separate impact due to IC, forcing, and OBC    **
** POSTERIOR_EOFS          if posterior analysis error covariance EOFS       **
** POSTERIOR_ERROR_F       if final posterior analysis error covariance      **
** POSTERIOR_ERROR_I       if initial posterior analysis error covariance    **
** PRIOR_BULK_FLUXES       if imposing prior NLM surface fluxes              **
** RECOMPUTE_4DVAR         if recomputing 4DVar in analysis algorithms       **
** RPCG                    if Restricted B-preconditioned Lanczos solver     **
** RPM_RELAXATION          if Picard iterations, Diffusive Relaxation of RPM **
** SKIP_NLM                to skip running NLM, reading NLM trajectory       **
** SO_SEMI_WHITE           to activate SO semi norm white/red noise processes**
** STOCH_OPT_WHITE         to activate SO white/red noise processes          **
** SPLINES_VCONV           to activate implicit splines vertical convolution **
** TIME_CONV               if weak-constraint 4D-Var time convolutions       **
** VCONVOLUTION            to add vertical correlation to 3D convolution     **
** VERIFICATION            if writing out solution at observation locations  **
** WEAK_NOINTERP           if not time interpolation in weak 4D-Var forcing  **
** ZETA_ELLIPTIC           if SSH elliptic Equation in balance operator      **
**                                                                           **
** OPTION for processing the full grid range (interior and boundary points)  **
** of the state vector in variational data assimilation and generalized      **
** stability theory analysis. Otherwise, only interior points are processed. **
**                                                                           **
** FULL_GRID               to consider both interior and boundary points     **
**                                                                           **
** Fennel et al. (2006) biology model OPTIONS:                               **
**                                                                           **
** BIO_FENNEL              if Fennel et al. (2006) nitrogen-based model      **
** BIO_SEDIMENT            to restore fallen material to the nutrient pool   **
** CARBON                  to add carbon constituents                        **
** DENITRIFICATION         to add denitrification processes                  **
** OCMIP_OXYGEN_SC         if O2 Schmidt number from Keeling et al. (1998)   **
** OXYGEN                  to add oxygen dynamics                            **
** PCO2AIR_DATA            if pCO2 climatology from Laurent et al. (2017)    **
** PCO2AIR_SECULAR         if pCO2 time-depedent evolution                   **
** RW14_C02_SC             if CO2 Schmidt number from Wanninkhof (2014)      **
** RW14_OXYGEN_SC          if O2  Schmidt number from Wanninkhof (2014)      **
** PO4                     if phytoplanckton growth limitef by Phosphorus    **
** RIVER_DON               if DON non-sinking source from rivers             **
** TALK_NONCONSERV         if nonconservative computation of alkalinity      **
**                                                                           **
** Hypoxia ecosysten model OPTIONS:                                          **
**                                                                           **
** HYPOXIA_SRM             if Hypoxia Simple Respiration Model               **
**                                                                           **
** Spectral Fennel biology model OPTIONS:                                    **
**                                                                           **
** SPECTRAL_LIGHT          for spectral light as part of Fennel. Gallegos    **
** CDOM_DEFAULT            for constant default CDOM for spectral light      **
** CDOM_VARIABLE           for variable CDOM for spectral light              **
** MOD_SWR_SPECT       modulate shortwave radiation with Gallegos spectrum   **
** MOD_SWR_HOMO        modulate shortwave radiation evenly through spectrum  **
** CHL_BACKSCAT        add chlorophyll backscatter in attenuation            **
**                                                                           **
** Vegetation growth and flow model OPTIONS:                                 **
**                                                                           **
** VEGETATION              to activate submerged/emergent vegetation effects **
** SAV                 activate submerged vegetation (seagrass) effects      **
** EMERGENT_VEG        activate emergent vegetation (marsh) effects          **
** SEAGRASS_SINK           bottom sink of N due to benthic seagrass proxy    **
** SEAGRASS_LIGHT      add N sink (seagrass proxy) as function of light      **
** SEAGRASS_LIGHT_CONST constant N sink (seagrass proxy) if light exceeded   **
**                                                                           **
** Bering Sea biology model OPTIONS:                                         **
**                                                                           **
** BEST_NPZ                if Gibson et al. Bering Sea model                 **
** STATIONARY              if ??                                             **
** BENTHIC                 if benthic components                             **
** ICE_BIO                 if ice algae                                      **
** JELLY                   if jellyfish                                      **
** CLIM_ICE_1D             if one-D with ice                                 **
**                                                                           **
** NPZD biology model OPTIONS:                                               **
**                                                                           **
** BIO_UMAINE              if Chai et al. (2002) CoSINE model                **
** NPZD_FRANKS             if NPZD Biology model, Franks et al. (1986)       **
** NPZD_IRON               if NPZD Biology model with iron limitation        **
** NPZD_POWELL             if NPZD Biology model, Powell et al. (2006)       **
** IRON_LIMIT              if Fe limitation on phytoplankton growth          **
** IRON_RELAX              if nudging Fe over the shelf, h <= FeHmin         **
**                                                                           **
** Bio-optical EcoSim model OPTIONS:                                         **
**                                                                           **
** ECOSIM                  if bio-optical EcoSim model                       **
** BIO_OPTICAL             to compute underwater spectral light properties   **
**                                                                           **
** Nemuro lower trophic level ecosystem model OPTIONS:                       **
**                                                                           **
**    Need to choose a zooplankton grazing option (HOLLING_GRAZING or        **
**    IVLEV_EXPLICIT). The default implicit IVLEV algorithm does not         **
**    work yet.                                                              **
**                                                                           **
** NEMURO                  if Nemuro ecosystem model.                        **
** NEMURO_SED1             if Nemuro sediment remineralization               **
** PRIMARY_PROD            if primary productivity output                    **
** BIO_SEDIMENT            to restore fallen material to the nutrient pool   **
** HOLLING_GRAZING         Holling-type s-shaped curve grazing (implicit)    **
** IVLEV_EXPLICIT          Ivlev explicit grazing algorithm                  **
**                                                                           **
** Red tide biological model OPTIONS:                                        **
**                                                                           **
** RED_TIDE                if red tide biological model.                     **
**                                                                           **
** Sediment transport model OPTIONS:                                         **
**                                                                           **
** SEDIMENT                to activate sediment transport model              **
** BEDLOAD_MPM             to activate Meyer-Peter-Mueller bed load          **
** BEDLOAD_SOULSBY         to activate Soulsby wave/current bed load         **
** SED_DENS                to activate sediment to affect equation of state  **
** SED_MORPH               to allow bottom model elevation to evolve         **
** SUSPLOAD                to activate suspended load transport              **
** SED_BIODIFF             to activate sediment biodiffusivity               **
** MIXED_BED               to activate mixed bed behavior                    **
** COHESIVE_BED            to activate cohesive bed model                    **
** NONCOHESIVE_BED1        original bed model of Warner et al 2008, default  **
** NONCOHESIVE_BED2        modified bed model of Sherwood et al, in press    **
** SED_FLOCS               flocculation model of Verney et al., 2011         **
** FLOC_TURB_DISS          dissipation for flocculation based on turbulence  **
** FLOC_BBL_DISS           dissipation for flocs from bottom boundary layer  **
** SED_DEFLOC              flocculation decomposition in sediment bed        **
** SED_TAU_CD_CONST        constant critical stress for deposition           **
** SED_TAU_CD_LIN          linear critical stress for deposition             **
**                                                                           **
** Wave effoct on currents (WEC) and shallow water OPTIONS:                  **
**                                                                           **
** WET_DRY             activate wetting and drying                           **
** WEC_MELLOR          activate radiation stress terms from Mellor 08.       **
** WEC_VF              activate wave-current stresses from Uchiyama et al.   **
** WDISS_THORGUZA      activate wave dissipation from Thorton/Guza.          **
** WDISS_CHURTHOR      activate wave dissipation from Church/Thorton.        **
** WDISS_WAVEMOD       activate wave dissipation from a wave model           **
** WDISS_GAMMA         activate wave dissipation when using InWave           **
** WDISS_ROELVINK      activate wave dissipation Roelvink when using InWave  **
** ROLLER_SVENDSEN     activate wave roller based on Svendsen                **
** ROLLER_MONO         activate wave roller for monchromatic waves           **
** ROLLER_RENIERS      activate wave roller based on Reniers                 **
** BOTTOM_STREAMING    activate wave enhanced bottom streaming               **
** SURFACE_STREAMING   activate wave enhanced surface streaming              **
** WAVE_MIXING         activate enhanced vertical viscosity mixing from waves**
**                                                                           **
** OPTIONS for grid nesting:                                                 **
**                                                                           **
** NESTING                 to activate grid nesting: composite/refinement    **
** NESTING_DEBUG           to check mass fluxes conservation in refinement   **
** NO_CORRECT_TRACER       to avoid two-way correction of boundary tracer    **
** ONE_WAY                 if one-way nesting in refinement grids            **
** TIME_INTERP_FLUX        time interpolate coarse mass flux instead persist **
**                                                                           **
** OPTIONS for coupling to other Earth System Models (ESM) via the Earth     **
** Modeling Framework (ESMF) or Modeling Coupling Toolkit (MCT) libraries.   **
** If coupling with ESMF library, it uses the National Unified Operational   **
** Prediction Capability (NUOPC) layer "cap" files to facilitate exchanges   **
** with other ESM components.                                                **
**                                                                           **
** ESMF_LIB                if coupling with the ESMF/NUOPC library           **
** MCT_LIB                 if Coupling with the MCT library                  **
**                                                                           **
** CICE_COUPLING           if coupling to CICE sea ice model                 **
** COAMPS_COUPLING         if coupling to COAMPS atmospheric model           **
** DATA_COUPLING           if coupling to DATA model                         **
** EXCLUDE_SPONGE          if excluding sponge point in export fields        **
** FRC_COUPLING            if forcing from Atmopheric or Data model          **
** REFDIF_COUPLING         if coupling to REFDIT wave model                  **
** REGCM_COUPLING          if coupling to RegCM atmospheric model            **
** SWAN_COUPLING           if coupling to SWAN wave model                    **
** TIME_INTERP             if importing snapshots for time interpolation     **
** WAM_COUPLING            if coupling to WAM wave model                     **
** WRF_COUPLING            if coupling to WRF atmospheric model              **
** WRF_TIMEAVG             if time-averaged fields over coupling interval    **
**                                                                           **
** Nearshore and shallow water model OPTIONS:                                **
**                                                                           **
** WET_DRY                 to activate wetting and drying                    **
** NEARSHORE_MELLOR05      to activate radiation stress terms (Mellor 2005). **
** NEARSHORE_MELLOR08      to activate radiation stress terms (Mellor 2008). **
**                                                                           **
** MPI communication OPTIONS:  The routines "mp_assemble" (used in nesting), **
**                             "mp_collect" (used in NetCDF I/O and 4D-Var), **
** and "mp_reduce" (used in global reductions) are coded in "distribution.F" **
** by either using low-level ("mpi_isend" and "mpi_irecv") or high-level     **
** ("mpi_allgather" and "mpi_allreduce") MPI calls. The default is to use    **
** the low-level MPI  calls. The options for routine "mp_boundary" (used to  **
** process lateral open boundary conditions is either "mpi_allgather" or     **
** "mpi_allreduce" (default).                                                **
**                                                                           **
** The user needs to be aware that the choice of these MPI communication     **
** routines it will affect performance issue. In some computers, the         **
** low-level are either slower or faster than the high-level MPI library     **
** calls. It depends on the computer (cluster) set-up. Some vendors provide  **
** native MPI libraries fine-tuned for the computer architecture. The        **
** user needs to find which function option performs better by carrying on   **
** benchmarks. We provides the following choices:                            **
**                                                                           **
** ASSEMBLE_ALLGATHER  use "mpi_allgather" in "mp_assemble"                  **
** ASSEMBLE_ALLREDUCE  use "mpi_allreduce" in "mp_assemble"                  **
**                                                                           **
** BOUNDARY_ALLGATHER  use "mpi_allgather" in "mp_boundary"                  **
**                                                                           **
** COLLECT_ALLGATHER   use "mpi_allgather" in "mp_collect"                   **
** COLLECT_ALLREDUCE   use "mpi_allreduce" in "mp_collect"                   **
**                                                                           **
** REDUCE_ALLGATHER    use "mpi_allgather" in "mp_reduce"                    **
** REDUCE_ALLREDUCE    use "mpi_allreduce" in "mp_reduce"                    **
**                                                                           **
** NetCDF input/output OPTIONS:                                              **
**                                                                           **
** DEFLATE                 to set compression NetCDF-4/HDF5 format files     **
** HDF5                    to create NetCDF-4/HDF5 format files              **
** NO_LBC_ATT              to not check NLM_LBC global attribute on restart  **
** NO_READ_GHOST           to not include ghost points during read/scatter   **
** NO_WRITE_GRID           if not writing grid arrays                        **
** PARALLEL_IO             if parallel I/O via HDF5 or pnetcdf libraries     **
** PERFECT_RESTART         to include perfect restart variables              **
** PNETCDF                 if parallel I/O with pnetcdf (classic format)     **
** POSITIVE_ZERO           to impose positive zero in ouput data             **
** READ_WATER              if only reading water points data                 **
** REGRID_SHAPIRO          to apply Shapiro Filter to regridded data         **
** WRITE_WATER             if only writing water points data                 **
** RST_SINGLE              if writing single precision restart fields        **
** OUT_DOUBLE              if writing double precision output fields         **
**                                                                           **
** OPTION to process 3D data by levels (2D slabs) to reduce memory needs in  **
** distributed-memory configurations. This option is convenient for large    **
** problems on nodes with limited memory.                                    **
**                                                                           **
** INLINE_2DIO             if processing 3D IO level by level                **
**                                                                           **
** OPTION to avoid writing current date and CPP options to NetCDF file       **
** headers. This is used to compare serial and parallel solutions where      **
** the UNIX command "diff" is used between NetCDF files. It will only        **
** tell us that the binary files are different or not. Finding the           **
** parallel bug is complete different story.                                 **
**                                                                           **
** DEBUGGING           use to activate parallel debugging switch             **
**                                                                           **
*******************************************************************************
*******************************************************************************
*******************************************************************************
**                                                                           **
** Idealized Test Problems:                                                  **
**                                                                           **
** BASIN               Big Bad Basin Example                                 **
** BENCHMARK           Benchmark Tests (small, Medium, big grids)            **
** BIO_TOY             One-dimension (vertical) Biology Toy                  **
** BL_TEST             Boundary Layers Test                                  **
** CHANNEL             Periodic channel, Optimal Perturbations Test          **
** CANYON              Coastal form stress Canyon Test                       **
** CHANNEL_NECK        Channel with a Constriction                           **
** COUPLING_TEST       Two-way Atmosphere-Ocean Coupling Test                **
** DOGBONE             Idealize nesting grids (Composite/Refinement) Test    **
** DOUBLE_GYRE         Idealized Double-gyre Example                         **
** ESTUARY_TEST        Test Estuary for Sediment                             **
** FLT_TEST            Float Tracking Example                                **
** GRAV_ADJ            Gravitational Adjustment Example                      **
** INLET_TEST          Test Inlet Application                                **
** KELVIN              Kelvin wave test                                      **
** LAB_CANYON          Lab Canyon, Polar Coordinates Example                 **
** LAKE_JERSEY         Lake Jersey Nesting Test Case                         **
** LAKE_SIGNELL        Lake Signell Sediment Test Case                       **
** LMD_TEST            Test for LMD and KPP                                  **
** OVERFLOW            Gravitational/Overflow Example                        **
** RIVERPLUME1         River Plume Example 1                                 **
** RIVERPLUME2         River plume Example 2 (Hyatt and Signell)             **
** SEAMOUNT            Seamount Example                                      **
** SED_TEST1           Suspended Sediment Test in a Channel                  **
** SED_TOY             One-dimension (vertical) Sediment Toy                 **
** SHOREFACE           Shore Face Planar Beach Test Case                     **
** SOLITON             Equatorial Rossby Wave Example                        **
** TEST_CHAN           Sediment Test Channel Case                            **
** TEST_HEAD           Sediment Test Headland Case                           **
** UPWELLING           Upwelling Example (default)                           **
** WEDDELL             Idealized Weddell Sea Shelf Application               **
** WINDBASIN           Linear Wind-driven Constant Coriolis Basin            **
**                                                                           **
** Climatological Applications: (See www.myroms.org/Datasets)                **
**                                                                           **
** DAMEE_4             North Atlantic DAMEE Application, 3/4 degree          **
**                                                                           **
** Selected Realistic Applications:                                          **
**                                                                           **
** ADRIA02             Adriatic Sea Application                              **
** NJ_BIGHT            New Jersey Bight Application                          **
** WC13                California Current System, 1/3 degree resolution      **
**                                                                           **
*******************************************************************************
*******************************************************************************
*******************************************************************************
**                                                                           **
**  The user needs to choose either a pre-defined application or his/her     **
**  own application. The application CPP flag to run is activated in the     **
**  makefile. For example, to activate the upwelling example (UPWELLING)     **
**  set:                                                                     **
**                                                                           **
**    ROMS_APPLICATION ?= UPWELLING                                          **
**                                                                           **
**  in the makefile. ROMS will include the associated header file located    **
**  in the ROMS/Include directory. The application header file name is the   **
**  lowercase value of ROMS_APPLICATION with the .h extension and passed     **
**  as ROMS_HEADER definition during  C-preprocessing.  For example, the     **
**  upwelling test problem includes the "upwelling.h" header file:           **
**                                                                           **
**    ROMS_HEADER="upwelling.h"                                              **
**                                                                           **
**  If building a new application, choose an unique CPP flag for it and      **
**  create its associated include file (*.h) to specify the appropriate      **
**  configuration options.                                                   **
**                                                                           **
*******************************************************************************
*/

#if defined ROMS_HEADER
# include ROMS_HEADER
#else
   CPPDEFS - Choose an appropriate ROMS application.
#endif

/*
**  Include internal CPP definitions.
*/

#include "globaldefs.h"

#ifdef INWAVE_MODEL
#include "../../InWave/Include/inwave.h"
#endif
