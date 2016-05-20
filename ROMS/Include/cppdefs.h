/*
** Include file "cppdefs.h"
**
** svn $Id: cppdefs.h 818 2008-11-02 14:56:47Z jcwarner $
********************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                               **
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
** UV_ADV              use to turn ON or OFF advection terms                 **
** UV_COR              use to turn ON or OFF Coriolis term                   **
** UV_U3ADV_SPLIT      use if 3rd-order upstream split momentum advection    **
** UV_C2ADVECTION      use to turn ON or OFF 2nd-order centered advection    **
** UV_C4ADVECTION      use to turn ON or OFF 4th-order centered advection    **
** UV_SADVECTION       use to turn ON or OFF splines vertical advection      **
** UV_VIS2             use to turn ON or OFF harmonic horizontal mixing      **
** UV_VIS4             use to turn ON or OFF biharmonic horizontal mixing    **
** UV_SMAGORINSKY      use to turn ON or OFF Smagorinsky-like viscosity      **
** UV_DRAG_GRID        use if spatially varying bottom friction parameters   **
** UV_LOGDRAG          use to turn ON or OFF logarithmic bottom friction     **
** UV_LDRAG            use to turn ON or OFF linear bottom friction          **
** UV_QDRAG            use to turn ON or OFF quadratic bottom friction       **
** UV_WAVEDRAG         use to turn ON or OFF extra linear bottom wave drag   **
** SPLINES_VVISC       use if splines reconstruction of vertical viscosity   **
**                                                                           **
** OPTIONS associated with tracers equations:                                **
**                                                                           **
**   The default horizontal and vertical advection is 4th-order centered.    **
**                                                                           **
**   The 3rd-order upstream split advection (TS_U3ADV_SPLIT) can be used     **
**   to correct for the spurious diapycnal diffusion of the advection        **
**   operator in terrain-following coordinates. If this is the case, the     **
**   advection operator is split in advective and diffusive components       **
**   and several internal flags are activated in "globaldefs.h".  Notice     **
**   that horizontal and vertical advection of tracer is 4th-order centered  **
**   plus biharmonic diffusion to correct for spurious diapycnal mixing.     **
**   The total time-dependent horizontal mixing coefficient are computed     **
**   in "hmixing.F". It is also recommended to use the rotated mixing        **
**   tensor along geopotentials (MIX_GEO_TS) for the biharmonic operator.    **
**                                                                           **
**   WARNING:  Use the splines vertical advection option (TS_SVADVECTION)    **
**             only in idealized, high vertical resolution applications.     **
**                                                                           **
** TS_U3ADV_SPLIT      use if 3rd-order upstream split tracer advection      **
** TS_A4HADVECTION     use if 4th-order Akima horizontal advection           **
** TS_C2HADVECTION     use if 2nd-order centered horizontal advection        **
** TS_C4HADVECTION     use if 4th-order centered horizontal advection        **
** TS_MPDATA           use if recursive MPDATA 3D advection                  **
** TS_U3HADVECTION     use if 3rd-order upstream horiz. advection            **
** TS_A4VADVECTION     use if 4th-order Akima vertical advection             **
** TS_C2VADVECTION     use if 2nd-order centered vertical advection          **
** TS_C4VADVECTION     use if 4th-order centered vertical advection          **
** TS_SVADVECTION      use if splines vertical advection                     **
** TS_DIF2             use to turn ON or OFF harmonic horizontal mixing      **
** TS_DIF4             use to turn ON or OFF biharmonic horizontal mixing    **
** TS_SMAGORINSKY      use to turn ON or OFF Smagorinsky-like diffusion      **
** TS_FIXED            use if diagnostic run, no evolution of tracers        **
** T_PASSIVE           use if inert passive tracers (dyes, etc)              **
** AGE_MEAN            use if computing Mean Age of inert passive tracers    **
** NONLIN_EOS          use if using nonlinear equation of state              **
** QCORRECTION         use if net heat flux correction                       **
** SALINITY            use if having salinity                                **
** SCORRECTION         use if freshwater flux correction                     **
** SSSC_THRESHOLD      limit on freshwater flux correction                   **
** SOLAR_SOURCE        use if solar radiation source term                    **
** SPLINES_VDIFF       use if splines reconstruction of vertical diffusion   **
** SRELAXATION         use if salinity relaxation as a freshwater flux       **
** TRC_PSOURCE         use if source of inert passive tracers (dyes, etc)    **
** ONE_TRACER_SOURCE   use if one value per tracer for all sources           **
** TWO_D_TRACER_SOURCE use if one value per tracer per source                **
** WTYPE_GRID          use to turn ON spatially varying Jerlov water type    **
**                                                                           **
** OPTION to suppress further surface cooling if the SST is at freezing      **
** point or below and the net surface heat flux is cooling:                  **
**                                                                           **
** LIMIT_STFLX_COOLING use to suppress SST cooling below freezing point      **
**                                                                           **
** OPTIONS for MPDATA 3D Advection:                                          **
**                                                                           **
** TS_MPDATA_LIMIT     use to limit upwind corrector fluxes for stability    **
**                                                                           **
** Tracer advection OPTIONS for adjoint-based algorithms:                    **
**                                                                           **
**   Some of the tracer advection algorithms are highly nonlinear and        **
**   may become unstable when running the tangent linear, representer,       **
**   and adjoint models. This may affect the convergence of the 4DVar        **
**   data assimilation algorithms. Therefore, it is possible to choose       **
**   a simpler (less nonlinear) horizontal and vertical tracer advection     **
**   scheme, if so desired, for the tangent linear, representer and          **
**   adjoint models. Notice that this strategy still allows us to use        **
**   highly nonlinear tracer advection schemes in the basic state upon       **
**   which the tangent linear and adjoint models are linearized. Also,       **
**   it allows us to use those schemes that have not been adjointed yet,     **
**   for example, TS_MPDATA.  Recall that basic state trajectory is          **
**   computed by running the nonlinear model.                                **
**                                                                           **
**   The flags below are optional. By default, the same options chosen       **
**   for the nonlinear model are selected for the tangent linear,            **
**   representer, and adjoint models.                                        **
**                                                                           **
** TS_A4HADVECTION_TL  use if 4th-order Akima horizontal advection           **
** TS_C2HADVECTION_TL  use if 2nd-order centered horizontal advection        **
** TS_C4HADVECTION_TL  use if 4th-order centered horizontal advection        **
** TS_U3HADVECTION_TL  use if 3rd-order upstream horiz. advection            **
**                                                                           **
** TS_A4VADVECTION_TL  use if 4th-order Akima vertical advection             **
** TS_C2VADVECTION_TL  use if 2nd-order centered vertical advection          **
** TS_C4VADVECTION_TL  use if 4th-order centered vertical advection          **
** TS_SVADVECTION_TL   use if splines vertical advection                     **
**                                                                           **
** Pressure gradient algorithm OPTIONS:                                      **
**                                                                           **
**   If no option is selected, the pressure gradient term is computed using  **
**   standard density Jacobian algorithm. Notice that there are two quartic  **
**   pressure Jacobian options. They differ on how the WENO reconciliation   **
**   step is done and in the monotonicity constraining algorithms.           **
**                                                                           **
** DJ_GRADPS           use if splines density Jacobian (Shchepetkin, 2000)   **
** PJ_GRADP            use if finite volume Pressure Jacobian (Lin,1997)     **
** PJ_GRADPQ2          use if quartic 2 Pressure Jacobian (Shchepetkin,2000) **
** PJ_GRADPQ4          use if quartic 4 Pressure Jacobian (Shchepetkin,2000) **
** WJ_GRADP            use if weighted density Jacobian (Song,1998)          **
**                                                                           **
** ATM_PRESS           use to impose atmospheric pressure onto sea surface   **
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
** BULK_FLUXES         use if bulk fluxes computation                        **
** CCSM_FLUXES         use if CCSM version of bulk fluxes computation        **
** NCEP_FLUXES         use if NCEP forcing files are used                    **
** NL_BULK_FLUXES      use bulk fluxes computed by nonlinear model           **
** COOL_SKIN           use if cool skin correction                           **
** LONGWAVE            use if computing net longwave radiation               **
** LONGWAVE_OUT        use if computing outgoing longwave radiation          **
** EMINUSP             use if computing E-P                                  **
** EMINUSP_SSH         use if computing changes in SSH due to E-P            **
** RUNOFF              use if adding runoff as a second rain field           **
** RUNOFF_SSH          use if adjusting zeta based on runoff field           **
**                                                                           **
** OPTIONS for wave roughness formulation in bulk fluxes:                    **
**                                                                           **
** COARE_TAYLOR_YELLAND  use Taylor and Yelland (2001) relation              **
** COARE_OOST            use Oost et al (2002) relation                      **
** DRENNAN               use Drennan (2003) relation                         **
** DEEPWATER_WAVES       use Deep water waves approximation                  **
**                                                                           **
** OPTIONS for shortwave radiation:                                          **
**                                                                           **
**   The shortwave radiation can be computed using the global albedo         **
**   equation with a cloud correction. Alternatively, input shortwave        **
**   radiation data computed from averaged data (with snapshots greater      **
**   or equal than 24 hours) can be modulated by the local diurnal cycle     **
**   which is a function longitude, latitude and day-of-year.                **
**                                                                           **
** ALBEDO_CLOUD        use if albedo equation for shortwave radiation        **
** ALBEDO_CSIM         use if albedo function from CSIM for ice              **
** ALBEDO_CURVE        use if albedo function of lat from Large and Yeager   **
** ALBEDO_FILE         use if albedo read from a file                        **
** DIURNAL_SRFLUX      use to impose shortwave radiation local diurnal cycle **
**                                                                           **
** Model configuration OPTIONS:                                              **
**                                                                           **
** SOLVE3D             use if solving 3D primitive equations                 **
** CURVGRID            use if curvilinear coordinates grid                   **
** MASKING             use if land/sea masking                               **
** BODYFORCE           use if applying stresses as bodyforces                **
** PROFILE             use if time profiling                                 **
** AVERAGES            use if writing out NLM time-averaged data             **
** AVERAGES2           use if writing out secondary time-averaged data       **
** HISTORY2            use if writing out secondary history data             **
** AVERAGES_DETIDE     use if writing out NLM time-averaged detided fields   **
** AD_AVERAGES         use if writing out ADM time-averaged data             **
** RP_AVERAGES         use if writing out TLM time-averaged data             **
** TL_AVERAGES         use if writing out ADM time-averaged data             **
** DIAGNOSTICS_BIO     use if writing out biological diagnostics             **
** DIAGNOSTICS_UV      use if writing out momentum diagnostics               **
** DIAGNOSTICS_TS      use if writing out tracer diagnostics                 **
** ICESHELF            use if including ice shelf cavities                   **
** SPHERICAL           use if analytical spherical grid                      **
** STATIONS            use if writing out station data                       **
** STATIONS_CGRID      use if extracting data at native C-grid               **
**                                                                           **
** OPTIONS for floating iceshelves:                                          **
**                                                                           **
** ICESHELF_3EQ        use if using 3 eqn. thermodynamic formulation         **
**                                                                           **
** OPTIONS for Lagrangian drifters:                                          **
**                                                                           **
** FLOATS              use to activate simulated Lagrangian drifters         **
** FLOAT_OYSTER        use to activate oyster model behavior in floats       **
** FLOAT_STICKY        use to reflect/stick floats that hit surface/bottom   **
** FLOAT_VWALK         use if vertical random walk                           **
** VWALK_FORWARD       use if forward time stepping vertical random walk     **
** DIAPAUSE            use to simulate diapause                              **
**                                                                           **
** OPTIONS for analytical fields configuration:                              **
**                                                                           **
**    Any of the analytical expressions are coded in "analytical.F".         **
**                                                                           **
** ANA_ALBEDO          use if analytical albedo values                       **
** ANA_BIOLOGY         use if analytical biology initial conditions          **
** ANA_BPFLUX          use if analytical bottom passive tracers fluxes       **
** ANA_BSFLUX          use if analytical bottom salinity flux                **
** ANA_BTFLUX          use if analytical bottom temperature flux             **
** ANA_CLOUD           use if analytical cloud fraction                      **
** ANA_DIAG            use if customized diagnostics                         **
** ANA_DQDSST          use if analytical surface heat flux sensitivity to SST**
** ANA_DRAG            use if analytical spatially varying drag parameters   **
** ANA_FSOBC           use if analytical free-surface boundary conditions    **
** ANA_GRID            use if analytical model grid set-up                   **
** ANA_HUMIDITY        use if analytical surface air humidity                **
** ANA_ICE             use for analytic ice initial conditions               **
** ANA_INITIAL         use if analytical initial conditions                  **
** ANA_LRFLUX          use if analytical surface longwave radiation flux     **
** ANA_M2CLIMA         use if analytical 2D momentum climatology             **
** ANA_M2OBC           use if analytical 2D momentum boundary conditions     **
** ANA_M3CLIMA         use if analytical 3D momentum climatology             **
** ANA_M3OBC           use if analytical 3D momentum boundary conditions     **
** ANA_MASK            use if analytical Land/Sea masking                    **
** ANA_NUDGCOEF        use if analytical climatology nudging coefficients    **
** ANA_PAIR            use if analytical surface air pressure                **
** ANA_PASSIVE         use if analytical inert tracers initial conditions    **
** ANA_PERTURB         use if analytical perturbation of initial conditions  **
** ANA_PSOURCE         use if analytical point Sources/Sinks                 **
** ANA_PTOBC           use if analytical passive tracers boundary conditions **
** ANA_RAIN            use if analytical rain fall rate                      **
** ANA_SEDIMENT        use if analytical sediment initial fields             **
** ANA_SMFLUX          use if analytical surface momentum stress             **
** ANA_SNOW            use for analytic snowfall rate                        **
** ANA_SPFLUX          use if analytical surface passive tracers fluxes      **
** ANA_SPINNING        use if analytical time-varying rotation force         **
** ANA_SPONGE          use if analytical enhanced viscosity/diffusion sponge **
** ANA_SRFLUX          use if analytical surface shortwave radiation flux    **
** ANA_SSFLUX          use if analytical surface salinity flux               **
** ANA_SSH             use if analytical sea surface height                  **
** ANA_SSS             use if analytical sea surface salinity                **
** ANA_SST             use if analytical sea surface temperature, SST        **
** ANA_STFLUX          use if analytical surface net heat flux               **
** ANA_TAIR            use if analytical surface air temperature             **
** ANA_TCLIMA          use if analytical tracers climatology                 **
** ANA_TOBC            use if analytical tracers boundary conditions         **
** ANA_TRC_PSOURCE     use if analytical point sources of inert tracers      **
** ANA_VMIX            use if analytical vertical mixing coefficients        **
** ANA_WINDS           use if analytical surface winds                       **
** ANA_WWAVE           use if analytical wind induced waves                  **
**                                                                           **
** OPTIONS for horizontal mixing of momentum:                                **
**                                                                           **
** VISC_GRID           use to scale viscosity coefficient by grid size       **
** MIX_S_UV            use if mixing along constant S-surfaces               **
** MIX_GEO_UV          use if mixing on geopotential (constant Z) surfaces   **
**                                                                           **
** OPTIONS for horizontal mixing of tracers:                                 **
**                                                                           **
** DIFF_GRID           use to scale diffusion coefficients by grid size      **
** MIX_S_TS            use if mixing along constant S-surfaces               **
** MIX_GEO_TS          use if mixing on geopotential (constant Z) surfaces   **
** MIX_ISO_TS          use if mixing on epineutral (constant RHO) surfaces   **
** TS_MIX_CLIMA        use if diffusion of tracer perturbation (t-tclm)      **
** TS_MIX_MAX_SLOPE    use if maximum slope in epineutral diffusion          **
** TS_MIX_MIN_STRAT    use if minimum stratification in epineutral diffusion **
** TS_MIX_STABILITY    use if weighting diffusion between two time levels    **
**                                                                           **
** OPTIONS for vertical turbulent mixing scheme of momentum and tracers      **
** (activate only one closure):                                              **
**                                                                           **
** BVF_MIXING          use if Brunt-Vaisala frequency mixing                 **
** GLS_MIXING          use if Generic Length-Scale mixing                    **
** MY25_MIXING         use if Mellor/Yamada Level-2.5 closure                **
** LMD_MIXING          use if Large et al. (1994) interior closure           **
**                                                                           **
** OPTIONS for the Generic Length-Scale closure (Warner et al., 2005):       **
**                                                                           **
**   The default horizontal advection is third-order upstream bias.  The     **
**   default vertical advection is 4th-order centered advection.             **
**                                                                           **
** CANUTO_A            use if Canuto A-stability function formulation        **
** CANUTO_B            use if Canuto B-stability function formulation        **
** CHARNOK             use if Charnok surface roughness from wind stress     **
** CRAIG_BANNER        use if Craig and Banner wave breaking surface flux    **
** KANTHA_CLAYSON      use if Kantha and Clayson stability function          **
** K_C2ADVECTION       use if 2nd-order centered advection                   **
** K_C4ADVECTION       use if 4th-order centered advection                   **
** N2S2_HORAVG         use if horizontal smoothing of buoyancy/shear         **
** RI_SPLINES          use if splines reconstruction for vertical sheer      **
** ZOS_HSIG            use if surface roughness from wave amplitude          **
** TKE_WAVEDISS        use if wave breaking surface flux from wave amplitude **
**                                                                           **
** OPTIONS for the Mellor/Yamada level 2.5 closure:                          **
**                                                                           **
**   The default horizontal advection is third-order upstream bias.  The     **
**   default vertical advection is 4th-order centered advection.             **
**                                                                           **
** N2S2_HORAVG         use if horizontal smoothing of buoyancy/shear         **
** KANTHA_CLAYSON      use if Kantha and Clayson stability function          **
** K_C2ADVECTION       use if 2nd-order centered advection                   **
** K_C4ADVECTION       use if 4th-order centered advection                   **
** RI_SPLINES          use if splines reconstruction for vertical sheer      **
**                                                                           **
** OPTIONS for the Large et al. (1994) K-profile parameterization mixing:    **
** mixing:                                                                   **
**                                                                           **
** LMD_BKPP            use if bottom boundary layer KPP mixing               **
** LMD_CONVEC          use to add convective mixing due to shear instability **
** LMD_DDMIX           use to add double-diffusive mixing                    **
** LMD_NONLOCAL        use if nonlocal transport                             **
** LMD_RIMIX           use to add diffusivity due to shear instability       **
** LMD_SHAPIRO         use if Shapiro filtering boundary layer depth         **
** LMD_SKPP            use if surface boundary layer KPP mixing              **
** M2TIDE_DIFF         use to add simulated tidal diffusion                  **
** RI_SPLINES          use if splines reconstruction for Richardson Number   **
**                                                                           **
** OPTIONS in the K-profile parameterization to activate smoothing of        **
** Richardson number, if RI_SPLINES is not activated:                        **
**                                                                           **
** RI_HORAVG           use if horizontal Richardson number smoothing         **
** RI_VERAVG           use if vertical   Richardson number smoothing         **
**                                                                           **
** OPTIONS for Meinte Blass bottom boundary layer closure:                   **
**                                                                           **
**   The Options MB_Z0BL and MB_Z0RIP should be activated concurrently.      **
**                                                                           **
** MB_BBL              use if Meinte Blaas BBL closure                       **
** MB_CALC_ZNOT        use if computing bottom roughness internally          **
** MB_CALC_UB          use if computing bottom orbital velocity internally   **
** MB_Z0BIO            use if biogenic bedform roughness for ripples         **
** MB_Z0BL             use if bedload roughness for ripples                  **
** MB_Z0RIP            use if bedform roughness for ripples                  **
**                                                                           **
** OPTIONS for Styles and Glenn (2000) bottom boundary layer closure:        **
**                                                                           **
** SG_BBL              use if Styles and Glenn (2000) BBL closure            **
** SG_CALC_ZNOT        use if computing bottom roughness internally          **
** SG_CALC_UB          use if computing bottom orbital velocity internally   **
** SG_LOGINT           use if logarithmic interpolation of (Ur,Vr)           **
**                                                                           **
** OPTIONS for the Sherwood/Signell/Warner bottom boundary layer closure:    **
**                                                                           **
** SSW_BBL             use if Sherwood et al. BBL closure                    **
** SSW_CALC_ZNOT       use if computing bottom roughness internally          **
** SSW_LOGINT          use if logarithmic interpolation of (Ur,Vr)           **
** SSW_CALC_UB         use if computing bottom orbital velocity internally   **
** SSW_FORM_DRAG_COR   use to activate form drag coefficient                 **
** SSW_ZOBIO           use if biogenic bedform roughness from ripples        **
** SSW_ZOBL            use if bedload roughness for ripples                  **
** SSW_ZORIP           use if bedform roughness from ripples                 **
**                                                                           **
** Lateral boundary conditions OPTIONS:                                      **
**                                                                           **
** RADIATION_2D        use if tangential phase speed in radiation conditions **
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
** SSH_TIDES           use if imposing tidal elevation                       **
** UV_TIDES            use if imposing tidal currents                        **
** POT_TIDES           use if imposing potential tides                       **
** RAMP_TIDES          use if ramping (over one day) tidal forcing           **
** FSOBC_REDUCED       use if SSH data and reduced physics conditions        **
** ADD_FSOBC           use to add tidal elevation to processed OBC data      **
** ADD_M2OBC           use to add tidal currents  to processed OBC data      **
**                                                                           **
** OPTIONS for reading and processing of climatological fields:              **
**                                                                           **
** OCLIMATOLOGY        use if processing 3D vertical momentum climatology    **
** AKTCLIMATOLOGY      use if processing 3D vertical salinity diffustion     **
**                                                                           **
** ROMS/TOMS driver OPTIONS:                                                 **
**                                                                           **
** ADM_DRIVER          use if generic adjoint model driver                   **
** AD_SENSITIVITY      use if adjoint sensitivity driver                     **
** AFT_EIGENMODES      use if adjoint finite time eingenmodes driver         **
** ARRAY_MODES         use if W4DVAR representer matrix array modes          **
** CLIPPING            use if W4DVAR representer matrix clipping analysis    **
** CORRELATION         use if background-error correlation model driver      **
** ENSEMBLE            use if ensemble prediction driver                     **
** FORCING_SV          use if forcing singular vectors driver                **
** FT_EIGENMODES       use if finite time eingenmodes driver: normal modes   **
** INNER_PRODUCT       use if tangent linear and adjoint inner product check **
** IS4DVAR             use if incremental 4DVar data assimilation            **
** IS4DVAR_SENSITIVITY use if I4DVar observations sensitivity driver         **
** HESSIAN_FSV         use if Hessian forcing singular vectors               **
** HESSIAN_SO          use if Hessian stochastic optimals                    **
** HESSIAN_SV          use if Hessian singular vectors                       **
** OPT_OBSERVATIONS    use if optimal observations driver                    **
** OPT_PERTURBATION    use if optimal perturbations driver, singular vectors **
** PICARD_TEST         use if representer tangent linear model test          **
** PSEUDOSPECTRA       use if pseudospectra of tangent linear resolvant      **
** R_SYMMETRY          use if representer matrix symmetry test               **
** RPM_DRIVER          use if generic representers model driver              **
** SANITY_CHECK        use if tangent linear and adjoint codes sanity check  **
** SO_SEMI             use if stochastic optimals driver, semi-norm          **
** SO_TRACE            use if stochastic optimals, randomized trace          **
** STOCHASTIC_OPT      use if stochastic optimals                            **
** TLM_CHECK           use if tangent linear model linearization check       **
** TLM_DRIVER          use if generic tangent linear model driver            **
** W4DPSAS             use if weak constraint 4DPSAS data assimilation       **
** W4DPSAS_SENSITIVITY use if weak constraint 4DPSAS observation sensitivity **
** W4DVAR              use if Weak constraint 4DVar data assimilation        **
** W4DVAR_SENSITIVITY  use if Weak constraint 4DVar observation sensitivity  **
**                                                                           **
** OPTIONS associated with tangent linear, representer and adjoint models:   **
**                                                                           **
** AD_IMPULSE          use to force adjoint model with intermittent impulses **
** ADJUST_BOUNDARY     use if including boundary conditions in 4DVar state   **
** ADJUST_STFLUX       use if including surface tracer flux in 4DVar state   **
** ADJUST_WSTRESS      use if including wind-stress in 4DVar state           **
** ARRAY_MODES_SPLIT   use to separate analysis due to IC, forcing, and OBC  **
** BALANCE_OPERATOR    use if error covariance multivariate balance term     **
** CELERITY_WRITE      use if writing radiation celerity in forward file     **
** CLIPPING_SPLIT      use to separate analysis due to IC, forcing, and OBC  **
** DATALESS_LOOPS      use if testing convergence of Picard iterations       **
** FORWARD_MIXING      use if processing forward vertical mixing coefficient **
** FORWARD_WRITE       use if writing out forward solution, basic state      **
** FORWARD_READ        use if reading in  forward solution, basic state      **
** FORWARD_RHS         use if processing forward right-hand-side terms       **
** IMPLICIT_VCONV      use if implicit vertical convolution algorithm        **
** IMPULSE             use if processing adjoint impulse forcing             **
** MINRES              use if Minimal Residual Method for 4DVar minimization **
** MULTIPLE_TLM        use if multiple TLM history files in 4DVAR            **
** NLM_OUTER           use if nonlinear model as basic state in outer loop   **
** OBS_IMPACT          use if observation impact to 4DVAR data assimilation  **
** OBS_IMPACT_SPLIT    use to separate impact due to IC, forcing, and OBC    **
** POSTERIOR_EOFS      use if posterior analysis error covariance EOFS       **
** POSTERIOR_ERROR_F   use if final posterior analysis error covariance      **
** POSTERIOR_ERROR_I   use if initial posterior analysis error covariance    **
** RECOMPUTE_4DVAR     use if recomputing 4DVar in analysis algorithms       **
** RPM_RELAXATION      use if Picard iterations, Diffusive Relaxation of RPM **
** SO_SEMI_WHITE       use to activate SO semi norm white/red noise processes**
** STOCH_OPT_WHITE     use to activate SO white/red noise processes          **
** SPLINES_VCONV       use to activate implicit splines vertical convolution **
** VCONVOLUTION        use to add vertical correlation to 3D convolution     **
** VERIFICATION        use if writing out solution at observation locations  **
** ZETA_ELLIPTIC       use if SSH elliptic Equation in balance operator      **
**                                                                           **
** OPTION for processing the full grid range (interior and boundary points)  **
** of the state vector in variational data assimilation and generalized      **
** stability theory analysis. Otherwise, only interior points are processed. **
**                                                                           **
** FULL_GRID           use to consider both interior and boundary points     **
**                                                                           **
** Fennel et al. (2006) biology model OPTIONS:                               **
**                                                                           **
** BIO_FENNEL          use if Fennel et al. (2006) nitrogen-based model      **
** BIO_SEDIMENT        use to restore fallen material to the nutrient pool   **
** CARBON              use to add carbon constituents                        **
** DENITRIFICATION     use to add denitrification processes                  **
** OXYGEN              use to add oxygen dynamics                            **
** OCMIP_OXYGEN_SC     use if Schmidt number from Keeling et al. (1998)      **
** TALK_NONCONSERV     use if nonconservative computation of alkalinity      **
**                                                                           **
** Spectral Fennel biology model OPTIONS:                                    **
**                                                                           **
** SPECTRAL_LIGHT      use for spectral light as part of Fennel. Gallegos    **
** CDOM_DEFAULT        use for constant default CDOM for spectral light      **
** CDOM_VARIABLE       use for variable CDOM for spectral light              **
** MOD_SWR_SPECT       modulate shortwave radiation with Gallegos spectrum   **
** MOD_SWR_HOMO        modulate shortwave radiation evenly through spectrum  **
** CHL_BACKSCAT        add chlorophyll backscatter in attenuation            **
**                                                                           **
** Vegetation growth and flow model OPTIONS:                                 **
**                                                                           **
** VEGETATION          use to activate submerged/emergent vegetation effects **
** SAV                 activate submerged vegetation (seagrass) effects      **
** EMERGENT_VEG        activate emergent vegetation (marsh) effects          **
** SEAGRASS_SINK       use bottom sink of N due to benthic seagrass proxy    **
** SEAGRASS_LIGHT      add N sink (seagrass proxy) as function of light      **
** SEAGRASS_LIGHT_CONST constant N sink (seagrass proxy) if light exceeded   **
**                                                                           **
** Bering Sea biology model OPTIONS:                                         **
**                                                                           **
** BEST_NPZ            use if Gibson et al. Bering Sea model                 **
** STATIONARY          use if ??                                             **
** BENTHIC             use if benthic components                             **
** ICE_BIO             use if ice algae                                      **
** JELLY               use if jellyfish                                      **
** CLIM_ICE_1D         use if one-D with ice                                 **
**                                                                           **
** NPZD biology model OPTIONS:                                               **
**                                                                           **
** BIO_UMAINE          use if Chai et al. (2002) CoSINE model                **
** NPZD_FRANKS         use if NPZD Biology model, Franks et al. (1986)       **
** NPZD_IRON           use if NPZD Biology model with iron limitation        **
** NPZD_POWELL         use if NPZD Biology model, Powell et al. (2006)       **
** IRON_LIMIT          use if Fe limitation on phytoplankton growth          **
** IRON_RELAX          use if nudging Fe over the shelf, h <= FeHmin         **
**                                                                           **
** Bio-optical EcoSim model OPTIONS:                                         **
**                                                                           **
** ECOSIM              use if bio-optical EcoSim model                       **
**                                                                           **
** Nemuro lower trophic level ecosystem model OPTIONS:                       **
**                                                                           **
**    Need to choose a zooplankton grazing option (HOLLING_GRAZING or        **
**    IVLEV_EXPLICIT). The default implicit IVLEV algorithm does not         **
**    work yet.                                                              **
**                                                                           **
** NEMURO              use if Nemuro ecosystem model.                        **
** NEMURO_SED1         use if Nemuro sediment remineralization               **
** PRIMARY_PROD        use if primary productivity output                    **
** BIO_SEDIMENT        use to restore fallen material to the nutrient pool   **
** HOLLING_GRAZING     use Holling-type s-shaped curve grazing (implicit)    **
** IVLEV_EXPLICIT      use Ivlev explicit grazing algorithm                  **
**                                                                           **
** Red tide biological model OPTIONS:                                        **
**                                                                           **
** RED_TIDE            use if red tide biological model.                     **
**                                                                           **
** Sediment transport model OPTIONS:                                         **
**                                                                           **
** SEDIMENT            use to activate sediment transport model              **
** BEDLOAD_MPM         use to activate Meyer-Peter-Mueller bed load          **
** BEDLOAD_SOULSBY     use to activate Soulsby wave/current bed load         **
** SED_DENS            use to activate sediment to affect equation of state  **
** SED_MORPH           use to allow bottom model elevation to evolve         **
** SUSPLOAD            use to activate suspended load transport              **
** SED_BIODIFF         use to activate sediment biodiffusivity               **
** MIXED_BED           use to activate mixed bed behavior                    **
** COHESIVE_BED        use to activate cohesive bed model                    **
** NONCOHESIVE_BED1    use original bed model of Warner et al 2008, default  **
** NONCOHESIVE_BED2    use modified bed model of Sherwood et al, in press    **
** SED_FLOCS           use flocculation model of Verney et al., 2011         **
** FLOC_TURB_DISS      use dissipation for flocculation based on turbulence  **
** FLOC_BBL_DISS       use dissipation for flocs from bottom boundary layer  **
** SED_DEFLOC          use flocculation decomposition in sediment bed        **
** SED_TAU_CD_CONST    use constant critical stress for deposition           **
** SED_TAU_CD_LIN      use linear critical stress for deposition             **
**                                                                           **
** Wave effoct on currents (WEC) and shallow water OPTIONS:                  **
**                                                                           **
** WET_DRY             activate wetting and drying                           **
** WEC_MELLOR          activate radiation stress terms from Mellor 08.       **
** WEC_VF              activate wave-current stresses from Uchiyama et al.   **
** WDISS_THORGUZA      activate wave dissipation from Thorton/Guza.          **
** WDISS_CHURTHOR      activate wave dissipation from Church/Thorton.        **
** WDISS_WAVEMOD       activate wave dissipation from a wave model           **
** WDISS_INWAVE        activate wave dissipation from a InWave model         **
** ROLLER_SVENDSEN     activate wave roller based on Svendsen                **
** ROLLER_MONO         activate wave roller for monchromatic waves           **
** ROLLER_RENIERS      activate wave roller based on Reniers                 **
** BOTTOM_STREAMING    activate wave enhanced bottom streaming               **
** SURFACE_STREAMING   activate wave enhanced surface streaming              **
** WAVE_MIXING         activate enhanced vertical viscosity mixing from waves**
**                                                                           **
** OPTIONS for grid nesting:                                                 **
**                                                                           **
** NESTING             use to activate grid nesting: composite/refinement    **
** NO_CORRECT_TRACER   use to avoid two-way correction of boundary tracer    **
** ONE_WAY             use if one-way nesting in refinement grids            **
** TIME_INTERP_FLUX    time interpolate coarse mass flux instead of persist  **
**                                                                           **
** NetCDF input/output OPTIONS:                                              **
**                                                                           **
** DEFLATE             use to set compression NetCDF-4/HDF5 format files     **
** HDF5                use to create NetCDF-4/HDF5 format files              **
** NO_LBC_ATT          use to not check NLM_LBC global attribute on restart  **
** NO_READ_GHOST       use to not include ghost points during read/scatter   **
** NO_WRITE_GRID       use if not writing grid arrays                        **
** PARALLEL_IN         use if parallel Input via HDF5 or pnetcdf libraries   **
** PARALLEL_OUT        use if parallel Output via HDF5 or pnetcdf libraries  **
** PERFECT_RESTART     use to include perfect restart variables              **
** PNETCDF             use if parallel I/O with pnetcdf (classic format)     **
** POSITIVE_ZERO       use to impose positive zero in ouput data             **
** READ_WATER          use if only reading water points data                 **
** WRITE_WATER         use if only writing water points data                 **
** RST_SINGLE          use if writing single precision restart fields        **
** OUT_DOUBLE          use if writing double precision output fields         **
**                                                                           **
** Coupling Library OPTIONS:                                                 **
**                                                                           **
** ESMF_LIB            use Earth System Modeling Framework Library           **
** MCT_LIB             use Model Coupling Toolkit Library                    **
**                                                                           **
** OPTION to process 3D data by levels (2D slabs) to reduce memory needs in  **
** distributed-memory configurations. This option is convenient for large    **
** problems on nodes with limited memory.                                    **
**                                                                           **
** INLINE_2DIO         use if processing 3D IO level by level                **
**                                                                           **
** Sea ice options
**
** ICE_MODEL           use for sea-ice model
** ICE_THERMO          use for thermodynamic component
** ICE_MK              use for Mellor-Kantha thermodynamics (no other choice)
** ICE_ALB_EC92        use for albedo computation from Ebert and Curry
** ICE_MOMENTUM        use for momentum component
** ICE_MOM_BULK        hmmm, some option for ice-water stress computation
** ICE_EVP             use for elastic-viscous-plastic rheology
** ICE_I_O             use for shortwave going to heat ice
** ICE_SHOREFAST       use for a simple shorefast-ice algorithm (Budgell)
** ICE_LANDFAST        use for a different shorefast-ice algorithm (Lemieux et al)
** ICE_ADVECT          use for advection of ice tracers
** ICE_SMOLAR          use for MPDATA advection scheme
** ICE_UPWIND          use for upwind advection scheme
** ICE_BULK_FLUXES     use for ice part of bulk flux computation
** ICE_CONVSNOW        use for conversion of flooded snow to ice
** ICE_STRENGTH_QUAD   use for ice strength a quadratic function of thickness
** INI_GLORYS_ICE      use for ice initial conditions from GLORYS
** NO_SCORRECTION_ICE  use for no scorrection under the ice
** OUTFLOW_MASK        use for Hibler style outflow cells
**
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
** DOGBONE             Idealize nesting grids (Composite and Refinement) Test**
** DOUBLE_GYRE         Idealized Double-gyre Example                         **
** ESTUARY_TEST        Test Estuary for Sediment                             **
** FLT_TEST            Float Tracking Example                                **
** GRAV_ADJ            Gravitational Adjustment Example                      **
** INLET_TEST          Test Inlet Application                                **
** ISOMIP              ISOMIP 2.01 Ice Shelf test case                       **
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
