/*
** Include file "globaldef.h"
**
** svn $Id: globaldefs.h 1054 2021-03-06 19:47:12Z arango $
********************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group     Alexander F. Shchepetkin  **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**                                                                           **
** WARNING: This  file  contains  a set of  predetermined macro definitions  **
** =======  which are inserted into the individual files by C-preprocessor.  **
** It is strongly recommended to NOT modify any of the definitions below.    **
**                                                                           **
*******************************************************************************
*/

/*
** Set assumed-shape array switch.  Imported arrays with dummy
** arguments that takes the shape of the actual argument passed
** to it.  If off, all the arrays are explicit-shape.  In some
** computer explicit-shape arrays slow down performacnce because
** the arrays are copied when passed by arguments.
*/

#if !((defined G95 && defined I686) || defined UNICOS_SN)
# define ASSUMED_SHAPE
#endif

/*
** Set switch for computer lacking 4-byte (32 bit) floating point
** representation, like some Crays.  This becomes important when
** defining attributes for 4-byte float variables in NetCDF files.
** We need to have the _FillValue attribute of the same type as
** as the NetCDF variable.
*/

#if defined UNICOS_SN
# define NO_4BYTE_REALS
#endif

/*
** Set internal distributed-memory switch.
*/

#if defined MPI
# define DISTRIBUTE
#endif

/*
** Make sure that either "mpi_allgather" or "mpi_allreduce" is used
** in mp_reduce.  Low-level routines give an error.
*/

#ifdef DISTRIBUTE
# if !(defined REDUCE_ALLGATHER || defined REDUCE_ALLREDUCE)
#  define REDUCE_ALLREDUCE
# endif
#endif

/*
** Turn ON/OFF time profiling.
*/

#define PROFILE

/*
** Set default time-averaging filter for barotropic fields.
**
*/

#ifdef SOLVE3D
# undef COSINE2
# define POWER_LAW
#endif

/*
** Turn ON/OFF switch to include/disregard the difference between
** rho0 and surface density in the computation of baroclinic pressure
** term.
*/

#define RHO_SURF

/*
** Turn ON/OFF double precision arithmetic in numerical kernel (default)
** and floating-point type variables and associated intrinsic functions.
*/

#ifdef SINGLE_PRECISION
# ifdef OUT_DOUBLE
#   undef OUT_DOUBLE
# endif
# ifndef RST_SINGLE
#   define RST_SINGLE
# endif
#else
# define DOUBLE_PRECISION
#endif

/*
** Turn ON masking when wetting and drying is activated.
*/

#if !defined MASKING && defined WET_DRY
# define MASKING
#endif

/*
** If wetting and drying, activate limiting of bottom stress.
** For COAWST, I use WET_DRY in place of the LIMIT option.
*/

#if !defined LIMIT_BSTRESS && defined WET_DRY
# define LIMIT_BSTRESS
#endif

/*
** Define macro for the first 2D time-step.
*/

#ifdef SOLVE3D
# define FIRST_2D_STEP iif(ng).eq.1
#else
# define FIRST_2D_STEP iic(ng).eq.ntfirst(ng)
#endif

/*
** Set switch for distributed-memory applications to gather and scatter
** I/O data in 2D slabs. This is necessary on some platforms to conserve
** memory.
*/

#if defined DISTRIBUTE
# if defined UNICOS_SN
#  define INLINE_2DIO
# endif
#endif

/*
** Remove OpenMP directives in serial and distributed memory
** Applications.  This definition will be used in conjunction with
** the pearl script "cpp_clean" to remove the full directive.
*/

#if !defined _OPENMP
# define $OMP !
#endif

/*
** Set tile range for distributed- or shared-memory configurations.
*/

#ifdef DISTRIBUTE
# define THREAD_RANGE MyRank,MyRank
#else
# define THREAD_RANGE 0,numthreads-1
#endif

/*
** Choice of double/single precision for real type variables and
** associated intrinsic functions.
*/

#if (defined CRAY || defined CRAYT3E) && !defined CRAYX1
# ifdef  DOUBLE_PRECISION
#  undef  DOUBLE_PRECISION
# endif
#endif

#ifdef DOUBLE_PRECISION
# ifdef DISTRIBUTE
#  define PDNAUPD pdnaupd
#  define PDNEUPD pdneupd
#  define PDSAUPD pdsaupd
#  define PDSEUPD pdseupd
#  define PDNORM2 pdnorm2
# else
#  define DNAUPD dnaupd
#  define DNEUPD dneupd
#  define DSAUPD dsaupd
#  define DSEUPD dseupd
#  define DNRM2 dnrm2
# endif
# define DAXPY daxpy
# define DSTEQR dsteqr
#else
# ifdef DISTRIBUTE
#  define PDNAUPD psnaupd
#  define PDNEUPD psneupd
#  define PDSAUPD pssaupd
#  define PDSEUPD psseupd
#  define PDNORM2 psnorm2
# else
#  define DNAUPD snaupd
#  define DNEUPD sneupd
#  define DSAUPD ssaupd
#  define DSEUPD sseupd
#  define DNRM2 snrm2
# endif
#  define DAXPY saxpy
#  define DSTEQR ssteqr
#endif

/*
** Single intrinsic Fortran functions.
*/

#ifdef SINGLE_PRECISION
# define DSIGN SIGN
#endif

/*
** Check for 4D-Var deprecated CPP options.
*/

#ifdef IS4DVAR
# ifndef I4DVAR
#  define I4DVAR
# endif
#endif

#ifdef IS4DVAR_SENSITIVITY
# ifndef I4DVAR_ANA_SENSITIVITY
#  define I4DVAR_ANA_SENSITIVITY
# endif
#endif

#ifdef W4DPSAS
# ifndef RBL4DVAR
#  define RBL4DVAR
# endif
#endif

#ifdef W4DPSAS_SENSITIVITY
# ifndef RBL4DVAR_ANA_SENSITIVITY
#  define RBL4DVAR_ANA_SENSITIVITY
# endif
#endif

#ifdef W4DPSAS_FCT_SENSITIVITY
# ifndef RBL4DVAR_FCT_SENSITIVITY
#  define RBL4DVAR_FCT_SENSITIVITY
# endif
#endif

#ifdef W4DVAR
# ifndef R4DVAR
#  define R4DVAR
# endif
#endif

#ifdef W4DVAR_SENSITIVITY
# ifndef R4DVAR_ANA_SENSITIVITY
#  define R4DVAR_ANA_SENSITIVITY
# endif
#endif

/*
** Set multiple excutables split 4D-Var.
*/

#if defined SPLIT_I4DVAR   || \
    defined SPLIT_RBL4DVAR || \
    defined SPLIT_R4DVAR   || \
    defined SPLIT_SP4DVAR
# define SPLIT_4DVAR
#endif

/*
** If split 4D-Var, activate the unsplit option since both share
** identical configuration to avoid too many directives changes.
*/

#if !defined I4DVAR && defined SPLIT_I4DVAR
# define I4DVAR
#endif

#if !defined RBL4DVAR && defined SPLIT_RBL4DVAR
# define RBL4DVAR
#endif

#if !defined R4DVAR && defined SPLIT_R4DVAR
# define R4DVAR
#endif

#if !defined SP4DVAR && defined SPLIT_SP4DVAR
# define SP4DVAR
#endif


/*
** Set 4D-Var sensitivity switch.
*/

#if defined RBL4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR_FCT_SENSITIVITY || \
    defined R4DVAR_ANA_SENSITIVITY
# define SENSITIVITY_4DVAR
#endif

#if (defined RBL4DVAR       && \
     defined OBS_SPACE)     || \
    (defined SPLIT_RBL4DVAR && \
     defined OBS_SPACE)
# undef OBS_SPACE
#endif

/*
** Set perturbation tangent linear, fine amplitude tangent linear,
** and adjoint model switches.
*/

#if defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined CORRELATION            || \
    defined FORCING_SV             || \
    defined FT_EIGENMODES          || \
    defined HESSIAN_FSV            || \
    defined HESSIAN_SO             || \
    defined HESSIAN_SV             || \
    defined INNER_PRODUCT          || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_PERTURBATION       || \
    defined OPT_OBSERVATIONS       || \
    defined PICARD_TEST            || \
    defined RBL4DVAR               || \
    defined RPM_DRIVER             || \
    defined R4DVAR                 || \
    defined R_SYMMETRY             || \
    defined SANITY_CHECK           || \
    defined SENSITIVITY_4DVAR      || \
    defined SPLIT_I4DVAR           || \
    defined SPLIT_RBL4DVAR         || \
    defined SPLIT_R4DVAR           || \
    defined SPLIT_SP4DVAR          || \
    defined SP4DVAR                || \
    defined STOCHASTIC_OPT         || \
    defined TLM_CHECK              || \
    defined TLM_DRIVER             || \
    defined TL_RBL4DVAR            || \
    defined TL_R4DVAR
# define TANGENT
#endif

#if defined AD_SENSITIVITY         || \
    defined ADM_DRIVER             || \
    defined AFT_EIGENMODES         || \
    defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined CORRELATION            || \
    defined FORCING_SV             || \
    defined HESSIAN_SO             || \
    defined HESSIAN_FSV            || \
    defined HESSIAN_SV             || \
    defined INNER_PRODUCT          || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_PERTURBATION       || \
    defined OPT_OBSERVATIONS       || \
    defined RBL4DVAR               || \
    defined R4DVAR                 || \
    defined R_SYMMETRY             || \
    defined SANITY_CHECK           || \
    defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI                || \
    defined SPLIT_I4DVAR           || \
    defined SPLIT_RBL4DVAR         || \
    defined SPLIT_R4DVAR           || \
    defined SPLIT_SP4DVAR          || \
    defined SP4DVAR                || \
    defined STOCHASTIC_OPT         || \
    defined TLM_CHECK              || \
    defined TL_RBL4DVAR            || \
    defined TL_R4DVAR
# define ADJOINT
#endif

#if defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined PICARD_TEST            || \
    defined RPM_DRIVER             || \
    defined TL_R4DVAR              || \
    defined R4DVAR                 || \
    defined R4DVAR_ANA_SENSITIVITY
# define TL_IOMS
#endif

#if !defined ANA_PERTURB     && \
    (defined CORRELATION     || \
     defined SANITY_CHECK    || \
     defined R_SYMMETRY)
# define ANA_PERTURB
#endif

/*
** Turn off nonlinear model switch.
*/

#define NONLINEAR
#if defined AD_SENSITIVITY   || \
    defined ADM_DRIVER       || \
    defined AFT_EIGENMODES   || \
    defined FORCING_SV       || \
    defined FT_EIGENMODES    || \
    defined HESSIAN_FSV      || \
    defined HESSIAN_SO       || \
    defined HESSIAN_SV       || \
    defined INNER_PRODUCT    || \
    defined OPT_OBSERVATIONS || \
    defined OPT_PERTURBATION || \
    defined PICARD_TEST      || \
    defined RPM_DRIVER       || \
    defined SANITY_CHECK     || \
    defined SO_SEMI          || \
    defined STOCHASTIC_OPT   || \
    defined TLM_DRIVER
# undef NONLINEAR
#endif

/*
** Activate bacroclinic pressure gradient response due to the
** perturbation of free-surface in the presence of stratification
** and bathymetry. This option does not pass the sanity check
** in adjoint and tangent linear applications.
*/

#ifdef SOLVE3D
# if !(defined ADJOINT || defined TANGENT)
#   define VAR_RHO_2D
# endif
#endif

/*
** Set output index for multi-time levels variables.
*/

#ifdef SOLVE3D
# if defined TANGENT || defined TL_IOMS
#  define TKOUT kstp(ng)
#  define KOUT kstp(ng)
#  define NOUT nrhs(ng)
# else
#  define KOUT kstp(ng)
#  define NOUT nrhs(ng)
# endif
#else
# if defined TANGENT || defined TL_IOMS
#  define TKOUT kstp(ng)
# endif
# define KOUT knew(ng)
#endif

/*
** Set internal switch for the need of a propagator driver.
*/

#if defined AFT_EIGENMODES   || \
    defined ENSEMBLE         || \
    defined FORCING_SV       || \
    defined FT_EIGENMODES    || \
    defined HESSIAN_FSV      || \
    defined HESSIAN_SO       || \
    defined HESSIAN_SV       || \
    defined OPT_PERTURBATION || \
    defined PSEUDOSPECTRA    || \
    defined SO_SEMI          || \
    defined SO_TRACE         || \
    defined STOCHASTIC_OPT
# define PROPAGATOR
#endif

/*
** Activate processing of forward vertical mixing.
**
*/

#if !defined FORWARD_MIXING  && \
    (defined TANGENT         || \
     defined TL_IOMS         || \
     defined ADJOINT)        && \
    (defined LMD_MIXING      || \
     defined GLS_MIXING      || \
     defined MY25_MIXING)
# define FORWARD_MIXING
#endif

/*
** Set internal switches for all the 4DVAR schemes.
*/

#if !defined WEAK_CONSTRAINT          && \
    (defined ARRAY_MODES              || \
     defined CLIPPING                 || \
     defined SP4DVAR                  || \
     defined RBL4DVAR                 || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined R4DVAR                   || \
     defined R4DVAR_ANA_SENSITIVITY   || \
     defined R_SYMMETRY               || \
     defined SPLIT_RBL4DVAR           || \
     defined SPLIT_R4DVAR             || \
     defined SPLIT_SP4DVAR            || \
     defined SP4DVAR                  || \
     defined TL_RBL4DVAR              || \
     defined TL_R4DVAR)
# define WEAK_CONSTRAINT
#endif

#if !defined WEAK_CONSTRAINT && defined RPM_RELAXATION
# undef RPM_RELAXATION
#endif

#if defined CORRELATION            || \
    defined HESSIAN_FSV            || \
    defined HESSIAN_SO             || \
    defined HESSIAN_SV             || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS       || \
    defined TLM_CHECK              || \
    defined WEAK_CONSTRAINT
# define FOUR_DVAR
#endif

#if defined I4DVAR
# define BACKGROUND
#endif

#if !(defined RBL4DVAR || defined R4DVAR) && defined POSTERIOR_EOFS
# undef POSTERIOR_EOFS
#endif

#if !(defined RBL4DVAR || defined R4DVAR) && defined POSTERIOR_ERROR_F
# undef POSTERIOR_ERROR_F
#endif

#if !(defined RBL4DVAR || defined R4DVAR) && defined POSTERIOR_ERROR_I
# undef POSTERIOR_ERROR_I
#endif

#if !(defined WEAK_CONSTRAINT || defined I4DVAR_ANA_SENSITIVITY) && \
      defined OBS_IMPACT
# undef OBS_IMPACT
#endif

#if !(defined OBS_IMPACT                && \
      (defined I4DVAR_ANA_SENSITIVITY   || \
       defined RBL4DVAR_ANA_SENSITIVITY || \
       defined RBL4DVAR_FCT_SENSITIVITY || \
       defined R4DVAR_ANA_SENSITIVITY))
# undef IMPACT_INNER
#endif

/*
** Activate internal switch to process 4DVAR observations.
*/

#if defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR               || \
    defined R4DVAR                 || \
    defined SENSITIVITY_4DVAR      || \
    defined SPLIT_I4DVAR           || \
    defined SPLIT_RBL4DVAR         || \
    defined SPLIT_R4DVAR           || \
    defined SPLIT_SP4DVAR          || \
    defined SP4DVAR                || \
    defined TLM_CHECK              || \
    defined TL_RBL4DVAR            || \
    defined TL_R4DVAR              || \
    defined VERIFICATION
# define OBSERVATIONS
#endif

#if defined ARRAY_MODES            || \
    defined CLIPPING               || \
    defined I4DVAR                 || \
    defined I4DVAR_ANA_SENSITIVITY || \
    defined RBL4DVAR               || \
    defined R4DVAR                 || \
    defined R_SYMMETRY             || \
    defined SENSITIVITY_4DVAR      || \
    defined SP4DVAR                || \
    defined TLM_CHECK              || \
    defined TL_RBL4DVAR            || \
    defined TL_R4DVAR
# define TLM_OBS
#endif

/*
** Activate reading and writting of the basic sate.
*/

#if !defined FORWARD_READ           && \
    (defined ARRAY_MODES            || \
     defined CLIPPING               || \
     defined I4DVAR                 || \
     defined I4DVAR_ANA_SENSITIVITY || \
     defined PROPAGATOR             || \
     defined RBL4DVAR               || \
     defined R4DVAR                 || \
     defined SENSITIVITY_4DVAR      || \
     defined SPLIT_I4DVAR           || \
     defined SPLIT_RBL4DVAR         || \
     defined SPLIT_R4DVAR           || \
     defined SPLIT_SP4DVAR          || \
     defined SP4DVAR                || \
     defined TL_RBL4DVAR            || \
     defined TL_R4DVAR)
# define FORWARD_READ
#endif

#if !defined FORWARD_WRITE          && \
    (defined ARRAY_MODES            || \
     defined CLIPPING               || \
     defined I4DVAR                 || \
     defined I4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR               || \
     defined R4DVAR                 || \
     defined SENSITIVITY_4DVAR      || \
     defined SPLIT_I4DVAR           || \
     defined SPLIT_RBL4DVAR         || \
     defined SPLIT_R4DVAR           || \
     defined SPLIT_SP4DVAR          || \
     defined SP4DVAR                || \
     defined TL_RBL4DVAR            || \
     defined TL_R4DVAR)
# define FORWARD_WRITE
#endif

/*
** Set internal weak constraint switches.
*/

#ifdef WEAK_CONSTRAINT
# define IMPULSE
#endif

/*
** Set in internal switch to activate computation of nonlinear
** equation of state expnasion polynomial T-derivatives.
*/

#if defined ADJOINT     || \
    defined BULK_FLUXES || \
    defined LMD_BKPP    || \
    defined LMD_SKPP    || \
    defined TANGENT     || \
    defined TL_IOMS
# define EOS_TDERIVATIVE
#endif

/*
** If splines, deactivate horizontal and vertical smoothing of
** Richardson number horizontally and/or vertically.
*/

#ifdef RI_SPLINES
# if defined LMD_MIXING
#  undef RI_HORAVG
#  undef RI_VERAVG
# endif
#endif

/*
** Activate internal switch for the computation of the Brunt-Vaisala
** frequency.
*/

#if defined BVF_MIXING  || \
    defined GLS_MIXING  || \
    defined LMD_MIXING  || \
    defined LMD_SKPP    || \
    defined LMD_BKPP    || \
    defined MY25_MIXING
# define BV_FREQUENCY
#endif

/*
** Activate internal switch for bottom boundary layer closure.
*/

#if defined MB_BBL  || \
    defined SG_BBL  || \
    defined SSW_BBL
# define BBL_MODEL
#endif

/*
** Check if spatially varying bottom friction parameters are needed.
*/

#if defined UV_DRAG_GRID && \
  !(defined BBL_MODEL    || \
    defined SEDIMENT     || \
    defined UV_LOGDRAG   || \
    defined UV_LDRAG     || \
    defined UV_QDRAG)
# undef UV_DRAG_GRID
#endif

#if defined ANA_DRAG     && \
   !defined UV_DRAG_GRID
# undef ANA_DRAG
#endif


/*
** Activate internal biology option when using any type of biological
** module.
*/

#if defined BIO_FENNEL  || defined ECOSIM      || \
    defined ESTUARYBGC  || defined NEMURO      || \
    defined HYPOXIA_SRM || defined NPZD_FRANKS || \
    defined NPZD_IRON   || defined NPZD_POWELL || \
    defined BIO_UMAINE  || defined BEST_NPZ    || \
    defined RED_TIDE
# define BIOLOGY
#endif

/*
** Use AVERAGES structures to write out filtered data
 */

#ifdef FILTERED
# define AVERAGES
# undef  FILTRIM            /* define for fewer time-filtered fields */
# define FILTERED_RST       /* define if use restart files for time-filtering */
#endif

/*
** Activate internal option for biological float behavior.
*/

#if defined FLOATS       && \
    defined FLOAT_OYSTER
# define FLOAT_BIOLOGY
#endif

/*
** Define internal option to couple to other models.
**
*/

/* this allows differentiation between Rutgers
   and coawst coupling methods */
#define COAWST_MODEL

#if defined ROMS_MODEL && (defined SWAN_MODEL || defined WRF_MODEL || \
                           defined WW3_MODEL  || defined WRFHYDRO_MODEL)
# define ROMS_COUPLING
#endif
#if defined SWAN_MODEL && (defined ROMS_MODEL || defined WRF_MODEL)
# define SWAN_COUPLING
#endif
#if defined WW3_MODEL && (defined ROMS_MODEL || defined WRF_MODEL)
# define WW3_COUPLING
#endif
#if defined WRF_MODEL && (defined SWAN_MODEL || defined WW3_MODEL || \
                          defined ROMS_MODEL)
# define WRF_COUPLING
#endif
#if defined WRFHYDRO_MODEL && (defined SWAN_MODEL || defined WRF_MODEL || \
                               defined WW3_MODEL  || defined ROMS_MODEL)
# define HYDRO_COUPLING
#endif

#if defined WRF_COUPLING && defined ROMS_COUPLING
# define AIR_OCEAN
# ifndef FRC_COUPLING
#  define FRC_COUPLING
# endif
#endif

#if defined HYDRO_COUPLING && defined ROMS_COUPLING
# define HYDRO_OCEAN
#endif

#if defined WRF_COUPLING && (defined SWAN_COUPLING || defined WW3_COUPLING)
# define AIR_WAVES
#endif

#if (defined WW3_COUPLING || defined SWAN_COUPLING) && \
     defined ROMS_COUPLING
# define WAVES_OCEAN
#endif

#if defined AIR_OCEAN || defined AIR_WAVES || defined WAVES_OCEAN || defined HYDRO_OCEAN
# define COAWST_COUPLING
#endif

#if defined COAWST_COUPLING && defined NESTING
# if defined WAVES_OCEAN
#  define MCT_INTERP_OC2WV
# endif
# if defined AIR_OCEAN
#  define MCT_INTERP_OC2AT
# endif
# ifdef AIR_WAVES
#  define MCT_INTERP_WV2AT
# endif
#endif

/* Start of ROMS model coupling cpps.
   COAWST coupling uses options above. These are here for
   consistenct between codes. */
#if defined COAMPS_COUPLING || defined REGCM_COUPLING || \
    defined WRF_COUPLING
# define ATM_COUPLING
# ifndef FRC_COUPLING
#  define FRC_COUPLING
# endif
#endif

#if defined CICE_COUPLING
# define ICE_COUPLING
#endif

#if defined REFDIF_COUPLING || \
    defined SWAN_COUPLING   || \
    defined WAM_COUPLING
# define WAV_COUPLING
#endif

#if defined ATM_COUPLING  || \
    defined DATA_COUPLING || \
    defined ICE_COUPLING  || \
    defined WAV_COUPLING
# define MODEL_COUPLING
#endif

#if defined MODEL_COUPLING && \
    defined ESMF_LIB
# define REGRESS_STARTCLOCK
#endif
/* end of ROMS coupling cpps */


/*
** Define internal option to process wave data.
*/

#if defined WEC_MELLOR || defined WEC_VF
#   define WEC
#endif

/*
** Activate internal switch for imposing REFDIF as a
** monochromatic wave driver.
*/

#if defined REFDIF_COUPLING && defined ROLLER_SVENDSEN
# define ROLLER_MONO
#endif

/*
** Activate internal switch for activating a roller.
*/

#if (defined ROLLER_SVENDSEN || defined ROLLER_MONO ||	\
     defined ROLLER_RENIERS) && defined WEC
# define WEC_ROLLER
#endif


#if defined BBL_MODEL   || defined WEC || \
    defined WAVES_OCEAN || defined INWAVE_MODEL
# define WAVES_DIR
# define WAVES_DIRP
#endif

#if  defined BBL_MODEL   && \
   !(defined SSW_CALC_UB || defined MB_CALC_UB ||  \
     defined SG_CALC_UB)
# define WAVES_UB
#endif

#if (defined BBL_MODEL        && !defined WAVES_UB) ||  \
     defined WEC              || \
     defined INWAVE_MODEL     || \
     defined ZOS_HSIG         || defined COARE_TAYLOR_YELLAND || \
     defined BEDLOAD_SOULSBY  || defined BEDLOAD_VANDERA || \
     defined WAVES_OCEAN      || defined DRENNAN
# define WAVES_HEIGHT
#endif

#if defined WEC || defined BEDLOAD_SOULSBY || defined BEDLOAD_VANDERA || \
    defined WAVES_OCEAN || defined INWAVE_MODEL
# define WAVES_LENGTH
#endif

#if (!defined DEEPWATER_WAVES      && \
     (defined COARE_TAYLOR_YELLAND || defined COARE_OOST)) || \
      defined DRENNAN
# define WAVES_LENGTHP
#endif

#if defined COARE_TAYLOR_YELLAND   || defined COARE_OOST || \
    defined DRENNAN || defined WAVES_OCEAN || defined WEC_VF
# define WAVES_TOP_PERIOD
#endif

#if defined BBL_MODEL || defined WAVES_OCEAN || defined INWAVE_MODEL
# define WAVES_BOT_PERIOD
#endif

#if (defined WDISS_GAMMA || defined WDISS_ROELVINK)
# define WDISS_INWAVE
#endif

#if (defined TKE_WAVEDISS || defined WEC_VF) && \
  (!defined WDISS_THORGUZA && \
   !defined WDISS_CHURTHOR && !defined WDISS_WAVEMOD \
   && !defined WDISS_INWAVE)
# define WAVES_DISS
#endif

#if !defined WAVES_OCEAN     && \
   ((defined BULK_FLUXES     && defined COARE_TAYLOR_YELLAND) || \
    (defined BULK_FLUXES     && defined COARE_OOST)           || \
     defined WEC_ROLLER      || defined WAVES_DISS            || \
     defined WAVES_DIR       || defined WAVES_BOT_PERIOD      || \
     defined WAVES_HEIGHT    || defined WAVES_TOP_PERIOD      || \
     defined WAVES_LENGTH    || defined WAVES_LENGTHP)
# define WAVE_DATA
#endif

/*
** Define internal option for bedload treatment.
*/

#if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY  || defined BEDLOAD_VANDERA
# define BEDLOAD
#endif

/*
** Define internal flag indicating processing of input forcing
** NetCDF file.
*/

#ifdef SOLVE3D
# ifdef BULK_FLUXES
#  ifdef ANA_SMFLUX
#   undef ANA_SMFLUX
#  endif
#  ifdef ANA_STFLUX
#   undef ANA_STFLUX
#  endif
# endif
# if !defined ANA_BTFLUX   || \
    (!defined FRC_COUPLING && \
     !defined BULK_FLUXES  && !defined ANA_SMFLUX)   || \
    (!defined FRC_COUPLING && \
     !defined BULK_FLUXES  && !defined ANA_STFLUX)   || \
    ( defined BIOLOGY      && !defined ANA_SPFLUX)   || \
    ( defined BIOLOGY      && !defined ANA_BPFLUX)   || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined LONGWAVE)     || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_PAIR)     || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_TAIR)     || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_HUMIDITY) || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_CLOUD)    || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_RAIN)     || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_WINDS)    || \
    (!defined FRC_COUPLING && \
      defined BULK_FLUXES  && !defined ANA_SRFLUX)   || \
    (!defined FRC_COUPLING && \
      defined LMD_SKPP     && !defined ANA_SRFLUX)   || \
    (!defined FRC_COUPLING && \
      defined RED_TIDE)    || \
    (!defined ANA_BSFLUX   && defined SALINITY)      || \
    (!defined FRC_COUPLING && \
      defined SALINITY     && !defined ANA_SSFLUX    && \
     (defined BULK_FLUXES  && !defined EMINUSP))     || \
    (!defined FRC_COUPLING && \
      defined SOLAR_SOURCE && !defined ANA_SRFLUX)   || \
    ( defined BBL_MODEL    && !defined ANA_WWAVE     && \
     !defined WAVES_OCEAN  && !defined INWAVE_MODEL) || \
    ( defined SEDIMENT     && !defined ANA_SPFLUX)   || \
    ( defined SEDIMENT     && !defined ANA_BPFLUX)   || \
    ( defined WAVE_DATA    && !defined ANA_WWAVE     && \
     !defined WAVES_OCEAN  && !defined INWAVE_MODEL)
#  define FRC_FILE
# endif
#else
# if(!defined FRC_COUPLING && !defined ANA_SMFLUX)
#  define FRC_FILE
# endif
#endif

/*
** Check analytical initial conditions options.
*/

#if defined ANA_BIOLOGY && \
   !defined BIOLOGY
# undef ANA_BIOLOGY
#endif

#if defined ANA_PASSIVE && \
   !defined T_PASSIVE
# undef ANA_PASSIVE
#endif

#if defined ANA_SEDIMENT && \
   !(defined SEDIMENT    || \
     defined BBL_MODEL)
# undef ANA_SEDIMENT
#endif
#if  !defined ANA_INITIAL || \
    ( defined BIOLOGY     && !defined ANA_BIOLOGY)  || \
    ( defined T_PASSIVE   && !defined ANA_PASSIVE)  || \
    ( defined SEDIMENT    && !defined ANA_SEDIMENT) || \
    ( defined BBL_MODEL   && !defined ANA_SEDIMENT) || \
    ( defined VEGETATION  && !defined ANA_VEGETATION)
# define INI_FILE
#endif

/*
** Define internal shortwave radiation option.  Undefine analytical
** shortwave option if not needed.
*/

#if defined LMD_SKPP     || defined SOLAR_SOURCE   || \
    defined BULK_FLUXES  || defined BIOLOGY        || \
    defined ATM2OCN_FLUXES
# define SHORTWAVE
#endif

#if !defined SHORTWAVE   && \
     defined ANA_SRFLUX
# undef ANA_SRFLUX
#endif

#if !defined SHORTWAVE      && \
     defined DIURNAL_SRFLUX
# undef DIURNAL_SRFLUX
#endif

/*
** Define internal clouds option.  Undefine analytical
** shortwave option if not needed.
*/

#if (defined ANA_SRFLUX   && \
     defined ALBEDO)      || \
    (defined BULK_FLUXES  && \
     defined LONGWAVE)    || \
    defined ECOSIM
# define CLOUDS
#endif

#if !defined CLOUDS    && \
     defined ANA_CLOUD
# undef ANA_CLOUD
#endif

/*
** Check if it is meaningful to write out momentum/tracer diagnostics
** and activate internal diagnostics option.
*/

#if !defined SOLVE3D  || \
     defined TS_FIXED
# if defined DIAGNOSTICS_TS
#   undef DIAGNOSTICS_TS
# endif
#endif

#if defined DIAGNOSTICS_BIO && \
  !(defined BIO_FENNEL      || \
    defined ECOSIM          || \
    defined HYPOXIA_SRM)
#  undef DIAGNOSTICS_BIO
#endif

#if defined DIAGNOSTICS_BIO || \
    defined DIAGNOSTICS_TS  || \
    defined DIAGNOSTICS_UV
# define DIAGNOSTICS
#endif

/*
** Check if any analytical expression is defined.
*/

#if defined ANA_BIOLOGY     || \
    defined ANA_BPFLUX      || \
    defined ANA_BSFLUX      || \
    defined ANA_BTFLUX      || \
    defined ANA_CLOUD       || \
    defined ANA_DIAG        || \
    defined ANA_DQDSST      || \
    defined ANA_DRAG        || \
    defined ANA_FSOBC       || \
    defined ANA_GRID        || \
    defined ANA_HUMIDITY    || \
    defined ANA_INITIAL     || \
    defined ANA_M2CLIMA     || \
    defined ANA_M2OBC       || \
    defined ANA_M3CLIMA     || \
    defined ANA_M3OBC       || \
    defined ANA_MASK        || \
    defined ANA_NUDGCOEF    || \
    defined ANA_PAIR        || \
    defined ANA_PASSIVE     || \
    defined ANA_PERTURB     || \
    defined ANA_PSOURCE     || \
    defined ANA_RAIN        || \
    defined ANA_RESPIRATION || \
    defined ANA_SEDIMENT    || \
    defined ANA_SMFLUX      || \
    defined ANA_SPFLUX      || \
    defined ANA_SPINNING    || \
    defined ANA_SPONGE      || \
    defined ANA_SRFLUX      || \
    defined ANA_SSFLUX      || \
    defined ANA_SSH         || \
    defined ANA_SSS         || \
    defined ANA_SST         || \
    defined ANA_STFLUX      || \
    defined ANA_TAIR        || \
    defined ANA_TCLIMA      || \
    defined ANA_TOBC        || \
    defined ANA_VMIX        || \
    defined ANA_WINDS       || \
    defined ANA_WWAVE       || \
    defined DIFF_GRID       || \
    defined VISC_GRID
# define ANALYTICAL
#endif

/*
** If splitting 3rd-order upstream bias horizontal advection of
** momentum, activate other needed flags.
*/

#ifdef UV_U3ADV_SPLIT
# define VISC_3DCOEF
# ifndef UV_C4ADVECTION
#  define UV_C4ADVECTION
# endif
# ifndef UV_VIS4
#  define UV_VIS4
# endif
# ifdef UV_VIS2
#  undef UV_VIS2
# endif
# ifdef UV_SMAGORINSKY
#  undef UV_SMAGORINSKY
# endif
#endif

/*
** Define internal switch for Smagorinsky-like mixing.
*/

#if !defined DIFF_3DCOEF && defined TS_SMAGORINSKY
# define DIFF_3DCOEF
#endif
#if !defined VISC_3DCOEF && defined UV_SMAGORINSKY
# define VISC_3DCOEF
#endif

/*
** Define internal switch for vegetation horizontal mixing
*/
#if defined VEG_HMIXING
# define VISC_3DCOEF
#endif
