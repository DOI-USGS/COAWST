/*
** Include file "globaldef.h"
**
** svn $Id: globaldefs.h 838 2008-11-17 04:22:18Z jcwarner $
********************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2010 The ROMS/TOMS Group     Alexander F. Shchepetkin  **
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
** If parallel I/O and applicable, turn on NetCDF-4 type files.
*/

#if !defined NETCDF4 && defined PARALLEL_IO
# define NETCDF4
#endif

/*
** Set internal distributed-memory switch.
*/

#if defined MPI
# define DISTRIBUTE
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
** Activate criteria for isopycnic diffusion of tracer as maximum
** density slope or minimum density stratification.  Choose only
** one option. If neither option is activated, the default criteria
** is used in the algorithm.
*/

#if defined MIX_ISO_TS && (defined TS_DIF2 || defined TS_DIF4)
# undef   MAX_SLOPE
# undef   MIN_STRAT
#endif

/*
** Turn ON/OFF double precision for real type variables and
** associated intrinsic functions.
*/

#define DOUBLE_PRECISION

/*
** Turn ON masking when wetting and drying is activated.
*/

#if !defined MASKING && defined WET_DRY
# define MASKING
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
** Set number of ghost-points in the halo region.
*/

#if defined EW_PERIODIC && defined REFINED_GRID
# define EW_PERIODIC_REFINED
#endif

#if defined TS_MPDATA || defined UV_VIS4 || defined COMPOSED_GRID || \
    defined REFINED_GRID
# define GHOST_POINTS 3
# if defined DISTRIBUTE || defined EW_PERIODIC || \
     defined NS_PERIODIC || defined COMPOSED_GRID || \
     defined REFINED_GRID
#   define THREE_GHOST
# endif
#else
# define GHOST_POINTS 2
#endif

/*
** Define global grid lower and upper bounds in the I- and
** J-directions. These values are a function of periodicity.
** They are used in both shared- and distributed-memory
** configurations.
*/

#if defined COMPOSED_GRID || defined REFINED_GRID
#  define LOWER_BOUND_I -GHOST_POINTS
/*#  define UPPER_BOUND_I Im(ng)+GHOST_POINTS*/
#  define UPPER_BOUND_I Lm(ng)+GHOST_POINTS
#  define LOWER_BOUND_J -GHOST_POINTS
/*#  define UPPER_BOUND_J Jm(ng)+GHOST_POINTS*/
#  define UPPER_BOUND_J Mm(ng)+GHOST_POINTS

#  define LOWER_BOUND_Ip -GHOST_POINTS
#  define UPPER_BOUND_Ip Lm(ngp)+GHOST_POINTS
#  define LOWER_BOUND_Jp -GHOST_POINTS
#  define UPPER_BOUND_Jp Mm(ngp)+GHOST_POINTS
#  define LOWER_BOUND_Ic -GHOST_POINTS
#  define UPPER_BOUND_Ic Lm(ngc)+GHOST_POINTS
#  define LOWER_BOUND_Jc -GHOST_POINTS
#  define UPPER_BOUND_Jc Mm(ngc)+GHOST_POINTS
#else
# ifdef EW_PERIODIC
#  ifdef NS_PERIODIC
#   define LOWER_BOUND_I -GHOST_POINTS
#   define UPPER_BOUND_I Im(ng)+GHOST_POINTS
#   define LOWER_BOUND_J -GHOST_POINTS
#   define UPPER_BOUND_J Jm(ng)+GHOST_POINTS
#  else
#   define LOWER_BOUND_I -GHOST_POINTS
#   define UPPER_BOUND_I Im(ng)+GHOST_POINTS
#   define LOWER_BOUND_J 0
#   define UPPER_BOUND_J Jm(ng)+1
#  endif
# else
#  ifdef NS_PERIODIC
#   define LOWER_BOUND_I 0
#   define UPPER_BOUND_I Im(ng)+1
#   define LOWER_BOUND_J -GHOST_POINTS
#   define UPPER_BOUND_J Jm(ng)+GHOST_POINTS
#  else
#   define LOWER_BOUND_I 0
#   define UPPER_BOUND_I Im(ng)+1
#   define LOWER_BOUND_J 0
#   define UPPER_BOUND_J Jm(ng)+1
#  endif
# endif
#endif
#define XI_DIM  LOWER_BOUND_I:UPPER_BOUND_I
#define ETA_DIM LOWER_BOUND_J:UPPER_BOUND_J
#define GLOBAL_2D_ARRAY XI_DIM,ETA_DIM
#ifdef REFINED_GRID
# define XI_DIMp  LOWER_BOUND_Ip:UPPER_BOUND_Ip
# define ETA_DIMp LOWER_BOUND_Jp:UPPER_BOUND_Jp
# define XI_DIMc  LOWER_BOUND_Ic:UPPER_BOUND_Ic
# define ETA_DIMc LOWER_BOUND_Jc:UPPER_BOUND_Jc
#endif

/* #define PRIVATE_1D_SCRATCH_ARRAY Istr-3:Iend+3
   #define PRIVATE_2D_SCRATCH_ARRAY Istr-3:Iend+3,Jstr-3:Jend+3 */
/*#ifdef REFINED_GRID
  # define PRIVATE_1D_SCRATCH_ARRAY Istr-4:Iend+3
  # define PRIVATE_2D_SCRATCH_ARRAY Istr-4:Iend+3,Jstr-4:Jend+3
  #else */
/*#define PRIVATE_1D_SCRATCH_ARRAY LBi-1:UBi+1
  #define PRIVATE_2D_SCRATCH_ARRAY LBi-1:UBi+1,LBj-1:UBj+1 */
#define PRIVATE_1D_SCRATCH_ARRAY IminS:ImaxS
#define PRIVATE_2D_SCRATCH_ARRAY IminS:ImaxS,JminS:JmaxS
#define PRIVATE_2D_SCRATCH_ARRAY_THETA IminS:ImaxS,0:ND(ng)+1
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
# define OMP !
#endif

/*
** Set tile variable for distributed- or shared-memory configurations.
*/

#ifdef DISTRIBUTE
# define TILE MyRank
#else
# define TILE tile
#endif

/*
** The following definitions contain fortran logical expressions
** equivalent to the question: ''Am I the thread working on a tile
** which is adjacent to the WESTERN, EASTERN, SOUTHERN, or NORTHERN
** edges of the model domain?'' These logical expressions are used to
** update domain boundaries and corners.
*/

#ifdef REFINED_GRID
# define WESTERN_EDGE_REF (Istr.eq.1).and.(ng.eq.1)
# define EASTERN_EDGE_REF (Iend.eq.Lm(ng)).and.(ng.eq.1)
# define SOUTHERN_EDGE_REF (Jstr.eq.1).and.(ng.eq.1)
# define NORTHERN_EDGE_REF (Jend.eq.Mm(ng)).and.(ng.eq.1)
# define SOUTH_WEST_CORNER_REF (Istr.eq.1).and.(Jstr.eq.1).and.(ng.eq.1)
# define NORTH_WEST_CORNER_REF (Istr.eq.1).and.(Jend.eq.Mm(ng)).and.(ng.eq.1)
# define SOUTH_EAST_CORNER_REF (Iend.eq.Lm(ng)).and.(Jstr.eq.1).and.(ng.eq.1)
# define NORTH_EAST_CORNER_REF (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng)).and.(ng.eq.1)
#else
# define WESTERN_EDGE_REF Istr.eq.1
# define EASTERN_EDGE_REF Iend.eq.Lm(ng)
# define SOUTHERN_EDGE_REF Jstr.eq.1
# define NORTHERN_EDGE_REF Jend.eq.Mm(ng)
# define SOUTH_WEST_CORNER_REF (Istr.eq.1).and.(Jstr.eq.1)
# define NORTH_WEST_CORNER_REF (Istr.eq.1).and.(Jend.eq.Mm(ng))
# define SOUTH_EAST_CORNER_REF (Iend.eq.Lm(ng)).and.(Jstr.eq.1)
# define NORTH_EAST_CORNER_REF (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))
#endif
# define WESTERN_EDGE Istr.eq.1
# define EASTERN_EDGE Iend.eq.Lm(ng)
# define SOUTHERN_EDGE Jstr.eq.1
# define NORTHERN_EDGE Jend.eq.Mm(ng)
#define SOUTH_WEST_CORNER (Istr.eq.1).and.(Jstr.eq.1)
#define NORTH_WEST_CORNER (Istr.eq.1).and.(Jend.eq.Mm(ng))
#define SOUTH_EAST_CORNER (Iend.eq.Lm(ng)).and.(Jstr.eq.1)
#define NORTH_EAST_CORNER (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))

/*
** The following definitions are fortran logical expressions use
** to update global variables while avoiding mutual overlap between
** threads in shared-memory configurations.
*/

#ifdef DISTRIBUTE
# define SOUTH_WEST_TEST .TRUE.
# define NORTH_WEST_TEST .TRUE.
# define SOUTH_EAST_TEST .TRUE.
# define NORTH_EAST_TEST .TRUE.
#else
# define SOUTH_WEST_TEST (Istr.eq.1).and.(Jstr.eq.1)
# define NORTH_WEST_TEST (Istr.eq.1).and.(Jend.eq.Mm(ng))
# define SOUTH_EAST_TEST (Iend.eq.Lm(ng)).and.(Jstr.eq.1)
# define NORTH_EAST_TEST (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))
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

#ifdef ICE_MODEL
# define IOUT linew(ng)
# define IUOUT liunw(ng)
# define IEOUT lienw(ng)
#endif

/*
** Set 4DVAR sensitivity switch.
*/

#if defined W4DPSAS_SENSITIVITY || \
    defined W4DVAR_SENSITIVITY
# define SENSITIVITY_4DVAR
#endif

/*
** Set tangent, tl_ioms and adjoint switches.
*/

#if defined CONVOLUTION         || defined CORRELATION        || \
    defined FT_EIGENMODES       || defined FORCING_SV         || \
    defined INNER_PRODUCT       || defined IS4DVAR            || \
    defined IS4DVAR_SENSITIVITY || defined OPT_PERTURBATION   || \
    defined OPT_OBSERVATIONS    || defined PICARD_TEST        || \
    defined R_SYMMETRY          || defined RPM_DRIVER         || \
    defined SANITY_CHECK        || defined SENSITIVITY_4DVAR  || \
    defined TLM_CHECK           || defined TLM_DRIVER         || \
    defined TL_W4DPSAS          || defined TL_W4DVAR          || \
    defined W4DPSAS             || defined W4DVAR
# define TANGENT
#endif
#if defined AD_SENSITIVITY      || defined ADM_DRIVER         || \
    defined AFT_EIGENMODES      || defined CONVOLUTION        || \
    defined CORRELATION         || defined FORCING_SV         || \
    defined INNER_PRODUCT       || defined IS4DVAR            || \
    defined IS4DVAR_SENSITIVITY || defined OPT_PERTURBATION   || \
    defined OPT_OBSERVATIONS    || defined R_SYMMETRY         || \
    defined SANITY_CHECK        || defined SENSITIVITY_4DVAR  || \
    defined SO_SEMI             || defined TLM_CHECK          || \
    defined TL_W4DPSAS          || defined TL_W4DVAR          || \
    defined W4DPSAS             || defined W4DVAR
# define ADJOINT
#endif
#if defined PICARD_TEST        || defined RPM_DRIVER         || \
    defined TL_W4DVAR          || defined W4DVAR             || \
    defined W4DVAR_SENSITIVITY
# define TL_IOMS
#endif
#if !defined ANA_PERTURB                                 && \
    (defined CORRELATION     || defined SANITY_CHECK     || \
     defined R_SYMMETRY)
# define ANA_PERTURB
#endif

/*
** Since some of the tracer advection alogorithms are highly nonlinear,
** it is possible to choose a simpler (less nonlinear) horizontal and
** vertical tracer advection option for the tangent linear, representer
** and adjoint routines. This is likely to improve the convergence of
** the 4DVar algorithms. Notice that this strategy still allows us to
** use highly nonlinear tracer advection schemes in the basic state
** when running the nonlinear model.
*/

#if defined TANGENT || defined TL_IOMS || defined ADJOINT
# if !defined TS_A4HADVECTION_TL       && \
     !defined TS_C2HADVECTION_TL       && \
     !defined TS_C4HADVECTION_TL       && \
     !defined TS_U3HADVECTION_TL
#  if defined TS_A4HADVECTION
#   define TS_A4HADVECTION_TL
#  elif defined TS_C2HADVECTION
#   define TS_C2HADVECTION_TL
#  elif defined TS_C4HADVECTION
#   define TS_C4HADVECTION_TL
#  elif defined TS_U3HADVECTION
#   define TS_U3HADVECTION_TL
#  endif
# endif

# if !defined TS_A4VADVECTION_TL       && \
     !defined TS_C2VADVECTION_TL       && \
     !defined TS_C4VADVECTION_TL       && \
     !defined TS_SVADVECTION_TL
#  if defined TS_A4VADVECTION
#   define TS_A4VADVECTION_TL
#  elif defined TS_C2VADVECTION
#   define TS_C2VADVECTION_TL
#  elif defined TS_C4VADVECTION
#   define TS_C4VADVECTION_TL
#  elif defined TS_SVADVECTION
#   define TS_SVADVECTION_TL
#  endif
# endif
#endif

/*
** Now, we need a switch the tracer advection schemes that
** have not adjointed yet.
*/

#if defined TANGENT || defined TL_IOMS || defined ADJOINT
# if defined TS_A4HADVECTION_TL || defined TS_C2HADVECTION_TL || \
     defined TS_C4HADVECTION_TL || defined TS_U3HADVECTION_TL
#  define TS_HADVECTION_TL
# endif

# if defined TS_A4VADVECTION_TL || defined TS_C2VADVECTION_TL || \
     defined TS_C4VADVECTION_TL || defined TS_SVADVECTION_TL
#  define TS_VADVECTION_TL
# endif
#endif

/*
** Turn off nonlinear model switch.
*/

#define NONLINEAR
#if defined AD_SENSITIVITY   || defined ADM_DRIVER       || \
    defined AFT_EIGENMODES   || defined FT_EIGENMODES    || \
    defined INNER_PRODUCT    || defined OPT_OBSERVATIONS || \
    defined OPT_PERTURBATION || defined PICARD_TEST      || \
    defined RPM_DRIVER       || defined SANITY_CHECK     || \
    defined SO_SEMI          || defined TLM_DRIVER
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
# ifdef ICE_MODEL
#  define IOUT linew(ng)
#  define IUOUT liunw(ng)
#  define IEOUT lienw(ng)
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
#  define KOUTP kstp(ngp)
#  define NOUT nrhs(ng)
#  define NOUTP nnew(ngp)
#  define KOUTC kstp(ngc)
#  define NOUTC nnew(ngc)
# endif
#else
# if defined TANGENT || defined TL_IOMS
#  define TKOUT kstp(ng)
# endif
# define KOUT knew(ng)
# define KOUTP knew(ngp)
# define KOUTC knew(ngc)
#endif

/*
** Set internal switch for the need of a propagator driver.
*/

#if defined AFT_EIGENMODES   || defined ENSEMBLE       || \
    defined FORCING_SV       || defined FT_EIGENMODES  || \
    defined OPT_PERTURBATION || defined PSEUDOSPECTRA  || \
    defined SO_SEMI          || defined SO_TRACE       || \
    defined STOCHASTIC_OPT
# define PROPAGATOR
#endif

/*
** Activate checkpointing switch for GST analysis.  This requires
** a modified ARPACK library for symmetric (*saupd, *seupd) and
** non-symmetric (*naupd, *neupd) drivers.
*/

#ifdef PROPAGATOR
# define CHECKPOINTING
#endif

/*
** Activate processing of forward vertical mixing.
**
*/

#if !defined FORWARD_MIXING  && \
    (defined TANGENT         || defined TL_IOMS    || \
     defined ADJOINT)        && \
    (defined LMD_MIXING      || defined GLS_MIXING || \
     defined MY25_MIXING)
# define FORWARD_MIXING
#endif

/*
** Set internal switches for all the 4DVAR schemes.
*/

#if !defined WEAK_CONSTRAINT     && \
    (defined CONVOLUTION         || defined R_SYMMETRY         || \
     defined TL_W4DPSAS          || defined TL_W4DVAR          || \
     defined W4DPSAS             || defined W4DVAR             || \
     defined W4DPSAS_SENSITIVITY || defined W4DVAR_SENSITIVITY)
# define WEAK_CONSTRAINT
#endif
#if !defined WEAK_CONSTRAINT     && defined RPM_RELAXATION
# undef RPM_RELAXATION
#endif
#if defined CONVOLUTION          || defined CORRELATION         || \
    defined IS4DVAR              || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS     || defined TLM_CHECK           || \
    defined WEAK_CONSTRAINT
# define FOUR_DVAR
#endif
#if !defined WEAK_CONSTRAINT && defined FOUR_DVAR
# define CONVOLVE
#endif
#if defined IS4DVAR
# define BACKGROUND
#endif
#if !(defined W4DPSAS || defined W4DVAR) && defined POSTERIOR_EOFS
# undef POSTERIOR_EOFS
#endif
#if !(defined W4DPSAS || defined W4DVAR) && defined POSTERIOR_ERROR_F
# undef POSTERIOR_ERROR_F
#endif
#if !(defined W4DPSAS || defined W4DVAR) && defined POSTERIOR_ERROR_I
# undef POSTERIOR_ERROR_I
#endif
#if !(defined WEAK_CONSTRAINT || defined IS4DVAR_SENSITIVITY) && \
      defined OBS_IMPACT
# undef OBS_IMPACT
#endif

/*
** Activate internal switch to process 4DVAR observations.
*/

#if defined IS4DVAR            || defined IS4DVAR_SENSITIVITY || \
    defined SENSITIVITY_4DVAR  || defined TLM_CHECK           || \
    defined TL_W4DPSAS         || defined TL_W4DVAR           || \
    defined VERIFICATION       || defined W4DPSAS             || \
    defined W4DVAR
# define OBSERVATIONS
#endif

#if defined IS4DVAR            || defined IS4DVAR_SENSITIVITY || \
    defined R_SYMMETRY         || defined SENSITIVITY_4DVAR   || \
    defined TLM_CHECK          || defined TL_W4DPSAS          || \
    defined TL_W4DVAR          || defined W4DPSAS             || \
    defined W4DVAR
# define TLM_OBS
#endif

/*
** Activate reading and writting of the basic sate.
*/

#if !defined FORWARD_READ      && \
    (defined IS4DVAR           || defined IS4DVAR_SENSITIVITY || \
     defined SENSITIVITY_4DVAR || defined TL_W4DPSAS          || \
     defined TL_W4DVAR         || defined W4DPSAS             || \
     defined W4DVAR)
# define FORWARD_READ
#endif
#if !defined FORWARD_WRITE     && \
    (defined IS4DVAR           || defined IS4DVAR_SENSITIVITY || \
     defined SENSITIVITY_4DVAR || defined TL_W4DPSAS          || \
     defined TL_W4DVAR         || defined W4DPSAS             || \
     defined W4DVAR)
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

#if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES || \
    defined TANGENT  || defined TL_IOMS  || defined ADJOINT
# define EOS_TDERIVATIVE
#endif

/*
** If splines, deactivate horizontal and vertical smoothing of
** Richardson number horizontally and/or vertically.
*/

#ifdef SPLINES
# if defined LMD_MIXING
#  undef RI_HORAVG
#  undef RI_VERAVG
# endif
#endif

/*
** Activate internal switch for the computation of the Brunt-Vaisala
** frequency.
*/

#if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_SKPP    || \
    defined LMD_BKPP   || defined GLS_MIXING  || defined MY25_MIXING
# define BV_FREQUENCY
#endif

/*
** Activate switch for processing climatology data.
*/

#if (defined  ZCLIMATOLOGY && !defined ANA_SSH)     || \
    (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
    (defined  TCLIMATOLOGY && !defined ANA_TCLIMA)  || \
    (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA) || \
    (defined CLIMA_TS_MIX  && defined SOLVE3D       && \
     (defined TS_DIF2      || defined TS_DIF4))
# define CLM_FILE
#endif
#if defined ZCLIMATOLOGY   || defined M2CLIMATOLOGY || \
    defined TCLIMATOLOGY   || defined M3CLIMATOLOGY || \
    defined ZCLM_NUDGING   || defined M2CLM_NUDGING || \
    defined TCLM_NUDGING   || defined M3CLM_NUDGING || \
    (defined CLIMA_TS_MIX  && defined SOLVE3D       && \
     (defined TS_DIF2      || defined TS_DIF4))
# define CLIMATOLOGY
#endif

/*
** Activate internal switch for bottom boundary layer closure.
*/

#if defined SSW_BBL || defined MB_BBL || defined SG_BBL
# define BBL_MODEL
#endif

/*
** Activate internal switch to set-up nudging coefficients.
*/

#if defined ZCLM_NUDGING    || defined M2CLM_NUDGING   || \
    defined TCLM_NUDGING    || defined M3CLM_NUDGING   || \
    defined WEST_FSNUDGING  || defined EAST_FSNUDGING  || \
    defined SOUTH_FSNUDGING || defined NORTH_FSNUDGING || \
    defined WEST_M2NUDGING  || defined EAST_M2NUDGING  || \
    defined SOUTH_M2NUDGING || defined NORTH_M2NUDGING || \
    defined WEST_TNUDGING   || defined EAST_TNUDGING   || \
    defined SOUTH_TNUDGING  || defined NORTH_TNUDGING  || \
    defined WEST_M3NUDGING  || defined EAST_M3NUDGING  || \
    defined SOUTH_M3NUDGING || defined NORTH_M3NUDGING
# define NUDGING_COFF
#endif

/*
** Internal switches to deactivate calling boundary conditions
** during initialization of 2D state variables. Basically,
** we need to apply only non-radiation type boundary conditions.
*/

#if defined WEST_M2RADIATION   || defined WEST_M2FLATHER  || \
    defined EAST_M2RADIATION   || defined EAST_M2FLATHER  || \
    defined SOUTH_M2RADIATION  || defined SOUTH_M2FLATHER || \
    defined NORTH_M2RADIATION  || defined NORTH_M2FLATHER
# define OBC_M2RADIATION
#endif

#if defined WEST_FSRADIATION   || defined WEST_FSCHAPMAN  || \
    defined EAST_FSRADIATION   || defined EAST_FSCHAPMAN  || \
    defined SOUTH_FSRADIATION  || defined SOUTH_FSCHAPMAN || \
    defined NORTH_FSRADIATION  || defined NORTH_FSCHAPMAN
# define OBC_FSRADIATION
#endif

/*
** Activate internal switches requiring open boundary data.
*/

#if (defined WEST_M2RADIATION  && defined WEST_M2NUDGING)  || \
     defined WEST_M2FLATHER    || defined WEST_M2CLAMPED
# define WEST_M2OBC
#endif
#if (defined EAST_M2RADIATION  && defined EAST_M2NUDGING)  || \
     defined EAST_M2FLATHER    || defined EAST_M2CLAMPED
# define EAST_M2OBC
#endif
#if (defined SOUTH_M2RADIATION && defined SOUTH_M2NUDGING) || \
     defined SOUTH_M2FLATHER   || defined SOUTH_M2CLAMPED
# define SOUTH_M2OBC
#endif
#if (defined NORTH_M2RADIATION && defined NORTH_M2NUDGING) || \
     defined NORTH_M2FLATHER   || defined NORTH_M2CLAMPED
# define NORTH_M2OBC
#endif

#if (defined WEST_FSRADIATION  && defined WEST_FSNUDGING)  || \
    (defined WEST_M2REDUCED    && defined FSOBC_REDUCED)   || \
     defined WEST_M2FLATHER    || defined WEST_FSCLAMPED
# define WEST_FSOBC
#endif
#if (defined EAST_FSRADIATION  && defined EAST_FSNUDGING)  || \
    (defined EAST_M2REDUCED    && defined FSOBC_REDUCED)   || \
     defined EAST_M2FLATHER    || defined EAST_FSCLAMPED
# define EAST_FSOBC
#endif
#if (defined SOUTH_FSRADIATION && defined SOUTH_FSNUDGING) || \
    (defined SOUTH_M2REDUCED   && defined FSOBC_REDUCED)   || \
     defined SOUTH_M2FLATHER   || defined SOUTH_FSCLAMPED
# define SOUTH_FSOBC
#endif
#if (defined NORTH_FSRADIATION && defined NORTH_FSNUDGING) || \
    (defined NORTH_M2REDUCED   && defined FSOBC_REDUCED)   || \
     defined NORTH_M2FLATHER   || defined NORTH_FSCLAMPED
# define NORTH_FSOBC
#endif

#if defined FSOBC_REDUCED && \
  !(defined WEST_M2REDUCED  || defined EAST_M2REDUCED  || \
    defined NORTH_M2REDUCED || defined SOUTH_M2REDUCED || \
    defined WEST_M2FLATHER  || defined EAST_M2FLATHER  || \
    defined NORTH_M2FLATHER || defined SOUTH_M2FLATHER)
# undef FSOBC_REDUCED
#endif

#if (defined WEST_M3RADIATION  && defined WEST_M3NUDGING)  || \
     defined WEST_M3CLAMPED
# define WEST_M3OBC
#endif
#if (defined EAST_M3RADIATION  && defined EAST_M3NUDGING)  || \
     defined EAST_M3CLAMPED
# define EAST_M3OBC
#endif
#if (defined SOUTH_M3RADIATION && defined SOUTH_M3NUDGING) || \
     defined SOUTH_M3CLAMPED
# define SOUTH_M3OBC
#endif
#if (defined NORTH_M3RADIATION && defined NORTH_M3NUDGING) || \
     defined NORTH_M3CLAMPED
# define NORTH_M3OBC
#endif

#if (defined WEST_TRADIATION   && defined WEST_TNUDGING)   || \
     defined WEST_TCLAMPED
# define WEST_TOBC
#endif
#if (defined EAST_TRADIATION   && defined EAST_TNUDGING)   || \
     defined EAST_TCLAMPED
# define EAST_TOBC
#endif
#if (defined SOUTH_TRADIATION  && defined SOUTH_TNUDGING)  || \
     defined SOUTH_TCLAMPED
# define SOUTH_TOBC
#endif
#if (defined NORTH_TRADIATION  && defined NORTH_TNUDGING)  || \
     defined NORTH_TCLAMPED
# define NORTH_TOBC
#endif

#ifdef ICE_MODEL
#if (defined WEST_AIRADIATION   && defined WEST_AINUDGING)   || \
     defined WEST_AICLAMPED
# define WEST_AIOBC
#endif
#if (defined EAST_AIRADIATION   && defined EAST_AINUDGING)   || \
     defined EAST_AICLAMPED
# define EAST_AIOBC
#endif
#if (defined SOUTH_AIRADIATION  && defined SOUTH_AINUDGING)  || \
     defined SOUTH_AICLAMPED
# define SOUTH_AIOBC
#endif
#if (defined NORTH_AIRADIATION  && defined NORTH_AINUDGING)  || \
     defined NORTH_AICLAMPED
# define NORTH_AIOBC
#endif

#if (defined WEST_HIRADIATION   && defined WEST_HINUDGING)   || \
     defined WEST_HICLAMPED
# define WEST_HIOBC
#endif
#if (defined EAST_HIRADIATION   && defined EAST_HINUDGING)   || \
     defined EAST_HICLAMPED
# define EAST_HIOBC
#endif
#if (defined SOUTH_HIRADIATION  && defined SOUTH_HINUDGING)  || \
     defined SOUTH_HICLAMPED
# define SOUTH_HIOBC
#endif
#if (defined NORTH_HIRADIATION  && defined NORTH_HINUDGING)  || \
     defined NORTH_HICLAMPED
# define NORTH_HIOBC
#endif

#if (defined WEST_HSNRADIATION   && defined WEST_HSNNUDGING)   || \
     defined WEST_HSNCLAMPED
# define WEST_HSNOBC
#endif
#if (defined EAST_HSNRADIATION   && defined EAST_HSNNUDGING)   || \
     defined EAST_HSNCLAMPED
# define EAST_HSNOBC
#endif
#if (defined SOUTH_HSNRADIATION  && defined SOUTH_HSNNUDGING)  || \
     defined SOUTH_HSNCLAMPED
# define SOUTH_HSNOBC
#endif
#if (defined NORTH_HSNRADIATION  && defined NORTH_HSNNUDGING)  || \
     defined NORTH_HSNCLAMPED
# define NORTH_HSNOBC
#endif

#if (defined WEST_TIRADIATION   && defined WEST_TINUDGING)   || \
     defined WEST_TICLAMPED
# define WEST_TIOBC
#endif
#if (defined EAST_TIRADIATION   && defined EAST_TINUDGING)   || \
     defined EAST_TICLAMPED
# define EAST_TIOBC
#endif
#if (defined SOUTH_TIRADIATION  && defined SOUTH_TINUDGING)  || \
     defined SOUTH_TICLAMPED
# define SOUTH_TIOBC
#endif
#if (defined NORTH_TIRADIATION  && defined NORTH_TINUDGING)  || \
     defined NORTH_TICLAMPED
# define NORTH_TIOBC
#endif

#if (defined WEST_SFWATRADIATION   && defined WEST_SFWATNUDGING) || \
     defined WEST_SFWATCLAMPED
# define WEST_SFWATOBC
#endif
#if (defined EAST_SFWATRADIATION   && defined EAST_SFWATNUDGING) || \
     defined EAST_SFWATCLAMPED
# define EAST_SFWATOBC
#endif
#if (defined SOUTH_SFWATRADIATION  && defined SOUTH_SFWATNUDGING) || \
     defined SOUTH_SFWATCLAMPED
# define SOUTH_SFWATOBC
#endif
#if (defined NORTH_SFWATRADIATION  && defined NORTH_SFWATNUDGING) || \
     defined NORTH_SFWATCLAMPED
# define NORTH_SFWATOBC
#endif

#if (defined WEST_AGEICERADIATION   && defined WEST_AGEICENUDGING) || \
     defined WEST_AGEICECLAMPED
# define WEST_AGEICEOBC
#endif
#if (defined EAST_AGEICERADIATION   && defined EAST_AGEICENUDGING) || \
     defined EAST_AGEICECLAMPED
# define EAST_AGEICEOBC
#endif
#if (defined SOUTH_AGEICERADIATION  && defined SOUTH_AGEICENUDGING) || \
     defined SOUTH_AGEICECLAMPED
# define SOUTH_AGEICEOBC
#endif
#if (defined NORTH_AGEICERADIATION  && defined NORTH_AGEICENUDGING) || \
     defined NORTH_AGEICECLAMPED
# define NORTH_AGEICEOBC
#endif

#if (defined WEST_SIG11RADIATION   && defined WEST_SIG11NUDGING) || \
     defined WEST_SIG11CLAMPED
# define WEST_SIG11OBC
#endif
#if (defined EAST_SIG11RADIATION   && defined EAST_SIG11NUDGING) || \
     defined EAST_SIG11CLAMPED
# define EAST_SIG11OBC
#endif
#if (defined SOUTH_SIG11RADIATION  && defined SOUTH_SIG11NUDGING) || \
     defined SOUTH_SIG11CLAMPED
# define SOUTH_SIG11OBC
#endif
#if (defined NORTH_SIG11RADIATION  && defined NORTH_SIG11NUDGING) || \
     defined NORTH_SIG11CLAMPED
# define NORTH_SIG11OBC
#endif

#if (defined WEST_SIG22RADIATION   && defined WEST_SIG22NUDGING) || \
     defined WEST_SIG22CLAMPED
# define WEST_SIG22OBC
#endif
#if (defined EAST_SIG22RADIATION   && defined EAST_SIG22NUDGING) || \
     defined EAST_SIG22CLAMPED
# define EAST_SIG22OBC
#endif
#if (defined SOUTH_SIG22RADIATION  && defined SOUTH_SIG22NUDGING) || \
     defined SOUTH_SIG22CLAMPED
# define SOUTH_SIG22OBC
#endif
#if (defined NORTH_SIG22RADIATION  && defined NORTH_SIG22NUDGING) || \
     defined NORTH_SIG22CLAMPED
# define NORTH_SIG22OBC
#endif

#if (defined WEST_SIG12RADIATION   && defined WEST_SIG12NUDGING) || \
     defined WEST_SIG12CLAMPED
# define WEST_SIG12OBC
#endif
#if (defined EAST_SIG12RADIATION   && defined EAST_SIG12NUDGING) || \
     defined EAST_SIG12CLAMPED
# define EAST_SIG12OBC
#endif
#if (defined SOUTH_SIG12RADIATION  && defined SOUTH_SIG12NUDGING) || \
     defined SOUTH_SIG12CLAMPED
# define SOUTH_SIG12OBC
#endif
#if (defined NORTH_SIG12RADIATION  && defined NORTH_SIG12NUDGING) || \
     defined NORTH_SIG12CLAMPED
# define NORTH_SIG12OBC
#endif

#if (defined WEST_MIRADIATION   && defined WEST_MINUDGING)   || \
     defined WEST_MICLAMPED
# define WEST_MIOBC
#endif
#if (defined EAST_MIRADIATION   && defined EAST_MINUDGING)   || \
     defined EAST_MICLAMPED
# define EAST_MIOBC
#endif
#if (defined SOUTH_MIRADIATION  && defined SOUTH_MINUDGING)  || \
     defined SOUTH_MICLAMPED
# define SOUTH_MIOBC
#endif
#if (defined NORTH_MIRADIATION  && defined NORTH_MINUDGING)  || \
     defined NORTH_MICLAMPED
# define NORTH_MIOBC
#endif
#endif

#ifdef SOLVE3D
# if defined WEST_FSOBC  || defined EAST_FSOBC  || \
     defined SOUTH_FSOBC || defined NORTH_FSOBC || \
     defined WEST_M2OBC  || defined EAST_M2OBC  || \
     defined SOUTH_M2OBC || defined NORTH_M2OBC || \
     defined WEST_M3OBC  || defined EAST_M3OBC  || \
     defined SOUTH_M3OBC || defined NORTH_M3OBC || \
     defined WEST_TOBC   || defined EAST_TOBC   || \
     defined SOUTH_TOBC  || defined NORTH_TOBC
#  define OBC
# endif
#else
# if defined WEST_FSOBC  || defined EAST_FSOBC  || \
     defined SOUTH_FSOBC || defined NORTH_FSOBC || \
     defined WEST_M2OBC  || defined EAST_M2OBC  || \
     defined SOUTH_M2OBC || defined NORTH_M2OBC
#  define OBC
# endif
#endif

/*
** Define internal flag indicating processing of input boundary
** NetCDF file.
*/

#if (!defined ANA_FSOBC && \
     (defined WEST_FSOBC  || defined EAST_FSOBC    || \
      defined SOUTH_FSOBC || defined NORTH_FSOBC)) || \
    (!defined ANA_M2OBC && \
     (defined WEST_M2OBC  || defined EAST_M2OBC    || \
      defined SOUTH_M2OBC || defined NORTH_M2OBC)) || \
    (!defined ANA_M3OBC && \
     (defined WEST_M3OBC  || defined EAST_M3OBC    || \
      defined SOUTH_M3OBC || defined NORTH_M3OBC)) || \
    (!defined ANA_TOBC && \
     (defined WEST_TOBC   || defined EAST_TOBC    || \
      defined SOUTH_TOBC  || defined NORTH_TOBC))
# define OBC_DATA
#endif

/*
** Activate internal switches for volume conservation at open boundary.
*/

#if defined WEST_VOLCONS  || defined EAST_VOLCONS  || \
    defined NORTH_VOLCONS || defined SOUTH_VOLCONS
# define OBC_VOLCONS
#endif

/*
** Activate assimilation switches.
*/

#if defined ASSIMILATION_SSH || defined ASSIMILATION_SST   || \
    defined ASSIMILATION_T   || defined ASSIMILATION_UVsur || \
    defined ASSIMILATION_UV
# define ASSIMILATION
#endif
#if defined NUDGING_SST   || defined NUDGING_T   || \
    defined NUDGING_UVsur || defined NUDGING_UV
# define NUDGING
#endif

/*
** Check if it is meaningful to write out time-averaged vertical
** mixing coefficients.
*/

#if !defined LMD_MIXING && !defined MY25_MIXING && !defined GLS_MIXING
# if defined AVERAGES
#  if defined AVERAGES_AKV
#    undef AVERAGES_AKV
#  endif
#  if defined AVERAGES_AKT
#    undef AVERAGES_AKT
#  endif
#  if defined AVERAGES_AKS && !defined SALINITY
#    undef AVERAGES_AKS
#  endif
# endif
#endif
#if defined AVERAGES_WEC && (!defined WEC_MELLOR || \
                             !defined WEC_VF)
# undef AVERAGES_WEC
#endif

/*
** Activate internal biology option when using any type of biological
** module.
*/

#if defined BIO_FENNEL  || defined ECOSIM      || \
    defined NEMURO      || defined NPZD_FRANKS || \
    defined NPZD_IRON   || defined NPZD_POWELL
# define BIOLOGY
#endif

/*
** Define internal option to couple to other models.
**
*/

#if defined ROMS_MODEL && (defined SWAN_MODEL || defined WRF_MODEL)
# define ROMS_COUPLING
#endif
#if defined SWAN_MODEL && (defined ROMS_MODEL || defined WRF_MODEL)
# define SWAN_COUPLING
#endif
#if defined WRF_MODEL && (defined SWAN_MODEL || defined ROMS_MODEL)
# define WRF_COUPLING
#endif

#if defined WRF_COUPLING && defined ROMS_COUPLING
# define AIR_OCEAN
#endif

#if defined WRF_COUPLING && defined SWAN_COUPLING
# define AIR_WAVES
#endif

#if (defined REFDIF_COUPLING || defined SWAN_COUPLING) && \
     defined ROMS_COUPLING
# define WAVES_OCEAN
#endif

#if defined AIR_OCEAN || defined AIR_WAVES || defined WAVES_OCEAN
# define MODEL_COUPLING
#endif

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
    defined WAVES_OCEAN
# define WAVES_DIR
#endif

#if  defined BBL_MODEL   && \
   !(defined SSW_CALC_UB || defined MB_CALC_UB ||  \
     defined SG_CALC_UB)
# define WAVES_UB
#endif

#if (defined BBL_MODEL        && !defined WAVES_UB) ||  \
     defined WEC              || \
     defined ZOS_HSIG         || defined COARE_TAYLOR_YELLAND || \
     defined BEDLOAD_SOULSBY  || defined WAVES_OCEAN || \
     defined DRENNAN
# define WAVES_HEIGHT
#endif

#if defined WEC || defined BEDLOAD_SOULSBY || \
    defined WAVES_OCEAN
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

#if defined BBL_MODEL || defined WAVES_OCEAN
# define WAVES_BOT_PERIOD
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

#if defined BEDLOAD_MPM || defined BEDLOAD_SOULSBY
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
    (!defined AIR_OCEAN    && \
     !defined BULK_FLUXES  && !defined ANA_SMFLUX)   || \
    (!defined AIR_OCEAN    && \
     !defined BULK_FLUXES  && !defined ANA_STFLUX)   || \
    ( defined SALINITY     && !defined ANA_SSFLUX)   || \
    ( defined BULK_FLUXES  && !defined LONGWAVE  && !defined AIR_OCEAN)      || \
    ( defined BULK_FLUXES  && (!defined ANA_PAIR && !defined AIR_OCEAN))     || \
    ( defined BULK_FLUXES  && (!defined ANA_TAIR && !defined AIR_OCEAN))     || \
    ( defined BULK_FLUXES  && (!defined ANA_HUMIDITY && !defined AIR_OCEAN)) || \
    ( defined BULK_FLUXES  && (!defined ANA_CLOUD && !defined AIR_OCEAN))    || \
    ( defined BULK_FLUXES  && (!defined ANA_RAIN && !defined AIR_OCEAN))     || \
    ( defined BULK_FLUXES  && (!defined ANA_WINDS && !defined AIR_OCEAN))    || \
    ( defined BULK_FLUXES  && (!defined ANA_SRFLUX && !defined AIR_OCEAN))   || \
    ( defined LMD_SKPP     && !defined ANA_SRFLUX)   || \
    ( defined SOLAR_SOURCE && !defined ANA_SRFLUX)   || \
    ( defined BBL_MODEL    && (!defined ANA_WWAVE    && \
     !defined WAVES_OCEAN))                          || \
    ( defined BIOLOGY      && !defined ANA_SPFLUX)   || \
    ( defined BIOLOGY      && !defined ANA_BPFLUX)   || \
    ( defined SEDIMENT     && !defined ANA_SPFLUX)   || \
    ( defined SEDIMENT     && !defined ANA_BPFLUX)   || \
    ( defined WAVE_DATA    && (!defined ANA_WWAVE    && \
     !defined WAVES_OCEAN  && !defined INWAVE_MODEL))
#  define FRC_FILE
# endif
#else
# if !defined ANA_SMFLUX
#  define FRC_FILE
# endif
#endif
#ifdef ANA_NCEP
# undef FRC_FILE
#endif

/*
** Check if processing timeless data.
*/

#if (!defined ANA_PSOURCE  && \
     (defined UV_PSOURCE   || defined TS_PSOURCE || \
      defined Q_PSOURCE))  || \
    (defined  SSH_TIDES    || defined UV_TIDES)
# define TIMELESS_DATA
#endif

/*
** Check analytical initial conditions options.
*/

#if defined ANA_BIOLOGY && !defined BIOLOGY
# undef ANA_BIOLOGY
#endif
#if defined ANA_PASSIVE && !defined T_PASSIVE
# undef ANA_PASSIVE
#endif
#if defined ANA_SEDIMENT && !(defined SEDIMENT || defined BBL_MODEL)
# undef ANA_SEDIMENT
#endif
#if  !defined ANA_INITIAL || \
    ( defined BIOLOGY     && !defined ANA_BIOLOGY)  || \
    ( defined T_PASSIVE   && !defined ANA_PASSIVE)  || \
    ( defined SEDIMENT    && !defined ANA_SEDIMENT) || \
    ( defined BBL_MODEL   && !defined ANA_SEDIMENT)
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
#if defined SHORTWAVE && defined NCEP_FLUXES
# undef SHORTWAVE
#endif
#if !defined SHORTWAVE   && defined ANA_SRFLUX
# undef ANA_SRFLUX
#endif
#if !defined SHORTWAVE   && defined DIURNAL_SRFLUX
# undef DIURNAL_SRFLUX
#endif

/*
** Define internal clouds option.  Undefine analytical
** shortwave option if not needed.
*/

#if (defined BULK_FLUXES && defined LONGWAVE) || defined ECOSIM || \
    (defined ANA_SRFLUX  && defined ALBEDO)
# define CLOUDS
#endif
#if !defined CLOUDS && defined ANA_CLOUD
# undef ANA_CLOUD
#endif

/*
** Check if it is meaningful to write out momentum/tracer diagnostics
** and activate internal diagnostics option.
*/

#if !defined SOLVE3D || defined TS_FIXED
# if defined DIAGNOSTICS_TS
#   undef DIAGNOSTICS_TS
# endif
#endif
#if !defined BIO_FENNEL && defined DIAGNOSTICS_BIO
#  undef DIAGNOSTICS_BIO
#endif
#if defined DIAGNOSTICS_BIO || defined DIAGNOSTICS_TS || \
    defined DIAGNOSTICS_UV
# define DIAGNOSTICS
#endif

/*
** Activate switch to modify MAIN3D to recompute depths and
** thicknesses using the new time filtered free-surface.  This
** call is moved from STEP2D to facilitate nesting.
** This strategy needs to be tested in the TLM, RPM, and ADM.
*/

#if !(defined ADJOINT || defined TANGENT || defined TL_IOMS)
# define MOVE_SET_DEPTH
#endif

/*
** Check if any analytical expression is defined.
*/

#if defined ANA_BIOLOGY    || defined ANA_BPFLUX     || \
    defined ANA_BSFLUX     || defined ANA_BTFLUX     || \
    defined ANA_CLOUD      || defined ANA_DIAG       || \
    defined ANA_FSOBC      || defined ANA_GRID       || \
    defined ANA_HUMIDITY   || defined ANA_INITIAL    || \
    defined ANA_M2CLIMA    || defined ANA_M2OBC      || \
    defined ANA_M3CLIMA    || defined ANA_M3OBC      || \
    defined ANA_MASK       || defined ANA_PAIR       || \
    defined ANA_PASSIVE    || defined ANA_PERTURB    || \
    defined ANA_PSOURCE    || defined ANA_RAIN       || \
    defined ANA_SEDIMENT   || defined ANA_SMFLUX     || \
    defined ANA_SPFLUX     || defined ANA_SPINNING   || \
    defined ANA_SRFLUX     || defined ANA_SSFLUX     || \
    defined ANA_SSH        || defined ANA_SSS        || \
    defined ANA_SST        || defined ANA_STFLUX     || \
    defined ANA_TAIR       || defined ANA_TCLIMA     || \
    defined ANA_TOBC       || defined ANA_VMIX       || \
    defined ANA_WINDS      || defined ANA_WWAVE      || \
    defined DIFF_GRID      || defined NUDGING_COFF   || \
    defined SPONGE         || defined VISC_GRID
# define ANALYTICAL
#endif

/*
** If splitting 3rd-order upstream bias horizontal advection of
** tracer, activate other needed flags.
*/

#ifdef TS_U3ADV_SPLIT
# define DIFF_3DCOEF
# ifdef TS_U3HADVECTION
#  undef TS_U3HADVECTION
# endif
# ifndef TS_C4HADVECTION
#  define TS_C4HADVECTION
# endif
# ifndef TS_C4VADVECTION
#  define TS_C4VADVECTION
# endif
# ifndef TS_DIF4
#  define TS_DIF4
# endif
# ifdef TS_DIF2
#  undef TS_DIF2
# endif
# ifdef TS_SMAGORINSKY
#  undef TS_SMAGORINSKY
# endif
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
