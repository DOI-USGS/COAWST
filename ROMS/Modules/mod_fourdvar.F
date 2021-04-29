#include "cppdefs.h"
      MODULE mod_fourdvar

#if defined FOUR_DVAR || defined VERIFICATION
!
!svn $Id: mod_fourdvar.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Variational data assimilation variables:                            !
!                                                                      !
!  ADmodVal      Adjoint model values at observation locations.        !
!  BackCost      Current background cost function misfit (mean squared !
!                  difference) between model and background state, Jb. !
!  CostFun       Current iteration total cost function (background     !
!                  plus observations), J.                              !
!  CostFunOld    Previous iteration total cost function (background    !
!                  plus observations), J.                              !
!  CostNorm      Total cost function normalization scales              !
!                  (minimization starting value, Jb+Jo).               !
!  Cost0         The total cost function at the beggining of the inner !
!                  loop (inner=0) for each outer loop.                 !
!  DTsizeH       Horizontal diffusion time-step size for spatial       !
!                  convolutions.                                       !
!  DTsizeV       Vertical diffusion time-step size for spatial         !
!                  convolutions.                                       !
!  GradErr       Upper bound on relatice error of the gradient.        !
!  HevecErr      Maximum error bound on Hessian eigenvectors.          !
!  KhMax         Maximum  horizontal diffusion coefficient.            !
!  KhMin         Minimum  horizontal diffusion coefficient.            !
!  KvMax         Maximum  vertical diffusion coefficient.              !
!  KvMin         Minimum  vertical diffusion coefficient.              !
!  Load_Zobs     Logical switch indicating that input Zobs is negative !
!                  and fractional vertical positions are computed in   !
!                  "extract_obs3d" and written to observation NetCDF   !
!                  file for latter use.                                !
!  LhessianEV    Logical switch to compute Hessian eigenvectors.       !
!  LhotStart     Logical switch to activate hot start of subsquent     !
!                  outer loops in the weak-constraint algorithms.      !
!  Lprecond      Logical switch to activate conjugate gradient         !
!                  preconditioning.                                    !
!  Lritz         Logical switch to activw Ritz limited-memory          !
!                  preconditioning.                                    !
# ifdef RPCG
!  LaugWeak      Logical switch for computing augmented model error    !
!                terms in the routine "error_covariance".              !
# endif
!  nConvRitz     Number of converged Ritz eigenvalues.                 !
!  Nimpact       Outer loop to consider for observations impact and    !
!                  observations sensitivity.                           !
!  NLmodVal      Nonlinear model values at observation locations.      !
!  NHsteps       Full number of horizontal diffusion steps for spatial !
!                  convolutions.                                       !
!  NLobsCost     Current NLM observation cost function misfit (mean    !
!                  squared difference) between NLM and observations.   !
!  NpostI        Number of Lanczos iterations used for the posterior   !
!                  analysis error covariance matrix estimation.        !
!  NritzEV       If preconditioning, number of eigenpairs to use.      !
!  NVsteps       Full number of vertical diffusion steps for spatial   !
!                  convolutions.                                       !
!  ObsAngler     Observation vector rotation angle at RHO-points.      !
!  ObsCost       Current observation cost function misfit (mean        !
!                  squared difference) between model and observations. !
!  ObsCount      Current observation counts per state variable.        !
!  ObsErr        Observation error.                                    !
!  ObsMeta       Observation meta values.                              !
!  ObsName       String describing the observation types.              !
!  ObsProv       Observation provenance flags.                         !
!  ObsReject     Current rejected observation counts per state         !
!                  variable.                                           !
!  ObsScale      Observation screening and quality control scale,      !
!                  0: reject 1: accept.                                !
!  ObsState2Type Mapping indices from state variable to observation    !
!                  type.                                               !
!  ObsType2State Mapping indices from observation type to state        !
!                  variable.                                           !
!  ObsType       Observation type identifier.                          !
!  ObsVal        Observation values.                                   !
!  ObsVetting    Processing flag used to reject (zero) or accept       !
!                  (unity) observations.                               !
!  Optimality    normalized, optimal cost function minimum.            !
!  Ritz          Ritz eigenvalues to compute approximated Hessian.     !
!  RitzMaxErr    Ritz values maximum error limit.                      !
!  TLmodVal      Tangent linear model values at observation locations. !
!  Vdecay        Covariance vertical decorrelation scale (m).          !
!  Tobs          Observations time (days).                             !
!  Xobs          Observations X-locations (grid coordinates).          !
!  Yobs          Observations Y-locations (grid coordinates).          !
!  Zobs          Observations Z-locations (grid coordinates).          !
!  cg_Gnorm      Initial gradient vector normalization factor.         !
!  cg_Greduc     Reduction in the gradient norm; excess cost function. !
!  cg_QG         Lanczos vector normalization factor.                  !
!  cg_Ritz       Eigenvalues of the Lanczos recurrence term Q(k)T(k).  !
!  cg_RitzErr    Eigenvalues relative error.                           !
!  cg_Tmatrix    Lanczos recurrence symmetric, tridiagonal matrix.     !
!  cg_alpha      Conjugate gradient coefficient.                       !
!  cg_beta       Conjugate gradient coefficient.                       !
!  cg_delta      Lanczos algorithm coefficient.                        !
!  cg_gamma      Lanczos algorithm coefficient.                        !
!  cg_tau        Conjugate gradient coefficient.                       !
!  cg_zu         Lanczos tridiagonal matrix, upper diagonal elements.  !
!  cg_zv         Eigenvectors of Lanczos recurrence term Q(k)T(k).     !
!  haveADmod     Logical switch indicating that there is representer   !
!                  coefficients available to process in DAV(ng)%name   !
!                  file.                                               !
!  haveNLmod     Logical switch indicating that there is nonlinear     !
!                  model data available to process in DAV(ng)%name     !
!                  file.                                               !
!  haveTLmod     Logical switch indicating that there is tangent       !
!                  model data available to process in DAV(ng)%name     !
!                  file.                                               !
!  haveObsMeta   Logical switch indicating loading and processing of   !
!                  observations meta field.                            !
!                                                                      !
# ifdef FOUR_DVAR
!  Variables in observation space needed Desroziers et al. (2005)      !
!  statistics at observation locations:                                !
!                                                                      !
!  BgErr         4D-Var Background error standard deviation.           !
!  NLincrement   4D-Var NLM increment, analysis minus background.      !
!  innovation    4D-Var innovation, observations minus background.     !
!  residual      4D-Var residual, observation minus analysis.          !
!                                                                      !
# endif
# ifdef BGQC
!  Background quality control of observations:                         !
!                                                                      !
!  BgThresh      Number of standard deviations threshold used in the   !
!                  background quality control of observations.         !
!  Jact_S        Saved actual cost function.                           !
# endif
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

# ifdef OBSERVATIONS
!
!  Define module structure.
!
        TYPE T_FOURDVAR

          integer , allocatable :: NobsSurvey(:)
          integer , allocatable :: ObsCount(:)
          integer , allocatable :: ObsReject(:)
#  ifdef FOUR_DVAR
          real(r8), allocatable :: BackCost(:)
          real(r8), allocatable :: Cost0(:)
          real(r8), allocatable :: CostFun(:)
          real(r8), allocatable :: CostFunOld(:)
          real(r8), allocatable :: CostNorm(:)
          real(r8), allocatable :: ObsCost(:)
          real(r8), allocatable :: DataPenalty(:)
#   if defined DATALESS_LOOPS    || defined RBL4DVAR  || \
       defined SENSITIVITY_4DVAR || defined SP4DVAR   || \
       defined TL_RBL4DVAR
          real(r8), allocatable :: NLPenalty(:)
#   endif
          real(r8), allocatable :: NLobsCost(:)

#   if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
          real(r8), allocatable :: bc_ak(:)
          real(r8), allocatable :: bc_bk(:)

          real(r8), allocatable :: zdf1(:)
          real(r8), allocatable :: zdf2(:)
          real(r8), allocatable :: zdf3(:)

          real(r8), allocatable :: pc_r2d(:,:)
          real(r8), allocatable :: rhs_r2d(:,:)
          real(r8), allocatable :: r2d_ref(:,:)
          real(r8), allocatable :: zeta_ref(:,:)

          real(r8), allocatable :: p_r2d(:,:,:)
          real(r8), allocatable :: r_r2d(:,:,:)
          real(r8), allocatable :: bp_r2d(:,:,:)
          real(r8), allocatable :: br_r2d(:,:,:)

          real(r8), allocatable :: tl_rhs_r2d(:,:)
          real(r8), allocatable :: tl_r2d_ref(:,:)
          real(r8), allocatable :: tl_zeta_ref(:,:)

          real(r8), allocatable :: ad_rhs_r2d(:,:)
          real(r8), allocatable :: ad_r2d_ref(:,:)
          real(r8), allocatable :: ad_zeta_ref(:,:)
#   endif
#  endif
          real(r8), allocatable :: SurveyTime(:)

#  ifdef RPCG
          real(r8), allocatable :: cg_pxout(:)
          real(r8), allocatable :: cg_pxsave(:)
          real(r8), allocatable :: cg_pxsum(:)
          real(r8), allocatable :: cg_Dpxsum(:)
          real(r8), allocatable :: tl_cg_pxsave(:)
          real(r8), allocatable :: ad_cg_pxsave(:)
#  endif

        END TYPE T_FOURDVAR

        TYPE (T_FOURDVAR), allocatable :: FOURDVAR(:)
!
!  Define other module arrays.
!
        integer,  allocatable :: ObsType(:)
        integer,  allocatable :: ObsState2Type(:)
        integer,  allocatable :: ObsType2State(:)
        integer,  allocatable :: ObsProv(:)

#  ifdef CURVGRID
        real(r8), allocatable :: ObsAngler(:)
#  endif
        real(r8), allocatable :: ObsErr(:)
        real(r8), allocatable :: ObsMeta(:)
        real(r8), allocatable :: ObsScale(:)
#  ifndef WEAK_CONSTRAINT
        real(r8), allocatable :: ObsScaleGlobal(:)
#  endif
        real(r8), allocatable :: ObsVetting(:)
        real(r8), allocatable :: ObsVal(:)
        real(dp), allocatable :: Tobs(:)
        real(r8), allocatable :: Xobs(:)
        real(r8), allocatable :: Yobs(:)
#  ifdef SOLVE3D
        real(r8), allocatable :: Zobs(:)
#  endif
        real(r8), allocatable :: ADmodVal(:)
        real(r8), allocatable :: NLmodVal(:)
#  ifdef SP4DVAR
        real(r8), allocatable :: ADmodVal_S(:,:)
        real(r8), allocatable :: harnoldi(:,:,:)
        real(r8), allocatable :: gmze(:,:)
        real(r8), allocatable :: cg_beta0(:)
        real(r8), allocatable :: Jobs(:,:)
#  endif
#  ifdef TLM_OBS
        real(r8), allocatable :: TLmodVal(:)
#   if defined ARRAY_MODES       || defined CLIPPING  || \
       defined RBL4DVAR          || defined R4DVAR    || \
       defined SENSITIVITY_4DVAR || defined SP4DVAR   || \
       defined TL_RBL4DVAR       || defined TL_R4DVAR
        real(r8), allocatable :: TLmodVal_S(:,:,:)
#    if defined RPCG              && \
       (defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
        defined TL_R4DVAR)
        real(r8), allocatable :: tl_TLmodVal_S(:,:,:)
        real(r8), allocatable :: ad_TLmodVal_S(:,:,:)
#    endif
#   endif
#   if defined ARRAY_MODES       || defined RBL4DVAR    || \
       defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR
        real(r8), allocatable :: BCKmodVal(:)

#    ifdef RPCG
        real(r8), allocatable :: Hbk(:,:)
#     if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
         defined TL_R4DVAR
        real(r8), allocatable :: tl_BCKmodVal(:)
        real(r8), allocatable :: tl_Hbk(:)
        real(r8), allocatable :: ad_BCKmodVal(:)
        real(r8), allocatable :: ad_Hbk(:)
#     endif
#    endif
#   endif
#  endif
#  ifdef WEAK_CONSTRAINT
        real(r8), allocatable :: zcglwk(:,:,:)
#   ifdef RPCG
        real(r8), allocatable :: vcglwk(:,:,:)
        real(r8), allocatable :: Jb0(:)
#   endif
        real(r8), allocatable :: vcglev(:,:,:)
        real(r8), allocatable :: zgrad0(:,:)
        real(r8), allocatable :: vgrad0(:)
        real(r8), allocatable :: vgrad0s(:)
        real(r8), allocatable :: gdgw(:)
        real(r8), allocatable :: vgnorm(:)
        real(r8), allocatable :: cg_innov(:)
        real(r8), allocatable :: cg_dla(:,:)
#   ifndef RPCG
        real(r8), allocatable :: cg_pxsave(:)
#   endif
#   if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
       defined TL_R4DVAR
        real(r8), allocatable :: ad_zcglwk(:,:)
        real(r8), allocatable :: ad_zgrad0(:)
#    if defined RPCG              && \
       (defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
        defined TL_R4DVAR)
        real(r8), allocatable :: ad_vcglwk(:,:)
        real(r8), allocatable :: ad_vgrad0(:)
        real(r8), allocatable :: ad_vgrad0s(:)
        real(r8), allocatable :: ad_Jb0(:)
        real(r8), allocatable :: ad_ObsErr(:)
        real(r8), allocatable :: ad_ObsScale(:)
#    endif
        real(r8), allocatable :: ad_cg_innov(:)
        real(r8), allocatable :: ad_cg_pxsave(:)
        real(r8), allocatable :: ad_cg_pxout(:)
#    ifdef IMPACT_INNER
        real(r8), allocatable :: ad_ObsVal(:,:)
#    else
        real(r8), allocatable :: ad_ObsVal(:)
#    endif

        real(r8), allocatable :: tl_zcglwk(:,:)
        real(r8), allocatable :: tl_zgrad0(:)
#    ifdef RPCG
        real(r8), allocatable :: tl_vcglwk(:,:)
        real(r8), allocatable :: tl_vgrad0(:)
        real(r8), allocatable :: tl_vgrad0s(:)
        real(r8), allocatable :: tl_Jb0(:)
        real(r8), allocatable :: tl_ObsErr(:)
        real(r8), allocatable :: tl_ObsScale(:)
#    endif
        real(r8), allocatable :: tl_cg_innov(:)
        real(r8), allocatable :: tl_cg_pxsave(:)
        real(r8), allocatable :: tl_cg_pxout(:)
        real(r8), allocatable :: tl_ObsVal(:)
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Observations parameters.
!-----------------------------------------------------------------------
!
!  Switches indicating that at input observations vertical position were
!  given in meters (Zobs < 0) so the fractional vertical level position
!  is computed during extraction and written to Observation NetCDF file.
!
        logical, allocatable :: Load_Zobs(:)
        logical, allocatable :: wrote_Zobs(:)
!
!  Maximum number of observations to process.
!
        integer :: Mobs
!
!  Number of model state variables to process.
!
        integer, allocatable :: NstateVar(:)
!
!  Number of observation types and names to process. If the observation
!  operators for a specific application have a 1-to-1 association with
!  the model state variables,  NobsVar = NstateVar.  However, if more
!  that one state variable is needed to evaluate a particular type of
!  observation (HF radials, travel time, pressure, etc.), a new class
!  name is added to the state variable list to assess and differentiate
!  its impact.  Then, NobsVar = NstateVar + NextraObs.
!
        integer, allocatable :: NobsVar(:)

        character(len=40), allocatable :: ObsName(:)
!
! Number of extra-observation classes, extra-observation type indices,
! and extra-observation type names. It is used in observation operators
! that require more than one state variable or not directly associated
! with state variables.
!
        integer :: NextraObs = 0

        integer, allocatable :: ExtraIndex(:)

        character(len=40), allocatable :: ExtraName(:)
!
!  Number of interpolation weights and (I,J,K) indices offsets.
!
# ifdef SOLVE3D
        integer, parameter :: Nweights = 8

        integer, parameter, dimension(Nweights) :: Ioffset =            &
     &                      (/ 0, 1, 0, 1, 0, 1, 0, 1 /)
        integer, parameter, dimension(Nweights) :: Joffset =            &
     &                      (/ 0, 0, 1, 1, 0, 0, 1, 1 /)
        integer, parameter, dimension(Nweights) :: Koffset =            &
     &                      (/ 0, 0, 0, 0, 1, 1, 1, 1 /)
# else
        integer, parameter :: Nweights = 4

        integer, parameter, dimension(Nweights) :: Ioffset =            &
     &                      (/ 0, 1, 0, 1 /)
        integer, parameter, dimension(Nweights) :: Joffset =            &
     &                      (/ 0, 0, 1, 1 /)
# endif
!
!  Size of observation NetCDF file unlimited dimension.
!
        integer, allocatable :: Ndatum(:)
!
!  Number of observations surveys available.
!
        integer, allocatable :: Nsurvey(:)
!
!  Observation surveys counter.
!
        integer, allocatable :: ObsSurvey(:)
!
!  Current number of observations processed.
!
        integer, allocatable :: Nobs(:)
!
!  Current starting and ending observation file index processed.
!
        integer, allocatable :: NstrObs(:)
        integer, allocatable :: NendObs(:)
!
!  Background error covariance normalization method:
!
!       [0] Exact, very expensive
!       [1] Approximated, randomization
!
        integer, allocatable :: Nmethod(:)
!
!  Random number generation scheme for randomization:
!
!       [0] Intrinsic F90 routine "randon_number"
!       [1] Gaussian distributed deviates, numerical recipes
!
        integer, allocatable :: Rscheme(:)
!
!  Number of iterations in the randomization ensemble used to
!  compute the background error covariance, B, normalization
!  factors. This factors ensure that the diagonal elements of
!  B are equal to unity (Fisher and Coutier, 1995).
!
        integer :: Nrandom = 1000
!
!  Number of Lanczos iterations used in the posterior analysis
!  error covariance matrix estimation.
!
        integer :: NpostI
!
!  Outer loop to consider in the computatio of the observations
!  impact or observations sensitivity, Nimpact =< Nouter.  This
!  facilitates the computations with multiple outer loop 4D-Var
!  applications.  The observation analysis needs to be computed
!  separately for each outer loop.  The full analysis for all
!  outer loops are combined offline.
!
        integer :: Nimpact
!
!  Parameter to either process the eigenvector of the stabilized
!  representer matrix to whe computing array modes (Nvct=1 less
!  important eigenvector, Nvct=Ninner most important eigenvector)
!  or cut-off eigenvector for the clipped analysis when (only
!  eigenvector Nvct:Ninner to will be processed, others will be
!  disgarded).
!
        integer :: Nvct
!
!  Optimal, normalized cost funtion minimum.  If the statistical
!  hypotheses between the background and observations errors
!  is consistemt, the cost function value at the minimum, Jmin,
!  is idealy equal to half the number of observations assimilated
!  (Optimality=1=2*Jmin/Nobs), for a linear system.
!
        real(r8), allocatable :: Optimality(:)
!
!  Switch to activate the processing of model state at the observation
!  locations.
!
        logical, allocatable :: ProcessObs(:)
#  ifdef SP4DVAR
!
!  Switch to (de)activate the assimilation of observations during
!  saddle point 4D-Var.
!
        logical, allocatable :: Lsadd(:)
#  endif
!
!  Switch to activate writting of nonlinear model values at
!  observations locations into observations NetCDF file.
!
        logical, allocatable :: wrtNLmod(:)
!
!  Sitch to process and write out observation accept/reject flag used
!  for screening and quality control.
!
        logical, allocatable :: wrtObsScale(:)
!
!  Switch to activate writting of representer model values at
!  observations locations into observations NetCDF file.
!
        logical, allocatable :: wrtRPmod(:)
!
!  Switch to activate writting of tangent linear model values at
!  observations locations into observations NetCDF file.
!
        logical, allocatable :: wrtTLmod(:)
!
!  Switch to activate writting of initial and final model-observation
!  misfit (innovation) vector into 4DVAR output NetCDF file.
!
        logical, allocatable :: wrtMisfit(:)

# if defined I4DVAR_ANA_SENSITIVITY
#  ifdef OBS_IMPACT
!
!  Switch to control the writing of the total observation impact for
!  incremental 4D-Var.
!
        logical, allocatable :: wrtIMPACT_TOT(:)
#  endif
#  ifdef OBS_IMPACT_SPLIT
!
!  Switch to control the writing of the observation impact for initial
!  conditions for incremental 4D-Var.
!
        logical, allocatable :: wrtIMPACT_IC(:)
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Switch to control the writing of the observation impact for surface
!  forcing for incremental 4D-Var.
!
        logical, allocatable :: wrtIMPACT_FC(:)
#   endif
#   if defined ADJUST_BOUNDARY
!
!  Switch to control the writing of the observation impact for surface
!  forcing for incremental 4D-Var.
!
        logical, allocatable :: wrtIMPACT_BC(:)
#   endif
#  endif
# endif
!
!  Swiches indicating that there is representer coeffiecients,
!  nonlinear, and tangent linear model data available to process
!  in 4DVAR NetCDF output file (DAV(ng)%name).  At the beeginning,
!  there is not data since this file has been just defined.
!
        logical, allocatable :: haveADmod(:)
        logical, allocatable :: haveNLmod(:)
        logical, allocatable :: haveTLmod(:)
!
!  Switch to indicate whether observations are present that
!  require the specific meta field be loaded and processed.
!
        logical, allocatable :: haveObsMeta(:)

# ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Variables in observation space needed Desroziers et al. (2005)
!  statistics at observation locations:
!-----------------------------------------------------------------------
!
!  4D-Var background error standard deviation [Mobs].
!
        real(r8), allocatable :: BgErr(:)
!
!  4D-Var NLM increment vector, analysis minus background {Mobs].
!
        real(r8), allocatable :: NLincrement(:)
!
!  4D-Var innovation vector, observations minus background [Mobs].
!
        real(r8), allocatable :: innovation(:)
!
!  4D-Var residual vector, observation minus analysis [Mobs].
!
        real(r8), allocatable :: residual(:)
# endif

# ifdef BGQC
!
!-----------------------------------------------------------------------
!  Background quality control of observations.
!-----------------------------------------------------------------------
!
!  Threshhold for innovations for quality control, [Mobs].
!
        real(r8), allocatable :: BgThresh(:)
!
!  Saved actual cost function, [0:Ninner].
!
        real(r8), allocatable :: Jact_S(:)
# endif
!
!-----------------------------------------------------------------------
!  Spatial convolutions parameters
!-----------------------------------------------------------------------
!
!  Initial conditions, surface forcing and model error covariance: Full
!  number of horizontal and vertical diffusion operator step for spatial
!  convolutions.
!
        integer, allocatable :: NHsteps(:,:)
        integer, allocatable :: NVsteps(:,:)
!
!  Boundary conditions error covariance: Full number of horizontal and
!  vertical diffusion operator step for spatial convolutions.
!
        integer, allocatable :: NHstepsB(:,:)
        integer, allocatable :: NVstepsB(:,:)
!
!  Initial conditions, surface forcing and model error covariance:
!  Horizontal and vertical diffusion operator time-step size for
!  spatial convolutions.
!
        real(r8), allocatable :: DTsizeH(:,:)
        real(r8), allocatable :: DTsizeV(:,:)
!
!  Boundary conditions error covariance: Horizontal and vertical
!  diffusion operator time-step size for spatial convolutions.
!
        real(r8), allocatable :: DTsizeHB(:,:)
        real(r8), allocatable :: DTsizeVB(:,:)
!
!  Minimum and maximum Horizontal and vertical diffusion coefficients
!  used in spatial convolutions.
!
        real(r8), allocatable :: KhMin(:)
        real(r8), allocatable :: KhMax(:)
        real(r8), allocatable :: KvMin(:)
        real(r8), allocatable :: KvMax(:)

# if defined I4DVAR || defined WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Conjugate gradient parameters.
!-----------------------------------------------------------------------
!
!  Conjugate gradient coefficients.
!
#  if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR  || \
      defined TL_R4DVAR
        real(r8), allocatable :: ad_cg_beta(:)

        real(r8), allocatable :: tl_cg_beta(:)
#  endif
        real(r8), allocatable :: cg_alpha(:,:)
        real(r8), allocatable :: cg_beta(:,:)
        real(r8), allocatable :: cg_tau(:,:)
!
!  Ritz values maximum error limit.
!
        real(r8) :: RitzMaxErr
!
!  Converged Ritz eigenvalues used to approximate Hessian matrix during
!  preconditioning.
!
        real(r8), allocatable :: Ritz(:)
!
!  Orthogonal eigenvectors of Lanczos recurrence term Q(k)T(k).
!
#  ifdef WEAK_CONSTRAINT
        real(r8), allocatable :: cg_zv(:,:,:)
#  else
        real(r8), allocatable :: cg_zv(:,:)
#  endif
!
!  Lanczos algorithm coefficients.
!
#  if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR  || \
      defined TL_R4DVAR
        real(r8), allocatable :: ad_cg_delta(:)
        real(r8), allocatable :: ad_cg_QG(:)
        real(r8), allocatable :: ad_cg_Gnorm(:)

        real(r8), allocatable :: tl_cg_delta(:)
        real(r8), allocatable :: tl_cg_QG(:)
#   ifdef RPCG
        real(r8), allocatable :: tl_cg_Gnorm_v(:)
        real(r8), allocatable :: ad_cg_Gnorm_v(:)
#   else
        real(r8), allocatable :: tl_cg_Gnorm(:)
#   endif
#  endif
        real(r8), allocatable :: cg_delta(:,:)
        real(r8), allocatable :: cg_gamma(:,:)
!
!  Initial gradient vector normalization factor.
!
        real(r8), allocatable :: cg_Gnorm(:)
        real(r8), allocatable :: cg_Gnorm_v(:)
        real(r8), allocatable :: cg_Gnorm_y(:)
!
!  Reduction in the gradient norm; excess cost function.
!
        real(r8), allocatable :: cg_Greduc(:,:)
!
!  Lanczos vector normalization factor.
!
        real(r8), allocatable :: cg_QG(:,:)
!
!  Lanczos recurrence symmetric, tridiagonal matrix, T(k).
!
        real(r8), allocatable :: cg_Tmatrix(:,:)
        real(r8), allocatable :: cg_zu(:,:)
!
!  Eigenvalues of the Lanczos recurrence term Q(k)T(k)
!  algorithm and their relative error (accuaracy).
!
        real(r8), allocatable :: cg_Ritz(:,:)
        real(r8), allocatable :: cg_RitzErr(:,:)
# endif
# if defined BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Balance operator variables.
!-----------------------------------------------------------------------
!
!  Logical switch to write balance operator reference free-surface.
!
        logical, allocatable :: wrtZetaRef(:)
# endif
# if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
     defined POSTERIOR_ERROR_F
!
!-----------------------------------------------------------------------
!  Posterior analysis error covariance matrix parameters.
!-----------------------------------------------------------------------
!
!  Lanczos algorithm parameters.
!
        real(r8), allocatable :: ae_delta(:,:)
        real(r8), allocatable :: ae_beta(:,:)
        real(r8), allocatable :: ae_zv(:,:)
        real(r8), allocatable :: ae_trace(:)
!
!  Initial gradient vector normalization factor.
!
        real(r8), allocatable :: ae_Gnorm(:)
!
!  Lanczos recurrence symmetric, tridiagonal matrix, T(k).
!
        real(r8), allocatable :: ae_Tmatrix(:,:)
        real(r8), allocatable :: ae_zu(:,:)
!
!  Eigenvalues of the Lanczos recurrence term Q(k)T(k)
!  algorithm and their relative error (accuaracy).
!
        real(r8), allocatable :: ae_Ritz(:,:)
        real(r8), allocatable :: ae_RitzErr(:,:)

#  if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
!
!  Initial or final full posterior error covariance Lanczos vectors
!  tridiagonal system.
!
        real(r8), allocatable :: zLanczos_coef(:,:)    ! coefficients
        real(r8), allocatable :: zLanczos_inv(:,:)     ! inverse
        real(r8), allocatable :: zLanczos_err(:,:)     ! error, identity
#  endif
# endif
# if defined HESSIAN_SV || defined HESSIAN_SO || defined HESSIAN_FSV
        real(r8), allocatable :: zLanczos_diag(:)    ! diagonal coefs
        real(r8), allocatable :: zLanczos_offdiag(:) ! off-diagonal coefs
#  ifdef LCZ_FINAL
        real(r8), allocatable :: GSmatrix(:,:) ! Gramm-Schmidt matrix
        real(r8), allocatable :: GSmatinv(:,:) ! Inverse Gramm-Schmidt matrix
#  endif
# endif
# if defined I4DVAR_ANA_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Lanczos algorithm coefficients.
!-----------------------------------------------------------------------
!
!  These coefficients are computed in the 4D-Var Lanczos data
!  assimilation algorithm and are used here to tangent linear
!  model initial conditions as weighted sum of the Lanczos
!  vectors.
!
        real(r8), allocatable :: cg_beta(:,:)
        real(r8), allocatable :: cg_delta(:,:)
        real(r8), allocatable :: cg_zu(:,:)
# endif
# if defined HESSIAN_SV || defined HESSIAN_SO || defined HESSIAN_FSV
!
!-----------------------------------------------------------------------
!  Lanczos algorithm coefficients.
!-----------------------------------------------------------------------
!
!  These coefficients are computed in the I4D-Var Lanczos data
!  assimilation algorithm and are used here to tangent linear
!  model initial conditions as weighted sum of the Lanczos
!  vectors.
!
        real(r8), allocatable :: cg_beta(:,:)
        real(r8), allocatable :: cg_delta(:,:)
# endif
!
!-----------------------------------------------------------------------
!  Descent algorithm parameters.
!-----------------------------------------------------------------------
!
!  Switch to compute Hessian approximated eigenvalues and eigenvectors.
!
        logical :: LhessianEV
!
!  Switch switch to activate hot start of subsquent outerloops in the
!  weak-constraint algorithms.
!
        logical :: LhotStart
!
!  Switch to activate conjugate gradient preconditioning.
!
        logical :: Lprecond
!
!  Switch to activate Ritz limited-memory preconditioning using the
!  eigenpairs approximation of the Hessian matrix (Tshimanga et al.,
!  2008).
!
        logical :: Lritz
# ifdef RPCG
!
!  Switch that controls computation of the augumented contributions
!  to the model error forcing terms in error_covariance.
!
        logical :: LaugWeak=.FALSE.
!
# endif
!
!  Number of converged Ritz eigenvalues. If input parameter NritzEV > 0,
!  then nConvRitz=NritzEV.  Therefore, only nConvRitz eigenpairs will
!  used for preconditioning and the RitzMaxErr threshold is ignored.
!
        integer :: NritzEV
        integer :: nConvRitz = 0
!
!  Weak contraint conjugate gradient norms.
!
        real(r8) :: cg_gammam1 = 0.0_r8
        real(r8) :: cg_sigmam1 = 0.0_r8
        real(r8) :: cg_rnorm = 0.0_r8
!
!  Upper bound on the relative error of the gradient when computing
!  Hessian eigenvectors.
!
        real(r8) :: GradErr
!
!  Maximum error bound on Hessian eigenvectors.  Note that even quite
!  innacurate eigenvectors are usefull for pre-condtioning purposes.
!
        real(r8) :: HevecErr

# ifdef FOUR_DVAR
!
!------------------------------------------------------------------------
!  Dot product parameters.
!------------------------------------------------------------------------
!
!  Dot product between tangent linear and adjoint vectors.
!
        real(r8) :: DotProduct
        real(r8) :: adDotProduct
!
!  Tangent linear model linearization check dot products.
!
        integer :: ig1count            ! counter for g1 dot product
        integer :: ig2count            ! counter for g2 dot product

        real(r8), dimension(1000) :: g1
        real(r8), dimension(1000) :: g2
# endif
# ifdef WEAK_CONSTRAINT
!
!------------------------------------------------------------------------
!  Weak constraint parameters.
!------------------------------------------------------------------------
!
!  Switch to activate writing of weak constraint forings
!  into the adjoint NetCDF file.
!
        logical, allocatable :: WRTforce(:)
#  ifdef RPM_RELAXATION
        logical, allocatable :: LweakRelax(:)
#  endif
#  ifdef RBL4DVAR_FCT_SENSITIVITY
        logical, allocatable :: LsenFCT(:)
#  endif
#  if defined RBL4DVAR_ANA_SENSITIVITY || \
      defined RBL4DVAR_FCT_SENSITIVITY || \
      defined R4DVAR_ANA_SENSITIVITY
!
!  Switch to activate reading of adjoint sensitivity initial conditions
!  for weak constraint 4DVAR observation sensitivity.
!
        logical, allocatable :: Lsen4DVAR(:)
#   ifdef OBS_SPACE
        logical, allocatable :: Lobspace(:)
#    ifndef OBS_IMPACT
        logical, allocatable :: LadjVAR(:)
#    endif
#   endif
#  endif
#  ifdef RBL4DVAR_FCT_SENSITIVITY
        integer, allocatable :: ntimes_ana(:)
        integer, allocatable :: ntimes_fct(:)
#  endif
!
!  Weak-constraint forcing time.  A new variable is required since this
!  forcing is delayed by nADJ time-steps.
!
        real(r8), allocatable :: ForceTime(:)
# endif

      CONTAINS

      SUBROUTINE allocate_fourdvar
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in module "mod_fourdvar"   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate module variables due to nested grids.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(Load_Zobs)) THEN
        allocate ( Load_Zobs(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(wrote_Zobs)) THEN
        allocate ( wrote_Zobs(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(NstateVar)) THEN
        allocate ( NstateVar(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(NobsVar)) THEN
        allocate ( NobsVar(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Ndatum)) THEN
        allocate ( Ndatum(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Nsurvey)) THEN
        allocate ( Nsurvey(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ObsSurvey)) THEN
        allocate ( ObsSurvey(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Nobs)) THEN
        allocate ( Nobs(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(NstrObs)) THEN
        allocate ( NstrObs(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(NendObs)) THEN
        allocate ( NendObs(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Nmethod)) THEN
        allocate ( Nmethod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(Rscheme)) THEN
        allocate ( Rscheme(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(Optimality)) THEN
        allocate ( Optimality(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ProcessObs)) THEN
        allocate ( ProcessObs(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

#  ifdef SP4DVAR
      IF (.not.allocated(Lsadd)) THEN
        allocate ( Lsadd(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#  endif

      IF (.not.allocated(wrtNLmod)) THEN
        allocate ( wrtNLmod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(wrtObsScale)) THEN
        allocate ( wrtObsScale(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(wrtRPmod)) THEN
        allocate ( wrtRPmod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(wrtTLmod)) THEN
        allocate ( wrtTLmod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(wrtMisfit)) THEN
        allocate ( wrtMisfit(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

# if defined I4DVAR_ANA_SENSITIVITY
#  ifdef OBS_IMPACT
      IF (.not.allocated(wrtIMPACT_TOT)) THEN
        allocate ( wrtIMPACT_TOT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#  endif
#  ifdef OBS_IMPACT_SPLIT
      IF (.not.allocated(wrtIMPACT_IC)) THEN
        allocate ( wrtIMPACT_IC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      IF (.not.allocated(wrtIMPACT_FC)) THEN
        allocate ( wrtIMPACT_FC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#   endif
#   if defined ADJUST_BOUNDARY
      IF (.not.allocated(wrtIMPACT_BC)) THEN
        allocate ( wrtIMPACT_BC(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#   endif
#  endif
# endif

      IF (.not.allocated(haveADmod)) THEN
        allocate ( haveADmod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(haveNLmod)) THEN
        allocate ( haveNLmod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(haveTLmod)) THEN
        allocate ( haveTLmod(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(haveObsMeta)) THEN
        allocate ( haveObsMeta(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

      IF (.not.allocated(KhMin)) THEN
        allocate ( KhMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(KhMax)) THEN
        allocate ( KhMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(KvMin)) THEN
        allocate ( KvMin(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(KvMax)) THEN
        allocate ( KvMax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF

# if defined BALANCE_OPERATOR
      IF (.not.allocated(wrtZetaRef)) THEN
        allocate ( wrtZetaRef(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
# endif

# ifdef WEAK_CONSTRAINT
      IF (.not.allocated(WRTforce)) THEN
        allocate ( WRTforce(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#  ifdef RPM_RELAXATION
      IF (.not.allocated(LweakRelax)) THEN
        allocate ( LweakRelax(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#  endif
#  ifdef RBL4DVAR_FCT_SENSITIVITY
      IF (.not.allocated(LsenFCT)) THEN
        allocate ( LsenFCT(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#  endif
#  if defined RBL4DVAR_ANA_SENSITIVITY || \
      defined RBL4DVAR_FCT_SENSITIVITY || \
      defined R4DVAR_ANA_SENSITIVITY
      IF (.not.allocated(Lsen4DVAR)) THEN
        allocate ( Lsen4DVAR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#   ifdef OBS_SPACE
      IF (.not.allocated(Lobspace)) THEN
        allocate ( Lobspace(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#    ifndef OBS_IMPACT
      IF (.not.allocated(LadjVAR)) THEN
        allocate ( LadjVAR(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#    endif
#   endif
#  endif
#  ifdef RBL4DVAR_FCT_SENSITIVITY
      IF (.not.allocated(ntimes_ana)) THEN
        allocate ( ntimes_ana(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
      IF (.not.allocated(ntimes_fct)) THEN
        allocate ( ntimes_fct(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
#  endif
      IF (.not.allocated(ForceTime)) THEN
        allocate ( ForceTime(Ngrids) )
        Dmem(1)=Dmem(1)+REAL(Ngrids,r8)
      END IF
# endif

      RETURN
      END SUBROUTINE allocate_fourdvar

# if defined RBL4DVAR_FCT_SENSITIVITY && defined OBS_SPACE
      SUBROUTINE deallocate_fourdvar
!
!=======================================================================
!                                                                      !
!  This routine deallocates several variables in module associated     !
!  with observation dimension parameters Ndatum, Nsurvey, and Mobs.    !
!  It is used in the observations impacts for the analysis-forecast    !
!  cycle.                                                              !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Deallocate structure variables for each nested grids.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (allocated(FOURDVAR(ng) % NobsSurvey)) THEN
          deallocate (FOURDVAR(ng) % NobsSurvey)
        END IF

        IF (allocated(FOURDVAR(ng) % SurveyTime)) THEN
          deallocate (FOURDVAR(ng) % SurveyTime)
        END IF

#  ifdef RPCG
        IF (allocated(FOURDVAR(ng) % cg_pxsave)) THEN
          deallocate (FOURDVAR(ng) % cg_pxsave)
        END IF

        IF (allocated(FOURDVAR(ng) % cg_pxout)) THEN
          deallocate (FOURDVAR(ng) % cg_pxout)
        END IF

        IF (allocated(FOURDVAR(ng) % cg_pxsum)) THEN
          deallocate (FOURDVAR(ng) % cg_pxsum)
        END IF

        IF (allocated(FOURDVAR(ng) % cg_Dpxsum)) THEN
          deallocate (FOURDVAR(ng) % cg_Dpxsum)
        END IF

        IF (allocated(FOURDVAR(ng) % tl_cg_pxsave)) THEN
          deallocate (FOURDVAR(ng) % tl_cg_pxsave)
        END IF

        IF (allocated(FOURDVAR(ng) % ad_cg_pxsave)) THEN
          deallocate (FOURDVAR(ng) % ad_cg_pxsave)
        END IF

#  endif
      END DO
!
!-----------------------------------------------------------------------
!  deallocate and initialize model and observation variables.
!-----------------------------------------------------------------------
!
#  ifdef CURVGRID
      IF (allocated(ObsAngler)) THEN
        deallocate (ObsAngler)
      END IF
#  endif

      IF (allocated(ObsType)) THEN
        deallocate (ObsType)
      END IF

      IF (allocated(ObsProv)) THEN
        deallocate (ObsProv)
      END IF

      IF (allocated(ObsErr)) THEN
        deallocate (ObsErr)
      END IF

      IF (allocated(ObsMeta)) THEN
        deallocate (ObsMeta)
      END IF

      IF (allocated(ObsScale)) THEN
        deallocate (ObsScale)
      END IF

#  ifndef WEAK_CONSTRAINT
      IF (allocated(ObsScaleGlobal)) THEN
        deallocate (ObsScaleGlobal)
      END IF
#  endif

      IF (allocated(ObsVetting)) THEN
        deallocate (ObsVetting)
      END IF

      IF (allocated(ObsVal)) THEN
        deallocate (ObsVal)
      END IF

      IF (allocated(Tobs)) THEN
        deallocate (Tobs)
      END IF

      IF (allocated(Xobs)) THEN
        deallocate (Xobs)
      END IF

      IF (allocated(Yobs)) THEN
        deallocate (Yobs)
      END IF

#  ifdef SOLVE3D
      IF (allocated(Zobs)) THEN
        deallocate (Zobs)
      END IF
#  endif

      IF (allocated(ADmodVal)) THEN
        deallocate (ADmodVal)
      END IF

      IF (allocated(NLmodVal)) THEN
        deallocate (NLmodVal)
      END IF

#  ifdef FOUR_DVAR
      IF (allocated(BgErr)) THEN
        deallocate (BgErr)
      END IF

      IF (allocated(NLincrement)) THEN
        deallocate (NLincrement)
      END IF

      IF (allocated(innovation)) THEN
        deallocate (innovation)
      END IF

      IF (allocated(residual)) THEN
        deallocate (residual)
      END IF
#  endif

#  ifdef BGQC
      IF (allocated(BgThresh)) THEN
        deallocate (BgThresh)
      END IF
#  endif

#  ifdef TLM_OBS
      IF (allocated(TLmodVal)) THEN
        deallocate (TLmodVal)
      END IF

#   if defined ARRAY_MODES       || defined CLIPPING    || \
       defined RBL4DVAR          || defined R4DVAR      || \
       defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
       defined TL_R4DVAR
      IF (allocated(TLmodVal_S)) THEN
        deallocate (TLmodVal_S)
      END IF

#    if defined RPCG              && \
       (defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
        defined TL_R4DVAR)

      IF (allocated(tl_TLmodVal_S)) THEN
        deallocate (tl_TLmodVal_S)
      END IF

      IF (allocated(ad_TLmodVal_S)) THEN
        deallocate (ad_TLmodVal_S)
      END IF
#    endif

#   endif

#   if defined SENSITIVITY_4DVAR || defined RBL4DVAR || \
       defined TL_RBL4DVAR

      IF (allocated(BCKmodVal)) THEN
        deallocate (BCKmodVal)
      END IF

#    ifdef RPCG
      IF (allocated(Hbk)) THEN
        deallocate (Hbk)
      END IF
#     if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
         defined TL_R4DVAR

      IF (allocated(tl_BCKmodVal)) THEN
        deallocate (tl_BCKmodVal)
      END IF
      IF (allocated(ad_BCKmodVal)) THEN
        deallocate (ad_BCKmodVal)
      END IF
      IF (allocated(tl_Hbk)) THEN
        deallocate (tl_Hbk)
      END IF
      IF (allocated(ad_Hbk)) THEN
        deallocate (ad_Hbk)
      END IF
#     endif
#    endif
#   endif
#  endif

#  ifdef WEAK_CONSTRAINT
!
!  Deallocate and initialize weak constraint conjugate gradient vectors.
!
#   if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
       defined TL_R4DVAR
      IF (allocated(ad_zcglwk)) THEN
        deallocate (ad_zcglwk)
      END IF

      IF (allocated(ad_zgrad0)) THEN
        deallocate (ad_zgrad0)
      END IF

      IF (allocated(ad_cg_innov)) THEN
        deallocate (ad_cg_innov)
      END IF

      IF (allocated(ad_cg_pxsave)) THEN
        deallocate (ad_cg_pxsave)
      END IF

      IF (allocated(tl_zcglwk)) THEN
        deallocate (tl_zcglwk)
      END IF

      IF (allocated(tl_zgrad0)) THEN
        deallocate (tl_zgrad0)
      END IF

      IF (allocated(tl_cg_innov)) THEN
        deallocate (tl_cg_innov)
      END IF

      IF (allocated(tl_cg_pxsave)) THEN
        deallocate (tl_cg_pxsave)
      END IF

      IF (allocated(tl_ObsVal)) THEN
        deallocate (tl_ObsVal)
      END IF

#    ifdef IMPACT_INNER
      IF (allocated(ad_ObsVal)) THEN
        deallocate (ad_ObsVal)
      END IF
#    else
      IF (allocated(ad_ObsVal)) THEN
        deallocate (ad_ObsVal)
      END IF
#    endif
#   endif

#   ifdef RPCG
      IF (allocated(zcglwk)) THEN
        deallocate (zcglwk)
      END IF

      IF (allocated(vcglwk)) THEN
        deallocate (vcglwk)
      END IF

      IF (allocated(vcglev)) THEN
        deallocate (vcglev)
      END IF

      IF (allocated(zgrad0)) THEN
        deallocate (zgrad0)
      END IF

      IF (allocated(vgrad0)) THEN
        deallocate (vgrad0)
      END IF

      IF (allocated(vgrad0s)) THEN
        deallocate (vgrad0s)
      END IF

      IF (allocated(cg_innov)) THEN
        deallocate (cg_innov)
      END IF

#    if defined RBL4DVAR_ANA_SENSITIVITY || \
        defined RBL4DVAR_FCT_SENSITIVITY || \
        defined R4DVAR_ANA_SENSITIVITY   || \
        defined TL_RBL4DVAR              || \
        defined TL_R4DVAR

      IF (allocated(ad_vcglwk)) THEN
        deallocate (ad_vcglwk)
      END IF

      IF (allocated(ad_vgrad0)) THEN
        deallocate (ad_vgrad0)
      END IF

      IF (allocated(ad_vgrad0s)) THEN
        deallocate (ad_vgrad0s)
      END IF

      IF (allocated(ad_Jb0)) THEN
        deallocate (ad_Jb0)
      END IF

      IF (allocated(ad_ObsErr)) THEN
        deallocate (ad_ObsErr)
      END IF

      IF (allocated(ad_ObsScale)) THEN
        deallocate (ad_ObsScale)
      END IF

      IF (allocated(ad_ObsVal)) THEN
        deallocate (ad_ObsVal)
      END IF

      IF (allocated(tl_vcglwk)) THEN
        deallocate (tl_vcglwk)
      END IF

      IF (allocated(tl_vgrad0)) THEN
        deallocate (tl_vgrad0)
      END IF

      IF (allocated(tl_vgrad0s)) THEN
        deallocate (tl_vgrad0s)
      END IF

      IF (allocated(tl_Jb0)) THEN
        deallocate (tl_Jb0)
      END IF

      IF (allocated(tl_ObsErr)) THEN
        deallocate (tl_ObsErr)
      END IF

      IF (allocated(tl_ObsScale)) THEN
        deallocate (tl_ObsScale)
      END IF

      IF (allocated(tl_ObsVal)) THEN
        deallocate (tl_ObsVal)
      END IF

      IF (allocated(tl_cg_Gnorm_v)) THEN
        deallocate (tl_cg_Gnorm_v)
      END IF

      IF (allocated(ad_cg_Gnorm_v)) THEN
        deallocate (ad_cg_Gnorm_v)
      END IF
#    endif

#   else

      IF (allocated(cg_pxsave)) THEN
        deallocate (cg_pxsave)
      END IF

      IF (allocated(zcglwk)) THEN
        deallocate (zcglwk)
      END IF

      IF (allocated(vcglev)) THEN
        deallocate (vcglev)
      END IF

      IF (allocated(zgrad0)) THEN
        deallocate (zgrad0)
      END IF

      IF (allocated(vgrad0)) THEN
        deallocate (vgrad0)
      END IF

      IF (allocated(vgrad0s)) THEN
        deallocate (vgrad0s)
      END IF

      IF (allocated(gdgw)) THEN
        deallocate (gdgw)
      END IF

      IF (allocated(cg_innov)) THEN
        deallocate (cg_innov)
      END IF

      IF (allocated(cg_pxsave)) THEN
        deallocate (cg_pxsave)
      END IF

#   endif
#  endif

      RETURN
      END SUBROUTINE deallocate_fourdvar
# endif

      SUBROUTINE initialize_fourdvar
!
!=======================================================================
!                                                                      !
!  This routine initializes several variables in module "mod_fourdvar" !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_scalars
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
!
      USE strings_mod, ONLY : FoundError
      USE strings_mod, ONLY : uppercase
!
!  Local variable declarations.
!
      integer :: my_NobsVar, i, icount, ng
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      integer :: LBi, UBi, LBj, UBj
      integer :: tile

      real(r8) :: size2d
# endif
      real(r8), parameter :: IniVal = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", initialize_fourdvar"

      SourceFile=MyFile

# ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Inquire observations NetCDF and determine the maximum dimension of
!  several observations arrays.
!-----------------------------------------------------------------------
!
!  Inquire about the size of the "datum" unlimitted dimension and
!  "survey" dimension.
!
      DO ng=1,Ngrids
        CALL netcdf_get_dim(ng, iNLM, TRIM(OBS(ng)%name))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Ndatum(ng)=0
        Nsurvey(ng)=0
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'datum') then
            Ndatum(ng)=dim_size(i)
          ELSE IF (TRIM(dim_name(i)).eq.'survey') THEN
            Nsurvey(ng)=dim_size(i)
          ELSE IF (TRIM(dim_name(i)).eq.'state_variable') THEN
            my_NobsVar=dim_size(i)
          END IF
        END DO
        IF (Ndatum(ng).eq.0) THEN
          WRITE (stdout,10) 'datum', TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (Nsurvey(ng).eq.0) THEN
          WRITE (stdout,10) 'survey', TRIM(OBS(ng)%name)
          exit_flag=2
          RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
#  if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
#   ifdef DISTRIBUTE
      tile=MyRank
#   else
      tile=0
#   endif
#  endif
!
!  Allocate structure (1:Ngrids).
!
      IF (.not.allocated( FOURDVAR)) THEN
        allocate ( FOURDVAR(Ngrids) )
      END IF
!
!  Allocate variables in structure.
!
      DO ng=1,Ngrids

#  if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
        LBi=BOUNDS(ng)%LBi(tile)
        UBi=BOUNDS(ng)%UBi(tile)
        LBj=BOUNDS(ng)%LBj(tile)
        UBj=BOUNDS(ng)%UBj(tile)

        size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
#  endif
#  ifdef SOLVE3D
        NstateVar(ng)=5+NT(ng)
#   ifdef ADJUST_STFLUX
        NstateVar(ng)=NstateVar(ng)+NT(ng)
#   endif
#  else
        NstateVar(ng)=3
#  endif
#  ifdef ADJUST_WSTRESS
        NstateVar(ng)=NstateVar(ng)+2
#  endif
        NobsVar(ng)=NstateVar(ng)+NextraObs

        IF (.not.allocated(FOURDVAR(ng) % NobsSurvey)) THEN
          allocate ( FOURDVAR(ng) % NobsSurvey(Nsurvey(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nsurvey(ng),r8)
          FOURDVAR(ng) % NobsSurvey = 0
        END IF

        IF (.not.allocated(FOURDVAR(ng) % SurveyTime)) THEN
          allocate ( FOURDVAR(ng) % SurveyTime(Nsurvey(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nsurvey(ng),r8)
          FOURDVAR(ng) % SurveyTime = IniVal
        END IF

#  ifdef FOUR_DVAR
        IF (.not.allocated(FOURDVAR(ng) % BackCost)) THEN
          allocate ( FOURDVAR(ng) % BackCost(0:NstateVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NstateVar(ng)+1,r8)
          FOURDVAR(ng) % BackCost = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % Cost0)) THEN
          allocate ( FOURDVAR(ng) % Cost0(Nouter) )
          Dmem(ng)=Dmem(ng)+REAL(Nouter,r8)
          FOURDVAR(ng) % Cost0 = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % CostFun)) THEN
          allocate ( FOURDVAR(ng) % CostFun(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % CostFun = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % CostFunOld)) THEN
          allocate ( FOURDVAR(ng) % CostFunOld(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % CostFunOld = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % CostNorm)) THEN
          allocate ( FOURDVAR(ng) % CostNorm(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % CostNorm = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % ObsCost)) THEN
          allocate ( FOURDVAR(ng) % ObsCost(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % ObsCost = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % DataPenalty)) THEN
          allocate ( FOURDVAR(ng) % DataPenalty(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % DataPenalty = IniVal
        END IF

#   if defined DATALESS_LOOPS    || defined RBL4DVAR || \
       defined SENSITIVITY_4DVAR || defined SP4DVAR  || \
       defined TL_RBL4DVAR

        IF (.not.allocated(FOURDVAR(ng) % NLPenalty)) THEN
          allocate ( FOURDVAR(ng) % NLPenalty(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % NLPenalty = IniVal
        END IF
#   endif

        IF (.not.allocated(FOURDVAR(ng) % NLobsCost)) THEN
          allocate ( FOURDVAR(ng) % NLobsCost(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % NLobsCost = IniVal
        END IF

#   if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
        IF (.not.allocated(FOURDVAR(ng) % bc_ak)) THEN
          allocate ( FOURDVAR(ng) % bc_ak(Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)
          FOURDVAR(ng) % bc_ak = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % bc_bk)) THEN
          allocate ( FOURDVAR(ng) % bc_bk(Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)
          FOURDVAR(ng) % bc_bk = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % zdf1)) THEN
          allocate ( FOURDVAR(ng) % zdf1(Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)
          FOURDVAR(ng) % zdf1 = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % zdf2)) THEN
          allocate ( FOURDVAR(ng) % zdf2(Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)
          FOURDVAR(ng) % zdf2 = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % zdf3)) THEN
          allocate ( FOURDVAR(ng) % zdf3(Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)
          FOURDVAR(ng) % zdf3 = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % pc_r2d)) THEN
          allocate ( FOURDVAR(ng) % pc_r2d(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % pc_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % rhs_r2d)) THEN
          allocate ( FOURDVAR(ng) % rhs_r2d(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % rhs_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % r2d_ref)) THEN
          allocate ( FOURDVAR(ng) % r2d_ref(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % r2d_ref = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % zeta_ref)) THEN
          allocate ( FOURDVAR(ng) % zeta_ref(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % zeta_ref = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % p_r2d)) THEN
          allocate ( FOURDVAR(ng) % p_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)*size2d
          FOURDVAR(ng) % p_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % r_r2d)) THEN
          allocate ( FOURDVAR(ng) % r_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)*size2d
          FOURDVAR(ng) % r_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % bp_r2d)) THEN
          allocate ( FOURDVAR(ng) % bp_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)*size2d
          FOURDVAR(ng) % bp_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % br_r2d)) THEN
          allocate ( FOURDVAR(ng) % br_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(Nbico(ng),r8)*size2d
          FOURDVAR(ng) % br_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % tl_rhs_r2d)) THEN
          allocate ( FOURDVAR(ng) % tl_rhs_r2d(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % tl_rhs_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % tl_r2d_ref)) THEN
          allocate ( FOURDVAR(ng) % tl_r2d_ref(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % tl_r2d_ref = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % tl_zeta_ref)) THEN
          allocate ( FOURDVAR(ng) % tl_zeta_ref(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % tl_zeta_ref = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % ad_rhs_r2d)) THEN
          allocate ( FOURDVAR(ng) % ad_rhs_r2d(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % ad_rhs_r2d = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % ad_r2d_ref)) THEN
          allocate ( FOURDVAR(ng) % ad_r2d_ref(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % ad_r2d_ref = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % ad_zeta_ref)) THEN
          allocate ( FOURDVAR(ng) % ad_zeta_ref(LBi:UBi,LBj:UBj) )
          Dmem(ng)=Dmem(ng)+size2d
          FOURDVAR(ng) % ad_zeta_ref = IniVal
        END IF
#   endif
#  endif

        IF (.not.allocated(FOURDVAR(ng) % ObsCount)) THEN
          allocate ( FOURDVAR(ng) % ObsCount(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % ObsCount = 0
        END IF

        IF (.not.allocated(FOURDVAR(ng) % ObsReject)) THEN
          allocate ( FOURDVAR(ng) % ObsReject(0:NobsVar(ng)) )
          Dmem(ng)=Dmem(ng)+REAL(NobsVar(ng)+1,r8)
          FOURDVAR(ng) % ObsReject = 0
        END IF

#  ifdef RPCG
        IF (.not.allocated(FOURDVAR(ng) % cg_pxout)) THEN
          allocate ( FOURDVAR(ng) % cg_pxout(Nouter) )
          Dmem(ng)=Dmem(ng)+REAL(Nouter)
          FOURDVAR(ng) % cg_pxout = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % cg_pxsave)) THEN
          allocate ( FOURDVAR(ng) % cg_pxsave(Ndatum(ng)+1) )
          Dmem(ng)=Dmem(ng)+REAL(Ndatum(ng)+1,r8)
          FOURDVAR(ng) % cg_pxsave = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % cg_pxsum)) THEN
          allocate ( FOURDVAR(ng) % cg_pxsum(Ndatum(ng)+1) )
          Dmem(ng)=Dmem(ng)+REAL(Ndatum(ng)+1,r8)
          FOURDVAR(ng) % cg_pxsum = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % cg_Dpxsum)) THEN
          allocate ( FOURDVAR(ng) % cg_Dpxsum(Ndatum(ng)+1) )
          Dmem(ng)=Dmem(ng)+REAL(Ndatum(ng)+1,r8)
          FOURDVAR(ng) % cg_Dpxsum = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % tl_cg_pxsave)) THEN
          allocate ( FOURDVAR(ng) % tl_cg_pxsave(Ndatum(ng)+1) )
          Dmem(ng)=Dmem(ng)+REAL(Ndatum(ng)+1,r8)
          FOURDVAR(ng) % tl_cg_pxsave = IniVal
        END IF

        IF (.not.allocated(FOURDVAR(ng) % ad_cg_pxsave)) THEN
          allocate ( FOURDVAR(ng) % ad_cg_pxsave(Ndatum(ng)+1) )
          Dmem(ng)=Dmem(ng)+REAL(Ndatum(ng)+1,r8)
          FOURDVAR(ng) % ad_cg_pxsave = IniVal
        END IF
#  endif

      END DO
!
!-----------------------------------------------------------------------
!  Read in number of observations available per survey at their times.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!
!  Read in number of observations available per survey.
!
        CALL netcdf_get_ivar (ng, iNLM, TRIM(OBS(ng)%name),             &
     &                        TRIM(Vname(1,idNobs)),                    &
     &                        FOURDVAR(ng)%NobsSurvey,                  &
     &                        start = (/1/),                            &
     &                        total = (/Nsurvey(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in the time of each observation survey.
!
        CALL netcdf_get_fvar (ng, iNLM, TRIM(OBS(ng)%name),             &
     &                        TRIM(Vname(1,idOday)),                    &
     &                        FOURDVAR(ng)%SurveyTime,                  &
     &                        start = (/1/),                            &
     &                        total = (/Nsurvey(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Determine maximum size of observation arrays.
!
#  ifdef WEAK_CONSTRAINT
        Mobs=MAXVAL(Ndatum)
#  else
        Mobs=0
        DO i=1,Nsurvey(ng)
          Mobs=MAX(Mobs, FOURDVAR(ng)%NobsSurvey(i))
        END DO
#  endif
      END DO
!
!-----------------------------------------------------------------------
!  Allocate and initialize model and observation variables.
!-----------------------------------------------------------------------
!
#  ifdef CURVGRID
      IF (.not.allocated(ObsAngler)) THEN
        allocate ( ObsAngler(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsAngler = IniVal
#  endif

      IF (.not.allocated(ObsType)) THEN
        allocate ( ObsType(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsType = 0

      IF (.not.allocated(ObsState2Type)) THEN
        allocate ( ObsState2Type(0:MAXVAL(NobsVar)) )
        Dmem(1)=Dmem(1)+REAL(MAXVAL(NobsVar)+1,r8)
      END IF
      ObsState2Type = 0

      IF (.not.allocated(ObsProv)) THEN
        allocate ( ObsProv(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsProv = 0

      IF (.not.allocated(ObsErr)) THEN
        allocate ( ObsErr(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsErr = IniVal

      IF (.not.allocated(ObsMeta)) THEN
        allocate ( ObsMeta(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsMeta = IniVal

      IF (.not.allocated(ObsScale)) THEN
        allocate ( ObsScale(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsScale = IniVal

#  ifndef WEAK_CONSTRAINT
      IF (.not.allocated(ObsScaleGlobal)) THEN
        allocate ( ObsScaleGlobal(MAXVAL(Ndatum)) )
        Dmem(1)=Dmem(1)+REAL(MAXVAL(Ndatum),r8)
      END IF
      ObsScaleGlobal = IniVal
#  endif

      IF (.not.allocated(ObsVetting)) THEN
        allocate ( ObsVetting(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsVetting = IniVal

      IF (.not.allocated(ObsVal)) THEN
        allocate ( ObsVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ObsVal = IniVal

      IF (.not.allocated(Tobs)) THEN
        allocate ( Tobs(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      Tobs = 0.0_dp

      IF (.not.allocated(Xobs)) THEN
        allocate ( Xobs(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      Xobs = IniVal

      IF (.not.allocated(Yobs)) THEN
        allocate ( Yobs(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      Yobs = IniVal

#  ifdef SOLVE3D
      IF (.not.allocated(Zobs)) THEN
        allocate ( Zobs(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      Zobs = IniVal
#  endif

      IF (.not.allocated(ADmodVal)) THEN
        allocate ( ADmodVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ADmodVal = IniVal

      IF (.not.allocated(NLmodVal)) THEN
        allocate ( NLmodVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      NLmodVal = IniVal

#  ifdef SP4DVAR
      IF (.not.allocated(ADmodVal_S)) THEN
        allocate ( ADmodVal_S(Mobs, Ninner+4) )
        Dmem(1)=Dmem(1)+REAL(Mobs*(Ninner+4),r8)
      END IF
      ADmodVal_S = IniVal

      IF (.not.allocated(harnoldi)) THEN
        allocate ( harnoldi(Ninner+1, Ninner, Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Ninner*Nouter,r8)
      END IF
      harnoldi = IniVal

      IF (.not.allocated(gmze)) THEN
        allocate ( gmze(Ninner, Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      gmze = IniVal

      IF (.not.allocated(cg_beta0)) THEN
        allocate ( cg_beta0(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      cg_beta0 = IniVal

      IF (.not.allocated(Jobs)) THEN
        allocate ( Jobs(0:Ninner,1:Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      Jobs = IniVal
#  endif

#  ifdef FOUR_DVAR
      IF (.not.allocated(BgErr)) THEN
        allocate ( BgErr(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      BgErr = IniVal

      IF (.not.allocated(NLincrement)) THEN
        allocate ( NLincrement(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      NLincrement = IniVal

      IF (.not.allocated(innovation)) THEN
        allocate ( innovation(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      innovation = IniVal

      IF (.not.allocated(residual)) THEN
        allocate ( residual(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      residual = IniVal
#  endif

#  ifdef BGQC
      IF (.not.allocated(BgThresh)) THEN
        allocate ( BgThresh(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      BgThresh = IniVal

      IF (.not.allocated(Jact_S)) THEN
        allocate ( Jact_S(0:Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner+1,r8)
      END IF
      Jact_S = IniVal
#  endif

#  ifdef TLM_OBS
      IF (.not.allocated(TLmodVal)) THEN
        allocate ( TLmodVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      TLmodVal = IniVal

#   if defined ARRAY_MODES       || defined CLIPPING  || \
       defined RBL4DVAR          || defined R4DVAR    || \
       defined SENSITIVITY_4DVAR || defined SP4DVAR   || \
       defined TL_RBL4DVAR       || defined TL_R4DVAR
!
!  The following arrays are only needed and used by the master node.
!  In order to avoid hogging memory penalties in the other nodes, they
!  are only allocated by the master node.
!
      IF (.not.allocated(TLmodVal_S)) THEN
        allocate ( TLmodVal_S(Mobs,Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Mobs*Ninner*Nouter,r8)
      END IF
      TLmodVal_S = IniVal

#    if defined RPCG              && \
       (defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
        defined TL_R4DVAR)
      IF (Master) THEN
        IF (.not.allocated(tl_TLmodVal_S)) THEN
          allocate ( tl_TLmodVal_S(Mobs,Ninner,Nouter) )
          Dmem(1)=Dmem(1)+REAL(Mobs*Ninner*Nouter,r8)
        END IF
        tl_TLmodVal_S = IniVal

        IF (.not.allocated(ad_TLmodVal_S)) THEN
          allocate ( ad_TLmodVal_S(Mobs,Ninner,Nouter) )
          Dmem(1)=Dmem(1)+REAL(Mobs*Ninner*Nouter,r8)
        END IF
        ad_TLmodVal_S = IniVal
      END IF
#    endif
#   endif

#   if defined ARRAY_MODES       || defined RBL4DVAR    || \
       defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR

      IF (.not.allocated(BCKmodVal)) THEN
        allocate ( BCKmodVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      BCKmodVal = IniVal

#    ifdef RPCG

      IF (.not.allocated(Hbk)) THEN
        allocate ( Hbk(Mobs,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Mobs*Nouter,r8)
      END IF
      Hbk = IniVal

#     if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
         defined TL_R4DVAR

      IF (.not.allocated(tl_BCKmodVal)) THEN
        allocate ( tl_BCKmodVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_BCKmodVal = IniVal

      IF (.not.allocated(tl_Hbk)) THEN
        allocate ( tl_Hbk(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_Hbk = IniVal

      IF (.not.allocated(ad_BCKmodVal)) THEN
        allocate ( ad_BCKmodVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_BCKmodVal = IniVal

      IF (.not.allocated(ad_Hbk)) THEN
        allocate ( ad_Hbk(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_Hbk = IniVal

#     endif
#    endif
#   endif
#  endif
!
!  Allocate and initialize observation types names and indices.
!  Notice that a mapping from state-to-type (ObsState2Type) and its
!  inverse type-to-state (ObsType2State) indices are needed because
!  the User is allowed to add extra-observation operators with
!  nonsequential type enumeration.
!
!  Both mapping arrays ObsState2Type and ObsType2State have a zero
!  array element to allow applications with no extra-observations
!  to work with their zero associated state index (as initialized
!  in mod_ncparam.F).  For example, if the index "isRadial" is not
!  redefined below, the following assignment
!
!        ObsState2Type(isRadial)=ObsState2Type(0)=0
!        ObsType2State(isRadial)=ObsType2State(0)=0
!
!  is still legal with isRadial=0. It avoids a Fortran segmentation
!  violation (i.e., subscript #1 of the array ObsState2Type has
!  value 0 which is less than the lower bound of 1).  Sorry for
!  the awkward logic but we need a generic way to specify extra-
!  observation operators.
!
      IF (.not.allocated(ObsName)) THEN
        allocate ( ObsName(MAXVAL(NobsVar)) )
        Dmem(1)=Dmem(1)+REAL(MAXVAL(NobsVar),r8)
      END IF
      icount=MAXVAL(NstateVar)
      ObsState2Type(0)=0
      DO i=1,icount                               ! 5+NT
        ObsState2Type(i)=i
        ObsName(i)=TRIM(Vname(1,idSvar(i)))
      END DO
      IF (NextraObs.gt.0) THEN
        DO i=1,NextraObs
          icount=icount+1
          ObsName(icount)=TRIM(ExtraName(i))
          SELECT CASE (TRIM(uppercase(ExtraName(i))))
            CASE ('RADIAL')
              isRadial=icount
              ObsState2Type(icount)=ExtraIndex(i)
          END SELECT
        END DO
      END IF
!
!  Set inverse mapping type-to-state.
!
      IF (.not.allocated(ObsType2State)) THEN
        allocate ( ObsType2State(0:MAXVAL(ObsState2Type)) )
        Dmem(1)=Dmem(1)+REAL(MAXVAL(ObsState2Type)+1,r8)
      END IF
      ObsType2State(0)=0
      ObsType2State(1:MAXVAL(ObsState2Type))=-1
      DO i=1,MAXVAL(NstateVar)
       ObsType2State(i)=i
      END DO
      IF (NextraObs.gt.0) THEN
        DO i=1,NextraObs
          icount=ExtraIndex(i)
          SELECT CASE (TRIM(uppercase(ExtraName(i))))
            CASE ('RADIAL')
              ObsType2State(icount)=isRadial
          END SELECT
        END DO
      END IF

#  ifdef WEAK_CONSTRAINT
!
!  Allocate and initialize weak constraint conjugate gradient vectors.
!
#   if defined RPCG
#    define MyMobs Mobs+1
#   else
#    define MyMobs Mobs
#   endif

      IF (.not.allocated(zcglwk)) THEN
        allocate ( zcglwk(MyMobs,Ninner+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*(Ninner+1)*Nouter,r8)
      END IF
      zcglwk = IniVal

#   ifdef RPCG
      IF (.not.allocated(vcglwk)) THEN
        allocate ( vcglwk(MyMobs,Ninner+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*(Ninner+1)*Nouter,r8)
      END IF
      vcglwk = IniVal

      IF (.not.allocated(Jb0)) THEN
        allocate ( Jb0(0:Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter+1,r8)
      END IF
      Jb0 = IniVal
#   endif

      IF (.not.allocated(vcglev)) THEN
        allocate ( vcglev(MyMobs,Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*Ninner*Nouter,r8)
      END IF
      vcglev = IniVal

      IF (.not.allocated(zgrad0)) THEN
        allocate ( zgrad0(MyMobs,Nouter) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*Nouter,r8)
      END IF
      zgrad0 = IniVal

      IF (.not.allocated(vgrad0)) THEN
        allocate ( vgrad0(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      vgrad0 = IniVal

     IF (.not.allocated(vgrad0s)) THEN
        allocate ( vgrad0s(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      vgrad0s = IniVal

      IF (.not.allocated(gdgw)) THEN
        allocate ( gdgw(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      gdgw = IniVal

      IF (.not.allocated(vgnorm)) THEN
        allocate ( vgnorm(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      vgnorm = IniVal

      IF (.not.allocated(cg_innov)) THEN
        allocate ( cg_innov(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      cg_innov = IniVal

      IF (.not.allocated(cg_dla)) THEN
        allocate ( cg_dla(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_dla = IniVal

#   ifndef RPCG

      IF (.not.allocated(cg_pxsave)) THEN
        allocate ( cg_pxsave(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      cg_pxsave = IniVal

#   endif

#   if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
       defined TL_R4DVAR

      IF (.not.allocated(ad_zcglwk)) THEN
        allocate ( ad_zcglwk(MyMobs,Ninner+1) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*(Ninner+1),r8)
      END IF
      ad_zcglwk = IniVal

      IF (.not.allocated(ad_zgrad0)) THEN
        allocate ( ad_zgrad0(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      ad_zgrad0 = IniVal

#    if defined RPCG              && \
       (defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR || \
        defined TL_R4DVAR)

      IF (.not.allocated(ad_vcglwk)) THEN
        allocate ( ad_vcglwk(MyMobs,Ninner+1) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*(Ninner+1),r8)
      END IF
      ad_vcglwk = IniVal

      IF (.not.allocated(ad_vgrad0)) THEN
        allocate ( ad_vgrad0(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      ad_vgrad0 = IniVal

      IF (.not.allocated(ad_vgrad0s)) THEN
        allocate ( ad_vgrad0s(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      ad_vgrad0s = IniVal

      IF (.not.allocated(ad_Jb0)) THEN
        allocate ( ad_Jb0(0:Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter+1,r8)
      END IF
      ad_Jb0 = IniVal

      IF (.not.allocated(ad_ObsErr)) THEN
        allocate ( ad_ObsErr(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_ObsErr = IniVal

      IF (.not.allocated(ad_ObsScale)) THEN
        allocate ( ad_ObsScale(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_ObsScale = IniVal

#    endif

      IF (.not.allocated(ad_cg_innov)) THEN
        allocate ( ad_cg_innov(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_cg_innov = IniVal

      IF (.not.allocated(ad_cg_pxsave)) THEN
        allocate ( ad_cg_pxsave(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_cg_pxsave = IniVal

      IF (.not.allocated(ad_cg_pxout)) THEN
        allocate ( ad_cg_pxout(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      ad_cg_pxout = IniVal

#    ifdef IMPACT_INNER

      IF (.not.allocated(ad_ObsVal)) THEN
        allocate ( ad_ObsVal(Mobs,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Mobs*Ninner,r8)
      END IF
      ad_ObsVal = IniVal

#    else

      IF (.not.allocated(ad_ObsVal)) THEN
        allocate ( ad_ObsVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      ad_ObsVal = IniVal

#    endif

      IF (.not.allocated(tl_zcglwk)) THEN
        allocate ( tl_zcglwk(MyMobs,Ninner+1) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*(Ninner+1),r8)
      END IF
      tl_zcglwk = IniVal

      IF (.not.allocated(tl_zgrad0)) THEN
        allocate ( tl_zgrad0(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      tl_zgrad0 = IniVal

#    ifdef RPCG

      IF (.not.allocated(tl_vcglwk)) THEN
        allocate ( tl_vcglwk(MyMobs,Ninner+1) )
        Dmem(1)=Dmem(1)+REAL((MyMobs)*(Ninner+1),r8)
      END IF
      tl_vcglwk = IniVal

      IF (.not.allocated(tl_vgrad0)) THEN
        allocate ( tl_vgrad0(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      tl_vgrad0 = IniVal

      IF (.not.allocated(tl_vgrad0s)) THEN
        allocate ( tl_vgrad0s(MyMobs) )
        Dmem(1)=Dmem(1)+REAL(MyMobs,r8)
      END IF
      tl_vgrad0s = IniVal

      IF (.not.allocated(tl_Jb0)) THEN
        allocate ( tl_Jb0(0:Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter+1,r8)
      END IF
      tl_Jb0 = IniVal

      IF (.not.allocated(tl_ObsErr)) THEN
        allocate ( tl_ObsErr(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_ObsErr = IniVal

      IF (.not.allocated(tl_ObsScale)) THEN
        allocate ( tl_ObsScale(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_ObsScale = IniVal

#    endif

      IF (.not.allocated(tl_cg_innov)) THEN
        allocate ( tl_cg_innov(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_cg_innov = IniVal

      IF (.not.allocated(tl_cg_pxsave)) THEN
        allocate ( tl_cg_pxsave(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_cg_pxsave = IniVal

      IF (.not.allocated(tl_cg_pxout)) THEN
        allocate ( tl_cg_pxout(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      tl_cg_pxout = IniVal

      IF (.not.allocated(tl_ObsVal)) THEN
        allocate ( tl_ObsVal(Mobs) )
        Dmem(1)=Dmem(1)+REAL(Mobs,r8)
      END IF
      tl_ObsVal = IniVal

#   endif
#  endif

# else
!
!-----------------------------------------------------------------------
!  Allocate other module variables.
!-----------------------------------------------------------------------
!
!  Set number of state variables.
!
      DO ng=1,Ngrids
#  ifdef SOLVE3D
        NstateVar(ng)=5+NT(ng)
#   ifdef ADJUST_STFLUX
        NstateVar(ng)=NstateVar(ng)+NT(ng)
#   endif
#  else
        NstateVar(ng)=3
#  endif
#  ifdef ADJUST_WSTRESS
        NstateVar(ng)=NstateVar(ng)+2
#  endif
      END DO
# endif
!
!  Allocate convolution parameters.
!
      IF (.not.allocated(NHsteps)) THEN
        allocate ( NHsteps(2,MstateVar) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar,r8)
      END IF
      NHsteps = 0

      IF (.not.allocated(NVsteps)) THEN
        allocate ( NVsteps(2,MstateVar) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar,r8)
      END IF
      NVsteps = 0

      IF (.not.allocated(DTsizeH)) THEN
        allocate ( DTsizeH(2,MstateVar) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar,r8)
      END IF
      DTsizeH = IniVal

      IF (.not.allocated(DTsizeV)) THEN
        allocate ( DTsizeV(2,MstateVar) )
        Dmem(1)=Dmem(1)+2.0_r8*REAL(MstateVar,r8)
      END IF
      DTsizeV = IniVal

      IF (.not.allocated(NVstepsB)) THEN
        allocate ( NVstepsB(4,MstateVar) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar,r8)
      END IF
      NVstepsB = 0

      IF (.not.allocated(NHstepsB)) THEN
        allocate ( NHstepsB(4,MstateVar) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar,r8)
      END IF
      NHstepsB = 0

      IF (.not.allocated(DTsizeHB)) THEN
        allocate ( DTsizeHB(4,MstateVar) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar,r8)
      END IF
      DTsizeHB = IniVal

      IF (.not.allocated(DTsizeVB)) THEN
        allocate ( DTsizeVB(4,MstateVar) )
        Dmem(1)=Dmem(1)+4.0_r8*REAL(MstateVar,r8)
      END IF
      DTsizeVB = IniVal


# if defined I4DVAR || defined WEAK_CONSTRAINT
!
!  Allocate conjugate gradient variables.
!
#  if defined SENSITIVITY_4DVAR || defined TL_RBL4DVAR  || \
      defined TL_R4DVAR
      IF (.not.allocated(ad_cg_beta)) THEN
        allocate ( ad_cg_beta(Ninner+1) )
        Dmem(1)=Dmem(1)+REAL(Ninner+1,r8)
      END IF
      ad_cg_beta = IniVal

      IF (.not.allocated(ad_cg_delta)) THEN
        allocate ( ad_cg_delta(Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner,r8)
      END IF
      ad_cg_delta = IniVal

      IF (.not.allocated(ad_cg_QG)) THEN
        allocate ( ad_cg_QG(Ninner+1) )
        Dmem(1)=Dmem(1)+REAL(Ninner+1,r8)
      END IF
      ad_cg_QG = IniVal

      IF (.not.allocated(ad_cg_Gnorm)) THEN
        allocate ( ad_cg_Gnorm(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      ad_cg_Gnorm = IniVal

      IF (.not.allocated(tl_cg_beta)) THEN
        allocate ( tl_cg_beta(Ninner+1) )
        Dmem(1)=Dmem(1)+REAL(Ninner+1,r8)
      END IF
      tl_cg_beta = IniVal

      IF (.not.allocated(tl_cg_delta)) THEN
        allocate ( tl_cg_delta(Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner,r8)
      END IF
      tl_cg_delta = IniVal

      IF (.not.allocated(tl_cg_QG)) THEN
        allocate ( tl_cg_QG(Ninner+1) )
        Dmem(1)=Dmem(1)+REAL(Ninner+1,r8)
      END IF
      tl_cg_QG = IniVal

#   ifdef RPCG
      IF (.not.allocated(tl_cg_Gnorm_v)) THEN
        allocate ( tl_cg_Gnorm_v(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      tl_cg_Gnorm_v = IniVal

      IF (.not.allocated(ad_cg_Gnorm_v)) THEN
        allocate ( ad_cg_Gnorm_v(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      ad_cg_Gnorm_v = IniVal

#   else
      IF (.not.allocated(tl_cg_Gnorm)) THEN
        allocate ( tl_cg_Gnorm(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      tl_cg_Gnorm = IniVal

#   endif
#  endif

      IF (.not.allocated(cg_alpha)) THEN
        allocate ( cg_alpha(0:Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      cg_alpha = IniVal

      IF (.not.allocated(cg_beta)) THEN
        allocate ( cg_beta(Ninner+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      cg_beta = IniVal

      IF (.not.allocated(cg_tau)) THEN
        allocate ( cg_tau(0:Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      cg_tau = IniVal

      IF (.not.allocated(Ritz)) THEN
        allocate ( Ritz(Ninner+1) )
        Dmem(1)=Dmem(1)+REAL(Ninner+1,r8)
      END IF
      Ritz = IniVal

#  if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
      defined POSTERIOR_ERROR_F

      IF (.not.allocated(ae_beta)) THEN
        allocate ( ae_beta(NpostI+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((NpostI+1)*Nouter,r8)
      END IF
      ae_beta = IniVal

      IF (.not.allocated(ae_zv)) THEN
        allocate ( ae_zv(NpostI,NpostI) )
        Dmem(1)=Dmem(1)+REAL(NpostI*NpostI,r8)
      END IF
      ae_zv = IniVal

      IF (.not.allocated(ae_delta)) THEN
        allocate ( ae_delta(NpostI,Nouter) )
        Dmem(1)=Dmem(1)+REAL(NpostI*Nouter,r8)
      END IF
      ae_delta = IniVal

      IF (.not.allocated(ae_trace)) THEN
        allocate ( ae_trace(NpostI+1) )
        Dmem(1)=Dmem(1)+REAL(NpostI+1,r8)
      END IF
      ae_trace = IniVal

      IF (.not.allocated(ae_Gnorm)) THEN
        allocate ( ae_Gnorm(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      ae_Gnorm = IniVal

      IF (.not.allocated(ae_Tmatrix)) THEN
        allocate ( ae_Tmatrix(NpostI,3) )
        Dmem(1)=Dmem(1)+3.0_r8*REAL(NpostI,r8)
      END IF
      ae_Tmatrix = IniVal

      IF (.not.allocated(ae_Ritz)) THEN
        allocate ( ae_Ritz(NpostI,Nouter) )
        Dmem(1)=Dmem(1)+REAL(NpostI*Nouter,r8)
      END IF
      ae_Ritz = IniVal

      IF (.not.allocated(ae_RitzErr)) THEN
        allocate ( ae_RitzErr(NpostI,Nouter) )
        Dmem(1)=Dmem(1)+REAL(NpostI*Nouter,r8)
      END IF
      ae_RitzErr = IniVal

#   if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      IF (.not.allocated(zLanczos_coef)) THEN
        allocate ( zLanczos_coef(Ninner,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner,r8)
      END IF
      zLanczos_coef = IniVal

      IF (.not.allocated(zLanczos_inv)) THEN
        allocate ( zLanczos_inv(Ninner,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner,r8)
      END IF
      zLanczos_inv = IniVal

      IF (.not.allocated(zLanczos_err)) THEN
        allocate ( zLanczos_err(Ninner,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner,r8)
      END IF
      zLanczos_err = IniVal
#   endif
#  endif

#  ifdef WEAK_CONSTRAINT
      IF (.not.allocated(cg_zv)) THEN
        allocate ( cg_zv(Ninner,Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner*Nouter,r8)
      END IF
#  else
      IF (.not.allocated(cg_zv)) THEN
        allocate ( cg_zv(Ninner,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner,r8)
      END IF
#  endif
      cg_zv = IniVal

      IF (.not.allocated(cg_delta)) THEN
        allocate ( cg_delta(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_delta = IniVal

      IF (.not.allocated(cg_gamma)) THEN
        allocate ( cg_gamma(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_gamma = IniVal

      IF (.not.allocated(cg_Gnorm)) THEN
        allocate ( cg_Gnorm(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      cg_Gnorm = IniVal

      IF (.not.allocated(cg_Gnorm_v)) THEN
        allocate ( cg_Gnorm_v(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      cg_Gnorm_v = IniVal

      IF (.not.allocated(cg_Gnorm_y)) THEN
        allocate ( cg_Gnorm_y(Nouter) )
        Dmem(1)=Dmem(1)+REAL(Nouter,r8)
      END IF
      cg_Gnorm_y = IniVal

      IF (.not.allocated(cg_Greduc)) THEN
        allocate ( cg_Greduc(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_Greduc = IniVal

      IF (.not.allocated(cg_QG)) THEN
        allocate ( cg_QG(Ninner+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      cg_QG = IniVal

      IF (.not.allocated(cg_Tmatrix)) THEN
        allocate ( cg_Tmatrix(Ninner,3) )
        Dmem(1)=Dmem(1)+3.0_r8*REAL(Ninner,r8)
      END IF
      cg_Tmatrix = IniVal

      IF (.not.allocated(cg_Ritz)) THEN
        allocate ( cg_Ritz(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_Ritz = IniVal

      IF (.not.allocated(cg_RitzErr)) THEN
        allocate ( cg_RitzErr(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_RitzErr = IniVal

      IF (.not.allocated(cg_zu)) THEN
        allocate ( cg_zu(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_zu = IniVal

# endif
# if defined I4DVAR_ANA_SENSITIVITY
!
!  Allocate Lanczos algorithm coefficients.
!
      IF (.not.allocated(cg_beta)) THEN
        allocate ( cg_beta(Ninner+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      cg_beta = IniVal

      IF (.not.allocated(cg_delta)) THEN
        allocate ( cg_delta(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_delta = IniVal

      IF (.not.allocated(cg_zu)) THEN
        allocate ( cg_zu(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_zu = IniVal
# endif
# if defined HESSIAN_SV || defined HESSIAN_SO || defined HESSIAN_FSV
!
!  Allocate Lanczos algorithm coefficients.
!
      IF (.not.allocated(cg_beta)) THEN
        allocate ( cg_beta(Ninner+1,Nouter) )
        Dmem(1)=Dmem(1)+REAL((Ninner+1)*Nouter,r8)
      END IF
      cg_beta = IniVal

      IF (.not.allocated(cg_delta)) THEN
        allocate ( cg_delta(Ninner,Nouter) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Nouter,r8)
      END IF
      cg_delta = IniVal

      IF (.not.allocated(zLanczos_diag)) THEN
        allocate ( zLanczos_diag(Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner,r8)
      END IF
      zLanczos_diag = IniVal

      IF (.not.allocated(zLanczos_offdiag)) THEN
        allocate ( zLanczos_offdiag(Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner,r8)
      END IF
      zLanczos_offdiag = IniVal
#  ifdef LCZ_FINAL

      IF (.not.allocated(GSmatrix)) THEN
        allocate ( GSmatrix(Ninner,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner,r8)
      END IF
      GSmatrix = IniVal

      IF (.not.allocated(GSmatinv)) THEN
        allocate ( GSmatinv(Ninner,Ninner) )
        Dmem(1)=Dmem(1)+REAL(Ninner*Ninner,r8)
      END IF
      GSmatinv = IniVal
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Initialize various variables.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        Load_Zobs(ng)=.FALSE.
        ProcessObs(ng)=.FALSE.
# ifdef SP4DVAR
        Lsadd(ng)=.FALSE.
# endif
        haveObsMeta(ng)=.FALSE.
        wrote_Zobs(ng)=.FALSE.
        wrtMisfit(ng)=.FALSE.
        wrtNLmod(ng)=.FALSE.
        wrtObsScale(ng)=.FALSE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
# ifdef BALANCE_OPERATOR
        wrtZetaRef(ng)=.FALSE.
# endif
# ifdef I4DVAR_ANA_SENSITIVITY
#  ifdef OBS_IMPACT
        wrtIMPACT_TOT(ng)=.FALSE.
#  endif
#  ifdef OBS_IMPACT_SPLIT
        wrtIMPACT_IC(ng)=.FALSE.
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        wrtIMPACT_FC(ng)=.FALSE.
#   endif
#   if defined ADJUST_BOUNDARY
        wrtIMPACT_BC(ng)=.FALSE.
#   endif
#  endif
# endif
        KhMin(ng)=1.0_r8
        KhMax(ng)=1.0_r8
        KvMin(ng)=1.0_r8
        KvMax(ng)=1.0_r8
        Optimality(ng)=0.0_r8
# ifdef WEAK_CONSTRAINT
        ForceTime(ng)=0.0_r8
# endif
      END DO

# ifdef OBSERVATIONS
!
 10   FORMAT (/,' MOD_FOURDVAR - error inquiring dimension: ',a,2x,     &
     &          ' in input NetCDF file: ',a)
# endif

      RETURN
      END SUBROUTINE initialize_fourdvar
#endif
      END MODULE mod_fourdvar
