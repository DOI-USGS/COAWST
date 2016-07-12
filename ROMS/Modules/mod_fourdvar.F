#include "cppdefs.h"
      MODULE mod_fourdvar

#if defined FOUR_DVAR || defined VERIFICATION
!
!svn $Id: mod_fourdvar.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Variational data assimilation variables:                            !
!                                                                      !
!  ADmodVal     Adjoint model values at observation locations.         !
!  BackCost     Current background cost function misfit (mean squared  !
!                 difference) between model and background state, Jb.  !
!  CostFun      Current iteration total cost function (background      !
!                 plus observations), J.                               !
!  CostFunOld   Previous iteration total cost function (background     !
!                 plus observations), J.                               !
!  CostNorm     Total cost function normalization scales (minimization !
!                 starting value, Jb+Jo).                              !
!  Cost0        The total cost function at the beggining of the inner  !
!                 loop (inner=0) for each outer loop.                  !
!  DTsizeH      Horizontal diffusion time-step size for spatial        !
!                 convolutions.                                        !
!  DTsizeV      Vertical diffusion time-step size for spatial          !
!                 convolutions.                                        !
!  GradErr      Upper bound on relatice error of the gradient.         !
!  HevecErr     Maximum error bound on Hessian eigenvectors.           !
!  KhMax        Maximum  horizontal diffusion coefficient.             !
!  KhMin        Minimum  horizontal diffusion coefficient.             !
!  KvMax        Maximum  vertical diffusion coefficient.               !
!  KvMin        Minimum  vertical diffusion coefficient.               !
!  Load_Zobs    Logical switch indicating that input Zobs is negative  !
!                 and fractional vertical positions are computed in    !
!                 "extract_obs3d" and written to observation NetCDF    !
!                 file for latter use.                                 !
!  LhessianEV   Logical switch to compute Hessian eigenvectors.        !
!  LhotStart    Logical switch to activate hot start of subsquent      !
!                 outer loops in the weak-constraint algorithms.       !
!  Lprecond     Logical switch for conjugate gradient preconditioning. !
!  Lritz        Logical switch for Ritz limited-memory preconditioning.!
#ifdef RPCG
!  LaugWeak     Logical switch for computing augmented model error     !
!               terms in the routine "error_covariance".               !
#endif
!  nConvRitz    Number of converged Ritz eigenvalues.                  !
!  NLmodVal     Nonlinear model values at observation locations.       !
!  NHsteps      Full number of horizontal diffusion steps for spatial  !
!                 convolutions.                                        !
!  NLobsCost    Current NLM observation cost function misfit (mean     !
!                 squared difference) between NLM and observations.    !
!  NpostI       Number of Lanczos iterations used for the posterior    !
!                 analysis error covariance matrix estimation.         !
!  NritzEV      If preconditioning, number of eigenpairs to use.       !
!  NVsteps      Full number of vertical diffusion steps for spatial    !
!                 convolutions.                                        !
!  ObsCost      Current observation cost function misfit (mean squared !
!                 difference) between model and observations, Jo.      !
!  ObsCount     Current observation counts per state variable.         !
!  ObsErr       Observation error.                                     !
!  ObsProv      Observation provenance flags.                          !
!  ObsReject    Current rejected observation counts per state variable.!
!  ObsScale     Observation scale used during screenning and/or        !
!                 normalization of the observations.                   !
!  ObsType      Observation type identifier.                           !
!  ObsVal       Observation values.                                    !
!  ObsVar       Global observation variance for each state variable.   !
!  Optimality   normalized, optimal cost function minimum.             !
!  Ritz         Ritz eigenvalues to compute approximated Hessian.      !
!  RitzMaxErr   Ritz values maximum error limit.                       !
!  TLmodVal     Tangent linear model values at observation locations.  !
!  Vdecay       Covariance vertical decorrelation scale (m).           !
!  Tobs         Observations time (days).                              !
!  Xobs         Observations X-locations (grid coordinates).           !
!  Yobs         Observations Y-locations (grid coordinates).           !
!  Zobs         Observations Z-locations (grid coordinates).           !
!  cg_Gnorm     Initial gradient vector normalization factor.          !
!  cg_Greduc    Reduction in the gradient norm; excess cost function.  !
!  cg_QG        Lanczos vector normalization factor.                   !
!  cg_Ritz      Eigenvalues of the Lanczos recurrence term Q(k)T(k).   !
!  cg_RitzErr   Eigenvalues relative error.                            !
!  cg_Tmatrix   Lanczos recurrence symmetric, tridiagonal matrix.      !
!  cg_alpha     Conjugate gradient coefficient.                        !
!  cg_beta      Conjugate gradient coefficient.                        !
!  cg_delta     Lanczos algorithm coefficient.                         !
!  cg_gamma     Lanczos algorithm coefficient.                         !
!  cg_tau       Conjugate gradient coefficient.                        !
!  cg_zu        Lanczos tridiagonal matrix, upper diagonal elements.   !
!  cg_zv        Eigenvectors of Lanczos recurrence term Q(k)T(k).      !
!  haveADmod    Logical switch indicating that there is representer    !
!                 coefficients available to process in DAV(ng)%name    !
!                 file.                                                !
!  haveNLmod    Logical switch indicating that there is nonlinear      !
!                 model data available to process in DAV(ng)%name      !
!                 file.                                                !
!  haveTLmod    Logical switch indicating that there is tangent        !
!                 model data available to process in DAV(ng)%name      !
!                 file.                                                !
# ifdef BGQC
!                                                                      !
!  Background quality control of observations:                         !
!                                                                      !
!  BgErr        Background error standard deviation values at          !
!                 observation locations.                               !
!  BgThresh     Number of standard deviations threshold used in the    !
!                 background quality control of observations.          !
!  Jact_S       Saved actual cost function.                            !
#endif
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

          integer , pointer :: NobsSurvey(:)
          integer , pointer :: ObsCount(:)
          integer , pointer :: ObsReject(:)
#  ifdef FOUR_DVAR
          real(r8), pointer :: BackCost(:)
          real(r8), pointer :: Cost0(:)
          real(r8), pointer :: CostFun(:)
          real(r8), pointer :: CostFunOld(:)
          real(r8), pointer :: CostNorm(:)
          real(r8), pointer :: ObsCost(:)
          real(r8), pointer :: DataPenalty(:)
#   if defined DATALESS_LOOPS || defined TL_W4DPSAS          || \
       defined W4DPSAS        || defined W4DPSAS_SENSITIVITY
          real(r8), pointer :: NLPenalty(:)
#   endif
          real(r8), pointer :: NLobsCost(:)

#   if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
          real(r8), pointer :: bc_ak(:)
          real(r8), pointer :: bc_bk(:)

          real(r8), pointer :: zdf1(:)
          real(r8), pointer :: zdf2(:)
          real(r8), pointer :: zdf3(:)

          real(r8), pointer :: pc_r2d(:,:)
          real(r8), pointer :: rhs_r2d(:,:)
          real(r8), pointer :: r2d_ref(:,:)
          real(r8), pointer :: zeta_ref(:,:)

          real(r8), pointer :: p_r2d(:,:,:)
          real(r8), pointer :: r_r2d(:,:,:)
          real(r8), pointer :: bp_r2d(:,:,:)
          real(r8), pointer :: br_r2d(:,:,:)

          real(r8), pointer :: tl_rhs_r2d(:,:)
          real(r8), pointer :: tl_r2d_ref(:,:)
          real(r8), pointer :: tl_zeta_ref(:,:)

          real(r8), pointer :: ad_rhs_r2d(:,:)
          real(r8), pointer :: ad_r2d_ref(:,:)
          real(r8), pointer :: ad_zeta_ref(:,:)
#   endif
#  endif
          real(r8), pointer :: ObsVar(:)
          real(r8), pointer :: SurveyTime(:)

#  ifdef RPCG
          real(r8), pointer :: cg_pxsave(:)
          real(r8), pointer :: cg_pxsum(:)
          real(r8), pointer :: cg_Dpxsum(:)
#  endif

        END TYPE T_FOURDVAR

        TYPE (T_FOURDVAR), allocatable :: FOURDVAR(:)
!
!  Define other module arrays.
!
        integer,  allocatable :: ObsType(:)
        integer,  allocatable :: ObsProv(:)

        real(r8), allocatable :: ObsErr(:)
        real(r8), allocatable :: ObsScale(:)
        real(r8), allocatable :: ObsVal(:)
        real(r8), allocatable :: Tobs(:)
        real(r8), allocatable :: Xobs(:)
        real(r8), allocatable :: Yobs(:)
#  ifdef SOLVE3D
        real(r8), allocatable :: Zobs(:)
#  endif
        real(r8), allocatable :: ADmodVal(:)
        real(r8), allocatable :: NLmodVal(:)
#  ifdef TLM_OBS
        real(r8), allocatable :: TLmodVal(:)
#   if defined SENSITIVITY_4DVAR    || defined TL_W4DPSAS  || \
       defined TL_W4DVAR            || defined W4DPSAS     || \
       defined W4DVAR               || defined ARRAY_MODES || \
       defined CLIPPING
        real(r8), allocatable :: TLmodVal_S(:,:,:)
#   endif
#   if defined TL_W4DPSAS           || defined W4DPSAS    || \
       defined W4DPSAS_SENSITIVITY

        real(r8), allocatable :: BCKmodVal(:)
        real(r8), allocatable :: Hbk(:,:)
#   endif
#  endif
#  ifdef WEAK_CONSTRAINT
        real(r8), allocatable :: zcglwk(:,:,:)
#   ifdef RPCG
        real(r8), allocatable :: vcglwk(:,:,:)
        real(r8), allocatable :: Jb0(:)
        real(r8), allocatable :: Hdx(:)
#   endif
        real(r8), allocatable :: vcglev(:,:,:)
        real(r8), allocatable :: zgrad0(:,:)
        real(r8), allocatable :: vgrad0(:)
        real(r8), allocatable :: vgrad0s(:)
        real(r8), allocatable :: gdgw(:)
        real(r8), allocatable :: vgnorm(:)
        real(r8), allocatable :: cg_innov(:)
        real(r8), allocatable :: cg_dla(:,:)
#     ifndef RPCG
        real(r8), allocatable :: cg_pxsave(:)
#     endif
#   if defined TL_W4DPSAS          || defined TL_W4DVAR          || \
       defined W4DPSAS_SENSITIVITY || defined W4DVAR_SENSITIVITY
        real(r8), allocatable :: ad_zcglwk(:,:)
        real(r8), allocatable :: ad_zgrad0(:)
        real(r8), allocatable :: ad_cg_innov(:)
        real(r8), allocatable :: ad_cg_pxsave(:)
        real(r8), allocatable :: ad_ObsVal(:)

        real(r8), allocatable :: tl_zcglwk(:,:)
        real(r8), allocatable :: tl_zgrad0(:)
        real(r8), allocatable :: tl_cg_innov(:)
        real(r8), allocatable :: tl_cg_pxsave(:)
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
!
!  Switch to activate writting of nonlinear model values at
!  observations locations into observations NetCDF file.
!
        logical, allocatable :: wrtNLmod(:)
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

# if defined IS4DVAR_SENSITIVITY
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

# ifdef BGQC
!
!-----------------------------------------------------------------------
!  Background quality control of observations.
!-----------------------------------------------------------------------
!
!  Background error standard deviation values at observation locations,
!  [Mobs].
!
        real(r8), allocatable :: BgErr(:)
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

# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!-----------------------------------------------------------------------
!  Conjugate gradient parameters.
!-----------------------------------------------------------------------
!
!  Conjugate gradient coefficients.
!
#  if defined SENSITIVITY_4DVAR || defined TL_W4DPSAS  || \
      defined TL_W4DVAR
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
#  if defined SENSITIVITY_4DVAR || defined TL_W4DPSAS  || \
      defined TL_W4DVAR
        real(r8), allocatable :: ad_cg_delta(:)
        real(r8), allocatable :: ad_cg_QG(:)
        real(r8), allocatable :: ad_cg_Gnorm(:)

        real(r8), allocatable :: tl_cg_delta(:)
        real(r8), allocatable :: tl_cg_QG(:)
        real(r8), allocatable :: tl_cg_Gnorm(:)
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
# if defined IS4DVAR_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Lanczos algorithm coefficients.
!-----------------------------------------------------------------------
!
!  These coefficients are computed in the 4DVAR Lanczos data
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
!  These coefficients are computed in the IS4DVAR Lanczos data
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
#ifdef RPCG
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
#  if defined W4DPSAS_SENSITIVITY  || defined W4DVAR_SENSITIVITY
!
!  Switch to activate reading of adjoint sensitivity initial conditions
!  for weak constraint 4DVAR observation sensitivity.
!
        logical, allocatable :: LsenPSAS(:)
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

# ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  Allocate module structure.
!-----------------------------------------------------------------------
!
      allocate ( FOURDVAR(Ngrids) )
# endif
!
!-----------------------------------------------------------------------
!  Allocate module variables due to nested grids.
!-----------------------------------------------------------------------
!
      allocate ( Load_Zobs(Ngrids) )
      allocate ( wrote_Zobs(Ngrids) )

      allocate ( NstateVar(Ngrids) )
      allocate ( Ndatum(Ngrids) )
      allocate ( Nsurvey(Ngrids) )
      allocate ( ObsSurvey(Ngrids) )
      allocate ( Nobs(Ngrids) )

      allocate ( NstrObs(Ngrids) )
      allocate ( NendObs(Ngrids) )

      allocate ( Nmethod(Ngrids) )
      allocate ( Rscheme(Ngrids) )

      allocate ( Optimality(Ngrids) )
      allocate ( ProcessObs(Ngrids) )

      allocate ( wrtNLmod(Ngrids) )
      allocate ( wrtRPmod(Ngrids) )
      allocate ( wrtTLmod(Ngrids) )
      allocate ( wrtMisfit(Ngrids) )

# if defined IS4DVAR_SENSITIVITY
#  ifdef OBS_IMPACT
      allocate ( wrtIMPACT_TOT(Ngrids) )
#  endif
#  ifdef OBS_IMPACT_SPLIT
      allocate ( wrtIMPACT_IC(Ngrids) )
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      allocate ( wrtIMPACT_FC(Ngrids) )
#   endif
#   if defined ADJUST_BOUNDARY
      allocate ( wrtIMPACT_BC(Ngrids) )
#   endif
#  endif
# endif

      allocate ( haveADmod(Ngrids) )
      allocate ( haveNLmod(Ngrids) )
      allocate ( haveTLmod(Ngrids) )

      allocate ( KhMin(Ngrids) )
      allocate ( KhMax(Ngrids) )
      allocate ( KvMin(Ngrids) )
      allocate ( KvMax(Ngrids) )

# if defined BALANCE_OPERATOR
      allocate ( wrtZetaRef(Ngrids) )
# endif

# ifdef WEAK_CONSTRAINT
      allocate ( WRTforce(Ngrids) )
#  ifdef RPM_RELAXATION
      allocate ( LweakRelax(Ngrids) )
#  endif
#  if defined W4DPSAS_SENSITIVITY  || defined W4DVAR_SENSITIVITY
      allocate ( LsenPSAS(Ngrids) )
#  endif
      allocate ( ForceTime(Ngrids) )
# endif

      RETURN
      END SUBROUTINE allocate_fourdvar

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
!  Local variable declarations.
!
      integer :: my_NstateVar, i, ng
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      integer :: LBi, UBi, LBj, UBj
      integer :: tile
# endif
      real(r8), parameter :: IniVal = 0.0_r8
!
      SourceFile='mod_fourdvar.F, initialize_fourdvar'

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
        IF (exit_flag.ne.NoError) RETURN
!
        Ndatum(ng)=0
        Nsurvey(ng)=0
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'datum') then
            Ndatum(ng)=dim_size(i)
          ELSE IF (TRIM(dim_name(i)).eq.'survey') THEN
            Nsurvey(ng)=dim_size(i)
          ELSE IF (TRIM(dim_name(i)).eq.'state_variable') THEN
            my_NstateVar=dim_size(i)
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

      DO ng=1,Ngrids

#  if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
        LBi=BOUNDS(ng)%LBi(tile)
        UBi=BOUNDS(ng)%UBi(tile)
        LBj=BOUNDS(ng)%LBj(tile)
        UBj=BOUNDS(ng)%UBj(tile)
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

        allocate ( FOURDVAR(ng) % NobsSurvey(Nsurvey(ng)) )
        FOURDVAR(ng) % NobsSurvey = 0

        allocate ( FOURDVAR(ng) % SurveyTime(Nsurvey(ng)) )
        FOURDVAR(ng) % SurveyTime = IniVal

#  ifdef FOUR_DVAR
        allocate ( FOURDVAR(ng) % BackCost(0:NstateVar(ng)) )
        FOURDVAR(ng) % BackCost = IniVal

        allocate ( FOURDVAR(ng) % Cost0(Nouter) )
        FOURDVAR(ng) % Cost0 = IniVal

        allocate ( FOURDVAR(ng) % CostFun(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostFun = IniVal

        allocate ( FOURDVAR(ng) % CostFunOld(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostFunOld = IniVal

        allocate ( FOURDVAR(ng) % CostNorm(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostNorm = IniVal

        allocate ( FOURDVAR(ng) % ObsCost(0:NstateVar(ng)) )
        FOURDVAR(ng) % ObsCost = IniVal

        allocate ( FOURDVAR(ng) % DataPenalty(0:NstateVar(ng)) )
        FOURDVAR(ng) % DataPenalty = IniVal

#   if defined DATALESS_LOOPS || defined TL_W4DPSAS          || \
       defined W4DPSAS        || defined W4DPSAS_SENSITIVITY
        allocate ( FOURDVAR(ng) % NLPenalty(0:NstateVar(ng)) )
        FOURDVAR(ng) % NLPenalty = IniVal
#   endif

        allocate ( FOURDVAR(ng) % NLobsCost(0:NstateVar(ng)) )
        FOURDVAR(ng) % NLobsCost = IniVal

#   if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
        allocate ( FOURDVAR(ng) % bc_ak(Nbico(ng)) )
        FOURDVAR(ng) % bc_ak = IniVal

        allocate ( FOURDVAR(ng) % bc_bk(Nbico(ng)) )
        FOURDVAR(ng) % bc_bk = IniVal

        allocate ( FOURDVAR(ng) % zdf1(Nbico(ng)) )
        FOURDVAR(ng) % zdf1 = IniVal

        allocate ( FOURDVAR(ng) % zdf2(Nbico(ng)) )
        FOURDVAR(ng) % zdf2 = IniVal

        allocate ( FOURDVAR(ng) % zdf3(Nbico(ng)) )
        FOURDVAR(ng) % zdf3 = IniVal

        allocate ( FOURDVAR(ng) % pc_r2d(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % pc_r2d = IniVal

        allocate ( FOURDVAR(ng) % rhs_r2d(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % rhs_r2d = IniVal

        allocate ( FOURDVAR(ng) % r2d_ref(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % r2d_ref = IniVal

        allocate ( FOURDVAR(ng) % zeta_ref(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % zeta_ref = IniVal

        allocate ( FOURDVAR(ng) % p_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
        FOURDVAR(ng) % p_r2d = IniVal

        allocate ( FOURDVAR(ng) % r_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
        FOURDVAR(ng) % r_r2d = IniVal

        allocate ( FOURDVAR(ng) % bp_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
        FOURDVAR(ng) % bp_r2d = IniVal

        allocate ( FOURDVAR(ng) % br_r2d(LBi:UBi,LBj:UBj,Nbico(ng)) )
        FOURDVAR(ng) % br_r2d = IniVal

        allocate ( FOURDVAR(ng) % tl_rhs_r2d(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % tl_rhs_r2d = IniVal

        allocate ( FOURDVAR(ng) % tl_r2d_ref(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % tl_r2d_ref = IniVal

        allocate ( FOURDVAR(ng) % tl_zeta_ref(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % tl_zeta_ref = IniVal

        allocate ( FOURDVAR(ng) % ad_rhs_r2d(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % ad_rhs_r2d = IniVal

        allocate ( FOURDVAR(ng) % ad_r2d_ref(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % ad_r2d_ref = IniVal

        allocate ( FOURDVAR(ng) % ad_zeta_ref(LBi:UBi,LBj:UBj) )
        FOURDVAR(ng) % ad_zeta_ref = IniVal
#   endif
#  endif

        allocate ( FOURDVAR(ng) % ObsVar(NstateVar(ng)) )
        FOURDVAR(ng) % ObsVar = 1.0_r8

        allocate ( FOURDVAR(ng) % ObsCount(0:NstateVar(ng)) )
        FOURDVAR(ng) % ObsCount = 0

        allocate ( FOURDVAR(ng) % ObsReject(0:NstateVar(ng)) )
        FOURDVAR(ng) % ObsReject = 0

#  ifdef RPCG
        allocate ( FOURDVAR(ng) % cg_pxsave(Ndatum(ng)+1) )
        FOURDVAR(ng) % cg_pxsave = 0

        allocate ( FOURDVAR(ng) % cg_pxsum(Ndatum(ng)+1) )
        FOURDVAR(ng) % cg_pxsum = 0

        allocate ( FOURDVAR(ng) % cg_Dpxsum(Ndatum(ng)+1) )
        FOURDVAR(ng) % cg_Dpxsum = 0
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in the time of each observation survey.
!
        CALL netcdf_get_fvar (ng, iNLM, TRIM(OBS(ng)%name),             &
     &                        TRIM(Vname(1,idOday)),                    &
     &                        FOURDVAR(ng)%SurveyTime,                  &
     &                        start = (/1/),                            &
     &                        total = (/Nsurvey(ng)/))
        IF (exit_flag.ne.NoError) RETURN
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
      allocate ( ObsType(Mobs) )
      ObsType = 0

      allocate ( ObsProv(Mobs) )
      ObsProv = 0

      allocate ( ObsErr(Mobs) )
      ObsErr = IniVal

      allocate ( ObsScale(Mobs) )
      ObsScale = IniVal

      allocate ( ObsVal(Mobs) )
      ObsVal = IniVal

      allocate ( Tobs(Mobs) )
      Tobs = IniVal

      allocate ( Xobs(Mobs) )
      Xobs = IniVal

      allocate ( Yobs(Mobs) )
      Yobs = IniVal

#  ifdef SOLVE3D
      allocate ( Zobs(Mobs) )
      Zobs = IniVal
#  endif

      allocate ( ADmodVal(Mobs) )
      ADmodVal = IniVal

      allocate ( NLmodVal(Mobs) )
      NLmodVal = IniVal

#  ifdef BGQC
      allocate ( BgErr(Mobs) )
      BgErr = IniVal

      allocate ( BgThresh(Mobs) )
      BgThresh = IniVal

      allocate ( Jact_S(0:Ninner) )
      Jact_S = IniVal
#  endif

#  ifdef TLM_OBS
      allocate ( TLmodVal(Mobs) )
      TLmodVal = IniVal
#   if defined SENSITIVITY_4DVAR    || defined TL_W4DPSAS  || \
       defined TL_W4DVAR            || defined W4DPSAS     || \
       defined W4DVAR               || defined ARRAY_MODES || \
       defined CLIPPING
      allocate ( TLmodVal_S(Mobs,Ninner,Nouter) )
      TLmodVal_S = IniVal
#   endif
#   if defined TL_W4DPSAS           || defined W4DPSAS    || \
       defined W4DPSAS_SENSITIVITY
      allocate ( BCKmodVal(Mobs) )
      BCKmodVal = IniVal
      allocate ( Hbk(Mobs,Nouter) )
      Hbk = IniVal
#   endif
#  endif

#  ifdef WEAK_CONSTRAINT
!
!  Allocate and initialize weak constraint conjugate gradient vectors.
!
#   if defined TL_W4DPSAS          || defined TL_W4DVAR          || \
       defined W4DPSAS_SENSITIVITY || defined W4DVAR_SENSITIVITY
      allocate ( ad_zcglwk(Mobs,Ninner+1) )
      ad_zcglwk = IniVal

      allocate ( ad_zgrad0(Mobs) )
      ad_zgrad0 = IniVal

      allocate ( ad_cg_innov(Mobs) )
      ad_cg_innov = IniVal

      allocate ( ad_cg_pxsave(Mobs) )
      ad_cg_pxsave = IniVal

      allocate ( tl_zcglwk(Mobs,Ninner+1) )
      tl_zcglwk = IniVal

      allocate ( tl_zgrad0(Mobs) )
      tl_zgrad0 = IniVal

      allocate ( tl_cg_innov(Mobs) )
      tl_cg_innov = IniVal

      allocate ( tl_cg_pxsave(Mobs) )
      tl_cg_pxsave = IniVal

      allocate ( tl_ObsVal(Mobs) )
      tl_ObsVal = IniVal

      allocate ( ad_ObsVal(Mobs) )
      ad_ObsVal = IniVal
#   endif
      allocate ( cg_dla(Ninner,Nouter) )
      cg_dla = IniVal

#   ifdef RPCG
      allocate ( Jb0(0:Nouter) )
      Jb0 = IniVal

      allocate ( Hdx(Mobs) )
      Hdx = IniVal

#    ifdef W4DPSAS_SENSITIVITY

      allocate ( zcglwk(Mobs,Ninner+1,Nouter) )
      zcglwk = IniVal

      allocate ( vcglwk(Mobs,Ninner+1,Nouter) )
      vcglwk = IniVal

      allocate ( vcglev(Mobs,Ninner,Nouter) )
      vcglev = IniVal

      allocate ( zgrad0(Mobs,Nouter) )
      zgrad0 = IniVal

#    else
      allocate ( zcglwk(Mobs+1,Ninner+1,Nouter) )
      zcglwk = IniVal

      allocate ( vcglwk(Mobs+1,Ninner+1,Nouter) )
      vcglwk = IniVal

      allocate ( vcglev(Mobs+1,Ninner,Nouter) )
      vcglev = IniVal

      allocate ( zgrad0(Mobs+1,Nouter) )
      zgrad0 = IniVal
#    endif

      allocate ( vgrad0(Mobs+1) )
      vgrad0 = IniVal

      allocate ( vgrad0s(Mobs+1) )
      vgrad0s = IniVal

      allocate ( cg_innov(Mobs+1) )
      cg_innov = IniVal

#     ifndef RPCG
      allocate ( cg_pxsave(Mobs) )
      cg_pxsave = IniVal
#     endif

#    else

      allocate ( zcglwk(Mobs,Ninner+1,Nouter) )
      zcglwk = IniVal

      allocate ( vcglev(Mobs,Ninner,Nouter) )
      vcglev = IniVal

      allocate ( zgrad0(Mobs,Nouter) )
      zgrad0 = IniVal

      allocate ( vgrad0(Mobs) )
      vgrad0 = IniVal

      allocate ( vgrad0s(Mobs) )
      vgrad0s = IniVal

      allocate ( gdgw(Mobs) )
      gdgw = IniVal

      allocate ( vgnorm(Nouter) )
      vgnorm = IniVal

      allocate ( cg_innov(Mobs) )
      cg_innov = IniVal

      allocate ( cg_pxsave(Mobs) )
      cg_pxsave = IniVal
#    endif
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
      allocate ( NHsteps(2,MstateVar) )
      NHsteps = 0

      allocate ( NVsteps(2,MstateVar) )
      NVsteps = 0

      allocate ( DTsizeH(2,MstateVar) )
      DTsizeH = IniVal

      allocate ( DTsizeV(2,MstateVar) )
      DTsizeV = IniVal

      allocate ( NVstepsB(4,MstateVar) )
      NVstepsB = 0

      allocate ( NHstepsB(4,MstateVar) )
      NHstepsB = 0

      allocate ( DTsizeHB(4,MstateVar) )
      DTsizeHB = IniVal

      allocate ( DTsizeVB(4,MstateVar) )
      DTsizeVB = IniVal


# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Allocate conjugate gradient variables.
!
#  if defined SENSITIVITY_4DVAR || defined TL_W4DPSAS  || \
      defined TL_W4DVAR
      allocate ( ad_cg_beta(Ninner+1) )
      ad_cg_beta = IniVal

      allocate ( ad_cg_delta(Ninner) )
      ad_cg_delta = IniVal

      allocate ( ad_cg_QG(Ninner+1) )
      ad_cg_QG = IniVal

      allocate ( ad_cg_Gnorm(Nouter) )
      ad_cg_Gnorm = IniVal

      allocate ( tl_cg_beta(Ninner+1) )
      tl_cg_beta = IniVal

      allocate ( tl_cg_delta(Ninner) )
      tl_cg_delta = IniVal

      allocate ( tl_cg_QG(Ninner+1) )
      tl_cg_QG = IniVal

      allocate ( tl_cg_Gnorm(Nouter) )
      tl_cg_Gnorm = IniVal
#  endif

      allocate ( cg_alpha(0:Ninner,Nouter) )
      cg_alpha = IniVal

      allocate ( cg_beta(Ninner+1,Nouter) )
      cg_beta = IniVal

      allocate ( cg_tau(0:Ninner,Nouter) )
      cg_tau = IniVal

      allocate ( Ritz(Ninner+1) )
      Ritz = IniVal

# if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
     defined POSTERIOR_ERROR_F

      allocate ( ae_beta(NpostI+1,Nouter) )
      ae_beta = IniVal

      allocate ( ae_zv(NpostI,NpostI) )
      ae_zv = IniVal

      allocate ( ae_delta(NpostI,Nouter) )
      ae_delta = IniVal

      allocate ( ae_trace(NpostI+1) )
      ae_trace = IniVal

      allocate ( ae_Gnorm(Nouter) )
      ae_Gnorm = IniVal

      allocate ( ae_Tmatrix(NpostI,3) )
      ae_Tmatrix = IniVal

      allocate ( ae_Ritz(NpostI,Nouter) )
      ae_Ritz = IniVal

      allocate ( ae_RitzErr(NpostI,Nouter) )
      ae_RitzErr = IniVal

#  if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      allocate ( zLanczos_coef(Ninner,Ninner) )
      zLanczos_coef = IniVal

      allocate ( zLanczos_inv(Ninner,Ninner) )
      zLanczos_inv = IniVal

      allocate ( zLanczos_err(Ninner,Ninner) )
      zLanczos_err = IniVal
#  endif
# endif

#  ifdef WEAK_CONSTRAINT
      allocate ( cg_zv(Ninner,Ninner,Nouter) )
#  else
      allocate ( cg_zv(Ninner,Ninner) )
#  endif
      cg_zv = IniVal

      allocate ( cg_delta(Ninner,Nouter) )
      cg_delta = IniVal

      allocate ( cg_gamma(Ninner,Nouter) )
      cg_gamma = IniVal

      allocate ( cg_Gnorm(Nouter) )
      cg_Gnorm = IniVal

      allocate ( cg_Gnorm_v(Nouter) )
      cg_Gnorm_v = IniVal

      allocate ( cg_Gnorm_y(Nouter) )
      cg_Gnorm_y = IniVal

      allocate ( cg_Greduc(Ninner,Nouter) )
      cg_Greduc = IniVal

      allocate ( cg_QG(Ninner+1,Nouter) )
      cg_QG = IniVal

      allocate ( cg_Tmatrix(Ninner,3) )
      cg_Tmatrix = IniVal

      allocate ( cg_Ritz(Ninner,Nouter) )
      cg_Ritz = IniVal

      allocate ( cg_RitzErr(Ninner,Nouter) )
      cg_RitzErr = IniVal

      allocate ( cg_zu(Ninner,Nouter) )
      cg_zu = IniVal

# endif
# if defined IS4DVAR_SENSITIVITY
!
!  Allocate Lanczos algorithm coefficients.
!
      allocate ( cg_beta(Ninner+1,Nouter) )
      cg_beta = IniVal

      allocate ( cg_delta(Ninner,Nouter) )
      cg_delta = IniVal

      allocate ( cg_zu(Ninner,Nouter) )
      cg_zu = IniVal
# endif
# if defined HESSIAN_SV || defined HESSIAN_SO || defined HESSIAN_FSV
!
!  Allocate Lanczos algorithm coefficients.
!
      allocate ( cg_beta(Ninner+1,Nouter) )
      cg_beta = IniVal

      allocate ( cg_delta(Ninner,Nouter) )
      cg_delta = IniVal

      allocate ( zLanczos_diag(Ninner) )
      zLanczos_diag = IniVal

      allocate ( zLanczos_offdiag(Ninner) )
      zLanczos_offdiag = IniVal
#  ifdef LCZ_FINAL

      allocate ( GSmatrix(Ninner,Ninner) )
      GSmatrix = IniVal

      allocate ( GSmatinv(Ninner,Ninner) )
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
        wrote_Zobs(ng)=.FALSE.
        wrtMisfit(ng)=.FALSE.
        wrtNLmod(ng)=.FALSE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
# ifdef BALANCE_OPERATOR
        wrtZetaRef(ng)=.FALSE.
# endif
# ifdef IS4DVAR_SENSITIVITY
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

# ifdef OBSERVATIONS
!
!  Read in observations global variance.
!
        CALL netcdf_get_fvar (ng, iNLM, TRIM(OBS(ng)%name),             &
     &                        TRIM(Vname(1,idOvar)),                    &
     &                        FOURDVAR(ng)%ObsVar,                      &
     &                        start = (/1/),                            &
     &                        total = (/my_NstateVar/))
        IF (exit_flag.ne.NoError) RETURN
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
