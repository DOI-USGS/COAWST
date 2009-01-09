#include "cppdefs.h"
      MODULE mod_fourdvar

#if defined FOUR_DVAR || defined VERIFICATION
!
!svn $Id: mod_fourdvar.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!  CostGrad     Cost function gradient norm.                           !
!  CostGradDot  Adjoint solution (gradient) dot product used to adjust !
!                 cost function values in the cgradient algorithm.     !
!  CostGradOld  Cost function gradient norm from previous iteration.   !
!  CostNorm     Total cost function normalization scales (minimization !
!                 starting value, Jb+Jo).                              !
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
!  Lprecond     Logical switch to pre-conditioning with Hessian        !
!                 eigenvectors.                                        !
!  nConvRitz    Number of converged Ritz eigenvalues.                  !
!  NLmodVal     Nonlinear model values at observation locations.       !
!  NritzEV      If preconditioning, number of eigenpairs to use.       !
!  NHsteps      Full number of horizontal diffusion steps for spatial  !
!                 convolutions.                                        !
!  NVsteps      Full number of vertical diffusion steps for spatial    !
!                 convolutions.                                        !
!  ObsCost      Current observation cost function misfit (mean squared !
!                 difference) between model and observations, Jo.      !
!  ObsCount     Current observation counts per state variable.         !
!  ObsErr       Observation error.                                     !
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
!                 coefficients available to process in MODname file.   !
!  haveNLmod    Logical switch indicating that there is nonlinear      !
!                 model data available to process in MODname file.     !
!  haveTLmod    Logical switch indicating that there is tangent        !
!                 model data available to process in MODname file.     !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

# ifdef OBSERVATIONS
        integer,  allocatable :: ObsType(:)

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
#  ifdef S4DVAR
        real(r8), allocatable :: NLmodVal(:,:)
#  else
        real(r8), allocatable :: NLmodVal(:)
#  endif
#  ifdef TLM_OBS
#   ifdef IS4DVAR_OLD
        real(r8), allocatable :: TLmodVal(:,:)
#   else
        real(r8), allocatable :: TLmodVal(:)
#   endif
#  endif
#  if defined WEAK_CONSTRAINT || defined IOM
        real(r8), allocatable :: zdelta(:)
        real(r8), allocatable :: zbeta(:)
        real(r8), allocatable :: zcglwk(:,:)
        real(r8), allocatable :: vcglev(:,:)
        real(r8), allocatable :: zgrad0(:)
        real(r8), allocatable :: zqg0(:)
        real(r8) :: zgnorm = 0.0_r8
#  endif
!
        TYPE T_FOURDVAR

          integer , pointer :: NobsSurvey(:)
          integer , pointer :: ObsCount(:)
          integer , pointer :: ObsReject(:)
#  ifdef FOUR_DVAR
          real(r8), pointer :: BackCost(:)
          real(r8), pointer :: CostFun(:)
          real(r8), pointer :: CostFunOld(:)
          real(r8), pointer :: CostGrad(:)
          real(r8), pointer :: CostGradOld(:)
          real(r8), pointer :: CostGradDot(:)
          real(r8), pointer :: CostNorm(:)
          real(r8), pointer :: ObsCost(:)
#  endif
          real(r8), pointer :: ObsVar(:)
          real(r8), pointer :: SurveyTime(:)
  
        END TYPE T_FOURDVAR

        TYPE (T_FOURDVAR), allocatable :: FOURDVAR(:)
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
        logical, dimension(Ngrids) :: Load_Zobs
        logical, dimension(Ngrids) :: wrote_Zobs
!
!  Maximum number of observations to process.
!
        integer :: Mobs
!
!  Number of model state variables to process.
!
        integer, dimension(Ngrids) :: NstateVar
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
        integer, dimension(Ngrids) :: Ndatum         
!
!  Number of observations surveys available.
!
        integer, dimension(Ngrids) :: Nsurvey
!
!  Observation surveys counter.
!
        integer, dimension(Ngrids) :: ObsSurvey
!
!  Current number of observations processed.
!
        integer, dimension(Ngrids) :: Nobs
!
!  Current starting and ending observation file index processed.
!
        integer, dimension(Ngrids) :: NstrObs
        integer, dimension(Ngrids) :: NendObs
!
!  Background error covariance normalization method:
!
!       [0] Exact, very expensive
!       [1] Approximated, randomization
!
        integer, dimension(Ngrids) :: Nmethod
!
!  Random number generation scheme for randomization:
!
!       [0] Intrinsic F90 routine "randon_number"
!       [1] Gaussian distributed deviates, numerical recipes
!
        integer, dimension(Ngrids) :: Rscheme
!
!  Number of iterations in the randomization ensemble used to
!  compute the background error covariance, B, normalization
!  factors. This factors ensure that the diagonal elements of
!  B are equal to unity (Fisher and Coutier, 1995).
!
        integer :: Nrandom = 1000
!
!  Optimal, normalized cost funtion minimum.  If the statistical
!  hypotheses between the background and observations errors
!  is consistemt, the cost function value at the minimum, Jmin,
!  is idealy equal to half the number of observations assimilated
!  (Optimality=1=2*Jmin/Nobs), for a linear system.
!
        real(r8), dimension(Ngrids) :: Optimality
!
!  Switch to activate the processing of model state at the observation
!  locations.
!
        logical, dimension(Ngrids) :: ProcessObs
!
!  Switch to activate writting of nonlinear model values at
!  observations locations into observations NetCDF file.
!
        logical, dimension(Ngrids) :: wrtNLmod
!
!  Switch to activate writting of representer model values at
!  observations locations into observations NetCDF file.
!
        logical, dimension(Ngrids) :: wrtRPmod
!
!  Switch to activate writting of tangent linear model values at
!  observations locations into observations NetCDF file.
!
        logical, dimension(Ngrids) :: wrtTLmod
!
!  Switch to activate writting of initial and final model-observation
!  misfit (innovation) vector into 4DVAR output NetCDF file.
!
        logical, dimension(Ngrids) :: wrtMisfit
!
!  Swiches indicating that there is representer coeffiecients,
!  nonlinear, and tangent linear model data available to process
!  in 4DVAR NetCDF output file (MODname).  At the beeginning,
!  there is not data since this file has been just defined.
!
        logical, dimension(Ngrids) :: haveADmod
        logical, dimension(Ngrids) :: haveNLmod
        logical, dimension(Ngrids) :: haveTLmod
!
!-----------------------------------------------------------------------
!  Spatial convolutions parameters
!-----------------------------------------------------------------------
!
!  Full number of horizontal and vertical diffusion operator step for
!  spatial convolutions.
!
        integer, allocatable :: NHsteps(:)
        integer, allocatable :: NVsteps(:)
!
!  Horizontal and vertical diffusion operator time-step size for
!  spatial convolutions.
!
        real(r8), allocatable :: DTsizeH(:)
        real(r8), allocatable :: DTsizeV(:)
!
!  Minimum and maximum Horizontal and vertical diffusion coefficients
!  used in spatial convolutions.
!
        real(r8), dimension(Ngrids) :: KhMin
        real(r8), dimension(Ngrids) :: KhMax
        real(r8), dimension(Ngrids) :: KvMin
        real(r8), dimension(Ngrids) :: KvMax

# if defined IS4DVAR
!
!-----------------------------------------------------------------------
!  Conjugate gradient parameters.
!-----------------------------------------------------------------------
!
!  Conjugate gradient coefficients.
!
        real(r8), allocatable :: cg_alpha(:,:)
        real(r8), allocatable :: cg_beta(:,:)
        real(r8), allocatable :: cg_tau(:,:)
!
!  Number of converged Ritz eigenvalues.
!
        integer :: nConvRitz
!
!  Ritz values maximum error limit.
!
        real(r8) :: RitzMaxErr
!
!  Converged Ritz eigenvalues used to approximate Hessian matrix during
!  preconditioning.
!
        real(r8), allocatable :: Ritz(:)

#  if defined LANCZOS
!
!  Lanczos algorithm coefficients.
!
        real(r8), allocatable :: cg_delta(:,:)
        real(r8), allocatable :: cg_gamma(:,:)
!
!  Initial gradient vector normalization factor.
!  
        real(r8), allocatable :: cg_Gnorm(:)
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
!
!  Orthogonal eigenvectors of Lanczos recurrence term Q(k)T(k).
!
        real(r8), allocatable :: cg_zv(:,:)
#  endif
# endif
# if defined OBS_SENSITIVITY
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
!  Switch to pre-conditioning with Hessian eigenpairs.
!
        logical :: Lprecond

!  Number of iteration between each steepest descent step.
!
        integer :: NiterSD
!
!  Iteration of the last steepest descent.
!
        integer :: IterSD
!
!  Conjugate gradient algorithm: [0] Fletcher-Reeves [1] Polak-Riviere
!
        integer :: ICG
!
!  Pass number into the descent algorithm within the same iteration.
!
        integer :: Ipass = 0
!
!  If preconditioning, number of Ritz eigenpairs to use.
!
        integer :: NritzEV
!
!  Conjugate gradient alogorithm step size.
!
        real(r8) :: CGstepI              ! initial step size
        real(r8) :: CGstepF              ! first guess step size
        real(r8) :: CGstepR              ! refined step size
!
!  Step size contributions from background and observations
!
        real(r8) :: StepBotBck = 0.0_r8  ! background denominator
        real(r8) :: StepBotObs = 0.0_r8  ! observation denominator
        real(r8) :: StepTopBck = 0.0_r8  ! background numerator
        real(r8) :: StepTopObs = 0.0_r8  ! observation numerator
!
!  Strong constraint conjugate gradient orthogonality tolerance: dot
!  product orthogonality between current and previous search.
!
        real(r8) :: CGtol
!
!  Weak constraint conjugate gradient algorithm convergence criteria.
!
        real(r8) :: CGeps
!
!  Weak contraint conjugate gradient norms.
!
        real(r8) :: cg_gammam1 = 0.0_r8
        real(r8) :: cg_sigmam1 = 0.0_r8
        real(r8) :: cg_rnorm = 0.0_r8
!
!  Cost function percentage change between succesive iterations for
!  convergence.
!
        real(r8) :: CostFunFac
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
!
!  Brunt-Vaisala expansion polynomial used for cost function
!  normalization.
!
        integer :: npN2
        real(r8), allocatable :: pcoef_N2(:)

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

      CONTAINS

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
      integer :: Nritz, i, ng

# ifdef OBSERVATIONS
      real(r8), parameter :: IniVal = 0.0_r8
# endif

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
        CALL netcdf_get_dim(ng, iNLM, TRIM(OBSname(ng)))
        IF (exit_flag.ne.NoError) RETURN
!
        Ndatum(ng)=0
        Nsurvey(ng)=0
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'datum') then
            Ndatum(ng)=dim_size(i)
          ELSE IF (TRIM(dim_name(i)).eq.'survey') THEN
            Nsurvey(ng)=dim_size(i)
          END IF
        END DO
        IF (Ndatum(ng).eq.0) THEN
          WRITE (stdout,10) 'datum', TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (Nsurvey(ng).eq.0) THEN
          WRITE (stdout,10) 'survey', TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Allocate module structure.
!-----------------------------------------------------------------------
!
      allocate ( FOURDVAR(Ngrids) )

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

        allocate ( FOURDVAR(ng) % NobsSurvey(Nsurvey(ng)) )
        FOURDVAR(ng) % NobsSurvey = 0

        allocate ( FOURDVAR(ng) % SurveyTime(Nsurvey(ng)) )
        FOURDVAR(ng) % SurveyTime = IniVal

#  ifdef FOUR_DVAR
        allocate ( FOURDVAR(ng) % BackCost(0:NstateVar(ng)) )
        FOURDVAR(ng) % BackCost = IniVal

        allocate ( FOURDVAR(ng) % CostFun(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostFun = IniVal

        allocate ( FOURDVAR(ng) % CostFunOld(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostFunOld = IniVal

        allocate ( FOURDVAR(ng) % CostGrad(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostGrad = IniVal

        allocate ( FOURDVAR(ng) % CostGradOld(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostGradOld = IniVal

        allocate ( FOURDVAR(ng) % CostGradDot(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostGradDot = IniVal

        allocate ( FOURDVAR(ng) % CostNorm(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostNorm = IniVal

        allocate ( FOURDVAR(ng) % ObsCost(0:NstateVar(ng)) )
        FOURDVAR(ng) % ObsCost = IniVal
#  endif
        allocate ( FOURDVAR(ng) % ObsVar(NstateVar(ng)) )
        FOURDVAR(ng) % ObsVar = 1.0_r8

        allocate ( FOURDVAR(ng) % ObsCount(0:NstateVar(ng)) )
        FOURDVAR(ng) % ObsCount = 0

        allocate ( FOURDVAR(ng) % ObsReject(0:NstateVar(ng)) )
        FOURDVAR(ng) % ObsReject = 0

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
        CALL netcdf_get_ivar (ng, iNLM, TRIM(OBSname(ng)),              &
     &                        TRIM(Vname(1,idNobs)),                    &
     &                        FOURDVAR(ng)%NobsSurvey)
        IF (exit_flag.ne.NoError) RETURN
!
!  Read in the time of each observation survey.
!
        CALL netcdf_get_fvar (ng, iNLM, TRIM(OBSname(ng)),              &
     &                        TRIM(Vname(1,idOday)),                    &
     &                        FOURDVAR(ng)%SurveyTime)
        IF (exit_flag.ne.NoError) RETURN
!
!  Determine maximum size of observation arrays.
!
#  if defined WEAK_CONSTRAINT || defined IOM
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

#  ifdef S4DVAR
      allocate ( NLmodVal(Mobs,2) )
#  else
      allocate ( NLmodVal(Mobs) )
#  endif
      NLmodVal = IniVal

#  ifdef TLM_OBS
#   ifdef IS4DVAR_OLD
      allocate ( TLmodVal(Mobs,2) )
#   else
      allocate ( TLmodVal(Mobs) )
#   endif
      TLmodVal = IniVal
#  endif

#  if defined WEAK_CONSTRAINT || defined IOM
!
!  Allocate and initialize weak constraint conjugate gradient vectors.
!
      allocate ( zdelta(Ninner) )
      zdelta = IniVal

      allocate ( zbeta(Ninner+1) )
      zbeta = IniVal

      allocate ( zcglwk(Mobs,Ninner+1) )
      zcglwk = IniVal

      allocate ( vcglev(Mobs,Ninner) )
      vcglev = IniVal

      allocate ( zgrad0(Mobs) )
      zgrad0 = IniVal

      allocate ( zqg0(Ninner+1) )
      zqg0 = IniVal
#  endif
# else
!
!-----------------------------------------------------------------------
!  Allocate module variables.
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
# ifdef IS4DVAR
!
!  Determine the number of Hessian eigenpairs available for
!  preconditioning.
!
      Nritz=Ninner+1
#  if !defined LANCZOS
      IF (Lprecond) THEN
        DO ng=1,Ngrids
          CALL netcdf_get_dim (ng, iNLM, TRIM(HSSname(ng)),             &
     &                         TRIM(Vname(1,idtime)),                   &
     &                         Nritz)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF
#  endif
#  ifndef LANCZOS
!
!  Set number of eigenpairs to use to available quantity.
!
      IF (Lprecond) THEN
        IF (NritzEV.eq.0) THEN
          NritzEV=Nritz
        END IF
      END IF
#  endif
# endif
!
!  Allocate convolution parameters.
!
      allocate ( NHsteps(MstateVar) )
      allocate ( NVsteps(MstateVar) )
      allocate ( DTsizeH(MstateVar) )
      allocate ( DTsizeV(MstateVar) )

# ifdef IS4DVAR
!
!  Allocate conjugate gradient variables.
!
      allocate ( cg_alpha(0:Ninner,Nouter) )
      cg_alpha = IniVal

#  ifdef LANCZOS
      allocate ( cg_beta(Ninner+1,Nouter) )
      cg_beta = IniVal
#  else
      allocate ( cg_beta(0:Ninner,Nouter) )
      cg_beta = IniVal
#  endif

      allocate ( cg_tau(0:Ninner,Nouter) )
      cg_tau = IniVal

      allocate ( Ritz(Nritz) )
      Ritz = IniVal

#  ifdef LANCZOS
      allocate ( cg_delta(Ninner,Nouter) )
      cg_delta = IniVal

      allocate ( cg_gamma(Ninner,Nouter) )
      cg_gamma = IniVal

      allocate ( cg_Gnorm(Nouter) )
      cg_Gnorm = IniVal

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

      allocate ( cg_zv(Ninner,Ninner) )
      cg_zv = IniVal

#  endif
# endif
# if defined OBS_SENSITIVITY
!
!  Allocate Lanczos algorithm coefficients.
!
      allocate ( cg_beta(Ninner+1,Nouter) )
      cg_beta = IniVal

      allocate ( cg_delta(Ninner,Nouter) )
      cg_delta = IniVal
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
        KhMin(ng)=1.0_r8
        KhMax(ng)=1.0_r8
        KvMin(ng)=1.0_r8
        KvMax(ng)=1.0_r8
        Optimality(ng)=0.0_r8

# ifdef OBSERVATIONS
!
!  Read in observations global variance.
!
        CALL netcdf_get_fvar (ng, iNLM, TRIM(OBSname(ng)),              &
     &                        TRIM(Vname(1,idOvar)),                    &
     &                        FOURDVAR(ng)%ObsVar)
        IF (exit_flag.ne.NoError) RETURN
# endif
      END DO

# ifdef OBSERVATIONS
!
 10   FORMAT (/,' MOD_FOURDVAR - error inquiring dimension: ',a,2x,      &
     &          ' in input NetCDF file: ',a)
# endif

      RETURN
      END SUBROUTINE initialize_fourdvar
#endif
      END MODULE mod_fourdvar
