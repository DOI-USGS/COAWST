      MODULE ocean_control_mod
!
!svn $Id: w4dpsas_ocean.h 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Weak Constraint 4-Dimensional Variational Data            !
!          Assimilation Driver: Physical-space Statistical Analysis    !
!          System (4D-PSAS). Dual formulation in observarion space.    !
!                                                                      !
!  This driver is used for strong/weak constraint 4D-Var where errors  !
!  may be considered in both model and observations.                   !
!                                                                      !
!  The routines in this driver control the initialization,  time-      !
!  stepping, and finalization of  ROMS/TOMS  model following ESMF      !
!  conventions:                                                        !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Moore, A.M., H.G. Arango, G. Broquet, B.S. Powell, A.T. Weaver,     !
!    and J. Zavala-Garay, 2011: The Regional Ocean Modeling System     !
!    (ROMS)  4-dimensional variational data assimilations systems,     !
!    Part I - System overview and formulation, Prog. Oceanogr., 91,    !
!    34-49, doi:10.1016/j.pocean.2011.05.004.                          !
!                                                                      !
!  Moore, A.M., H.G. Arango, G. Broquet, C. Edward, M. Veneziani,      !
!    B. Powell, D. Foley, J.D. Doyle, D. Costa, and P. Robinson,       !
!    2011: The Regional Ocean Modeling System (ROMS) 4-dimensional     !
!    variational data assimilations systems, Part II - Performance     !
!    and application to the California Current System, Prog.           !
!    Oceanogr., 91, 50-73, doi:10.1016/j.pocean.2011.05.003.           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS/TOMS state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#ifdef MCT_LIB
!
# ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_ocn2wav_coupling
# endif
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: STDrec, Tindex
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(mpiCOMM)) THEN
        OCN_COMM_WORLD=mpiCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel control switches. These scalars switches are
!  independent from standard input parameters.
!
        CALL initialize_parallel
!
!  Read in model tunable parameters from standard input. Allocate and
!  initialize variables in several modules after the number of nested
!  grids and dimension parameters are known.
!
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
!$OMP PARALLEL
#if defined _OPENMP
      MyThread=my_threadnum()
#elif defined DISTRIBUTE
      MyThread=MyRank
#else
      MyThread=0
#endif
      DO ng=1,Ngrids
        chunk_size=(NtileX(ng)*NtileE(ng)+numthreads-1)/numthreads
        first_tile(ng)=MyThread*chunk_size
        last_tile (ng)=first_tile(ng)+chunk_size-1
      END DO
!$OMP END PARALLEL
!
!  Initialize internal wall clocks. Notice that the timings does not
!  includes processing standard input because several parameters are
!  needed to allocate clock variables.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Process Information:',/)
        END IF
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iNLM, 0)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL mod_arrays (allocate_vars)
!$OMP END PARALLEL
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar

      END IF

#if defined MCT_LIB && (defined AIR_OCEAN || defined WAVES_OCEAN)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef AIR_OCEAN
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Read in standard deviation factors for error covariance.
!-----------------------------------------------------------------------
!
!  Initial conditions standard deviation. They are loaded in Tindex=1
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 6, 6, STD(1,ng)%name, STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Model error standard deviation. They are loaded in Tindex=2
!  of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      IF (NSA.eq.2) THEN
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STD(2,ng)%name, STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Open boundary conditions standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 8, 8, STD(3,ng)%name, STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Surface forcing standard deviation.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 9, 9, STD(4,ng)%name, STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Create 4D-Var analysis file that used as initial conditions for the
!  next data assimilation cycle.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        LdefDAI(ng)=.TRUE.
        CALL def_dai (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS nonlinear, tangent linear and     !
!  adjoint models.                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE convolve_mod, ONLY : convolve
      USE convolve_mod, ONLY : error_covariance
      USE ini_adjust_mod, ONLY : load_TLtoAD
#ifdef ADJUST_BOUNDARY
      USE mod_boundary, ONLY : initialize_boundary
#endif
      USE mod_forces, ONLY : initialize_forces
      USE mod_ocean, ONLY : initialize_ocean
      USE normalization_mod, ONLY : normalization
#if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
    defined POSTERIOR_ERROR_F
      USE posterior_mod, ONLY : posterior
      USE random_ic_mod, ONLY : random_ic
#endif
#if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
      USE posterior_var_mod, ONLY : posterior_var
#endif
#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY: balance_ref, biconj
#endif
#ifdef RPCG
      USE ini_adjust_mod, ONLY : ini_adjust
      USE sum_grad_mod, ONLY : sum_grad
      USE sum_imp_mod, ONLY : sum_imp
      USE comp_Jb0_mod, ONLY : comp_Jb0, aug_oper
#endif
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: Lcgini, Linner, Lposterior, add
#ifdef POSTERIOR_EOFS
      logical :: Ltrace
#endif
      integer :: my_inner, my_outer
      integer :: Lbck, Lini, Rec, Rec1, Rec2
      integer :: i, ng, status, tile
      integer :: Fcount, NRMrec
#ifdef RPCG
      integer :: ADrec, nADrec, nLAST
      integer :: irec, jrec, jrec1, jrec2
      integer :: Rec3, Rec4, Rec5,  LTLM1, LTLM2, LiNL
#endif

      integer, dimension(Ngrids) :: indxSave
      integer, dimension(Ngrids) :: Nrec

      character (len=7 ) :: driver
      character (len=20) :: string
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
!  Initialize relevant parameters.
!
      DO ng=1,Ngrids
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        Lfinp(ng)=1         ! forcing index for input
        Lfout(ng)=1         ! forcing index for output history files
#endif
#ifdef ADJUST_BOUNDARY
        Lbinp(ng)=1         ! boundary index for input
        Lbout(ng)=1         ! boundary index for output history files
#endif
        Lold(ng)=1          ! old minimization time index
        Lnew(ng)=2          ! new minimization time index
      END DO
      Lini=1                ! NLM initial conditions record in INI
      Lbck=2                ! background record in INI
      Rec1=1
      Rec2=2
#ifdef RPCG
      Rec3=3
      Rec4=4
      Rec5=5
      LTLM1=1
      LTLM2=2
#endif
      Nrun=1
      outer=0
      inner=0
      ERstr=1
      ERend=Nouter
      driver='w4dpsas'
!
!-----------------------------------------------------------------------
!  Configure weak constraint 4DVAR algorithm: PSAS Approach.
!-----------------------------------------------------------------------
!
!  Initialize the switch to gather weak constraint forcing.
!
      DO ng=1,Ngrids
        WRTforce(ng)=.FALSE.
      END DO
!
!  Clear nonlinear mixing arrays.
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_mixing (ng, tile, iTLM)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Initialize and set nonlinear model initial conditions.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%Fcount
        RST(ng)%Nrec(Fcount)=0
      END DO

!$OMP PARALLEL
      CALL initial
!$OMP END PARALLEL
      IF (exit_flag.ne.NoError) RETURN
!
!  Save nonlinear initial conditions (currently in time index 1,
!  background) into record "Lbck" of INI(ng)%name NetCDF file. The
!  record "Lbck" becomes the background state record and the record
!  "Lini" becomes current nonlinear initial conditions.
!
      DO ng=1,Ngrids
        INI(ng)%Rindex=1
        Fcount=INI(ng)%Fcount
        INI(ng)%Nrec(Fcount)=1
        CALL wrt_ini (ng, 1)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Set nonlinear output history file as the initial basic state
!  trajectory.
!
      DO ng=1,Ngrids
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%base), outer
      END DO

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
      DO ng=1,Ngrids
        BLK(ng)%name=HIS(ng)%name
      END DO
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Model-error covariance normalization and stardard deviation factors.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Compute or read in the error covariance normalization factors.
!  If computing, write out factors to NetCDF. This is an expensive
!  computation that needs to be computed only once for a particular
!  application grid and decorrelation scales.
!
      DO ng=1,Ngrids
        IF (ANY(LwrtNRM(:,ng))) THEN
          CALL def_norm (ng, iNLM, 1)
          IF (exit_flag.ne.NoError) RETURN

          IF (NSA.eq.2) THEN
            CALL def_norm (ng, iNLM, 2)
          IF (exit_flag.ne.NoError) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
          CALL def_norm (ng, iNLM, 3)
          IF (exit_flag.ne.NoError) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL def_norm (ng, iNLM, 4)
          IF (exit_flag.ne.NoError) RETURN
#endif

!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL normalization (ng, tile, 2)
          END DO
!$OMP END PARALLEL
          LdefNRM(1:4,ng)=.FALSE.
          LwrtNRM(1:4,ng)=.FALSE.
        ELSE
          NRMrec=1
          CALL get_state (ng, 5, 5, NRM(1,ng)%name, NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN

          IF (NSA.eq.2) THEN
            CALL get_state (ng, 5, 5, NRM(2,ng)%name, NRMrec, 2)
            IF (exit_flag.ne.NoError) RETURN
          END IF

#ifdef ADJUST_BOUNDARY
          CALL get_state (ng, 10, 10, NRM(3,ng)%name, NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN
#endif

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          CALL get_state (ng, 11, 11, NRM(4,ng)%name, NRMrec, 1)
          IF (exit_flag.ne.NoError) RETURN
#endif
        END IF
      END DO
!
!  Define tangent linear initial conditions file.
!
      DO ng=1,Ngrids
        LdefITL(ng)=.TRUE.
        CALL tl_def_ini (ng)
        LdefITL(ng)=.FALSE.
        IF (exit_flag.ne.NoError) RETURN
      END DO
#ifdef RPCG
!
!  Initialize all records of the ITL file to zero.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, Rec1, Rec1)
          IF (exit_flag.ne.NoError) RETURN
          CALL tl_wrt_ini (ng, Rec1, Rec2)
          IF (exit_flag.ne.NoError) RETURN
          CALL tl_wrt_ini (ng, Rec1, Rec3)
          IF (exit_flag.ne.NoError) RETURN
          CALL tl_wrt_ini (ng, Rec1, Rec4)
          IF (exit_flag.ne.NoError) RETURN
          CALL tl_wrt_ini (ng, Rec1, Rec5)
          IF (exit_flag.ne.NoError) RETURN
          nADrec=0
          IF (nADJ(ng).lt.ntimes(ng)) THEN
            nLAST=Rec5
            nADrec=2*(1+ntimes(ng)/nADJ(ng))
            DO irec=1,nADrec
              CALL tl_wrt_ini (ng, Rec1, nLAST+irec)
              IF (exit_flag.ne.NoError) RETURN
            END DO
          END IF
        END DO
#endif
!
!  Define impulse forcing NetCDF file.
!
      DO ng=1,Ngrids
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Define output 4DVAR NetCDF file containing all processed data
!  at observation locations.
!
      DO ng=1,Ngrids
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO

#if defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_I || \
    defined POSTERIOR_ERROR_F
!
!  Define output Hessian NetCDF file that will eventually contain
!  the intermediate posterior analysis error covariance matrix
!  fields or its EOFs.
!
      DO ng=1,Ngrids
        LdefHSS(ng)=.TRUE.
        CALL def_hessian (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
#if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
!
!  Define output initial or final full posterior error covariance
!  (diagonal) matrix NetCDF.
!
      DO ng=1,Ngrids
        LdefERR (ng)=.TRUE.
        CALL def_error (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute background state trajectory, X_n-1(t)
!  and the background values at the observation points and times.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      DO ng=1,Ngrids
#ifdef AVERAGES
        LdefAVG(ng)=.TRUE.
        LwrtAVG(ng)=.TRUE.
        WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%base), outer
#endif
#ifdef DIAGNOSTICS
        LdefDIA(ng)=.TRUE.
        LwrtDIA(ng)=.TRUE.
        WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%base), outer
#endif
        wrtMisfit(ng)=.TRUE.
        SporadicImpulse(ng)=.FALSE.
        FrequentImpulse(ng)=.FALSE.
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
        END IF
      END DO

!$OMP PARALLEL
#ifdef SOLVE3D
      CALL main3d (RunInterval)
#else
      CALL main2d (RunInterval)
#endif
!$OMP END PARALLEL
      IF (exit_flag.ne.NoError) RETURN

      DO ng=1,Ngrids
#ifdef AVERAGES
        LdefAVG(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
#endif
#ifdef DIAGNOSTICS
        LdefDIA(ng)=.FALSE.
        LwrtDIA(ng)=.FALSE.
#endif
        wrtNLmod(ng)=.FALSE.
      END DO
!
!  Report data penalty function.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          DO i=0,NstateVar(ng)
            IF (i.eq.0) THEN
              string='Total'
            ELSE
              string=Vname(1,idSvar(i))
            END IF
            IF (FOURDVAR(ng)%NLPenalty(i).ne.0.0_r8) THEN
              WRITE (stdout,30) outer, inner, 'NLM',                    &
     &                          FOURDVAR(ng)%NLPenalty(i),              &
     &                          TRIM(string)
            END IF
          END DO
        END IF
!
!  Write out initial data penalty function to NetCDF file.
!
        SourceFile='w4dpsas_ocean.h, ROMS_run'

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                   &
     &                        'NL_iDataPenalty',                        &
     &                        FOURDVAR(ng)%NLPenalty(0:),               &
     &                        (/1/), (/NstateVar(ng)+1/),               &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Clean penalty array before next run of NL model.
!
        FOURDVAR(ng)%NLPenalty=0.0_r8
      END DO
!
!  Set forward basic state NetCDF ID to nonlinear model trajectory to
!  avoid the inquiring stage.
!
      DO ng=1,Ngrids
        FWD(ng)%ncid=HIS(ng)%ncid
      END DO
!
!-----------------------------------------------------------------------
!  Solve the system (following Courtier, 1997):
!
!              (H M_n B (M_n)' H' + Cobs) * w_n = d_n
!
!              d_n = yo - H * Xb_n
!
!  where M_n is the tangent linear model matrix, Cobs is the
!  observation-error covariance, B is the background error covariance
!  and dx_n=B M' H' w_n is the analysis increment so that Xa=Xb+dx_n.
!  d_n is the misfit between observations (yo) and model (H * Xb_n),
!  and H is the linearized observation operator.
!
!  Here, _n denotes a sequence of outer-loop estimates.
!
!  The system does not need to be solved explicitly by inverting the
!  symmetric matrix, P_n:
!
!              P_n = H M_n B (M_n)' H' + Cobs
!
!  but by computing the action of P_n on any vector PSI, such that
!
!              P_n * PSI =  H M_n B (M_n)' H' * PSI + Cobs * PSI
!
!  The (H M_n B (M_n)' H') matrix is not explicitly computed but
!  evaluated by one integration backward of the adjoint model and
!  one integration forward of the tangent linear model for any
!  forcing vector PSI.
!
!  A preconditioned conjugate gradient algorithm is used to compute
!  an approximation PSI for w_n.
!
!-----------------------------------------------------------------------
!
      OUTER_LOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0
!
!  Set basic state trajectory (X_n-1) file to previous outer loop file
!  (outer-1).
!
        DO ng=1,Ngrids
          WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer-1
        END DO
!
!  Clear tangent linear forcing arrays before entering inner-loop.
!  This is very important since these arrays are non-zero and must
!  be zero when running the tangent linear model.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_forces (ng, tile, iTLM)
#ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#endif
          END DO
!$OMP END PARALLEL
        END DO

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
        IF (balance(isFsur)) THEN
          DO ng=1,Ngrids
            CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
            IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL balance_ref (ng, tile, Lini)
              CALL biconj (ng, tile, iNLM, Lini)
            END DO
!$OMP END PARALLEL
            wrtZetaRef(ng)=.TRUE.
          END DO
        END IF
#endif
!
        INNER_LOOP : DO my_inner=0,Ninner
          inner=my_inner
#ifdef RPCG
          IF (inner.ne.Ninner) THEN
            Linner=.TRUE.
          ELSE
            Linner=.FALSE.
          END IF
!
!  Retrieve TLmodVal and NLmodVal when inner=0 and outer>1 for use as Hbk
!  and BCKmodVal respectively.
!
          IF (inner.eq.0.and.outer.gt.1) THEN
            DO ng=1,Ngrids
              CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,             &
     &                               'TLmodel_value', TLmodVal)
              IF (exit_flag.ne. NoError) RETURN
              CALL netcdf_get_fvar (ng, iTLM, DAV(ng)%name,             &
     &                              'NLmodel_value', NLmodVal,          &
     &                              ncid = DAV(ng)%ncid,                &
     &                              start = (/1,outer-1/),              &
     &                              total = (/Ndatum(ng),1/))
              IF (exit_flag.ne. NoError) RETURN
            END DO
          END IF
!
          IF (inner.eq.0) Lcgini=.TRUE.
          DO ng=1,Ngrids
            CALL rpcg_lanczos (ng, iRPM, outer, inner, Ninner, Lcgini)
          END DO
#else
!
!  Initialize conjugate gradient algorithm depending on hot start or
!  outer loop index.
!
          IF (inner.eq.0) THEN
            Lcgini=.TRUE.
            DO ng=1,Ngrids
              CALL congrad (ng, iRPM, outer, inner, Ninner, Lcgini)
            END DO
          END IF
!
!  If initialization step, skip the inner-loop computations.
!
          Linner=.FALSE.
          IF ((inner.ne.0).or.(Nrun.ne.1)) THEN
            IF (((inner.eq.0).and.LhotStart).or.(inner.ne.0)) THEN
              Linner=.TRUE.
            END IF
          END IF
#endif
!
!  Start inner loop computations.
!
          INNER_COMPUTE : IF (Linner) THEN
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate adjoint model forced with any vector PSI at the observation
!  locations and generate adjoint trajectory, Lambda_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize the adjoint model from rest.
!
            DO ng=1,Ngrids
!$OMP PARALLEL
              CALL ad_initial (ng)
!$OMP END PARALLEL
              IF (exit_flag.ne.NoError) RETURN
              wrtMisfit(ng)=.FALSE.
            END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file only once to avoid opening too many files.
!
            DO ng=1,Ngrids
              WRTforce(ng)=.TRUE.
              IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
              Fcount=ADM(ng)%Fcount
              ADM(ng)%Nrec(Fcount)=0
              ADM(ng)%Rindex=0
            END DO
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
            DO ng=1,Ngrids
              IF (Master) THEN
                WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
              END IF
            END DO

!$OMP PARALLEL
#ifdef SOLVE3D
            CALL ad_main3d (RunInterval)
#else
            CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
            DO ng=1,Ngrids
              CALL ad_wrt_his (ng)
              IF (exit_flag.ne.NoError) RETURN
            END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
            DO ng=1,Ngrids
              WRTforce(ng)=.FALSE.
              CALL ad_wrt_his (ng)
              IF (exit_flag.ne.NoError) RETURN
            END DO
!
!  Convolve adjoint trajectory with error covariances.
!
#ifdef POSTERIOR_ERROR_I
            Lposterior=.TRUE.
#else
            Lposterior=.FALSE.
#endif
#ifdef RPCG
!
!  Set the flag that controls the augmentation of the model error
!  forcing terms. This is ONLY done in the outer-loop so
!  LaugWeak=.FALSE. here.
!
            LaugWeak=.FALSE.
#endif
            CALL error_covariance (iTLM, driver, outer, inner,          &
     &                             Lbck, Lini, Lold, Lnew,              &
     &                             Rec1, Rec2, Lposterior)
            IF (exit_flag.ne.NoError) RETURN
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
            IF (Master) THEN
              WRITE (stdout,40) outer, inner
            END IF
            DO ng=1,Ngrids
              TLF(ng)%Rindex=0
#ifdef DISTRIBUTE
              tile=MyRank
#else
              tile=-1
#endif
              CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
              IF (exit_flag.ne.NoError) RETURN
            END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Initialize tangent linear model from initial impulse which is now
!  stored in file ITL(ng)%name.
!
            DO ng=1,Ngrids
              wrtNLmod(ng)=.FALSE.
              wrtTLmod(ng)=.TRUE.
            END DO
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
            DO ng=1,Ngrids
              IF (FrcRec(ng).gt.3) THEN
                FrequentImpulse(ng)=.TRUE.
              END IF
            END DO
!
!  Initialize tangent linear model from ITL(ng)%name, record Rec1.
!
            DO ng=1,Ngrids
              ITL(ng)%Rindex=Rec1
!$OMP PARALLEL
              CALL tl_initial (ng)
!$OMP END PARALLEL
              IF (exit_flag.ne.NoError) RETURN
            END DO
!
!  Activate switch to write out initial misfit between model and
!  observations.
!
            IF ((outer.eq.1).and.(inner.eq.1)) THEN
              DO ng=1,Ngrids
                wrtMisfit(ng)=.TRUE.
              END DO
            END IF
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute (H M B M' H')_n * PSI at observation
!  points which are used in the conjugate gradient algorithm.
!
            DO ng=1,Ngrids
              IF (Master) THEN
                WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
              END IF
            END DO

!$OMP PARALLEL
#ifdef SOLVE3D
            CALL tl_main3d (RunInterval)
#else
            CALL tl_main2d (RunInterval)
#endif
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN

            DO ng=1,Ngrids
              wrtNLmod(ng)=.FALSE.
              wrtTLmod(ng)=.FALSE.
            END DO

#ifdef POSTERIOR_ERROR_F
!
!  Copy the final time tl_var(Lold) into ad_var(Lold) so that it can be
!  written to the Hessian NetCDF file.
!
            add=.FALSE.
            DO ng=1,Ngrids
!$OMP PARALLEL
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL load_TLtoAD (ng, tile, Lold(ng), Lold(ng), add)
              END DO
!$OMP END PARALLEL
            END DO
!
!  Write evolved tangent solution into hessian netcdf file for use
!  later.
!
            IF (inner.ne.0) THEN
              DO ng=1,Ngrids
                CALL wrt_hessian (ng, Lold(ng), Lold(ng))
                IF (exit_flag.ne.NoERRor) RETURN
              END DO
            END IF
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Use conjugate gradient algorithm to find a better approximation
!  PSI to coefficients Beta_n.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
            Nrun=Nrun+1
#ifndef RPCG
            DO ng=1,Ngrids
              Lcgini=.FALSE.
              CALL congrad (ng, iTLM, outer, inner, Ninner, Lcgini)
              IF (exit_flag.ne.NoError) RETURN
            END DO
#else
            Lcgini=.FALSE.
#endif

          END IF INNER_COMPUTE

        END DO INNER_LOOP
!
!-----------------------------------------------------------------------
!  Once the w_n, have been approximated with sufficient accuracy,
!  compute estimates of Lambda_n and Xhat_n by carrying out one
!  backward intergration of the adjoint model and one forward
!  itegration of the nonlinear model.
!-----------------------------------------------------------------------
!
!  Initialize the adjoint model always from rest.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          CALL ad_initial (ng)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Set adjoint history NetCDF parameters.  Define adjoint history
!  file one to avoid opening to many files.
!
        DO ng=1,Ngrids
          WRTforce(ng)=.TRUE.
          IF (Nrun.gt.1) LdefADJ(ng)=.FALSE.
          Fcount=ADM(ng)%Fcount
          ADM(ng)%Nrec(Fcount)=0
          ADM(ng)%Rindex=0
        END DO
!
!  Time-step adjoint model backwards forced with estimated coefficients,
!  Beta_n.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            WRITE (stdout,20) 'AD', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL ad_main3d (RunInterval)
#else
        CALL ad_main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN
!
!  Write out last weak-constraint forcing (WRTforce is still .TRUE.)
!  record into the adjoint history file.  Note that the weak-constraint
!  forcing is delayed by nADJ time-steps.
!
        DO ng=1,Ngrids
          CALL ad_wrt_his (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Write out adjoint initial condition record into the adjoint
!  history file.
!
        DO ng=1,Ngrids
          WRTforce(ng)=.FALSE.
          CALL ad_wrt_his (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#ifdef RPCG
!
!  Get number of records in adjoint NetCDF.
!  We need to do this here since ADM(ng)%Nrec is reset to zero in
!  error_covariance.
!
          DO ng=1,Ngrids
            Fcount=ADM(ng)%Fcount
            Nrec(ng)=ADM(ng)%Nrec(Fcount)
          END DO
#endif
!
!  Convolve adjoint trajectory with error covariances.
!
        Lposterior=.FALSE.
#ifdef RPCG
!
!  Set the flag that controls the augmentation of the model error
!  forcing terms. This is ONLY done in the outer-loop so
!  LaugWeak=.TRUE. here.
!
        LaugWeak=.TRUE.
        LiNL=outer+1
        CALL error_covariance (iNLM, driver, outer, inner,              &
     &                         LiNL, Lini, Lold, Lnew,                  &
     &                         Rec1, Rec2, Lposterior)
        IF (exit_flag.ne.NoError) RETURN
#else
        CALL error_covariance (iNLM, driver, outer, inner,              &
     &                         Lbck, Lini, Lold, Lnew,                  &
     &                         Rec1, Rec2, Lposterior)
        IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef RPCG
!
!  Augmented solver:
!
!  NOTES: The ITL file contains 5 records -
!         Rec2 = the new TL initial condition
!         Rec3 = the sum of the TL initial conditions
!         Rec4 = B^-1(xb-xk)=sum_j^k-1 G_j lambda_j
!         Rec5 = the augmented correction to Rec2.
!         Rec5+1 to Rec5+nADrec/2 = sum of the TLF forcing increments
!         Rec5+nADrec/2+1 to Rec5+nADrec = B^-1 of the sum of the TLF
!                                         forcing increments.
!  Reset the flag LaugWeak flag.
!
          LaugWeak=.FALSE.

#else
!
!  Write out nonlinear model initial conditions into INIname, record
!  INI(ng)%Rindex.
!
          DO ng=1,Ngrids
            CALL wrt_ini (ng, Lnew(ng))
            IF (exit_flag.ne.NoError) RETURN
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
            CALL wrt_frc_AD (ng, Lold(ng), INI(ng)%Rindex)
            IF (exit_flag.ne.NoError) RETURN
          END DO
#  endif
# endif
#ifdef RPCG
!
! Compute the augmented correction to the adjoint propagator.
! We need to use sum (x(k)-x(k-1)) before it is updated. This is
! in record 3 of the ITL file.
!
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, Rec3, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL aug_oper (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO
!
! Save this augmented piece in record 5 of the ITL file.
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, LTLM2, Rec5)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
! Complete the computation of the TL initial condition by adding the
! contribution from the augmented adjoint propagator.
!
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, Rec2, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, Rec5, LTLM2)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO
!
! Write the final TL increment to Rec2 of the ITL file.
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, LTLM2, Rec2)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
! Now update the non-linear model initial conditions.
!
          LiNL=outer+1
          DO ng=1,Ngrids
            CALL get_state (ng, iNLM, 9, INI(ng)%name, LiNL, Lnew(ng))
            CALL get_state (ng, iADM, 4, ITL(ng)%name, Rec2, LTLM1)
          END DO
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL ini_adjust (ng, tile, LTLM1, Lnew(ng))
            END DO
!$OMP END PARALLEL
          END DO
!
!  Write out nonlinear model initial conditions into INIname, record
!  INI(ng)%Rindex.
!
          DO ng=1,Ngrids
            CALL wrt_ini (ng, Lnew(ng))
            IF (exit_flag.ne.NoError) RETURN
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
            CALL wrt_frc_AD (ng, LTLM1, INI(ng)%Rindex)
            IF (exit_flag.ne.NoError) RETURN
# endif
          END DO
!
! Compute the new B^-1(x(k)-x(k-1)) term.
! Gather the final adjoint solutions increments, sum and
! save in record 4 of the ITL file. Use the tl arrays as temporary
! storage.
!
! First add the augmented piece which is computed from the previous
! sum.
!
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, Rec4, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL aug_oper (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO

!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO
!
          DO ng=1,Ngrids
            ADrec=Nrec(ng)
            CALL get_state (ng, iTLM, 4, ADM(ng)%name, ADrec, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO
!
!  Write the current sum of adjoint solutions into record 4 of the ITL
!  file.
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, LTLM2, Rec4)
            IF (exit_flag.ne.NoError) RETURN
          END DO
#endif
!
!  Convert the current adjoint solution in ADM(ng)%name to impulse
!  forcing. Write out impulse forcing into TLF(ng)%name NetCDF file.
!  To facilitate the forcing to the TLM and RPM, the forcing is
!  processed and written in increasing time coordinates (recall that
!  the adjoint solution in ADM(ng)%name is backwards in time).
!
        IF (Master) THEN
          WRITE (stdout,40) outer, inner
        END IF
        DO ng=1,Ngrids
          TLF(ng)%Rindex=0
#ifdef DISTRIBUTE
          tile=MyRank
#else
          tile=-1
#endif
          CALL wrt_impulse (ng, tile, iADM, ADM(ng)%name)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#ifdef RPCG
!
! Now Compute the augmented corrections to the weak constraint
! forcing terms. The sums to far are in records 6 to
! 6+nADrec/2.
!
        DO i=1,nADrec/2
          irec=i
          jrec=Rec5+i
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, jrec, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL aug_oper (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO
!
! Complete the computation of the TLF forcing term by adding the
! contribution from the augmented adjoint propagator. Specify
! the value of 7 for the model variable since this the special
! case in get_state for reading the impulse forcing.
!
          DO ng=1,Ngrids
!! TEST     CALL get_state (ng, iTLM, 4, TLF(ng)%name, irec, LTLM1)
            CALL get_state (ng, 7, 4, TLF(ng)%name, irec, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
!! TEST       CALL sum_grad (ng, tile, LTLM1, LTLM2)
              CALL sum_imp (ng, tile, LTLM2)
            END DO
!$OMP END PARALLEL
          END DO
!
! Write the final forcing increment to he TLF file.
! Note the original TLF file is overwritten at this point.
!
          DO ng=1,Ngrids
!AMM DO WE NEED TO Reset Rindex HERE????
            TLF(ng)%Rindex=0
#ifdef DISTRIBUTE
            tile=MyRank
#else
            tile=-1
#endif
            CALL wrt_aug_imp (ng, tile, iTLM, LTLM2, i, TLF(ng)%name)
          END DO
!
        END DO
!
!
! Gather the increments from the final inner-loop and
! save in record 3 of the ITL file. The current increment
! is in record 2 and the sum so far is in record 3.
!
        DO ng=1,Ngrids
          CALL get_state (ng, iTLM, 4, ITL(ng)%name, Rec2, LTLM1)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_state (ng, iTLM, 4, ITL(ng)%name, Rec3, LTLM2)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL sum_grad (ng, tile, LTLM1, LTLM2)
          END DO
!$OMP END PARALLEL
        END DO
!
! Write the current sum into record 3 of the ITL file.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, LTLM2, Rec3)
          IF (exit_flag.ne.NoError) RETURN
         END DO
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
    defined ADJUST_BOUNDARY
!
!  Write the current sum of the forcing and boundary increments
!  into the INI file into index INI(ng)%Rindex-1 (i.e. we are
!  overwriting the previous fields. Read the current sum into the
!  adjoint arrays since this is what wrt_frc_AD uses.
!
        DO ng=1,Ngrids
          CALL get_state (ng, iADM, 4, ITL(ng)%name, Rec3, LTLM1)
          IF (exit_flag.ne.NoError) RETURN
        END DO
         DO ng=1,Ngrids
           CALL wrt_frc_AD (ng, LTLM1, INI(ng)%Rindex)
           IF (exit_flag.ne.NoError) RETURN
         END DO
#endif
!
! Gather the model error forcing increments and update the
! sume in records 6 to 6+nADrec/2.
!
        DO i=1,nADrec/2
          irec=i
          jrec=Rec5+i
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, jrec, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
!! TEST     CALL get_state (ng, iTLM, 4, TLF(ng)%name, irec, LTLM2)
            CALL get_state (ng, 7, 4, TLF(ng)%name, irec, LTLM2)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
!!TEST        CALL sum_grad (ng, tile, LTLM1, LTLM2)
              CALL sum_imp (ng, tile, LTLM1)
            END DO
!$OMP END PARALLEL
          END DO
!
! Write the current sum into record jrec of the ITL file.
!
          DO ng=1,Ngrids
!!TEST      CALL tl_wrt_ini (ng, LTLM2, jrec)
            CALL tl_wrt_ini (ng, LTLM1, jrec)
            IF (exit_flag.ne.NoError) RETURN
          END DO

        END DO
!
!  Now compute the background cost function Jb0.
!  First compute the contribution from the increments in the
!  initial conditions, surface forcing, and boundary conditions.
!
        Jb0(outer)=0.0_r8
        DO ng=1,Ngrids
          CALL get_state (ng, iADM, 8, ITL(ng)%name, Rec4, LTLM1)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec3, LTLM1)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL comp_Jb0 (ng, tile, iTLM, outer, LTLM1, LTLM1)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Now compute the contribution of model error terms to the
!  background cost function Jb0.
!
!  NOTE: I THINK WE NEED A BETTER WAY OF COMPUTING THE OBS ERROR
!  CONTRIBUTION TO Jb. THE FOLLOWING CALCULATION IS DANGEROUS BECAUSE
!  IT RELIES ON THE FACT THAT THE SURFACING FORCING AND OBC TL ARRAYS
!  READ FROM jrec2 ARE ZERO. THIS MEANS THAT ONLY THE MODEL STATE
!  CONTRIBUTES TO Jb FOR THE MODEL ERROR TERM, WHICH OF COURSE IS
!  WHAT SHOULD BE THE CASE. THE SURFACE FORCING AND OBC CONTRIBUTIONS
!  ARE COMPUTED IN THE PREVIOUS CALL TO COMP_JB0.
!
        DO irec=1,nADrec/2
          jrec1=Rec5+irec
          jrec2=Rec5+nADrec/2+irec
          DO ng=1,Ngrids
            CALL get_state (ng, iADM, 8, ITL(ng)%name, jrec1, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, jrec2, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
          END DO

          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
                CALL comp_Jb0 (ng, tile, iTLM, outer, LTLM1, LTLM1)
            END DO
!$OMP END PARALLEL
          END DO
        END DO
!
!  Overwrite the TFL netcdf file with the sum of the model error
!  forcing increments - required for the next run of the NLM and
!  TLM.
!
        DO irec=1,nADrec/2
          jrec=Rec5+irec
          DO ng=1,Ngrids
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, jrec, LTLM1)
          END DO
          DO ng=1,Ngrids
!!AMM DO WE NEED TO Reset Rindex HERE????
            TLF(ng)%Rindex=0
#ifdef DISTRIBUTE
            tile=MyRank
#else
            tile=-1
#endif
            CALL wrt_aug_imp (ng, tile, iTLM, LTLM1, irec, TLF(ng)%name)
          END DO
       END DO
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Run nonlinear model and compute a "new estimate" of the state
!  trajectory, X_n(t).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Set new basic state trajectory for next outer loop.
!
        DO ng=1,Ngrids
          LdefHIS(ng)=.TRUE.
          LwrtHIS(ng)=.TRUE.
          wrtNLmod(ng)=.TRUE.
          wrtTLmod(ng)=.FALSE.
          WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%base), outer
        END DO
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
        DO ng=1,Ngrids
          IF (FrcRec(ng).gt.3) THEN
            FrequentImpulse(ng)=.TRUE.
          END IF
        END DO
!
!  Clear tangent arrays and the nonlinear model mixing arrays
!  before running nonlinear model (important).
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
            CALL initialize_forces (ng, tile, iTLM)
#ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#endif
            CALL initialize_mixing (ng, tile, iNLM)
          END DO
!$OMF END PARALLEL
        END DO
!
!  Initialize nonlinear model INI(ng)%name file, record outer+2.
!  Notice that NetCDF record index counter is saved because this
!  counter is used to write initial conditions.
!
        DO ng=1,Ngrids
          indxSave(ng)=INI(ng)%Rindex
          INI(ng)%Rindex=outer+2
          RST(ng)%Rindex=0
          Fcount=RST(ng)%Fcount
          RST(ng)%Nrec(Fcount)=0
        END DO

!$OMP PARALLEL
        CALL initial
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN

        DO ng=1,Ngrids
          INI(ng)%Rindex=indxSave(ng)
        END DO
!
!  Activate switch to write out final misfit between model and
!  observations.
!
        IF (outer.eq.Nouter) THEN
          DO ng=1,Ngrids
            wrtMisfit(ng)=.TRUE.
          END DO
        END IF
!
!  Run nonlinear forced by convolved adjoint trajectory impulses and
!  compute new basic state trajectory X_n.
!
        DO ng=1,Ngrids
#ifdef AVERAGES
          LdefAVG(ng)=.TRUE.
          LwrtAVG(ng)=.TRUE.
          WRITE (AVG(ng)%name,10) TRIM(AVG(ng)%base), outer
#endif
#ifdef DIAGNOSTICS
          LdefDIA(ng)=.TRUE.
          LwrtDIA(ng)=.TRUE.
          WRITE (DIA(ng)%name,10) TRIM(DIA(ng)%base), outer
#endif
          IF (Master) THEN
            WRITE (stdout,20) 'NL', ng, ntstart(ng), ntend(ng)
          END IF
        END DO

!$OMP PARALLEL
#ifdef SOLVE3D
        CALL main3d (RunInterval)
#else
        CALL main2d (RunInterval)
#endif
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN

        DO ng=1,Ngrids
#ifdef AVERAGES
          LdefAVG(ng)=.FALSE.
          LwrtAVG(ng)=.FALSE.
#endif
#ifdef DIAGNOSTICS
          LdefDIA(ng)=.FALSE.
          LwrtDIA(ng)=.FALSE.
#endif
          wrtNLmod(ng)=.FALSE.
          wrtTLmod(ng)=.FALSE.
        END DO
!
!  Report data penalty function.
!
        DO ng=1,Ngrids
          IF (Master) THEN
            DO i=0,NstateVar(ng)
              IF (i.eq.0) THEN
                string='Total'
              ELSE
                string=Vname(1,idSvar(i))
              END IF
              IF (FOURDVAR(ng)%NLPenalty(i).ne.0.0_r8) THEN
                WRITE (stdout,30) outer, inner, 'NLM',                  &
     &                            FOURDVAR(ng)%NLPenalty(i),            &
     &                            TRIM(string)
              END IF
            END DO
          END IF
!
!  Write out final data penalty function to NetCDF file.
!
          SourceFile='w4dpsas_ocean.h, ROMS_run'

          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'NL_fDataPenalty',                      &
     &                          FOURDVAR(ng)%NLPenalty(0:),             &
     &                          (/1,outer/), (/NstateVar(ng)+1,1/),     &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
!
!  Clean penalty array before next run of NL model.
!
          FOURDVAR(ng)%NLPenalty=0.0_r8
        END DO
!
!  Close current forward NetCDF file.
!
        DO ng=1,Ngrids
          CALL netcdf_close (ng, iNLM, FWD(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END DO

#ifdef RPCG
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model again to evaluate G(x(k)-x(k-1)).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
        IF (outer.ne.Nouter) THEN
!
!  Clear tangent linear forcing arrays.
!  This is very important since these arrays are non-zero and must
!  be zero when running the tangent linear model.
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
              CALL initialize_boundary (ng, tile, iTLM)
# endif
            END DO
!$OMP END PARALLEL
          END DO
!
!  Set FWDname to the new basic state just computed.
!
          DO ng=1,Ngrids
            WRITE (FWD(ng)%name,10) TRIM(FWD(ng)%base), outer
          END DO
!
!  Initialize tangent linear model from initial impulse which is now
!  stored in file ITL(ng)%name.
!
          DO ng=1,Ngrids
            wrtNLmod(ng)=.FALSE.
            wrtTLmod(ng)=.TRUE.
          END DO
!
!  If weak constraint, the impulses are time-interpolated at each
!  time-steps.
!
          DO ng=1,Ngrids
            IF (FrcRec(ng).gt.3) THEN
              FrequentImpulse(ng)=.TRUE.
            END IF
          END DO
!
!  Initialize tangent linear model from ITLname, record Rec3.
!  The sum of the initial condition increments is in record 3.
!
           DO ng=1,Ngrids
             ITL(ng)%Rindex=Rec3
!$OMP PARALLEL
             CALL tl_initial (ng)
!$OMP END PARALLEL
             IF (exit_flag.ne.NoError) RETURN
           END DO
!
!  Run tangent linear model.
!
           DO ng=1,Ngrids
             IF (Master) THEN
               WRITE (stdout,20) 'TL', ng, ntstart(ng), ntend(ng)
             END IF
           END DO

!$OMP PARALLEL
# ifdef SOLVE3D
           CALL tl_main3d (RunInterval)
# else
           CALL tl_main2d (RunInterval)
# endif
!$OMP END PARALLEL
           IF (exit_flag.ne.NoError) RETURN

           DO ng=1,Ngrids
             wrtNLmod(ng)=.FALSE.
             wrtTLmod(ng)=.FALSE.
           END DO
!
         END IF
!
#endif

      END DO OUTER_LOOP

#if defined POSTERIOR_ERROR_I || defined POSTERIOR_ERROR_F
!
!-----------------------------------------------------------------------
!  Compute full (diagonal) posterior analysis error covariance matrix.
!
!  NOTE: Currently, this code only works for a single outer-loop.
!-----------------------------------------------------------------------
!
!  Clear tangent and adjoint arrays because they are used as
!  work arrays below.
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
          CALL initialize_ocean (ng, tile, iTLM)
          CALL initialize_forces (ng, tile, iADM)
          CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iADM)
          CALL initialize_boundary (ng, tile, iTLM)
# endif
        END DO
!$OMF END PARALLEL
      END DO
!
!  Compute the diagonal of the posterior/analysis error covariance
!  matrix. The result is written to record 2 of the ITL netcdf file.
!
      VAR_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL posterior_var (ng, tile, iTLM, outer)
          END DO
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END DO VAR_OLOOP
!
!  Write out the diagonal of the posterior/analysis covariance matrix
!  which is in tl_var(Rec1) to 4DVar error NetCDF file.
!
      DO ng=1,Ngrids
        CALL wrt_error (ng, Rec1, Rec1)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Clear tangent and adjoint arrays because they are used as
!  work arrays below.
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
          CALL initialize_ocean (ng, tile, iTLM)
          CALL initialize_forces (ng, tile, iADM)
          CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iADM)
          CALL initialize_boundary (ng, tile, iTLM)
# endif
        END DO
!$OMF END PARALLEL
      END DO
#endif

#ifdef POSTERIOR_EOFS
!
!-----------------------------------------------------------------------
!  Compute the posterior analysis error covariance matrix EOFs using a
!  Lanczos algorithm.
!
!  NOTE: Currently, this code only works for a single outer-loop.
!-----------------------------------------------------------------------
!
      IF (Master) WRITE (stdout,50)
!
!  Estimate first the trace of the posterior analysis error
!  covariance matrix since the evolved and convolved Lanczos
!  vectors stored in the Hessian NetCDF file will be destroyed
!  later.
!
      Ltrace=.TRUE.

      TRACE_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0

        TRACE_ILOOP : DO my_inner=1,NpostI
          inner=my_inner
!
!  Initialize the tangent linear variables with a random vector
!  comprised of +1 and -1 elements randomly chosen.
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL random_ic (ng, tile, iTLM, inner, outer,             &
     &                        Lold(ng), Ltrace)
            END DO
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  Apply horizontal convolution.
!
          CALL convolve (driver, Lini, Lold, Lnew)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute Lanczos vector and eigenvectors of the posterior analysis
!  error covariance matrix.
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL posterior (ng, tile, iTLM, inner, outer, Ltrace)
            END DO
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO

        END DO TRACE_ILOOP

      END DO TRACE_OLOOP
!
!  Estimate posterior analysis error covariance matrix.
!
      Ltrace=.FALSE.

      POST_OLOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0
!
!  The Lanczos algorithm requires to save all the Lanczos vectors.
!  They are used to compute the posterior EOFs.
!
        DO ng=1,Ngrids
          ADM(ng)%Rindex=0
          Fcount=ADM(ng)%Fcount
          ADM(ng)%Nrec(Fcount)=0
        END DO

        POST_ILOOP : DO my_inner=0,NpostI
          inner=my_inner
!
!  Read first record of ITL file and apply convolutions.
!
!  NOTE: If inner=0, we would like to use a random starting vector.
!        For now we can use what ever is in record 1.
!
          IF (inner.ne.0) THEN
            DO ng=1,Ngrids
              Rec=1
              CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec, Lold(ng))
              IF (exit_flag.ne.NoError) RETURN
            END DO
          ELSE

            DO ng=1,Ngrids
!$OMP PARALLEL
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL random_ic (ng, tile, iTLM, inner, outer,           &
     &                          Lold(ng), Ltrace)
              END DO
!$OMP END PARALLEL
              IF (exit_flag.ne.NoError) RETURN
            END DO
          END IF
!
!  Apply horizontal convolution.
!
          CALL convolve (driver, Lini, Lold, Lnew)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute Lanczos vector and eigenvectors of the posterior analysis
!  error covariance matrix.
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL posterior (ng, tile, iTLM, inner, outer, Ltrace)
            END DO
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!   Write the Lanczos vectors of the posterior error covariance
!   to the adjoint NetCDF file.
!
          DO ng=1,Ngrids
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            Lfout(ng)=Lnew(ng)
# endif
# ifdef ADJUST_BOUNDARY
            Lbout(ng)=Lnew(ng)
# endif
            kstp(ng)=Lnew(ng)
# ifdef SOLVE3D
            nstp(ng)=Lnew(ng)
# endif
            LwrtState2d(ng)=.TRUE.
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
            LwrtState2d(ng)=.FALSE.
          END DO
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing adjustments for next inner
!  loop into ITL(ng)%name (record Rec1).
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, Lnew(ng), Rec1)
            IF (exit_flag.ne.NoError) RETURN
          END DO

        END DO POST_ILOOP

      END DO POST_OLOOP

#endif /* POSTERIOR_EOFS */
!
!  Done.  Set history file ID to closed state since we manipulated
!  its indices with the forward file ID which was closed above.
!
      DO ng=1,Ngrids
        HIS(ng)%ncid=-1
      END DO
!
 10   FORMAT (a,'_',i3.3,'.nc')
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (' (',i3.3,',',i3.3,'): ',a,' data penalty, Jdata = ',     &
     &        1p,e17.10,0p,t68,a)
 40   FORMAT (/,' Converting Convolved Adjoint Trajectory to',          &
     &          ' Impulses: Outer = ',i3.3,' Inner = ',i3.3,/)
#ifdef POSTERIOR_EOFS
 50   FORMAT (/,' <<<< Posterior Analysis Error Covariance Matrix',     &
     &          ' Estimation >>>>',/)
#endif

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear, tangent linear, and    !
!  adjoint models execution.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Local variable declarations.
!
      integer :: Fcount, ng, tile, thread
!
!-----------------------------------------------------------------------
!  Write out 4D-Var analysis fields that used as initial conditions for
!  the next data assimilation cycle.
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
      tile=MyRank
#else
      tile=-1
#endif
!
      IF (exit_flag.eq.NoError) THEN
        DO ng=1,Ngrids
          CALL wrt_dai (ng, tile)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Compute and report model-observation comparison statistics.
!-----------------------------------------------------------------------
!
      IF (exit_flag.eq.NoError) THEN
        DO ng=1,Ngrids
          CALL stats_modobs (ng)
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into record 3.
!
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
            Fcount=RST(ng)%Fcount
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
            CALL wrt_rst (ng)
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Close IO files.
!
      CALL close_out

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
