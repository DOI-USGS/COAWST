#include "cppdefs.h"
      MODULE ocean_control_mod
!
!svn $Id: is4dvar_ocean.h 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Strong Constraint 4-Dimensional Variational Data          !
!            Assimilation Driver, Incremental Approach (I4D-Var).      !
!            Primal formulation in model space.                        !
!                                                                      !
!  This driver is used for strong constraint 4D-Var where the only     !
!  errors considered are those for the observations.  The model is     !
!  assumed to be perfect.  This is the  incremental method and the     !
!  nonlinear, tangent linear and adjoint models are needed.            !
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
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE back_cost_mod, ONLY : back_cost
      USE cgradient_mod, ONLY : cgradient
      USE cost_grad_mod, ONLY : cost_grad
      USE ini_adjust_mod, ONLY : ini_adjust
      USE ini_fields_mod, ONLY : ini_fields
#ifdef ADJUST_BOUNDARY
      USE mod_boundary, ONLY : initialize_boundary
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces, ONLY : initialize_forces
#endif
      USE mod_ocean, ONLY : initialize_ocean
      USE normalization_mod, ONLY : normalization
      USE sum_grad_mod, ONLY : sum_grad
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      USE zeta_balance_mod, ONLY: balance_ref, biconj
#endif
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: converged
      logical :: Lweak = .FALSE.

      integer :: my_inner, my_outer
      integer :: AdjRec, Lbck, Lini, Lsav, Rec1, Rec2, Rec3, Rec4
      integer :: i, my_iic, ng, tile
      integer :: Lcon, LTLM1, LTLM2, LTLM3, LADJ1, LADJ2
      integer :: Fcount, NRMrec
      integer :: status

      real(r8) :: rate
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
!  Initialize relevant parameters.
!
      DO ng=1,Ngrids
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
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
      LTLM1=1               ! trial x-space TLM IC record in ITL
      LTLM2=2               ! previous v-space TLM IC record in ITL
      LTLM3=3               ! trial v-space TLM IC record in ITL
      LADJ1=1               ! initial cost gradient
      LADJ2=2               ! new cost gradient (not normalized)
      Lini=1                ! NLM initial conditions record in INI
      Lbck=2                ! background record in INI
      Rec1=1
      Rec2=2
      Rec3=3
      Rec4=4
      Nrun=1
      ERstr=1
      ERend=Nouter
!
!-----------------------------------------------------------------------
!  OUTER LOOP: time-step nonlinear model.
!-----------------------------------------------------------------------
!
      OUTER_LOOP : DO my_outer=1,Nouter
        outer=my_outer
        inner=0
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop.
!
        DO ng=1,Ngrids
          LdefHIS(ng)=.TRUE.
          LwrtHIS(ng)=.TRUE.
          WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%base), outer-1
        END DO

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
        IF (outer.eq.1) THEN
          DO ng=1,Ngrids
            BLK(ng)%name=HIS(ng)%name
          END DO
        END IF
#endif
!
!  Clear nonlinear mixing arrays.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_mixing (ng, tile, iNLM)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Initialize nonlinear model. If outer=1, the model is initialized
!  with the background or reference state. Otherwise, the model is
!  initialized with the estimated initial conditions from previous
!  iteration, X(0) = X(0) + deltaX(0).
!
        DO ng=1,Ngrids
          wrtNLmod(ng)=.TRUE.
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
!  If first pass, save nonlinear initial conditions (currently in time
!  index 1, background) into next record (Lbck) of INI(ng)%name NetCDF
!  file. The record "Lbck" becomes the background state record and the
!  record "Lini" becomes current nonlinear initial conditions.  Both
!  records are used in the algorithm below.
!
        IF (Nrun.eq.1) THEN
          DO ng=1,Ngrids
            INI(ng)%Rindex=1
            Fcount=INI(ng)%Fcount
            INI(ng)%Nrec(Fcount)=1
            CALL wrt_ini (ng, 1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
        END IF

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
        IF (balance(isFsur)) THEN
          DO ng=1,Ngrids
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
!  If first pass, compute or read in background-error covariance
!  normalization factors. If computing, write out factors to
!  NetCDF. This is an expensive computation that needs to be
!  computed only once for a particular application grid.
!
        IF (Nrun.eq.1) THEN
          DO ng=1,Ngrids
            IF (ANY(LwrtNRM(:,ng))) THEN
              CALL def_norm (ng, iNLM, 1)
              IF (exit_flag.ne.NoError) RETURN

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
        END IF
!
!  If first pass, define output 4DVAR NetCDF file containing all
!  processed data at observation locations.
!
        IF (Nrun.eq.1) THEN
          DO ng=1,Ngrids
            LdefMOD(ng)=.TRUE.
            CALL def_mod (ng)
            IF (exit_flag.ne.NoError) RETURN
          END DO
        END IF
!
!  Run nonlinear model. Save nonlinear tracjectory needed by the
!  adjoint and tangent linear models. Interpolate nonlinear model
!  to boservation locations (compute and save H x).
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
          wrtTLmod(ng)=.TRUE.
        END DO

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
!
!  Write out initial and background surface forcing into initial
!  INI(ng)%name NetCDF file for latter use.
!
        DO ng=1,Ngrids
          CALL wrt_frc (ng, Lfout(ng), Lini)
          IF (exit_flag.ne.NoError) RETURN
          IF (Nrun.eq.1) THEN
            CALL wrt_frc (ng, Lfout(ng), Lbck)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#endif
!
!  Write out nonlinear model misfit cost function into DAV(ng)%name
!  NetCDF file.
!
        SourceFile='is4dvar_ocean.h, ROMS_run'
        DO ng=1,Ngrids
          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name,                 &
     &                          'NLcost_function',                      &
     &                          FOURDVAR(ng)%NLobsCost(0:),             &
     &                          (/1,outer/), (/NstateVar(ng)+1,1/),     &
     &                          ncid = DAV(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!-----------------------------------------------------------------------
!  INNER LOOP: iterate using tangent linear model increments.
!-----------------------------------------------------------------------
!
!  The minimization algorithm requires to save all the gradient
!  solutions for each inner loop iteration.  They are used for
!  orthogonalization in the conjugate gradient algorithm.  Thus,
!  we need to reset adjoint file record indices.
!
        DO ng=1,Ngrids
          ADM(ng)%Rindex=0
          Fcount=ADM(ng)%Fcount
          ADM(ng)%Nrec(Fcount)=0
        END DO
!
!  An adjoint NetCDF is created for each outer loop.
!
        DO ng=1,Ngrids
          LdefADJ(ng)=.TRUE.
          WRITE (ADM(ng)%name,10) TRIM(ADM(ng)%base), outer
        END DO
!
!  Define output Hessian NetCDF file containing the eigenvectors
!  approximation to the Hessian matrix computed from the Lanczos
!  algorithm. Notice that the file name is a function of the
!  outer loop. That is, a file is created for each outer loop.
!
        DO ng=1,Ngrids
          WRITE (HSS(ng)%name,10) TRIM(HSS(ng)%base), outer
          LdefHSS(ng)=.TRUE.
          CALL def_hessian (ng)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Notice that inner loop iteration start from zero. This is needed to
!  compute the minimization initial increment deltaX(0), its associated
!  gradient G(0), and descent direction d(0) used in the conjugate
!  gradient algorithm.
!
        INNER_LOOP : DO my_inner=0,Ninner
          inner=my_inner
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time-step tangent linear model: compute cost function.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  If first pass inner=0, initialize tangent linear state (increments,
!  deltaX) from rest. Otherwise, use trial initial conditions estimated
!  by the conjugate gradient algorithm in previous inner loop. The TLM
!  initial conditions are read from ITL(ng)%name, record 1.
!
          DO ng=1,Ngrids
            ITL(ng)%Rindex=1
!$OMP PARALLEL
            CALL tl_initial (ng)
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  On first pass, initialize records 2, 3 and 4 of the ITL file to zero.
!
          IF (inner.eq.0.and.outer.eq.1) THEN
            DO ng=1,Ngrids
              CALL tl_wrt_ini (ng, LTLM1, Rec2)
              IF (exit_flag.ne.NoError) RETURN
              CALL tl_wrt_ini (ng, LTLM1, Rec3)
              IF (exit_flag.ne.NoError) RETURN
              CALL tl_wrt_ini (ng, LTLM1, Rec4)
              IF (exit_flag.ne.NoError) RETURN
            END DO
          END IF

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, activate writing and determine
!  output file name. The multiple file option is use to perturb initial
!  state and create ensembles.  The TLM final trajectory is written for
!  each inner loop on separated NetCDF files.
!
          DO ng=1,Ngrids
            LdefTLM(ng)=.TRUE.
            LwrtTLM(ng)=.TRUE.
            WRITE (TLM(ng)%name,10) TRIM(TLM(ng)%base), Nrun
          END DO
#endif
!
!  Activate switch to write out initial and final misfit between
!  model and observations.
!
          DO ng=1,Ngrids
            wrtMisfit(ng)=.FALSE.
            IF (((outer.eq.1).and.(inner.eq.0)).or.                     &
     &          ((outer.eq.Nouter).and.(inner.eq.Ninner))) THEN
              wrtMisfit(ng)=.TRUE.
            END IF
          END DO
!
!  Run tangent linear model. Compute misfit observation cost function,
!  Jo.
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

#ifdef EVOLVED_LCZ
!
!  Write evolved tangent Lanczos vector into hessian netcdf file for use
!  later.
!
!  NOTE: When using this option, it is important to set LhessianEV and
!  Lprecond to FALSE in s4dvar.in, otherwise the evolved Lanczos vectors
!  with be overwritten by the Hessian eigenvectors. The fix to this is to
!  define a new netcdf file that contains the evolved Lanczos vectors.
!
          IF (inner.ne.0) THEN
            DO ng=1,Ngrids
              CALL wrt_evolved (ng, kstp(ng), nrhs(ng))
              IF (exit_flag.ne.NoERRor) RETURN
            END DO
          END IF
#endif

#ifdef MULTIPLE_TLM
!
!  If multiple TLM history NetCDF files, close current NetCDF file.
!
          DO ng=1,Ngrids
            IF (TLM(ng)%ncid.ne.-1) THEN
              SourceFile='is4dvar_ocean.h, ROMS_run'

              CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
              IF (exit_flag.ne.NoError) RETURN
            END IF
          END DO
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Time step adjoint model backwards: compute cost function gradient.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
!  Time-step adjoint model backwards. The adjoint model is forced with
!  the adjoint of the observation misfit (Jo) term.
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
!  Clear adjoint arrays.  Is it needed?
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL initialize_ocean (ng, tile, iADM)
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
              CALL initialize_forces (ng, tile, iADM)
#endif
#ifdef ADJUST_BOUNDARY
              CALL initialize_boundary (ng, tile, iADM)
#endif
            END DO
!$OMP END PARALLEL
          END DO
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Descent algorithm.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Read TLM v-space initial conditions, record 3 in ITL(ng)%name, and
!  load it into time index LTLM1. This is needed to compute background
!  cost function. Also read in new (x-space) gradient vector, GRADx(Jo),
!  from adjoint history file ADM(ng)%name.  Read in the sum of all the
!  previous outer-loop increments which are always in record 4 of
!  the ITL file.
!
          DO ng=1,Ngrids
            IF (inner.eq.0) THEN
              CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec1, LTLM1)
              IF (exit_flag.ne.NoError) RETURN
            ELSE
              CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec3, LTLM1)
              IF (exit_flag.ne.NoError) RETURN
            END IF
            CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec4, LTLM2)
            IF (exit_flag.ne.NoError) RETURN
            CALL get_state (ng, iADM, 4, ADM(ng)%name, ADM(ng)%Rindex,  &
     &                      LADJ2)
            IF (exit_flag.ne.NoError) RETURN
#ifdef BALANCE_OPERATOR
            CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
            IF (exit_flag.ne.NoError) RETURN
            nrhs(ng)=Lini
#endif
          END DO
!
!  Convert observation cost function gradient, GRADx(Jo), from model
!  space (x-space) to minimization space (v-space):
!
!     GRADv(Jo) = B^(T/2) GRADx(Jo),  operator: S G L^(T/2) W^(-1/2)
!
!  First, multiply the adjoint solution, GRADx(Jo), by the background-
!  error standard deviations, S.  Second, convolve result with the
!  adjoint diffusion operator, G L^(T/2) W^(-1/2). Then, backgound
!  cost function contribution (BackCost) and cost function gradient
!  (v-space) by adding background and observation contributions:
!
!     GRADv(J) = GRADv(Jb) + GRADv(Jo) = deltaV + GRADv(Jo)
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
#ifdef BALANCE_OPERATOR
              CALL ad_balance (ng, tile, Lini, LADJ2)
#endif
              CALL ad_variability (ng, tile, LADJ2, Lweak)
              CALL ad_convolution (ng, tile, LADJ2, Lweak, 2)
              CALL cost_grad (ng, tile, LTLM1, LTLM2, LADJ2)
            END DO
!$OMP END PARALLEL
          END DO
!
!  Compute current total cost function.
!
          DO ng=1,Ngrids
            IF (Nrun.eq.1) THEN
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
                FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
              END DO
            ELSE
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
              END DO
            END IF
          END DO
!
!  Prepare for background cost function (Jb) calculation:
!
!  Read the convolved gradient from inner=0 (which is permanently
!  saved in record 1 of the adjoint file)  ALWAYS into record 1.
!
          IF (inner.gt.0) THEN
            DO ng=1,Ngrids
              CALL get_state (ng, iADM, 3, ADM(ng)%name, LADJ1,         &
     &                        LADJ1)
              IF (exit_flag.ne.NoError) RETURN
            END DO
          END IF
!
!  Compute background cost function (Jb) for inner=0:
!
!  If first pass of inner loop, read in the sum of previous v-space
!  gradients from record 4 of ITL file using the TLM model variables
!  as temporary storage. Also add background cost function to Cost0.
!
          IF (inner.eq.0) THEN
            DO ng=1,Ngrids
              CALL get_state (ng, iTLM, 2, ITL(ng)%name, Rec4, LTLM2)
              IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL back_cost (ng, tile, LTLM2)
              END DO
!$OMP END PARALLEL
!
              FOURDVAR(ng)%Cost0(outer)=FOURDVAR(ng)%Cost0(outer)+      &
     &                                  FOURDVAR(ng)%BackCost(0)
            END DO
          END IF
!
!  Compute current total cost function.
!
          DO ng=1,Ngrids
            IF (Nrun.eq.1) THEN
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostNorm(i)=FOURDVAR(ng)%CostNorm(i)+      &
     &                                   FOURDVAR(ng)%BackCost(i)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostNorm(i)
                FOURDVAR(ng)%CostFun(i)=FOURDVAR(ng)%CostNorm(i)
              END DO
            ELSE
              DO i=0,NstateVar(ng)
                FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
              END DO
            END IF
          END DO
!
!  Determine the descent direction in which the quadractic total cost
!  function decreases. Then, determine the TLM initial conditions,
!  deltaV(LTLM1), and its gradient, GRADv{J(Lnew)} at the new
!  direction.  Also, Compute TLM v-space trial initial conditions for
!  next inner loop, deltaV(LTLM2). The new gradient minimize the
!  quadratic cost function spanned by current and previous inner loop
!  iterations.  This is achieved by orthogonalizing (Gramm-Schmidt
!  algorithm) against all previous inner loop gradients.
!
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL cgradient (ng, tile, iTLM, inner, outer)
            END DO
!$OMP END PARALLEL
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  Report background (Jb) and observations (Jo) cost function values
!  normalized by their first minimization value. It also reports the
!  percentage change on total cost function value with respect to
!  previous iteration. Compute the optimality of the minimization to
!  check the statistical hypotheses between the background and
!  observations errors: the cost function value at the minimum, Jmin,
!  is idealy equal to half the number of observations assimilated
!  (Optimality=1=2*Jmin/Nobs), for a linear system.
!
          IF (Master) THEN
            DO ng=1,Ngrids
              IF (Nrun.gt.1) THEN
                rate=100.0_r8*ABS(FOURDVAR(ng)%CostFun(0)-              &
     &                            FOURDVAR(ng)%CostFunOld(0))/          &
     &                        FOURDVAR(ng)%CostFunOld(0)
              ELSE
                rate=0.0_r8
              END IF
              Optimality(ng)=2.0_r8*FOURDVAR(ng)%CostFun(0)/            &
     &                       (FOURDVAR(ng)%ObsCount(0)-                 &
     &                        FOURDVAR(ng)%ObsReject(0))
              WRITE (stdout,30) outer, inner,                           &
     &                          FOURDVAR(ng)%BackCost(0)/               &
     &                          FOURDVAR(ng)%CostNorm(0),               &
     &                          FOURDVAR(ng)%ObsCost(0)/                &
     &                          FOURDVAR(ng)%CostNorm(0),               &
     &                          rate
              IF (inner.eq.0) THEN
                DO i=0,NstateVar(ng)
                  IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
                    IF (i.eq.0) THEN
                      WRITE (stdout,40) outer, inner,                   &
     &                                  FOURDVAR(ng)%NLobsCost(i)/      &
     &                                  FOURDVAR(ng)%CostNorm(i)
                    ELSE
                      WRITE (stdout,50) outer, inner,                   &
     &                                  FOURDVAR(ng)%NLobsCost(i)/      &
     &                                  FOURDVAR(ng)%CostNorm(i),       &
     &                                  TRIM(Vname(1,idSvar(i)))
                    END IF
                  END IF
                  FOURDVAR(ng)%NLobsCost(i)=0.0
                END DO
              END IF
              WRITE (stdout,60) outer, inner, Optimality(ng)
            END DO
          END IF
!
!  Save total v-space cost function gradient, GRADv{J(Lnew)}, into
!  ADM(ng)%name history NetCDF file. Noticed that the lastest adjoint
!  solution record is over-written in the NetCDF file for future use.
!  The switch "LwrtState2d" is activated to write out state arrays
!  instead ad_*_sol arrays.
!
          DO ng=1,Ngrids
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            Lfout(ng)=LADJ2
#endif
#ifdef ADJUST_BOUNDARY
            Lbout(ng)=LADJ2
#endif
            kstp(ng)=LADJ2
#ifdef SOLVE3D
            nstp(ng)=LADJ2
#endif
            ADM(ng)%Rindex=ADM(ng)%Rindex-1
            LwrtState2d(ng)=.TRUE.
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
            LwrtState2d(ng)=.FALSE.
          END DO
!
!  Write out trial v-space TLM initial conditions, currently in time
!  index LTM2, into record 3 of ITL(ng)%name NetCDF file.
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, LTLM2, Rec3)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  Read current outer loop nonlinear model initial conditions and
!  background state vectors.
!
          DO ng=1,Ngrids
            CALL get_state (ng, iNLM, 2, INI(ng)%name, Lini, Lini)
            IF (exit_flag.ne.NoError) RETURN
            CALL get_state (ng, iNLM, 9, INI(ng)%name, Lbck, Lbck)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!  Convert increment vector, deltaV, from minimization space (v-space)
!  to model space (x-space):
!
!     deltaX = B^(1/2) deltaV
!  or
!     deltaX = W^(1/2) L^(1/2) G S
!
!  First, convolve estimated increment vector (v-space) by with the
!  tangent linear diffusion operator, W^(1/2) L^(1/2) G.  Second,
!  multiply result by the background-error standard deviation, S.
!
          Lcon=LTLM2
          DO ng=1,Ngrids
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_convolution (ng, tile, Lcon, Lweak, 2)
              CALL tl_variability (ng, tile, Lcon, Lweak)
#ifdef BALANCE_OPERATOR
              CALL tl_balance (ng, tile, Lini, Lcon)
#endif
            END DO
!$OMP END PARALLEL
          END DO
!
!  Write out trial x-space (convolved) TLM initial conditions, currently
!  in time index Lcon, into record 1 of ITL(ng)%name NetCDF file.
!
          DO ng=1,Ngrids
            CALL tl_wrt_ini (ng, Lcon, Rec1)
            IF (exit_flag.ne.NoError) RETURN
          END DO
!
!-----------------------------------------------------------------------
!  Update counters.
!-----------------------------------------------------------------------
!
          DO ng=1,Ngrids
            Lsav=Lnew(ng)
            Lnew(ng)=Lold(ng)
            Lold(ng)=Lsav
            Nrun=Nrun+1
          END DO

        END DO INNER_LOOP
!
!  Close adjoint NetCDF file.
!
        DO ng=1,Ngrids
          IF (ADM(ng)%ncid.ne.-1) THEN
            SourceFile='is4dvar_ocean.h, ROMS_run'

            CALL netcdf_close (ng, iADM, ADM(ng)%ncid)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Close Hessian NetCDF file.
!
        DO ng=1,Ngrids
          IF (HSS(ng)%ncid.ne.-1) THEN
            SourceFile='is4dvar_ocean.h, ROMS_run'

            CALL netcdf_close (ng, iADM, HSS(ng)%ncid)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Clear nonlinear state variables.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iNLM)
          END DO
!$OMP END PARALLEL
        END DO
!
!-----------------------------------------------------------------------
!  Compute new nonlinear initial conditions by adding minimization
!  increment to previous outer loop initial conditions:
!
!         Xi(outer+1) = Xi(outer) + deltaX(Lcon)
!
!-----------------------------------------------------------------------
!
!  Notice that "ini_fields" is called here for output purposes only.
!  It computes the vertically integrated momentum in 3D applications.
!  In order to use the correct fields, the model time indices are set
!  to Lini.
!  The appropriate tl correction for the NL model resides in record 1
!  of the ITL file.
!
        DO ng=1,Ngrids
          kstp(ng)=Lini
#ifdef SOLVE3D
          nstp(ng)=Lini
#endif
          CALL get_state (ng, iNLM, 1, INI(ng)%name, Lini, Lini)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_state (ng, iTLM, 1, ITL(ng)%name, LTLM1, LTLM1)
          IF (exit_flag.ne.NoError) RETURN

!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL ini_adjust (ng, tile, LTLM1, Lini)
            CALL ini_fields (ng, tile, iNLM)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Write out new nonlinear model initial conditions into record Lini
!  of INI(ng)%name.
!
        DO ng=1,Ngrids
          INI(ng)%Rindex=0
          Fcount=INI(ng)%Fcount
          INI(ng)%Nrec(Fcount)=1
          CALL wrt_ini (ng, Lini)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
! Gather the v-space increments from the final inner-loop and
! save in record 4 of the ITL file. The current v-space increment
! is in record 3 and the sum so far is in record 4.
!
        DO ng=1,Ngrids
          CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec3, LTLM1)
          IF (exit_flag.ne.NoError) RETURN
          CALL get_state (ng, iTLM, 8, ITL(ng)%name, Rec4, LTLM2)
          IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL sum_grad (ng, tile, LTLM1, LTLM2)
          END DO
!$OMP END PARALLEL
        END DO
!
! Write the current sum into record 4 of the ITL file.
!
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, LTLM2, Rec4)
          IF (exit_flag.ne.NoError) RETURN
        END DO

#if defined ADJUST_STFLUX   || defined ADJUST_WSTRESS || \
    defined ADJUST_BOUNDARY
!
!  Set index containing the surface forcing increments used the run
!  the nonlinear model in the outer loop and read the forcing
!  increments. For bulk fluxes, we read Rec1 because the stress
!  fluxes change by virtue of the changing initial conditions.
!  When not using bulk fluxes, we read Rec4 because the background
!  stress and flux is prescribed by input files which are not
!  overwritten so we need to correct the background using the
!  sum of the increments from all previous outer-loops.
!  If using Rec4 we need to convert from v-space to x-space
!  by applying the convolution.
!  Note that Lfinp=Lbinp so the the forcing and boundary
!  adjustments are both processsed correctly.
# ifdef BALANCE_OPERATOR
!  Currently, We do not need the call to tl_balance below, but we
!  might later if we impose a balance constraint on the wind stress
!  corrections.
# endif
!
!  AMM: CHECK WHAT HAPPENS WITH SECONDARY PRECONDITIONING.
!
        DO ng=1,Ngrids
          Lfinp(ng)=LTLM1
# if defined BULK_FLUXES && !defined NL_BULK_FLUXES
          CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec1, Lfinp(ng))
# endif
# if defined NL_BULK_FLUXES || !defined BULK_FLUXES
          CALL get_state (ng, iTLM, 1, ITL(ng)%name, Rec4, Lfinp(ng))
          Lcon=Lfinp(ng)
!
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tl_convolution (ng, tile, Lcon, Lweak, 2)
            CALL tl_variability (ng, tile, Lcon, Lweak)
# ifdef BALANCE_OPERATOR
!!          CALL tl_balance (ng, tile, Lini, Lcon)
# endif
          END DO
!$OMP END PARALLEL
# endif
          IF (exit_flag.ne.NoError) RETURN
        END DO
#endif
!
!-----------------------------------------------------------------------
!  Clear tangent linear state variables.
!-----------------------------------------------------------------------
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, iTLM)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Close current forward NetCDF file.
!
        SourceFile='is4dvar_ocean.h, ROMS_run'
        DO ng=1,Ngrids
          CALL netcdf_close (ng, iNLM, FWD(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
        END DO

      END DO OUTER_LOOP
!
!-----------------------------------------------------------------------
!  Done with data assimilation. Initialize the nonlinear model with
!  estimated initial conditions. Save nonlinear solution at observation
!  points for posterior analysis.
!-----------------------------------------------------------------------
!
!  Set nonlinear output history file name. Create a basic state file
!  for each outher loop.
!
      DO ng=1,Ngrids
        LdefHIS(ng)=.TRUE.
        LwrtHIS(ng)=.TRUE.
        HIS(ng)%Rindex=0
        Fcount=HIS(ng)%Fcount
        HIS(ng)%Nrec(Fcount)=0
        WRITE (HIS(ng)%name,10) TRIM(FWD(ng)%base), Nouter
      END DO
!
!  Clear nonlinear mixing arrays.
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_mixing (ng, tile, iNLM)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Initialize nonlinear model with estimated initial conditions.
!
      DO ng=1,Ngrids
        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.
        wrtMisfit(ng)=.FALSE.
        RST(ng)%Rindex=0
        Fcount=RST(ng)%Fcount
        RST(ng)%Nrec(Fcount)=0
      END DO

!$OMP PARALLEL
      CALL initial
!$OMP END PARALLEL
      IF (exit_flag.ne.NoError) RETURN
!
! Clear NLobsCost.
!
      DO ng=1,Ngrids
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%NLobsCost(i)=0.0_r8
        END DO
      END DO
!
!  Run nonlinear model. Interpolate nonlinear model to observation
!  locations.
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
!
!  Write out nonlinear model final misfit cost function into
!  DAV(ng)%name NetCDF file. Notice that it is written in the
!  Nouter+1 record.
!
      SourceFile='is4dvar_ocean.h, ROMS_run'
      DO ng=1,Ngrids
        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'NLcost_function',&
     &                        FOURDVAR(ng)%NLobsCost(0:),               &
     &                        (/1,Nouter+1/), (/NstateVar(ng)+1,1/),    &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Report the final value of the nonlinear model misfit cost function.
!
      IF (Master) THEN
        DO ng=1,Ngrids
          DO i=0,NstateVar(ng)
            IF (FOURDVAR(ng)%NLobsCost(i).ne.0.0_r8) THEN
              IF (i.eq.0) THEN
                WRITE (stdout,40) outer, inner,                         &
     &                            FOURDVAR(ng)%NLobsCost(i)/            &
     &                            FOURDVAR(ng)%CostNorm(i)
              ELSE
                WRITE (stdout,50) outer, inner,                         &
     &                            FOURDVAR(ng)%NLobsCost(i)/            &
     &                            FOURDVAR(ng)%CostNorm(i),             &
     &                            TRIM(Vname(1,idSvar(i)))
              END IF
            END IF
          END DO
        END DO
      END IF
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
 30   FORMAT (/,' (',i3.3,',',i3.3,'): TLM Cost Jb, J  = ',             &
     &        1p,e17.10,0p,1x,1p,e17.10,0p,t68,1p,e11.4,' %')
 40   FORMAT (/,'>(',i3.3,',',i3.3,'): NLM Cost     J  = ',             &
     &        18x,1p,e17.10,0p)
 50   FORMAT (' (',i3.3,',',i3.3,'): NLM Cost     J  = ',               &
     &        18x,1p,e17.10,0p,t69,a)
 60   FORMAT (/,1x,'(',i3.3,',',i3.3,'): Optimality (2*J/Nobs) = ',     &
     &        1p,e17.10,/)

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
