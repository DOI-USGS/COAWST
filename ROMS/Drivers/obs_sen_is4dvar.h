      MODULE ocean_control_mod
!
!svn $Id: obs_sen_is4dvar.h 431 2009-12-26 20:36:20Z arango $
!=================================================== Andrew M. Moore ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS 4DVAR Observation Sensitivity Analysis Driver:            !
!                                                                      !
!  This driver evaluates the impact of each observation in the         !
!  4DVAR analysis increment by measuring their sensitivity over        !
!  a specified circulation functional index, J, similar to the         !
!  adjoint sensitivity driver. This is equivalent to taking the        !
!  adjoint of the 4DVAR algorithm.                                     !
!                                                                      !
!  Algorithm Outline:                                                  !
!                                                                      !
!                                                                      !
!    |----------------------------|-----------------------------|      !
!   t=t0                         t=t1                          t=t2    !
!         Assimilation window             Forecast period              !
!                                                                      !
!    |----------------------------------------------------------> NLM  !
!                                                                      !
!    <----------------------------------------------------------| ADM  !
!                                                                      !
!    |----------------------------> TLM                                !
!                                                                      !
!  (1) We begin by running an IS4DVAR Lanczos calculation using k      !
!      inner-loops and 1 outer-loop for the period t=t0 to t1. We      !
!      will denote by xb(0) the background initial condition, and      !
!      the observations vector by y.  The resulting Lanczos vectors    !
!      that we save in the adjoint NetCDF file will be denoted by      !
!      q_i, where i=1,2,...,k.                                         !
!                                                                      !
!  (2) Next we run the NLM for the combined assimilation+forecast      !
!      period t=t0 to t2, where t2>t1. This represents the final       !
!      sweep of the NLM for the period t=t0 to t1 after exiting the    !
!      inner-loop in the IS4DVAR plus the forecast period t=t1 to t2.  !
!      The initial condition for the NLM at t=t0 is xb(0) and not      !
!      the new estimated initial conditions. We save the basic state   !
!      trajectory, xb(t), of this NLM run for use in the adjoint       !
!      sensitivity calculation next, and for use in the TLM run        !
!      later in step (7).                                              !
!                                                                      !
!      Depending on time for which the sensitivity functional J(t)     !
!      is defined, this will dictate t2. For example, if J(t) is a     !
!      functional defined during the forecast interval t1<t<t2, then   !
!      t2>t1 for this run of the NLM. However, if J(t) is defined      !
!      during the assimilation interval t0<t<t1, then t2=t1. That is,  !
!      the definition of t2 should be flexible depending on the        !
!      choice of J.                                                    !
!                                                                      !
!  (3) The next step involves an adjoint sensitivity calculation       !
!      for the combined assimilation+forecast period t=t0 to t2. The   !
!      basic state trajectory for this calculation will be that from   !
!      the NLM run in step (2).                                        !
!                                                                      !
!  (4) After running the regular adjoint sensitivity calculation in    !
!      (3), we will have a full 3D-adjoint state vector at time t=t0.  !
!      Let's call this vector x(0). The next thing we want to do is    !
!      to compute the dot-product of x(0) with each of the Lanczos     !
!      vectors from the previous IS4DVAR run. So if we ran IS4DVAR     !
!      with k inner-loops we will have k Lanczos vectors which we      !
!      denote as q_i where i=1,2,...,k. So we will compute a_i=x'q_i   !
!      where x' is the transpose of the vector x(0), and a_i for       !
!      i=1,2,...,k are scalars, so there will be k of them.            !
!                                                                      !
!  (5) The next step is to invert the tridiagonal matrix associated    !
!      with the Lanczos vectors. Let's denote this matrix as T. So     !
!      what we want to solve T*b=a, where a is the k by 1 vector of    !
!      scalars a_i from step (4), and b is the k by 1 vector that we   !
!      want to find. So we solve for b by using a tridiagonal solver.  !
!                                                                      !
!  (6) The next step is to compute a weighted sum of the Lanczos       !
!      vectors. Let's call this z, where:                              !
!                                                                      !
!      z = SUM_i (b_i * q_i)    and i=1,2,...,k                        !
!                                                                      !
!      The b_i are obtained from solving the tridiagonal equation in   !
!      (5), and the q_i are the Lanczos vectors. The vector z is a     !
!      full-state vector and be used as an initial condition for the   !
!      TLM in step (8).                                                !
!                                                                      !
!  (7) Finally, we run the TLM from t=t0 to t=t1 using z from (6) as   !
!      the TLM initial condition. During this run of the TLM, we       !
!      need to read and process the observations that we used in       !
!      the IS4DVAR of step (1) and write the TLM solution at the       !
!      observation points and times to the MODname NetCDF file.        !
!      The values that we write into this MODname are actually the     !
!      TLM values multiplied by error variance assigned to each        !
!      observation during the IS4DVAR in step (1).                     !
!                                                                      !
!  These routines control the initialization,  time-stepping,  and     !
!  finalization of  ROMS/TOMS  model following ESMF conventions:       !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize

      CONTAINS

      SUBROUTINE ROMS_initialize (first, MyCOMM)
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
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
#ifdef AIR_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_atmos_coupling
#endif
#ifdef WAVES_OCEAN
      USE ocean_coupler_mod, ONLY : initialize_waves_coupling
#endif
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

      integer :: NRMrec, STDrec, Tindex, ng, thread

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (MPI) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(MyCOMM)) THEN
        OCN_COMM_WORLD=MyCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
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
!  Initialize parallel parameters.
!
        CALL initialize_parallel
!
!  Initialize wall clocks.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (' Process Information:',/)
        END IF
        DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
          DO thread=0,numthreads-1
            CALL wclock_on (ng, iADM, 0)
          END DO
!$OMP END PARALLEL DO
        END DO

#if defined AIR_OCEAN || defined WAVES_OCEAN
!
!  Initialize coupling streams between model(s).
!
        DO ng=1,Ngrids
# ifdef AIR_OCEAN
          CALL initialize_atmos_coupling (ng, MyRank)
# endif
# ifdef WAVES_OCEAN
          CALL initialize_waves_coupling (ng, MyRank)
# endif
        END DO
#endif
!
!  Read in model tunable parameters from standard input. Initialize
!  "mod_param", "mod_ncparam" and "mod_scalar" modules.
!
        CALL inp_par (iADM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
      END IF
!
!-----------------------------------------------------------------------
!  Read in Lanczos algorithm coefficients (cg_beta, cg_delta) from
!  file LCZname NetCDF (IS4DVAR adjoint file), as computed in the
!  IS4DVAR Lanczos data assimilation algorithm for the first outer
!  loop.  They are needed here, in routine "ini_lanczos", to compute
!  the tangent linear model initial conditions as the weighted sum
!  of the Lanczos vectors. The weighting coefficient are computed
!  by solving a tri-diagonal system that uses cg_beta and cg_gamma.
!-----------------------------------------------------------------------
!
      SourceFile='obs_sen_ocean.h, ROMS_initialize'

      DO ng=1,Ngrids
        CALL netcdf_get_fvar (ng, iADM, LCZname(ng), 'cg_beta',         &
     &                        cg_beta)
        IF (exit_flag.ne. NoError) RETURN
        CALL netcdf_get_fvar (ng, iADM, LCZname(ng), 'cg_delta',        &
     &                        cg_delta)
        IF (exit_flag.ne. NoError) RETURN
      END DO
!
!-----------------------------------------------------------------------
!  Read in error covariance standard deviation factors used in the
!  spatial convolutions.
!-----------------------------------------------------------------------
!
!  Standard deviation factors for initial conditions error covariance.
!  They are loaded in Tindex=1 of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 6, 6, STDname(1,ng), STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Standard deviation factors for model error covariance. They are
!  loaded in Tindex=2 of the e_var(...,Tindex) state variables.
!
      STDrec=1
      Tindex=2
      DO ng=1,Ngrids
        IF (NSA.eq.2) THEN
          CALL get_state (ng, 6, 6, STDname(2,ng), STDrec, Tindex)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO

#ifdef ADJUST_BOUNDARY
!
!  Standard deviation factors for boundary conditions error covariance.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 8, 8, STDname(3,ng), STDrec, Tindex)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Standard deviation factors for boundary conditions error covariance.
!
      STDrec=1
      Tindex=1
      DO ng=1,Ngrids
        CALL get_state (ng, 9, 9, STDname(4,ng), STDrec, 1)
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Read in initial conditions, model, boundary conditions, and surface
!  forcing error covariance covariance normalization factors.
!-----------------------------------------------------------------------
!
      NRMrec=1
      DO ng=1,Ngrids
        CALL get_state (ng, 5, 5, NRMname(1,ng), NRMrec, 1)
        IF (exit_flag.ne.NoError) RETURN

        IF (NSA.eq.2) THEN
          CALL get_state (ng, 5, 5, NRMname(2,ng), NRMrec, 2)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#ifdef ADJUST_BOUNDARY
        CALL get_state (ng, 10, 10, NRMname(3,ng), NRMrec, 1)
        IF (exit_flag.ne.NoError) RETURN
#endif
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        CALL get_state (ng, 11, 11, NRMname(4,ng), NRMrec, 1)
        IF (exit_flag.ne.NoError) RETURN
#endif
      END DO

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine computes the adjoint sensitivity analysis, dJ/dS,      !
!  to the specified functional J.  The sensitivity masking arrays      !
!  Rscope, Uscope, and Vscope are used to evaluate the functional      !
!  in the desired spatial area.                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
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
#ifdef ADJUST_BOUNDARY
      USE mod_boundary, ONLY : initialize_boundary
#endif
#if defined OBS_IMPACT && defined OBS_IMPACT_SPLIT
      USE mod_forces, ONLY : initialize_forces
      USE mod_ocean, ONLY : initialize ocean
#endif
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
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: Lweak = .FALSE.

      integer :: i, my_iic, ng,  subs, tile, thread
      integer :: Lbck, Lini, Litl, Rec

      real (r8) :: str_day, end_day
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!  Initialize relevant parameters.
!
        Lini=1        ! 4DVAR initial conditions record in INIname
        Lbck=2        ! First guess initial conditions record in INIname
        Litl=1        ! TLM initial conditions record
        Lnew(ng)=1
!
!  Initialize nonlinear model with the estimated initial conditions
!  from the IS4DVAR Lanczos algorithm.
!
        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
        tRSTindx(ng)=0
        NrecRST(ng)=0
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN

#if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Compute the reference zeta and biconjugate gradient arrays
!  required for the balance of free surface.
!
          IF (balance(isFsur)) THEN
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL balance_ref (ng, TILE, Lini)
                CALL biconj (ng, TILE, iNLM, Lini)
              END DO
            END DO
!$OMP END PARALLEL DO
            wrtZetaRef(ng)=.TRUE.
          END IF
#endif
!
!  Run nonlinear model for the combined assimilation plus forecast
!  period, t=t0 to t2. Save nonlinear (basic state) tracjectory, xb(t),
!  needed by the adjoint model.
!
        IF (Master) THEN
          WRITE (stdout,10) 'NL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        NL_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL main3d (ng)
#else
          CALL main2d (ng)
#endif
          IF (exit_flag.ne.NoError) RETURN

        END DO NL_LOOP1
!
!  Initialize adjoint model and define sensitivity functional.
!
        Lstiffness=.FALSE.
        CALL ad_initial (ng)
        IF (exit_flag.ne.NoError) RETURN

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
        BLKname(ng)=HISname(ng)
#endif
!
!  Activate adjoint output.
!
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LcycleADJ(ng)=.FALSE.
!
!  Time-step adjoint model for the combined plus forecast period,
!  t=t2 to t0. Compute the gradient or index, dJ/dS, of the
!  sensitivity functional.
!
        str_day=time(ng)*sec2day
        end_day=str_day-ntimes(ng)*dt(ng)*sec2day
        IF ((DstrS(ng).eq.0.0_r8).and.(DendS(ng).eq.0.0_r8)) THEN
          DstrS(ng)=end_day
          DendS(ng)=str_day
        END IF
        IF (Master) THEN
          WRITE (stdout,20) 'AD', ntstart(ng), ntend(ng),               &
     &                      DendS(ng), DstrS(ng)
        END IF
#ifndef OBS_IMPACT
        IF ((DstrS(ng).gt.str_day).or.(DstrS(ng).lt.end_day)) THEN
          IF (Master)  WRITE (stdout,30) 'DstrS = ', DstrS(ng),         &
     &                                   end_day, str_day
          exit_flag=7
          RETURN
        END IF
        IF ((DendS(ng).gt.str_day).or.(DendS(ng).lt.end_day)) THEN
          IF (Master)  WRITE (stdout,30) 'DendS = ', DendS(ng),         &
     &                                   end_day, str_day
          exit_flag=7
          RETURN
        END IF
#endif

!
!   If DstrS=DendS=dstart, skip the adjoint model and read adjoint solution
!   from ADS netcdf file, record 1.
!
        IF ((DstrS(ng).eq.DendS(ng)).and.(DstrS(ng).eq.dstart)) THEN
          Rec=1
          CALL get_state (ng, iADM, 4, ADSname(ng), Rec, Lnew(ng))
        ELSE
          time(ng)=time(ng)+dt(ng)

          AD_LOOP : DO my_iic=ntstart(ng),ntend(ng),-1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL ad_main3d (ng)
#else
            CALL ad_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO AD_LOOP
!
!  Load full adjoint sensitivity vector, x(0), for t=t0 into adjoint
!  state arrays at index Lnew.
!
          CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),       &
     &                    Lnew(ng))
        END IF
#ifdef BALANCE_OPERATOR
        CALL get_state (ng, iNLM, 2, INIname(ng), Lini, Lini)
        IF (exit_flag.ne.NoError) RETURN
        nrhs(ng)=Lini
#endif
!
!  Convert adjoint solution to v-space since the Lanczos vectors
!  are in v-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
#ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lini, Lnew(ng))
#endif
            CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
            CALL ad_convolution (ng, TILE, Lnew(ng), Lweak, 2)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Check Lanczos vector input file and determine t=t1. That is, the
!  time to run the tangent linear model.  This time must be the same
!  as the IS4DVAR Lanczos algorithm.
!
        SourceFile='obs_sen_ocean.h, ROMS_run'

        CALL netcdf_get_ivar (ng, iADM, LCZname(ng), 'ntimes',          &
     &                        ntimes(ng))
        IF (exit_flag.ne. NoError) RETURN

#ifndef OBS_IMPACT
!
!  Initialize nonlinear model with the same initial conditions, xb(0),
!  Lbck record in INIname. This is the first guess NLM initial
!  conditions used to start the IS4DVAR Lanczos algorithm.
!
        LdefINI(ng)=.FALSE.
        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
        tRSTindx(ng)=0
        tINIindx(ng)=Lbck
        NrecRST(ng)=0
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Run nonlinear model for the combined assimilation plus forecast
!  period, t=t0 to t2. Save nonlinear (basic state) tracjectory, xb(t),
!  needed by the tangent linear model.
!
        IF (Master) THEN
          WRITE (stdout,10) 'NL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        NL_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
# ifdef SOLVE3D
          CALL main3d (ng)
# else
          CALL main2d (ng)
# endif
          IF (exit_flag.ne.NoError) RETURN

        END DO NL_LOOP2
#endif
!
!  Initialize tangent linear model with the weighted sum of the
!  Lanczos vectors, steps (4) to (6) from the algorithm summary
!  above.
!
        CALL tl_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        wrtTLmod(ng)=.TRUE.
!
!  Convert TL initial condition from v-space to x-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_convolution (ng, TILE, Litl, Lweak, 2)
            CALL tl_variability (ng, TILE, Litl, Lweak)
#ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, TILE, Lini, Litl)
#endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Define output 4DVAR NetCDF file containing the sensitivity at the
!  observation locations.
!
        LdefMOD(ng)=.TRUE.
        CALL def_mod (ng)
        IF (exit_flag.ne.NoError) RETURN
        wrtIMPACT_TOT(ng)=.TRUE.
#ifdef OBS_IMPACT_SPLIT
        wrtIMPACT_IC(ng)=.FALSE.
#endif
!
!  Run tangent linear model for the assimilation period, t=t0 to t1.
!  Read and process the 4DVAR observations.
!
        IF (Master) THEN
          WRITE (stdout,10) 'TL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        TL_LOOP1 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL tl_main3d (ng)
#else
          CALL tl_main2d (ng)
#endif
          IF (exit_flag.ne.NoError) RETURN

        END DO TL_LOOP1

#if defined OBS_IMPACT && defined OBS_IMPACT_SPLIT
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with initial condition increments
!  only to compute the observation impact associated with the initial
!  conditions.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Load full adjoint sensitivity vector, x(0), for t=t0 into adjoint
!  state arrays at index Lnew.
!
        CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),         &
     &                    Lnew(ng))
# ifdef BALANCE_OPERATOR
        CALL get_state (ng, iNLM, 2, INIname(ng), Lini, Lini)
        IF (exit_flag.ne.NoError) RETURN
        nrhs(ng)=Lini
# endif
!
!  Clear the adjoint forcing and boundary condition increment arrays.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL initialize_forces (ng, TILE, iADM)
# ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, TILE, iADM)
# endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Convert adjoint solution to v-space since the Lanczos vectors
!  are in v-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
# ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lini, Lnew(ng))
# endif
            CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
            CALL ad_convolution (ng, TILE, Lnew(ng), Lweak, 2)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Check Lanczos vector input file and determine t=t1. That is, the
!  time to run the tangent linear model.  This time must be the same
!  as the IS4DVAR Lanczos algorithm.
!
        SourceFile='obs_sen_ocean.h, ROMS_run'

        CALL netcdf_get_ivar (ng, iADM, LCZname(ng), 'ntimes',          &
     &                        ntimes(ng))
        IF (exit_flag.ne. NoError) RETURN
!
!  Initialize tangent linear model with the weighted sum of the
!  Lanczos vectors, steps (4) to (6) from the algorithm summary
!  above.
!
        CALL tl_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        wrtTLmod(ng)=.TRUE.
        wrtIMPACT_TOT(ng)=.FALSE.
        wrtIMPACT_IC(ng)=.TRUE.
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        wrtIMPACT_FC(ng)=.FALSE.
# endif
# if defined ADJUST_BOUNDARY
        wrtIMPACT_BC(ng)=.FALSE.
# endif
!
!  Convert TL initial condition from v-space to x-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_convolution (ng, TILE, Litl, Lweak, 2)
            CALL tl_variability (ng, TILE, Litl, Lweak)
# ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, TILE, Lini, Litl)
# endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Run tangent linear model for the assimilation period, t=t0 to t1.
!  Read and process the 4DVAR observations.
!
        IF (Master) THEN
          WRITE (stdout,10) 'TL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        TL_LOOP2 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
# ifdef SOLVE3D
          CALL tl_main3d (ng)
# else
          CALL tl_main2d (ng)
# endif
          IF (exit_flag.ne.NoError) RETURN

        END DO TL_LOOP2

# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with surface forcing increments only
!  to compute the observation impact associated with the surface forcing.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Load full adjoint sensitivity vector, x(0), for t=t0 into adjoint
!  state arrays at index Lnew.
!
        CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),         &
     &                    Lnew(ng))
#  ifdef BALANCE_OPERATOR
        CALL get_state (ng, iNLM, 2, INIname(ng), Lini, Lini)
        IF (exit_flag.ne.NoError) RETURN
        nrhs(ng)=Lini
#  endif
!
!  Clear the adjoint initial condition and boundary condition increment
!  arrays.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL initialize_ocean (ng, TILE, iADM)
#  ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, TILE, iADM)
#  endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Convert adjoint solution to v-space since the Lanczos vectors
!  are in v-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
#  ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lini, Lnew(ng))
#  endif
            CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
            CALL ad_convolution (ng, TILE, Lnew(ng), Lweak, 2)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Check Lanczos vector input file and determine t=t1. That is, the
!  time to run the tangent linear model.  This time must be the same
!  as the IS4DVAR Lanczos algorithm.
!
        SourceFile='obs_sen_ocean.h, ROMS_run'

        CALL netcdf_get_ivar (ng, iADM, LCZname(ng), 'ntimes',          &
     &                        ntimes(ng))
        IF (exit_flag.ne. NoError) RETURN
!
!  Initialize tangent linear model with the weighted sum of the
!  Lanczos vectors, steps (4) to (6) from the algorithm summary
!  above.
!
        CALL tl_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        wrtTLmod(ng)=.TRUE.
        wrtIMPACT_TOT(ng)=.FALSE.
        wrtIMPACT_IC(ng)=.FALSE.
        wrtIMPACT_FC(ng)=.TRUE.
#  if defined ADJUST_BOUNDARY
        wrtIMPACT_BC(ng)=.FALSE.
#  endif
!
!  Convert TL initial condition from v-space to x-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_convolution (ng, TILE, Litl, Lweak, 2)
            CALL tl_variability (ng, TILE, Litl, Lweak)
#  ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, TILE, Lini, Litl)
#  endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Run tangent linear model for the assimilation period, t=t0 to t1.
!  Read and process the 4DVAR observations.
!
        IF (Master) THEN
          WRITE (stdout,10) 'TL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        TL_LOOP3 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#  ifdef SOLVE3D
          CALL tl_main3d (ng)
#  else
          CALL tl_main2d (ng)
#  endif
          IF (exit_flag.ne.NoError) RETURN

        END DO TL_LOOP3
# endif

# if defined ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Integrate tangent linear model with boundary increments only
!  to compute the observation impact associated with the boundaries.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  Load full adjoint sensitivity vector, x(0), for t=t0 into adjoint
!  state arrays at index Lnew.
!
        CALL get_state (ng, iADM, 4, ADJname(ng), tADJindx(ng),         &
     &                    Lnew(ng))
#  ifdef BALANCE_OPERATOR
        CALL get_state (ng, iNLM, 2, INIname(ng), Lini, Lini)
        IF (exit_flag.ne.NoError) RETURN
        nrhs(ng)=Lini
#  endif
!
!  Clear the adjoint increment initial condition and forcing arrays.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL initialize_ocean (ng, TILE, iADM)
            CALL initialize_forces (ng, TILE, iADM)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Convert adjoint solution to v-space since the Lanczos vectors
!  are in v-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
#  ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, TILE, Lini, Lnew(ng))
#  endif
            CALL ad_variability (ng, TILE, Lnew(ng), Lweak)
            CALL ad_convolution (ng, TILE, Lnew(ng), Lweak, 2)
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Check Lanczos vector input file and determine t=t1. That is, the
!  time to run the tangent linear model.  This time must be the same
!  as the IS4DVAR Lanczos algorithm.
!
        SourceFile='obs_sen_ocean.h, ROMS_run'

        CALL netcdf_get_ivar (ng, iADM, LCZname(ng), 'ntimes',          &
     &                        ntimes(ng))
        IF (exit_flag.ne. NoError) RETURN
!
!  Initialize tangent linear model with the weighted sum of the
!  Lanczos vectors, steps (4) to (6) from the algorithm summary
!  above.
!
        CALL tl_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        wrtTLmod(ng)=.TRUE.
        wrtIMPACT_TOT(ng)=.FALSE.
        wrtIMPACT_IC(ng)=.FALSE.
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
        wrtIMPACT_FC(ng)=.FALSE.
#  endif
        wrtIMPACT_BC(ng)=.TRUE.
!
!  Convert TL initial condition from v-space to x-space.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile,Lini) SHARED(numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1,+1
            CALL tl_convolution (ng, TILE, Litl, Lweak, 2)
            CALL tl_variability (ng, TILE, Litl, Lweak)
#  ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, TILE, Lini, Litl)
#  endif
          END DO
        END DO
!$OMP END PARALLEL DO
!
!  Run tangent linear model for the assimilation period, t=t0 to t1.
!  Read and process the 4DVAR observations.
!
        IF (Master) THEN
          WRITE (stdout,10) 'TL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        TL_LOOP4 : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#  ifdef SOLVE3D
          CALL tl_main3d (ng)
#  else
          CALL tl_main2d (ng)
#  endif
          IF (exit_flag.ne.NoError) RETURN

        END DO TL_LOOP4
# endif

#endif

      END DO NEST_LOOP
!
 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/,15x,                &
     &        'adjoint forcing time range: ',f12.4,' - ',f12.4 ,/)
 30   FORMAT (/,' Out of range adjoint forcing time, ',a,f12.4,/,       &
     &        ' It must be between ',f12.4,' and ',f12.4)

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear and adjoint models      !
!  execution.                                                          !
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
      integer :: ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      DO ng=1,Ngrids
        IF (LwrtRST(ng).and.(exit_flag.eq.1)) THEN
          IF (Master) WRITE (stdout,10)
 10       FORMAT (/,' Blowing-up: Saving latest model state into ',     &
     &              ' RESTART file',/)
          IF (LcycleRST(ng).and.(NrecRST(ng).ge.2)) THEN
            tRSTindx(ng)=2
            LcycleRST(ng)=.FALSE.
          END IF
          blowup=exit_flag
          exit_flag=NoError
          CALL wrt_rst (ng)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks.  Close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL DO PRIVATE(thread) SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          CALL wclock_off (ng, iADM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE ROMS_finalize

      END MODULE ocean_control_mod
