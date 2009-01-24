      MODULE ocean_control_mod
!
!svn $Id: so_semi_ocean.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Stochastic Optimals, Seminorm Estimation Driver:          !
!                                                                      !
!  This driver computes the eigenvectors of the stochastic optimals    !
!  operator with respect the seminorm of the chosen functional. The    !
!  stochastic  optimals  provide information about the influence of    !
!  stochastic  variations  (biases) in  ocean forcing.  They can be    !
!  used to build forecast ensembles.                                   !
!                                                                      !
!  These routines control the initialization,  time-stepping,  and     !
!  finalization of  ROMS/TOMS  model following ESMF conventions:       !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Moore, A.M. et al., 2004: A comprehensive ocean prediction and    !
!      analysis system based on the tangent linear and adjoint of a    !
!      regional ocean model, Ocean Modelling, 7, 227-258.              !
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
      USE mod_iounits
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

      integer :: ng, thread

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
            CALL wclock_on (ng, iNLM, 0)
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
        CALL inp_par (iNLM)
        IF (exit_flag.ne.NoError) RETURN
!
!  Allocate and initialize modules variables.
!
        CALL mod_arrays (allocate_vars)

      END IF

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
!
!=======================================================================
!                                                                      !
!  This routine computes the eigenvectors of the stochastic optimal    !
!  matrix, defined as:                                                 !
!                                                                      !
!       S = INTEGRAL [ R'(t,To) X R(t,To) dt ]   from To to T          !
!                                                                      !
!  where T is the  forecast interval from the initial state at To,     !
!  R(t,To) is  the  forward  tangent  propagator of the linearized     !
!  dynamical model,  R'(t,To) is the adjoint of  R(t,To),  and the     !
!  matrix X defines the norm of interest. Here, X is a seminorm of     !
!  the chosen functional.                                              !
!                                                                      !
!  The above operator is symmetric and the  ARPACK library is used     !
!  to select eigenvectors and eigenvalues:                             !
!                                                                      !
!  Lehoucq, R.B., D.C. Sorensen, and C. Yang, 1997:  ARPACK user's     !
!    guide:  solution  of  large  scale  eigenvalue  problems with     !
!    implicit restarted Arnoldi Methods, Rice University, 140p.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_storage
!
      USE propagator_mod
      USE packing_mod, ONLY : so_unpack
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: ITERATE, LwrtGST

      integer :: NconvRitz, i, iter, ng, status, varid
      integer :: thread, subs, tile
      integer :: start(4), total(4)

#ifdef DISTRIBUTE
      real(r8), external :: pdnorm2
#else
      real(r8), external :: dnrm2
#endif

      character (len=55) :: string
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
!  Initialize tangent linear for all grids first in order to compute
!  the size of the state vector, Nstate.  This size is computed in
!  routine "wpoints".
!
      DO ng=1,Ngrids
        CALL ad_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Currently, only non-nested applications are considered.  Otherwise,
!  a different structure for mod_storage is needed. 
!
      NEST_LOOP : DO ng=1,Ngrids

        IF (ng.eq.1) THEN
          CALL allocate_storage (ng)
        END IF
!
!  Initialize various parameters.
!
        Nrun=0

        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LwrtPER(ng)=.FALSE.
        LwrtGST=.TRUE.
        LcycleTLM(ng)=.FALSE.
        LcycleADJ(ng)=.FALSE.
!
!-----------------------------------------------------------------------
!  Implicit Restarted Arnoldi Method (IRAM) for the computation of
!  optimal perturbation Ritz eigenfunctions.
!-----------------------------------------------------------------------
!
        ITERATE=.TRUE.

        Lrvec=.TRUE.              ! Compute Ritz vectors
        ido=0                     ! reverse communication flag
        bmat='I'                  ! standard eigenvalue problem
        which='LM'                ! compute NEV largest eigenvalues
        howmany='A'               ! compute NEV Ritz vectors
        info=0                    ! random initial residual vector
        iparam(1)=1               ! exact shifts
        iparam(3)=MaxIterGST      ! maximum number of Arnoldi iterations
        iparam(4)=1               ! block size in the recurrence
        iparam(7)=1               ! type of eigenproblem being solved
!
!  ARPACK debugging parameters.
!
        logfil=stdout             ! output logical unit
        ndigit=-3                 ! number of decimal digits
        msaupd=1                  ! iterations, timings, Ritz
        msaup2=1                  ! norms, Ritz values
        msaitr=0
        mseigt=0
        msapps=0
        msgets=0
        mseupd=0
!
!  Determine size of the eigenproblem (Nsize) and size of work space
!  array SworkL (LworkL).
!
        Nsize=Nend(ng)-Nstr(ng)+1
        Nconv=0

#ifdef CHECKPOINTING
!
!  If restart, read in checkpointing data GST restart NetCDF file.
!  Otherwise, create checkpointing restart NetCDF file.
!
        IF (LrstGST) THEN
          CALL get_gst (ng, iTLM)
          ido=-2
          laup2(1)=.FALSE.        ! cnorm
          laup2(2)=.FALSE.        ! getv0
          laup2(3)=.FALSE.        ! initv
          laup2(4)=.FALSE.        ! update
          laup2(5)=.TRUE.         ! ushift
        ELSE
          CALL def_gst (ng, iTLM)
        END IF
        IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Iterate until either convergence or maximum iterations has been
!  exceeded.
!
        iter=0
!
        ITER_LOOP : DO WHILE (ITERATE)
          iter=iter+1
!
!  Reverse communication interface.
!
#ifdef PROFILE
          CALL wclock_on (ng, iTLM, 38)
#endif
#ifdef DISTRIBUTE
          CALL pdsaupd (OCN_COMM_WORLD,                                 &
     &                  ido, bmat, Nsize, which, NEV, Ritz_tol,         &
     &                  resid(Nstr(ng)), NCV, Bvec(Nstr(ng),1),         &
     &                  Nsize, iparam, ipntr,                           &
     &                  SworkD, SworkL, LworkL, info)
#else
          CALL dsaupd (ido, bmat, Nsize, which, NEV, Ritz_tol,          &
     &                 resid, NCV, Bvec, Nsize, iparam, ipntr,          &
     &                 SworkD, SworkL, LworkL, info)
#endif
#ifdef PROFILE
          CALL wclock_off (ng, iTLM, 38)
#endif
#ifdef CHECKPOINTING2
!
!  If appropriate, write out check point data into GST restart NetCDF
!  file. Notice that the restart data is always saved if MaxIterGST
!  is reached without convergence. It is also saved when convergence
!  is achieved (ido=99).
!
          IF ((MOD(iter,nGST).eq.0).or.(iter.ge.MaxIterGST).or.         &
              (ido.eq.99)) THEN
            CALL wrt_gst (ng, iTLM)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#endif
!
!  Terminate computations if maximum number of iterations is reached.
!  This will faciliate splitting the analysis in several computational
!  cycles using the restart option.
!
          IF ((iter.ge.MaxIterGST).and.(ido.ne.99)) THEN
            ITERATE=.FALSE.
            EXIT ITER_LOOP
          END IF
!
!  Perform matrix-vector operation:  R`(t,0)XR(0,t)u
!
          IF (ABS(ido).eq.1) THEN
            NrecADJ(ng)=0
            NrecTLM(ng)=0
            tADJindx(ng)=0
            tTLMindx(ng)=0
            Nconv=iaup2(4)
            CALL propagator (ng, Nstr(ng), Nend(ng),                    &
     &                       SworkD(ipntr(1):), SworkD(ipntr(2):))
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            IF (info.ne.0) THEN
              IF (Master) THEN
                CALL IRAM_error (info, string)
                WRITE (stdout,10) 'DSAUPD', TRIM(string),               &
     &                            ', info = ', info
              END IF
              RETURN
            ELSE
!
!  Compute Ritz vectors. (The only choice left is IDO=99).
!
              IF (Master) THEN
                WRITE (stdout,20) 'Number of converged Ritz values:',   &
     &                            iparam(5)
                WRITE (stdout,20) 'Number of Arnoldi iterations taken:',&
     &                            iparam(3)
              END IF
#ifdef PROFILE
              CALL wclock_on (ng, iTLM, 38)
#endif
#ifdef DISTRIBUTE
              CALL pdseupd (OCN_COMM_WORLD,                             &
     &                      Lrvec, howmany, select,                     &
     &                      RvalueR, Rvector(Nstr(ng),1), Nsize,        &
     &                      sigmaR, bmat, Nsize, which, NEV, Ritz_tol,  &
     &                      resid(Nstr(ng)), NCV, Bvec(Nstr(ng),1),     &
     &                      Nsize, iparam, ipntr,                       &
     &                      SworkD, SworkL, LworkL, info)
#else
              CALL dseupd (Lrvec, howmany, select,                      &
     &                     RvalueR, Rvector, Nsize,                     &
     &                     sigmaR, bmat, Nsize, which, NEV, Ritz_tol,   &
     &                     resid, NCV, Bvec, Nsize, iparam,             &
     &                     ipntr, SworkD, SworkL, LworkL, info)
#endif
#ifdef PROFILE
              CALL wclock_off (ng, iTLM, 38)
#endif
              IF (info.ne.0) THEN
                IF (Master) THEN
                  CALL IRAM_error (info, string)
                  WRITE (stdout,10) 'DSEUPD', TRIM(string),             &
     &                              ', info = ', info
                END IF
                RETURN
              ELSE
!
!  The converged stochastic vectors are written into the nonlinear
!  history NetCDF file. Notice that only surface forcing variables
!  are defined.
!
                NrecHIS(ng)=0
                tHISindx(ng)=0
                ndefHIS(ng)=0
                LwrtHIS(ng)=.TRUE.
                LdefHIS(ng)=.TRUE.
                DO i=1,NV
                  Hout(i,ng)=.FALSE.
                END DO
                IF (SCALARS(ng)%SOstate(isUstr)) Hout(idUsms,ng)=.TRUE.
                IF (SCALARS(ng)%SOstate(isVstr)) Hout(idVsms,ng)=.TRUE.
#ifdef SOLVE3D
                DO i=1,NT(ng)
                  IF (SCALARS(ng)%SOstate(isTsur(i))) THEN
                    Hout(idTsur(i),ng)=.TRUE.
                  END IF
                END DO
#endif
                HISname(ng)=TLMname(ng)
                CALL def_his (ng, LdefHIS(ng))
                IF (exit_flag.ne.NoError) RETURN
!
!  Clear forcing arrays.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile) SHARED(ng,numthreads)
                DO thread=0,numthreads-1
                  subs=NtileX(ng)*NtileE(ng)/numthreads
                  DO tile=subs*thread,subs*(thread+1)-1
                    CALL initialize_forces (ng, TILE, 0)
                  END DO
                END DO
!$OMP END PARALLEL DO
!
!  Check residuals (Euclidean norm) and Ritz values.
!
                NconvRitz=iparam(5)
                DO i=1,NconvRitz
                  CALL propagator (ng, Nstr(ng), Nend(ng),              &
     &                             Rvector(Nstr(ng):,i), SworkD)
                  IF (exit_flag.ne.NoError) RETURN
!
!  Unpack surface forcing eigenvectors from Rvector and write into
!  nonlinear history file.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads,Nstr,Nend,state,ad_state)
                  DO thread=0,numthreads-1
                    subs=NtileX(ng)*NtileE(ng)/numthreads
                    DO tile=subs*thread,subs*(thread+1)-1,+1
                      CALL so_unpack (ng, TILE, Nstr(ng), Nend(ng),     &
     &                                Rvector(Nstr(ng):,i))
                    END DO
                  END DO
!$OMP END PARALLEL DO
                  CALL wrt_his (ng)
                  IF (exit_flag.ne.NoError) RETURN
!
!  Compute Euclidean norm.
!
                  CALL daxpy (Nsize, -RvalueR(i), Rvector(Nstr(ng):,i), &
     &                        1, SworkD, 1)
#ifdef DISTRIBUTE
                  norm(i)=pdnorm2(OCN_COMM_WORLD, Nsize, SworkD, 1)
#else
                  norm(i)=dnrm2(Nstate(ng), SworkD, 1)
#endif
                  IF (Master) THEN
                    WRITE (stdout,30) i, norm(i), RvalueR(i)
                  END IF
!
!  Write out Ritz eigenvalues and Ritz eigenvector Euclidean norm to
!  NetCDF file.
!
                  IF (OutThread) THEN
                    start(1)=i
                    total(1)=1
                    status=nf90_inq_varid(ncHISid(ng), 'Ritz_rvalue',   &
     &                                    varid)
                    status=nf90_put_var(ncHISid(ng), varid, RvalueR(i:),&
     &                                  start, total)
                    status=nf90_inq_varid(ncHISid(ng), 'Ritz_norm',     &
     &                                    varid)
                    status=nf90_put_var(ncHISid(ng), varid, norm(i:),   &
     &                                  start, total)
                    IF (i.eq.1) THEN
                      status=nf90_inq_varid(ncHISid(ng), 'SO_trace',    &
     &                                      varid)
                      status=nf90_put_var(ncHISid(ng), varid,           &
     &                                    TRnorm(ng))
                    END IF
                  END IF
                END DO
              END IF
            END IF
            ITERATE=.FALSE.
          END IF

#ifdef CHECKPOINTING
!
!  If appropriate, write out check point data into GST restart NetCDF
!  file. Notice that the restart data is always saved if MaxIterGST
!  is reached without convergence. It is also saved when convergence
!  is achieved (ido=99).
!
          IF ((MOD(iter,nGST).eq.0).or.(iter.ge.MaxIterGST).or.         &
              ((ido.eq.99).and.LwrtGST)) THEN
            CALL wrt_gst (ng, iTLM)
            IF (exit_flag.ne.NoError) RETURN
            IF (ido.eq.99) LwrtGST=.FALSE.
          END IF
#endif

        END DO ITER_LOOP

      END DO NEST_LOOP
!
 10   FORMAT (/,1x,'Error in ',a,1x,a,a,1x,i5,/)
 20   FORMAT (/,a,1x,i2,/)
 30   FORMAT (4x,i4.4,'-th residual ',1p,e14.6,0p,                      &
     &        ' Corresponding to Ritz value ',1pe14.6)

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
          CALL wclock_off (ng, iNLM, 0)
        END DO
!$OMP END PARALLEL DO
      END DO
!
!  Close IO files.
!
      CALL close_io

      RETURN
      END SUBROUTINE ROMS_finalize

      SUBROUTINE IRAM_error (info, string)
!
!=======================================================================
!                                                                      !
!  This routine decodes internal error messages from the Implicit      !
!  Restarted Arnoldi Method (IRAM) for the computation of optimal      !
!  perturbation Ritz eigenfunctions.                                   !
!                                                                      !
!=======================================================================
!
!
!  imported variable declarations.
!
      integer, intent(in) :: info

      character (len=*), intent(out) :: string
!
!-----------------------------------------------------------------------
!  Decode error message from IRAM.
!-----------------------------------------------------------------------
!
      IF (info.eq.0)  THEN
        string='Normal exit                                            '
      ELSE IF (info.eq.1) THEN
        string='Maximum number of iterations taken                     '
      ELSE IF (info.eq.3) THEN
        string='No shifts could be applied during an IRAM cycle        '
      ELSE IF (info.eq.-1) THEN
        string='Nstate must be positive                                '
      ELSE IF (info.eq.-2) THEN
        string='NEV must be positive                                   '
      ELSE IF (info.eq.-3) THEN
        string='NCV must be greater NEV and less than or equal Nstate  '
      ELSE IF (info.eq.-4) THEN
        string='Maximum number of iterations must be greater than zero '
      ELSE IF (info.eq.-5) THEN
        string='WHICH must be one of LM, SM, LA, SA or BE              '
      ELSE IF (info.eq.-6) THEN
        string='BMAT must be one of I or G                             '
      ELSE IF (info.eq.-7) THEN
        string='Length of private work array SworkL is not sufficient  '
      ELSE IF (info.eq.-8) THEN
        string='Error in DSTEQR in the eigenvalue calculation          '
      ELSE IF (info.eq.-9) THEN
        string='Starting vector is zero                                '
      ELSE IF (info.eq.-10) THEN
        string='IPARAM(7) must be 1, 2, 3, 4, 5                        '
      ELSE IF (info.eq.-11) THEN
        string='IPARAM(7) = 1 and BMAT = G are incompatable            '
      ELSE IF (info.eq.-12) THEN
        string='IPARAM(1) must be equal to 0 or 1                      '
      ELSE IF (info.eq.-13) THEN
        string='NEV and WHICH = BE are incompatable                    '
      ELSE IF (info.eq.-14) THEN
        string='Did not find any eigenvalues to sufficient accuaracy   '
      ELSE IF (info.eq.-15) THEN
        string='HOWMANY must be one of A or S if RVEC = .TRUE.         '
      ELSE IF (info.eq.-16) THEN
        string='HOWMANY = S not yet implemented                        '
      ELSE IF (info.eq.-17) THEN
        string='Different count of converge Ritz values in DSEUPD      '
      ELSE IF (info.eq.-9999) THEN
        string='Could not build and Arnoldi factorization              '
      END IF

      RETURN
      END SUBROUTINE IRAM_error

      END MODULE ocean_control_mod
