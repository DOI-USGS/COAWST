      MODULE ocean_control_mod
!
!svn $Id: fsv_ocean.h 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Forcing Singular Vectors (FSV) Driver:                    !
!                                                                      !
!  This driver computes the forcing singular vectors of the propagator !
!  R(0,t) when the forcing is constant in time.  The solution is then: !
!                                                                      !
!      s(t) = M(t) * f                                                 !
!                                                                      !
!  where                                                               !
!                                                                      !
!      M(t) = integral[R(t',t) dt']   from t'=0 to t'=t                !
!                                                                      !
!  and f is the stochastic forcing constant in time.  The eigenvectors !
!  of transpose(M)XM  are the forcing singular vectors and can be used !
!  to generate ensembles of forecasts  associated  with the  different !
!  possible realizations of systematic errors in surface forcing.      !
!                                                                      !
!  These  routines  control  the  initialization,  time-stepping,  and !
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
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_storage

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
        CALL inp_par (iTLM)
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
            CALL wclock_on (ng, iTLM, 0)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL mod_arrays (allocate_vars)
!$OMP END PARALLEL

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
!  Initialize tangent linear for all grids first in order to compute
!  the size of the state vector, Nstate.  This size is computed in
!  routine "wpoints".
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#if defined BULK_FLUXES && defined NL_BULK_FLUXES
        BLK(ng)%name=FWD(ng)%name
#endif
!$OMP PARALLEL
        CALL tl_initial (ng)
!$OMP END PARALLEL
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Allocate arrays associated with Generalized Stability Theory (GST)
!  analysis.
!
      CALL allocate_storage
!
!  Initialize various IO flags.
!
      Nrun=0
      DO ng=1,Ngrids
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LwrtPER(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        LcycleADJ(ng)=.FALSE.
        nADJ(ng)=ntimes(ng)
        nTLM(ng)=ntimes(ng)
      END DO
!
!  Initialize ARPACK parameters.
!
      Lrvec=.TRUE.                ! Compute Ritz vectors
      bmat='I'                    ! standard eigenvalue problem
      which='LM'                  ! compute NEV largest eigenvalues
      howmany='A'                 ! compute NEV Ritz vectors
      DO ng=1,Ngrids
        ido(ng)=0                 ! reverse communication flag
        info(ng)=0                ! random initial residual vector
        iparam(1,ng)=1            ! exact shifts
        iparam(3,ng)=MaxIterGST   ! maximum number of Arnoldi iterations
        iparam(4,ng)=1            ! block size in the recurrence
        iparam(7,ng)=1            ! type of eigenproblem being solved
      END DO
!
!  ARPACK debugging parameters.
!
      logfil=stdout               ! output logical unit
      ndigit=-3                   ! number of decimal digits
      msaupd=1                    ! iterations, timings, Ritz
      msaup2=1                    ! norms, Ritz values
      msaitr=0
      mseigt=0
      msapps=0
      msgets=0
      mseupd=0
!
!  Determine size of the eigenproblem (Nsize) and size of work space
!  array SworkL (LworkL).
!
      DO ng=1,Ngrids
        Nconv(ng)=0
        Nsize(ng)=Nend(ng)-Nstr(ng)+1
      END DO

#ifdef CHECKPOINTING
!
!  If restart, read in check pointing data GST restart NetCDF file.
!  Otherwise, create check pointing restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (LrstGST) THEN
          CALL get_gst (ng, iTLM)
          ido(ng)=-2
        ELSE
          CALL def_gst (ng, iTLM)
        END IF
        IF (exit_flag.ne.NoError) RETURN
      END DO
#endif

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the forcing singular vectors of R(0,t) by a   !
!  singel integration of a perturbation  "u" forward intime with the   !
!  tangent linear model  over [0,t],  multiplication  of the  result   !
!  by "X",  followed by an  integration of the  result  backwards in   !
!  time with the  adjoint model over [t,0].  This  is  equivalmet to   !
!  the matrix-vector operation:                                        !
!                                                                      !
!       transpose[R(t,0)] X R(0,t) u                                   !
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
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE mod_storage
!
      USE propagator_mod
      USE packing_mod, ONLY : r_norm2
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: ITERATE
#ifdef CHECKPOINTING
      logical :: LwrtGST
#endif

      integer :: Fcount, Is, Ie, i, iter, ng
      integer :: NconvRitz(Ngrids)

      real(r8) :: Enorm

      real(r8), dimension(2) :: my_norm, my_value

      TYPE (T_GST), allocatable :: ad_state(:)
      TYPE (T_GST), allocatable :: state(:)

      character (len=55) :: string
!
!-----------------------------------------------------------------------
!  Implicit Restarted Arnoldi Method (IRAM) for the computation of
!  optimal perturbation Ritz eigenfunctions.
!-----------------------------------------------------------------------
!
!  Allocate nested grid pointers for state vectors.
!
      IF (.not.allocated(ad_state)) THEN
        allocate ( ad_state(Ngrids) )
      END IF
      IF (.not.allocated(state)) THEN
        allocate ( state(Ngrids) )
      END IF
!
!  Iterate until either convergence or maximum iterations has been
!  exceeded.
!
      iter=0
      ITERATE=.TRUE.
#ifdef CHECKPOINTING
      LwrtGST=.TRUE.
#endif
!
      ITER_LOOP : DO WHILE (ITERATE)
        iter=iter+1
!
!  Reverse communication interface.
!
        DO ng=1,Ngrids
#ifdef PROFILE
          CALL wclock_on (ng, iTLM, 38)
#endif
#ifdef DISTRIBUTE
          CALL pdsaupd (OCN_COMM_WORLD,                                 &
     &                  ido(ng), bmat, Nsize(ng), which, NEV,           &
     &                  Ritz_tol,                                       &
     &                  STORAGE(ng)%resid(Nstr(ng)), NCV,               &
     &                  STORAGE(ng)%Bvec(Nstr(ng),1), Nsize(ng),        &
     &                  iparam(1,ng), ipntr(1,ng),                      &
     &                  STORAGE(ng)%SworkD,                             &
     &                  SworkL(1,ng), LworkL, info(ng))
#else
          CALL dsaupd (ido(ng), bmat, Nsize(ng), which, NEV,            &
     &                 Ritz_tol,                                        &
     &                 STORAGE(ng)%resid, NCV,                          &
     &                 STORAGE(ng)%Bvec, Nsize(ng),                     &
     &                 iparam(1,ng), ipntr(1,ng),                       &
     &                 STORAGE(ng)%SworkD,                              &
     &                 SworkL(1,ng), LworkL, info(ng))
#endif
          Nconv(ng)=iaup2(4)
#ifdef PROFILE
          CALL wclock_off (ng, iTLM, 38)
#endif
#ifdef CHECKPOINTING
!
!  If appropriate, write out check point data into GST restart NetCDF
!  file. Notice that the restart data is always saved if MaxIterGST
!  is reached without convergence. It is also saved when convergence
!  is achieved (ido=99).
!
          IF ((MOD(iter,nGST).eq.0).or.(iter.ge.MaxIterGST).or.         &
     &        (ANY(ido.eq.99))) THEN
            CALL wrt_gst (ng, iTLM)
            IF (exit_flag.ne.NoError) RETURN
          END IF
#endif
        END DO
!
!  Terminate computations if maximum number of iterations is reached.
!  This will faciliate splitting the analysis in several computational
!  cycles using the restart option.
!
        IF ((iter.ge.MaxIterGST).and.ANY(ido.ne.99)) THEN
          ITERATE=.FALSE.
          EXIT ITER_LOOP
        END IF
!
!  Perform matrix-vector operation:  R`(t,0)XR(0,t)u
!
        IF (ANY(ABS(ido).eq.1)) THEN
          DO ng=1,Ngrids
            Fcount=ADM(ng)%Fcount
            ADM(ng)%Nrec(Fcount)=0
            Fcount=TLM(ng)%Fcount
            TLM(ng)%Nrec(Fcount)=0
            ADM(ng)%Rindex=0
            TLM(ng)%Rindex=0
          END DO
!
!  Set state vectors to process by the propagator via pointer
!  equivalence.
!
          DO ng=1,Ngrids
            IF (ASSOCIATED(state(ng)%vector)) THEN
              nullify (state(ng)%vector)
            END IF
            Is=ipntr(1,ng)
            Ie=Is+Nsize(ng)-1
            state(ng)%vector => STORAGE(ng)%SworkD(Is:Ie)

            IF (ASSOCIATED(ad_state(ng)%vector)) THEN
              nullify (ad_state(ng)%vector)
            END IF
            Is=ipntr(2,ng)
            Ie=Is+Nsize(ng)-1
            ad_state(ng)%vector => STORAGE(ng)%SworkD(Is:Ie)
          END DO

!$OMP PARALLEL
          CALL propagator (RunInterval, state, ad_state)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          IF (ANY(info.ne.0)) THEN
            DO ng=1,Ngrids
              IF (info(ng).ne.0) THEN
                IF (Master) THEN
                  CALL IRAM_error (info(ng), string)
                  WRITE (stdout,10) 'DSAUPD', TRIM(string),             &
     &                              ', info = ', info(ng)
                END IF
                RETURN
              END IF
            END DO
          ELSE
!
!  Compute Ritz vectors (the only choice left is IDO=99).  They are
!  generated in ARPACK in decreasing magnitude of its eigenvalue.
!  The most significant is first.
!
            DO ng=1,Ngrids
              NconvRitz(ng)=iparam(5,ng)
              IF (Master) THEN
                WRITE (stdout,20) 'Number of converged Ritz values:',   &
     &                            iparam(5,ng)
                WRITE (stdout,20) 'Number of Arnoldi iterations:',      &
     &                            iparam(3,ng)
              END IF
#ifdef PROFILE
              CALL wclock_on (ng, iTLM, 38)
#endif
#ifdef DISTRIBUTE
              CALL pdseupd (OCN_COMM_WORLD,                             &
     &                      Lrvec, howmany, select(1,ng),               &
     &                      RvalueR(1,ng),                              &
     &                      STORAGE(ng)%Rvector(Nstr(ng),1),            &
     &                      Nsize(ng), sigmaR,                          &
     &                      bmat, Nsize(ng), which, NEV, Ritz_tol,      &
     &                      STORAGE(ng)%resid(Nstr(ng)), NCV,           &
     &                      STORAGE(ng)%Bvec(Nstr(ng),1), Nsize(ng),    &
     &                      iparam(1,ng), ipntr(1,ng),                  &
     &                      STORAGE(ng)%SworkD,                         &
     &                      SworkL(1,ng), LworkL, info(ng))
#else
              CALL dseupd (Lrvec, howmany, select(1,ng),                &
     &                     RvalueR(1,ng),                               &
     &                     STORAGE(ng)%Rvector, Nsize(ng),              &
     &                     sigmaR, bmat, Nsize, which, NEV, Ritz_tol,   &
     &                     STORAGE(ng)%resid, NCV,                      &
     &                     STORAGE(ng)%Bvec, Nsize(ng),                 &
     &                     iparam(1,ng), ipntr(1,ng),                   &
     &                     STORAGE(ng)%SworkD,                          &
     &                     SworkL(1,ng), LworkL, info(ng))
#endif
#ifdef PROFILE
              CALL wclock_off (ng, iTLM, 38)
#endif
            END DO

            IF (ANY(info.ne.0)) THEN
              DO ng=1,Ngrids
                IF (info(ng).ne.0) THEN
                  IF (Master) THEN
                    CALL IRAM_error (info(ng), string)
                    WRITE (stdout,10) 'DSEUPD', TRIM(string),           &
     &                                ', info = ', info(ng)
                  END IF
                  RETURN
                END IF
              END DO
            ELSE
!
!  Activate writing of each eigenvector into the adjoint and tangent
!  linear history NetCDF files. The "ocean_time" is the eigenvector
!  number.
!
              Nrun=0

              DO i=1,MAXVAL(NconvRitz)
                DO ng=1,Ngrids
                  IF ((i.eq.1).or.LmultiGST) THEN
                    Fcount=ADM(ng)%Fcount
                    ADM(ng)%Nrec(Fcount)=0
                    Fcount=TLM(ng)%Fcount
                    TLM(ng)%Nrec(Fcount)=0
                    ADM(ng)%Rindex=0
                    TLM(ng)%Rindex=0
                  END IF
                  IF (LmultiGST) THEN
                    LdefADJ(ng)=.TRUE.
                    LdefTLM(ng)=.TRUE.
                    WRITE (ADM(ng)%name,30) TRIM(ADM(ng)%base), i
                    WRITE (TLM(ng)%name,30) TRIM(TLM(ng)%base), i
                  END IF
                END DO
!
!  Compute and write Ritz eigenvectors.
!
                DO ng=1,Ngrids
                  Is=Nstr(ng)
                  Ie=Nend(ng)
                  IF (ASSOCIATED(state(ng)%vector)) THEN
                    nullify (state(ng)%vector)
                  END IF

                  IF (ASSOCIATED(ad_state(ng)%vector)) THEN
                    nullify (ad_state(ng)%vector)
                  END IF
                  state(ng)%vector => STORAGE(ng)%Rvector(Is:Ie,i)
                  ad_state(ng)%vector => SworkR(Is:Ie)
                END DO

!$OMP PARALLEL
                CALL propagator (RunInterval, state, ad_state)
!$OMP END PARALLEL
                IF (exit_flag.ne.NoError) RETURN
!
                DO ng=1,Ngrids
                  CALL r_norm2 (ng, iTLM, Nstr(ng), Nend(ng),           &
     &                          -RvalueR(i,ng),                         &
     &                          state(ng)%vector,                       &
     &                          ad_state(ng)%vector, Enorm)
                  norm(i,ng)=Enorm
                  IF (Master) THEN
                    WRITE (stdout,40) i, norm(i,ng), RvalueR(i,ng), i
                  END IF
                END DO
!
!  Write out Ritz eigenvalues and Ritz eigenvector Euclidean norm
!  (residual) to NetCDF file(s).  Notice that we write the same value
!  twice in the TLM file for the initial and final perturbation of
!  the eigenvector.
!
                DO ng=1,Ngrids
                  SourceFile='fsv_ocean.h, ROMS_run'
                  my_norm(1)=norm(i,ng)
                  my_norm(2)=my_norm(1)
                  my_value(1)=RvalueR(i,ng)
                  my_value(2)=my_value(1)

                  IF (LwrtTLM(ng)) THEN
                    CALL netcdf_put_fvar (ng, iTLM, TLM(ng)%name,       &
     &                                    'Ritz_rvalue',                &
     &                                    my_value,                     &
     &                                    start = (/TLM(ng)%Rindex/),   &
     &                                    total = (/1/),                &
     &                                    ncid = TLM(ng)%ncid)
                    IF (exit_flag.ne. NoError) RETURN

                    CALL netcdf_put_fvar (ng, iTLM, TLM(ng)%name,       &
     &                                    'Ritz_norm',                  &
     &                                    my_norm,                      &
     &                                    start = (/TLM(ng)%Rindex/),   &
     &                                    total = (/1/),                &
     &                                    ncid = TLM(ng)%ncid)
                    IF (exit_flag.ne. NoError) RETURN

                    IF (LmultiGST) THEN
                      CALL netcdf_close (ng, iTLM, TLM(ng)%ncid,        &
     &                                   TLM(ng)%name)
                      IF (exit_flag.ne.NoError) RETURN
                    END IF
                  END IF
                  IF (LwrtADJ(ng)) THEN
                    CALL netcdf_put_fvar (ng, iADM, ADM(ng)%name,       &
     &                                    'Ritz_rvalue',                &
     &                                    RvalueR(i:,ng),               &
     &                                    start = (/ADM(ng)%Rindex/),   &
     &                                    total = (/1/),                &
     &                                    ncid = ADM(ng)%ncid)
                    IF (exit_flag.ne. NoError) RETURN

                    CALL netcdf_put_fvar (ng, iADM, ADM(ng)%name,       &
     &                                    'Ritz_norm', norm(i:,ng),     &
     &                                    start = (/ADM(ng)%Rindex/),   &
     &                                    total = (/1/),                &
     &                                    ncid = ADM(ng)%ncid)
                    IF (exit_flag.ne. NoError) RETURN

                    IF (LmultiGST) THEN
                      CALL netcdf_close (ng, iADM, ADM(ng)%ncid,        &
     &                                   ADM(ng)%name)
                      IF (exit_flag.ne.NoError) RETURN
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
          ITERATE=.FALSE.
        END IF

      END DO ITER_LOOP
!
 10   FORMAT (/,1x,'Error in ',a,1x,a,a,1x,i5,/)
 20   FORMAT (/,a,1x,i2,/)
 30   FORMAT (a,'_',i3.3,'.nc')
 40   FORMAT (1x,i4.4,'-th residual',1p,e14.6,0p,                       &
     &        '  Ritz value',1pe14.6,0p,2x,i4.4)

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
      integer :: Fcount, ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
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
 20     FORMAT (/,' Elapsed CPU time (seconds):',/)
      END IF

      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 0)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Close IO files.
!
      CALL close_out

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
