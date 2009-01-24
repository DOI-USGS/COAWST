      MODULE ocean_control_mod
!
!svn $Id: symmetry.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Representer Driver:                                       !
!                                                                      !
!  This driver executes ROMS/TOMS weak constraint 4DVAR inner loop and !
!  checks the symmetry of the H R R' H' operator,  where R' and H' are !
!  transpose of R and H, respectively. The  R' H'  term is computed by !
!  integrating the adjoint model backwards while the  H R  is computed !
!  integrating forward the tangent linear model.                       !
!                                                                      !
!  It controls the initialization,  time-stepping, and finalization of !
!  the model execution following ESMF conventions:                     !
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
      USE mod_scalars
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first

      integer, intent(in), optional :: MyCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.

      integer :: STDrec, ng, thread

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
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
!
!  Read in background/model error standard deviation factors and
!  spatial convolution diffusion coefficients.
!  
        STDrec=1
        DO ng=1,Ngrids
          CALL get_state (ng, 6, 6, STDname(ng), STDrec, 1)
          IF (exit_flag.ne.NoError) RETURN
        END DO

      END IF

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Ted)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS weak constraint 4DVAR inner loop. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE ad_variability_mod, ONLY : ad_variability
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
      USE impulse_mod, ONLY : impulse
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
      USE normalization_mod, ONLY : normalization
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod, ONLY: tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      logical :: BOUNDED_TL, add
      logical :: Lweak, outer_impulse

      integer :: i, j, my_iic, ng, subs, tile, thread
      integer :: ADrec, Lstate, Nrec, rec
      integer :: IperAD, JperAD, KperAD, ivarAD
      integer :: IoutTL, JoutTL, KoutTL, ivarTL
#ifdef DISTRIBUTE
      integer :: Istr, Iend, Jstr, Jend
#endif
#ifdef BALANCE_OPERATOR
      integer :: Lbck = 1
#endif
      integer, allocatable :: StateVar(:)

      real(r8), allocatable :: R(:,:), Rerr(:,:)

      character (len=20) :: frmt
!
!-----------------------------------------------------------------------
!  Run model for all nested grids, if any.
!-----------------------------------------------------------------------
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  Run nonlinear model and compute basic state tracjectory.
!-----------------------------------------------------------------------
!
!  Initialize and set nonlinear model initial conditions.
!
        wrtNLmod(ng)=.FALSE.
        wrtRPmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Run nonlinear model and compute basic state trajectory.
!
        IF (Master) THEN
          WRITE (stdout,10) 'NL', ntstart(ng), ntend(ng)
        END IF

        time(ng)=time(ng)-dt(ng)

        NL_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

          iic(ng)=my_iic
#ifdef SOLVE3D
          CALL main3d (ng)
#else
          CALL main2d (ng)
#endif
          IF (exit_flag.ne.NoError) RETURN

        END DO NL_LOOP
!
!  Set forward file trajectory.
!
        FWDname(ng)=HISname(ng)
        ncFWDid(ng)=ncHISid(ng)
!
!-----------------------------------------------------------------------
!  Compute or read in background-error correlations normalization
!  factors.
!-----------------------------------------------------------------------
!
!  If computing, write out factors to NetCDF. This is an expensive
!  computation and needs to be computed once for a particular
!  application grid.
!
        IF (LwrtNRM(ng)) THEN
          CALL def_norm (ng)
          IF (exit_flag.ne.NoError) RETURN
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL normalization (ng, TILE, 2)
            END DO
          END DO
!$OMP END PARALLEL DO
          LdefNRM(ng)=.FALSE.
          LwrtNRM(ng)=.FALSE.
        ELSE
          tNRMindx(ng)=1
          CALL get_state (ng, 5, 5, NRMname(ng), tNRMindx(ng), 1)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define TLM impulse forcing NetCDF file.
!
        LdefTLF(ng)=.TRUE.
        CALL def_impulse (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!=======================================================================
!  Perturb the adjoint state variable, one at the time, with a delta
!  function at the specified perturbation point. Integrate the adjoint
!  model backwards and the force the tangent linear model with the
!  result. Then, store tangent linear solution at the requested point
!  for each state variable.
!=======================================================================
!  
!  Determine variables to perturb.
!
        Lold(ng)=1
        Ipass=1
        ERstr=1
#ifdef SOLVE3D
        ERend=NstateVar(ng)-2
        allocate ( StateVar(NstateVar(ng)-2) )
        StateVar(1)=isFsur
        StateVar(2)=isUvel
        StateVar(3)=isVvel
        DO i=1,NT(ng)
          StateVar(3+i)=isTvar(i)
        END DO
        Lstate=MAXVAL(StateVar)
#else
        ERend=NstateVar(ng)
        allocate ( StateVar(NstateVar(ng)) )
        StateVar(1)=isFsur
        StateVar(2)=isUbar
        StateVar(3)=isVbar
        Lstate=3
#endif
        allocate ( R(Lstate,Lstate) )
        allocate ( Rerr(Lstate,Lstate) )
!
!  Initialize sample representer matrix to report.
!
        R(1:Lstate,1:Lstate)=0.0_r8
!
!  Set point to perturb and point to report.
!
        ivarTL=INT(user(1))
        ivarAD=INT(user(2))
        IoutTL=INT(user(3))
        IperAD=INT(user(4))
        JoutTL=INT(user(5))
        JperAD=INT(user(6))
        KoutTL=INT(user(7))
        KperAD=INT(user(8))
!
        PERT_LOOP : DO Nrun=ERstr,ERend
          user(2)=StateVar(Nrun)          
!
!-----------------------------------------------------------------------
!  Time-step the adjoint model.
!-----------------------------------------------------------------------
!
!  Perturb the adjoint model at specified point.
!
          ADmodel=.TRUE.
          TLmodel=.FALSE.
          CALL ad_initial (ng)
          IF (exit_flag.ne.NoError) RETURN
          NrecADJ(ng)=0
          tADJindx(ng)=0
!
!  Time-step adjoint model backwards forced with current PSI vector.
!
          IF (Master) THEN
            WRITE (stdout,10) 'AD', ntstart(ng), ntend(ng)
          END IF

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
!-----------------------------------------------------------------------
!  Convolve adjoint trajectory with model-error covariance and convert
!  to impulse forcing.
!-----------------------------------------------------------------------
!
          Nrec=NrecADJ(ng)
          NrecADJ(ng)=0
          tADJindx(ng)=0
          LwrtState2d(ng)=.TRUE.
!
!  Clear adjoint state arrays.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile), SHARED(numthreads)
          DO thread=0,numthreads-1
            subs=NtileX(ng)*NtileE(ng)/numthreads
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, iADM)
            END DO
          END DO
!$OMF END PARALLEL DO

#ifdef BALANCE_OPERATOR
!
!  Read background state.
!
          CALL get_state (ng, iNLM, 9, FWDname(ng), Lbck, Lbck)
          IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Proccess each time record of current adjoint solution in ADJname.
!
          DO rec=1,Nrec
!
!  Set switch to scale model error covariace with background error
!  covariance factor Cfscale(:).
!
          IF (rec.eq.Nrec) THEN
            Lweak=.FALSE.
          ELSE
            Lweak=.TRUE.
          END IF
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
            ADrec=rec
            CALL get_state (ng, iTLM, 4, ADJname(ng), ADrec, Lold(ng))
            IF (exit_flag.ne.NoError) RETURN
!
!  Load interior solution, read above, into adjoint state arrays. 
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the 
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
            add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1
                CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
#ifdef BALANCE_OPERATOR
                CALL ad_balance (ng, TILE, Lbck, Lold(ng))
#endif
                CALL ad_variability (ng, TILE, Lold(ng), Lweak)
                CALL ad_convolution (ng, TILE, Lold(ng), 2)
                CALL initialize_ocean (ng, TILE, iTLM)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for output purposes.
!
            add=.FALSE.
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile)                          &
!$OMP&            SHARED(inner,add,numthreads)
            DO thread=0,numthreads-1
              subs=NtileX(ng)*NtileE(ng)/numthreads
              DO tile=subs*thread,subs*(thread+1)-1,+1
                CALL load_ADtoTL (ng, TILE, Lold(ng), Lold(ng), add)
                CALL tl_convolution (ng, TILE, Lold(ng), 2)
                CALL tl_variability (ng, TILE, Lold(ng), Lweak)
#ifdef BALANCE_OPERATOR
                CALL tl_balance (ng, TILE, Lbck, Lold(ng))
#endif
                CALL load_TLtoAD (ng, TILE, Lold(ng), Lold(ng), add)
              END DO
            END DO
!$OMP END PARALLEL DO
!
!  Overwrite ADJname history NetCDF file with convolved adjoint
!  solution.
!
            kstp(ng)=Lold(ng)
#ifdef SOLVE3D
            nstp(ng)=Lold(ng)
#endif
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
!
!  Convert convolved adjoint solution to impulse forcing. Write out
!  impulse forcing into TLFname NetCDF file. To facilitate the forcing
!  by the TLM and RPM, the forcing is process and written in
!  increasing time coordinates.
!
          tTLFindx(ng)=0
          outer_impulse=.FALSE.
#ifdef DISTRIBUTE
          tile=MyRank
#else
          tile=-1
#endif
          CALL impulse (ng, tile, iADM, outer_impulse, ADJname(ng))
!
!-----------------------------------------------------------------------
!  Integrate tangent linear model forced by the convolved adjoint
!  trajectory (impulse forcing) to compute R_n * PSI at observation
!  points.
!-----------------------------------------------------------------------
!
!  Initialize tangent linear model from rest. The initial contribution
!  from the adjoint model will be added as impulse forcing in
!  "tl_forcing".
!
          wrtTLmod(ng)=.FALSE.
          ADmodel=.FALSE.
          TLmodel=.FALSE.
          CALL tl_initial (ng)
          IF (exit_flag.ne.NoError) RETURN
!
!  Run tangent linear model forward and force with convolved adjoint
!  trajectory impulses. Compute R_n * PSI at observation points which
!  are used in the conjugate gradient algorithm.
!
          IF (Master) THEN
            WRITE (stdout,10) 'TL', ntstart(ng), ntend(ng)
          END IF

          time(ng)=time(ng)-dt(ng)

          TL_LOOP : DO my_iic=ntstart(ng),ntend(ng)+1

            iic(ng)=my_iic
#ifdef SOLVE3D
            CALL tl_main3d (ng)
#else
            CALL tl_main2d (ng)
#endif
            IF (exit_flag.ne.NoError) RETURN

          END DO TL_LOOP
!
!  Extract solution at requested points.
!
#ifdef DISTRIBUTE
          Istr=BOUNDS(ng)%Istr(MyRank)
          Iend=BOUNDS(ng)%Iend(MyRank)
          Jstr=BOUNDS(ng)%Jstr(MyRank)
          Jend=BOUNDS(ng)%Jend(MyRank)
          BOUNDED_TL=((Istr.le.IoutTL).and.(IoutTL.le.Iend)).and.       &
     &               ((Jstr.le.JoutTL).and.(JoutTL.le.Jend))
#else
          BOUNDED_TL=.TRUE.
#endif
          R(1,Nrun)=OCEAN(ng)%tl_zeta(IoutTL,JoutTL,knew(ng))
#ifndef SOLVE3D
          R(2,Nrun)=OCEAN(ng)%tl_ubar(IoutTL,JoutTL,knew(ng))
          R(3,Nrun)=OCEAN(ng)%tl_vbar(IoutTL,JoutTL,knew(ng))
#else
          R(2,Nrun)=OCEAN(ng)%tl_u(IoutTL,JoutTL,KoutTL,nstp(ng))
          R(3,Nrun)=OCEAN(ng)%tl_v(IoutTL,JoutTL,KoutTL,nstp(ng))
          DO i=1,NT(ng)
            R(i+3,Nrun)=OCEAN(ng)%tl_t(IoutTL,JoutTL,KoutTL,nstp(ng),i)
          END DO           
#endif

        END DO PERT_LOOP
!
!  Report sampled representer matrix and report symmetry.
!
#ifdef DISTRIBUTE
        DO i=1,NstateVar(ng)
          CALL mp_bcastf (ng, iTLM, R(:,i), NstateVar(ng))
        END DO
#endif
        IF (Master) THEN
          WRITE (stdout,20) 'Representer Matrix Symmetry Test: ',       &
     &                      'Perturbing Point: ',                       &
     &                      IperAD, JperAD, KperAD,                     &
     &                      'Sampling   Point: ',                       &
     &                      IoutTL, JoutTL, KoutTL
          WRITE (stdout,30) 'Sampled Representer Matrix: '
          IF (Lstate.lt.10) THEN
            WRITE (frmt,'(i1,a)') Lstate, '(1x,1p,e14.7,0p)'
          ELSE
            WRITE (frmt,'(i2,a)') Lstate, '(1x,1p,e14.7,0p)'
          END IF
          DO i=1,Lstate
            DO j=1,Lstate
              Rerr(i,j)=R(i,j)-R(j,i)
            END DO
          END DO
          DO i=1,Lstate
            WRITE (stdout,frmt) (R(i,j),j=1,Lstate)
          END DO
          WRITE (stdout,30) 'Representer Matrix Symmetry Error: '
          DO i=1,Lstate
            WRITE (stdout,frmt) (Rerr(i,j),j=1,Lstate)
          END DO
        END IF        

      END DO NEST_LOOP

 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 20   FORMAT (/,1x,a,/,/,3x,a,' i = ',i4.4,' j = ',i4.4,' k = ',i4.4,   &
     &                 /,3x,a,' i = ',i4.4,' j = ',i4.4,' k = ',i4.4)
 30   FORMAT (/,1x,a,/)

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS nonlinear model execution.        !
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

      END MODULE ocean_control_mod
