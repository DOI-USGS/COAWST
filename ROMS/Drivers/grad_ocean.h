      MODULE ocean_control_mod
!
!svn $Id: grad_ocean.h 652 2008-07-24 23:20:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Tangent Linear and Adjoint Models Gradient Test:          !
!                                                                      !
!  This driver is used to check the structure for variational data     !
!  assimilation using the gradient test:                               !
!                                                                      !
!  Denote the state vector as "s".  The cost function is the given     !
!  by J(s). Suppose we perturb "s" by "pds" where "p" is a scalar.     !
!  Then, using Taylor expansion to first-order we have:                !
!                                                                      !
!             J(s+pds) = J(s) + Transpose[grad(J)](pds)                !
!                                                                      !
!  Consider the functions:                                             !
!                                                                      !
!             g(p) = [J(s+pds) - J(s)] / Transpose[grad(J)](pds)       !
!                                                                      !
!  and                                                                 !
!                                                                      !
!             h(p) = [g(p) - 1] / p                                    !
!                                                                      !
!  As "p" goes to zero, we require g(p) to go to unity and h(p) to     !
!  go to a constant. These consitions will be satisfied if we have     !
!  the correct tangent linear and adjoint models.  Practically, we     !
!  define "ds" to be a steepest descent direction given by grad(J).    !
!  We get J(s) from the nonlinear model,  grad(J) from the adjoint     !
!  and J(s+pds) from the tangent linear model.                         !
!                                                                      !
!  The subroutines in this driver control the initialization, time-    !
!  stepping,  and  finalization of ROMS/TOMS  model following ESMF     !
!  conventions:                                                        !
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
!
!  Allocate and initialize observation arrays.
!
        CALL initialize_fourdvar
      END IF

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (Tstr, Tend)
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
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
# ifdef GRADIENT_CHECK
      USE dotproduct_mod, ONLY : ad_dotproduct
# endif
!
!  Imported variable declarations
!
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      integer :: IniRec, i, ic, my_iic, ng, subs, tile, thread

      real(r8) :: gp, hp, p

      character (len=80), dimension(MT+4) :: Pvars
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids
!
!-----------------------------------------------------------------------
!  Run nonlinear and adjoint models.
!-----------------------------------------------------------------------
!
!  Initialize relevant parameters.
!
        Lold(ng)=1
        Lnew(ng)=1
        Nrun=1
        Ipass=1
        ERstr=1
#ifdef SOLVE3D
        ERend=NstateVar(ng)-1
#else
        ERend=NstateVar(ng)+1
#endif
        ig1count=0
        IniRec=1
!
!  Initialize nonlinear model with first guess initial conditions.
!
        CALL initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Run nonlinear model. Extact and store nonlinear model values at
!  observation locations.
!
        IF (Master) THEN
          WRITE (stdout,20) 'NL', ntstart(ng), ntend(ng)
        END IF

        wrtNLmod(ng)=.TRUE.
        wrtTLmod(ng)=.FALSE.

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

        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.TRUE.
!
!  Save and Report cost function between nonlinear model and
!  observations.
!
        DO i=0,NstateVar(ng)
          FOURDVAR(ng)%CostFunOld(i)=FOURDVAR(ng)%CostFun(i)
        END DO
        IF (Master) THEN        
          WRITE (stdout,30) FOURDVAR(ng)%CostFunOld(0)
          DO i=1,NstateVar(ng)
            IF (FOURDVAR(ng)%CostFunOld(i).gt.0.0_r8) THEN
              IF (i.eq.1) THEN
                WRITE (stdout,40) FOURDVAR(ng)%CostFunOld(i),           &
     &                            TRIM(Vname(1,idSvar(i)))
              ELSE
                WRITE (stdout,50) FOURDVAR(ng)%CostFunOld(i),           &
     &                            TRIM(Vname(1,idSvar(i)))
              END IF
            END IF
          END DO
        END IF
!
!  Initialize the adjoint model from rest.
!
        CALL ad_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Time-step adjoint model: Compute model state gradient, GRAD(J).
!  Force the adjoint model with the adjoint misfit between nonlinear
!  model and observations.
!
        IF (Master) THEN
          WRITE (stdout,20) 'AD', ntstart(ng), ntend(ng)
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
!  Perturb each tangent linear state variable using the steepest decent
!  direction of GRAD(J).
!-----------------------------------------------------------------------
!
!  Load adjoint solution.
!
        CALL get_state (ng, iADM, 3, ADJname(ng), tADJindx(ng),         &
     &                  Lnew(ng))
!
!  Compute adjoint solution dot product for scaling purposes.
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads)
        DO thread=0,numthreads-1
          subs=NtileX(ng)*NtileE(ng)/numthreads
          DO tile=subs*thread,subs*(thread+1)-1
            CALL ad_dotproduct (ng, TILE, Lnew(ng))
          END DO
        END DO
!$OMP END PARALLEL DO

#ifdef SOLVE3D
!
!  OUTER LOOP: First, Perturb all state variables (zeta, u, v, t) at
!  ==========  once. Then, perturb one state variable at the time.
!              Notice, that ubar and vbar are not perturbed. They are
!              computed by vertically integarting u and v.
!
#else
!
!  OUTER LOOP: First, perturb all state variables (zeta, ubar, vbar) at
!  ==========  once. Then, perturb one state variable at the time.
!
#endif
        OUTER_LOOP : DO outer=ERstr,ERend

          CALL get_state (ng, iNLM, 1, FWDname(ng), IniRec, Lnew(ng))
!
!  INNER LOOP: scale perturbation amplitude by selecting "p" scalar,
!  ==========  such that:
!                              p = 10 ** FLOAT(-inner) 
!
          INNER_LOOP : DO inner=1,Ninner
!
!  Initialize tangent linear with the steepest decent direction
!  (adjoint state, GRAD(J)) times the perturbation amplitude "p".
! 
            CALL tl_initial (ng)
            IF (exit_flag.ne.NoError) RETURN
!
!  Time-step tangent linear model:  Compute misfit cost function
!  between model (nonlinear + tangent linear) and observations.
!
            IF (Master) THEN
              WRITE (stdout,20) 'TL', ntstart(ng), ntend(ng)
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
!  Report cost function between model (nonlinear + tangent linear)
!  and observations.
!
            IF (Master) THEN
              WRITE (stdout,60) outer, inner, FOURDVAR(ng)%CostFun(0)
              DO i=1,NstateVar(ng)
                IF (FOURDVAR(ng)%CostFun(i).gt.0.0_r8) THEN
                  IF (i.eq.1) THEN
                    WRITE (stdout,70) outer, inner,                     &
     &                                FOURDVAR(ng)%CostFun(i),          &
     &                                TRIM(Vname(1,idSvar(i)))
                  ELSE
                    WRITE (stdout,80) outer, inner,                     &
     &                                FOURDVAR(ng)%CostFun(i),          &
     &                                TRIM(Vname(1,idSvar(i)))
                  END IF
                END IF
              END DO
            END IF
!
!  Compute gradient check functions.
!
            p=10.0_r8**REAL(-inner,r8)
            gp=(FOURDVAR(ng)%CostFun(0)-                                &
     &          FOURDVAR(ng)%CostFunOld(0))/DotProduct
            hp=(gp-1.0_r8)/p
            IF (Master) THEN
              WRITE (stdout,90) outer, inner, p, gp, hp, DotProduct
              ig1count=ig1count+1
              g1(ig1count)=p
              g2(ig1count)=gp
            END IF

            Nrun=Nrun+1

          END DO INNER_LOOP

        END DO OUTER_LOOP
!
!  Report gradient funtions summary
!
        IF (Master) THEN
#ifdef SOLVE3D
          Pvars(1)='zeta, u, v'
          Pvars(2)='zeta'
          Pvars(3)='u'
          Pvars(4)='v'
          DO i=1,NT(ng)
            WRITE(Pvars(1),100) TRIM(Pvars(1)),                         &
     &                          TRIM(Vname(1,idSvar(isTvar(i))))
            WRITE(Pvars(i+4),110) TRIM(Vname(1,idSvar(isTvar(i))))
          END DO
#else
          Pvars(1)='zeta, ubar, vbar'
          Pvars(2)='zeta'
          Pvars(3)='ubar'
          Pvars(4)='vbar'
#endif 
          WRITE (stdout,120)                                            &
     &      'ADM Test - Gradient Functions Summary: p, g, h, error'
          ic=1
          WRITE (stdout,130)
          WRITE (stdout,140) TRIM(Pvars(ic))
          DO i=1,ig1count
            WRITE (stdout,150) i, g1(i), g2(i), (g2(i)-1.0_r8)/g1(i),   &
     &                         100.0_r8*ABS(g2(i)-1.0_r8)
            IF (MOD(i,Ninner).eq.0) THEN
              ic=ic+1
              WRITE (stdout,130)
              IF (ic.le.ERend) THEN
                WRITE (stdout,140) TRIM(Pvars(ic))
              END IF
            END IF
          END DO      
        END IF

      END DO NEST_LOOP
!
 20   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 30   FORMAT (/,' Nonlinear Model Cost Function = ',1p,e21.14)
 40   FORMAT (' --------------- ','cost function = ',1p,e21.14,2x,a)
 50   FORMAT (17x,'cost function = ',1p,e21.14,2x,a)
 60   FORMAT (/,' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,           &
     &        'Cost Function = ',1p,e21.14)
 70   FORMAT (' --------------- ','(',i4.4,',',i4.4,')',3x,             &
     &        'cost function = ',1p,e21.14,2x,a)
 80   FORMAT (17x,'(',i4.4,',',i4.4,')',3x,'cost function = ',          &
     &        1p,e21.14,2x,a)
 90   FORMAT (/,' (Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,           &
     &        'Gradient Function, p = ',1p,e21.14,0p,/,31x,             &
     &        'Gradient Function, g = ',1p,e21.14,0p,/,31x,             &
     &        'Gradient Function, h = ',1p,e21.14,0p,/,31x,             &
     &        'Gradient dot product = ',1p,e21.14,0p)
100   FORMAT (a,', ',a)
110   FORMAT (a)
120   FORMAT (/,a,/)
130   FORMAT (79('-'))
140   FORMAT (3x,'Perturbed state variable(s): ',a,/)
150   FORMAT (i4,2x,1p,e8.1,2(2x,1p,e21.14,0p),2x,f15.10,' %')

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
      integer :: ng, thread
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into record 3.
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
