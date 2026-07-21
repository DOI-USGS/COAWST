#include "cppdefs.h"
      MODULE convolve_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Multi-scale, Background-Error Covariace modeling:                   !
!                                                                      !
!  This module simulates the propagation effects of the background-    !
!  error covariance matrix (B) for each variable in the data           !
!  assimilation control vector by convolving generalized implicit      !
!  pseudo-diffusion tangent-linear and adjoint operators.              !
!                                                                      !
!  In the multiscale approach, B is expressed as a linear combination  !
!  of distinct spatial scales, from large to small. This formulation   !
!  allows for the representation of both broad structures and fine-    !
!  scale features, and reduces scale aliasing in the data assimilation !
!  cost function  (Weaver et al., 2016). It uses implicit horizontal   !
!  pseudo-diffusion operators implemented via Conjugate Gradient (CG)  !
!  and Chebyshev Iterations (CI). Error correlations are considered    !
!  separable in the horizontal and vertical directions. The vertical   !
!  diffusion operator is also implicit.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     driver       Calling driver (string)                             !
!     model        Initial conditions model to process (iTLM or iRPM)  !
!     outLoop      Current 4D-Var outer loop                           !
!     innLoop      Current 4D-Var inner loop                           !
!     Rbck         NLM background record to process                    !
!     Rini         NLM initial condition record to process             !
!     Rold         State vector old minimization time index            !
!     Rnew         State vector new minimization time index            !
!     Rec1         TLM initial conditions record to write              !
!     Rec2         RPM initial conditions record to write              !
!     Lposterior   Switch to process posterior error covariance        !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Weaver, A. and P. Courtier, 2001: Correlation modeling on the     !
!      sphere using a generalized diffusion equation, Q.J.R. Meteorol. !
!      Soc, 127, 1815-1846, doi:10.1002/qj.49712757518.                !
!                                                                      !
!    Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation  !
!      operators based on an implicitly formulated diffusion equation  !
!      solved with the Chebyshev iteration, Q.J.R. Meteorol. Soc.,     !
!      142, 455-471, doi:10.1002/qj.2664.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef RPCG
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
#ifdef ADJUST_BOUNDARY
      USE mod_boundary,       ONLY : initialize_boundary
#endif
      USE mod_forces,         ONLY : initialize_forces
      USE mod_ocean,          ONLY : initialize_ocean
!
#ifdef BALANCE_OPERATOR
      USE ad_balance_mod,     ONLY : ad_balance
#endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE sum_multi_B_mod
      USE ad_variability_mod, ONLY : ad_variability
      USE ad_wrt_his_mod,     ONLY : ad_wrt_his
#ifdef RPCG
      USE comp_Jb0_mod,       ONLY : aug_oper
#endif
      USE get_state_mod,      ONLY : get_state
      USE ini_adjust_mod,     ONLY : load_ADtoTL
      USE ini_adjust_mod,     ONLY : load_TLtoAD
#if defined ARRAY_MODES            || defined R4DVAR || \
     defined R4DVAR_ANA_SENSITIVITY
      USE ini_adjust_mod,     ONLY : rp_ini_adjust
#elif defined RBL4DVAR  || defined RBL4DVAR_ANA_SENSITIVITY
      USE ini_adjust_mod,     ONLY : ini_adjust
#endif
#if defined ARRAY_MODES            || defined R4DVAR || \
     defined R4DVAR_ANA_SENSITIVITY
      USE rp_wrt_ini_mod,     ONLY : rp_wrt_ini
#endif
      USE strings_mod,        ONLY : FoundError
#ifdef RPCG
      USE sum_grad_mod,       ONLY : sum_grad
#endif
#ifdef TIME_CONV
      USE time_corr_mod,      ONLY : time_corr
#endif
#ifdef BALANCE_OPERATOR
      USE tl_balance_mod,     ONLY : tl_balance
#endif
      USE tl_convolution_mod, ONLY : tl_convolution
      USE tl_variability_mod, ONLY : tl_variability
      USE tl_wrt_ini_mod,     ONLY : tl_wrt_ini
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
      USE wrt_ini_mod,        ONLY : wrt_frc_AD
#endif
      USE wrt_ini_mod,        ONLY : wrt_ini
#ifdef POSTERIOR_ERROR_I
      USE wrt_hessian_mod,    ONLY : wrt_hessian
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined ADJUST_BOUNDARY
      USE wrt_ini_mod,        ONLY : wrt_frc_AD
#endif
      USE wrt_ini_mod,        ONLY : wrt_ini
!
      implicit none
!
      PUBLIC  :: convolve
      PUBLIC  :: error_covariance
#ifdef SP4DVAR
      PUBLIC  :: saddlec
#endif
      PRIVATE
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE convolve (driver, Rini, Rold, Rnew)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Rini
      integer, intent(in) :: Rold(Ngrids)
      integer, intent(in) :: Rnew(Ngrids)
!
      character (len=*), intent(in) :: driver
!
!  Local variable declarations.
!
      logical :: Lweak, add
!
      integer, parameter :: ifac = 2             ! squared-root operator
      integer :: IniRec
      integer :: ng, tile, ns
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Specify spatial error covariance on the TLM control vector, at time
!  index "Rold", via convolutions of a pseudo-diffusion operator.  The
!  convolved control vector is loaded into time index "Rnew" and index
!  "Rold" is preserved.
!-----------------------------------------------------------------------

#ifdef BALANCE_OPERATOR
!
!  Read NLM initial condition in readiness for the balance operator.
!
      IniRec=Rini
      DO ng=1,Ngrids
        CALL get_state (ng, iNLM, 2, INI(ng), IniRec, IniRec)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        nrhs(ng)=Rini
      END DO
#endif
!
!  Load "Rold" time index TLM control vector into ADM state arrays.
!  Apply balance operator, if activated. Then, scale control vector
!  with error covariance standard deviations. Next, convolve resulting
!  vector with the squared-root adjoint diffusion operator. Notice
!  that the spatial convolution is only done for half of the diffusion
!  steps (squared-root filter, ifac=2).
!
      Lweak=.FALSE.
      add=.FALSE.
!
      NESTED_GRID_LOOP : DO ng=1,Ngrids
        MULTISCALE_LOOP : DO ns=1,Nscale(ng)
#ifdef PROFILE
          CALL wclock_on (ng, iADM, 82, __LINE__, MyFile)
#endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
#ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, tile, Rini, Rold(ng))
#endif
#ifndef DIRAC
            CALL ad_variability (ng, tile, Rold(ng), Lweak)
#endif
            CALL ad_convolution (ng, tile, ns, Rold(ng), Lweak, ifac)
          END DO
#ifdef PROFILE
          CALL wclock_off (ng, iADM, 82, __LINE__, MyFile)
#endif
!
!  Since we wish to preserve what is in tl_var(Rold), load resulting
!  filtered solution from above, ad_var(Rold), into tl_var(Rnew).
!  Then, convolve with the squared-root (half of steps, ifac=2) tangent
!  linear diffusion operator. Next, scale results with error
!  covariance standard deviations. Apply balance operator, if
!  activated.
!
          add=.FALSE.
#ifdef PROFILE
          CALL wclock_on (ng, iTLM, 82, __LINE__, MyFile)
#endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_ADtoTL (ng, tile, Rold(ng), Rnew(ng), add)
            CALL tl_convolution (ng, tile, ns, Rnew(ng), Lweak, ifac)
#ifndef DIRAC
            CALL tl_variability (ng, tile, Rnew(ng), Lweak)
#endif
#ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, tile, Rini, Rnew(ng))
#endif
          END DO
#ifdef PROFILE
          CALL wclock_off (ng, iTLM, 82, __LINE__, MyFile)
#endif
!
!  Accumulate the implicit diffusion operator solutions which are in
!  tl_var(Rnew) in ad_var(Rnew).
!
#ifdef PROFILE
          CALL wclock_on (ng, iTLM, 82, __LINE__, MyFile)
#endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL sum_multi_B (ng, tile, ns, Rnew(ng), Rnew(ng))
          END DO
        END DO MULTISCALE_LOOP
      END DO NESTED_GRID_LOOP
!
!  Copy accumulated ad_var(Rnew) to tl_var(Rnew).
!
      DO ng=1,Ngrids
#ifdef PROFILE
        CALL wclock_on (ng, iTLM, 82, __LINE__, MyFile)
#endif
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_ADtoTL (ng, tile, Rnew(ng), Rnew(ng), add)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE convolve
!
!***********************************************************************
      SUBROUTINE error_covariance (model, driver, outLoop, innLoop,     &
     &                             Rbck, Rini, Rold, Rnew,              &
     &                             Rec1, Rec2, Lposterior)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lposterior
!
      integer, intent(in) :: model, outLoop, innLoop
      integer, intent(in) :: Rbck, Rini, Rec1, Rec2
      integer, intent(in) :: Rold(Ngrids)
      integer, intent(in) :: Rnew(Ngrids)
!
      character (len=*), intent(in) :: driver
!
!  Local variable declarations.
!
      logical :: Lweak, add
!
      integer, parameter :: ifac = 2             ! squared-root operator
      integer :: ADrec, BckRec, Fcount, IniRec
      integer :: i, irec, ng, tile, ns
#ifdef RPCG
      integer :: LTLM1, LTLM2, Rec5, jrec, nADrec
#endif
      integer, dimension(Ngrids) :: Nrec
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", error_covariance"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Convolve adjoint trajectory.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        IF ((outLoop.lt.0).and.(innLoop.lt.0)) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Convolving Adjoint Trajectory...')
        ELSE
          WRITE (stdout,20) outLoop, innLoop
 20       FORMAT (/,' Convolving Adjoint Trajectory: Outer = ',i3.3,    &
     &            ' Inner = ',i3.3)
        END IF
      END IF
!
!  Clear adjoint state arrays.
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
        END DO
      END DO
!
!  Get number of records in adjoint NetCDF and initialize indices.
!
      DO ng=1,Ngrids
        Fcount=ADM(ng)%Fcount
        Nrec(ng)=ADM(ng)%Nrec(Fcount)
        ADM(ng)%Nrec(Fcount)=0
        ADM(ng)%Rindex=0
        LwrtState2d(ng)=.TRUE.
        LwrtTime(ng)=.FALSE.
      END DO
!
!  Read in adjoint control vector to convolve from NetCDF file
!  ADM(ng)%name for time record Nrec, which correspond to the
!  initial conditions time. If adjusting open boundaries and/or
!  surface forcing, only record Nrec is read since it is the
!  only record for which adjoint forcing arrays are complete.
!
!  Notice that since routine "get_state" loads data into the
!  ghost points, the adjoint contol vector is read into the
!  tangent linear state arrays by using iTLM instead of iADM
!  in the calling arguments.
!
      NESTED_GRID_LOOP : DO ng=1,Ngrids
        MULTISCALE_LOOP : DO ns=1,Nscale(ng)

          ADrec=Nrec(ng)
          FrcRec(ng)=Nrec(ng)
          CALL get_state (ng, iTLM, 4, ADM(ng), ADrec, Rold(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef BALANCE_OPERATOR
!
!  Read NLM initial condition in readiness for the balance
!  operator.
!
          IniRec=Rini
          CALL get_state (ng, iNLM, 2, INI(ng), IniRec, IniRec)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          nrhs(ng)=Rini
#endif
!
!  Load ADM control vector read above from TLM into ADM state arrays.
!  Then, multiply control vector by the error covariance standard
!  deviations. Next, convolve resulting adjoint solution with the
!  squared-root adjoint diffusion operator to impose apriori error
!  hypothesis. Notice that the spatial convolution are only done
!  for half of the diffusion steps (squared-root filter, ifac=2).
!  Clear tangent linear state arrays when done.
!
          Lweak=.FALSE.
          add=.FALSE.
!
#ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
!
#endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
#ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, tile, Rini, Rold(ng))
#endif
            CALL ad_variability (ng, tile, Rold(ng), Lweak)
            CALL ad_convolution (ng, tile, ns, Rold(ng), Lweak, ifac)
            CALL initialize_ocean (ng, tile, iTLM)
            CALL initialize_forces (ng, tile, iTLM)
#ifdef ADJUST_BOUNDARY
            CALL initialize_boundary (ng, tile, iTLM)
#endif
          END DO
#ifdef PROFILE
!
          CALL wclock_off (ng, model, 82, __LINE__, MyFile)
#endif
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps, ifac=2) tangent
!  diffusion operator. Then, multiply result with its corresponding
!  linear error covariance standard deviations. Since the convolved
!  solution is in the adjoint state arrays, first copy to tangent
!  linear state arrays including the ghosts points.
!
          Lweak=.FALSE.
          add=.FALSE.
!
#ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
!
#endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_ADtoTL (ng, tile, Rold(ng), Rold(ng), add)
            CALL tl_convolution (ng, tile, ns, Rold(ng), Lweak, ifac)
            CALL tl_variability (ng, tile, Rold(ng), Lweak)
#ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, tile, Rini, Rold(ng))
#endif
          END DO
!
!  Accumulate the implicit diffusion operator solutions which are in
!  tl_var(Rold) in ad_var(Rnew).
!
#ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
#endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL sum_multi_B (ng, tile, ns, Rold(ng), Rnew(ng))
          END DO

#ifdef POSTERIOR_ERROR_I
!
!  If computing the analysis error covariance matrix, copy TLM back
!  into ADM so that it can be written to Hessian NetCDF file.
!
          DO ng=1,Ngrids
            IF (Lposterior) THEN
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
              END DO
            END IF
          END DO
#endif
#ifdef PROFILE
!
          CALL wclock_off (ng, model, 82, __LINE__, MyFile)
#endif

        END DO MULTISCALE_LOOP
      END DO NESTED_GRID_LOOP
!
!  Copy accumulated ad_var(Rnew) to tl_var(Rold).
!
      DO ng=1,Ngrids
#ifdef PROFILE
        CALL wclock_on (ng, model, 82, __LINE__, MyFile)
#endif
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_ADtoTL (ng, tile, Rnew(ng), Rold(ng), add)
        END DO
      END DO

#if defined ARRAY_MODES            || defined R4DVAR || \
     defined R4DVAR_ANA_SENSITIVITY
!
!  Copy back to adjoint state arrays when done with the convolution.
!  Compute representer model initial conditions by adding convolved
!  adjoint solution to the reference nonlinear state (INI(ng)%name,
!  record Rbck).
!
      IF ((model.eq.iRPM).and.                                          &
     &    (INDEX(driver,'r4dvar').ne.0)) THEN
        BckRec=Rbck
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 9, INI(ng), BckRec, Rnew(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
        add=.FALSE.
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            CALL rp_ini_adjust (ng, tile, Rnew(ng), Rold(ng))
          END DO
        END DO
      END IF

#elif defined RBL4DVAR || defined RBL4DVAR_ANA_SENSITIVITY
# ifndef RPCG
!
!  Copy back to adjoint state arrays when done with the convolution.
!  Compute nonlinear model initial conditions by adding convolved
!  adjoint solution to the reference nonlinear state (INI(ng)%name,
!  record Rbck).
!
      IF ((model.eq.iNLM).and.                                          &
     &    (INDEX(driver,'rbl4dvar').ne.0)) THEN
        BckRec=Rbck
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 9, INI(ng), BckRec, Rnew(ng))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
!
        add=.FALSE.
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            CALL ini_adjust (ng, tile, Rold(ng), Rnew(ng))
          END DO
        END DO
      END IF
# endif
#endif
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing adjustments for next inner loop into
!  ITL(ng)%name (record Rec1). The tangent model initial
!  conditions are set to the convolved adjoint solution.
!
      IF (model.eq.iTLM) THEN
        DO ng=1,Ngrids
#ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rold(ng), Rec1)
#else
          CALL tl_wrt_ini (ng, -1, Rold(ng), Rec1)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
#if defined ARRAY_MODES            || defined R4DVAR || \
     defined R4DVAR_ANA_SENSITIVITY
      ELSE IF (model.eq.iRPM) THEN
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          CALL rp_wrt_ini (ng, MyRank, Rold(ng), Rec2)
# else
          CALL rp_wrt_ini (ng, -1, Rold(ng), Rec2)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
#elif (defined RBL4DVAR || defined RBL4DVAR_ANA_SENSITIVITY) && \
       !defined RPCG
      ELSE IF (model.eq.iNLM) THEN
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          CALL wrt_ini (ng, MyRank, Rnew(ng))
# else
          CALL wrt_ini (ng, -1, Rnew(ng))
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
      defined ADJUST_BOUNDARY
#  ifdef DISTRIBUTE
          CALL wrt_frc_AD (ng, MyRank, Rold(ng), INI(ng)%Rindex)
#  else
          CALL wrt_frc_AD (ng, -1, Rold(ng), INI(ng)%Rindex)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
        END DO
#endif
      END IF

#ifdef RPCG
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing and obc adjustments for next outer
!  loop into ITLname (record Rec2). The tangent model initial
!  conditions are set to the convolved adjoint solution.
!
      IF (model.eq.iNLM) THEN
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          CALL tl_wrt_ini (ng, MyRank, Rold(ng), Rec2)
# else
          CALL tl_wrt_ini (ng, -1, Rold(ng), Rec2)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#endif

#ifdef POSTERIOR_ERROR_I
!
!  If weak constraint, convolve records 2:Nrec in ADM(ng)%name and
!  Write convolved adjoint solution into Hessian NetCDF file for use
!  later.
!
      IF (Lposterior.and.(inner.ne.0)) THEN
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          CALL wrt_hessian (ng, MyRank, Rold(ng), Rold(ng))
# else
          CALL wrt_hessian (ng, -1, Rold(ng), Rold(ng))
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
#endif

#ifdef RPCG
!
!  Now compute the new B^-1(x(k)-x(k-1)) term for the weak constraint
!  forcing terms. Records 6+nADrec/2 to nADrec+5 contain the sums so
!  far. This is ONLY done in the outer-loop not the inner-loops, and
!  is controlled by the flag "LaugWeak".
!
      IF (LaugWeak) THEN
        LTLM1=1
        LTLM2=2
        Rec5=5
        DO ng=1,Ngrids
          IF (nADJ(ng).lt.ntimes(ng)) THEN
            nADrec=2*(1+ntimes(ng)/nADJ(ng))
          ELSE
            nADrec=0
          END IF
          DO irec=1,nADrec/2
            jrec=Rec5+nADrec/2+irec
!
!  First add the augmented piece which is computed from the previous
!  sum.
!
            CALL get_state (ng, iTLM, 4, ITL(ng), jrec, LTLM1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL aug_oper (ng, tile, LTLM1, LTLM2)
            END DO
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!
!  Need to read adjoint netcdf file in reverse.
!
            ADrec=(Nrec(ng)-1)-(irec-1)
            CALL get_state (ng, iTLM, 4, ADM(ng), ADrec, LTLM1)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!
!  Write the current sum of adjoint solutions into record jrec of the
!  ITL file.
!
# ifdef DISTRIBUTE
            CALL tl_wrt_ini (ng, MyRank, LTLM2, jrec)
# else
            CALL tl_wrt_ini (ng, -1, LTLM2, jrec)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          END DO
        END DO
      END IF
#endif

#ifdef TIME_CONV_NOT_YET
!
!  The multiscale option IS NOT CODED YET FOR TIME CONVOLUTION.
!
!  Apply the factorized space-time model error covariance
!  matrix of the form ATA' where A is the square-root of the
!  model error covariance matrix, and T is the time-correlation
!  matrix.
!
!  If weak constraint, convolve records 2-Nrec in ADM(ng)%name
!  and impose model error covariance. NOTE: We will not use the
!  convolved forcing increments generated here since these arrays
!  do not contain the complete solution and are redundant.
!  AMM: We might want to get rid of these unwanted records to
!  avoid any confusion in the future.
!
      DO ng=1,Ngrids
        IF (Nrec(ng).gt.3) THEN
# ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
# endif
          DO irec=1,Nrec(ng)-1
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
            ADrec=irec
            CALL get_state (ng, iTLM, 4, ADM(ng), ADrec, Rold(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Load interior solution, read above, into adjoint state arrays.
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter, ifac=2).
!  Clear tangent linear state arrays when done.
!
            Lweak=.TRUE.
            add=.FALSE.
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
              CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
              CALL ad_variability (ng, tile, Rold(ng), Lweak)
              CALL ad_convolution (ng, tile, Rold(ng), Lweak, ifac)
              CALL initialize_ocean (ng, tile, iTLM)
              CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
              CALL initialize_boundary (ng, tile, iTLM)
# endif
            END DO
!
!  Overwrite ADM(ng)%name history NetCDF file with convolved adjoint
!  solution.
!
            kstp(ng)=Rold(ng)
# ifdef SOLVE3D
            nstp(ng)=Rold(ng)
# endif
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
          LwrtTime(ng)=.TRUE.
# ifdef PROFILE
        CALL wclock_off (ng, model, 82, __LINE__, MyFile)
# endif
        END IF
      END DO
!
!  Apply the time correlation matrix to the newly convolved adjoint
!  fields.
!
      DO ng=1,Ngrids
        TLF(ng)%Rindex=0
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL time_corr (ng, tile, iADM, ADM(ng)%IOtype, ADM(ng)%name)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END DO
!
      DO ng=1,Ngrids
        IF (Nrec(ng).gt.3) THEN
# ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
# endif
          ADM(ng)%Rindex=0
          LwrtState2d(ng)=.TRUE.
          DO irec=Nrec(ng)-1,1,-1
!
!  Read the latest TLF solution. We need to read the TLF file in reverse
!  order since the final solution is written into the adjoint netcdf
!  file which will be rewritten in ascending order of time in the TLF
!  file by wrt_impulse.
!  NOTE: model=6 here so as to only read the state variables.
!
            ADrec=irec
            CALL get_state (ng, 6, 6, TLF(ng), ADrec, Rold(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
            Lweak=.TRUE.
            add=.FALSE.
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_convolution (ng, tile, Rold(ng), Lweak, ifac)
              CALL tl_variability (ng, tile, Rold(ng), Lweak)
# ifdef BALANCE_OPERATOR
              CALL tl_balance (ng, tile, Rini, Rold(ng))
# endif
              CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            END DO
!
!  Overwrite ADM(ng)%name history NetCDF file with convolved adjoint
!  solution.
!
            kstp(ng)=Rold(ng)
# ifdef SOLVE3D
            nstp(ng)=Rold(ng)
# endif
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
          LwrtTime(ng)=.TRUE.
# ifdef PROFILE
          CALL wclock_off (ng, model, 82, __LINE__, MyFile)
# endif
        END IF
      END DO

#else

!
!  If weak constraint, convolve records 2-Nrec in ADM(ng)%name
!  and impose model error covariance. NOTE: We will not use the
!  convolved forcing increments generated here since these arrays
!  do not contain the complete solution and are redundant.
!  AMM: We might want to get rid of these unwanted records to
!  avoid any confusion in the future.
!
      WEAK_NESTED_LOOP : DO ng=1,Ngrids
        WEAK_SWITCH : IF (Nrec(ng).gt.3) THEN

# ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
# endif

          WEAK_REC_LOOP : DO irec=1,Nrec(ng)-1
            WEAK_MULTISCALE_LOOP : DO ns=1,Nscale(ng)

!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
              ADrec=irec
              CALL get_state (ng, iTLM, 4, ADM(ng), ADrec, Rold(ng))
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
!
!  Load interior solution, read above, into adjoint state arrays.
!  Then, multiply adjoint solution by the background-error standard
!  deviations. Next, convolve resulting adjoint solution with the
!  squared-root adjoint diffusion operator which impose the model-error
!  spatial correlations. Notice that the spatial convolution is only
!  done for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
              Lweak=.TRUE.
              add=.FALSE.
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
                CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
                CALL ad_variability (ng, tile, Rold(ng), Lweak)
                CALL ad_convolution (ng, tile, ns, Rold(ng),Lweak, ifac)
                CALL initialize_ocean (ng, tile, iTLM)
                CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
                CALL initialize_boundary (ng, tile, iTLM)
# endif
              END DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  background-error standard deviations.  Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points. Copy back to adjoint state
!  arrays when done with the convolution for output purposes.
!
              Lweak=.TRUE.
              add=.FALSE.
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL load_ADtoTL (ng, tile, Rold(ng), Rold(ng), add)
                CALL tl_convolution (ng, tile, ns, Rold(ng),Lweak, ifac)
                CALL tl_variability (ng, tile, Rold(ng), Lweak)
# ifdef BALANCE_OPERATOR
                CALL tl_balance (ng, tile, Rini, Rold(ng))
# endif
              END DO
!
!  Accumulate the implicit diffusion operator solutions which are in
!  tl_var(Rnew) in ad_var(Rnew).
!
# ifdef PROFILE
              CALL wclock_on (ng, model, 82, __LINE__, MyFile)
# endif
              DO tile=first_tile(ng),last_tile(ng),+1
                CALL sum_multi_B (ng, tile, ns, Rold(ng), Rnew(ng))
              END DO

            END DO WEAK_MULTISCALE_LOOP
!
!  Overwrite ADM(ng)%name history NetCDF file with final solution.
!
            kstp(ng)=Rnew(ng)
# ifdef SOLVE3D
            nstp(ng)=Rnew(ng)
# endif
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          END DO WEAK_REC_LOOP
!
          LwrtState2d(ng)=.FALSE.
          LwrtTime(ng)=.TRUE.
# ifdef PROFILE
          CALL wclock_off (ng, model, 82, __LINE__, MyFile)
# endif
        END IF WEAK_SWITCH

      END DO WEAK_NESTED_LOOP
#endif
!
      RETURN
      END SUBROUTINE error_covariance

#ifdef SP4DVAR
!
!***********************************************************************
      SUBROUTINE saddlec (driver, Lselect, Rini, Rold, Rnew)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Rini
      integer, intent(in) :: Rold(Ngrids)
      integer, intent(in) :: Rnew(Ngrids)
!
      logical, intent(in) :: Lselect
!
      character (len=*), intent(in) :: driver
!
!  Local variable declarations.
!
      logical :: Lweak, add
!
      integer, parameter :: ifac = 2             ! squared-root operator
      integer :: IniRec
      integer :: ng, tile, ns
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", saddlec"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Specify spatial error covariance on the TLM control vector, at time
!  index "Rold", via convolutions of a pseudo-diffusion operator.  The
!  convolved control vector is loaded into time index "Rnew" and index
!  "Rold" is preserved.
!-----------------------------------------------------------------------

# ifdef BALANCE_OPERATOR
!
!  Read NLM initial condition in readiness for the balance operator.
!
      IniRec=Rini
      DO ng=1,Ngrids
        CALL get_state (ng, iNLM, 2, INI(ng), IniRec, IniRec)
        IF (exit_flag.ne.NoError) RETURN
        nrhs(ng)=Rini
      END DO
# endif
!
!  Load "Rold" time index TLM control vector into ADM state arrays.
!  Apply balance operator, if activated. Then, scale control vector
!  with error covariance standard deviations. Next, convolve resulting
!  vector with the squared-root adjoint diffusion operator. Notice
!  that the spatial convolution is only done for half of the diffusion
!  steps (squared-root filter, ifac=2).
!
      Lweak=Lselect
      add=.FALSE.
!
      NESTED_GRID_LOOP : DO ng=1,Ngrids
        MULTISCALE_LOOP : DO ns=1,Nscale(ng)

# ifdef PROFILE
          CALL wclock_on (ng, iADM, 82, __LINE__, MyFile)
# endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
            CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
            CALL ad_variability (ng, tile, Rold(ng), Lweak)
            CALL ad_convolution (ng, tile, ns, Rold(ng), Lweak, ifac)
            CALL initialize_ocean (ng, tile, iTLM)       ! AMM important
            CALL initialize_forces (ng, tile, iTLM)      ! AMM important
          END DO
# ifdef PROFILE
          CALL wclock_off (ng, iADM, 82, __LINE__, MyFile)
# endif
!
!  Since we wish to preserve what is in tl_var(Rold), load resulting
!  filtered solution from above, ad_var(Rold), into tl_var(Rnew).
!  Then, convolve with the squared-root (half of steps, ifac=2) tangent
!  linear diffusion operator. Next, scale results with error
!  covariance standard deviations. Apply balance operator, if
!  activated.
!
          add=.FALSE.
# ifdef PROFILE
          CALL wclock_on (ng, iTLM, 82, __LINE__, MyFile)
# endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_ADtoTL (ng, tile, Rold(ng), Rnew(ng), add)
            CALL tl_convolution (ng, tile, ns, Rnew(ng), Lweak, ifac)
            CALL tl_variability (ng, tile, Rnew(ng), Lweak)
# ifdef BALANCE_OPERATOR
            CALL tl_balance (ng, tile, Rini, Rnew(ng))
# endif
          END DO
# ifdef PROFILE
          CALL wclock_off (ng, iTLM, 82, __LINE__, MyFile)
# endif
!
!  Accumulate the implicit diffusion operator solutions which are in
!  tl_var(Rnew) in ad_var(Rnew).
!
# ifdef PROFILE
          CALL wclock_on (ng, model, 82, __LINE__, MyFile)
# endif
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL sum_multi_B (ng, tile, ns, Rnew(ng), Rnew(ng))
          END DO

        END DO MULTISCALE_LOOP
      END DO NESTED_GRID_LOOP
!
!  Copy accumulated ad_var(Rnew) to tl_var(Rnew).
!
      DO ng=1,Ngrids
# ifdef PROFILE
        CALL wclock_on (ng, model, 82, __LINE__, MyFile)
# endif
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_ADtoTL (ng, tile, Rnew(ng), Rnew(ng), add)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE saddlec
#endif
      END MODULE convolve_mod
