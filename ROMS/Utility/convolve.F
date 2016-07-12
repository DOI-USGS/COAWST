#include "cppdefs.h"
      MODULE convolve_mod

#ifdef FOUR_DVAR
!
!svn $Id: convolve.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!=======================================================================
!                                                                      !
!  These routines impose the error covariance matrices for the 4D-Var  !
!  control variables. The  error covariances matrices  are modeled by  !
!  the spatial convolutions of pseudo-diffusion equations.             !
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
!      sphere using a generalized diffusion equation, Q.J.R. Meteo.    !
!      Soc, 127, 1815-1846.                                            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: convolve
      PUBLIC  :: error_covariance

      CONTAINS
!
!***********************************************************************
      SUBROUTINE convolve (driver, Rini, Rold, Rnew)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
!
# ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
      USE tl_balance_mod, ONLY: tl_balance
# endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE tl_convolution_mod, ONLY : tl_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE tl_variability_mod, ONLY : tl_variability
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
!
!  Imported variable declarations.
!
      integer, intent(in) :: Rini
      integer, intent(in) :: Rold(Ngrids)
      integer, intent(in) :: Rnew(Ngrids)

      character (len=*), intent(in) :: driver
!
!  Local variable declarations.
!
      logical :: Lweak, add
      integer :: IniRec
      integer :: ng, tile
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
        CALL get_state (ng, iNLM, 2, INI(ng)%name, IniRec, IniRec)
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
!  steps (squared-root filter).
!
      Lweak=.FALSE.
      add=.FALSE.
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
          CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
          CALL ad_variability (ng, tile, Rold(ng), Lweak)
          CALL ad_convolution (ng, tile, Rold(ng), Lweak, 2)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Since we wish to preserve what is in tl_var(Rold), load resulting
!  filtered solution from above, ad_var(Rold), into tl_var(Rnew).
!  Then, convolve with the squared-root (half of steps) tangent
!  linear diffusion operator. Next, scale results with error
!  covariance standard deviations. Apply balance operator, if
!  activated.
!
      add=.FALSE.
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_ADtoTL (ng, tile, Rold(ng), Rnew(ng), add)
          CALL tl_convolution (ng, tile, Rnew(ng), Lweak, 2)
          CALL tl_variability (ng, tile, Rnew(ng), Lweak)
# ifdef BALANCE_OPERATOR
          CALL tl_balance (ng, tile, Rini, Rnew(ng))
# endif
        END DO
!$OMP END PARALLEL
      END DO

      RETURN
      END SUBROUTINE convolve

!
!***********************************************************************
      SUBROUTINE error_covariance (model, driver, outLoop, innLoop,     &
     &                             Rbck, Rini, Rold, Rnew,              &
     &                             Rec1, Rec2, Lposterior)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_stepping
# ifdef RPCG
      USE mod_fourdvar
# endif
!
# ifdef BALANCE_OPERATOR
      USE ad_balance_mod, ONLY: ad_balance
      USE tl_balance_mod, ONLY: tl_balance
# endif
      USE ad_convolution_mod, ONLY : ad_convolution
      USE tl_convolution_mod, ONLY : tl_convolution
      USE ad_variability_mod, ONLY : ad_variability
      USE tl_variability_mod, ONLY : tl_variability
      USE ini_adjust_mod, ONLY : load_ADtoTL
      USE ini_adjust_mod, ONLY : load_TLtoAD
# if defined ARRAY_MODES        || defined W4DVAR || \
     defined W4DVAR_SENSITIVITY
      USE ini_adjust_mod, ONLY : rp_ini_adjust
# elif defined W4DPSAS  || defined W4DPSAS_SENSITIVITY
      USE ini_adjust_mod, ONLY : ini_adjust
# endif
      USE mod_forces, ONLY : initialize_forces
      USE mod_ocean, ONLY : initialize_ocean
      USE mod_boundary, ONLY : initialize_boundary
# ifdef RPCG
      USE sum_grad_mod, ONLY : sum_grad
      USE comp_Jb0_mod, ONLY : aug_oper
# endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lposterior
      integer, intent(in) :: model, outLoop, innLoop
      integer, intent(in) :: Rbck, Rini, Rec1, Rec2
      integer, intent(in) :: Rold(Ngrids)
      integer, intent(in) :: Rnew(Ngrids)

      character (len=*), intent(in) :: driver
!
!  Local variable declarations.
!
      logical :: Lweak, add
      integer :: ADrec, BckRec, Fcount, IniRec
      integer :: i, irec, ng, tile
# ifdef RPCG
      integer :: LTLM1, LTLM2, Rec5, jrec, nADrec
# endif

      integer, dimension(Ngrids) :: Nrec
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
 20       FORMAT (/,' Convolving Adjoint Trajectory: Outer = ',i3.3,     &
     &            ' Inner = ',i3.3)
        END IF
      END IF
!
!  Clear adjoint state arrays.
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL initialize_ocean (ng, tile, iADM)
        END DO
!$OMP END PARALLEL
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
      DO ng=1,Ngrids
        ADrec=Nrec(ng)
        FrcRec(ng)=Nrec(ng)
        CALL get_state (ng, iTLM, 4, ADM(ng)%name, ADrec, Rold(ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO

# ifdef BALANCE_OPERATOR
!
!  Read NLM initial condition in readiness for the balance
!  operator.
!
      IniRec=Rini
      DO ng=1,Ngrids
        CALL get_state (ng, iNLM, 2, INI(ng)%name, IniRec, IniRec)
        IF (exit_flag.ne.NoError) RETURN
        nrhs(ng)=Rini
      END DO
# endif
!
!  Load ADM control read above from TLM into ADM state arrays.
!  Then, multiply control vector by the error covariance standard
!  deviations. Next, convolve resulting adjoint solution with the
!  squared-root adjoint diffusion operator to impose apriori error
!  hypothesis. Notice that the spatial convolution are only done
!  for half of the diffusion steps (squared-root filter). Clear
!  tangent linear state arrays when done.
!
      Lweak=.FALSE.
      add=.FALSE.
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
          CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
          CALL ad_variability (ng, tile, Rold(ng), Lweak)
          CALL ad_convolution (ng, tile, Rold(ng), Lweak, 2)
          CALL initialize_ocean (ng, tile, iTLM)
          CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
          CALL initialize_boundary (ng, tile, iTLM)
# endif
        END DO
!$OMP END PARALLEL
      END DO
!
!  To insure symmetry, convolve resulting filtered adjoint solution
!  from above with the squared-root (half of steps) tangent linear
!  diffusion operator. Then, multiply result with its corresponding
!  error covariance standard deviations. Since the convolved solution
!  is in the adjoint state arrays, first copy to tangent linear state
!  arrays including the ghosts points.
!
      Lweak=.FALSE.
      add=.FALSE.
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL load_ADtoTL (ng, tile, Rold(ng), Rold(ng), add)
          CALL tl_convolution (ng, tile, Rold(ng), Lweak, 2)
          CALL tl_variability (ng, tile, Rold(ng), Lweak)
# ifdef BALANCE_OPERATOR
          CALL tl_balance (ng, tile, Rini, Rold(ng))
# endif
        END DO
!$OMP END PARALLEL

# ifdef POSTERIOR_ERROR_I
!
!  If computing the analysis error covariance matrix, copy TLM back
!  into ADM so that it can be written to Hessian NetCDF file.
!
        IF (Lposterior) THEN
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
          END DO
!$OMP END PARALLEL
        END IF
# endif
      END DO

# if defined ARRAY_MODES        || defined W4DVAR || \
     defined W4DVAR_SENSITIVITY
!
!  Copy back to adjoint state arrays when done with the convolution.
!  Compute representer model initial conditions by adding convolved
!  adjoint solution to the reference nonlinear state (INI(ng)%name,
!  record Rbck).
!
      IF ((model.eq.iRPM).and.                                          &
     &    (INDEX(driver,'w4dvar').ne.0)) THEN
        BckRec=Rbck
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 9, INI(ng)%name, BckRec, Rnew(ng))
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
        add=.FALSE.
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            CALL rp_ini_adjust (ng, tile, Rnew(ng), Rold(ng))
          END DO
!$OMP END PARALLEL
        END DO
      END IF

# elif defined W4DPSAS || defined W4DPSAS_SENSITIVITY
#  ifndef RPCG
!
!  Copy back to adjoint state arrays when done with the convolution.
!  Compute nonlinear model initial conditions by adding convolved
!  adjoint solution to the reference nonlinear state (INI(ng)%name,
!  record Rbck).
!
      IF ((model.eq.iNLM).and.                                          &
     &    (INDEX(driver,'w4dpsas').ne.0)) THEN
        BckRec=Rbck
        DO ng=1,Ngrids
          CALL get_state (ng, iNLM, 9, INI(ng)%name, BckRec, Rnew(ng))
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
        add=.FALSE.
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            CALL ini_adjust (ng, tile, Rold(ng), Rnew(ng))
          END DO
!$OMP END PARALLEL
        END DO
      END IF
#  endif
# endif
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing adjustments for next inner loop into
!  ITL(ng)%name (record Rec1). The tangent model initial
!  conditions are set to the convolved adjoint solution.
!
      IF (model.eq.iTLM) THEN
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, Rold(ng), Rec1)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# if defined ARRAY_MODES        || defined W4DVAR || \
     defined W4DVAR_SENSITIVITY
      ELSE IF (model.eq.iRPM) THEN
        DO ng=1,Ngrids
          CALL rp_wrt_ini (ng, Rold(ng), Rec2)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# elif (defined W4DPSAS || defined W4DPSAS_SENSITIVITY) && \
       !defined RPCG
      ELSE IF (model.eq.iNLM) THEN
        DO ng=1,Ngrids
          CALL wrt_ini (ng, Rnew(ng))
          IF (exit_flag.ne.NoError) RETURN
#  if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
      defined ADJUST_BOUNDARY
          CALL wrt_frc_AD (ng, Rold(ng), INI(ng)%Rindex)
          IF (exit_flag.ne.NoError) RETURN
#  endif
        END DO
# endif
      END IF
# ifdef RPCG
!
!  Write out tangent linear model initial conditions and tangent
!  linear surface forcing and obc adjustments for next outer
!  loop into ITLname (record Rec2). The tangent model initial
!  conditions are set to the convolved adjoint solution.
!
      IF (model.eq.iNLM) THEN
        DO ng=1,Ngrids
          CALL tl_wrt_ini (ng, Rold(ng), Rec2)
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF
# endif

# ifdef POSTERIOR_ERROR_I
!
!  If weak constraint, convolve records 2:Nrec in ADM(ng)%name and
!  Write convolved adjoint solution into Hessian NetCDF file for use
!  later.
!
      IF (Lposterior.and.(inner.ne.0)) THEN
        DO ng=1,Ngrids
          CALL wrt_hessian (ng, Rold(ng), Rold(ng))
          IF (exit_flag.ne.NoERRor) RETURN
        END DO
      END IF
# endif
# ifdef RPCG
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
            CALL get_state (ng, iTLM, 4, ITL(ng)%name, jrec, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL aug_oper (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
!
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
!
!  Need to read adjoint netcdf file in reverse.
!
            ADrec=(Nrec(ng)-1)-(irec-1)
            CALL get_state (ng, iTLM, 4, ADM(ng)%name, ADrec, LTLM1)
            IF (exit_flag.ne.NoError) RETURN
!
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL sum_grad (ng, tile, LTLM1, LTLM2)
            END DO
!$OMP END PARALLEL
!
!  Write the current sum of adjoint solutions into record jrec of the
!  ITL file.
!
            CALL tl_wrt_ini (ng, LTLM2, jrec)
            IF (exit_flag.ne.NoError) RETURN
!
          END DO
        END DO
      END IF
# endif
# ifdef TIME_CONV
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
          DO irec=1,Nrec(ng)-1
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
            ADrec=irec
            CALL get_state (ng, iTLM, 4, ADM(ng)%name, ADrec, Rold(ng))
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
            Lweak=.TRUE.
            add=.FALSE.
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
              CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
              CALL ad_variability (ng, tile, Rold(ng), Lweak)
              CALL ad_convolution (ng, tile, Rold(ng), Lweak, 2)
              CALL initialize_ocean (ng, tile, iTLM)
              CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
              CALL initialize_boundary (ng, tile, iTLM)
# endif
            END DO
!$OMP END PARALLEL
!
!  Overwrite ADM(ng)%name history NetCDF file with convolved adjoint
!  solution.
!
            kstp(ng)=Rold(ng)
# ifdef SOLVE3D
            nstp(ng)=Rold(ng)
# endif
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
          LwrtTime(ng)=.TRUE.
        END IF
      END DO
!
!  Apply the time correlation matrix to the newly convolved adjoint
!  fields.
!
      DO ng=1,Ngrids
        TLF(ng)%Rindex=0
!$OMP PARALLEL
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL time_corr (ng, tile, iADM, ADM(ng)%name)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!$OMP END PARALLEL
      END DO
!
      DO ng=1,Ngrids
        IF (Nrec(ng).gt.3) THEN
          ADM(ng)%Rindex=0
          LwrtState2d(ng)=.TRUE.
          DO irec=Nrec(ng)-1,1,-1
!
!  Read the latest TLF solution. We need to read the TLF file in reverse order
!  since the final solution is written into the adjoint netcdf file which
!  will be rewritten in ascending order of time in the TLF file by wrt_impulse.
!  NOTE: model=14 here so as to only read the state variables.
!
            ADrec=irec
            CALL get_state (ng, 14, 4, TLF(ng)%name, ADrec, Rold(ng))
            IF (exit_flag.ne.NoError) RETURN
!
            Lweak=.TRUE.
            add=.FALSE.
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL tl_convolution (ng, tile, Rold(ng), Lweak, 2)
              CALL tl_variability (ng, tile, Rold(ng), Lweak)
# ifdef BALANCE_OPERATOR
              CALL tl_balance (ng, tile, Rini, Rold(ng))
# endif
              CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            END DO
!$OMP END PARALLEL
!
!  Overwrite ADM(ng)%name history NetCDF file with convolved adjoint
!  solution.
!
            kstp(ng)=Rold(ng)
# ifdef SOLVE3D
            nstp(ng)=Rold(ng)
# endif
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
          LwrtTime(ng)=.TRUE.
        END IF
      END DO

# else
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
          DO irec=1,Nrec(ng)-1
!
!  Read adjoint solution. Since routine "get_state" loads data into the
!  ghost points, the adjoint solution is read in the tangent linear
!  state arrays by using iTLM instead of iADM in the calling arguments.
!
            ADrec=irec
            CALL get_state (ng, iTLM, 4, ADM(ng)%name, ADrec, Rold(ng))
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
            Lweak=.TRUE.
            add=.FALSE.
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
# ifdef BALANCE_OPERATOR
              CALL ad_balance (ng, tile, Rini, Rold(ng))
# endif
              CALL ad_variability (ng, tile, Rold(ng), Lweak)
              CALL ad_convolution (ng, tile, Rold(ng), Lweak, 2)
              CALL initialize_ocean (ng, tile, iTLM)
              CALL initialize_forces (ng, tile, iTLM)
# ifdef ADJUST_BOUNDARY
              CALL initialize_boundary (ng, tile, iTLM)
# endif
            END DO
!$OMP END PARALLEL
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
!$OMP PARALLEL
            DO tile=first_tile(ng),last_tile(ng),+1
              CALL load_ADtoTL (ng, tile, Rold(ng), Rold(ng), add)
              CALL tl_convolution (ng, tile, Rold(ng), Lweak, 2)
              CALL tl_variability (ng, tile, Rold(ng), Lweak)
# ifdef BALANCE_OPERATOR
              CALL tl_balance (ng, tile, Rini, Rold(ng))
# endif
              CALL load_TLtoAD (ng, tile, Rold(ng), Rold(ng), add)
            END DO
!$OMP END PARALLEL
!
!  Overwrite ADM(ng)%name history NetCDF file with final solution.
!
            kstp(ng)=Rold(ng)
# ifdef SOLVE3D
            nstp(ng)=Rold(ng)
# endif
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END DO
          LwrtState2d(ng)=.FALSE.
          LwrtTime(ng)=.TRUE.
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE error_covariance
#endif
      END MODULE convolve_mod
