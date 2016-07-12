      MODULE ocean_control_mod
!
!svn $Id: pert_ocean.h 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ROMS/TOMS Tangent Linear and Adjoint Models Sanity Test:            !
!                                                                      !
!  This driver is used to test the tangent linear model and adjoint    !
!  models for ROMS.  It is used to check whether or not both models    !
!  are correct. This driver is usually run in "ensemble" mode where    !
!  each member represents a perturbed interior point in the tangent    !
!  linear and adjoint models.  There is a maximum of Lm*Mm ensemble    !
!  members. Each member is ran for few time steps.                     !
!                                                                      !
!  If each interior point is perturbed at one time,  the  resulting    !
!  tangent linear (T) and adjoint (A) M-by-N matrices yield:           !
!                                                                      !
!                T - tranpose(A) = 0    within round off               !
!                                                                      !
!  That is, their inner product give a symmetric matrix. Here, M is    !
!  the number of state points and N is the number of perturbations.    !
!  This option is activated with the INNER_PRODUCT test.               !
!                                                                      !
!  In realistic applications, it is awkward to perturb all interior    !
!  points for each state variable. Alternatively, random check at a    !
!  specified point is inexpensive. This option is activate with the    !
!  SANITY_CHECK switch.  The standard input "User" array is used to    !
!  specify such point:                                                 !
!                                                                      !
!     INT(user(1)) => state tangent variable to perturb                !
!     INT(user(2)) => state adjoint variable to perturb                !
!     INT(user(3)) => I-index of tangent variable to perturb           !
!     INT(user(4)) => I-index of adjoint variable to perturb           !
!     INT(user(5)) => J-index of tangent variable to perturb           !
!     INT(user(6)) => J-index of adjoint variable to perturb           !
!                                                                      !
!  In 3D applications:                                                 !
!                                                                      !
!     INT(user(7)) => J-index of tangent variable to perturb           !
!     INT(user(8)) => J-index of adjoint variable to perturb           !
!                                                                      !
!  The subroutines in this driver control the initialization, time-    !
!  stepping,  and  finalization of  ROMS/TOMS  model following ESMF    !
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

      RETURN
      END SUBROUTINE ROMS_initialize

      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine time-steps ROMS/TOMS tangent linear and adjoint        !
!  models.                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
#endif
!
!  Imported variable declarations
!
      real(r8), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      integer :: ng, tile
#ifdef SANITY_CHECK
      logical :: BOUNDED_AD, BOUNDED_TL, SAME_VAR
# ifdef DISTRIBUTE
      integer :: Istr, Iend, Jstr, Jend
# endif
      integer :: IperAD, JperAD, KperAD, ivarAD
      integer :: IperTL, JperTL, KperTL, ivarTL
      integer :: i

      real(r8) :: IniVal = 0.0_r8

      real(r8), dimension(4,Ngrids) :: val
#endif
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
#ifdef INNER_PRODUCT
!
!  Set end of pertubation loos as the total number of state variables
!  interior points. Is not possible to run the inner product test over
!  nested grids. If nested grids, run each grid separately.
!
      IF (Ngrids.gt.1) THEN
        WRITE (stdout,10) 'Nested grids are not allowed, Ngrids = ',    &
                          Ngrids
        STOP
      END IF
      ng=1
# ifdef SOLVE3D
      ERend=(Lm(ng)-1)*Mm(ng)+                                          &
     &       Lm(ng)*(Mm(ng)-1)+                                         &
     &       Lm(ng)*Mm(ng)+                                             &
     &       (Lm(ng)-1)*Mm(ng)*N(ng)+                                   &
     &       Lm(ng)*(Mm(ng)-1)*N(ng)+                                   &
     &       Lm(ng)*Mm(ng)*N(ng)*NAT
# else
      ERend=(Lm(ng)-1)*Mm(ng)+                                          &
     &       Lm(ng)*(Mm(ng)-1)+                                         &
     &       Lm(ng)*Mm(ng)
# endif
#endif
#ifdef SANITY_CHECK
!
!  Set tangent and adjoint variable and random point to perturb.
!
      ivarTL=INT(user(1))
      ivarAD=INT(user(2))
      IperTL=INT(user(3))
      IperAD=INT(user(4))
      JperTL=INT(user(5))
      JperAD=INT(user(6))
# ifdef SOLVE3D
      KperTL=INT(user(7))
      KperAD=INT(user(8))
      SAME_VAR=(ivarTL.eq.ivarAD).and.                                  &
     &         (IperTL.eq.IperAD).and.                                  &
     &         (JperTL.eq.JperAD).and.                                  &
     &         (KperTL.eq.KperAD)
# else
      SAME_VAR=(ivarTL.eq.ivarAD).and.                                  &
     &         (IperTL.eq.IperAD).and.                                  &
     &         (JperTL.eq.JperAD)
# endif
#endif
!
!  Set relevant IO switches.
!
      DO ng=1,Ngrids
        IF (nTLM(ng).gt.0) LdefTLM(ng)=.TRUE.
        IF (nADJ(ng).gt.0) LdefADJ(ng)=.TRUE.
      END DO
      Lstiffness=.FALSE.

#if defined BULK_FLUXES || defined NL_BULK_FLUXES

!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
      DO ng=1,Ngrids
        BLK(ng)%name=FWD(ng)%name
      END DO
#endif
!
!=======================================================================
!  Perturbation loop.
!=======================================================================
!
      PERT_LOOP : DO Nrun=ERstr,ERend
!
!-----------------------------------------------------------------------
!  Time step tangent linear model.
!-----------------------------------------------------------------------
!
        TLmodel=.TRUE.
        ADmodel=.FALSE.
        DO ng=1,Ngrids
!$OMP PARALLEL
          CALL tl_initial (ng)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN

          IF (Master) THEN
            WRITE (stdout,10) 'TL', ng, ntstart(ng), ntend(ng)
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

#ifdef SANITY_CHECK
!
!  Get check value for tangent linear state variable.
!
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          Istr=BOUNDS(ng)%Istr(MyRank)
          Iend=BOUNDS(ng)%Iend(MyRank)
          Jstr=BOUNDS(ng)%Jstr(MyRank)
          Jend=BOUNDS(ng)%Jend(MyRank)
          BOUNDED_TL=((Istr.le.IperTL).and.(IperTL.le.Iend)).and.       &
     &               ((Jstr.le.JperTL).and.(JperTL.le.Jend))
          BOUNDED_AD=((Istr.le.IperAD).and.(IperAD.le.Iend)).and.       &
     &               ((Jstr.le.JperAD).and.(JperAD.le.Jend))
# else
          BOUNDED_TL=.TRUE.
          BOUNDED_AD=.TRUE.
# endif
          DO i=1,4
            val(i,ng)=IniVal
          END DO
          IF (BOUNDED_TL) THEN
# ifdef SOLVE3D
            IF (ivarTL.eq.isUbar) THEN
              val(1,ng)=OCEAN(ng)%tl_ubar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(1,ng)=OCEAN(ng)%tl_vbar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(1,ng)=OCEAN(ng)%tl_zeta(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isUvel) THEN
              val(1,ng)=OCEAN(ng)%tl_u(IperTL,JperTL,KperTL,nstp(ng))
            ELSE IF (ivarTL.eq.isVvel) THEN
              val(1,ng)=OCEAN(ng)%tl_v(IperTL,JperTL,KperTL,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarTL.eq.isTvar(i)) THEN
                  val(1,ng)=OCEAN(ng)%tl_t(IperTL,JperTL,KperTL,        &
     &                                     nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarTL.eq.isUbar) THEN
              val(1,ng)=OCEAN(ng)%tl_ubar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(1,ng)=OCEAN(ng)%tl_vbar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(1,ng)=OCEAN(ng)%tl_zeta(IperTL,JperTL,knew(ng))
            END IF
# endif
          END IF

          IF (BOUNDED_AD) THEN
# ifdef SOLVE3D
            IF (ivarAD.eq.isUbar) THEN
              val(3,ng)=OCEAN(ng)%tl_ubar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(3,ng)=OCEAN(ng)%tl_vbar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(3,ng)=OCEAN(ng)%tl_zeta(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isUvel) THEN
              val(3,ng)=OCEAN(ng)%tl_u(IperAD,JperAD,KperAD,nstp(ng))
            ELSE IF (ivarAD.eq.isVvel) THEN
              val(3,ng)=OCEAN(ng)%tl_v(IperAD,JperAD,KperAD,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarAD.eq.isTvar(i)) THEN
                  val(3,ng)=OCEAN(ng)%tl_t(IperAD,JperAD,KperAD,        &
     &                                     nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarAD.eq.isUbar) THEN
              val(3,ng)=OCEAN(ng)%tl_ubar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(3,ng)=OCEAN(ng)%tl_vbar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(3,ng)=OCEAN(ng)%tl_zeta(IperAD,JperAD,knew(ng))
            END IF
# endif
          END IF
        END DO
#endif
!
!  Clear all model state arrays arrays.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, 0)
          END DO
!$OMP END PARALLEL
        END DO
!
!-----------------------------------------------------------------------
!  Time step adjoint model backwards.
!-----------------------------------------------------------------------
!
        TLmodel=.FALSE.
        ADmodel=.TRUE.
        DO ng=1,Ngrids
!$OMP PARALLEL
          CALL ad_initial (ng)
!$OMP END PARALLEL
          IF (exit_flag.ne.NoError) RETURN

          IF (Master) THEN
            WRITE (stdout,10) 'AD', ng, ntstart(ng), ntend(ng)
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

#ifdef SANITY_CHECK
!
!  Get check value for adjoint state variable.
!
        DO ng=1,Ngrids
# ifdef DISTRIBUTE
          Istr=BOUNDS(ng)%Istr(MyRank)
          Iend=BOUNDS(ng)%Iend(MyRank)
          Jstr=BOUNDS(ng)%Jstr(MyRank)
          Jend=BOUNDS(ng)%Jend(MyRank)
          BOUNDED_AD=((Istr.le.IperAD).and.(IperAD.le.Iend)).and.       &
     &               ((Jstr.le.JperAD).and.(JperAD.le.Jend))
          BOUNDED_TL=((Istr.le.IperTL).and.(IperTL.le.Iend)).and.       &
     &               ((Jstr.le.JperTL).and.(JperTL.le.Jend))
# else
          BOUNDED_AD=.TRUE.
          BOUNDED_TL=.TRUE.
# endif
          IF (BOUNDED_AD) THEN
# ifdef SOLVE3D
            IF (ivarAD.eq.isUbar) THEN
              val(2,ng)=OCEAN(ng)%ad_ubar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(2,ng)=OCEAN(ng)%ad_vbar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(2,ng)=OCEAN(ng)%ad_zeta(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isUvel) THEN
              val(2,ng)=OCEAN(ng)%ad_u(IperAD,JperAD,KperAD,nstp(ng))
            ELSE IF (ivarAD.eq.isVvel) THEN
              val(2,ng)=OCEAN(ng)%ad_v(IperAD,JperAD,KperAD,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarAD.eq.isTvar(i)) THEN
                  val(2,ng)=OCEAN(ng)%ad_t(IperAD,JperAD,KperAD,        &
     &                                     nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarAD.eq.isUbar) THEN
              val(2,ng)=OCEAN(ng)%ad_ubar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(2,ng)=OCEAN(ng)%ad_vbar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(2,ng)=OCEAN(ng)%ad_zeta(IperAD,JperAD,kstp(ng))
            END IF
# endif
          END IF

          IF (BOUNDED_TL) THEN
# ifdef SOLVE3D
            IF (ivarTL.eq.isUbar) THEN
              val(4,ng)=OCEAN(ng)%ad_ubar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(4,ng)=OCEAN(ng)%ad_vbar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(4,ng)=OCEAN(ng)%ad_zeta(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isUvel) THEN
              val(4,ng)=OCEAN(ng)%ad_u(IperTL,JperTL,KperTL,nstp(ng))
            ELSE IF (ivarTL.eq.isVvel) THEN
              val(4,ng)=OCEAN(ng)%ad_v(IperTL,JperTL,KperTL,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarTL.eq.isTvar(i)) THEN
                  val(4,ng)=OCEAN(ng)%ad_t(IperTL,JperTL,KperTL,        &
     &                                     nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarTL.eq.isUbar) THEN
              val(4,ng)=OCEAN(ng)%ad_ubar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(4,ng)=OCEAN(ng)%ad_vbar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(4,ng)=OCEAN(ng)%ad_zeta(IperTL,JperTL,kstp(ng))
            END IF
# endif
          END IF

!
!  Report sanity check values.
!
# ifdef DISTRIBUTE
          CALL mp_collect (ng, iTLM, 4, IniVal, val(1,ng))
# endif
          IF (Master) THEN
            IF (SAME_VAR) THEN
              WRITE (stdout,20) 'Perturbing',                           &
     &                          TRIM(Vname(1,idSvar(ivarTL)))
              IF (ivarTL.le.3) THEN
                WRITE (stdout,30) 'Tangent:    ', val(1,ng),            &
     &                                            IperTL,JperTL
                WRITE (stdout,30) 'Adjoint:    ', val(2,ng),            &
     &                                            IperAD, JperAD
                WRITE (stdout,30) 'Difference: ', val(2,ng)-val(1,ng),  &
     &                                            IperTL,JperTL
              ELSE
                WRITE (stdout,40) 'Tangent:    ', val(1,ng),            &
     &                                            IperTL,JperTL,KperTL
                WRITE (stdout,40) 'Adjoint:    ', val(2,ng),            &
     &                                            IperAD,JperAD,KperAD
                WRITE (stdout,40) 'Difference: ', val(2,ng)-val(1,ng),  &
     &                                            IperTL,JperTL,KperTL
              END IF
            ELSE
              IF (ivarTL.le.3) THEN
                WRITE (stdout,50) 'Tangent, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            IperTL,JperTL
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(1,ng),IperTL,JperTL
              ELSE
                WRITE (stdout,70) 'Tangent, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            IperTL,JperTL,KperTL
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(1,ng),IperTL,JperTL,KperTL
              END IF
              IF (ivarAD.le.3) THEN
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(3,ng),IperAD,JperAD
              ELSE
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(3,ng),IperAD,JperAD,KperAD
              END IF
!
              IF (ivarAD.le.3) THEN
                WRITE (stdout,50) 'Adjoint, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            IperAD,JperAD
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(2,ng),IperAD,JperAD
              ELSE
                WRITE (stdout,70) 'Adjoint, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            IperAD,JperAD,KperAD
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(2,ng),IperAD,JperAD,KperAD
              END IF
              IF (ivarTL.le.3) THEN
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(4,ng),IperTL,JperTL
              ELSE
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(4,ng),IperTL,JperTL,KperTL
              END IF
              WRITE (stdout,90) val(3,ng)-val(4,ng)
            END IF
          END IF
        END DO
#endif
!
!  Clear model state arrays arrays.
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL initialize_ocean (ng, tile, 0)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Close current forward NetCDF file.
!
        SourceFile='pert_ocean.h, ROMS_run'

        DO ng=1,Ngrids
          CALL netcdf_close (ng, iTLM, FWD(ng)%ncid, Lupdate = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END DO

      END DO PERT_LOOP
!
 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        ' (Grid: ',i2.2,' TimeSteps: ',i8.8,' - ',i8.8,')',/)
#ifdef SANITY_CHECK
 20   FORMAT (/,' Sanity Check - ',a,' variable: ',a,t60)
 30   FORMAT (' Sanity Check - ', a, 1p,e19.12,                         &
     &        3x, 'at (i,j)   ',2i4)
 40   FORMAT (' Sanity Check - ', a, 1p,e19.12,                         &
     &        3x, 'at (i,j,k) ',3i4)
 50   FORMAT (/,' Sanity Check - ',a, a, t52, 'at (i,j)   ', 2i4)
 60   FORMAT (' Sanity Check - ', a, ' =', t30, 1p,e19.12,              &
     &        t52, 'at (i,j)   ',2i4)
 70   FORMAT (/,' Sanity Check - ',a, a, t52, 'at (i,j,k) ', 3i4)
 80   FORMAT (' Sanity Check - ', a, ' =', t30, 1p,e19.12,              &
     &        t52, 'at (i,j,k) ',3i4)
 90   FORMAT (/,' Sanity Check - Difference = ', 1p,e19.12)
#endif

      RETURN
      END SUBROUTINE ROMS_run

      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS/TOMS tangent linear and adjoint        !
!  models execution.                                                   !
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
