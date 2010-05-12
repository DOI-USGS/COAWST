      MODULE ocean_control_mod
!
!svn $Id: pert_ocean.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
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
      integer, dimension(Ngrids) :: Tstr
      integer, dimension(Ngrids) :: Tend
!
!  Local variable declarations.
!
      integer :: my_iic, ng, subs, tile, thread
#ifdef SANITY_CHECK
      logical :: BOUNDED_AD, BOUNDED_TL, SAME_VAR
# ifdef DISTRIBUTE
      integer :: Istr, Iend, Jstr, Jend
# endif
      integer :: IperAD, JperAD, KperAD, ivarAD
      integer :: IperTL, JperTL, KperTL, ivarTL
      integer :: i

      real(r8) :: IniVal = 0.0_r8

      real(r8), dimension(4) :: val
#endif
!
!=======================================================================
!  Run model for all nested grids, if any.
!=======================================================================
!
      NEST_LOOP : DO ng=1,Ngrids

#ifdef INNER_PRODUCT
!
!  Set end of pertubation loos as the total number of state variables
!  interior points.
!
# ifdef SOLVE3D
        ERend=(Lm(ng)-1)*Mm(ng)+                                        &
     &        Lm(ng)*(Mm(ng)-1)+                                        &
     &        Lm(ng)*Mm(ng)+                                            &
     &        (Lm(ng)-1)*Mm(ng)*N(ng)+                                  &
     &        Lm(ng)*(Mm(ng)-1)*N(ng)+                                  &
     &        Lm(ng)*Mm(ng)*N(ng)*NAT
# else
        ERend=(Lm(ng)-1)*Mm(ng)+                                        &
     &        Lm(ng)*(Mm(ng)-1)+                                        &
     &        Lm(ng)*Mm(ng)
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
        SAME_VAR=(ivarTL.eq.ivarAD).and.                                &
     &           (IperTL.eq.IperAD).and.                                &
     &           (JperTL.eq.JperAD).and.                                &
     &           (KperTL.eq.KperAD)
# else
        SAME_VAR=(ivarTL.eq.ivarAD).and.                                &
     &           (IperTL.eq.IperAD).and.                                &
     &           (JperTL.eq.JperAD)
# endif
#endif
!
!  Set relevant IO switches.
!
        IF (nTLM(ng).gt.0) LdefTLM(ng)=.TRUE.
        IF (nADJ(ng).gt.0) LdefADJ(ng)=.TRUE.
        Lstiffness=.FALSE.

#if defined BULK_FLUXES || defined NL_BULK_FLUXES

!  Set file name containing the nonlinear model bulk fluxes to be read
!  and processed by other algorithms.
!
        BLKname(ng)=FWDname(ng)
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
          CALL tl_initial (ng)
          IF (exit_flag.ne.NoError) RETURN

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

#ifdef SANITY_CHECK
!
!  Get check value for tangent linear state variable.
!
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
            val(i)=IniVal
          END DO
          IF (BOUNDED_TL) THEN
# ifdef SOLVE3D
            IF (ivarTL.eq.isUbar) THEN
              val(1)=OCEAN(ng)%tl_ubar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(1)=OCEAN(ng)%tl_vbar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(1)=OCEAN(ng)%tl_zeta(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isUvel) THEN
              val(1)=OCEAN(ng)%tl_u(IperTL,JperTL,KperTL,nstp(ng))
            ELSE IF (ivarTL.eq.isVvel) THEN
              val(1)=OCEAN(ng)%tl_v(IperTL,JperTL,KperTL,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarTL.eq.isTvar(i)) THEN
                  val(1)=OCEAN(ng)%tl_t(IperTL,JperTL,KperTL,nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarTL.eq.isUbar) THEN
              val(1)=OCEAN(ng)%tl_ubar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(1)=OCEAN(ng)%tl_vbar(IperTL,JperTL,knew(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(1)=OCEAN(ng)%tl_zeta(IperTL,JperTL,knew(ng))
            END IF
# endif
          END IF

          IF (BOUNDED_AD) THEN
# ifdef SOLVE3D
            IF (ivarAD.eq.isUbar) THEN
              val(3)=OCEAN(ng)%tl_ubar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(3)=OCEAN(ng)%tl_vbar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(3)=OCEAN(ng)%tl_zeta(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isUvel) THEN
              val(3)=OCEAN(ng)%tl_u(IperAD,JperAD,KperAD,nstp(ng))
            ELSE IF (ivarAD.eq.isVvel) THEN
              val(3)=OCEAN(ng)%tl_v(IperAD,JperAD,KperAD,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarAD.eq.isTvar(i)) THEN
                  val(3)=OCEAN(ng)%tl_t(IperAD,JperAD,KperAD,nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarAD.eq.isUbar) THEN
              val(3)=OCEAN(ng)%tl_ubar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(3)=OCEAN(ng)%tl_vbar(IperAD,JperAD,knew(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(3)=OCEAN(ng)%tl_zeta(IperAD,JperAD,knew(ng))
            END IF
# endif
          END IF

#endif
!
!  Clear all model state arrays arrays.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
#if defined _OPENMP || defined DISTRIBUTE
            subs=NtileX(ng)*NtileE(ng)/numthreads
#else
            subs=1
#endif
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, 0)
            END DO
          END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Time step adjoint model backwards.
!-----------------------------------------------------------------------
!
          TLmodel=.FALSE.
          ADmodel=.TRUE.
          CALL ad_initial (ng)
          IF (exit_flag.ne.NoError) RETURN

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

#ifdef SANITY_CHECK
!
!  Get check value for adjoint state variable.
!
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
              val(2)=OCEAN(ng)%ad_ubar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(2)=OCEAN(ng)%ad_vbar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(2)=OCEAN(ng)%ad_zeta(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isUvel) THEN
              val(2)=OCEAN(ng)%ad_u(IperAD,JperAD,KperAD,nstp(ng))
            ELSE IF (ivarAD.eq.isVvel) THEN
              val(2)=OCEAN(ng)%ad_v(IperAD,JperAD,KperAD,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarAD.eq.isTvar(i)) THEN
                  val(2)=OCEAN(ng)%ad_t(IperAD,JperAD,KperAD,nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarAD.eq.isUbar) THEN
              val(2)=OCEAN(ng)%ad_ubar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isVbar) THEN
              val(2)=OCEAN(ng)%ad_vbar(IperAD,JperAD,kstp(ng))
            ELSE IF (ivarAD.eq.isFsur) THEN
              val(2)=OCEAN(ng)%ad_zeta(IperAD,JperAD,kstp(ng))
            END IF
# endif
          END IF

          IF (BOUNDED_TL) THEN
# ifdef SOLVE3D
            IF (ivarTL.eq.isUbar) THEN
              val(4)=OCEAN(ng)%ad_ubar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(4)=OCEAN(ng)%ad_vbar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(4)=OCEAN(ng)%ad_zeta(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isUvel) THEN
              val(4)=OCEAN(ng)%ad_u(IperTL,JperTL,KperTL,nstp(ng))
            ELSE IF (ivarTL.eq.isVvel) THEN
              val(4)=OCEAN(ng)%ad_v(IperTL,JperTL,KperTL,nstp(ng))
            ELSE
              DO i=1,NT(ng)
                IF (ivarTL.eq.isTvar(i)) THEN
                  val(4)=OCEAN(ng)%ad_t(IperTL,JperTL,KperTL,nstp(ng),i)
                END IF
              END DO
            END IF
# else
            IF (ivarTL.eq.isUbar) THEN
              val(4)=OCEAN(ng)%ad_ubar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isVbar) THEN
              val(4)=OCEAN(ng)%ad_vbar(IperTL,JperTL,kstp(ng))
            ELSE IF (ivarTL.eq.isFsur) THEN
              val(4)=OCEAN(ng)%ad_zeta(IperTL,JperTL,kstp(ng))
            END IF
# endif
          END IF

!
!  Report sanity check values.
!
# ifdef DISTRIBUTE
          CALL mp_collect (ng, iTLM, 4, IniVal, val)
# endif
          IF (Master) THEN
            IF (SAME_VAR) THEN
              WRITE (stdout,20) 'Perturbing',                           &
     &                          TRIM(Vname(1,idSvar(ivarTL)))
              IF (ivarTL.le.3) THEN
                WRITE (stdout,30) 'Tangent:    ', val(1),IperTL,JperTL
                WRITE (stdout,30) 'Adjoint:    ', val(2),IperAD,JperAD
                WRITE (stdout,30) 'Difference: ', val(2)-val(1),        &
     &                                            IperTL,JperTL
              ELSE
                WRITE (stdout,40) 'Tangent:    ', val(1),IperTL,JperTL, &
     &                                                 KperTL
                WRITE (stdout,40) 'Adjoint:    ', val(2),IperAD,JperAD, &
     &                                                   KperAD
                WRITE (stdout,40) 'Difference: ', val(2)-val(1),        &
     &                                            IperTL,JperTL,KperTL
              END IF
            ELSE
              IF (ivarTL.le.3) THEN
                WRITE (stdout,50) 'Tangent, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            IperTL,JperTL
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(1),IperTL,JperTL
              ELSE
                WRITE (stdout,70) 'Tangent, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            IperTL,JperTL,KperTL
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(1),IperTL,JperTL,KperTL
              END IF
              IF (ivarAD.le.3) THEN
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(3),IperAD,JperAD
              ELSE
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(3),IperAD,JperAD,KperAD
              END IF
!
              IF (ivarAD.le.3) THEN
                WRITE (stdout,50) 'Adjoint, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            IperAD,JperAD
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(2),IperAD,JperAD
              ELSE
                WRITE (stdout,70) 'Adjoint, Perturbing: ',              &
     &                            TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            IperAD,JperAD,KperAD
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarAD))),        &
     &                            val(2),IperAD,JperAD,KperAD
              END IF
              IF (ivarTL.le.3) THEN
                WRITE (stdout,60) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(4),IperTL,JperTL
              ELSE
                WRITE (stdout,80) TRIM(Vname(1,idSvar(ivarTL))),        &
     &                            val(4),IperTL,JperTL,KperTL
              END IF
              WRITE (stdout,90) val(3)-val(4)
            END IF
          END IF
#endif
!
!  Clear model state arrays arrays.
!
!$OMP PARALLEL DO PRIVATE(ng,thread,subs,tile) SHARED(numthreads)
          DO thread=0,numthreads-1
#if defined _OPENMP || defined DISTRIBUTE
            subs=NtileX(ng)*NtileE(ng)/numthreads
#else
            subs=1
#endif
            DO tile=subs*thread,subs*(thread+1)-1
              CALL initialize_ocean (ng, TILE, 0)
            END DO
          END DO
!$OMP END PARALLEL DO

        END DO PERT_LOOP
!
!  Close current forward NetCDF file.
!
        SourceFile='pert_ocean.h, ROMS_run'

        CALL netcdf_close (ng, iTLM, ncFWDid(ng), Lupdate = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

      END DO NEST_LOOP
!
 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
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
