      SUBROUTINE propagator (ng, Nstr, Nend, state, ad_state)
!
!svn $Id: propagator_so_semi.h 652 2008-07-24 23:20:53Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Stochastic Optimals, Seminorm Estimation:                           !
!                                                                      !
!  This routine is used during the computation of the eigenvectors of  !
!  the stochastic optimals operator  with respect the seminorm of the  !
!  chosen functional.                                                  !
!                                                                      !
!   Reference:                                                         !
!                                                                      !
!     Moore, A.M. et al., 2004: A comprehensive ocean prediction and   !
!       analysis system based on the tangent linear and adjoint of a   !
!       regional ocean model, Ocean Modelling, 7, 227-258.             !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
#ifdef SO_SEMI_WHITE
      USE packing_mod, ONLY : so_semi_white
#else
      USE packing_mod, ONLY : so_semi_red
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Nstr, Nend

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Nstr:)
      real(r8), intent(out) :: ad_state(Nstr:)
#else
      real(r8), intent(in) :: state(Nstr:Nend)
      real(r8), intent(out) :: ad_state(Nstr:Nend)
#endif
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      logical :: FirstPass = .TRUE.
#endif
      integer :: my_iic, subs, tile, thread
!
!=======================================================================
!  Backward integration of adjoint model forced with the seminorm of
!  the chosen functional. The adjoint model is run only only once in
!  the first iteration.
!=======================================================================
!
      Nrun=Nrun+1
!
      FIRST_PASS : IF (Nrun.eq.1) THEN
!
!  Initialize the adjoint model always from rest.
!
        CALL ad_initial (ng)
        IF (exit_flag.ne.NoError) RETURN
!
!  Activate adjoint output.
!
        LdefADJ(ng)=.TRUE.
        LwrtADJ(ng)=.TRUE.
        LcycleADJ(ng)=.FALSE.
!
!  Time-step adjoint model forced with chosen functional at initial
!  time only.
!
        DstrS(ng)=time(ng)*sec2day
        DendS(ng)=DstrS(ng)

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

      END IF FIRST_PASS
!
!-----------------------------------------------------------------------
!  Compute new packed adjoint state vector containing surface forcing
!  variables.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE(thread,subs,tile)                             &
!$OMP&            SHARED(ng,numthreads,Nstr,Nend,state,ad_state)
      DO thread=0,numthreads-1
        subs=NtileX(ng)*NtileE(ng)/numthreads
        DO tile=subs*thread,subs*(thread+1)-1,+1
# ifdef SO_SEMI_WHITE
          CALL so_semi_white (ng, TILE, Nstr, Nend, state, ad_state)
# else
          CALL so_semi_red (ng, TILE, Nstr, Nend, state, ad_state)
# endif
        END DO
      END DO
!$OMP END PARALLEL DO
!
!-----------------------------------------------------------------------
!  Report iteration and trace or stochastic optimals matrix.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        WRITE (stdout,20) ' PROPAGATOR - Iteration Run: ', Nrun,        &
     &                    ',  number converged RITZ values: ', Nconv,   &
     &                    'TRnorm = ', TRnorm(ng)
      END IF

 10   FORMAT (/,1x,a,1x,'ROMS/TOMS: started time-stepping:',            &
     &        '( TimeSteps: ',i8.8,' - ',i8.8,')',/)
 20   FORMAT (/,a,i3,a,i3,/,35x,a,1p,e15.8)

      RETURN
      END SUBROUTINE propagator
