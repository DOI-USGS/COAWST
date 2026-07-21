/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore  **
**    Licensed under a MIT/X style license                            **
**    See License_ROMS.md                                             **
**                                                                    **
************************************************************************
**                                                                    **
** Implicit Multiscale Background-error Covariance Solver:            **
**                                                                    **
** Routines to invert the implicit diffusion equation using Chebyshev **
** Iterations (CI). It requires estimates of the maximum and minimum  **
** eigenvalues of the K-Laplacian operator.                           **
**                                                                    **
** Reference:                                                         **
**                                                                    **
** Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation   **
**   operators based on an implicitly formulated diffusion equation   **
**   solved with the Chebyshev iteration, Q. J. R. Meteorol. Soc.,    **
**   142, 455-471, doi:10.1002/qj.2664.                               **
**                                                                    **
************************************************************************
**
**
**  <><> CLASS MULTISCALE:  Chebyshev Iterations Solver <><><><><><><><>
**
*/

!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for 2D variables in the control
!  vector. It operates in conjunction with its adjoint to enforce
!  symmetry-preserving correlation functions. The algorithm requires
!  the extrema eigenvalues of matrix A, which are obtained from the
!  Conjugate Gradient (CG) solver. Additionally, initial values for
!  the estimate x(0) and its residual r(0) are required.
!
      SUBROUTINE multiscale_CI_2d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
      integer                           :: itrc
#ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#endif
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCI)               :: ci_alpha
      real (r8), dimension(0:NiterCI+1)             :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_2d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ci2d_r=0.0_r8
      self%ci2d_p=0.0_r8
      self%ci2d_q=0.0_r8
      self%ci2d_x=0.0_r8
      tl_scale=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extrema
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%zeta_eigen(:,ms,1)
          eigMax => self%zeta_eigen(:,ms,2)
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%ubar_eigen(:,ms,1)
          eigMax => self%ubar_eigen(:,ms,2)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%vbar_eigen(:,ms,1)
          eigMax => self%vbar_eigen(:,ms,2)
#ifdef ADJUST_WSTRESS
        CASE ('sustr')                         ! surface U-stress
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%sustr_eigen(:,ms,1)
          eigMax => self%sustr_eigen(:,ms,2)
        CASE ('svstr')                         ! surface V-stress
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%svstr_eigen(:,ms,1)
          eigMax => self%svstr_eigen(:,ms,2)
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
        CASE ('shflux', 'ssflux')              ! surface trace flux
          itrc = tracer_index(TRIM(StateVarName(ifield)))
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%stflux_eigen(:,ms,1,itrc)
          eigMax => self%stflux_eigen(:,ms,2,itrc)
#endif
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j)=tl_scale(i,j)*tl_A(i,j)
            self%ci2d_x(i,j)=0.0_r8                     ! step 9
            self%ci2d_r(i,j)=-tl_A(i,j)                 ! step 10
            self%ci2d_p(i,j)=tl_A(i,j)                  ! step 11
          END DO
        END DO

#ifdef MULTI_SCALE_DEBUG
!
        dotn = dot_prod_2d (ng, tile, model, ctype,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      self%ci2d_r, self%ci2d_r)
#endif
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)/tl_scale(i,j)
            END DO
          END DO
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r2dvar)
              CALL self%tl_Klap_r2d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               tl_A)
            CASE (u2dvar)
              CALL self%tl_Klap_u2d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               tl_A)
            CASE (v2dvar)
              CALL self%tl_Klap_v2d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               tl_A)
          END SELECT
!                                                        step 13
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)*tl_scale(i,j)
              self%ci2d_q(i,j)=tl_A(i,j)
            END DO
          END DO

#ifdef MULTI_SCALE_DEBUG
!
          dotr = dot_prod_2d (ng, tile, model, ctype,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        self%ci2d_r, self%ci2d_r)
          deps=SQRT(dotr/dotn)
!
          status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,       &
     &                            iterCI, NiterCI, 'CI_2d', deps)
#endif
!                                                        steps 14, 15
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci2d_x(i,j)=self%ci2d_x(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_p(i,j)
              self%ci2d_r(i,j)=self%ci2d_r(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_q(i,j)
            END DO
          END DO
!                                                        step 16
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci2d_p(i,j)=-self%ci2d_r(i,j)+                       &
     &                         ci_beta(iterCI+1)*self%ci2d_p(i,j)
              tl_A(i,j)=self%ci2d_p(i,j)
            END DO
          END DO

        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%ci2d_p(i,j)=0.0_r8
            self%ci2d_r(i,j)=0.0_r8
            self%ci2d_q(i,j)=0.0_r8
          END DO
        END DO
!
!  Update control vector 2D variable with the implicit diffusion
!  solution iterate.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j)=self%ci2d_x(i,j)/tl_scale(i,j)
          END DO
        END DO
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_2d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, utilizing Chebyshev Iterations
!  (CI) to model the multiscale propagation of background-error
!  covariance for 2D variables in the control vector. It is needed to
!  enforce symmetry-preserving correlation functions. The algorithm
!  requires the extrema eigenvalues of matrix A, which are obtained
!  from the Conjugate Gradient (CG) solver.
!
      SUBROUTINE multiscale_CI_2d_ad (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
      integer                           :: itrc
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_scale
      real (r8), dimension(0:NiterCI)               :: ci_alpha
      real (r8), dimension(0:NiterCI+1)             :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_2d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      ad_scale=0.0_r8
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ci2d_r=0.0_r8
      self%ci2d_p=0.0_r8
      self%ci2d_q=0.0_r8
      self%ci2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extrema
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%zeta_eigen(:,ms,1)
          eigMax => self%zeta_eigen(:,ms,2)
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%ubar_eigen(:,ms,1)
          eigMax => self%ubar_eigen(:,ms,2)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%vbar_eigen(:,ms,1)
          eigMax => self%vbar_eigen(:,ms,2)
#ifdef ADJUST_WSTRESS
        CASE ('sustr')                         ! surface U-stress
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%sustr_eigen(:,ms,1)
          eigMax => self%sustr_eigen(:,ms,2)
        CASE ('svstr')                         ! surface V-stress
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%svstr_eigen(:,ms,1)
          eigMax => self%svstr_eigen(:,ms,2)
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
        CASE ('shflux', 'ssflux')              ! surface tracer flux
          itrc = tracer_index(TRIM(StateVarName(ifield)))
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%stflux_eigen(:,ms,1,itrc)
          eigMax => self%stflux_eigen(:,ms,2,itrc)
#endif
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
!        DO j=Jmin,Jmax
!          DO i=Imin,Imax
!            self%ci2d_p(i,j)=0.0_r8
!            self%ci2d_r(i,j)=0.0_r8
!            self%ci2d_q(i,j)=0.0_r8
!            self%ci2d_x(i,j)=0.0_r8
!          END DO
!        END DO
!
!  Adjoint of invert the implicit diffusion equation using Chebyshev
!  Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
!  Adjoint of update control vector 2D variable with the implicit
!  diffusion solution iterate.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          tl_A(i,j)=self%ci2d_x(i,j)/tl_scale(i,j)
!^
            self%ci2d_x(i,j)=ad_A(i,j)/ad_scale(i,j)
            ad_A(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of reset iterative conjugate gradient arrays.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          self%ci2d_q(i,j)=0.0_r8
!^
            self%ci2d_q(i,j)=0.0_r8
!^          self%ci2d_r(i,j)=0.0_r8
!^
            self%ci2d_r(i,j)=0.0_r8
!^          self%ci2d_p(i,j)=0.0_r8
!^
            self%ci2d_p(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          DO j=Jmin,Jmax                                ! step 16
            DO i=Imin,Imax
!^            tl_A(i,j)=self%ci2d_p(i,j)
!^
              self%ci2d_p(i,j)=self%ci2d_p(i,j)+ad_A(i,j)
              ad_A(i,j)=0.0_r8
!^            self%ci2d_p(i,j)=-self%ci2d_r(i,j)+                       &
!^   &                         ci_beta(iterCI+1)*self%ci2d_p(i,j)
!^
              self%ci2d_r(i,j)=self%ci2d_r(i,j)-self%ci2d_p(i,j)
              self%ci2d_p(i,j)=ci_beta(iterCI+1)*self%ci2d_p(i,j)
            END DO
          END DO
!                                                        steps 15, 14
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci2d_r(i,j)=self%ci2d_r(i,j)+                        &
!^   &                         ci_alpha(iterCI)*self%ci2d_q(i,j)
!^
              self%ci2d_q(i,j)=self%ci2d_q(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_r(i,j)
!^            self%ci2d_x(i,j)=self%ci2d_x(i,j)+                        &
!^   &                         ci_alpha(iterCI)*self%ci2d_p(i,j)
!^
              self%ci2d_p(i,j)=self%ci2d_p(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_x(i,j)
            END DO
          END DO
!                                                        step 13
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci2d_q(i,j)=tl_A(i,j)
!^
              ad_A(i,j)=ad_A(i,j)+self%ci2d_q(i,j)
              self%ci2d_q(i,j)=0.0_r8
!^            tl_A(i,j)=tl_A(i,j)*tl_scale(i,j)
!^
              ad_A(i,j)=ad_A(i,j)*ad_scale(i,j)
            END DO
          END DO
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r2dvar)
              CALL self%ad_Klap_r2d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               ad_A)
            CASE (u2dvar)
              CALL self%ad_Klap_u2d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               ad_A)
            CASE (v2dvar)
              CALL self%ad_Klap_v2d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               ad_A)
          END SELECT
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            tl_A(i,j)=tl_A(i,j)/tl_scale(i,j)
!^
              ad_A(i,j)=ad_A(i,j)/ad_scale(i,j)
            END DO
          END DO

        END DO KLAP_ITER
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          self%cg2d_p(i,j)=tl_A(i,j)                    ! step 11
!^
            ad_A(i,j)=ad_A(i,j)+self%ci2d_p(i,j)
            self%ci2d_p(i,j)=0.0_r8
!^          self%cg2d_r(i,j)=-tl_A(i,j)                   ! step 10
!^
            ad_A(i,j)=ad_A(i,j)-self%ci2d_r(i,j)
            self%ci2d_r(i,j)=0.0_r8
!^          self%cg2d_x(i,j)=0.0_r8                       ! step 9
!^
            self%ci2d_x(i,j)=0.0_r8
!^          tl_A(i,j)=tl_scale(i,j)*tl_A(i,j)
!^
            ad_A(i,j)=ad_scale(i,j)*ad_A(i,j)
          END DO
        END DO

      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_2d_ad

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for 3D variables in the control
!  vector. It operates in conjunction with its adjoint to enforce
!  symmetry-preserving correlation functions. The algorithm requires
!  the extrema eigenvalues of matrix A, which are obtained from the
!  Conjugate Gradient (CG) solver. Additionally, initial values for
!  the estimate x(0) and its residual r(0) are required.
!
      SUBROUTINE multiscale_CI_3d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
      integer                           :: itrc
# ifdef MULTI_SCALE_DEBUG
      integer                           :: status
# endif
!
      real (r8)                         :: dotr, deps
      real (r8), dimension(N(ng))       :: dotn
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCI,N(ng))         :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng))       :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_3d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ci3d_r=0.0_r8
      self%ci3d_p=0.0_r8
      self%ci3d_q=0.0_r8
      self%ci3d_x=0.0_r8
      tl_scale=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extrema
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%u_eigen(:,:,ms,1)
          eigMax => self%u_eigen(:,:,ms,2)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%v_eigen(:,:,ms,1)
          eigMax => self%v_eigen(:,:,ms,2)
        CASE ('temp', 'salt')                  ! tracers
          itrc = tracer_index(TRIM(StateVarName(ifield)))
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%t_eigen(:,:,ms,1,itrc)
          eigMax => self%t_eigen(:,:,ms,2,itrc)
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!                                            why r = -tl_A & p = +tl_A?
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j,k)=tl_scale(i,j)*tl_A(i,j,k)
              self%ci3d_x(i,j,k)=0.0_r8                 ! step 9
              self%ci3d_r(i,j,k)=-tl_A(i,j,k)           ! step 10
              self%ci3d_p(i,j,k)=tl_A(i,j,k)            ! step 11
            END DO
          END DO

# ifdef MULTI_SCALE_DEBUG
!
          dotn(k) = dot_prod_2d (ng, tile, model, ctype,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           self%ci3d_r(:,:,k),                    &
     &                           self%ci3d_r(:,:,k))
# endif
        END DO LEVEL_LOOP1
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)/tl_scale(i,j)
              END DO
            END DO
          END DO LEVEL_LOOP2
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r3dvar)
              CALL self%tl_Klap_r3d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               tl_A)
            CASE (u3dvar)
              CALL self%tl_Klap_u3d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               tl_A)
            CASE (v3dvar)
              CALL self%tl_Klap_v3d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               tl_A)
          END SELECT
!                                                        step 13
          LEVEL_LOOP3 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)*tl_scale(i,j)
                self%ci3d_q(i,j,k)=tl_A(i,j,k)
              END DO
            END DO

# ifdef MULTI_SCALE_DEBUG
!
            dotr = dot_prod_2d (ng, tile, model, ctype,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          self%ci3d_r(:,:,k),                     &
     &                          self%ci3d_r(:,:,k))
            deps=sqrt(dotr/dotn(k))
!
            status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,     &
     &                              iterCI, NiterCI, 'CI_3d', deps,     &
     &                              level = k)
# endif
!                                                        steps 14, 15
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%ci3d_x(i,j,k)=self%ci3d_x(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_p(i,j,k)
                self%ci3d_r(i,j,k)=self%ci3d_r(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_q(i,j,k)
              END DO
            END DO
!                                                        step 16
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%ci3d_p(i,j,k)=-self%ci3d_r(i,j,k)+                 &
     &                             ci_beta(iterCI+1,k)*                 &
     &                             self%ci3d_p(i,j,k)
                tl_A(i,j,k)=self%ci3d_p(i,j,k)
              END DO
            END DO
          END DO LEVEL_LOOP3

        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci3d_p(i,j,k)=0.0_r8
              self%ci3d_r(i,j,k)=0.0_r8
              self%ci3d_q(i,j,k)=0.0_r8
            END DO
          END DO
!
!  Update control vector 3D variable with the implicit diffusion
!  solution iterate.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j,k)=self%ci3d_x(i,j,k)/tl_scale(i,j)
            END DO
          END DO
        END DO LEVEL_LOOP4
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_3d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, utilizing Chebyshev Iterations
!  (CI) to model the multiscale propagation of background-error
!  covariance for 3D variables in the control vector. It is needed to
!  enforce symmetry-preserving correlation functions. The algorithm
!  requires the extrema eigenvalues of matrix A, which are obtained
!  from the Conjugate Gradient (CG) solver.
!
      SUBROUTINE multiscale_CI_3d_ad (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
      integer                           :: itrc
!
      real (r8)                         :: dotr, deps
      real (r8), dimension(N(ng))       :: dotn
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_scale
      real (r8), dimension(0:NiterCI,N(ng))         :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng))       :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_3d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      ad_scale=0.0_r8
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ci3d_r=0.0_r8
      self%ci3d_p=0.0_r8
      self%ci3d_q=0.0_r8
      self%ci3d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extrema
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%u_eigen(:,:,ms,1)
          eigMax => self%u_eigen(:,:,ms,2)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%v_eigen(:,:,ms,1)
          eigMax => self%v_eigen(:,:,ms,2)
        CASE ('temp', 'salt')                  ! tracers
          itrc = tracer_index(TRIM(StateVarName(ifield)))
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%t_eigen(:,:,ms,1,itrc)
          eigMax => self%t_eigen(:,:,ms,2,itrc)
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
!  Initialize implicit diffusion operator for the Chebyshev Iterations
!  solver. Non-adjointable steps.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci3d_r(i,j,k)=0.0_r8
              self%ci3d_p(i,j,k)=0.0_r8
              self%ci3d_q(i,j,k)=0.0_r8
              self%ci3d_x(i,j,k)=0.0_r8
            END DO
          END DO
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!
!  Adjoint of update control vector 3D variable with the implicit
!  diffusion solution iterate.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            tl_A(i,j,k)=self%ci3d_x(i,j,k)/tl_scale(i,j)
!^
              self%ci3d_x(i,j,k)=ad_A(i,j,k)/ad_scale(i,j)
              ad_A(i,j,k)=0.0_r8
            END DO
          END DO
!
!  Reset iterative conjugate gradient arrays.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci3d_q(i,j,k)=0.0_r8
!^
              self%ci3d_q(i,j,k)=0.0_r8
!^            self%ci3d_r(i,j,k)=0.0_r8
!^
              self%ci3d_r(i,j,k)=0.0_r8
!^            self%ci3d_p(i,j,k)=0.0_r8
!^
              self%ci3d_p(i,j,k)=0.0_r8
            END DO
          END DO
        END DO LEVEL_LOOP4
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          LEVEL_LOOP3 : DO k=1,N(ng)
!                                                        step 16
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              tl_A(i,j,k)=self%ci3d_p(i,j,k)
!^
                self%ci3d_p(i,j,k)=self%ci3d_p(i,j,k)+ad_A(i,j,k)
                ad_A(i,j,k)=0.0_r8
!^              self%ci3d_p(i,j,k)=-self%ci3d_r(i,j,k)+                 &
!^   &                             ci_beta(iterCI+1,k)*                 &
!^   &                             self%ci3d_p(i,j,k)
!^
                self%ci3d_r(i,j,k)=self%ci3d_r(i,j,k)-                  &
     &                             self%ci3d_p(i,j,k)
                self%ci3d_p(i,j,k)=ci_beta(iterCI+1,k)*                 &
     &                             self%ci3d_p(i,j,k)
              END DO
            END DO
!                                                        steps 15, 14
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              self%ci3d_r(i,j,k)=self%ci3d_r(i,j,k)+                  &
!^   &                             ci_alpha(iterCI,k)*                  &
!^   &                             self%ci3d_q(i,j,k)
!^
                self%ci3d_q(i,j,k)=self%ci3d_q(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_r(i,j,k)
!^              self%ci3d_x(i,j,k)=self%ci3d_x(i,j,k)+                  &
!^   &                             ci_alpha(iterCI,k)*                  &
!^   &                             self%ci3d_p(i,j,k)
!^
                self%ci3d_p(i,j,k)=self%ci3d_p(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_x(i,j,k)
              END DO
            END DO
!                                                        step 13
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              self%ci3d_q(i,j,k)=tl_A(i,j,k)
!^
                ad_A(i,j,k)=ad_A(i,j,k)+self%ci3d_q(i,j,k)
                self%ci3d_q(i,j,k)=0.0_r8
!^              tl_A(i,j,k)=tl_A(i,j,k)*tl_scale(i,j)
!^
                ad_A(i,j,k)=ad_A(i,j,k)*ad_scale(i,j)
              END DO
            END DO

          END DO LEVEL_LOOP3
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r3dvar)
              CALL self%ad_Klap_r3d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               ad_A)
            CASE (u3dvar)
              CALL self%ad_Klap_u3d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               ad_A)
            CASE (v3dvar)
              CALL self%ad_Klap_v3d (ng, tile, model,                   &
     &                               ifield, ctype, ms, Lweak,          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               ad_A)
          END SELECT
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              tl_A(i,j,k)=tl_A(i,j,k)/tl_scale(i,j)
!^
                ad_A(i,j,k)=ad_A(i,j,k)/ad_scale(i,j)
              END DO
            END DO
          END DO LEVEL_LOOP2

        END DO KLAP_ITER
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci3d_p(i,j,k)=tl_A(i,j,k)            ! step 11
!^
              ad_A(i,j,k)=ad_A(i,j,k)+self%ci3d_p(i,j,k)
              self%ci3d_p(i,j,k)=0.0_r8
!^            self%ci3d_r(i,j,k)=-tl_A(i,j,k)           ! step 10
!^
              ad_A(i,j,k)=ad_A(i,j,k)-self%ci3d_r(i,j,k)
              self%ci3d_r(i,j,k)=0.0_r8
!^            self%ci3d_x(i,j,k)=0.0_r8                 ! step 9
!^
              self%ci3d_x(i,j,k)=0.0_r8
!^            tl_A(i,j,k)=tl_scale(i,j)*tl_A(i,j,k)
!^
              ad_A(i,j,k)=ad_scale(i,j)*ad_A(i,j,k)
            END DO
          END DO
        END DO LEVEL_LOOP1
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_3d_ad
#endif /* SOLVE3D */

#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for the boundary adjustment of 2D
!  variables in the control vector. It operates in conjunction with
!  its adjoint to enforce symmetry-preserving correlation functions.
!  The algorithm requires the extrema eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b1d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI, ifac,  &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
# ifdef MULTI_SCALE_DEBUG
      integer                           :: status
# endif
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(LBij:UBij)   :: tl_scale
      real (r8), dimension(0:NiterCI)   :: ci_alpha
      real (r8), dimension(0:NiterCI+1) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b1d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ciB1d_r=0.0_r8
      self%ciB1d_p=0.0_r8
      self%ciB1d_q=0.0_r8
      self%ciB1d_x=0.0_r8
      tl_scale=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%zeta_obc_eigen(:,ibry,ms,1)
          eigMax => self%zeta_obc_eigen(:,ibry,ms,2)
        CASE ('ubar','ubar_eastward')          ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%ubar_obc_eigen(:,ibry,ms,1)
          eigMax => self%ubar_obc_eigen(:,ibry,ms,2)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%vbar_obc_eigen(:,ibry,ms,1)
          eigMax => self%vbar_obc_eigen(:,ibry,ms,2)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              tl_A(j)=tl_scale(j)*tl_A(j)
              self%ciB1d_x(j)=0.0_r8                    ! step 9
              self%ciB1d_r(j)=-tl_A(j)                  ! step 10
              self%ciB1d_p(j)=tl_A(j)                   ! step 11
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              tl_A(i)=tl_scale(i)*tl_A(i)
              self%ciB1d_x(i)=0.0_r8                    ! step 9
              self%ciB1d_r(i)=-tl_A(i)                  ! step 10
              self%ciB1d_p(i)=tl_A(i)                   ! step 11
            END DO
        END SELECT

# ifdef MULTI_SCALE_DEBUG
!
        dotn = dot_prod_1d (ng, tile, model, ctype, ibry,               &
     &                      LBij, UBij,                                 &
     &                      self%ciB1d_r,                               &
     &                      self%ciB1d_r)
# endif
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                tl_A(j)=tl_A(j)/tl_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                tl_A(i)=tl_A(i)/tl_scale(i)
              END DO
          END SELECT
!
          CALL self%tl_Klap_b1d (ng, tile, model,                       &
     &                           ifield, ibry, ctype, ms,               &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_A)                  ! step 12
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                tl_A(j)=tl_A(j)*tl_scale(j)
                self%ciB1d_q(j)=tl_A(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                tl_A(i)=tl_A(i)*tl_scale(i)
                self%ciB1d_q(i)=tl_A(i)
              END DO
          END SELECT

# ifdef MULTI_SCALE_DEBUG
!
          dotr = dot_prod_1d (ng, tile, model, ctype, ibry,             &
     &                        LBij, UBij,                               &
     &                        self%ciB1d_r,                             &
     &                        self%ciB1d_r)
          deps=SQRT(dotr/dotn)
!
          status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,       &
     &                            iterCI, NiterCI, 'CI_B1d', deps,      &
     &                            boundary = ibry)
# endif
!                                                         steps 14, 15
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%ciB1d_x(j)=self%ciB1d_x(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_p(j)
                self%ciB1d_r(j)=self%ciB1d_r(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_q(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%ciB1d_x(i)=self%ciB1d_x(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_p(i)
                self%ciB1d_r(i)=self%ciB1d_r(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_q(i)
              END DO
          END SELECT
!                                                         step 16
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%ciB1d_p(j)=-self%ciB1d_r(j)+                       &
     &                          ci_beta(iterCI+1)*self%ciB1d_p(j)
                tl_A(j)=self%ciB1d_p(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%ciB1d_p(i)=-self%ciB1d_r(i)+                       &
     &                          ci_beta(iterCI+1)*self%ciB1d_p(i)
                tl_A(i)=self%ciB1d_p(i)
              END DO
          END SELECT
!
        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              self%ciB1d_p(j)=0.0_r8
              self%ciB1d_r(j)=0.0_r8
              self%ciB1d_q(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              self%ciB1d_p(i)=0.0_r8
              self%ciB1d_r(i)=0.0_r8
              self%ciB1d_q(i)=0.0_r8
            END DO
        END SELECT
!
!  Update control vector lateral boundary variable with the implicit
!  diffusion solution iterate.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              tl_A(j)=self%ciB1d_x(j)/tl_scale(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              tl_A(i)=self%ciB1d_x(i)/tl_scale(i)
            END DO
        END SELECT
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_b1d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, using Chebyshev Iterations (CI)
!  to model the multiscale propagation of background-error covariance
!  for the boundary adjustment of 2D variables in the control vector.
!  The algorithm requires the extrema eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b1d_ad (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI, ifac,  &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: adfac, cff
      real (r8)                         :: ci_alpham1, ci_delta, ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(LBij:UBij)   :: ad_scale
      real (r8), dimension(0:NiterCI)   :: ci_alpha
      real (r8), dimension(0:NiterCI+1) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b1d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      ad_scale=0.0_r8
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ciB1d_r=0.0_r8
      self%ciB1d_p=0.0_r8
      self%ciB1d_q=0.0_r8
      self%ciB1d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%zeta_obc_eigen(:,ibry,ms,1)
          eigMax => self%zeta_obc_eigen(:,ibry,ms,2)
        CASE ('ubar','ubar_eastward')          ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%ubar_obc_eigen(:,ibry,ms,1)
          eigMax => self%ubar_obc_eigen(:,ibry,ms,2)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%vbar_obc_eigen(:,ibry,ms,1)
          eigMax => self%vbar_obc_eigen(:,ibry,ms,2)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              self%ciB1d_p(j)=0.0_r8
              self%ciB1d_r(j)=0.0_r8
              self%ciB1d_q(j)=0.0_r8
              self%ciB1d_x(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              self%ciB1d_p(i)=0.0_r8
              self%ciB1d_r(i)=0.0_r8
              self%ciB1d_q(i)=0.0_r8
              self%ciB1d_x(i)=0.0_r8
            END DO
        END SELECT
!
!  Adjoint of invert the implicit diffusion equation using Chebyshev
!  Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
!  Adjoint of update control vector lateral boundary variable with the
!  implicit diffusion solution iterate.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
!^            tl_A(j)=self%ciB1d_x(j)/tl_scale(j)
!^
              self%ciB1d_x(j)=ad_A(j)/ad_scale(j)
              ad_A(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
!^            tl_A(i)=self%ciB1d_x(i)/tl_scale(i)
!^
              self%ciB1d_x(i)=ad_A(i)/ad_scale(i)
              ad_A(i)=0.0_r8
            END DO
        END SELECT
!
!  Adjoint of reset iterative conjugate gradient arrays.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
!^            self%ciB1d_q(j)=0.0_r8
!^
              self%ciB1d_q(j)=0.0_r8
!^            self%ciB1d_r(j)=0.0_r8
!^
              self%ciB1d_r(j)=0.0_r8
!^            self%ciB1d_p(j)=0.0_r8
!^
              self%ciB1d_p(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
!^            self%ciB1d_q(i)=0.0_r8
!^
              self%ciB1d_q(i)=0.0_r8
!^            self%ciB1d_r(i)=0.0_r8
!^
              self%ciB1d_r(i)=0.0_r8
!^            self%ciB1d_p(i)=0.0_r8
!^
              self%ciB1d_p(i)=0.0_r8
            END DO
        END SELECT
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          SELECT CASE (ibry)                            ! step 16
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              tl_A(j)=self%ciB1d_p(j)
!^
                self%ciB1d_p(j)=self%ciB1d_p(j)+ad_A(j)
                ad_A(j)=0.0_r8
!^              self%ciB1d_p(j)=-self%ciB1d_r(j)+                       &
!^   &                          ci_beta(iterCI+1)*self%ciB1d_p(j)
!^
                self%ciB1d_r(j)=self%ciB1d_r(j)-self%ciB1d_p(j)
                self%ciB1d_p(j)=ci_beta(iterCI+1)*self%ciB1d_p(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              tl_A(i)=self%ciB1d_p(i)
!^
                self%ciB1d_p(i)=self%ciB1d_p(i)+ad_A(i)
                ad_A(i)=0.0_r8
!^              self%ciB1d_p(i)=-self%ciB1d_r(i)+                       &
!^   &                          ci_beta(iterCI+1)*self%ciB1d_p(i)
!^
                self%ciB1d_r(i)=self%ciB1d_r(i)-self%ciB1d_p(i)
                self%ciB1d_p(i)=ci_beta(iterCI+1)*self%ciB1d_p(i)
              END DO
          END SELECT
!                                                         steps 15, 14
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              self%ciB1d_r(j)=self%ciB1d_r(j)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_q(j)
!^
                self%ciB1d_q(j)=self%ciB1d_q(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_r(j)
!^              self%ciB1d_x(j)=self%ciB1d_x(j)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_p(j)
!^
                self%ciB1d_p(j)=self%ciB1d_p(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_x(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              self%ciB1d_r(i)=self%ciB1d_r(i)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_q(i)
!^
                self%ciB1d_q(i)=self%ciB1d_q(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_r(i)
!^              self%ciB1d_x(i)=self%ciB1d_x(i)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_p(i)
!^
                self%ciB1d_p(i)=self%ciB1d_p(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_x(i)
              END DO
          END SELECT
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              self%ciB1d_q(j)=tl_A(j)
!^
                ad_A(j)=ad_A(j)+self%ciB1d_q(j)
                self%ciB1d_q(j)=0.0_r8
!^              tl_A(j)=tl_A(j)*tl_scale(j)
!^
                ad_A(j)=ad_A(j)*ad_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              self%ciB1d_q(i)=tl_A(i)
!^
                ad_A(i)=ad_A(i)+self%ciB1d_q(i)
                self%ciB1d_q(i)=0.0_r8
!^              tl_A(i)=tl_A(i)*tl_scale(i)
!^
                ad_A(i)=ad_A(i)*ad_scale(i)
              END DO
          END SELECT
!
          CALL self%ad_Klap_b1d (ng, tile, model,                       &
     &                           ifield, ibry, ctype, ms,               &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_A)                  ! step 12
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              tl_A(j)=tl_A(j)/tl_scale(j)
!^
                ad_A(j)=ad_A(j)/ad_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              tl_A(i)=tl_A(i)/tl_scale(i)
!^
                ad_A(i)=ad_A(i)/ad_scale(i)
              END DO
          END SELECT
!
        END DO KLAP_ITER
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
!^            self%ciB1d_p(j)=tl_A(j)                   ! step 11
!^
              ad_A(j)=ad_A(j)+self%ciB1d_p(j)
              self%ciB1d_p(j)=0.0_r8
!^            self%ciB1d_r(j)=-tl_A(j)                  ! step 10
!^
              ad_A(j)=ad_A(j)-self%ciB1d_r(j)
!^            self%ciB1d_x(j)=0.0_r8                    ! step 9
!^
              self%ciB1d_x(j)=0.0_r8
!^            tl_A(j)=tl_scale(j)*tl_A(j)
!^
              ad_A(j)=ad_scale(j)*ad_A(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
!^            self%ciB1d_p(i)=tl_A(i)                   ! step 11
!^
              ad_A(i)=ad_A(i)+self%ciB1d_p(i)
              self%ciB1d_p(i)=0.0_r8
!^            self%ciB1d_r(i)=-tl_A(i)                  ! step 10
!^
              ad_A(i)=ad_A(i)-self%ciB1d_r(i)
!^            self%ciB1d_x(i)=0.0_r8                    ! step 9
!^
              self%ciB1d_x(i)=0.0_r8
!^            tl_A(i)=tl_scale(i)*tl_A(i)
!^
              ad_A(i)=ad_scale(i)*ad_A(i)
            END DO
        END SELECT
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_b1d_ad

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for the boundary adjustment of 3D
!  variables in the control vector. It operates in conjunction with
!  its adjoint to enforce symmetry-preserving correlation functions.
!  The algorithm requires the extrema eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b2d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI, ifac,  &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
      integer                           :: itrc
#  ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#  endif
!
      real (r8)                         :: dotr, deps
      real (r8), dimension(N(ng))       :: dotn
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(LBij:UBij)         :: tl_scale
      real (r8), dimension(0:NiterCI,N(ng))   :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng)) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b2d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ciB2d_r=0.0_r8
      self%ciB2d_p=0.0_r8
      self%ciB2d_q=0.0_r8
      self%ciB2d_x=0.0_r8
      tl_scale=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%u_obc_eigen(:,:,ibry,ms,1)
          eigMax => self%u_obc_eigen(:,:,ibry,ms,2)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%v_obc_eigen(:,:,ibry,ms,1)
          eigMax => self%v_obc_eigen(:,:,ibry,ms,2)
        CASE ('temp', 'salt')                  ! tracers
          itrc = tracer_index(TRIM(StateVarName(ifield)))
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%t_obc_eigen(:,:,ibry,ms,1,itrc)
          eigMax => self%t_obc_eigen(:,:,ibry,ms,2,itrc)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 :  DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                tl_A(j,k)=tl_scale(j)*tl_A(j,k)
                self%ciB2d_x(j,k)=0.0_r8               ! step 9
                self%ciB2d_r(j,k)=-tl_A(j,k)           ! step 10
                self%ciB2d_p(j,k)=tl_A(j,k)            ! step 11
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                tl_A(i,k)=tl_scale(i)*tl_A(i,k)
                self%ciB2d_x(i,k)=0.0_r8               ! step 9
                self%ciB2d_r(i,k)=-tl_A(i,k)           ! step 10
                self%ciB2d_p(i,k)=tl_A(i,k)            ! step 11
              END DO
          END SELECT

#  ifdef MULTI_SCALE_DEBUG
!
          dotn(k) = dot_prod_1d (ng, tile, model, ctype, ibry,          &
     &                           LBij, UBij,                            &
     &                           self%ciB2d_r(:,k),                     &
     &                           self%ciB2d_r(:,k))
#  endif
        END DO LEVEL_LOOP1
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  tl_A(j,k)=tl_A(j,k)/tl_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  tl_A(i,k)=tl_A(i,k)/tl_scale(i)
                END DO
            END SELECT
          END DO LEVEL_LOOP2
!
          CALL self%tl_Klap_b2d (ng, tile, model,                       &
     &                           ifield, ibry, ctype, ms,               &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_A)                  ! step 12
!
          LEVEL_LOOP3 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  tl_A(j,k)=tl_A(j,k)*tl_scale(j)
                  self%ciB2d_q(j,k)=tl_A(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  tl_A(i,k)=tl_A(i,k)*tl_scale(i)
                  self%ciB2d_q(i,k)=tl_A(i,k)
                END DO
            END SELECT

#  ifdef MULTI_SCALE_DEBUG
!
            dotr = dot_prod_1d (ng, tile, model, ctype, ibry,           &
     &                          LBij, UBij,                             &
     &                          self%ciB2d_r(:,k),                      &
     &                          self%ciB2d_r(:,k))
            deps=SQRT(dotr/dotn(k))
!
            status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,     &
     &                              iterCI, NiterCI, 'CI_B2d', deps,    &
     &                              level = k, boundary = ibry)
#  endif
!                                                         steps 14, 15
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  self%ciB2d_x(j,k)=self%ciB2d_x(j,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_p(j,k)
                  self%ciB2d_r(j,k)=self%ciB2d_r(j,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_q(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  self%ciB2d_x(i,k)=self%ciB2d_x(i,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_p(i,k)
                  self%ciB2d_r(i,k)=self%ciB2d_r(i,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_q(i,k)
                END DO
            END SELECT
!                                                         step 16
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  self%ciB2d_p(j,k)=-self%ciB2d_r(j,k)+                 &
     &                              ci_beta(iterCI+1,k)*                &
     &                              self%ciB2d_p(j,k)
                  tl_A(j,k)=self%ciB2d_p(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  self%ciB2d_p(i,k)=-self%ciB2d_r(i,k)+                 &
     &                              ci_beta(iterCI+1,k)*                &
     &                              self%ciB2d_p(i,k)
                    tl_A(i,k)=self%ciB2d_p(i,k)
                END DO
            END SELECT
          END DO LEVEL_LOOP3
!
        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%ciB2d_p(j,k)=0.0_r8
                self%ciB2d_r(j,k)=0.0_r8
                self%ciB2d_q(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%ciB2d_p(i,k)=0.0_r8
                self%ciB2d_r(i,k)=0.0_r8
                self%ciB2d_q(i,k)=0.0_r8
              END DO
          END SELECT
!
!  Update control vector lateral boundary variable with the implicit
!  diffusion solution iterate.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                tl_A(j,k)=self%ciB2d_x(j,k)/tl_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                tl_A(i,k)=self%ciB2d_x(i,k)/tl_scale(i)
              END DO
          END SELECT
        END DO LEVEL_LOOP4
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_b2d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, using Chebyshev Iterations (CI)
!  to model the multiscale propagation of background-error covariance
!  for the boundary adjustment of 3D variables in the control vector.
!  The algorithm requires the extrema eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b2d_ad (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI, ifac,  &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
      integer                           :: itrc
!
      real (r8)                         :: dotr, deps
      real (r8)                         :: adfac, cff
      real (r8)                         :: ci_alpham1, ci_delta, ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(LBij:UBij)         :: ad_scale
      real (r8), dimension(0:NiterCI,N(ng))   :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng)) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b2d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      ad_scale=0.0_r8
      ci_alpha=0.0_r8
      ci_beta=0.0_r8
      self%ciB2d_r=0.0_r8
      self%ciB2d_p=0.0_r8
      self%ciB2d_q=0.0_r8
      self%ciB2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
!  Note that Mlap/ifac iterations are used in this approach. When
!  ifac=2, the pseudo-diffusion operator is applied for only half
!  of the iterations, resulting in a square-root smoothing filter.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%u_obc_eigen(:,:,ibry,ms,1)
          eigMax => self%u_obc_eigen(:,:,ibry,ms,2)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%v_obc_eigen(:,:,ibry,ms,1)
          eigMax => self%v_obc_eigen(:,:,ibry,ms,2)
        CASE ('temp', 'salt')                  ! tracers
          itrc = tracer_index(TRIM(StateVarName(ifield)))
          Mlap=self%Mlap(ifield,ms)/ifac
          eigMin => self%t_obc_eigen(:,:,ibry,ms,1,itrc)
          eigMax => self%t_obc_eigen(:,:,ibry,ms,2,itrc)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
        LEVEL_LOOP4 :  DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%ciB2d_p(j,k)=0.0_r8
                self%ciB2d_r(j,k)=0.0_r8
                self%ciB2d_q(j,k)=0.0_r8
                self%ciB2d_x(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%ciB2d_p(i,k)=0.0_r8
                self%ciB2d_r(i,k)=0.0_r8
                self%ciB2d_q(i,k)=0.0_r8
                self%ciB2d_x(i,k)=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of invert the implicit diffusion equation using Chebyshev
!  Iterations.
!
!     A x = b  linear system
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!
!  Adjoint of update control vector lateral boundary variable with the
!  implicit diffusion solution iterate.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              tl_A(j,k)=self%ciB2d_x(j,k)/tl_scale(j)
!^
                self%ciB2d_x(j,k)=ad_A(j,k)/ad_scale(j)
                ad_A(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              tl_A(i,k)=self%ciB2d_x(i,k)/tl_scale(i)
!^
                self%ciB2d_x(i,k)=ad_A(i,k)/ad_scale(i)
                ad_A(i,k)=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of reset iterative conjugate gradient arrays.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              self%ciB2d_q(j,k)=0.0_r8
!^
                self%ciB2d_q(j,k)=0.0_r8
!^              self%ciB2d_r(j,k)=0.0_r8
!^
                self%ciB2d_r(j,k)=0.0_r8
!^              self%ciB2d_p(j,k)=0.0_r8
!^
                self%ciB2d_p(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              self%ciB2d_q(i,k)=0.0_r8
!^
                self%ciB2d_q(i,k)=0.0_r8
!^              self%ciB2d_r(i,k)=0.0_r8
!^
                self%ciB2d_r(i,k)=0.0_r8
!^              self%ciB2d_p(i,k)=0.0_r8
!^
                self%ciB2d_p(i,k)=0.0_r8
              END DO
          END SELECT
        END DO LEVEL_LOOP4
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          LEVEL_LOOP3 : DO k=1,N(ng)
!                                                         step 16
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
!^                tl_A(j,k)=self%ciB2d_p(j,k)
!^
                  self%ciB2d_p(j,k)=self%ciB2d_p(j,k)+ad_A(j,k)
                  ad_A(j,k)=0.0_r8
!^                self%ciB2d_p(j,k)=-self%ciB2d_r(j,k)+                 &
!^   &                              ci_beta(iterCI+1,k)*                &
!&   &                              self%ciB2d_p(j,k)
!^
                  self%ciB2d_r(j,k)=self%ciB2d_r(j,k)-                  &
     &                              self%ciB2d_p(j,k)
                  self%ciB2d_p(j,k)=ci_beta(iterCI+1,k)*                &
     &                              self%ciB2d_p(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
!^                tl_A(i,k)=self%ciB2d_p(i,k)
!^
                  self%ciB2d_p(i,k)=self%ciB2d_p(i,k)+ad_A(i,k)
                  ad_A(i,k)=0.0_r8
!^                self%ciB2d_p(i,k)=-self%ciB2d_r(i,k)+                 &
!^   &                              ci_beta(iterCI+1,k)*                &
!^   &                              self%ciB2d_p(i,k)
!^
                  self%ciB2d_r(i,k)=self%ciB2d_r(i,k)-                  &
     &                              self%ciB2d_p(i,k)
                  self%ciB2d_p(i,k)=ci_beta(iterCI+1,k)*                &
     &                              self%ciB2d_p(i,k)
                END DO
            END SELECT
!                                                         steps 15, 14
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
!^                self%ciB2d_r(j,k)=self%ciB2d_r(j,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_q(j,k)
!^
                  self%ciB2d_q(j,k)=self%ciB2d_q(j,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_r(j,k)
!^                self%ciB2d_x(j,k)=self%ciB2d_x(j,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_p(j,k)
!^
                  self%ciB2d_p(j,k)=self%ciB2d_p(j,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_x(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
!^                self%ciB2d_r(i,k)=self%ciB2d_r(i,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_q(i,k)
!^
                  self%ciB2d_q(i,k)=self%ciB2d_q(i,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_r(i,k)
!^                self%ciB2d_x(i,k)=self%ciB2d_x(i,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_p(i,k)
!^
                  self%ciB2d_p(i,k)=self%ciB2d_p(i,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_x(i,k)
                END DO
            END SELECT
!
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
!^                self%ciB2d_q(j,k)=tl_A(j,k)
!^
                  ad_A(j,k)=ad_A(j,k)+self%ciB2d_q(j,k)
                  self%ciB2d_q(j,k)=0.0_r8
!^                tl_A(j,k)=tl_A(j,k)*tl_scale(j)
!^
                  ad_A(j,k)=ad_A(j,k)*ad_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
!^                self%ciB2d_q(i,k)=tl_A(i,k)
!^
                  ad_A(i,k)=ad_A(i,k)+self%ciB2d_q(i,k)
                  self%ciB2d_q(i,k)=0.0_r8
!^                tl_A(i,k)=tl_A(i,k)*tl_scale(i)
!^
                  ad_A(i,k)=ad_A(i,k)*ad_scale(i)
                END DO
            END SELECT
!
          END DO LEVEL_LOOP3
!
          CALL self%ad_Klap_b2d (ng, tile, model,                       &
     &                           ifield, ibry, ctype, ms,               &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           ad_A)                  ! step 12
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
!^                tl_A(j,k)=tl_A(j,k)/tl_scale(j)
!^
                  ad_A(j,k)=ad_A(j,k)/ad_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
!^                tl_A(i,k)=tl_A(i,k)/tl_scale(i)
!^
                  ad_A(i,k)=ad_A(i,k)/ad_scale(i)
                END DO
            END SELECT
          END DO LEVEL_LOOP2
!
        END DO KLAP_ITER
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
!^              self%ciB2d_p(j,k)=tl_A(j,k)             ! step 11
!^
                ad_A(j,k)=ad_A(j,k)+self%ciB2d_p(j,k)
                self%ciB2d_p(j,k)=0.0_r8
!^              self%ciB2d_r(j,k)=-tl_A(j,k)            ! step 10
!^
                ad_A(j,k)=ad_A(j,k)-self%ciB2d_r(j,k)
!^              self%ciB2d_x(j,k)=0.0_r8                ! step 9
!^
                self%ciB2d_x(j,k)=0.0_r8
!^              tl_A(j,k)=tl_scale(j)*tl_A(j,k)
!^
                ad_A(j,k)=ad_scale(j)*ad_A(j,k)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
!^              self%ciB2d_p(i,k)=tl_A(i,k)             ! step 11
!^
                ad_A(i,k)=ad_A(i,k)+self%ciB2d_p(i,k)
                self%ciB2d_p(i,k)=0.0_r8
!^              self%ciB2d_r(i,k)=-tl_A(i,k)            ! step 10
!^
                ad_A(i,k)=ad_A(i,k)-self%ciB2d_r(i,k)
!^              self%ciB2d_x(i,k)=0.0_r8                ! step 9
!^
                self%ciB2d_x(i,k)=0.0_r8
!^              tl_A(i,k)=tl_scale(i)*tl_A(i,k)
!^
                ad_A(i,k)=ad_scale(i)*ad_A(i,k)
              END DO
          END SELECT
        END DO LEVEL_LOOP1
!
      END DO KDIFF_ITER
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CI_b2d_ad

# endif /* SOLVE3D */
#endif /* ADJUST_BOUNDARY */
