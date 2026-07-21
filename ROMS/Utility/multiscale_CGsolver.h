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
** Routines to invert the implicit diffusion equation using Conjugate **
** Gradient (CG) iterations to determine their extrema eigenvalues.   **
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
**  <><> CLASS MULTISCALE:  Conjugate Gradient Solver <><><><><><><><><>
**
*/

!  This routine computes the extrema eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for 2D variables in the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator are
!  computed using a Conjugate Gradient (CG) approach to solve the
!  associated linear system, Ax = b. The eigenvalue spectrum of this
!  operator remains invariant for a fixed application grid and a given
!  value of K. Consequently, estimates can be precomputed via the
!  Lanczos formulation of the Conjugate Gradient method, initialized
!  with random vectors, as implemented in the multiscale_eigen.F
!  module.
!
      SUBROUTINE multiscale_CG_2d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCG, ifac, Lweak,  &
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
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j
      integer                           :: itrc
#ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#endif
!
      real (r8)                         :: dotn, dotpq, dotr, dotrd
#ifdef MULTI_SCALE_DEBUG
      real (r8)                         :: deps
#endif
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCG)               :: cg_a
      real (r8), dimension(0:NiterCG+1)             :: cg_b
      real (r8), dimension(NiterCG+1)               :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)         :: work
      real (r8), dimension(NiterCG)                 :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1)     :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_2d_tl"
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
      self%cg2d_r=0.0_r8
      self%cg2d_p=0.0_r8
      self%cg2d_q=0.0_r8
      self%cg2d_x=0.0_r8
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
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
#ifdef MASKING
              tl_A(i,j)=GRID(ng)%rmask(i,j)*tl_A(i,j)
#endif
            END DO
          END DO
        CASE (u2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
#ifdef MASKING
              tl_A(i,j)=GRID(ng)%umask(i,j)*tl_A(i,j)
#endif
            END DO
          END DO
        CASE (v2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
#ifdef MASKING
              tl_A(i,j)=GRID(ng)%vmask(i,j)*tl_A(i,j)
#endif
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion 2D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%cg2d_x(i,j)=0.0_r8                       ! step 1
            self%cg2d_r(i,j)=-tl_A(i,j)                   ! step 2
            self%cg2d_p(i,j)=tl_A(i,j)                    ! step 3
          END DO
        END DO
        cg_b(0)=0.0_r8                                    ! step 4
!
        dotn = dot_prod_2d (ng, tile, model, ctype,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      self%cg2d_r,                                &
     &                      self%cg2d_r)
!
!  Invert the implicit diffusion equation using CG iterationsc.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)/tl_scale(i,j)
            END DO
          END DO
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          SELECT CASE (ctype)
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
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)*tl_scale(i,j)
              self%cg2d_q(i,j)=tl_A(i,j)                  ! step 7
            END DO
          END DO
!
          dotr = dot_prod_2d (ng, tile, model, ctype,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        self%cg2d_r,                              &
     &                        self%cg2d_r)
          dotpq = dot_prod_2d (ng, tile, model, ctype,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         self%cg2d_p,                             &
     &                         self%cg2d_q)
!
          cg_a(iterCG)=dotr/dotpq                         ! step 8

#ifdef MULTI_SCALE_DEBUG
!
          deps=SQRT(dotr/dotn)
          status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,       &
     &                            iterCG, NiterCG, 'CG_2d', deps)
#endif
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg2d_x(i,j)=self%cg2d_x(i,j)+                        &
     &                         cg_a(iterCG)*                            &
     &                         self%cg2d_p(i,j)           ! step 9
              self%cg2d_r(i,j)=self%cg2d_r(i,j)+                        &
     &                         cg_a(iterCG)*                            &
     &                         self%cg2d_q(i,j)           ! step 10
            END DO
          END DO
!
!  Reorthogonalize r(k+1) with respect to r(0:k).           step 11
!
          dotrd=dotr
          dotr = dot_prod_2d (ng, tile, model, ctype,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        self%cg2d_r,                              &
     &                        self%cg2d_r)
!
          cg_b(iterCG+1)=dotr/dotrd                       ! step 12
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg2d_p(i,j)=-self%cg2d_r(i,j)+                       &
     &                         cg_b(iterCG+1)*                          &
     &                         self%cg2d_p(i,j)           ! step 13
              tl_A(i,j)=self%cg2d_p(i,j)
            END DO
          END DO

        END DO KLAP_ITER
!
!  Update control vector 2D variable with the implicit diffusion
!  solution iterate.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j)=self%cg2d_x(i,j)
          END DO
        END DO
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and 
!  Ritz values of the tridiagonal matrix.
!
        DO i=1,NiterCG+1                                  ! step 14
          cg_Rv(i)=1.0_r8/cg_a(i-1)
          IF (i.gt.1) THEN
            cg_Rv(i)=cg_Rv(i)+cg_b(i-1)/cg_a(i-2)
          END IF
        END DO
        DO i=1,NiterCG
          zwork(i)=-SQRT(cg_b(i))/cg_a(i-1)
        END DO
!
        CALL DSTEQR ('I', NiterCG+1, cg_Rv, zwork,                      &
     &               zgv, NiterCG+1, work, info)
        IF (info.ne.0) THEN
          WRITE (stdout,*) ' multiscale_CG_2d_tl - Error in DSTEQR:',   &
     &                     ' info = ', info
          exit_flag=8
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        eigMax(iterDiff)=MAXVAL(cg_Rv)     ! maximum Ritz eigenvalue
        eigMin(iterDiff)=MINVAL(cg_Rv)     ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%cg2d_p(i,j)=0.0_r8
            self%cg2d_r(i,j)=0.0_r8
            self%cg2d_q(i,j)=0.0_r8
          END DO
        END DO
!
      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          tl_A(i,j)=self%cg2d_x(i,j)/tl_scale(i,j)
        END DO
      END DO
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CG_2d_tl

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  This routine computes the extrema eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for 3D variables in the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system.
!
      SUBROUTINE multiscale_CG_3d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCG, ifac, Lweak,  &
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
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j, k
      integer                           :: itrc
# ifdef MULTI_SCALE_DEBUG
      integer                           :: status
# endif
!
      real (r8)                         :: dotpq, dotr, dotrd, deps
      real (r8), dimension(N(ng))       :: dotn
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCG,N(ng))         :: cg_a
      real (r8), dimension(0:NiterCG+1,N(ng))       :: cg_b
      real (r8), dimension(NiterCG+1)               :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)         :: work
      real (r8), dimension(NiterCG)                 :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1)     :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_3d_tl"
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
      self%cg3d_r=0.0_r8
      self%cg3d_p=0.0_r8
      self%cg3d_q=0.0_r8
      self%cg3d_x=0.0_r8
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
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*                    &
     &                           GRID(ng)%on_r(i,j))
# ifdef MASKING
              DO k=1,N(ng)
                tl_A(i,j,k)=tl_A(i,j,k)*GRID(ng)%rmask(i,j)
              END DO
# endif
            END DO
          END DO
        CASE (u3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*                    &
     &                           GRID(ng)%on_u(i,j))
# ifdef MASKING
              DO k=1,N(ng)
                tl_A(i,j,k)=tl_A(i,j,k)*GRID(ng)%umask(i,j)
              END DO
# endif
            END DO
          END DO
        CASE (v3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*                    &
     &                           GRID(ng)%on_v(i,j))
# ifdef MASKING
              DO k=1,N(ng)
                tl_A(i,j,k)=tl_A(i,j,k)*GRID(ng)%vmask(i,j)
              END DO
# endif
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion 3D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg3d_x(i,j,k)=0.0_r8                    ! step 1
              self%cg3d_r(i,j,k)=-tl_A(i,j,k)              ! step 2
              self%cg3d_p(i,j,k)=tl_A(i,j,k)               ! step 3
            END DO
          END DO
          cg_b(0,k)=0.0_r8                                 ! step 4
!
          dotn(k) = dot_prod_2d (ng, tile, model, ctype,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           self%cg3d_r(:,:,k),                    &
     &                           self%cg3d_r(:,:,k))
        END DO LEVEL_LOOP1
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)/tl_scale(i,j)
              END DO
            END DO
          END DO LEVEL_LOOP2
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          SELECT CASE (ctype)
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
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          LEVEL_LOOP3 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)*tl_scale(i,j)
                self%cg3d_q(i,j,k)=tl_A(i,j,k)            ! step 7
              END DO
            END DO
!
            dotr = dot_prod_2d (ng, tile, model, ctype,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          self%cg3d_r(:,:,k),                     &
     &                          self%cg3d_r(:,:,k))
            dotpq = dot_prod_2d (ng, tile, model, ctype,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           self%cg3d_p(:,:,k),                    &
     &                           self%cg3d_q(:,:,k))
!
            cg_a(iterCG,k)=dotr/dotpq                     ! step 8

# ifdef MULTI_SCALE_DEBUG
!
            deps=SQRT(dotr/dotn(k))
            status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,     &
     &                              iterCG, NiterCG, 'CG_3d', deps,     &
     &                              level = k)
# endif
!
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%cg3d_x(i,j,k)=self%cg3d_x(i,j,k)+                  &
     &                             cg_a(iterCG,k)*                      &
     &                             self%cg3d_p(i,j,k)     ! step 9
                self%cg3d_r(i,j,k)=self%cg3d_r(i,j,k)+                  &
     &                             cg_a(iterCG,k)*                      &
     &                             self%cg3d_q(i,j,k)     ! step 10
              END DO
            END DO
!
!  Reorthogonalize r(k+1) with respect to r(0:k).
!
            dotrd=dotr
            dotr = dot_prod_2d (ng, tile, model, ctype,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          self%cg3d_r(:,:,k),                     &
     &                          self%cg3d_r(:,:,k))
!
            cg_b(iterCG+1,k)=dotr/dotrd                     ! step 12
!
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%cg3d_p(i,j,k)=-self%cg3d_r(i,j,k)+                 &
     &                             cg_b(iterCG+1,k)*                    &
     &                             self%cg3d_p(i,j,k)       ! step 13
                tl_A(i,j,k)=self%cg3d_p(i,j,k)
              END DO
            END DO
          END DO LEVEL_LOOP3

        END DO KLAP_ITER
!
!  Update control vector 3D variable with the implicit diffusion
!  solution iterate.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j,k)=self%cg3d_x(i,j,k)
            END DO
          END DO
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and
!  Ritz values of the tridiagonal matrix.
!
          DO i=1,NiterCG+1
            cg_Rv(i)=1.0_r8/cg_a(i-1,k)
            IF (i.gt.1) THEN
              cg_Rv(i)=cg_Rv(i)+cg_b(i-1,k)/cg_a(i-2,k)
            END IF
          END DO
          DO i=1,NiterCG
            zwork(i)=-SQRT(cg_b(i,k))/cg_a(i-1,k)
          END DO
!
          CALL DSTEQR ('I', NiterCG+1, cg_Rv, zwork,                    &
     &                 zgv, NiterCG+1, work, info)
          IF (info.ne.0) THEN
            WRITE (stdout,*) ' multiscale_CG_2d_tl - Error in DSTEQR:', &
     &                       ' info = ', info
            exit_flag=8
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          eigMax(iterDiff,k)=MAXVAL(cg_Rv) ! maximum Ritz eigenvalue
          eigMin(iterDiff,k)=MINVAL(cg_Rv) ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg3d_p(i,j,k)=0.0_r8
              self%cg3d_r(i,j,k)=0.0_r8
              self%cg3d_q(i,j,k)=0.0_r8
            END DO
          END DO
        END DO LEVEL_LOOP4

      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      DO k=1,N(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j,k)=self%cg3d_x(i,j,k)/tl_scale(i,j)
          END DO
        END DO
      END DO
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CG_3d_tl
#endif /* SOLVE3D */

#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  This routine computes the extrema eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for the boundary adjustments of 2D variables in
!  the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system, A x = b.
!
      SUBROUTINE multiscale_CG_b1d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCG, ifac,  &
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
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j, k
# ifdef MULTI_SCALE_DEBUG
      integer                           :: status
# endif
!
      real (r8)                         :: dotn, dotpq, dotr, dotrd
# ifdef MULTI_SCALE_DEBUG
      real (r8)                         :: deps
# endif
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(LBij:UBij)           :: tl_scale
      real (r8), dimension(0:NiterCG)           :: cg_a
      real (r8), dimension(0:NiterCG+1)         :: cg_b
      real (r8), dimension(NiterCG+1)           :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)     :: work
      real (r8), dimension(NiterCG)             :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1) :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_b1d_tl"
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
      self%cgB1d_r=0.0_r8
      self%cgB1d_p=0.0_r8
      self%cgB1d_q=0.0_r8
      self%cgB1d_x=0.0_r8
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
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
# ifdef MASKING
                 tl_A(j)=GRID(ng)%rmask(i,j)*tl_A(j)
# endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
# ifdef MASKING
                 tl_A(i)=GRID(ng)%rmask(i,j)*tl_A(i)
# endif
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
# ifdef MASKING
                 tl_A(j)=GRID(ng)%umask(i,j)*tl_A(j)
# endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
# ifdef MASKING
                 tl_A(i)=GRID(ng)%umask(i,j)*tl_A(i)
# endif
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
# ifdef MASKING
                 tl_A(j)=GRID(ng)%vmask(i,j)*tl_A(j)
# endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
# ifdef MASKING
                 tl_A(i)=GRID(ng)%vmask(i,j)*tl_A(i)
# endif
               END DO
             END IF
        END SELECT
      END IF

# ifdef DISTRIBUTE
!
      CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
     &                        LBij, UBij,                               &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        tl_A)
#  endif
!
!  Advance in K-space implicit pseudo-diffusion 2D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              self%cgB1d_x(j)=0.0_r8                      ! step 1
              self%cgB1d_r(j)=-tl_A(j)                    ! step 2
              self%cgB1d_p(j)=tl_A(j)                     ! step 3
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              self%cgB1d_x(i)=0.0_r8                      ! step 1
              self%cgB1d_r(i)=-tl_A(i)                    ! step 2
              self%cgB1d_p(i)=tl_A(i)                     ! step 3
            END DO
        END SELECT
        cg_b(0)=0.0_r8                                    ! step 4
!
        dotn = dot_prod_1d (ng, tile, model, ctype, ibry,               &
     &                      LBij, UBij,                                 &
     &                      self%cgB1d_r,                               &
     &                      self%cgB1d_r)
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
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
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          CALL self%tl_Klap_b1d (ng, tile, model,                       &
     &                           ifield, ibry, ctype, ms,               &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_A)
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                tl_A(j)=tl_A(j)*tl_scale(j)
                self%cgB1d_q(j)=tl_A(j)                   ! step 7
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                tl_A(i)=tl_A(i)*tl_scale(i)
                self%cgB1d_q(i)=tl_A(i)                   ! step 7
              END DO
          END SELECT
!
          dotr = dot_prod_1d (ng, tile, model, ctype, ibry,             &
     &                        LBij, UBij,                               &
     &                        self%cgB1d_r,                             &
     &                        self%cgB1d_r)
          dotpq = dot_prod_1d (ng, tile, model, ctype, ibry,            &
     &                         LBij, UBij,                              &
     &                         self%cgB1d_p,                            &
     &                         self%cgB1d_q)
!
          cg_a(iterCG)=dotr/dotpq                         ! step 8

# ifdef MULTI_SCALE_DEBUG
!
          deps=SQRT(dotr/dotn)
          status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,       &
     &                            iterCG, NiterCG, 'CG_b1d', deps,      &
     &                            boundary = ibry)
# endif
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%cgB1d_x(j)=self%cgB1d_x(j)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_p(j)           ! step 9
                self%cgB1d_r(j)=self%cgB1d_r(j)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_q(j)           ! step 10
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%cgB1d_x(i)=self%cgB1d_x(i)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_p(i)           ! step 9
                self%cgB1d_r(i)=self%cgB1d_r(i)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_q(i)           ! step 10
              END DO
          END SELECT
!
!  Reorthogonalize r(k+1) with respect to r(0:k).           step 11
!
          dotrd=dotr
          dotr = dot_prod_1d (ng, tile, model, ctype, ibry,             &
     &                        LBij, UBij,                               &
     &                        self%cgB1d_r,                             &
     &                        self%cgB1d_r)
!
          cg_b(iterCG+1)=dotr/dotrd                       ! step 12
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%cgB1d_p(j)=-self%cgB1d_r(j)+                       &
     &                          cg_b(iterCG+1)*                         &
     &                          self%cgB1d_p(j)           ! step 13
                tl_A(j)=self%cgB1d_p(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%cgB1d_p(i)=-self%cgB1d_r(i)+                       &
     &                          cg_b(iterCG+1)*                         &
     &                          self%cgB1d_p(i)           ! step 13
                tl_A(i)=self%cgB1d_p(i)
              END DO
          END SELECT
!
        END DO KLAP_ITER
!
!  Update control vector 2D variable lateral boundary consition with the
!  implicit diffusion solution iterate.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              tl_A(j)=self%cgB1d_x(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              tl_A(i)=self%cgB1d_x(i)
            END DO
        END SELECT
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and 
!  Ritz values of the tridiagonal matrix.
!
        DO i=1,NiterCG+1                                  ! step 14
          cg_Rv(i)=1.0_r8/cg_a(i-1)
          IF (i.gt.1) THEN
            cg_Rv(i)=cg_Rv(i)+cg_b(i-1)/cg_a(i-2)
          END IF
        END DO
        DO i=1,NiterCG
          zwork(i)=-SQRT(cg_b(i))/cg_a(i-1)
        END DO
!
        CALL DSTEQR ('I', NiterCG+1, cg_Rv, zwork,                      &
     &               zgv, NiterCG+1, work, info)
        IF (info.ne.0) THEN
          WRITE (stdout,*) ' multiscale_CG_b1d_tl - Error in DSTEQR:',  &
     &                     ' info = ', info
          exit_flag=8
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        eigMax(iterDiff)=MAXVAL(cg_Rv)  ! maximum Ritz eigenvalue
        eigMin(iterDiff)=MINVAL(cg_Rv)  ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              self%cgB1d_p(j)=0.0_r8
              self%cgB1d_r(j)=0.0_r8
              self%cgB1d_q(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              self%cgB1d_p(i)=0.0_r8
              self%cgB1d_r(i)=0.0_r8
              self%cgB1d_q(i)=0.0_r8
            END DO
        END SELECT
!
      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO j=Jmin,Jmax
            tl_A(j)=self%cgB1d_x(j)/tl_scale(j)
          END DO
        CASE (isouth, inorth)
          DO i=Imin,Imax
            tl_A(i)=self%cgB1d_x(i)/tl_scale(i)
          END DO
      END SELECT
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CG_b1d_tl

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  This routine computes the extrema eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for the boundary adjustments of 3D variables in
!  the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system, A x = b.
!
      SUBROUTINE multiscale_CG_b2d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCG, ifac,  &
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
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      integer,            intent(in   ) :: ifac      ! iteraction factor
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j, k
      integer                           :: itrc
#  ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#  endif
!
      real (r8)                         :: dotpq, dotr, dotrd, deps
      real (r8), dimension(N(ng))       :: dotn
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(LBij:UBij)           :: tl_scale
      real (r8), dimension(0:NiterCG,N(ng))     :: cg_a
      real (r8), dimension(0:NiterCG+1,N(ng))   :: cg_b
      real (r8), dimension(NiterCG+1)           :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)     :: work
      real (r8), dimension(NiterCG)             :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1) :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_b2d_tl"
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
      self%cgB2d_r=0.0_r8
      self%cgB2d_p=0.0_r8
      self%cgB2d_q=0.0_r8
      self%cgB2d_x=0.0_r8
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
#  ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(j,k)=GRID(ng)%rmask(i,j)*tl_A(j,k)
                 END DO
#  endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
#  ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(i,k)=GRID(ng)%rmask(i,j)*tl_A(i,k)
                 END DO
#  endif
               END DO
             END IF
           CASE (u3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
#  ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(j,k)=GRID(ng)%umask(i,j)*tl_A(j,k)
                 END DO
#  endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
#  ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(i,k)=GRID(ng)%umask(i,j)*tl_A(i,k)
                 END DO
#  endif
               END DO
             END IF
           CASE (v3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
#  ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(j,k)=GRID(ng)%vmask(i,j)*tl_A(j,k)
                 END DO
#  endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
#  ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(i,k)=GRID(ng)%vmask(i,j)*tl_A(i,k)
                 END DO
#  endif
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion 2D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%cgB2d_x(j,k)=0.0_r8                  ! step 1
                self%cgB2d_r(j,k)=-tl_A(j,k)              ! step 2
                self%cgB2d_p(j,k)=tl_A(j,k)               ! step 3
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%cgB2d_x(i,k)=0.0_r8                  ! step 1
                self%cgB2d_r(i,k)=-tl_A(i,k)              ! step 2
                self%cgB2d_p(i,k)=tl_A(i,k)               ! step 3
              END DO
          END SELECT
          cg_b(0,k)=0.0_r8                                ! step 4
!
          dotn(k) = dot_prod_1d (ng, tile, model, ctype, ibry,          &
     &                           LBij, UBij,                            &
     &                           self%cgB2d_r(:,k),                     &
     &                           self%cgB2d_r(:,k))
        END DO LEVEL_LOOP1
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
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
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          CALL self%tl_Klap_b2d (ng, tile, model,                       &
     &                           ifield, ibry, ctype, ms,               &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           tl_A)
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          LEVEL_LOOP3 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  tl_A(j,k)=tl_A(j,k)*tl_scale(j)
                  self%cgB2d_q(j,k)=tl_A(j,k)             ! step 7
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  tl_A(i,k)=tl_A(i,k)*tl_scale(i)
                  self%cgB2d_q(i,k)=tl_A(i,k)             ! step 7
                END DO
            END SELECT
!
            dotr = dot_prod_1d (ng, tile, model, ctype, ibry,           &
     &                          LBij, UBij,                             &
     &                          self%cgB2d_r(:,k),                      &
     &                          self%cgB2d_r(:,k))
            dotpq = dot_prod_1d (ng, tile, model, ctype, ibry,          &
     &                           LBij, UBij,                            &
     &                           self%cgB2d_p(:,k),                     &
     &                           self%cgB2d_q(:,k))
!
            cg_a(iterCG,k)=dotr/dotpq                     ! step 8

#  ifdef MULTI_SCALE_DEBUG
!
            deps=SQRT(dotr/dotn(k))
            status=multiscale_print(ng, ifield, ms, iterDiff, Mlap,     &
     &                              iterCG, NiterCG, 'CG_b2d', deps,    &
     &                              level = k, boundary = ibry)
#  endif
!
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  self%cgB2d_x(j,k)=self%cgB2d_x(j,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_p(j,k)     ! step 9
                  self%cgB2d_r(j,k)=self%cgB2d_r(j,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_q(j,k)     ! step 10
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  self%cgB2d_x(i,k)=self%cgB2d_x(i,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_p(i,k)     ! step 9
                  self%cgB2d_r(i,k)=self%cgB2d_r(i,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_q(i,k)     ! step 10
                END DO
            END SELECT
!
!  Reorthogonalize r(k+1) with respect to r(0:k).           step 11
!
            dotrd=dotr
            dotr = dot_prod_1d (ng, tile, model, ctype, ibry,           &
     &                          LBij, UBij,                             &
     &                          self%cgB2d_r(:,k),                      &
     &                          self%cgB2d_r(:,k))
!
            cg_b(iterCG+1,k)=dotr/dotrd                   ! step 12
!
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmax
                  self%cgB2d_p(j,k)=-self%cgB2d_r(j,k)+                 &
     &                              cg_b(iterCG+1,k)*                   &
     &                              self%cgB2d_p(j,k)     ! step 13
                  tl_A(j,k)=self%cgB2d_p(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imax
                  self%cgB2d_p(i,k)=-self%cgB2d_r(i,k)+                 &
     &                              cg_b(iterCG+1,k)*                   &
     &                              self%cgB2d_p(i,k)     ! step 13
                  tl_A(i,k)=self%cgB2d_p(i,k)
                END DO
            END SELECT
          END DO LEVEL_LOOP3
!
        END DO KLAP_ITER
!
!  Update control vector 2D variable with the implicit diffusion
!  solution iterate.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                tl_A(j,k)=self%cgB2d_x(j,k)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                tl_A(i,k)=self%cgB2d_x(i,k)
              END DO
          END SELECT
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and 
!  Ritz values of the tridiagonal matrix.
!
          DO i=1,NiterCG+1                                ! step 14
            cg_Rv(i)=1.0_r8/cg_a(i-1,k)
            IF (i.gt.1) THEN
              cg_Rv(i)=cg_Rv(i)+cg_b(i-1,k)/cg_a(i-2,k)
            END IF
          END DO
          DO i=1,NiterCG
            zwork(i)=-SQRT(cg_b(i,k))/cg_a(i-1,k)
          END DO
!
          CALL DSTEQR ('I', NiterCG+1, cg_Rv, zwork,                    &
     &                 zgv, NiterCG+1, work, info)
          IF (info.ne.0) THEN
            WRITE (stdout,*) ' multiscale_CG_b2d_tl - Error in DSTEQR:',&
     &                       ' info = ', info
            exit_flag=8
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          eigMax(iterDiff,k)=MAXVAL(cg_Rv)  ! maximum Ritz eigenvalue
          eigMin(iterDiff,k)=MINVAL(cg_Rv)  ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmax
                self%cgB2d_p(j,k)=0.0_r8
                self%cgB2d_r(j,k)=0.0_r8
                self%cgB2d_q(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imax
                self%cgB2d_p(i,k)=0.0_r8
                self%cgB2d_r(i,k)=0.0_r8
                self%cgB2d_q(i,k)=0.0_r8
              END DO
          END SELECT
        END DO LEVEL_LOOP4
!
      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              tl_A(j,k)=self%cgB2d_x(j,k)/tl_scale(j)
            END DO
          END DO
        CASE (isouth, inorth)
          DO k=1,N(ng)
            DO i=Imin,Imax
              tl_A(i,k)=self%cgB2d_x(i,k)/tl_scale(i)
            END DO
          END DO
      END SELECT
!
!  Nullify local pointers.
!
      nullify (eigMin)
      nullify (eigMax)
!
      RETURN
      END SUBROUTINE multiscale_CG_b2d_tl

# endif /* SOLVE3D */
#endif /* ADJUST_BOUNDARY */
