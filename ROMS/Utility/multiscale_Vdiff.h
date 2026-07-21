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
** Routines to advance the implicit vertical diffusion operator.      **
** Error correlations are assumed to be separable in the horizontal   **
** and vertical directions.                                           **
**                                                                    **
************************************************************************
**
**
**  <><> CLASS MULTISCALE:  Vertical Diffusion Solver <><><><><><><><><>
**
*/

!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix.
!  Note that error correlations are considered separable in the
!  horizontal and vertical directions.
!
      SUBROUTINE multiscale_Vdiff_r3d_tl (self, ng, tile, model,        &
     &                                    ifield, NVsteps,              &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV, Kv, tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: cff
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: tl_Awrk
#if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
#else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: tl_DC
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Nold=1
      Nnew=2
!
      FC=0.0_r8
      tl_Awrk=0.0_r8
      tl_DC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
          END DO
#else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-        &
     &                                    GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#endif
        END DO
      END DO
!
!  Set operator initial conditions.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    Nghostpoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Nghostpoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
!
#endif
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!^          Awrk(i,j,k,Nold)=A(i,j,k)
!^
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO

#ifdef SPLINES_VCONV
!
! <><><> Integrate vertical diffusion equation implicitly using
!        parabolic splines <><><>
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=Jstr,Jend
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                        &
     &                DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                        &
     &                DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
!^          DC(i,0)=0.0_r8
!^
            tl_DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i,j,k+1))+   &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
!^            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
!^   &                     FC(i,k)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
!^          DC(i,N(ng))=0.0_r8
!^
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!^            DC(i,k)=DC(i,k)*Kv(i,j,k)
!^
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!^   &                         DTsizeV*oHz(i,j,k)*                      &
!^   &                         (DC(i,k)-DC(i,k-1))
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

#else
!
! <><><> Integrate vertical diffusion equation implicitly <><><>
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tangent diffusion equation into tl_DC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,j,k)-FC(i,j,k-1)
!^            DC(i,k)=Awrk(i,j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*GRID(ng)%Hz(i,j,k)
            END DO
          END DO
!
!  Solve the tangent linear tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!^          DC(i,1)=cff*DC(i,1)
!^
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!^            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
!^          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!^   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!^   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
# ifdef MASKING
!^          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*                  &
!^   &                           GRID(ng)%rmask(i,j)
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*            &
     &                              GRID(ng)%rmask(i,j)
# endif
          END DO
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^            Awrk(i,j,k,Nnew)=DC(i,k)
!^
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
# ifdef MASKING
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*                        &
!^   &                         GRID(ng)%rmask(i,j)
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*                  &
     &                            GRID(ng)%rmask(i,j)
# endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#endif
!
!  Load convolved solution.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          A(i,j,k)=Awrk(i,j,k,Nold)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

#ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_r3d_tl
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion adjoint operator using an
!  implicit solver using parabolic splines or inverting a tridiagonal
!  matrix. Note that error correlations are considered separable in the
!  horizontal and vertical directions.
!
      SUBROUTINE multiscale_Vdiff_r3d_ad (self, ng, tile, model,        &
     &                                    ifield, NVsteps,              &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV, Kv, ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: adfac, cff
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: ad_Awrk
#if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
#else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: ad_DC
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Nold=1
      Nnew=2
!
      ad_Awrk=0.0_r8
      ad_DC=0.0_r8
      FC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
          END DO
#else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-        &
     &                                    GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#endif
        END DO
      END DO
!
!  Adjoint of load convolved solution.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#endif
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
!^
            ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

#ifdef SPLINES_VCONV
!
! <><><> Integrate adjoint vertical diffusion equation implicitly
!        using parabolic splines <><><>
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
        DO j=Jstr,Jend
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives. Then, time step adjoint vertical diffusion
!  equation implicitly.
!
!  Compute BASIC STATE time-invariant coefficients.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                        &
     &                DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                        &
     &                DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
          END DO
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i,j,k+1))+   &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
            END DO
          END DO
!
!  Adjoint of backward substitution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
!^   &                            DTsizeV*oHz(i,j,k)*                   &
!^   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!^
              adfac=DTsizeV*oHz(i,j,k)*ad_Awrk(i,j,k,Nnew)
              ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
              ad_DC(i,k  )=ad_DC(i,k  )+adfac
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+                  &
     &                            ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^
              ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
            END DO
          END DO
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
          DO i=Istr,Iend
!^          tl_DC(i,N(ng))=0.0_r8
!^
            ad_DC(i,N(ng))=0.0_r8
          END DO
!
!  Adjoint of LU decomposition and forward substitution.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
!^   &                        tl_Awrk(i,j,k  ,Nold)-                    &
!^   &                        FC(i,k)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_Awrk(i,j,k  ,Nold)=ad_Awrk(i,j,k  ,Nold)-adfac
              ad_Awrk(i,j,k+1,Nold)=ad_Awrk(i,j,k+1,Nold)+adfac
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
!
          DO i=Istr,Iend
!^          tl_DC(i,0)=0.0_r8
!^
            ad_DC(i,0)=0.0_r8
          END DO
        END DO
      END DO

#else

!
! <><><> Integrate adjoint vertical diffusion equation implicitly <><><>
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
        DO j=Jstr,Jend
!
!  Adjoint of compute diagonal matrix coefficients BC.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,j,k)-FC(i,j,k-1)
            END DO
          END DO
!
!  Adjoint of compute new solution by back substitution.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
            END DO
          END DO
!
!^        DO k=N(ng)-1,1,-1
!^
          DO k=1,N(ng)-1
            DO i=Istr,Iend
# ifdef MASKING
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*                  &
!^   &                            GRID(ng)%rmask(i,j)
!^
              ad_Awrk(i,j,k,Nnew)=ad_Awrk(i,j,k,Nnew)*                  &
     &                            GRID(ng)%rmask(i,j)
# endif
!^            tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
!^
              ad_DC(i,k)=ad_DC(i,k)+                                    &
     &                   ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
!
          DO i=Istr,Iend
# ifdef MASKING
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*            &
!^   &                              GRID(ng)%rmask(i,j)
!^
            ad_Awrk(i,j,N(ng),Nnew)=ad_Awrk(i,j,N(ng),Nnew)*            &
     &                              GRID(ng)%rmask(i,j)
# endif
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
!^
            ad_DC(i,N(ng))=ad_DC(i,N(ng))+                              &
     &                     ad_Awrk(i,j,N(ng),Nnew)
            ad_Awrk(i,j,N(ng),Nnew)=0.0_r8
!^          tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
!^   &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!^   &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            adfac=ad_DC(i,N(ng))/                                       &
     &            (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-FC(i,j,N(ng)-1)*adfac
            ad_DC(i,N(ng)  )=adfac
          END DO
!
!  Adjoint of solve tridiagonal system.
!
!^        DO k=2,N(ng)-1
!^
          DO k=N(ng)-1,2,-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,j,k-1)*adfac
              ad_DC(i,k  )=adfac
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
!^          tl_DC(i,1)=cff*tl_DC(i,1)
!^
            ad_DC(i,1)=cff*ad_DC(i,1)
          END DO
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!^            tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+                  &
     &                            GRID(ng)%Hz(i,j,k)*ad_DC(i,k)
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
        END DO
      END DO
#endif
!
!  Adjoint of set operator initial conditions.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k,Nold)
            ad_Awrk(i,j,k,Nold)=0.0_r8
          END DO
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    Nghostpoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    ad_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       Nghostpoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_r3d_ad
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix.
!  Note that error correlations are considered separable in the
!  horizontal and vertical directions.
! 
      SUBROUTINE multiscale_Vdiff_u3d_tl (self, ng, tile, model,        &
     &                                    ifield, NVsteps,              &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV, Kv, tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: IstrU, Iend, Jstr, Jend
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: cff
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: tl_Awrk
#if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
      real(r8), dimension(IminS:ImaxS,N(ng))               :: Hzk
#else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: tl_DC
!
!  Initialize.
!
      IstrU=BOUNDS(ng)%IstrU(tile)   ! tile computational range
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Nold=1
      Nnew=2
!
      FC=0.0_r8
      tl_Awrk=0.0_r8
      tl_DC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
#ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=2.0_r8/(GRID(ng)%Hz(i-1,j,k)+                    &
     &                         GRID(ng)%Hz(i  ,j,k))
          END DO
#else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                 &
     &                 (GRID(ng)%z_r(i-1,j,k+1)+GRID(ng)%z_r(i,j,k+1)-  &
     &                  GRID(ng)%z_r(i-1,j,k  )-GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#endif
        END DO
      END DO
!
!  Set operator initial conditions.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    Nghostpoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Nghostpoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
!
#endif
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
!^          Awrk(i,j,k,Nold)=A(i,j,k)
!^
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO

#ifdef SPLINES_VCONV
!
! <><><> Integrate vertical diffusion equation implicitly using
!        parabolic splines <><><>
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              Hzk(i,k)=0.5_r8*(GRID(ng)%Hz(i-1,j,k)+                    &
     &                         GRID(ng)%Hz(i  ,j,k))
            END DO
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              FC(i,k)=cff1*Hzk(i,k  )-DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=IstrU,Iend
            CF(i,0)=0.0_r8
!^          DC(i,0)=0.0_r8
!^
            tl_DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
!^            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
!^   &                     FC(i,k)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=IstrU,Iend
!^          DC(i,N(ng))=0.0_r8
!^
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
!^            DC(i,k)=DC(i,k)*Kv(i,j,k)
!^
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!^   &                         DTsizeV*oHz(i,j,k)*                      &
!^   &                         (DC(i,k)-DC(i,k-1))
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

#else

!
! <><><> Integrate vertical diffusion equation implicitly <><><>
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tangent diffusion equation into tl_DC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff=0.5_r8*(GRID(ng)%Hz(i-1,j,k)+GRID(ng)%Hz(i,j,k))
              BC(i,k)=cff-FC(i,j,k)-FC(i,j,k-1)
!^            DC(i,k)=Awrk(i,j,k,Nold)*cff
!^
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*cff
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=IstrU,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!^          DC(i,1)=cff*DC(i,1)
!^
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=IstrU,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!^            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=IstrU,Iend
!^          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!^   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!^   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
# ifdef MASKING
!^          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*                  &
!^   &                           GRID(ng)%umask(i,j)
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*            &
     &                              GRID(ng)%umask(i,j)
# endif
          END DO
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^            Awrk(i,j,k,Nnew)=DC(i,k)
!^
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
# ifdef MASKING
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*                        &
!^   &                         GRID(ng)%umask(i,j)
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*                  &
     &                            GRID(ng)%umask(i,j)
# endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#endif
!
!  Load convolved solution.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!^          A(i,j,k)=Awrk(i,j,k,Nold)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_u3d_tl
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix.
!  Note that error correlations are considered separable in the
!  horizontal and vertical directions.
! 
      SUBROUTINE multiscale_Vdiff_u3d_ad (self, ng, tile, model,        &
     &                                    ifield, NVsteps,              &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV, Kv, ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: IstrU, Iend, Jstr, Jend
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: adfac, cff
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: ad_Awrk
#if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
      real(r8), dimension(IminS:ImaxS,N(ng))               :: Hzk
#else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: ad_DC
!
!  Initialize.
!
      IstrU=BOUNDS(ng)%IstrU(tile)   ! tile computational range
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Nold=1
      Nnew=2
!
      ad_Awrk=0.0_r8
      ad_DC=0.0_r8
      FC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
#ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=2.0_r8/(GRID(ng)%Hz(i-1,j,k)+                    &
     &                         GRID(ng)%Hz(i  ,j,k))
          END DO
#else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*(Kv(i-1,j,k)+Kv(i,j,k))/                 &
     &                 (GRID(ng)%z_r(i-1,j,k+1)+GRID(ng)%z_r(i,j,k+1)-  &
     &                  GRID(ng)%z_r(i-1,j,k  )-GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#endif
        END DO
      END DO
!
!  Adjoint of load convolved solution.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#endif
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!^          tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
!^
            ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

#ifdef SPLINES_VCONV
!
! <><><> Integrate adjoint vertical diffusion equation implicitly
!        using parabolic splines <><><>
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute BASIC STATE time-invariant coefficients.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              Hzk(i,k)=0.5_r8*(GRID(ng)%Hz(i-1,j,k)+                    &
     &                         GRID(ng)%Hz(i  ,j,k))
            END DO
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              FC(i,k)=cff1*Hzk(i,k  )-DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=IstrU,Iend
            CF(i,0)=0.0_r8
          END DO
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
            END DO
          END DO
!
!  Adjoint of backward substitution.
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
!^   &                            DTsizeV*oHz(i,j,k)*                   &
!^   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!^
              adfac=DTsizeV*oHz(i,j,k)*ad_Awrk(i,j,k,Nnew)
              ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
              ad_DC(i,k  )=ad_DC(i,k  )+adfac
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+                  &
     &                            ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^
              ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
            END DO
          END DO
!
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
          DO i=IstrU,Iend
!^          tl_DC(i,N(ng))=0.0_r8
!^
            ad_DC(i,N(ng))=0.0_r8
          END DO
!
!  Adjoint of LU decomposition and forward substitution.
!
          DO k=N(ng)-1,1,-1
            DO i=IstrU,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
!^   &                        tl_Awrk(i,j,k  ,Nold)-                    &
!^   &                        FC(i,k)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_Awrk(i,j,k  ,Nold)=ad_Awrk(i,j,k  ,Nold)-adfac
              ad_Awrk(i,j,k+1,Nold)=ad_Awrk(i,j,k+1,Nold)+adfac
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
          DO i=IstrU,Iend
!^          tl_DC(i,0)=0.0_r8
!^
            ad_DC(i,0)=0.0_r8
          END DO
        END DO
      END DO

#else

!
! <><><> Integrate adjoint vertical diffusion equation implicitly <><><>
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Adjoint of compute diagonal matrix coefficients BC.
!
        DO j=Jstr,Jend
          DO k=1,N(ng)
            DO i=IstrU,Iend
              BC(i,k)=0.5_r8*(GRID(ng)%Hz(i-1,j,k)+                     &
     &                        GRID(ng)%Hz(i  ,j,k))-                    &
     &                FC(i,j,k)-FC(i,j,k-1)
            END DO
          END DO
!
!  Adjoint of compute new solution by back substitution.
!
          DO i=IstrU,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
          END DO
          DO k=2,N(ng)-1
            DO i=IstrU,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
            END DO
          END DO
!
!^        DO k=N(ng)-1,1,-1
!^
          DO k=1,N(ng)-1
            DO i=IstrU,Iend
# ifdef MASKING
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*                  &
!^   &                            GRID(ng)%umask(i,j)
!^
              ad_Awrk(i,j,k,Nnew)=ad_Awrk(i,j,k,Nnew)*                  &
     &                            GRID(ng)%umask(i,j)
# endif
!^            tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
!^
              ad_DC(i,k)=ad_DC(i,k)+                                    &
     &                   ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
!
          DO i=IstrU,Iend
# ifdef MASKING
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*            &
!^   &                              GRID(ng)%umask(i,j)
!^
            ad_Awrk(i,j,N(ng),Nnew)=ad_Awrk(i,j,N(ng),Nnew)*            &
     &                              GRID(ng)%umask(i,j)
# endif
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
!^
            ad_DC(i,N(ng))=ad_DC(i,N(ng))+                              &
     &                     ad_Awrk(i,j,N(ng),Nnew)
            ad_Awrk(i,j,N(ng),Nnew)=0.0_r8
!^          tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
!^   &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!^   &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            adfac=ad_DC(i,N(ng))/                                       &
     &            (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-FC(i,j,N(ng)-1)*adfac
            ad_DC(i,N(ng)  )=adfac
          END DO
!
!  Adjoint of solve tridiagonal system.
!
!^        DO k=2,N(ng)-1
!^
          DO k=N(ng)-1,2,-1
            DO i=IstrU,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,j,k-1)*adfac
              ad_DC(i,k  )=adfac
            END DO
          END DO
          DO i=IstrU,Iend
            cff=1.0_r8/BC(i,1)
!^          tl_DC(i,1)=cff*tl_DC(i,1)
!^
            ad_DC(i,1)=cff*ad_DC(i,1)
          END DO
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          DO k=1,N(ng)
            DO i=IstrU,Iend
              cff=0.5*(GRID(ng)%Hz(i-1,j,k)+GRID(ng)%Hz(i,j,k))
!^            tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*cff
!^
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+cff*ad_DC(i,k)
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
        END DO
      END DO
#endif
!
!  Adjoint of set operator initial conditions.
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
!^          tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k,Nold)
            ad_Awrk(i,j,k,Nold)=0.0_r8
          END DO
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    Nghostpoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       Nghostpoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_u3d_ad
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix.
!  Note that error correlations are considered separable in the
!  horizontal and vertical directions.
! 
      SUBROUTINE multiscale_Vdiff_v3d_tl (self, ng, tile, model,        &
     &                                    ifield, NVsteps,              &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV, Kv, tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, Iend, JstrV, Jend
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: cff
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: tl_Awrk
#if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
      real(r8), dimension(IminS:ImaxS,N(ng))               :: Hzk
#else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: tl_DC
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Nold=1
      Nnew=2
!
      FC=0.0_r8
      tl_Awrk=0.0_r8
      tl_DC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=2.0_r8/(GRID(ng)%Hz(i,j-1,k)+                    &
     &                         GRID(ng)%Hz(i,j  ,k))
          END DO
#else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                 &
     &                 (GRID(ng)%z_r(i,j-1,k+1)+GRID(ng)%z_r(i,j,k+1)-  &
     &                  GRID(ng)%z_r(i,j-1,k  )-GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#endif
        END DO
      END DO
!
!  Set operator initial conditions.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    Nghostpoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Nghostpoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
!
#endif
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
!^          Awrk(i,j,k,Nold)=A(i,j,k)
!^
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO

#ifdef SPLINES_VCONV
!
! <><><> Integrate vertical diffusion equation implicitly using
!        parabolic splines <><><>
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              Hzk(i,k)=0.5_r8*(GRID(ng)%Hz(i,j-1,k)+                    &
     &                         GRID(ng)%Hz(i,j  ,k))
            END DO
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
!^          DC(i,0)=0.0_r8
!^
            tl_DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
!^            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-                          &
!^   &                     Awrk(i,j,k  ,Nold)-                          &
!^   &                     FC(i,k)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Istr,Iend
!^          DC(i,N(ng))=0.0_r8
!^
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!^            DC(i,k)=DC(i,k)*Kv(i,j,k)
!^
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!^   &                         DTsizeV*oHz(i,j,k)*                      &
!^   &                         (DC(i,k)-DC(i,k-1))
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

#else

!
! <><><> Integrate vertical diffusion equation implicitly <><><>
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tangent linear diffusion equation into tl_DC.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=0.5_r8*(GRID(ng)%Hz(i,j-1,k)+GRID(ng)%Hz(i,j,k))
              BC(i,k)=cff-FC(i,j,k)-FC(i,j,k-1)
!^            DC(i,k)=Awrk(i,j,k,Nold)*cff
!^
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*cff
            END DO
          END DO
!
!  Solve the tridiagonal system.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!^          DC(i,1)=cff*DC(i,1)
!^
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!^            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
!^          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!^   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!^   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
# ifdef MASKING
!^          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*                  &
!^   &                           GRID(ng)%vmask(i,j)
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*            &
     &                              GRID(ng)%vmask(i,j)
# endif
          END DO
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^            Awrk(i,j,k,Nnew)=DC(i,k)
!^
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
# ifdef MASKING
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*                        &
!^   &                         GRID(ng)%vmask(i,j)
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*                  &
     &                            GRID(ng)%vmask(i,j)
# endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#endif
!
!  Load convolved solution.
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          A(i,j,k)=Awrk(i,j,k,Nold)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_v3d_tl
!
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix.
!  Note that error correlations are considered separable in the
!  horizontal and vertical directions.
! 
      SUBROUTINE multiscale_Vdiff_v3d_ad (self, ng, tile, model,        &
     &                                    ifield, NVsteps,              &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV, Kv, ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Istr, Iend, JstrV, Jend
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: adfac, cff
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: ad_Awrk
#if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
      real(r8), dimension(IminS:ImaxS,N(ng))               :: Hzk
#else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: ad_DC
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Nold=1
      Nnew=2
!
      ad_Awrk=0.0_r8
      ad_DC=0.0_r8
      FC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=2.0_r8/(GRID(ng)%Hz(i,j-1,k)+                    &
     &                         GRID(ng)%Hz(i,j  ,k))
          END DO
#else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*(Kv(i,j-1,k)+Kv(i,j,k))/                 &
     &                 (GRID(ng)%z_r(i,j-1,k+1)+GRID(ng)%z_r(i,j,k+1)-  &
     &                  GRID(ng)%z_r(i,j-1,k  )-GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#endif
        END DO
      END DO
!
!  Adjoint of load convolved solution.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#endif
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
!^
            ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

#ifdef SPLINES_VCONV
!
! <><><> Integrate adjoint vertical diffusion equation implicitly
!        using parabolic splines <><><>
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute BASIC STATE time-invariant coefficients.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              Hzk(i,k)=0.5_r8*(GRID(ng)%Hz(i,j-1,k)+                    &
     &                         GRID(ng)%Hz(i,j  ,k))
            END DO
          END DO
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              FC(i,k)=cff1*Hzk(i,k  )-DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*Hzk(i,k+1)-DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
          END DO
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Istr,Iend
              BC(i,k)=cff1*(Hzk(i,k)+Hzk(i,k+1))+                       &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
            END DO
          END DO
!
!  Adjoint of backward substitution.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
!^   &                            DTsizeV*oHz(i,j,k)*                   &
!^   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!^
              adfac=DTsizeV*oHz(i,j,k)*ad_Awrk(i,j,k,Nnew)
              ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
              ad_DC(i,k  )=ad_DC(i,k  )+adfac
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+                  &
     &                            ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^
              ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
            END DO
          END DO
!
          DO k=1,N(ng)-1
            DO i=Istr,Iend
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
          DO i=Istr,Iend
!^          tl_DC(i,N(ng))=0.0_r8
!^
            ad_DC(i,N(ng))=0.0_r8
          END DO
!
!  Adjoint of LU decomposition and forward substitution.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
!^   &                        tl_Awrk(i,j,k  ,Nold)-                    &
!^   &                        FC(i,k)*tl_DC(i,k-1))
              adfac=cff*ad_DC(i,k)
              ad_Awrk(i,j,k  ,Nold)=ad_Awrk(i,j,k  ,Nold)-adfac
              ad_Awrk(i,j,k+1,Nold)=ad_Awrk(i,j,k+1,Nold)+adfac
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
          DO i=Istr,Iend
!^          tl_DC(i,0)=0.0_r8
!^
            ad_DC(i,0)=0.0_r8
          END DO
        END DO
      END DO

#else

!
! <><><> Integrate adjoint vertical diffusion equation implicitly <><><>
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Compute diagonal matrix coefficients BC.
!
        DO j=JstrV,Jend
          DO k=1,N(ng)
            DO i=Istr,Iend
              BC(i,k)=0.5_r8*(GRID(ng)%Hz(i,j-1,k)+                     &
     &                        GRID(ng)%Hz(i,j  ,k))-                    &
     &                FC(i,j,k)-FC(i,j,k-1)
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
            END DO
          END DO
!
!^        DO k=N(ng)-1,1,-1
!^
          DO k=1,N(ng)-1
            DO i=Istr,Iend
# ifdef MASKING
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*                  &
!^   &                            GRID(ng)%vmask(i,j)
!^
              ad_Awrk(i,j,k,Nnew)=ad_Awrk(i,j,k,Nnew)*                  &
     &                            GRID(ng)%vmask(i,j)
# endif
!^            tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
!^
              ad_DC(i,k)=ad_DC(i,k)+                                    &
     &                   ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
!
          DO i=Istr,Iend
# ifdef MASKING
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*            &
!^   &                              GRID(ng)%vmask(i,j)
!^
            ad_Awrk(i,j,N(ng),Nnew)=ad_Awrk(i,j,N(ng),Nnew)*            &
     &                              GRID(ng)%vmask(i,j)
# endif
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
!^
            ad_DC(i,N(ng))=ad_DC(i,N(ng))+                              &
     &                     ad_Awrk(i,j,N(ng),Nnew)
            ad_Awrk(i,j,N(ng),Nnew)=0.0_r8
!^          tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
!^   &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!^   &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            adfac=ad_DC(i,N(ng))/                                       &
     &            (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-FC(i,j,N(ng)-1)*adfac
            ad_DC(i,N(ng)  )=adfac
          END DO
!
!  Adjoint of solve tridiagonal system.
!
!^        DO k=2,N(ng)-1
!^
          DO k=N(ng)-1,2,-1
            DO i=Istr,Iend
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,j,k-1)*adfac
              ad_DC(i,k  )=adfac
            END DO
          END DO
          DO i=Istr,Iend
            cff=1.0_r8/BC(i,1)
!^          tl_DC(i,1)=cff*tl_DC(i,1)
!^
            ad_DC(i,1)=cff*ad_DC(i,1)
          END DO
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=0.5_r8*(GRID(ng)%Hz(i,j-1,k)+GRID(ng)%Hz(i,j,k))
!^            tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*cff
!^
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+cff*ad_DC(i,k)
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
        END DO
      END DO
#endif
!
!  Adjoint of set operator initial conditions.
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k,Nold)
            ad_Awrk(i,j,k,Nold)=0.0_r8
          END DO
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    Nghostpoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       Nghostpoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_v3d_ad

#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix for
!  adjusting lateral boundaries in the control vector.
!
      SUBROUTINE multiscale_bry_Vdiff_tl (self, ng, tile, model,        &
     &                                    ifield, ibry, ctype, NVsteps, &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    LBij, UBij,                   &
     &                                    DTsizeV, Kv, tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: LBij, UBij
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
      integer                           :: bindex
!
      real (r8)                         :: cff, cff1
!
      real (r8), pointer                :: mask(:,:) => NULL()
!
      real(r8), dimension(LBij:UBij,N(ng),2) :: tl_Awrk
# if defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng))   :: oHz
# endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: tl_DC
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_bry_Vdiff_tl"
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
        CASE (r3dvar)
          bindex=r2dvar
# ifdef MASKING
          mask => GRID(ng)%rmask
# endif
        CASE (u3dvar)
          bindex=u2dvar
          Imin=IstrU
# ifdef MASKING
          mask => GRID(ng)%umask
# endif
        CASE (v3dvar)
          bindex=v2dvar
          Jmin=JstrV
# ifdef MASKING
          mask => GRID(ng)%vmask
# endif
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      Nold=1
      Nnew=2
!
      FC=0.0_r8
      tl_Awrk=0.0_r8
      tl_DC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,bindex)
            DO j=Jmin,Jmax
# ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
# else
              FC(j,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(j,0)=0.0_r8
# endif
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,bindex)
            DO i=Imin,Imax
# ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(i,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
# else
              FC(i,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(i,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(i,0)=0.0_r8
# endif
            END DO
        END SELECT
      END IF
!
!  Set initial conditions.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              Awrk(j,k,Nold)=A(j,k)
!^
                tl_Awrk(j,k,Nold)=tl_A(j,k)
              END DO
            END DO
          CASE (isouth, inorth)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              Awrk(i,k,Nold)=A(i,k)
!^
                tl_Awrk(i,k,Nold)=tl_A(i,k)
              END DO
            END DO
        END SELECT
      END IF

# ifdef SPLINES_VCONV
!
!  Integrate vertical diffusion equation implicitly using parabolic
!  splines.
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/6.0_r8
              DO j=Jmin,Jmax
                DO k=1,N(ng)-1
                  FC(j,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                  CF(j,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
                END DO
                CF(j,0)=0.0_r8
!^              DC(j,0)=0.0_r8
!^
                tl_DC(j,0)=0.0_r8
              END DO
            CASE (isouth, inorth)
              i=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/6.0_r8
              DO i=Imin,Imax
                DO k=1,N(ng)-1
                  FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                  CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
                END DO
                CF(i,0)=0.0_r8
!^              DC(i,0)=0.0_r8
!^
                tl_DC(i,0)=0.0_r8
              END DO
          END SELECT
!
!  LU decomposition and forward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
                  BC(j,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(j,k)+oHz(j,k+1))
                  cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                  CF(j,k)=cff*CF(j,k)
!^                DC(j,k)=cff*(Awrk(j,k+1,Nold)-                        &
!^   &                         Awrk(j,k  ,Nold)-                        &
!^   &                         FC(j,k)*DC(j,k-1))
!^
                  tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                  &
     &                            tl_Awrk(j,k  ,Nold)-                  &
     &                            FC(j,k)*tl_DC(j,k-1))
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO i=Imin,Imax
                  BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(i,k)+oHz(i,k+1))
                  cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                  CF(i,k)=cff*CF(i,k)
!^                DC(i,k)=cff*(Awrk(i,k+1,Nold)-                        &
!^   &                         Awrk(i,k  ,Nold)-                        &
!^   &                         FC(i,k)*DC(i,k-1))
!^
                  tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                  &
     &                            tl_Awrk(i,k  ,Nold)-                  &
     &                            FC(i,k)*tl_DC(i,k-1))
                END DO
              END DO
          END SELECT
!
!  Backward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO j=Jmin,Jmax
!^              DC(j,N(ng))=0.0_r8
!^
                tl_DC(j,N(ng))=0.0_r8
              END DO
              DO k=N(ng)-1,1,-1
                DO j=Jmin,Jmax
!^                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
!^
                  tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
                END DO
              END DO
              DO k=1,N(ng)
                DO j=Jmin,Jmax
!^                DC(j,k)=DC(j,k)*Kv(i,j,k)
!^
                  tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!^                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                        &
!^   &                           DTsizeV*oHz(j,k)*                      &
!^   &                           (DC(j,k)-DC(j,k-1))
!^
                  tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                  &
     &                              DTsizeV*oHz(j,k)*                   &
     &                              (tl_DC(j,k)-tl_DC(j,k-1))
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO i=Imin,Imax
!^              DC(i,N(ng))=0.0_r8
!^
                tl_DC(i,N(ng))=0.0_r8
              END DO
              DO k=N(ng)-1,1,-1
                DO i=Imin,Imax
!^                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
                  tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
                END DO
              END DO
              DO k=1,N(ng)
                DO i=Imin,Imax
!^                DC(i,k)=DC(i,k)*Kv(i,j,k)
!^
                  tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                        &
!^   &                           DTsizeV*oHz(i,k)*                      &
!^   &                           (DC(i,k)-DC(i,k-1))
!^
                  tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                  &
     &                              DTsizeV*oHz(i,k)*                   &
     &                              (tl_DC(i,k)-tl_DC(i,k-1))
                END DO
              END DO
          END SELECT
        END IF
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

# else
!
!  Integrate vertical diffusion equation implicitly.
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the diffusion equation into DC.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  BC(j,k)=GRID(ng)%Hz(i,j,k)-FC(j,k)-FC(j,k-1)
!^                DC(j,k)=Awrk(j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  tl_DC(j,k)=tl_Awrk(j,k,Nold)*GRID(ng)%Hz(i,j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO i=Imin,Imax
                  BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,k)-FC(i,k-1)
!^                DC(i,k)=Awrk(i,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  tl_DC(i,k)=tl_Awrk(i,k,Nold)*GRID(ng)%Hz(i,j,k)
                END DO
              END DO
          END SELECT
!
!  Solve the tridiagonal system.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO j=Jmin,Jmax
                cff=1.0_r8/BC(j,1)
                CF(j,1)=cff*FC(j,1)
!^              DC(j,1)=cff*DC(j,1)
!^
                tl_DC(j,1)=cff*tl_DC(j,1)
              END DO
              DO k=2,N(ng)-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                  CF(j,k)=cff*FC(j,k)
!^                DC(j,k)=cff*(DC(j,k)-FC(j,k-1)*DC(j,k-1))
!^
                  tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO i=Imin,Imax
                cff=1.0_r8/BC(i,1)
                CF(i,1)=cff*FC(i,1)
!^              DC(i,1)=cff*DC(i,1)
!^
                tl_DC(i,1)=cff*tl_DC(i,1)
              END DO
              DO k=2,N(ng)-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                  CF(i,k)=cff*FC(i,k)
!^                DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
!^
                  tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
                END DO
              END DO
          END SELECT
!
!  Compute new solution by back substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO j=Jmin,Jmax
!^              DC(j,N(ng))=(DC(j,N(ng))-                               &
!^   &                       FC(j,N(ng)-1)*DC(j,N(ng)-1))/              &
!^   &                      (BC(j,N(ng))-                               &
!^   &                       FC(j,N(ng)-1)*CF(j,N(ng)-1))
!^
                tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                         &
     &                          FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/        &
     &                         (BC(j,N(ng))-                            &
     &                          FC(j,N(ng)-1)*CF(j,N(ng)-1))
!^              Awrk(j,N(ng),Nnew)=DC(j,N(ng))
!^
                tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
#  ifdef MASKING
!^              Awrk(j,N(ng),Nnew)=Awrk(j,N(ng),Nnew)*mask(i,j)
!^
                tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*mask(i,j)
#  endif
              END DO
              DO k=N(ng)-1,1,-1
                DO j=Jmin,Jmax
!^                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
!^
                  tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!^                Awrk(j,k,Nnew)=DC(j,k)
!^
                  tl_Awrk(j,k,Nnew)=tl_DC(j,k)
#  ifdef MASKING
!^                Awrk(j,k,Nnew)=Awrk(j,k,Nnew)*mask(i,j)
!^
                  tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*mask(i,j)
#  endif
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO i=Imin,Imax
!^              DC(i,N(ng))=(DC(i,N(ng))-                               &
!^   &                       FC(i,N(ng)-1)*DC(i,N(ng)-1))/              &
!^   &                      (BC(i,N(ng))-                               &
!^   &                       FC(i,N(ng)-1)*CF(i,N(ng)-1))
!^
                tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                         &
     &                          FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/        &
     &                         (BC(i,N(ng))-                            &
     &                          FC(i,N(ng)-1)*CF(i,N(ng)-1))
!^              Awrk(i,N(ng),Nnew)=DC(i,N(ng))
!^
                tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
#  ifdef MASKING
!^              Awrk(i,N(ng),Nnew)=Awrk(i,N(ng),Nnew)*mask(i,j)
!^
                tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*mask(i,j)
#  endif
              END DO
              DO k=N(ng)-1,1,-1
                DO i=Imin,Imax
!^                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
                  tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^                Awrk(i,k,Nnew)=DC(i,k)
!^
                  tl_Awrk(i,k,Nnew)=tl_DC(i,k)
#  ifdef MASKING
!^                Awrk(i,k,Nnew)=Awrk(i,k,Nnew)*mask(i,j)
!^
                  tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*mask(i,j)
#  endif
                END DO
              END DO
          END SELECT
        END IF
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
# endif
!
!  Load convolved solution.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,bindex)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              A(j,k)=Awrk(j,k,Nold)
!^
                tl_A(j,k)=tl_Awrk(j,k,Nold)
              END DO
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,bindex)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              A(i,k)=Awrk(i,k,Nold)
!^
                tl_A(i,k)=tl_Awrk(i,k,Nold)
              END DO
            END DO
        END SELECT
      END IF

# ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        A)
!^
      CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        tl_A)
# endif
!
      RETURN
      END SUBROUTINE multiscale_bry_Vdiff_tl
!
!-----------------------------------------------------------------------
!  It advances the adjoint vertical pseudo-difusion operator using an
!  implicit solver using parabolic splines or inverting a tridiagonal
!  matrix for adjusting lateral boundaries in the control vector.
!
      SUBROUTINE multiscale_bry_Vdiff_ad (self, ng, tile, model,        &
     &                                    ifield, ibry, ctype, NVsteps, &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    LBij, UBij,                   &
     &                                    DTsizeV, Kv, ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: LBij, UBij
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: ad_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
      integer                           :: bindex
!
      real (r8)                         :: adfac, cff, cff1
!
      real (r8), pointer                :: mask(:,:) => NULL()
!
      real(r8), dimension(LBij:UBij,N(ng),2) :: ad_Awrk
# if defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng))   :: oHz
# endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_DC
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_bry_Vdiff_ad"
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
        CASE (r3dvar)
          bindex=r2dvar
# ifdef MASKING
          mask => GRID(ng)%rmask
# endif
        CASE (u3dvar)
          bindex=u2dvar
          Imin=IstrU
# ifdef MASKING
          mask => GRID(ng)%umask
# endif
        CASE (v3dvar)
          bindex=v2dvar
          Jmin=JstrV
# ifdef MASKING
          mask => GRID(ng)%vmask
# endif
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      Nold=1
      Nnew=2
!
      ad_Awrk=0.0_r8
      ad_DC=0.0_r8
      FC=0.0_r8
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,bindex)
            DO j=Jmin,Jmax
# ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
# else
              FC(j,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(j,0)=0.0_r8
# endif
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,bindex)
            DO i=Imin,Imax
# ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(i,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
# else
              FC(i,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(i,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(i,0)=0.0_r8
# endif
            END DO
        END SELECT
      END IF
!
!  Adjoint of load convolved solution.
!
# ifdef DISTRIBUTE
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_A
!^
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif

      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,bindex)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              tl_A(j,k)=tl_Awrk(j,k,Nold)
!^
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_A(j,k)
                ad_A(j,k)=0.0_r8
              END DO
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,bindex)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              tl_A(i,k)=tl_Awrk(i,k,Nold)
!^
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_A(i,k)
                ad_A(i,k)=0.0_r8
              END DO
            END DO
        END SELECT
      END IF

# ifdef SPLINES_VCONV
!
!  Integrate adjoint vertical diffusion equation implicitly using
!  parabolic splines.
!
      STEP_LOOP : DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute basic state time-invariant coefficients.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/6.0_r8
              DO j=Jmin,Jmax
                DO k=1,N(ng)-1
                  FC(j,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                  CF(j,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
                END DO
                CF(j,0)=0.0_r8
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/6.0_r8
              DO i=Imin,Imax
                DO k=1,N(ng)-1
                  FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                  CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
                END DO
                CF(i,0)=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of LU decomposition and forward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
                  BC(j,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(j,k)+oHz(j,k+1))
                  cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                  CF(j,k)=cff*CF(j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO i=Imin,Imax
                  BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(i,k)+oHz(i,k+1))
                  cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                  CF(i,k)=cff*CF(i,k)
                END DO
              END DO
          END SELECT
!
!  Adjoint backward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
!^                tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                  &
!^   &                              DTsizeV*oHz(j,k)*                   &
!^   &                              (tl_DC(j,k)-tl_DC(j,k-1))
!^
                  adfac=DTsizeV*oHz(j,k)*ad_Awrk(j,k,Nnew)
                  ad_DC(j,k-1)=ad_DC(j,k-1)-adfac
                  ad_DC(j,k  )=ad_DC(j,k  )+adfac
                  ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                  &
     &                              ad_Awrk(j,k,Nnew)
                  ad_Awrk(j,k,Nnew)=0.0_r8
!^                tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!^
                  ad_DC(j,k)=ad_DC(j,k)*Kv(i,j,k)
                END DO
              END DO
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
!^                tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!^
                  ad_DC(j,k+1)=ad_DC(j,k+1)-CF(j,k)*ad_DC(j,k)
                END DO
              END DO
              DO j=Jmin,Jmax
!^              tl_DC(j,N(ng))=0.0_r8
!^
                ad_DC(j,N(ng))=0.0_r8
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO i=Imin,Imax
!^                tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                  &
!^   &                              DTsizeV*oHz(i,k)*                   &
!^   &                              (tl_DC(i,k)-tl_DC(i,k-1))
!^
                  adfac=DTsizeV*oHz(i,k)*ad_Awrk(i,k,Nnew)
                  ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
                  ad_DC(i,k  )=ad_DC(i,k  )+adfac
                  ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                  &
     &                              ad_Awrk(i,k,Nnew)
                  ad_Awrk(i,k,Nnew)=0.0_r8
!^                tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^
                  ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
                END DO
              END DO
              DO k=1,N(ng)-1
                DO i=Imin,Imax
!^                tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
                  ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
                END DO
              END DO
              DO i=Imin,Imax
!^              tl_DC(i,N(ng))=0.0_r8
!^
                ad_DC(i,N(ng))=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of LU decomposition and forward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO k=N(ng)-1,1,-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
!^                tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                  &
!^   &                            tl_Awrk(j,k  ,Nold)-                  &
!^   &                            FC(j,k)*tl_DC(j,k-1))
!^
                  adfac=cff*ad_DC(j,k)
                  ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                  ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                  ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k)*adfac
                  ad_DC(j,k)=0.0_r8
                END DO
              END DO
              DO j=Jmin,Jmax
!^              tl_DC(j,0)=0.0_r8
!^
                ad_DC(j,0)=0.0_r8
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO k=N(ng)-1,1,-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!^                tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                  &
!^   &                            tl_Awrk(i,k  ,Nold)-                  &
!^   &                            FC(i,k)*tl_DC(i,k-1))
!^
                  adfac=cff*ad_DC(i,k)
                  ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                  ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                  ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
                  ad_DC(i,k)=0.0_r8
                END DO
              END DO
              DO i=Imin,Imax
!^              tl_DC(i,0)=0.0_r8
!^
                ad_DC(i,0)=0.0_r8
              END DO
          END SELECT
        END IF
      END DO STEP_LOOP

# else
!
!  Integrate adjoint vertical diffusion equation implicitly.
!
      STEP_LOOP : DO step=1,NVsteps
!
!  Adjoint of update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Compute diagonal matrix coefficients BC and CF.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  BC(j,k)=GRID(ng)%Hz(i,j,k)-FC(j,k)-FC(j,k-1)
                END DO
              END DO
              DO j=Jmin,Jmax
                cff=1.0_r8/BC(j,1)
                CF(j,1)=cff*FC(j,1)
              END DO
              DO k=2,N(ng)-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                  CF(j,k)=cff*FC(j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO i=Imin,Imax
                  BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,k)-FC(i,k-1)
                END DO
              END DO
              DO i=Imin,Imax
                cff=1.0_r8/BC(i,1)
                CF(i,1)=cff*FC(i,1)
              END DO
              DO k=2,N(ng)-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                  CF(i,k)=cff*FC(i,k)
                END DO
              END DO
          END SELECT
!
!  Adjoint of compute new solution by back substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
!^            DO k=N(ng)-1,1,-1
!^
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
#  ifdef MASKING
!^                tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*mask(i,j)
!^
                  ad_Awrk(j,k,Nnew)=ad_Awrk(j,k,Nnew)*mask(i,j)
#  endif
!^                tl_Awrk(j,k,Nnew)=tl_DC(j,k)
!^
                  ad_DC(j,k)=ad_DC(j,k)+                                &
     &                       ad_Awrk(j,k,Nnew)
                  ad_Awrk(j,k,Nnew)=0.0_r8
!^                tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!^
                  ad_DC(j,k+1)=-CF(j,k)*ad_DC(j,k)
                END DO
              END DO
              DO j=Jmin,Jmax
#  ifdef MASKING
!^              tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*mask(i,j)
!^
                ad_Awrk(j,N(ng),Nnew)=ad_Awrk(j,N(ng),Nnew)*mask(i,j)

#  endif
!^              tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
!^
                ad_DC(j,N(ng))=ad_DC(j,N(ng))+                          &
     &                         ad_Awrk(j,N(ng),Nnew)
                ad_Awrk(j,N(ng),Nnew)=0.0_r8
!^              tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                         &
!^   &                          FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/        &
!^   &                         (BC(j,N(ng))-                            &
!^   &                          FC(j,N(ng)-1)*CF(j,N(ng)-1))
!^
                adfac=ad_DC(j,N(ng))/                                   &
     &                (BC(j,N(ng))-                                     &
     &                 FC(j,N(ng)-1)*CF(j,N(ng)-1))
                ad_DC(j,N(ng)-1)=ad_DC(j,N(ng)-1)-                      &
     &                           FC(j,N(ng)-1)*adfac
                ad_DC(j,N(ng)  )=adfac
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
!^            DO k=N(ng)-1,1,-1
!^
              DO k=1,N(ng)-1
                DO i=Imin,Imax
#  ifdef MASKING
!^                tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*mask(i,j)
!^
                  ad_Awrk(i,k,Nnew)=ad_Awrk(i,k,Nnew)*mask(i,j)
#  endif
!^                tl_Awrk(i,k,Nnew)=tl_DC(i,k)
!^
                  ad_DC(i,k)=ad_DC(i,k)+                                &
     &                       ad_Awrk(i,k,Nnew)
                  ad_Awrk(i,k,Nnew)=0.0_r8
!^                tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
                  ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
                END DO
              END DO
              DO i=Imin,Imax
#  ifdef MASKING
!^              tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*mask(i,j)
!^
                ad_Awrk(i,N(ng),Nnew)=ad_Awrk(i,N(ng),Nnew)*mask(i,j)
#  endif
!^              tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
!^
                ad_DC(i,N(ng))=ad_DC(i,N(ng))+                          &
     &                         ad_Awrk(i,N(ng),Nnew)
                ad_Awrk(i,N(ng),Nnew)=0.0_r8
!^              tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                         &
!^   &                          FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/        &
!^   &                         (BC(i,N(ng))-                            &
!^   &                          FC(i,N(ng)-1)*CF(i,N(ng)-1))
!^
                adfac=ad_DC(i,N(ng))/                                   &
     &                (BC(i,N(ng))-                                     &
     &                 FC(i,N(ng)-1)*CF(i,N(ng)-1))
                ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-                      &
     &                           FC(i,N(ng)-1)*adfac
                ad_DC(i,N(ng)  )=adfac
              END DO
            END SELECT
!
!  Solve the adjoint tridiagonal system.
!
           SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
!^            DO k=2,N(ng)-1
!^
              DO k=N(ng)-1,2,-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
!^                tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
!^
                  adfac=cff*ad_DC(j,k)
                  ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k-1)*adfac
                  ad_DC(j,k  )=adfac
                END DO
              END DO
              DO j=Jmin,Jmax
                cff=1.0_r8/BC(j,1)
!^              tl_DC(j,1)=cff*tl_DC(j,1)
!^
                ad_DC(j,1)=cff*ad_DC(j,1)
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
!^            DO k=2,N(ng)-1
!^
              DO k=N(ng)-1,2,-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!^                tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
!^
                  adfac=cff*ad_DC(i,k)
                  ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
                  ad_DC(i,k  )=adfac
                END DO
              END DO
              DO i=Imin,Imax
                cff=1.0_r8/BC(i,1)
!^              tl_DC(i,1)=cff*tl_DC(i,1)
!^
                ad_DC(i,1)=cff*ad_DC(i,1)
              END DO
          END SELECT
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
!^                tl_DC(j,k)=tl_Awrk(j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                  &
     &                              GRID(ng)%Hz(i,j,k)*ad_DC(j,k)
                  ad_DC(j,k)=0.0_r8
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,bindex)
              DO k=1,N(ng)
                DO i=Imin,Imax
!^                tl_DC(i,k)=tl_Awrk(i,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                  &
     &                              GRID(ng)%Hz(i,j,k)*ad_DC(i,k)
                  ad_DC(i,k)=0.0_r8
                END DO
              END DO
          END SELECT
        END IF
      END DO STEP_LOOP
# endif
!
!  Set adjoint initial conditions.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,bindex)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              tl_Awrk(j,k,Nold)=tl_A(j,k)
!^
                ad_A(j,k)=ad_A(j,k)+ad_Awrk(j,k,Nold)
                ad_Awrk(j,k,Nold)=0.0_r8
              END DO
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,bindex)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              tl_Awrk(i,k,Nold)=tl_A(i,k)
!^
                ad_A(i,k)=ad_A(i,k)+ad_Awrk(i,k,Nold)
                ad_Awrk(i,k,Nold)=0.0_r8
              END DO
            END DO
        END SELECT
      END IF

# ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_A)
!^
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!
      RETURN
      END SUBROUTINE multiscale_bry_Vdiff_ad
#endif /* ADJUST_BOUNDARY */
