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
** These routines computes the K-Laplacian operator, [1 + Del(K*Del)],**
** which is subsequently used to determine the minimum and maximum    **
** eigenvalues required by the implicit diffusion iterative algorithm **
** for modeling multiscale background-error covariance matrix effects **
** (Weaver et al., 2016).                                             **
**                                                                    **
** The horizontal diffusion coefficients (Khx, Khy) are expressed in  **
** units of correlation length squared. These coefficients are        **
** derived from Equation (44) of Weaver and Mirouze (2013) for d=2.   **
** Thus, kappa=D*D/(2*(M-2)). Here, D denotes the Daley length scale, **
** and d represents the spatial dimension.                            **
**                                                                    **
** References:                                                        **
**                                                                    **
** Weaver, A.T. and I. Mirouze, 2013: On the diffusion equation and   **
**   its application to isotropic and anisotropic correlation         **
**   modeling in variational assimilation, Q. J. R. Meteorol. Soc.,   **
**   139, 242-260, doi:10.1002/qj.1955.                               **
**                                                                    **
** Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation   **
**   operators based on an implicitly formulated diffusion equation   **
**   solved with the Chebyshev iteration, Q. J. R. Meteorol. Soc.,    **
**   142, 455-471, doi:10.1002/qj.2664.                               **
**                                                                    **
************************************************************************
**
**
**  <><> CLASS MULTISCALE:  K-Laplacian Operators <><><><><><><><><><><>
**
*/

!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of a
!  2D control variable at RHO-points.
!
      SUBROUTINE multiscale_Klap_r2d_tl (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Mlap, i, j, rec
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
      integer                           :: itrc
!
      real (r8)                         :: cffx, cffy

#ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: tl_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
!
!  Initialize.
!
      Istr=BOUNDS(ng)%Istr(tile)     ! tile computational indices
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')
#ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%zeta_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%zeta_Bcorr(is:ie,js:je,2,ms)
#endif
          Mlap=self%Mlap(ifield,ms)
#if defined ADJUST_STFLUX && defined SOLVE3D
        CASE ('shflux', 'ssflux')
          itrc = tracer_index(TRIM(StateVarName(ifield)))
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%stflux_Bcorr(is:ie,js:je,1,ms,itrc)
          BscaleY(is:,js:) => self%stflux_Bcorr(is:ie,js:je,2,ms,itrc)
# endif
          Mlap=self%Mlap(ifield,ms)
#endif
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
#else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
#endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Set operator initial conditions.
!
!^    CALL dabc_r2d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    A)
!^
      CALL dabc_r2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_A)
#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!^        Awrk(i,j)=A(i,j)
!^
          tl_Awrk(i,j)=tl_A(i,j)
        END DO
      END DO
!
!  Compute XI- and ETA-components of diffusive flux.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend+1
!^        FX(i,j)=GRID(ng)%pmon_u(i,j)*                                 &
!^   &            0.5_r8*(Khx(i-1,j)+Khx(i,j))*                         &
!^   &            (Awrk(i,j)-Awrk(i-1,j))
!^
          tl_FX(i,j)=GRID(ng)%pmon_u(i,j)*                              &
     &               0.5_r8*(Khx(i-1,j)+Khx(i,j))*                      &
     &               (tl_Awrk(i,j)-tl_Awrk(i-1,j))
#ifdef MASKING
!^        FX(i,j)=FX(i,j)*GRID(ng)%umask(i,j)
!^
          tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%umask(i,j)
#endif
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
!^        FE(i,j)=GRID(ng)%pnom_v(i,j)*                                 &
!^                0.5_r8*(Khy(i,j-1)+Khy(i,j))*                         &
!^   &              (Awrk(i,j)-Awrk(i,j-1))
!^
          tl_FE(i,j)=GRID(ng)%pnom_v(i,j)*                              &
     &               0.5_r8*(Khy(i,j-1)+Khy(i,j))*                      &
     &               (tl_Awrk(i,j)-tl_Awrk(i,j-1))
#ifdef MASKING
!^        FE(i,j)=FE(i,j)*GRID(ng)%vmask(i,j)
!^
          tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%vmask(i,j)
#endif
        END DO
      END DO
!
!  Compute K-laplacian.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!^        Awrk(i,j)=Awrk(i,j)-                                          &
!^   &              Hfac(i,j)*(FX(i+1,j)-FX(i,j)+                       &
!^   &                         FE(i,j+1)-FE(i,j))
!^
          tl_Awrk(i,j)=tl_Awrk(i,j)-                                    &
     &                 Hfac(i,j)*(tl_FX(i+1,j)-tl_FX(i,j)+              &
     &                            tl_FE(i,j+1)-tl_FE(i,j))
        END DO
      END DO
!
!  Load K-Laplacian solution.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!^        A(i,j)=Awrk(i,j)
!^
          tl_A(i,j)=tl_Awrk(i,j)
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif

#ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_r2d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of a
!  2D control variable at RHO-points.
!
      SUBROUTINE multiscale_Klap_r2d_ad (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:)
!
      integer                           :: Mlap, i, j, rec
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
      integer                           :: itrc
!
      real (r8)                         :: adfac, cffx, cffy

#ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: ad_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
!
!  Initialize.
!
      Istr=BOUNDS(ng)%Istr(tile)     ! tile computational indices
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')
#ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%zeta_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%zeta_Bcorr(is:ie,js:je,2,ms)
#endif
          Mlap=self%Mlap(ifield,ms)
#if defined ADJUST_STFLUX && defined SOLVE3D
        CASE ('shflux', 'ssflux')
          itrc = tracer_index(TRIM(StateVarName(ifield)))
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%stflux_Bcorr(is:ie,js:je,1,ms,itrc)
          BscaleY(is:,js:) => self%stflux_Bcorr(is:ie,js:je,2,ms,itrc)
# endif
          Mlap=self%Mlap(ifield,ms)
#endif
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
#else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
#endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Adjoint of load K-Laplacian solution.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#endif
      DO j=Jstr,Jend
        DO i=Istr,Iend
!^        tl_A(i,j)=tl_Awrk(i,j)
!^
          ad_Awrk(i,j)=ad_Awrk(i,j)+ad_A(i,j)
          ad_A(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of compute K-Laplacian.
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
!^        tl_Awrk(i,j)=tl_Awrk(i,j)-                                    &
!^   &                 Hfac(i,j)*(tl_FX(i+1,j)-tl_FX(i,j)+              &
!^   &                            tl_FE(i,j+1)-tl_FE(i,j))
!^
          adfac=-Hfac(i,j)*ad_Awrk(i,j)
          ad_FE(i,j  )=ad_FE(i,j  )-adfac
          ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
          ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
          ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
        END DO
      END DO
!
!  Adjoint of compute XI- and ETA-components of diffusive flux.
!
      DO j=Jstr,Jend+1
        DO i=Istr,Iend
#ifdef MASKING
!^        tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%vmask(i,j)
!^
          ad_FE(i,j)=ad_FE(i,j)*GRID(ng)%vmask(i,j)
#endif
!^        tl_FE(i,j)=GRID(ng)%pnom_v(i,j)*                              &
!^   &               0.5_r8*(Khy(i,j-1)+Khy(i,j))*                      &
!^   &               (tl_Awrk(i,j)-tl_Awrk(i,j-1))
!^
          adfac=GRID(ng)%pnom_v(i,j)*                                   &
     &          0.5_r8*(Khy(i,j-1)+Khy(i,j))*ad_FE(i,j)
          ad_Awrk(i,j-1)=ad_Awrk(i,j-1)-adfac
          ad_Awrk(i,j  )=ad_Awrk(i,j  )+adfac
          ad_FE(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=Istr,Iend+1
#ifdef MASKING
!^        tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%umask(i,j)
!^
          ad_FX(i,j)=ad_FX(i,j)*GRID(ng)%umask(i,j)
#endif
!^        tl_FX(i,j)=GRID(ng)%pmon_u(i,j)*                              &
!^   &               0.5_r8*(Khx(i-1,j)+Khx(i,j))*                      &
!^   &               (tl_Awrk(i,j)-tl_Awrk(i-1,j))
!^
          adfac=GRID(ng)%pmon_u(i,j)*                                   &
     &          0.5_r8*(Khx(i-1,j)+Khx(i,j))*ad_FX(i,j)
          ad_Awrk(i-1,j)=ad_Awrk(i-1,j)-adfac
          ad_Awrk(i  ,j)=ad_Awrk(i  ,j)+adfac
          ad_FX(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of set operator initial conditions.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
!^        tl_Awrk(i,j)=tl_A(i,j)
!^
          ad_A(i,j)=ad_A(i,j)+ad_Awrk(i,j)
          ad_Awrk(i,j)=0.0_r8
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^                        EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#endif
!^    CALL dabc_r2d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    tl_A)
!^
      CALL ad_dabc_r2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)

#ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_r2d_ad
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of a
!  2D control variable at U-points.
!
      SUBROUTINE multiscale_Klap_u2d_tl (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Mlap, i, j, rec
      integer                           :: IstrU, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: cffx, cffy

#ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: tl_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
!
!  Initialize.
!
      IstrU=BOUNDS(ng)%IstrU(tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('ubar', 'ubar_eastward')
#ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%ubar_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%ubar_Bcorr(is:ie,js:je,2,ms)
#endif
          Mlap=self%Mlap(ifield,ms)
#ifdef ADJUST_WSTRESS
        CASE ('sustr')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%sustr_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%sustr_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
#endif
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
#ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
#else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
#endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Set operator initial conditions.
!
!^    CALL dabc_u2d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    A)
!^
      CALL dabc_u2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_A)
#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
!^        Awrk(i,j)=A(i,j)
!^
          tl_Awrk(i,j)=tl_A(i,j)
        END DO
      END DO
!
!  Compute XI- and ETA-components of diffusive flux.
!
      DO j=Jstr,Jend
        DO i=IstrU-1,Iend
!^        FX(i,j)=GRID(ng)%pmon_r(i,j)*                                 &
!^   &            Khx(i,j)*(Awrk(i+1,j)-Awrk(i,j))
!^
          tl_FX(i,j)=GRID(ng)%pmon_r(i,j)*                              &
     &               Khx(i,j)*(tl_Awrk(i+1,j)-tl_Awrk(i,j))
        END DO
      END DO
!
      DO j=Jstr,Jend+1
        DO i=IstrU,Iend
!^        FE(i,j)=GRID(ng)%pnom_p(i,j)*                                 &
!^   &            0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                     &
!^   &                     Khy(i-1,j-1)+Khy(i,j-1))*                    &
!^   &            (Awrk(i,j)-Awrk(i,j-1))
!^
          tl_FE(i,j)=GRID(ng)%pnom_p(i,j)*                              &
     &               0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                  &
     &                        Khy(i-1,j-1)+Khy(i,j-1))*                 &
     &               (tl_Awrk(i,j)-tl_Awrk(i,j-1))
#ifdef MASKING
!^        FE(i,j)=FE(i,j)*GRID(ng)%pmask(i,j)
!^
          tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%pmask(i,j)
#endif
        END DO
      END DO
!
!  Compute K-Laplacian operator.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!^        Awrk(i,j)=Awrk(i,j)-                                          &
!^   &              Hfac(i,j)*(FX(i,j)-FX(i-1,j)+                       &
!^   &                         FE(i,j+1)-FE(i,j))
!^
          tl_Awrk(i,j)=tl_Awrk(i,j)-                                    &
     &                 Hfac(i,j)*(tl_FX(i,j)-tl_FX(i-1,j)+              &
     &                            tl_FE(i,j+1)-tl_FE(i,j))
        END DO
      END DO
!
!  Load K-Laplacian solution.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!^        A(i,j)=Awrk(i,j)
!^
          tl_A(i,j)=tl_Awrk(i,j)
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif

#ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_u2d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of a
!  2D control variable at U-points.
!
      SUBROUTINE multiscale_Klap_u2d_ad (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:)
!
      integer                           :: Mlap, i, j, rec
      integer                           :: IstrU, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: adfac, cffx, cffy

#ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: ad_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
!
!  Initialize.
!
      IstrU=BOUNDS(ng)%IstrU(tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('ubar', 'ubar_eastward')
#ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%ubar_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%ubar_Bcorr(is:ie,js:je,2,ms)
#endif
          Mlap=self%Mlap(ifield,ms)
#ifdef ADJUST_WSTRESS
        CASE ('sustr')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%sustr_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%sustr_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
#endif
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
#ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
#else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
#endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Adjoint of load K-Laplacian solution.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_Awrk)
!
#endif
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!^        tl_A(i,j)=tl_Awrk(i,j)
!^
          ad_Awrk(i,j)=ad_Awrk(i,j)+ad_A(i,j)
          ad_A(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of compute K-Laplacian.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
!^        tl_Awrk(i,j)=tl_Awrk(i,j)-                                    &
!^   &                 Hfac(i,j)*(tl_FX(i,j)-tl_FX(i-1,j)+              &
!^   &                            tl_FE(i,j+1)-tl_FE(i,j))
!^
          adfac=-Hfac(i,j)*ad_Awrk(i,j)
          ad_FE(i,j  )=ad_FE(i,j  )-adfac
          ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
          ad_FX(i-1,j)=ad_FX(i-1,j)-adfac
          ad_FX(i  ,j)=ad_FX(i  ,j)+adfac
        END DO
      END DO
!
!  Adjoint of compute XI- and ETA-components of diffusive flux.
!
      DO j=Jstr,Jend+1
        DO i=IstrU,Iend
#ifdef MASKING
!^        tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%pmask(i,j)
!^
          ad_FE(i,j)=ad_FE(i,j)*GRID(ng)%pmask(i,j)
#endif
!^        tl_FE(i,j)=GRID(ng)%pnom_p(i,j)*                              &
!^   &               0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                  &
!^   &                        Khy(i-1,j-1)+Khy(i,j-1))*                 &
!^   &               (tl_Awrk(i,j)-tl_Awrk(i,j-1))
!^
          adfac=GRID(ng)%pnom_p(i,j)*                                   &
     &          0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                       &
     &                   Khy(i-1,j-1)+Khy(i,j-1))*ad_FE(i,j)
          ad_Awrk(i,j-1)=ad_Awrk(i,j-1)-adfac
          ad_Awrk(i,j  )=ad_Awrk(i,j  )+adfac
          ad_FE(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=Jstr,Jend
        DO i=IstrU-1,Iend
!^        tl_FX(i,j)=GRID(ng)%pmon_r(i,j)*                              &
!^   &               Khx(i,j)*(tl_Awrk(i+1,j)-tl_Awrk(i,j))
!^
          adfac=GRID(ng)%pmon_r(i,j)*Khx(i,j)*ad_FX(i,j)
          ad_Awrk(i  ,j)=ad_Awrk(i  ,j)-adfac
          ad_Awrk(i+1,j)=ad_Awrk(i+1,j)+adfac
          ad_FX(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of set operator initial conditions.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
!^        tl_Awrk(i,j)=tl_A(i,j)
!^
          ad_A(i,j)=ad_A(i,j)+ad_Awrk(i,j)
          ad_Awrk(i,j)=0.0_r8
        END DO
      END DO
#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
#endif
!
!^    CALL dabc_u2d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    tl_A)
!^
      CALL ad_dabc_u2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)

#ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_u2d_ad
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of a
!  2D control variable at V-points.
!
      SUBROUTINE multiscale_Klap_v2d_tl (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Mlap, i, j, rec
      integer                           :: Istr, Iend, JstrV, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: cffx, cffy

#ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: tl_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('vbar', 'vbar_northward')
#ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%vbar_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%vbar_Bcorr(is:ie,js:je,2,ms)
#endif
          Mlap=self%Mlap(ifield,ms)
#ifdef ADJUST_WSTRESS
        CASE ('svstr')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%svstr_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%svstr_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
#endif
      END SELECT
!
!  Compute metrics factor.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
#else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
#endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Set operator initial conditions.
!
!^    CALL dabc_v2d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    A)
!^
      CALL dabc_v2d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    tl_A)

#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#endif
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
!^        Awrk(i,j)=A(i,j)
!^
          tl_Awrk(i,j)=tl_A(i,j)
        END DO
      END DO
!
!  Compute XI- and ETA-components of diffusive flux.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend+1
!^        FX(i,j)=GRID(ng)%pmon_p(i,j)*                                 &
!^   &            0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                     &
!^   &                     Khx(i-1,j-1)+Khx(i,j-1))*                    &
!^   &            (Awrk(i,j)-Awrk(i-1,j))
!^
          tl_FX(i,j)=GRID(ng)%pmon_p(i,j)*                              &
     &               0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                  &
     &                        Khx(i-1,j-1)+Khx(i,j-1))*                 &
     &               (tl_Awrk(i,j)-tl_Awrk(i-1,j))
#ifdef MASKING
!^        FX(i,j)=FX(i,j)*GRID(ng)%pmask(i,j)
!^
          tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%pmask(i,j)
#endif
        END DO
      END DO
!
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
!^        FE(i,j)=GRID(ng)%pnom_r(i,j)*                                 &
!^   &            Khy(i,j)*(Awrk(i,j+1)-Awrk(i,j))
!^
          tl_FE(i,j)=GRID(ng)%pnom_r(i,j)*                              &
     &               Khy(i,j)*(tl_Awrk(i,j+1)-tl_Awrk(i,j))
        END DO
      END DO
!
!  Compute K-Laplacian operator.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
!^        Awrk(i,j)=Awrk(i,j)-                                          &
!^   &              Hfac(i,j)*(FX(i+1,j)-FX(i,j)+                       &
!^   &                         FE(i,j)-FE(i,j-1))
!^
          tl_Awrk(i,j)=tl_Awrk(i,j)-                                    &
     &                 Hfac(i,j)*(tl_FX(i+1,j)-tl_FX(i,j)+              &
     &                            tl_FE(i,j)-tl_FE(i,j-1))
        END DO
      END DO
!
!  Load K-Laplacian solution.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
!^        A(i,j)=Awrk(i,j)
!^
          tl_A(i,j)=tl_Awrk(i,j)
        END DO
      END DO

#ifdef DISTRIBUTE
!
!^      CALL mp_exchange2d (ng, tile, model, 1,                         &
!^   &                      LBi, UBi, LBj, UBj,                         &
!^   &                      NghostPoints,                               &
!^   &                      EWperiodic(ng), NSperiodic(ng),             &
!^   &                      A)
!^
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      tl_A)
#endif

#ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_v2d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of a
!  2D control variable at V-points.
!
      SUBROUTINE multiscale_Klap_v2d_ad (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:)
!
      integer                           :: Mlap, i, j, rec
      integer                           :: Istr, Iend, JstrV, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: adfac, cffx, cffy

#ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj)         :: ad_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('vbar', 'vbar_northward')
#ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%vbar_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%vbar_Bcorr(is:ie,js:je,2,ms)
#endif
          Mlap=self%Mlap(ifield,ms)
#ifdef ADJUST_WSTRESS
        CASE ('svstr')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%svstr_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%svstr_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
#endif
      END SELECT
!
!  Compute metrics factor.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
#ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
#else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
#endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Adjoint of load K-Laplacian solution.
!
#ifdef DISTRIBUTE
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#endif
      DO j=JstrV,Jend
        DO i=Istr,Iend
!^        tl_A(i,j)=tl_Awrk(i,j)
!^
          ad_Awrk(i,j)=ad_Awrk(i,j)+ad_A(i,j)
          ad_A(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of compute K-Laplacian operator.
!
      DO j=JstrV,Jend
        DO i=Istr,Iend
!^        tl_Awrk(i,j)=tl_Awrk(i,j)-                                    &
!^   &                 Hfac(i,j)*(tl_FX(i+1,j)-tl_FX(i,j)+              &
!^   &                            tl_FE(i,j)-tl_FE(i,j-1))
!^
          adfac=-Hfac(i,j)*ad_Awrk(i,j)
          ad_FE(i,j-1)=ad_FE(i,j-1)-adfac
          ad_FE(i,j  )=ad_FE(i,j  )+adfac
          ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
          ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
        END DO
      END DO
!
!  Adjoint of compute XI- and ETA-components of diffusive flux.
!
      DO j=JstrV-1,Jend
        DO i=Istr,Iend
!^        tl_FE(i,j)=GRID(ng)%pnom_r(i,j)*                              &
!^   &               Khy(i,j)*(tl_Awrk(i,j+1)-tl_Awrk(i,j))
!^
          adfac=GRID(ng)%pnom_r(i,j)*Khy(i,j)*ad_FE(i,j)
          ad_Awrk(i,j  )=ad_Awrk(i,j  )-adfac
          ad_Awrk(i,j+1)=ad_Awrk(i,j+1)+adfac
          ad_FE(i,j)=0.0_r8
        END DO
      END DO
!
      DO j=JstrV,Jend
        DO i=Istr,Iend+1
#ifdef MASKING
!^        tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%pmask(i,j)
!^
          ad_FX(i,j)=ad_FX(i,j)*GRID(ng)%pmask(i,j)
#endif
!^        tl_FX(i,j)=GRID(ng)%pmon_p(i,j)*                              &
!^   &               0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                  &
!^   &                        Khx(i-1,j-1)+Khx(i,j-1))*                 &
!^   &               (tl_Awrk(i,j)-tl_Awrk(i-1,j))
!^
          adfac=GRID(ng)%pmon_p(i,j)*                                   &
     &          0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                       &
     &                   Khx(i-1,j-1)+Khx(i,j-1))*ad_FX(i,j)
          ad_Awrk(i-1,j)=ad_Awrk(i-1,j)-adfac
          ad_Awrk(i  ,j)=ad_Awrk(i  ,j)+adfac
          ad_FX(i,j)=0.0_r8
        END DO
      END DO
!
!  Adjoint of set operator initial conditions.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
!^        tl_Awrk(i,j)=tl_A(i,j)
!^
          ad_A(i,j)=ad_A(i,j)+ad_Awrk(i,j)
          ad_Awrk(i,j)=0.0_r8
        END DO
      END DO
#ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange2d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
#endif
!
!^    CALL dabc_v2d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj,                           &
!^   &                    tl_A)
!^
      CALL ad_dabc_v2d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ad_A)

#ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
#endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_v2d_ad

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of a
!  3D control variable at RHO-points.
!
      SUBROUTINE multiscale_Klap_r3d_tl (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Mlap, i, j, k, k1, k2, rec
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
      integer                           :: itrc
!
      real (r8)                         :: cff, cff1, cff2, cff3, cff4
      real (r8)                         :: cffx, cffy

# ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: tl_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAde
# endif
!
!  Initialize.
!
      Istr=BOUNDS(ng)%Istr(tile)     ! tile computational indices
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('temp', 'salt')
          itrc = tracer_index(TRIM(StateVarName(ifield)))
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%t_Bcorr(is:ie,js:je,1,ms,itrc)
          BscaleY(is:,js:) => self%t_Bcorr(is:ie,js:je,2,ms,itrc)
# endif
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
# ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
# else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
# endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Set operator initial conditions.
!
!^    CALL dabc_r3d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    A)4
!^
      CALL dabc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    tl_A)
# ifdef DISTRIBUTE
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
# endif
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!^          Awrk(i,j,k)=A(i,j,k)
!^
            tl_Awrk(i,j,k)=tl_A(i,j,k)
          END DO
        END DO
      END DO

# ifdef GEOPOTENTIAL_HCONV
!
!  Diffusion along geopotential surfaces.
!  =====================================
!
!  Compute horizontal and vertical gradients. Notice the recursive
!  blocking sequence. The vertical placement of the gradients is:
!
!        dAdx,dAde(:,:,k1) k     rho-points
!        dAdx,dAde(:,:,k2) k+1   rho-points
!          FZ,dAdz(:,:,k1) k-1/2   W-points
!          FZ,dAdz(:,:,k2) k+1/2   W-points
!
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%umask(i,j)
#  endif
              dZdx(i,j,k2)=cff*(GRID(ng)%z_r(i  ,j,k+1)-                &
     &                          GRID(ng)%z_r(i-1,j,k+1))
#  ifdef MASKING
!^            dAdx(i,j,k2)=cff*                                         &
!^   &                     (Awrk(i  ,j,k+1)*GRID(ng)%rmask(i  ,j)-      &
!^   &                      Awrk(i-1,j,k+1)*GRID(ng)%rmask(i-1,j))
!^
              tl_dAdx(i,j,k2)=cff*                                      &
     &                        (tl_Awrk(i  ,j,k+1)*                      &
     &                            GRID(ng)%rmask(i  ,j)-                &
     &                         tl_Awrk(i-1,j,k+1)*                      &
     &                            GRID(ng)%rmask(i-1,j))
#  else
!^            dAdx(i,j,k2)=cff*(Awrk(i  ,j,k+1)-                        &
!^   &                          Awrk(i-1,j,k+1))
!^
              tl_dAdx(i,j,k2)=cff*(tl_Awrk(i  ,j,k+1)-                  &
     &                             tl_Awrk(i-1,j,k+1))
#  endif
            END DO
          END DO
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%vmask(i,j)
#  endif
              dZde(i,j,k2)=cff*(GRID(ng)%z_r(i,j  ,k+1)-                &
     &                          GRID(ng)%z_r(i,j-1,k+1))
#  ifdef MASKING
!^            dAde(i,j,k2)=cff*
!^   &                     (Awrk(i,j  ,k+1)*GRID(ng)%rmask(i,j  )-      &
!^   &                      Awrk(i,j-1,k+1)*GRID(ng)%rmask(i,j-1))
!^
              tl_dAde(i,j,k2)=cff*                                      &
     &                        (tl_Awrk(i,j  ,k+1)*                      &
     &                            GRID(ng)%rmask(i,j  )-                &
     &                         tl_Awrk(i,j-1,k+1)*                      &
     &                            GRID(ng)%rmask(i,j-1))
#  else
!^            dAde(i,j,k2)=cff*(Awrk(i,j  ,k+1)-                        &
!^   &                          Awrk(i,j-1,k+1))
!^
              tl_dAde(i,j,k2)=cff*(tl_Awrk(i,j  ,k+1)-                  &
     &                             tl_Awrk(i,j-1,k+1))
#  endif
            END DO
          END DO
        END IF
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=Istr-1,Iend+1
!^            dAdz(i,j,k2)=0.0_r8
!^
              tl_dAdz(i,j,k2)=0.0_r8
!^            FZ(i,j,k2)=0.0_r8
!^
              tl_FZ(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(GRID(ng)%z_r(i,j,k+1)-GRID(ng)%z_r(i,j,k))
!^            dAdz(i,j,k2)=cff*(Awrk(i,j,k+1)-                          &
!^   &                          Awrk(i,j,k))
!^
              tl_dAdz(i,j,k2)=cff*(tl_Awrk(i,j,k+1)-                    &
     &                             tl_Awrk(i,j,k))
#  ifdef MASKING
!^            dAdz(i,j,k2)=dAdz(i,j,k2)*GRID(ng)%rmask(i,j)
!^
              tl_dAdz(i,j,k2)=tl_dAdz(i,j,k2)*GRID(ng)%rmask(i,j)
#  endif
            END DO
          END DO
        END IF
!
!  Compute components of the rotated A flux along geopotential
!  surfaces.
!
        IF (k.gt.0) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.5_r8*(Khx(i,j)+Khx(i-1,j))*GRID(ng)%on_u(i,j)
              cff1=MIN(dZdx(i,j,k1),0.0_r8)
              cff2=MAX(dZdx(i,j,k1),0.0_r8)
!^            FX(i,j)=cff*                                              &
!^   &                (dAdx(i,j,k1)-                                    &
!^   &                 0.5_r8*(cff1*(dAdz(i-1,j,k1)+                    &
!^   &                               dAdz(i  ,j,k2))+                   &
!^   &                         cff2*(dAdz(i-1,j,k2)+                    &
!^   &                               dAdz(i  ,j,k1))))
!^
              tl_FX(i,j)=cff*                                           &
     &                   (tl_dAdx(i,j,k1)-                              &
     &                    0.5_r8*(cff1*(tl_dAdz(i-1,j,k1)+              &
     &                                  tl_dAdz(i  ,j,k2))+             &
     &                            cff2*(tl_dAdz(i-1,j,k2)+              &
     &                                  tl_dAdz(i  ,j,k1))))
            END DO
          END DO
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.5_r8*(Khy(i,j-1)+Khy(i,j))*GRID(ng)%om_v(i,j)
              cff1=MIN(dZde(i,j,k1),0.0_r8)
              cff2=MAX(dZde(i,j,k1),0.0_r8)
!^            FE(i,j)=cff*                                              &
!^   &                (dAde(i,j,k1)-                                    &
!^   &                 0.5_r8*(cff1*(dAdz(i,j-1,k1)+                    &
!^   &                               dAdz(i,j  ,k2))+                   &
!^   &                         cff2*(dAdz(i,j-1,k2)+                    &
!^   &                               dAdz(i,j  ,k1))))
!^
              tl_FE(i,j)=cff*                                           &
     &                   (tl_dAde(i,j,k1)-                              &
     &                    0.5_r8*(cff1*(tl_dAdz(i,j-1,k1)+              &
     &                                  tl_dAdz(i,j  ,k2))+             &
     &                            cff2*(tl_dAdz(i,j-1,k2)+              &
     &                                  tl_dAdz(i,j  ,k1))))
            END DO
          END DO
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff=0.5_r8*Khx(i,j)
                cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
!^              FZ(i,j,k2)=cff*                                         &
!^   &                     (cff1*(cff1*dAdz(i,j,k2)-dAdx(i  ,j,k1))+    &
!^   &                      cff2*(cff2*dAdz(i,j,k2)-dAdx(i+1,j,k2))+    &
!^   &                      cff3*(cff3*dAdz(i,j,k2)-dAdx(i  ,j,k2))+    &
!^   &                      cff4*(cff4*dAdz(i,j,k2)-dAdx(i+1,j,k1)))
!^
                tl_FZ(i,j,k2)=cff*                                      &
     &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i  ,j,k1))+                &
     &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i+1,j,k2))+                &
     &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i  ,j,k2))+                &
     &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i+1,j,k1)))
!
                cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde(i,j+1,k1),0.0_r8)
                cff=0.5_r8*Khy(i,j)
!^              FZ(i,j,k2)=FZ(i,j,k2)+                                  &
!^   &                     cff*                                         &
!^   &                     (cff1*(cff1*dAdz(i,j,k2)-dAde(i,j  ,k1))+    &
!^   &                      cff2*(cff2*dAdz(i,j,k2)-dAde(i,j+1,k2))+    &
!^   &                      cff3*(cff3*dAdz(i,j,k2)-dAde(i,j  ,k2))+    &
!^   &                      cff4*(cff4*dAdz(i,j,k2)-dAde(i,j+1,k1)))
!^
                tl_FZ(i,j,k2)=tl_FZ(i,j,k2)+                            &
     &                        cff*                                      &
     &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j  ,k1))+                &
     &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j+1,k2))+                &
     &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j  ,k2))+                &
     &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j+1,k1)))
              END DO
            END DO
          END IF
!
!  Compute horizontal K-Laplacian operator.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
!^            Awrk(i,j,k)=Awrk(i,j,k)-                                  &
!^   &                    Hfac(i,j)*                                    &
!^   &                    (FX(i+1,j  )-FX(i,j)+                         &
!^   &                     FE(i  ,j+1)-FE(i,j))-                        &
!^   &                    (FZ(i,j,k2)-FZ(i,j,k1))/GRID(ng)%Hz(i,j,k)
!^
              tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                            &
     &                       Hfac(i,j)*                                 &
     &                       (tl_FX(i+1,j  )-tl_FX(i,j)+                &
     &                        tl_FE(i  ,j+1)-tl_FE(i,j))-               &
     &                       (tl_FZ(i,j,k2)-tl_FZ(i,j,k1))/             &
     &                       GRID(ng)%Hz(i,j,k)
            END DO
          END DO
        END IF
      END DO K_LOOP

#else

!
!  Diffusion along S-coordinates.
!  =============================
!
!  Compute XI- and ETA-components of diffusive flux.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
!^          FX(i,j)=GRID(ng)%pmon_u(i,j)*                               &
!^   &              0.5_r8*(Khx(i-1,j)+Khx(i,j))*                       &
!^   &              (Awrk(i,j,k)-Awrk(i-1,j,k))
!^
            tl_FX(i,j)=GRID(ng)%pmon_u(i,j)*                            &
     &                 0.5_r8*(Khx(i-1,j)+Khx(i,j))*                    &
     &                 (tl_Awrk(i,j,k)-tl_Awrk(i-1,j,k))
#  ifdef MASKING
!^          FX(i,j)=FX(i,j)*GRID(ng)%umask(i,j)
!^
            tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%umask(i,j)
#  endif
          END DO
        END DO
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
!^          FE(i,j)=GRID(ng)%pnom_v(i,j)*                               &
!^                  0.5_r8*(Khy(i,j-1)+Khy(i,j))*                       &
!^   &              (Awrk(i,j,k)-Awrk(i,j-1,k))
!^
            tl_FE(i,j)=GRID(ng)%pnom_v(i,j)*                            &
     &                 0.5_r8*(Khy(i,j-1)+Khy(i,j))*                    &
     &                 (tl_Awrk(i,j,k)-tl_Awrk(i,j-1,k))
#  ifdef MASKING
!^          FE(i,j)=FE(i,j)*GRID(ng)%vmask(i,j)
!^
            tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%vmask(i,j)
#  endif
          END DO
        END DO
!
!  Compute horizontal K-Laplacian operator.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          Awrk(i,j,k)=Awrk(i,j,k)-                                    &
!^   &                  Hfac(i,j)*                                      &
!^   &                  (FX(i+1,j)-FX(i,j)+                             &
!^   &                   FE(i,j+1)-FE(i,j))
!^
            tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                              &
     &                     Hfac(i,j)*                                   &
     &                     (tl_FX(i+1,j)-tl_FX(i,j)+                    &
     &                      tl_FE(i,j+1)-tl_FE(i,j))
          END DO
        END DO
      END DO
# endif
!
!  Load K-Laplacian solution.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          A(i,j,k)=Awrk(i,j,k)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k)
          END DO
        END DO
      END DO

# ifdef DISTRIBUTE
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
# endif

# ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_r3d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of a
!  3D control variable at RHO-points.
!
      SUBROUTINE multiscale_Klap_r3d_ad (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Mlap, i, j, k, rec
      integer                           :: kk, kt, k1, k1b, k2, k2b
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
      integer                           :: itrc
!
      real (r8)                         :: adfac, adfac1, adfac2
      real (r8)                         :: cff, cff1, cff2, cff3, cff4
      real (r8)                         :: cffx, cffy

# ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: ad_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAde
# endif
!
!  Initialize.
!
      Istr=BOUNDS(ng)%Istr(tile)     ! tile computational indices
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
# ifdef GEOPOTENTIAL_HCONV
      ad_FZ=0.0_r8
      ad_dAdz=0.0_r8
      ad_dAdx=0.0_r8
      ad_dAde=0.0_r8
# endif
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('temp', 'salt')
          itrc = tracer_index(TRIM(StateVarName(ifield)))
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%t_Bcorr(is:ie,js:je,1,ms,itrc)
          BscaleY(is:,js:) => self%t_Bcorr(is:ie,js:je,2,ms,itrc)
# endif
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
# ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
# else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
# endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Adjoint of load K-Laplacian solution.
!
# ifdef DISTRIBUTE
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
# endif
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          tl_A(i,j,k)=tl_Awrk(i,j,k)
!^
            ad_Awrk(i,j,k)=ad_Awrk(i,j,k)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

# ifdef GEOPOTENTIAL_HCONV
!
!  Adjoint of diffusion along geopotential surfaces.
!  ================================================
!
!  Compute horizontal and vertical gradients. Notice the recursive
!  blocking sequence. The vertical placement of the gradients is:
!
!        dAdx,dAde(:,:,k1) k     rho-points
!        dAdx,dAde(:,:,k2) k+1   rho-points
!          FZ,dAdz(:,:,k1) k-1/2   W-points
!          FZ,dAdz(:,:,k2) k+1/2   W-points
!
!  Compute adjoint of starting values of k1 and k2.
!
      k1=2
      k2=1
      DO k=0,N(ng)
!
!  Note: The following code is equivalent to
!
!       kt=k1
!       k1=k2
!       k2=kt
!
!  We use the adjoint of above code.
!
        k1=k2
        k2=3-k1
      END DO
!
      K_LOOP : DO k=N(ng),0,-1
!
!  Compute required BASIC STATE fields. Need to look forward in
!  recursive kk index.
!
        k2b=1
        DO kk=0,k
          k1b=k2b
          k2b=3-k1b
!
!  Compute components of the rotated tracer flux (A m3/s) along
!  geopotential surfaces (required BASIC STATE fields).
!
          IF (kk.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
                cff=cff*GRID(ng)%umask(i,j)
#  endif
                dZdx(i,j,k2b)=cff*(GRID(ng)%z_r(i  ,j,kk+1)-            &
     &                             GRID(ng)%z_r(i-1,j,kk+1))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend+1
                  dZdx(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
            DO j=Jstr,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
                cff=cff*GRID(ng)%vmask(i,j)
#  endif
                dZde(i,j,k2b)=cff*(GRID(ng)%z_r(i,j  ,kk+1)-            &
     &                             GRID(ng)%z_r(i,j-1,kk+1))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=Jstr,Jend+1
                DO i=Istr,Iend
                  dZde(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END DO
!
        IF (k.gt.0) THEN
!
!  Adjoint of compute K-Laplacian geopotential operator.
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
!^            tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                            &
!^   &                       Hfac(i,j)*                                 &
!^   &                       (tl_FX(i+1,j  )-tl_FX(i,j)+                &
!^   &                        tl_FE(i  ,j+1)-tl_FE(i,j))-               &
!^   &                       (tl_FZ(i,j,k2)-tl_FZ(i,j,k1))/             &
!^   &                       GRID(ng)%Hz(i,j,k)
!^
              adfac1=-Hfac(i,j)*ad_Awrk(i,j,k)
              adfac2=-ad_Awrk(i,j,k)/Hz(i,j,k)
              ad_FE(i,j  )=ad_FE(i,j  )-adfac1
              ad_FE(i,j+1)=ad_FE(i,j+1)+adfac1
              ad_FX(i  ,j)=ad_FX(i  ,j)-adfac1
              ad_FX(i+1,j)=ad_FX(i+1,j)+adfac1
              ad_FZ(i,j,k1)=ad_FZ(i,j,k1)-adfac2
              ad_FZ(i,j,k2)=ad_FZ(i,j,k2)+adfac2
            END DO
          END DO
!
!  Adjoint of compute components of the rotated A flux (A m3/s) along
!  geopotential surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=Istr,Iend
                cff=0.5_r8*Khy(i,j)
                cff1=MIN(dZde(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde(i,j+1,k1),0.0_r8)
!^              tl_FZ(i,j,k2)=tl_FZ(i,j,k2)+                            &
!^   &                        cff*                                      &
!^   &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j  ,k1))+                &
!^   &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j+1,k2))+                &
!^   &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j  ,k2))+                &
!^   &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j+1,k1)))
!^
                adfac=cff*ad_FZ(i,j,k2)
                ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)+                        &
     &                          (cff1*cff1+                             &
     &                           cff2*cff2+                             &
     &                           cff3*cff3+                             &
     &                           cff4*cff4)*adfac
                ad_dAde(i,j  ,k1)=ad_dAde(i,j  ,k1)-cff1*adfac
                ad_dAde(i,j+1,k2)=ad_dAde(i,j+1,k2)-cff2*adfac
                ad_dAde(i,j  ,k2)=ad_dAde(i,j  ,k2)-cff3*adfac
                ad_dAde(i,j+1,k1)=ad_dade(i,j+1,k1)-cff4*adfac
!
                cff1=MIN(dZdx(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx(i+1,j,k1),0.0_r8)
                cff=0.5_r8*Khx(i,j)
!
!^              tl_FZ(i,j,k2)=cff*                                      &
!^   &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAdx(i  ,j,k1))+                &
!^   &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAdx(i+1,j,k2))+                &
!^   &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAdx(i  ,j,k2))+                &
!^   &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAdx(i+1,j,k1)))
!^
                ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)+                        &
     &                          (cff1*cff1+                             &
     &                           cff2*cff2+                             &
     &                           cff3*cff3+                             &
     &                           cff4*cff4)*adfac
                ad_dAdx(i  ,j,k1)=ad_dAdx(i  ,j,k1)-cff1*adfac
                ad_dAdx(i+1,j,k2)=ad_dAdx(i+1,j,k2)-cff2*adfac
                ad_dAdx(i  ,j,k2)=ad_dAdx(i  ,j,k2)-cff3*adfac
                ad_dAdx(i+1,j,k1)=ad_dAdx(i+1,j,k1)-cff4*adfac
                ad_FZ(i,j,k2)=0.0_r8
              END DO
            END DO
          END IF
!
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.5_r8*(Khy(i,j-1)+Khy(i,j))*GRID(ng)%om_v(i,j)
              cff1=MIN(dZde(i,j,k1),0.0_r8)
              cff2=MAX(dZde(i,j,k1),0.0_r8)
!^            tl_FE(i,j)=cff*                                           &
!^   &                   (tl_dAde(i,j,k1)-                              &
!^   &                    0.5_r8*(cff1*(tl_dAdz(i,j-1,k1)+              &
!^   &                                  tl_dAdz(i,j  ,k2))+             &
!^   &                            cff2*(tl_dAdz(i,j-1,k2)+              &
!^   &                                  tl_dAdz(i,j  ,k1))))
!^
              adfac=cff*ad_FE(i,j)
              adfac1=adfac*0.5_r8*cff1
              adfac2=adfac*0.5_r8*cff2
              ad_dAde(i,j,k1)=ad_dAde(i,j,k1)+adfac
              ad_dAdz(i,j-1,k1)=ad_dAdz(i,j-1,k1)-adfac1
              ad_dAdz(i,j  ,k2)=ad_dAdz(i,j  ,k2)-adfac1
              ad_dAdz(i,j-1,k2)=ad_dAdz(i,j-1,k2)-adfac2
              ad_dAdz(i,j  ,k1)=ad_dAdz(i,j  ,k1)-adfac2
              ad_FE(i,j)=0.0_r8
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.5_r8*(Khx(i,j)+Khx(i-1,j))*GRID(ng)%on_u(i,j)
              cff1=MIN(dZdx(i,j,k1),0.0_r8)
              cff2=MAX(dZdx(i,j,k1),0.0_r8)
!^            tl_FX(i,j)=cff*                                           &
!^   &                   (tl_dAdx(i,j,k1)-                              &
!^   &                    0.5_r8*(cff1*(tl_dAdz(i-1,j,k1)+              &
!^   &                                  tl_dAdz(i  ,j,k2))+             &
!^   &                            cff2*(tl_dAdz(i-1,j,k2)+              &
!^   &                                  tl_dAdz(i  ,j,k1))))
!^
              adfac=cff*ad_FX(i,j)
              adfac1=adfac*0.5_r8*cff1
              adfac2=adfac*0.5_r8*cff2
              ad_dAdx(i,j,k1)=ad_dAdx(i,j,k1)+adfac
              ad_dAdz(i-1,j,k1)=ad_dAdz(i-1,j,k1)-adfac1
              ad_dAdz(i  ,j,k2)=ad_dAdz(i  ,j,k2)-adfac1
              ad_dAdz(i-1,j,k2)=ad_dAdz(i-1,j,k2)-adfac2
              ad_dAdz(i  ,j,k1)=ad_dAdz(i  ,j,k1)-adfac2
              ad_FX(i,j)=0.0_r8
            END DO
          END DO
        END IF
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=Istr-1,Iend+1
!^            tl_FZ(i,j,k2)=0.0_r8
!^
              ad_FZ(i,j,k2)=0.0_r8
!^            tl_dAdz(i,j,k2)=0.0_r8
!^
              ad_dAdz(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(GRID(ng)%z_r(i,j,k+1)-GRID(ng)%z_r(i,j,k))
#  ifdef MASKING
!^            tl_dAdz(i,j,k2)=tl_dAdz(i,j,k2)*GRID(ng)%rmask(i,j)
!^
              ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)*GRID(ng)%rmask(i,j)
#  endif
!^            tl_dAdz(i,j,k2)=cff*(tl_Awrk(i,j,k+1)-                    &
!^   &                             tl_Awrk(i,j,k  ))
!^
              adfac=cff*ad_dAdz(i,j,k2)
              ad_Awrk(i,j,k)=ad_Awrk(i,j,k)-adfac
              ad_Awrk(i,j,k+1)=ad_Awrk(i,j,k+1)+adfac
              ad_dAdz(i,j,k2)=0.0_r8
            END DO
          END DO
        END IF
!
        IF (k.lt.N(ng)) THEN
          DO j=Jstr,Jend+1
            DO i=Istr,Iend
              cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%vmask(i,j)
!^            tl_dAde(i,j,k2)=cff*
!^   &                        (tl_Awrk(i,j  ,k+1)*                      &
!^   &                            GRID(ng)%rmask(i,j  )-                &
!^   &                         tl_Awrk(i,j-1,k+1)*                      &
!^   &                            GRID(ng)%rmask(i,j-1))
!^
              adfac=cff*ad_dAde(i,j,k2)
              ad_Awrk(i,j-1,k+1)=ad_Awrk(i,j-1,k+1)-                    &
     &                           GRID(ng)%rmask(i,j-1)*adfac
              ad_Awrk(i,j  ,k+1)=ad_Awrk(i,j  ,k+1)+                    &
     &                           GRID(ng)%rmask(i,j  )*adfac
              ad_dAde(i,j,k2)=0.0_r8
#  else
!^            tl_dAde(i,j,k2)=cff*(tl_Awrk(i,j  ,k+1)-                  &
!^   &                             tl_Awrk(i,j-1,k+1))
!^
              adfac=cff*ad_dAde(i,j,k2)
              ad_Awrk(i,j-1,k+1)=ad_Awrk(i,j-1,k+1)-adfac
              ad_Awrk(i,j  ,k+1)=ad_Awrk(i,j  ,k+1)+adfac
              ad_dAde(i,j,k2)=0.0_r8
#  endif
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=Istr,Iend+1
              cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
              cff=cff*umask(i,j)
!^            tl_dAdx(i,j,k2)=cff*                                      &
!^   &                        (tl_Awrk(i  ,j,k+1)*                      &
!^   &                            GRID(ng)%rmask(i  ,j)-                &
!^   &                         tl_Awrk(i-1,j,k+1)*                      &
!^   &                            GRID(ng)%rmask(i-1,j))
!^
              adfac=cff*ad_dAdx(i,j,k2)
              ad_Awrk(i-1,j,k+1)=ad_Awrk(i-1,j,k+1)-                    &
     &                           GRID(ng)%rmask(i-1,j)*adfac
              ad_Awrk(i  ,j,k+1)=ad_Awrk(i  ,j,k+1)+                    &
     &                           GRID(ng)%rmask(i  ,j)*adfac
              ad_dAdx(i,j,k2)=0.0_r8
#  else
!^            tl_dAdx(i,j,k2)=cff*(tl_Awrk(i  ,j,k+1)-                  &
!^   &                             tl_Awrk(i-1,j,k+1))
!^
              adfac=cff*ad_dAdx(i,j,k2)
              ad_Awrk(i-1,j,k+1)=ad_Awrk(i-1,j,k+1)-adfac
              ad_Awrk(i  ,j,k+1)=ad_Awrk(i  ,j,k+1)+adfac
              ad_dAdx(i,j,k2)=0.0_r8
#  endif
            END DO
          END DO
        END IF
!
!  Compute new storage recursive indices.
!
        kt=k2
        k2=k1
        k1=kt
      END DO K_LOOP

# else
!
!  Adjoint of diffusion along S-coordinates.
!  ========================================
!
!  Adjoint of compute horizontal K-Laplacian operator.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=Istr,Iend
!^          tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                              &
!^   &                     Hfac(i,j)*                                   &
!^   &                     (tl_FX(i+1,j)-tl_FX(i,j)+                    &
!^   &                      tl_FE(i,j+1)-tl_FE(i,j))
!^
            adfac=-Hfac(i,j)*ad_Awrk(i,j,k)
            ad_FE(i,j  )=ad_FE(i,j  )-adfac
            ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
            ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
            ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
          END DO
        END DO
!
!  Adjoint of compute XI- and ETA-components of diffusive flux.
!
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
#  ifdef MASKING
!^          tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%vmask(i,j)
!^
            ad_FE(i,j)=ad_FE(i,j)*GRID(ng)%vmask(i,j)
#  endif
!^          tl_FE(i,j)=GRID(ng)%pnom_v(i,j)*                            &
!^   &                 0.5_r8*(Khy(i,j-1)+Khy(i,j))*                    &
!^   &                 (tl_Awrk(i,j,k)-tl_Awrk(i,j-1,k))
!^
            adfac=GRID(ng)%pnom_v(i,j)*                                 &
     &            0.5_r8*(Khy(i,j-1)+Khy(i,j))*ad_FE(i,j)
            ad_Awrk(i,j-1,k)=ad_Awrk(i,j-1,k)-adfac
            ad_Awrk(i,j  ,k)=ad_Awrk(i,j  ,k)+adfac
            ad_FE(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
#  ifdef MASKING
!^          tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%umask(i,j)
!^
            ad_FX(i,j)=ad_FX(i,j)*GRID(ng)%umask(i,j)
#  endif
!^          tl_FX(i,j)=GRID(ng)%pmon_u(i,j)*                            &
!^   &                 0.5_r8*(Khx(i-1,j)+Khx(i,j))*                    &
!^   &                 (tl_Awrk(i,j,k)-tl_Awrk(i-1,j,k))
!^
            adfac=GRID(ng)%pmon_u(i,j)*                                 &
     &            0.5_r8*(Khx(i-1,j)+Khx(i,j))*ad_FX(i,j)
            ad_Awrk(i-1,j,k)=ad_Awrk(i-1,j,k)-adfac
            ad_Awrk(i  ,j,k)=ad_Awrk(i  ,j,k)+adfac
            ad_FX(i,j)=0.0_r8
          END DO
        END DO
      END DO
# endif
!
!  Adjoint of set operator initial conditions.
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!^          tl_Awrk(i,j,k)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k)
            ad_Awrk(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

# ifdef DISTRIBUTE
!
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
# endif
!
!^    CALL dabc_r3d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    tl_A)
!^
      CALL ad_dabc_r3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       ad_A)

# ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_r3d_ad
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of a
!  3D control variable at U-points.
!
      SUBROUTINE multiscale_Klap_u3d_tl (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Mlap, i, j, k, k1, k2, rec
      integer                           :: IstrU, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: cff, cff1, cff2, cff3, cff4
      real (r8)                         :: cffx, cffy

# ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: tl_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS)   :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS)   :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAde
# endif
!
!  Initialize.
!
      IstrU=BOUNDS(ng)%IstrU(tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%u_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%u_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
# ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
# else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
# endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Set operator initial conditions.
!
!^    CALL dabc_u3d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    A)
!^
      CALL dabc_u3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    tl_A)
# ifdef DISTRIBUTE
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
# endif
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
!^          Awrk(i,j,k)=A(i,j,k)
!^
            tl_Awrk(i,j,k)=tl_A(i,j,k)
          END DO
        END DO
      END DO

# ifdef GEOPOTENTIAL_HCONV
!
!  Diffusion along geopotential surfaces.
!  =====================================
!
!  Compute horizontal and vertical gradients. Notice the recursive
!  blocking sequence. The vertical placement of the gradients is:
!
!        dZdx_r,dAdx,dAde(:,:,k1) k     rho-points
!        dZdx_r,dAdx,dAde(:,:,k2) k+1   rho-points
!                  dZde_p(:,:,k1) k     psi-points
!                  dZde_p(:,:,k2) k+1   psi-points
!                 FZ,dAdz(:,:,k1) k-1/2   W-points
!                 FZ,dAdz(:,:,k2) k+1/2   W-points
!
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend+1
              cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%umask(i,j)
#  endif
              dZdx(i,j)=cff*(GRID(ng)%z_r(i  ,j,k+1)-                   &
     &                       GRID(ng)%z_r(i-1,j,k+1))
            END DO
          END DO
!
          DO j=Jstr,Jend+1
            DO i=IstrU-1,Iend
              cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%vmask(i,j)
#  endif
              dZde(i,j)=cff*(GRID(ng)%z_r(i,j  ,k+1)-                   &
     &                       GRID(ng)%z_r(i,j-1,k+1))
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
#  ifdef MASKING
!^            dAdx(i,j,k2)=GRID(ng)%pm(i,j)*                            &
!^   &                     (Awrk(i+1,j,k+1)*GRID(ng)%umask(i+1,j)-      &
!^   &                      Awrk(i  ,j,k+1)*GRID(ng)%umask(i  ,j))
!^
              tl_dAdx(i,j,k2)=GRID(ng)%pm(i,j)*                         &
     &                        (tl_Awrk(i+1,j,k+1)*                      &
     &                            GRID(ng)%umask(i+1,j)-                &
     &                         tl_Awrk(i  ,j,k+1)*                      &
     &                            GRID(ng)%umask(i  ,j))
!^            dAdx(i,j,k2)=dAdx(i,j,k2)*GRID(ng)%rmask(i,j)
!^
              tl_dAdx(i,j,k2)=tl_dAdx(i,j,k2)*GRID(ng)%rmask(i,j)
#  else
!^            dAdx(i,j,k2)=GRID(ng)%pm(i,j)*(Awrk(i+1,j,k+1)-           &
!^   &                                       Awrk(i  ,j,k+1))
!^
              tl_dAdx(i,j,k2)=GRID(ng)%pm(i,j)*(tl_Awrk(i+1,j,k+1)-     &
     &                                          tl_Awrk(i  ,j,k+1))
#  endif
              dZdx_r(i,j,k2)=0.5_r8*(dZdx(i  ,j)+                       &
     &                               dZdx(i+1,j))
            END DO
          END DO
!
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              cff=0.25_r8*(GRID(ng)%pn(i-1,j  )+GRID(ng)%pn(i,j  )+     &
     &                     GRID(ng)%pn(i-1,j-1)+GRID(ng)%pn(i,j-1))
#  ifdef MASKING
!^            dAde(i,j,k2)=cff*                                         &
!^   &                     (Awrk(i,j  ,k+1)*GRID(ng)%umask(i,j  )-      &
!^   &                      Awrk(i,j-1,k+1)*GRID(ng)%umask(i,j-1))
!^
              tl_dAde(i,j,k2)=cff*                                      &
     &                        (tl_Awrk(i,j  ,k+1)*                      &
     &                            GRID(ng)%umask(i,j  )-                &
     &                         tl_Awrk(i,j-1,k+1)*                      &
     &                            GRID(ng)%umask(i,j-1))
!^            dAde(i,j,k2)=dAde(i,j,k2)*GRID(ng)%pmask(i,j)
!^
              tl_dAde(i,j,k2)=tl_dAde(i,j,k2)*GRID(ng)%pmask(i,j)
#  else
!^            dAde(i,j,k2)=cff*(Awrk(i,j  ,k+1)-                        &
!^   &                          Awrk(i,j-1,k+1))
!^
              tl_dAde(i,j,k2)=cff*(tl_Awrk(i,j  ,k+1)-                  &
     &                             tl_Awrk(i,j-1,k+1))
#  endif
              dZde_p(i,j,k2)=0.5_r8*(dZde(i-1,j)+                       &
     &                               dZde(i  ,j))
            END DO
          END DO
        END IF
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
!^            dAdz(i,j,k2)=0.0_r8
!^
              tl_dAdz(i,j,k2)=0.0_r8
!^            FZ(i,j,k2)=0.0_r8
!^
              tl_FZ(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(GRID(ng)%z_r(i,j,k+1)-GRID(ng)%z_r(i,j,k))
!^            dAdz(i,j,k2)=cff*(Awrk(i,j,k+1)-                          &
!^   &                          Awrk(i,j,k))
!^
              tl_dAdz(i,j,k2)=cff*(tl_Awrk(i,j,k+1)-                    &
     &                             tl_Awrk(i,j,k  ))
#  ifdef MASKING
!^            dAdz(i,j,k2)=dAdz(i,j,k2)*GRID(ng)%umask(i,j)
!^
              tl_dAdz(i,j,k2)=tl_dAdz(i,j,k2)*GRID(ng)%umask(i,j)
#  endif
            END DO
          END DO
        END IF
!
!  Compute components of the rotated A flux along geopotential
!  surfaces.
!
        IF (k.gt.0) THEN
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
              cff=Khx(i,j)*GRID(ng)%on_r(i,j)
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
!^            FX(i,j)=cff*                                              &
!^   &                (dAdx(i,j,k1)-                                    &
!^   &                 0.5_r8*(cff1*(dAdz(i  ,j,k1)+                    &
!^   &                               dAdz(i+1,j,k2))+                   &
!^   &                         cff2*(dAdz(i  ,j,k2)+                    &
!^   &                               dAdz(i+1,j,k1))))
!^
              tl_FX(i,j)=cff*                                           &
     &                   (tl_dAdx(i,j,k1)-                              &
     &                    0.5_r8*(cff1*(tl_dAdz(i  ,j,k1)+              &
     &                                  tl_dAdz(i+1,j,k2))+             &
     &                            cff2*(tl_dAdz(i  ,j,k2)+              &
     &                                  tl_dAdz(i+1,j,k1))))
            END DO
          END DO
!
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              cff=0.25_r8*(Khy(i-1,j-1)+Khy(i-1,j)+                     &
     &                     Khy(i  ,j-1)+Khy(i  ,j))*GRID(ng)%om_p(i,j)
              cff1=MIN(dZde_p(i,j,k1),0.0_r8)
              cff2=MAX(dZde_p(i,j,k1),0.0_r8)
!^            FE(i,j)=cff*                                              &
!^   &                (dAde(i,j,k1)-                                    &
!^   &                 0.5_r8*(cff1*(dAdz(i,j-1,k1)+                    &
!^   &                               dAdz(i,j  ,k2))+                   &
!^   &                         cff2*(dAdz(i,j-1,k2)+                    &
!^   &                               dAdz(i,j  ,k1))))
!^
              tl_FE(i,j)=cff*                                           &
     &                   (tl_dAde(i,j,k1)-                              &
     &                    0.5_r8*(cff1*(tl_dAdz(i,j-1,k1)+              &
     &                                  tl_dAdz(i,j  ,k2))+             &
     &                            cff2*(tl_dAdz(i,j-1,k2)+              &
     &                                  tl_dAdz(i,j  ,k1))))
            END DO
          END DO
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                cff=0.25_r8*(Khx(i-1,j)+Khx(i,j))
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
!^              FZ(i,j,k2)=cff*                                         &
!^   &                     (cff1*(cff1*dAdz(i,j,k2)-dAdx(i-1,j,k1))+    &
!^   &                      cff2*(cff2*dAdz(i,j,k2)-dAdx(i  ,j,k2))+    &
!^   &                      cff3*(cff3*dAdz(i,j,k2)-dAdx(i-1,j,k2))+    &
!^   &                      cff4*(cff4*dAdz(i,j,k2)-dAdx(i  ,j,k1)))
!^
                tl_FZ(i,j,k2)=cff*                                      &
     &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i-1,j,k1))+                &
     &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i  ,j,k2))+                &
     &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i-1,j,k2))+                &
     &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i  ,j,k1)))
!
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
                cff=0.25_r8*(Khy(i-1,j)+Khy(i,j))
!
!^              FZ(i,j,k2)=FZ(i,j,k2)+                                  &
!^   &                     cff*                                         &
!^   &                     (cff1*(cff1*dAdz(i,j,k2)-dAde(i,j  ,k1))+    &
!^   &                      cff2*(cff2*dAdz(i,j,k2)-dAde(i,j+1,k2))+    &
!^   &                      cff3*(cff3*dAdz(i,j,k2)-dAde(i,j  ,k2))+    &
!^   &                      cff4*(cff4*dAdz(i,j,k2)-dAde(i,j+1,k1)))
!^
                tl_FZ(i,j,k2)=tl_FZ(i,j,k2)+                            &
     &                        cff*                                      &
     &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j  ,k1))+                &
     &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j+1,k2))+                &
     &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j  ,k2))+                &
     &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j+1,k1)))
              END DO
            END DO
          END IF
!
!  Compute horizontal K-Laplacian operator.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!^            Awrk(i,j,k)=Awrk(i,j,k)-                                  &
!^   &                    Hfac(i,j)*                                    &
!^   &                    (FX(i,j  )-FX(i-1,j)+                         &
!^   &                     FE(i,j+1)-FE(i  ,j))-                        &
!^   &                    2.0_r8*(FZ(i,j,k2)-FZ(i,j,k1))/               &
!^   &                           (GRID(ng)%Hz(i  ,j,k)+                 &
!^   &                            GRID(ng)%Hz(i-1,j,k))
!^
              tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                            &
     &                       Hfac(i,j)*                                 &
     &                       (tl_FX(i,j  )-tl_FX(i-1,j)+                &
     &                        tl_FE(i,j+1)-tl_FE(i  ,j))-               &
     &                       2.0_r8*(tl_FZ(i,j,k2)-tl_FZ(i,j,k1))/      &
     &                              (GRID(ng)%Hz(i  ,j,k)+              &
                                     GRID(ng)%Hz(i-1,j,k))
            END DO
          END DO
        END IF
      END DO K_LOOP

# else

!
!  Diffusion along S-coordinates.
!  =============================
!
!  Compute XI- and ETA-components of diffusive flux.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!^          FX(i,j)=GRID(ng)%pmon_r(i,j)*                               &
!^   &              Khx(i,j)*(Awrk(i+1,j,k)-Awrk(i,j,k))
!^
            tl_FX(i,j)=GRID(ng)%pmon_r(i,j)*                            &
     &                 Khx(i,j)*(tl_Awrk(i+1,j,k)-tl_Awrk(i,j,k))
          END DO
        END DO
!
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
!^          FE(i,j)=GRID(ng)%pnom_p(i,j)*
!^   &              0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                   &
!^   &                       Khy(i-1,j-1)+Khy(i,j-1))*                  &
!^   &              (Awrk(i,j,k)-Awrk(i,j-1,k))
!^
            tl_FE(i,j)=GRID(ng)%pnom_p(i,j)*                            &
     &                 0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                &
     &                          Khy(i-1,j-1)+Khy(i,j-1))*               &
     &                 (tl_Awrk(i,j,k)-tl_Awrk(i,j-1,k))
#  ifdef MASKING
!^          FE(i,j)=FE(i,j)*GRID(ng)%pmask(i,j)
!^
            tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%pmask(i,j)
#  endif
          END DO
        END DO
!
!  Compute horizontal K-Laplacian operator.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!^          Awrk(i,j,k)=Awrk(i,j,k)-                                    &
!^   &                  Hfac(i,j)*                                      &
!^   &                  (FX(i,j)-FX(i-1,j)+                             &
!^   &                   FE(i,j+1)-FE(i,j))
!^
            tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                              &
     &                     Hfac(i,j)*                                   &
     &                     (tl_FX(i,j)-tl_FX(i-1,j)+                    &
     &                      tl_FE(i,j+1)-tl_FE(i,j))
          END DO
        END DO
      END DO
# endif
!
!  Load K-Laplacian solution.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!^          A(i,j,k)=Awrk(i,j,k)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k)
          END DO
        END DO
      END DO

# ifdef DISTRIBUTE
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
# endif

# ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_u3d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of a
!  3D control variable at U-points.
!
      SUBROUTINE multiscale_Klap_u3d_ad (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Mlap, i, j, k, rec
      integer                           :: kk, kt, k1, k1b, k2, k2b
      integer                           :: IstrU, Iend, Jstr, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: adfac, adfac1, adfac2
      real (r8)                         :: cff, cff1, cff2, cff3, cff4
      real (r8)                         :: cffx, cffy

# ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: ad_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAde
# endif
!
!  Initialize.
!
      IstrU=BOUNDS(ng)%IstrU(tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
# ifdef GEOPOTENTIAL_HCONV
      ad_FZ=0.0_r8
      ad_dAdz=0.0_r8
      ad_dAdx=0.0_r8
      ad_dAde=0.0_r8
# endif
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%u_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%u_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=Jstr-1,Jend+1
        DO i=IstrU-1,Iend+1
# ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
# else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
# endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Adjoint of load K-Laplacian solution.
!
# ifdef DISTRIBUTE
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
# endif
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!^          tl_A(i,j,k)=tl_Awrk(i,j,k)
!^
            ad_Awrk(i,j,k)=ad_Awrk(i,j,k)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

# ifdef GEOPOTENTIAL_HCONV
!
!  Adjoint of diffusion along geopotential surfaces.
!  ================================================
!
!  Compute horizontal and vertical gradients. Notice the recursive
!  blocking sequence. The vertical placement of the gradients is:
!
!        dZdx_r,dAdx,dAde(:,:,k1) k     rho-points
!        dZdx_r,dAdx,dAde(:,:,k2) k+1   rho-points
!                  dZde_p(:,:,k1) k     psi-points
!                  dZde_p(:,:,k2) k+1   psi-points
!                 FZ,dAdz(:,:,k1) k-1/2   W-points
!                 FZ,dAdz(:,:,k2) k+1/2   W-points
!
!  Compute adjoint of starting values of k1 and k2.
!
      k1=2
      k2=1
      DO k=0,N(ng)
!
!  Note: The following code is equivalent to
!
!       kt=k
!       k1=k2
!       k2=kt
!
!  We use the adjoint of above code.
!
        k1=k2
        k2=3-k1
      END DO
!
      K_LOOP : DO k=N(ng),0,-1

!  Compute required BASIC STATE fields. Need to look forward in
!  recursive kk index.
!
        k2b=1
        DO kk=0,k
          k1b=k2b
          k2b=3-k1b
!
!  Compute components of the rotated tracer flux (A m3/s) along
!  geopotential surfaces (required BASIC STATE fields).
!
          IF (kk.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU-1,Iend+1
                cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
                cff=cff*GRID(ng)%umask(i,j)
#  endif
                dZdx(i,j)=cff*(GRID(ng)%z_r(i  ,j,kk+1)-                &
     &                         GRID(ng)%z_r(i-1,j,kk+1))
              END DO
            END DO
            DO j=Jstr,Jend
              DO i=IstrU-1,Iend
                dZdx_r(i,j,k2b)=0.5_r8*(dZdx(i  ,j)+                    &
     &                                  dZdx(i+1,j))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=Jstr,Jend
                DO i=IstrU-1,Iend
                  dZdx_r(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
!
            DO j=Jstr,Jend+1
              DO i=IstrU-1,Iend
                cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
                  cff=cff*GRID(ng)%vmask(i,j)
#  endif
                  dZde(i,j)=cff*(GRID(ng)%z_r(i,j  ,kk+1)-              &
     &                           GRID(ng)%z_r(i,j-1,kk+1))
              END DO
            END DO
            DO j=Jstr,Jend+1
              DO i=IstrU,Iend
                dZde_p(i,j,k2b)=0.5_r8*(dZde(i-1,j)+                    &
     &                                  dZde(i  ,j))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=Jstr,Jend+1
                DO i=IstrU,Iend
                  dZde_p(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END DO
!
        IF (k.gt.0) THEN
!
!  Adjoint of compute geopotential K-Laplacian operator.
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
!^            tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                            &
!^   &                       Hfac(i,j)*                                 &
!^   &                       (tl_FX(i,j  )-tl_FX(i-1,j)+                &
!^   &                        tl_FE(i,j+1)-tl_FE(i  ,j))-               &
!^   &                       2.0_r8*(tl_FZ(i,j,k2)-tl_FZ(i,j,k1))/      &
!^   &                              (GRID(ng)%Hz(i  ,j,k)+              &
!^   &                               GRID(ng)%Hz(i-1,j,k))
!^
              adfac1=-Hfac(i,j)*ad_Awrk(i,j,k)
              adfac2=-2.0_r8*ad_Awrk(i,j,k)/                            &
     &               (GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i-1,j,k))
              ad_FE(i,j  )=ad_FE(i,j  )-adfac1
              ad_FE(i,j+1)=ad_FE(i,j+1)+adfac1
              ad_FX(i-1,j)=ad_FX(i-1,j)-adfac1
              ad_FX(i  ,j)=ad_FX(i  ,j)+adfac1
              ad_FZ(i,j,k1)=ad_FZ(i,j,k1)-adfac2
              ad_FZ(i,j,k2)=ad_FZ(i,j,k2)+adfac2
            END DO
          END DO
!
!  Adjoint of compute components of rotated A flux (A m3/s) along
!  geopotential surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                cff=0.25_r8*(Khy(i-1,j)+Khy(i,j))
                cff1=MIN(dZde_p(i,j  ,k1),0.0_r8)
                cff2=MIN(dZde_p(i,j+1,k2),0.0_r8)
                cff3=MAX(dZde_p(i,j  ,k2),0.0_r8)
                cff4=MAX(dZde_p(i,j+1,k1),0.0_r8)
!^              tl_FZ(i,j,k2)=tl_FZ(i,j,k2)+                            &
!^   &                        cff*                                      &
!^   &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j  ,k1))+                &
!^   &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j+1,k2))+                &
!^   &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j  ,k2))+                &
!^   &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j+1,k1)))
!^
                adfac=cff*ad_FZ(i,j,k2)
                ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)+                        &
     &                          (cff1*cff1+                             &
     &                           cff2*cff2+                             &
     &                           cff3*cff3+                             &
     &                           cff4*cff4)*adfac
                ad_dAde(i,j  ,k1)=ad_dAde(i,j  ,k1)-cff1*adfac
                ad_dAde(i,j+1,k2)=ad_dAde(i,j+1,k2)-cff2*adfac
                ad_dAde(i,j  ,k2)=ad_dAde(i,j  ,k2)-cff3*adfac
                ad_dAde(i,j+1,k1)=ad_dade(i,j+1,k1)-cff4*adfac
!
                cff1=MIN(dZdx_r(i-1,j,k1),0.0_r8)
                cff2=MIN(dZdx_r(i  ,j,k2),0.0_r8)
                cff3=MAX(dZdx_r(i-1,j,k2),0.0_r8)
                cff4=MAX(dZdx_r(i  ,j,k1),0.0_r8)
                cff=0.25_r8*(Khx(i-1,j)+Khx(i,j))
!^              tl_FZ(i,j,k2)=cff*                                      &
!^   &                          (cff1*(cff1*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i-1,j,k1))+              &
!^   &                           cff2*(cff2*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i  ,j,k2))+              &
!^   &                           cff3*(cff3*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i-1,j,k2))+              &
!^   &                           cff4*(cff4*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i  ,j,k1)))
!^
                ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)+                        &
     &                          (cff1*cff1+                             &
     &                           cff2*cff2+                             &
     &                           cff3*cff3+                             &
     &                           cff4*cff4)*adfac
                ad_dAdx(i-1,j,k1)=ad_dAdx(i-1,j,k1)-cff1*adfac
                ad_dAdx(i  ,j,k2)=ad_dAdx(i  ,j,k2)-cff2*adfac
                ad_dAdx(i-1,j,k2)=ad_dAdx(i-1,j,k2)-cff3*adfac
                ad_dAdx(i  ,j,k1)=ad_dAdx(i  ,j,k1)-cff4*adfac
                ad_FZ(i,j,k2)=0.0_r8
              END DO
            END DO
          END IF
!
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              cff=0.25_r8*(Khy(i-1,j-1)+Khy(i-1,j)+                     &
     &                     Khy(i  ,j-1)+Khy(i  ,j))*GRID(ng)%om_p(i,j)
              cff1=MIN(dZde_p(i,j,k1),0.0_r8)
              cff2=MAX(dZde_p(i,j,k1),0.0_r8)
!^            tl_FE(i,j)=cff*                                           &
!^   &                   (tl_dAde(i,j,k1)-                              &
!^   &                    0.5_r8*(cff1*(tl_dAdz(i,j-1,k1)+              &
!^   &                                  tl_dAdz(i,j  ,k2))+             &
!^   &                            cff2*(tl_dAdz(i,j-1,k2)+              &
!^   &                                  tl_dAdz(i,j  ,k1))))
!^
              adfac=cff*ad_FE(i,j)
              adfac1=adfac*0.5_r8*cff1
              adfac2=adfac*0.5_r8*cff2
              ad_dAde(i,j,k1)=ad_dAde(i,j,k1)+adfac
              ad_dAdz(i,j-1,k1)=ad_dAdz(i,j-1,k1)-adfac1
              ad_dAdz(i,j  ,k2)=ad_dAdz(i,j  ,k2)-adfac1
              ad_dAdz(i,j-1,k2)=ad_dAdz(i,j-1,k2)-adfac2
              ad_dAdz(i,j  ,k1)=ad_dAdz(i,j  ,k1)-adfac2
              ad_FE(i,j)=0.0_r8
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
              cff=Khx(i,j)*GRID(ng)%on_r(i,j)
              cff1=MIN(dZdx_r(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_r(i,j,k1),0.0_r8)
!^            tl_FX(i,j)=cff*                                           &
!^   &                   (tl_dAdx(i,j,k1)-                              &
!^   &                    0.5_r8*(cff1*(tl_dAdz(i  ,j,k1)+              &
!^   &                                  tl_dAdz(i+1,j,k2))+             &
!^   &                            cff2*(tl_dAdz(i  ,j,k2)+              &
!^   &                                  tl_dAdz(i+1,j,k1))))
!^
              adfac=cff*ad_FX(i,j)
              adfac1=adfac*0.5_r8*cff1
              adfac2=adfac*0.5_r8*cff2
              ad_dAdx(i,j,k1)=ad_dAdx(i,j,k1)+adfac
              ad_dAdz(i  ,j,k1)=ad_dAdz(i  ,j,k1)-adfac1
              ad_dAdz(i+1,j,k2)=ad_dAdz(i+1,j,k2)-adfac1
              ad_dAdz(i  ,j,k2)=ad_dAdz(i  ,j,k2)-adfac2
              ad_dAdz(i+1,j,k1)=ad_dAdz(i+1,j,k1)-adfac2
              ad_FX(i,j)=0.0_r8
            END DO
          END DO
        END IF
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
!^            tl_FZ(i,j,k2)=0.0_r8
!^
              ad_FZ(i,j,k2)=0.0_r8
!^            tl_dAdz(i,j,k2)=0.0_r8
!^
              ad_dAdz(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              cff=1.0_r8/(GRID(ng)%z_r(i,j,k+1)-GRID(ng)%z_r(i,j,k))
#  ifdef MASKING
!^            tl_dAdz(i,j,k2)=tl_dAdz(i,j,k2)*GRID(ng)%umask(i,j)
!^
              ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)*GRID(ng)%umask(i,j)
#  endif
!^            tl_dAdz(i,j,k2)=cff*(tl_Awrk(i,j,k+1)-                    &
!^   &                             tl_Awrk(i,j,k))
!^
              adfac=cff*ad_dAdz(i,j,k2)
              ad_Awrk(i,j,k)=ad_Awrk(i,j,k)-adfac
              ad_Awrk(i,j,k+1)=ad_Awrk(i,j,k+1)+adfac
              ad_dAdz(i,j,k2)=0.0_r8
            END DO
          END DO
        END IF
!
        IF (k.lt.N(ng)) THEN
          DO j=Jstr,Jend+1
            DO i=IstrU,Iend
              cff=0.25_r8*(GRID(ng)%pn(i-1,j  )+GRID(ng)%pn(i,j  )+     &
     &                     GRID(ng)%pn(i-1,j-1)+GRID(ng)%pn(i,j-1))
#  ifdef MASKING
!^            tl_dAde(i,j,k2)=tl_dAde(i,j,k2)*GRID(ng)%pmask(i,j)
!^
              ad_dAde(i,j,k2)=ad_dAde(i,j,k2)*GRID(ng)%pmask(i,j)
!^            tl_dAde(i,j,k2)=cff*                                      &
!^   &                        (tl_Awrk(i,j  ,k+1)*                      &
!^   &                            GRID(ng)%umask(i,j  )-                &
!^   &                         tl_Awrk(i,j-1,k+1)*
!^   &                            GRID(ng)%umask(i,j-1))
!^
              adfac=cff*ad_dAde(i,j,k2)
              ad_Awrk(i,j  ,k+1)=ad_Awrk(i,j  ,k+1)+                    &
     &                           GRID(ng)%umask(i,j  )*adfac
              ad_Awrk(i,j-1,k+1)=ad_Awrk(i,j-1,k+1)-                    &
     &                           GRID(ng)%umask(i,j-1)*adfac
              ad_dAde(i,j,k2)=0.0_r8
#  else
!^            tl_dAde(i,j,k2)=cff*(tl_Awrk(i,j  ,k+1)-                  &
!^   &                             tl_Awrk(i,j-1,k+1))
!^
              adfac=cff*ad_dAde(i,j,k2)
              ad_Awrk(i,j  ,k+1)=ad_Awrk(i,j  ,k+1)+adfac
              ad_Awrk(i,j-1,k+1)=ad_Awrk(i,j-1,k+1)-adfac
              ad_dAde(i,j,k2)=0.0_r8
#  endif
            END DO
          END DO
!
          DO j=Jstr,Jend
            DO i=IstrU-1,Iend
#  ifdef MASKING
!^            tl_dAdx(i,j,k2)=tl_dAdx(i,j,k2)*GRID(ng)%rmask(i,j)
!^
              ad_dAdx(i,j,k2)=ad_dAdx(i,j,k2)*GRID(ng)%rmask(i,j)
!^            tl_dAdx(i,j,k2)=GRID(ng)%pm(i,j)*                         &
!^   &                        (tl_Awrk(i+1,j,k+1)*                      &
!^   &                            GRID(ng)%umask(i+1,j)-                &
!^   &                         tl_Awrk(i  ,j,k+1)*                      &
!^   &                            GRID(ng)%umask(i  ,j))
!^
              adfac=GRID(ng)%pm(i,j)*ad_dAdx(i,j,k2)
              ad_Awrk(i  ,j,k+1)=ad_Awrk(i  ,j,k+1)-                    &
     &                           GRID(ng)%umask(i  ,j)*adfac
              ad_Awrk(i+1,j,k+1)=ad_Awrk(i+1,j,k+1)+                    &
     &                           GRID(ng)%umask(i+1,j)*adfac
              ad_dAdx(i,j,k2)=0.0_r8
#  else
!^            tl_dAdx(i,j,k2)=GRID(ng)%pm(i,j)*                         &
!^   &                        (tl_Awrk(i+1,j,k+1)-                      &
!^   &                         tl_Awrk(i  ,j,k+1))
!^
              adfac=GRID(ng)%pm(i,j)*ad_dAdx(i,j,k2)
              ad_Awrk(i  ,j,k+1)=ad_Awrk(i  ,j,k+1)-adfac
              ad_Awrk(i+1,j,k+1)=ad_Awrk(i+1,j,k+1)+adfac
              ad_dAdx(i,j,k2)=0.0_r8
#  endif
            END DO
          END DO
        END IF
!
!  Compute new storage recursive indices.
!
        kt=k2
        k2=k1
        k1=kt
      END DO K_LOOP

# else

!
!  Adjoint of diffusion along S-coordinates.
!  ========================================
!
!  Adjoint of compute horizontal K-Laplacian operator.
!
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!^          tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                              &
!^   &                     Hfac(i,j)*                                   &
!^   &                     (tl_FX(i,j)-tl_FX(i-1,j)+                    &
!^   &                      tl_FE(i,j+1)-tl_FE(i,j))
!^
            adfac=-Hfac(i,j)*ad_Awrk(i,j,k)
            ad_FE(i,j  )=ad_FE(i,j  )-adfac
            ad_FE(i,j+1)=ad_FE(i,j+1)+adfac
            ad_FX(i-1,j)=ad_FX(i-1,j)-adfac
            ad_FX(i  ,j)=ad_FX(i  ,j)+adfac
          END DO
        END DO
!
!  Adjoint of compute XI- and ETA-components adjoint diffusive flux.
!
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
#  ifdef MASKING
!^          tl_FE(i,j)=tl_FE(i,j)*GRID(ng)%pmask(i,j)
!^
            ad_FE(i,j)=ad_FE(i,j)*GRID(ng)%pmask(i,j)
#  endif
!^          tl_FE(i,j)=GRID(ng)%pnom_p(i,j)*                            &
!^   &                 0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                &
!^   &                          Khy(i-1,j-1)+Khy(i,j-1))*               &
!^   &                   (tl_Awrk(i,j,k)-tl_Awrk(i,j-1,k))
!^
            adfac=GRID(ng)%pnom_p(i,j)*                                 &
     &            0.25_r8*(Khy(i-1,j  )+Khy(i,j  )+                     &
     &                     Khy(i-1,j-1)+Khy(i,j-1))*ad_FE(i,j)
            ad_Awrk(i,j-1,k)=ad_Awrk(i,j-1,k)-adfac
            ad_Awrk(i,j  ,k)=ad_Awrk(i,j  ,k)+adfac
            ad_FE(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
!^          tl_FX(i,j)=GRID(ng)%pmon_r(i,j)*                            &
!^   &                 Khx(i,j)*(tl_Awrk(i+1,j,k)-tl_Awrk(i,j,k))
!^
            adfac=GRID(ng)%pmon_r(i,j)*Khx(i,j)*ad_FX(i,j)
            ad_Awrk(i  ,j,k)=ad_Awrk(i  ,j,k)-adfac
            ad_Awrk(i+1,j,k)=ad_Awrk(i+1,j,k)+adfac
            ad_FX(i,j)=0.0_r8
          END DO
        END DO
      END DO
# endif
!
!  Adjoint of set operator initial conditions.
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
!^          tl_Awrk(i,j,k)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k)
            ad_Awrk(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

# ifdef DISTRIBUTE
!
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
# endif
!
!^    CALL dabc_u3d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    tl_A)
!^
      CALL ad_dabc_u3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       ad_A)

# ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_u3d_ad
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of a
!  3D control variable at V-points.
!
      SUBROUTINE multiscale_Klap_v3d_tl (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Mlap, i, j, k, k1, k2, rec
      integer                           :: Istr, Iend, JstrV, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: cff, cff1, cff2, cff3, cff4
      real (r8)                         :: cffx, cffy

# ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: tl_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_FX
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS)   :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS)   :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: tl_dAde
# endif
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('v', 'v_northward')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%v_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%v_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
# ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
# else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
# endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Set operator initial conditions.
!
!^    CALL dabc_v3d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    A)
!^
      CALL dabc_v3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    tl_A)
# ifdef DISTRIBUTE
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
# endif
!
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
!^          Awrk(i,j,k)=A(i,j,k)
!^
            tl_Awrk(i,j,k)=tl_A(i,j,k)
          END DO
        END DO
      END DO

# ifdef GEOPOTENTIAL_HCONV
!
!  Diffusion along geopotential surfaces: Compute horizontal and
!  vertical gradients.  Notice the recursive blocking sequence.  The
!  vertical placement of the gradients is:
!
!        dZde_r,dAdx,dAde(:,:,k1) k     rho-points
!        dZde_r,dAdx,dAde(:,:,k2) k+1   rho-points
!                  dZdx_p(:,:,k1) k     psi-points
!                  dZdx_p(:,:,k2) k+1   psi-points
!                 FZ,dAdz(:,:,k1) k-1/2   W-points
!                 FZ,dAdz(:,:,k2) k+1/2   W-points
!
      k2=1
      K_LOOP : DO k=0,N(ng)
        k1=k2
        k2=3-k1
        IF (k.lt.N(ng)) THEN
          DO j=JstrV-1,Jend
            DO i=Istr,Iend+1
              cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%umask(i,j)
#  endif
              dZdx(i,j)=cff*(GRID(ng)%z_r(i  ,j,k+1)-                   &
     &                       GRID(ng)%z_r(i-1,j,k+1))
            END DO
          END DO
!
          DO j=JstrV-1,Jend+1
            DO i=Istr,Iend
              cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
              cff=cff*GRID(ng)%vmask(i,j)
#  endif
              dZde(i,j)=cff*(GRID(ng)%z_r(i,j  ,k+1)-                   &
     &                       GRID(ng)%z_r(i,j-1,k+1))
            END DO
          END DO
!
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(GRID(ng)%pm(i-1,j-1)+GRID(ng)%pm(i-1,j)+     &
     &                     GRID(ng)%pm(i  ,j-1)+GRID(ng)%pm(i  ,j))
#  ifdef MASKING
!^            dAdx(i,j,k2)=cff*                                         &
!^   &                     (Awrk(i  ,j,k+1)*GRID(ng)%vmask(i  ,j)-      &
!^   &                      Awrk(i-1,j,k+1)*GRID(ng)%vmask(i-1,j))
!^
              tl_dAdx(i,j,k2)=cff*                                      &
     &                        (tl_Awrk(i  ,j,k+1)*                      &
     &                            GRID(ng)%vmask(i  ,j)-                &
     &                         tl_Awrk(i-1,j,k+1)*
     7                            GRID(ng)%vmask(i-1,j))
!^            dAdx(i,j,k2)=dAdx(i,j,k2)*GRID(ng)%pmask(i,j)
!^
              tl_dAdx(i,j,k2)=tl_dAdx(i,j,k2)*GRID(ng)%pmask(i,j)
#  else
!^            dAdx(i,j,k2)=cff*(Awrk(i  ,j,k+1)-                        &
!^   &                          Awrk(i-1,j,k+1))
!^
              tl_dAdx(i,j,k2)=cff*(tl_Awrk(i  ,j,k+1)-                  &
     &                             tl_Awrk(i-1,j,k+1))
#  endif
              dZdx_p(i,j,k2)=0.5_r8*(dZdx(i,j-1)+                       &
     &                               dZdx(i,j  ))
            END DO
          END DO
!
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
#  ifdef MASKING
!^            dAde(i,j,k2)=GRID(ng)%pn(i,j)*                            &
!^   &                     (Awrk(i,j+1,k+1)*GRID(ng)%vmask(i,j+1)-      &
!^   &                      Awrk(i,j  ,k+1)*GRID(ng)%vmask(i,j  ))
!^
              tl_dAde(i,j,k2)=GRID(ng)%pn(i,j)*                         &
     &                        (tl_Awrk(i,j+1,k+1)*                      &
     &                            GRID(ng)%vmask(i,j+1)-                &
     &                         tl_Awrk(i,j  ,k+1)*
     &                            GRID(ng)%vmask(i,j  ))
!^            dAde(i,j,k2)=dAde(i,j,k2)*GRID(ng)%rmask(i,j)
!^
              tl_dAde(i,j,k2)=tl_dAde(i,j,k2)*GRID(ng)%rmask(i,j)
#  else
!^            dAde(i,j,k2)=GRID(ng)%pn(i,j)*(Awrk(i,j+1,k+1)-           &
!^   &                                       Awrk(i,j  ,k+1))
!^
              tl_dAde(i,j,k2)=GRID(ng)%pn(i,j)*(tl_Awrk(i,j+1,k+1)-     &
     &                                          tl_Awrk(i,j  ,k+1))
#  endif
              dZde_r(i,j,k2)=0.5_r8*(dZde(i,j  )+                       &
     &                               dZde(i,j+1))
            END DO
          END DO
        END IF
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
!^            dAdz(i,j,k2)=0.0_r8
!^
              tl_dAdz(i,j,k2)=0.0_r8
!^            FZ(i,j,k2)=0.0_r8
!^
              tl_FZ(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(GRID(ng)%z_r(i,j,k+1)-GRID(ng)%z_r(i,j,k))
!^            dAdz(i,j,k2)=cff*(Awrk(i,j,k+1)-                          &
!^   &                          Awrk(i,j,k))
!^
              tl_dAdz(i,j,k2)=cff*(tl_Awrk(i,j,k+1)-                    &
     &                             tl_Awrk(i,j,k))
#  ifdef MASKING
!^            dAdz(i,j,k2)=dAdz(i,j,k2)*GRID(ng)%vmask(i,j)
!^
              tl_dAdz(i,j,k2)=tl_dAdz(i,j,k2)*GRID(ng)%vmask(i,j)
#  endif
            END DO
          END DO
        END IF
!
!  Compute components of the rotated A flux along geopotential
!  surfaces.
!
        IF (k.gt.0) THEN
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(Khx(i-1,j-1)+Khx(i-1,j)+                     &
     &                     Khx(i  ,j-1)+Khx(i  ,j))*GRID(ng)%on_p(i,j)
              cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
!^            FX(i,j)=cff*                                              &
!^   &                (dAdx(i,j,k1)-                                    &
!^   &                 0.5_r8*(cff1*(dAdz(i-1,j,k1)+                    &
!^   &                               dAdz(i  ,j,k2))+                   &
!^   &                         cff2*(dAdz(i-1,j,k2)+                    &
!^   &                               dAdz(i  ,j,k1))))
!^
              tl_FX(i,j)=cff*                                           &
     &                   (tl_dAdx(i,j,k1)-                              &
     &                    0.5_r8*(cff1*(tl_dAdz(i-1,j,k1)+              &
     &                                  tl_dAdz(i  ,j,k2))+             &
     &                            cff2*(tl_dAdz(i-1,j,k2)+              &
     &                                  tl_dAdz(i  ,j,k1))))
            END DO
          END DO
!
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
              cff=Khy(i,j)*GRID(ng)%om_r(i,j)
              cff1=MIN(dZde_r(i,j,k1),0.0_r8)
              cff2=MAX(dZde_r(i,j,k1),0.0_r8)
!^            FE(i,j)=cff*                                              &
!^   &                (dAde(i,j,k1)-                                    &
!^   &                 0.5_r8*(cff1*(dAdz(i,j  ,k1)+                    &
!^   &                               dAdz(i,j+1,k2))+                   &
!^   &                         cff2*(dAdz(i,j  ,k2)+                    &
!^   &                               dAdz(i,j+1,k1))))
!^
              tl_FE(i,j)=cff*                                           &
     &                   (tl_dAde(i,j,k1)-                              &
     &                    0.5_r8*(cff1*(tl_dAdz(i,j  ,k1)+              &
     &                                  tl_dAdz(i,j+1,k2))+             &
     &                            cff2*(tl_dAdz(i,j  ,k2)+              &
     &                                  tl_dAdz(i,j+1,k1))))
            END DO
          END DO
!
          IF (k.lt.N(ng)) THEN
            DO j=JstrV,Jend
              DO i=Istr,Iend
                cff=0.5_r8*(Khx(i,j-1)+Khx(i,j))
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
!^              FZ(i,j,k2)=cff*                                         &
!^   &                     (cff1*(cff1*dAdz(i,j,k2)-dAdx(i  ,j,k1))+    &
!^   &                      cff2*(cff2*dAdz(i,j,k2)-dAdx(i+1,j,k2))+    &
!^   &                      cff3*(cff3*dAdz(i,j,k2)-dAdx(i  ,j,k2))+    &
!^   &                      cff4*(cff4*dAdz(i,j,k2)-dAdx(i+1,j,k1)))
!^
                tl_FZ(i,j,k2)=cff*                                      &
     &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i  ,j,k1))+                &
     &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i+1,j,k2))+                &
     &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i  ,j,k2))+                &
     &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
     &                               tl_dAdx(i+1,j,k1)))
!
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
                cff=0.5_r8*(Khy(i,j-1)+Khy(i,j))
!^              FZ(i,j,k2)=FZ(i,j,k2)+                                  &
!^   &                     cff*                                         &
!^   &                     (cff1*(cff1*dAdz(i,j,k2)-dAde(i,j-1,k1))+    &
!^   &                      cff2*(cff2*dAdz(i,j,k2)-dAde(i,j  ,k2))+    &
!^   &                      cff3*(cff3*dAdz(i,j,k2)-dAde(i,j-1,k2))+    &
!^   &                      cff4*(cff4*dAdz(i,j,k2)-dAde(i,j  ,k1)))
!^
                tl_FZ(i,j,k2)=tl_FZ(i,j,k2)+                            &
     &                        cff*                                      &
     &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j-1,k1))+                &
     &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j  ,k2))+                &
     &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j-1,k2))+                &
     &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
     &                               tl_dAde(i,j  ,k1)))
              END DO
            END DO
          END IF
!
!  Compute horizontal K-Laplacian operator.
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
!^            Awrk(i,j,k)=Awrk(i,j,k)-                                  &
!^   &                    Hfac(i,j)*                                    &
!^   &                    (FX(i+1,j)-FX(i,j  )+                         &
!^   &                     FE(i  ,j)-FE(i,j-1))-                        &
!^   &                    2.0_r8*(FZ(i,j,k2)-FZ(i,j,k1))/               &
!^   &                           (GRID(ng)%Hz(i,j  ,k)+                 &
!^   &                            GRID(ng)%Hz(i,j-1,k))
!^
              tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                            &
     &                       Hfac(i,j)*                                 &
     &                       (tl_FX(i+1,j)-tl_FX(i,j  )+                &
     &                        tl_FE(i  ,j)-tl_FE(i,j-1))-               &
     &                       2.0_r8*(tl_FZ(i,j,k2)-tl_FZ(i,j,k1))/      &
     &                              (GRID(ng)%Hz(i,j  ,k)+              &
     &                               GRID(ng)%Hz(i,j-1,k))
            END DO
          END DO
        END IF
      END DO K_LOOP

# else

!
!  Compute XI- and ETA-components of diffusive flux.
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
!^          FX(i,j)=GRID(ng)%pmon_p(i,j)*                               &
!^   &              0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                   &
!^   &                       Khx(i-1,j-1)+Khx(i,j-1))*                  &
!^   &              (Awrk(i,j,k)-Awrk(i-1,j,k))
!^
            tl_FX(i,j)=GRID(ng)%pmon_p(i,j)*                            &
     &                 0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                &
     &                          Khx(i-1,j-1)+Khx(i,j-1))*               &
     &                 (tl_Awrk(i,j,k)-tl_Awrk(i-1,j,k))
#  ifdef MASKING
!^          FX(i,j)=FX(i,j)*GRID(ng)%pmask(i,j)
!^
            tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%pmask(i,j)
#  endif
          END DO
        END DO
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
!^          FE(i,j)=GRID(ng)%pnom_r(i,j)*                               &
!^   *              Khy(i,j)*(Awrk(i,j+1,k)-Awrk(i,j,k))
!^
              tl_FE(i,j)=GRID(ng)%pnom_r(i,j)*                          &
     &                   Khy(i,j)*(tl_Awrk(i,j+1,k)-tl_Awrk(i,j,k))
          END DO
        END DO
!
!  Compute horizontal K-Laplacian operator.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          Awrk(i,j,k)=Awrk(i,j,k)-                                    &
!^   &                  Hfac(i,j)*                                      &
!^   &                  FX(i+1,j)-FX(i,j)+                              &
!^   &                  FE(i,j)-FE(i,j-1))
!^
            tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                              &
     &                     Hfac(i,j)*                                   &
     &                     (tl_FX(i+1,j)-tl_FX(i,j)+                    &
     &                      tl_FE(i,j)-tl_FE(i,j-1))
          END DO
        END DO
      END DO
# endif
!
!  Load K-Laplacian solution.
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          A(i,j,k)=Awrk(i,j,k)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k)
          END DO
        END DO
      END DO

# ifdef DISTRIBUTE
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
# endif

# ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_v3d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of a
!  3D control variable at V-points.
!
      SUBROUTINE multiscale_Klap_v3d_ad (self, ng, tile, model,         &
     &                                   ifield, ctype, ms, Lweak,      &
     &                                   LBi, UBi, LBj, UBj,            &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Mlap, i, j, k, k1, k2, rec
      integer                           :: Istr, Iend, JstrV, Jend
      integer                           :: is, ie, js, je
!
      real (r8)                         :: cff, cff1, cff2, cff3, cff4
      real (r8)                         :: adfac, cffx, cffy

# ifdef NONUNIFORM_SCALES
!
      real (r8), pointer                :: BscaleX(:,:) => NULL()
      real (r8), pointer                :: BscaleY(:,:) => NULL()
# endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: ad_Awrk
!
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hfac
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Khy
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_FX
# ifdef GEOPOTENTIAL_HCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS)   :: dZdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS)   :: dZde
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZdx_p
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: dZde_r
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_FZ
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAdz
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAdx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,2) :: ad_dAde
# endif
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      IF (Lweak) THEN
        rec=2                        ! weak constraint correlations
      ELSE
        rec=1                        ! strong constraint correlations
      END IF
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
# ifdef GEOPOTENTIAL_HCONV
      ad_FZ=0.0_r8
      ad_dAdz=0.0_r8
      ad_dAdx=0.0_r8
      ad_dAde=0.0_r8
# endif
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('v', 'v_northward')
# ifdef NONUNIFORM_SCALES
          BscaleX(is:,js:) => self%v_Bcorr(is:ie,js:je,1,ms)
          BscaleY(is:,js:) => self%v_Bcorr(is:ie,js:je,2,ms)
# endif
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
          Hfac(i,j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))*      &
     &                      (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
        END DO
      END DO
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=2, kappa=D*D/(2*(M-2)), where D is the Daley length scale and
!  d is the space dimension.
!
      DO j=JstrV-1,Jend+1
        DO i=Istr-1,Iend+1
# ifdef NONUNIFORM_SCALES
          cffx=BscaleX(i,j)*BscaleX(i,j)   ! spatially varying
          cffy=BscaleY(i,j)*BscaleY(i,j)
# else
          cffx=HdecayX(rec,ifield,ms,ng)*HdecayX(rec,ifield,ms,ng)
          cffy=HdecayY(rec,ifield,ms,ng)*HdecayY(rec,ifield,ms,ng)
# endif
          Khx(i,j)=0.5_r8*cffx/REAL(Mlap-2,r8)
          Khy(i,j)=0.5_r8*cffy/REAL(Mlap-2,r8)
        END DO
      END DO
!
!  Adjoint of load K-Laplacian solution.
!
# ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_Awrk)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_Awrk)
!
# endif
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          tl_A(i,j,k)=tl_Awrk(i,j,k)
!^
            ad_Awrk(i,j,k)=ad_Awrk(i,j,k)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

# ifdef GEOPOTENTIAL_HCONV
!
!  Adjoint of diffusion along geopotential surfaces.
!  ================================================
!
!  Compute horizontal and vertical gradients. Notice the recursive
!  blocking sequence. The vertical placement of the gradients is:
!
!        dZde_r,dAdx,dAde(:,:,k1) k     rho-points
!        dZde_r,dAdx,dAde(:,:,k2) k+1   rho-points
!                  dZdx_p(:,:,k1) k     psi-points
!                  dZdx_p(:,:,k2) k+1   psi-points
!                 FZ,dAdz(:,:,k1) k-1/2   W-points
!                 FZ,dAdz(:,:,k2) k+1/2   W-points
!
!  Compute adjoint of starting values of k1 and k2.
!
      k1=2
      k2=1
      DO k=0,N(ng)
!
!  Note: The following code is equivalent to
!
!       kt=k1
!       k1=k2
!       k2=kt
!
!  We use the adjoint of above code.
!
        k1=k2
        k2=3-k1
      END DO
!
      K_LOOP : DO k=N(ng),0,-1

!  Compute required BASIC STATE fields. Need to look forward in
!  recursive kk index.
!
        k2b=1
        DO kk=0,k
          k1b=k2b
          k2b=3-k1b
!
!  Compute components of the rotated tracer flux (A m3/s) along
!  geopotential surfaces (required BASIC STATE fields).
!
          IF (kk.lt.N(ng)) THEN
            DO j=JstrV-1,Jend
              DO i=Istr,Iend+1
                cff=0.5_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))
#  ifdef MASKING
                cff=cff*GRID(ng)%umask(i,j)
#  endif
                dZdx(i,j)=cff*(GRID(ng)%z_r(i  ,j,kk+1)-                &
     &                         GRID(ng)%z_r(i-1,j,kk+1))
              END DO
            END DO
            DO j=JstrV,Jend
              DO i=Istr,Iend+1
                dZdx_p(i,j,k2b)=0.5_r8*(dZdx(i,j-1)+                    &
     &                                  dZdx(i,j  ))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=JstrV,Jend
                DO i=Istr,Iend+1
                  dZdx_p(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
!
            DO j=JstrV-1,Jend+1
              DO i=Istr,Iend
                cff=0.5_r8*(GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
#  ifdef MASKING
                cff=cff*GRID(ng)%vmask(i,j)
#  endif
                dZde(i,j)=cff*(GRID(ng)%z_r(i,j  ,kk+1)-                &
     &                         GRID(ng)%z_r(i,j-1,kk+1))
              END DO
            END DO
            DO j=JstrV-1,Jend
              DO i=Istr,Iend
                dZde_r(i,j,k2b)=0.5_r8*(dZde(i,j  )+                    &
     &                                  dZde(i,j+1))
              END DO
            END DO
            IF (kk.eq.0) THEN
              DO j=JstrV-1,Jend
                DO i=Istr,Iend
                  dZde_r(i,j,k1b)=0.0_r8
                END DO
              END DO
            END IF
          END IF
        END DO
!
        IF (k.gt.0) THEN
!
!  Adjoint of compute geopotential K-Laplacian operator.
!
          DO j=JstrV,Jend
            DO i=Istr,Iend
!^            tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                            &
!^   &                       Hfac(i,j)*                                 &
!^   &                       (tl_FX(i+1,j)-tl_FX(i,j  )+                &
!^   &                        tl_FE(i  ,j)-tl_FE(i,j-1))-               &
!^   &                       2.0_r8*(tl_FZ(i,j,k2)-tl_FZ(i,j,k1))/
!^   &                       (GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i,j-1,k))
!^
              adfac1=-Hfac(i,j)*ad_Awrk(i,j,k)
              adfac2=-2.0_r8*ad_Awrk(i,j,k)/                            &
     &               (GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i,j-1,k))
              ad_FE(i,j-1)=ad_FE(i,j-1)-adfac1
              ad_FE(i,j  )=ad_FE(i,  j)+adfac1
              ad_FX(i  ,j)=ad_FX(i  ,j)-adfac1
              ad_FX(i+1,j)=ad_FX(i+1,j)+adfac1
              ad_FZ(i,j,k1)=ad_FZ(i,j,k1)-adfac2
              ad_FZ(i,j,k2)=ad_FZ(i,j,k2)+adfac2
            END DO
          END DO
!
!  Adjoint of compute components of rotated A flux (A m3/s) along
!  geopotential surfaces.
!
          IF (k.lt.N(ng)) THEN
            DO j=JstrV,Jend
              DO i=Istr,Iend
                cff=0.5_r8*(Khy(i,j-1)+Khy(i,j))
                cff1=MIN(dZde_r(i,j-1,k1),0.0_r8)
                cff2=MIN(dZde_r(i,j  ,k2),0.0_r8)
                cff3=MAX(dZde_r(i,j-1,k2),0.0_r8)
                cff4=MAX(dZde_r(i,j  ,k1),0.0_r8)
!^              tl_FZ(i,j,k2)=tl_FZ(i,j,k2)+                            &
!^   &                        cff*                                      &
!^   &                        (cff1*(cff1*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j-1,k1))+                &
!^   &                         cff2*(cff2*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j  ,k2))+                &
!^   &                         cff3*(cff3*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j-1,k2))+                &
!^   &                         cff4*(cff4*tl_dAdz(i,j,k2)-              &
!^   &                               tl_dAde(i,j  ,k1)))
!^
                adfac=cff*ad_FZ(i,j,k2)
                ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)+                        &
     &                          (cff1*cff1+                             &
     &                           cff2*cff2+                             &
     &                           cff3*cff3+                             &
     &                           cff4*cff4)*adfac
                ad_dAde(i,j-1,k1)=ad_dAde(i,j-1,k1)-cff1*adfac
                ad_dAde(i,j  ,k2)=ad_dAde(i,j  ,k2)-cff2*adfac
                ad_dAde(i,j-1,k2)=ad_dAde(i,j-1,k2)-cff3*adfac
                ad_dAde(i,j  ,k1)=ad_dade(i,j  ,k1)-cff4*adfac
!
                cff1=MIN(dZdx_p(i  ,j,k1),0.0_r8)
                cff2=MIN(dZdx_p(i+1,j,k2),0.0_r8)
                cff3=MAX(dZdx_p(i  ,j,k2),0.0_r8)
                cff4=MAX(dZdx_p(i+1,j,k1),0.0_r8)
                cff=0.5_r8*(Khx(i,j-1)+Khx(i,j))
!^              tl_FZ(i,j,k2)=cff*                                      &
!^   &                          (cff1*(cff1*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i  ,j,k1))+              &
!^   &                           cff2*(cff2*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i+1,j,k2))+              &
!^   &                           cff3*(cff3*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i  ,j,k2))+              &
!^   &                           cff4*(cff4*tl_dAdz(i,j,k2)-            &
!^   &                                 tl_dAdx(i+1,j,k1)))
!^
                ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)+                        &
     &                          (cff1*cff1+                             &
     &                           cff2*cff2+                             &
     &                           cff3*cff3+                             &
     &                           cff4*cff4)*adfac
                ad_dAdx(i  ,j,k1)=ad_dAdx(i  ,j,k1)-cff1*adfac
                ad_dAdx(i+1,j,k2)=ad_dAdx(i+1,j,k2)-cff2*adfac
                ad_dAdx(i  ,j,k2)=ad_dAdx(i  ,j,k2)-cff3*adfac
                ad_dAdx(i+1,j,k1)=ad_dAdx(i+1,j,k1)-cff4*adfac
                ad_FZ(i,j,k2)=0.0_r8
              END DO
            END DO
          END IF
!
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
              cff=Khy(i,j)*GRID(ng)%om_r(i,j)
              cff1=MIN(dZde_r(i,j,k1),0.0_r8)
              cff2=MAX(dZde_r(i,j,k1),0.0_r8)
!^            tl_FE(i,j)=cff*                                           &
!^   &                   (tl_dAde(i,j,k1)-                              &
!^   &                    0.5_r8*(cff1*(tl_dAdz(i,j  ,k1)+              &
!^   &                                  tl_dAdz(i,j+1,k2))+             &
!^   &                            cff2*(tl_dAdz(i,j  ,k2)+              &
!^   &                                  tl_dAdz(i,j+1,k1))))
!^
              adfac=cff*ad_FE(i,j)
              adfac1=adfac*0.5_r8*cff1
              adfac2=adfac*0.5_r8*cff2
              ad_dAde(i,j,k1)=ad_dAde(i,j,k1)+adfac
              ad_dAdz(i,j  ,k1)=ad_dAdz(i,j  ,k1)-adfac1
              ad_dAdz(i,j+1,k2)=ad_dAdz(i,j+1,k2)-adfac1
              ad_dAdz(i,j  ,k2)=ad_dAdz(i,j  ,k2)-adfac2
              ad_dAdz(i,j+1,k1)=ad_dAdz(i,j+1,k1)-adfac2
              ad_FE(i,j)=0.0_r8
            END DO
          END DO
!
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(Khx(i-1,j-1)+Khx(i-1,j)+                     &
     &                     Khx(i  ,j-1)+Khx(i  ,j))*on_p(i,j)
              cff1=MIN(dZdx_p(i,j,k1),0.0_r8)
              cff2=MAX(dZdx_p(i,j,k1),0.0_r8)
!^            tl_FX(i,j)=cff*                                           &
!^   &                   (tl_dAdx(i,j,k1)-                              &
!^   &                    0.5_r8*(cff1*(tl_dAdz(i-1,j,k1)+              &
!^   &                                  tl_dAdz(i  ,j,k2))+             &
!^   &                            cff2*(tl_dAdz(i-1,j,k2)+              &
!^   &                                  tl_dAdz(i  ,j,k1))))
!^
              adfac=cff*ad_FX(i,j)
              adfac1=adfac*0.5_r8*cff1
              adfac2=adfac*0.5_r8*cff2
              ad_dAdx(i,j,k1)=ad_dAdx(i,j,k1)+adfac
              ad_dAdz(i-1,j,k1)=ad_dAdz(i-1,j,k1)-adfac1
              ad_dAdz(i  ,j,k2)=ad_dAdz(i  ,j,k2)-adfac1
              ad_dAdz(i-1,j,k2)=ad_dAdz(i-1,j,k2)-adfac2
              ad_dAdz(i  ,j,k1)=ad_dAdz(i  ,j,k1)-adfac2
              ad_FX(i,j)=0.0_r8
            END DO
          END DO
        END IF
!
        IF ((k.eq.0).or.(k.eq.N(ng))) THEN
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
!^            tl_FZ(i,j,k2)=0.0_r8
!^
              ad_FZ(i,j,k2)=0.0_r8
!^            tl_dAdz(i,j,k2)=0.0_r8
!^
              ad_dAdz(i,j,k2)=0.0_r8
            END DO
          END DO
        ELSE
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              cff=1.0_r8/(GRID(ng)%z_r(i,j,k+1)-GRID(ng)%z_r(i,j,k))
#  ifdef MASKING
!^            tl_dAdz(i,j,k2)=tl_dAdz(i,j,k2)*GRID(ng)%vmask(i,j)
!^
              ad_dAdz(i,j,k2)=ad_dAdz(i,j,k2)*GRID(ng)%vmask(i,j)
#  endif
!^            tl_dAdz(i,j,k2)=cff*(tl_Awrk(i,j,k+1)-                    &
!^   &                             tl_Awrk(i,j,k))
!^
              adfac=cff*ad_dAdz(i,j,k2)
              ad_Awrk(i,j,k)=ad_Awrk(i,j,k)-adfac
              ad_Awrk(i,j,k+1)=ad_Awrk(i,j,k+1)+adfac
              ad_dAdz(i,j,k2)=0.0_r8
            END DO
          END DO
        END IF
!
        IF (k.lt.N(ng)) THEN
          DO j=JstrV-1,Jend
            DO i=Istr,Iend
#  ifdef MASKING
!^            tl_dAde(i,j,k2)=tl_dAde(i,j,k2)*GRID(ng)%rmask(i,j)
!^
              ad_dAde(i,j,k2)=ad_dAde(i,j,k2)*GRID(ng)%rmask(i,j)
!^            tl_dAde(i,j,k2)=GRID(ng)%pn(i,j)*                         &
!^   &                        (tl_Awrk(i,j+1,k+1)*                      &
!^   &                            GRID(ng)%vmask(i,j+1)-                &
!^   &                         tl_Awrk(i,j  ,k+1)*                      &
!^   &                            GRID(ng)%vmask(i,j  ))
!^
              adfac=GRID(ng)%pn(i,j)*ad_dAde(i,j,k2)
              ad_Awrk(i,j  ,k+1)=ad_Awrk(i,j  ,k+1)-                    &
     &                           GRID(ng)%vmask(i,j  )*adfac
              ad_Awrk(i,j+1,k+1)=ad_Awrk(i,j+1,k+1)+                    &
     &                           GRID(ng)%vmask(i,j+1)*adfac
              ad_dAde(i,j,k2)=0.0_r8
#  else
!^            tl_dAde(i,j,k2)=GRID(ng)%pn(i,j)*                         &
!^   &                        (tl_Awrk(i,j+1,k+1)-                      &
!^   &                         tl_Awrk(i,j  ,k+1))
!^
              adfac=GRID(ng)%pn(i,j)*ad_dAde(i,j,k2)
              ad_Awrk(i,j  ,k+1)=ad_Awrk(i,j  ,k+1)-adfac
              ad_Awrk(i,j+1,k+1)=ad_Awrk(i,j+1,k+1)+adfac
              ad_dAde(i,j,k2)=0.0_r8
#  endif
            END DO
          END DO
!
          DO j=JstrV,Jend
            DO i=Istr,Iend+1
              cff=0.25_r8*(GRID(ng)%pm(i-1,j-1)+GRID(ng)%pm(i-1,j)+     &
     &                     GRID(ng)%pm(i  ,j-1)+GRID(ng)%pm(i  ,j))
#  ifdef MASKING
!^            tl_dAdx(i,j,k2)=tl_dAdx(i,j,k2)*GRID(ng)%pmask(i,j)
!^
              ad_dAdx(i,j,k2)=ad_dAdx(i,j,k2)*GRID(ng)%pmask(i,j)
!^            tl_dAdx(i,j,k2)=cff*                                      &
!^   &                        (tl_Awrk(i  ,j,k+1)*                      &
!^   &                            GRID(ng)%vmask(i  ,j)-                &
!^   &                         tl_Awrk(i-1,j,k+1)*                      &
!^   &                            GRID(ng)%vmask(i-1,j))
!^
              adfac=cff*ad_dAdx(i,j,k2)
              ad_Awrk(i-1,j,k+1)=ad_Awrk(i-1,j,k+1)-                    &
     &                           GRID(ng)%vmask(i-1,j)*adfac
              ad_Awrk(i  ,j,k+1)=ad_Awrk(i  ,j,k+1)+                    &
     &                           GRID(ng)%vmask(i  ,j)*adfac
              ad_dAdx(i,j,k2)=0.0_r8
#  else
!^            tl_dAdx(i,j,k2)=cff*(tl_Awrk(i  ,j,k+1)-                  &
!^   &                             tl_Awrk(i-1,j,k+1))
!^
              adfac=cff*ad_dAdx(i,j,k2)
              ad_Awrk(i-1,j,k+1)=ad_Awrk(i-1,j,k+1)-adfac
              ad_Awrk(i  ,j,k+1)=ad_Awrk(i  ,j,k+1)+adfac
              ad_dAdx(i,j,k2)=0.0_r8
#  endif
            END DO
          END DO
        END IF
!
!  Compute new storage recursive indices.
!
        kt=k2
        k2=k1
        k1=kt
      END DO K_LOOP

# else

!
!  Adjoint of diffusion along S-coordinates.
!  ========================================
!
!  Adjoint of compute K-Laplacian operator.
!
      DO k=1,N(ng)
        DO j=JstrV,Jend
          DO i=Istr,Iend
!^          tl_Awrk(i,j,k)=tl_Awrk(i,j,k)-                              &
!^   &                     Hfac(i,j)*                                   &
!^   &                     (tl_FX(i+1,j)-tl_FX(i,j)+                    &
!^   &                      tl_FE(i,j)-tl_FE(i,j-1))
!^
            adfac=-Hfac(i,j)*ad_Awrk(i,j,k)
            ad_FE(i,j-1)=ad_FE(i,j-1)-adfac
            ad_FE(i,j  )=ad_FE(i,j  )+adfac
            ad_FX(i  ,j)=ad_FX(i  ,j)-adfac
            ad_FX(i+1,j)=ad_FX(i+1,j)+adfac
          END DO
        END DO
!
!  Adjoint of compute XI- and ETA-components of diffusive flux.
!
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
!^          tl_FE(i,j)=GRID(ng)%pnom_r(i,j)*                            &
!^   &                 Khy(i,j)*(tl_Awrk(i,j+1,k)-tl_Awrk(i,j,k))
!^
            adfac=GRID(ng)%pnom_r(i,j)*Khy(i,j)*ad_FE(i,j)
            ad_Awrk(i,j  ,k)=ad_Awrk(i,j  ,k)-adfac
            ad_Awrk(i,j+1,k)=ad_Awrk(i,j+1,k)+adfac
            ad_FE(i,j)=0.0_r8
          END DO
        END DO
!
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
#  ifdef MASKING
!^          tl_FX(i,j)=tl_FX(i,j)*GRID(ng)%pmask(i,j)
!^
            ad_FX(i,j)=ad_FX(i,j)*GRID(ng)%pmask(i,j)
#  endif
!^          tl_FX(i,j)=GRID(ng)%pmon_p(i,j)*                            &
!^   &                 0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                &
!^   &                          Khx(i-1,j-1)+Khx(i,j-1))*               &
!^   &                 (tl_Awrk(i,j,k)-tl_Awrk(i-1,j,k))
!^
            adfac=GRID(ng)%pmon_p(i,j)*                                 &
     &            0.25_r8*(Khx(i-1,j  )+Khx(i,j  )+                     &
     &                     Khx(i-1,j-1)+Khx(i,j-1))*ad_FX(i,j)
            ad_Awrk(i-1,j,k)=ad_Awrk(i-1,j,k)-adfac
            ad_Awrk(i  ,j,k)=ad_Awrk(i  ,j,k)+adfac
            ad_FX(i,j)=0.0_r8
          END DO
        END DO
      END DO
# endif
!
!  Adjoint of set operator initial conditions.
!
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
!^          tl_Awrk(i,j,k)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k)
            ad_Awrk(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

# ifdef DISTRIBUTE
!
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
# endif
!
!^    CALL dabc_v3d_tile (ng, tile,                                     &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    tl_A)
!^
      CALL ad_dabc_v3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       ad_A)

# ifdef NONUNIFORM_SCALES
!
!  Nullify local pointers.
!
      nullify (BscaleX)
      nullify (BscaleY)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_v3d_ad
#endif /* SOLVE3D */

#ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of
!  lateral boundary conditions for a 2D control variable.
!
      SUBROUTINE multiscale_Klap_b1d_tl (self, ng, tile, model,         &
     &                                   ifield, ibry, ctype, ms,       &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Mlap, i, j
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
!
      real (r8)                         :: cffx, cffy
!
      real (r8), dimension(LBij:UBij)   :: tl_Awrk
!
      real (r8), dimension(LBij:UBij)   :: Hfac
      real (r8), dimension(IminS:ImaxS) :: Khx
      real (r8), dimension(JminS:JmaxS) :: Khy
      real (r8), dimension(JminS:JmaxS) :: tl_FE
      real (r8), dimension(IminS:ImaxS) :: tl_FX
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
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
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta', 'ubar', 'ubar_eastward', 'vbar', 'vbar_northward')
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jstr-1,Jend+1
                 Hfac(j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Istr-1,Iend+1
                 Hfac(i)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jstr-1,Jend+1
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=IstrU-1,Iend+1
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=JstrV-1,Jend+1
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Istr-1,Iend+1
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=1, kappa=D*D/(2*M-d-2), where D is the Daley length scale
!  and d is the space dimension.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO j=JminS,JmaxS
            cffy=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khy(j)=cffy/REAL(2*Mlap-3,r8)
          END DO
        CASE (isouth, inorth)
          DO i=IminS,ImaxS
            cffx=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khx(i)=cffx/REAL(2*Mlap-3,r8)
          END DO
      END SELECT
!
!  Set operator initial conditions.
!
      SELECT CASE (ctype)
        CASE (r2dvar)
!^        CALL bc_r2d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij,                             &
!^   &                          A)
!^
          CALL bc_r2d_bry_tile (ng, tile, ibry,                         &
     &                          LBij, UBij,                             &
     &                          tl_A)
        CASE (u2dvar)
!^        CALL bc_u2d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij,                             &
!^   &                          A)
!^
          CALL bc_u2d_bry_tile (ng, tile, ibry,                         &
     &                          LBij, UBij,                             &
     &                          tl_A)
        CASE (v2dvar)
!^        CALL bc_u2d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij,                             &
!^   &                          A)
!^
          CALL bc_v2d_bry_tile (ng, tile, ibry,                         &
     &                          LBij, UBij,                             &
     &                          tl_A)
      END SELECT

# ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij,                               &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        A)
!^
      CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
     &                        LBij, UBij,                               &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        tl_A)
# endif
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin-1,Jmax+1
!^            Awrk(j)=A(j)
!^
              tl_Awrk(j)=tl_A(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin-1,Imax+1
!^            Awrk(i)=A(i)
!^
              tl_Awrk(i)=tl_A(i)
            END DO
        END SELECT
      END IF
!
!  Compute XI- or ETA-components of diffusive flux.
!
      SELECT CASE (ctype)
!
        CASE (r2dvar)                           ! 2D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,r2dvar)
              DO j=Jstr,Jend+1
!^              FE(j)=GRID(ng)%pnom_v(i,j)*                             &
!^   &                0.5_r8*(Khy(j-1)+Khy(j))*                         &
!^   &                (Awrk(j)-Awrk(j-1))
!^
                tl_FE(j)=GRID(ng)%pnom_v(i,j)*                          &
     &                   0.5_r8*(Khy(j-1)+Khy(j))*                      &
     &                   (tl_Awrk(j)-tl_Awrk(j-1))
# ifdef MASKING
!^              FE(j)=FE(j)*GRID(ng)%vmask(i,j)
!^
                tl_FE(j)=tl_FE(j)*GRID(ng)%vmask(i,j)
# endif
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,r2dvar)
              DO i=Istr,Iend+1
!^              FX(i)=GRID(ng)%pmon_u(i,j)*                             &
!^                    0.5_r8*(Khx(i-1)+Khx(i))*                         &
!^   &                (Awrk(i)-Awrk(i-1))
!^
                tl_FX(i)=GRID(ng)%pmon_u(i,j)*                          &
     &                   0.5_r8*(Khx(i-1)+Khx(i))*                      &
     &                   (tl_Awrk(i)-tl_Awrk(i-1))
# ifdef MASKING
!^              FX(i)=FX(i)*GRID(ng)%umask(i,j)
!^
                tl_FX(i)=tl_FX(i)*GRID(ng)%umask(i,j)
# endif
              END DO
            END IF
          END IF
!
        CASE (u2dvar)                           ! 2D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,u2dvar)
              DO j=Jstr,Jend+1
!^              FE(j)=GRID(ng)%pnom_p(i,j)*                             &
!^   &                0.5_r8*(Khy(j)+Khy(j-1))*                         &
!^   &                (Awrk(j)-Awrk(j-1))
!^
                tl_FE(j)=GRID(ng)%pnom_p(i,j)*                          &
     &                   0.5_r8*(Khy(j)+Khy(j-1))*                      &
     &                   (tl_Awrk(j)-tl_Awrk(j-1))
# ifdef MASKING
!^              FE(j)=FE(j)*GRID(ng)%pmask(i,j)
!^
                tl_FE(j)=tl_FE(j)*GRID(ng)%pmask(i,j)
# endif
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,u2dvar)
              DO i=IstrU-1,Iend
!^              FX(i)=GRID(ng)%pmon_r(i,j)*                             &
!^                    Khx(i)*(Awrk(i+1)-Awrk(i))
!^
                tl_FX(i)=GRID(ng)%pmon_r(i,j)*                          &
     &                   Khx(i)*(tl_Awrk(i+1)-tl_Awrk(i))
              END DO
            END IF
          END IF
!
        CASE (v2dvar)                           ! 2D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,v2dvar)
              DO j=JstrV-1,Jend
!^              FE(j)=GRID(ng)%pnom_r(i,j)*                             &
!^   &                Khy(j)*(Awrk(j+1)-Awrk(j))
!^
                tl_FE(j)=GRID(ng)%pnom_r(i,j)*                          &
     &                   Khy(j)*(tl_Awrk(j+1)-tl_Awrk(j))
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,v2dvar)
              DO i=Istr,Iend+1
!^              FX(i)=GRID(ng)%pmon_p(i,j)*                             &
!^   &                0.5_r8*(Khx(i-1)+Khx(i))*                         &
!^   &                (Awrk(i)-Awrk(i-1))
!^
                tl_FX(i)=GRID(ng)%pmon_p(i,j)*                          &
     &                   0.5_r8*(Khx(i-1)+Khx(i))*                      &
     &                   (tl_Awrk(i)-tl_Awrk(i-1))
# ifdef MASKING
!^              FX(i)=FX(i)*GRID(ng)%pmask(i,j)
!^
                tl_FX(i)=tl_FX(i)*GRID(ng)%pmask(i,j)
# endif
            END DO
          END IF
        END IF
!
      END SELECT
!
!  Compute lateral boundary edge K-Laplacian.
!
      SELECT CASE (ctype)
!
        CASE (r2dvar)                           ! 2D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO j=Jstr,Jend
!^              Awrk(j)=Awrk(j)-                                        &
!^   &                  Hfac(j)*(FE(j+1)-FE(j))
!^
                tl_Awrk(j)=tl_Awrk(j)-                                  &
     &                     Hfac(j)*(tl_FE(j+1)-tl_FE(j))
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO i=Istr,Iend
!^              Awrk(i)=Awrk(i)-                                        &
!^   &                  Hfac(i)*(FX(i+1)-FX(i))
!^
                tl_Awrk(i)=tl_Awrk(i)-                                  &
     &                     Hfac(i)*(tl_FX(i+1)-tl_FX(i))
              END DO
            END IF
          END IF
!
        CASE (u2dvar)                           ! 2D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO j=Jstr,Jend
!^              Awrk(j)=Awrk(j)-                                        &
!^   &                  Hfac(j)*(FE(j+1)-FE(j))
!^
                tl_Awrk(j)=tl_Awrk(j)-                                  &
     &                     Hfac(j)*(tl_FE(j+1)-tl_FE(j))
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO i=IstrU,Iend
!^              Awrk(i)=Awrk(i)-                                        &
!^   &                  Hfac(i)*(FX(i)-FX(i-1))
!^
                tl_Awrk(i)=tl_Awrk(i)-                                  &
     &                     Hfac(i)*(tl_FX(i)-tl_FX(i-1))
              END DO
            END IF
          END IF
!
        CASE (v2dvar)                           ! 2D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO j=JstrV,Jend
!^              Awrk(j)=Awrk(j)-                                        &
!^   &                  Hfac(j)*(FE(j)-FE(j-1))
!^
                tl_Awrk(j)=tl_Awrk(j)-                                  &
     &                     Hfac(j)*(tl_FE(j)-tl_FE(j-1))
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO i=Istr,Iend
!^              Awrk(i)=Awrk(i)-                                        &
!^   &                  Hfac(i)*(FX(i+1)-FX(i))
!^
                tl_Awrk(i)=tl_Awrk(i)-                                  &
     &                     Hfac(i)*(tl_FX(i+1)-tl_FX(i))
              END DO
            END IF
          END IF
!
      END SELECT
!
!  Load K-Laplacian solution.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
!^            A(j)=Awrk(j)
!^
              tl_A(j)=tl_Awrk(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
!^            A(i)=Awrk(i)
!^
              tl_A(i)=tl_Awrk(i)
            END DO
        END SELECT
      END IF

# ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij,                               &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        A)
!^
      CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
     &                        LBij, UBij,                               &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        tl_A)
# endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_b1d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of lateral
!  boundary conditions for a 2D control variable.
!
      SUBROUTINE multiscale_Klap_b1d_ad (self, ng, tile, model,         &
     &                                   ifield, ibry, ctype, ms,       &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Mlap, i, j
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
!
      real (r8)                         :: adfac, cffx, cffy
!
      real (r8), dimension(LBij:UBij)   :: ad_Awrk
!
      real (r8), dimension(LBij:UBij)   :: Hfac
      real (r8), dimension(IminS:ImaxS) :: Khx
      real (r8), dimension(JminS:JmaxS) :: Khy
      real (r8), dimension(JminS:JmaxS) :: ad_FE
      real (r8), dimension(IminS:ImaxS) :: ad_FX

!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
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
!  Initialize private adjoint variables.
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8     
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta', 'ubar', 'ubar_eastward', 'vbar', 'vbar_northward')
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jstr-1,Jend+1
                 Hfac(j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Istr-1,Iend+1
                 Hfac(i)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jstr-1,Jend+1
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=IstrU-1,Iend+1
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=JstrV-1,Jend+1
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Istr-1,Iend+1
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=1, kappa=D*D/(2*M-d-2), where D is the Daley length scale
!  and d is the space dimension.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO j=JminS,JmaxS
            cffy=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khy(j)=cffy/REAL(2*Mlap-3,r8)
          END DO
        CASE (isouth, inorth)
          DO i=IminS,ImaxS
            cffx=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khx(i)=cffx/REAL(2*Mlap-3,r8)
          END DO
      END SELECT
!
!  Adjoint of load K-Laplacian solution.
!
# ifdef DISTRIBUTE
!^    CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij,                               &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_A)
!^
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij,                            &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
!
# endif
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
!^            tl_A(j)=tl_Awrk(j)
!^
              ad_Awrk(j)=ad_Awrk(j)+ad_A(j)
              ad_A(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
!^            tl_A(i)=tl_Awrk(i)
!^
              ad_Awrk(i)=ad_Awrk(i)+ad_A(i)
              ad_A(i)=0.0_r8
            END DO
        END SELECT
      END IF
!
!  Adjoint of compute lateral boundary edge K-Laplacian.
!
      SELECT CASE (ctype)
!
        CASE (r2dvar)                           ! 2D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO j=Jstr,Jend
!^              tl_Awrk(j)=tl_Awrk(j)-                                  &
!^   &                     Hfac(j)*(tl_FE(j+1)-tl_FE(j))
!^
                adfac=-Hfac(j)*ad_Awrk(j)
                ad_FE(j  )=ad_FE(j  )-adfac
                ad_FE(j+1)=ad_FE(j+1)+adfac
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO i=Istr,Iend
!^              tl_Awrk(i)=tl_Awrk(i)-                                  &
!^   &                     Hfac(i)*(tl_FX(i+1)-tl_FX(i))
!^
                adfac=-Hfac(i)*ad_Awrk(i)
                ad_FX(i  )=ad_FX(i  )-adfac
                ad_FX(i+1)=ad_FX(i+1)+adfac
              END DO
            END IF
          END IF
!
        CASE (u2dvar)                           ! 2D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO j=Jstr,Jend
!^              tl_Awrk(j)=tl_Awrk(j)-                                  &
!^   &                     Hfac(j)*(tl_FE(j+1)-tl_FE(j))
!^
                adfac=-Hfac(j)*ad_Awrk(j)
                ad_FE(j  )=ad_FE(j  )-adfac
                ad_FE(j+1)=ad_FE(j+1)+adfac
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO i=IstrU,Iend
!^              tl_Awrk(i)=tl_Awrk(i)-                                  &
!^   &                     Hfac(i)*(tl_FX(i)-tl_FX(i-1))
!^
                adfac=-Hfac(i)*ad_Awrk(i)
                ad_FX(i-1)=ad_FX(i-1)-adfac
                ad_FX(i  )=ad_FX(i  )+adfac
              END DO
            END IF
          END IF
!
        CASE (v2dvar)                           ! 2D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO j=JstrV,Jend
!^              tl_Awrk(j)=tl_Awrk(j)-                                  &
!^   &                     Hfac(j)*(tl_FE(j)-tl_FE(j-1))
!^
                adfac=-Hfac(j)*ad_Awrk(j)
                ad_FE(j-1)=ad_FE(j-1)-adfac
                ad_FE(j  )=ad_FE(j  )+adfac
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO i=Istr,Iend
!^              tl_Awrk(i)=tl_Awrk(i)-                                  &
!^   &                     Hfac(i)*(tl_FX(i+1)-tl_FX(i))
!^
                adfac=-Hfac(i)*ad_Awrk(i)
                ad_FX(i  )=ad_FX(i  )-adfac
                ad_FX(i+1)=ad_FX(i+1)+adfac
              END DO
            END IF
          END IF
!
      END SELECT
!
!  Adjoint of compute XI- or ETA-components of diffusive flux.
!
      SELECT CASE (ctype)
!
        CASE (r2dvar)                           ! 2D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,r2dvar)
              DO j=Jstr,Jend+1
# ifdef MASKING
!^              tl_FE(j)=tl_FE(j)*GRID(ng)%vmask(i,j)
!^
                ad_FE(j)=ad_FE(j)*GRID(ng)%vmask(i,j)
# endif
!^              tl_FE(j)=GRID(ng)%pnom_v(i,j)*                          &
!^   &                   0.5_r8*(Khy(j-1)+Khy(j))*                      &
!^   &                   (tl_Awrk(j)-tl_Awrk(j-1))
!^
                adfac=GRID(ng)%pnom_v(i,j)*                             &
     &                0.5_r8*(Khy(j-1)+Khy(j))*ad_FE(j)
                ad_Awrk(j-1)=ad_Awrk(j-1)-adfac
                ad_Awrk(j  )=ad_Awrk(j  )+adfac
                ad_FE(j)=0.0_r8
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,r2dvar)
              DO i=Istr,Iend+1
# ifdef MASKING
!^              tl_FX(i)=tl_FX(i)*GRID(ng)%umask(i,j)
!^
                ad_FX(i)=ad_FX(i)*GRID(ng)%umask(i,j)
# endif
!^              tl_FX(i)=GRID(ng)%pmon_u(i,j)*                          &
!^   &                   0.5_r8*(Khx(i-1)+Khx(i))*                      &
!^   &                   (tl_Awrk(i)-tl_Awrk(i-1))
!^
                adfac=GRID(ng)%pmon_u(i,j)*                             &
     &                0.5_r8*(Khx(i-1)+Khx(i))*ad_FX(i)
                ad_Awrk(i-1)=ad_Awrk(i-1)-adfac
                ad_Awrk(i  )=ad_Awrk(i  )+adfac
                ad_FX(i)=0.0_r8
              END DO
            END IF
          END IF
!
        CASE (u2dvar)                           ! 2D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,u2dvar)
              DO j=Jstr,Jend+1
# ifdef MASKING
!^              tl_FE(j)=tl_FE(j)*GRID(ng)%pmask(i,j)
!^
                ad_FE(j)=ad_FE(j)*GRID(ng)%pmask(i,j)
# endif
!^              tl_FE(j)=GRID(ng)%pnom_p(i,j)*                          &
!^   &                   0.5_r8*(Khy(j)+Khy(j-1))*                      &
!^   &                   (tl_Awrk(j)-tl_Awrk(j-1))
!^
                adfac=GRID(ng)%pnom_p(i,j)*                             &
     &                0.5_r8*(Khy(j)+Khy(j-1))*ad_FE(j)
                ad_Awrk(j-1)=ad_Awrk(j-1)-adfac
                ad_Awrk(j  )=ad_Awrk(j  )+adfac
                ad_FE(j)=0.0_r8
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,u2dvar)
              DO i=IstrU-1,Iend
!^              tl_FX(i)=GRID(ng)%pmon_r(i,j)*                          &
!^   &                   Khx(i)*(tl_Awrk(i+1)-tl_Awrk(i))
!^
                adfac=GRID(ng)%pmon_r(i,j)*Khx(i)*ad_FX(i)
                ad_Awrk(i  )=ad_Awrk(i  )-adfac
                ad_Awrk(i+1)=ad_Awrk(i+1)+adfac
                ad_FX(i)=0.0_r8
              END DO
            END IF
          END IF
!
        CASE (v2dvar)                           ! 2D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,v2dvar)
              DO j=JstrV-1,Jend
!^              tl_FE(j)=GRID(ng)%pnom_r(i,j)*                          &
!^   &                   Khy(j)*(tl_Awrk(j+1)-tl_Awrk(j))
!^
                adfac=GRID(ng)%pnom_r(i,j)*Khy(j)*ad_FE(j)
                ad_Awrk(j  )=ad_Awrk(j  )-adfac
                ad_Awrk(j+1)=ad_Awrk(j+1)+adfac
                ad_FE(j)=0.0_r8
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,v2dvar)
              DO i=Istr,Iend+1
# ifdef MASKING
!^              tl_FX(i)=tl_FX(i)*GRID(ng)%pmask(i,j)
!^
                ad_FX(i)=ad_FX(i)*GRID(ng)%pmask(i,j)
# endif
!^              tl_FX(i)=GRID(ng)%pmon_p(i,j)*                          &
!^   &                   0.5_r8*(Khx(i-1)+Khx(i))*                      &
!^   &                   (tl_Awrk(i)-tl_Awrk(i-1))
!^
                adfac=GRID(ng)%pmon_p(i,j)*                             &
     &                0.5_r8*(Khx(i-1)+Khx(i))*ad_FX(i)
                ad_Awrk(i-1)=ad_Awrk(i-1)-adfac
                ad_Awrk(i  )=ad_Awrk(i  )+adfac
                ad_FX(i)=0.0_r8
              END DO
            END IF
          END IF
!
      END SELECT
!
!  Set adjoint operator initial conditions.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin-1,Jmax+1
!^            Awrk(j)=A(j)
!^
              ad_Awrk(j)=ad_A(j)
              ad_Awrk(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin-1,Imax+1
!^            Awrk(i)=A(i)
!^
              ad_Awrk(i)=ad_A(i)
              ad_Awrk(i)=0.0_r8
            END DO
        END SELECT
      END IF

# ifdef DISTRIBUTE
!
!^    CALL mp_exchange2d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij,                               &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_A)
!^
      CALL ad_mp_exchange2d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij,                            &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
# endif
!
      SELECT CASE (ctype)
        CASE (r2dvar)
!^        CALL bc_r2d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij,                             &
!^   &                          tl_A)
!^
          CALL ad_bc_r2d_bry_tile (ng, tile, ibry,                      &
     &                             LBij, UBij,                          &
     &                             ad_A)
        CASE (u2dvar)
!^        CALL bc_u2d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij,                             &
!^   &                          tl_A)
!^
          CALL ad_bc_u2d_bry_tile (ng, tile, ibry,                      &
     &                             LBij, UBij,                          &
     &                             ad_A)
        CASE (v2dvar)
!^        CALL bc_v2d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij,                             &
!^   &                          tl_A)
!^
          CALL ad_bc_v2d_bry_tile (ng, tile, ibry,                      &
     &                             LBij, UBij,                          &
     &                             ad_A)
      END SELECT
!
      RETURN
      END SUBROUTINE multiscale_Klap_b1d_ad

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  It computes the tangent linear K-Laplacian, [1 + Del(K*Del)], of
!  lateral boundary conditions for a 3D control variable.
!
      SUBROUTINE multiscale_Klap_b2d_tl (self, ng, tile, model,         &
     &                                   ifield, ibry, ctype, ms,       &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Mlap, i, j, k
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
!
      real (r8)                         :: cffx, cffy
!
      real (r8), dimension(LBij:UBij,N(ng))   :: tl_Awrk
!
      real (r8), dimension(LBij:UBij)         :: Hfac
      real (r8), dimension(IminS:ImaxS)       :: Khx
      real (r8), dimension(JminS:JmaxS)       :: Khy
      real (r8), dimension(JminS:JmaxS,N(ng)) :: tl_FE
      real (r8), dimension(IminS:ImaxS,N(ng)) :: tl_FX
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
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
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8
      tl_Awrk=0.0_r8
      tl_FE=0.0_r8
      tl_FX=0.0_r8
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward', 'v', 'v_northward')
          Mlap=self%Mlap(ifield,ms)
        CASE ('temp', 'salt')
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               DO j=Jstr-1,Jend+1
                 i=BOUNDS(ng)%edge(ibry,r2dvar)
                 Hfac(j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               DO i=Istr-1,Iend+1
                 j=BOUNDS(ng)%edge(ibry,r2dvar)
                 Hfac(i)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             END IF
           CASE (u3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               DO j=Jstr-1,Jend+1
                 i=BOUNDS(ng)%edge(ibry,u2dvar)
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               DO i=IstrU-1,Iend+1
                 j=BOUNDS(ng)%edge(ibry,u2dvar)
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             END IF
           CASE (v3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               DO j=JstrV-1,Jend+1
                 i=BOUNDS(ng)%edge(ibry,v2dvar)
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               DO i=Istr-1,Iend+1
                 j=BOUNDS(ng)%edge(ibry,v2dvar)
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=1, kappa=D*D/(2*M-d-2), where D is the Daley length scale
!  and d is the space dimension.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO j=Jmin-1,Jmax+1
            cffy=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khy(j)=cffy/REAL(2*Mlap-3,r8)
          END DO
        CASE (isouth, inorth)
          DO i=Imin-1,Imax+1
            cffx=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khx(i)=cffx/REAL(2*Mlap-3,r8)
          END DO
      END SELECT
!
!  Set operator initial conditions.
!
      SELECT CASE (ctype)
        CASE (r3dvar)
!^        CALL bc_r3d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij, 1, N(ng),                   &
!^   &                          A)
!^
          CALL bc_r3d_bry_tile (ng, tile, ibry,                         &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          tl_A)
        CASE (u3dvar)
!^        CALL bc_u3d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij, 1, N(ng),                   &
!^   &                          A)
!^
          CALL bc_u3d_bry_tile (ng, tile, ibry,                         &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          tl_A)
        CASE (v3dvar)
!^        CALL bc_u3d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij, 1, N(ng),                   &
!^   &                          A)
!^
          CALL bc_v3d_bry_tile (ng, tile, ibry,                         &
     &                          LBij, UBij, 1, N(ng),                   &
     &                          tl_A)
      END SELECT

#  ifdef DISTRIBUTE
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
#  endif
!
      IF (Lboundary(ibry)) THEN
        IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN 
          DO k=1,N(ng)
           DO j=Jmin-1,Jmax+1
!^            Awrk(j,k)=A(j,k)
!^
              tl_Awrk(j,k)=tl_A(j,k)
            END DO
          END DO
        ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Imin-1,Imax+1
!^            Awrk(i,k)=A(i,k)
!^
              tl_Awrk(i,k)=tl_A(i,k)
            END DO
          END DO
        END IF
      END IF
!
!  Compute XI- or ETA-components of diffusive flux.
!
      SELECT CASE (ctype)
!
        CASE (r3dvar)                           ! 3D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,r2dvar)
              DO k=1,N(ng)
                DO j=Jstr,Jend+1
!^                FE(j,k)=GRID(ng)%pnom_v(i,j)*                         &
!^   &                    0.5_r8*(Khy(j-1)+Khy(j))*                     &
!^   &                    (Awrk(j,k)-Awrk(j-1,k))
!^
                  tl_FE(j,k)=GRID(ng)%pnom_v(i,j)*                      &
     &                       0.5_r8*(Khy(j-1)+Khy(j))*                  &
     &                       (tl_Awrk(j,k)-tl_Awrk(j-1,k))
#  ifdef MASKING
!^                FE(j,k)=FE(j,k)*GRID(ng)%vmask(i,j)
!^
                  tl_FE(j,k)=tl_FE(j,k)*GRID(ng)%vmask(i,j)
#  endif
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,r2dvar)
              DO k=1,N(ng)
                DO i=Istr,Iend+1
!^                FX(i,k)=GRID(ng)%pmon_u(i,j)*                         &
!^                        0.5_r8*(Khx(i-1)+Khx(i))*                     &
!^   &                    (Awrk(i,k)-Awrk(i-1,k))
!^
                  tl_FX(i,k)=GRID(ng)%pmon_u(i,j)*                      &
     &                       0.5_r8*(Khx(i-1)+Khx(i))*                  &
     &                       (tl_Awrk(i,k)-tl_Awrk(i-1,k))
#  ifdef MASKING
!^                FX(i,k)=FX(i,k)*GRID(ng)%umask(i,j)
!^
                  tl_FX(i,k)=tl_FX(i,k)*GRID(ng)%umask(i,j)
#  endif
                END DO
              END DO
            END IF
          END IF
!
        CASE (u3dvar)                           ! 3D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,u2dvar)
              DO k=1,N(ng)
                DO j=Jstr,Jend+1
!^                FE(j,k)=GRID(ng)%pnom_p(i,j)*                         &
!^   &                    0.5_r8*(Khy(j)+Khy(j-1))*                     &
!^   &                    (Awrk(j,k)-Awrk(j-1,k))
!^
                  tl_FE(j,k)=GRID(ng)%pnom_p(i,j)*                      &
     &                       0.5_r8*(Khy(j)+Khy(j-1))*                  &
     &                       (tl_Awrk(j,k)-tl_Awrk(j-1,k))
#  ifdef MASKING
!^                FE(j,k)=FE(j,k)*GRID(ng)%pmask(i,j)
!^
                  tl_FE(j,k)=tl_FE(j,k)*GRID(ng)%pmask(i,j)
#  endif
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,u2dvar)
              DO k=1,N(ng)
                DO i=IstrU-1,Iend
!^                FX(i,k)=GRID(ng)%pmon_r(i,j)*                         &
!^   &                    Khx(i)*(Awrk(i+1,k)-Awrk(i,k))
!^
                  tl_FX(i,k)=GRID(ng)%pmon_r(i,j)*                      &
     &                       Khx(i)*(tl_Awrk(i+1,k)-tl_Awrk(i,k))
                END DO
              END DO
            END IF
          END IF
!
        CASE (v3dvar)                           ! 3D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,v2dvar)
              DO k=1,N(ng)
                DO j=JstrV-1,Jend
!^                FE(j,k)=GRID(ng)%pnom_r(i,j)*                         &
!^   &                    Khy(j)*(Awrk(j+1,k)-Awrk(j,k))
!^
                  tl_FE(j,k)=GRID(ng)%pnom_r(i,j)*                      &
     &                       Khy(j)*(tl_Awrk(j+1,k)-tl_Awrk(j,k))
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,v2dvar)
              DO k=1,N(ng)
                DO i=Istr,Iend+1
!^                FX(i,k)=GRID(ng)%pmon_p(i,j)*                         &
!^   &                    0.5_r8*(Khx(i-1)+Khx(i))*                     &
!^   &                    (Awrk(i,k)-Awrk(i-1,k))
!^
                  tl_FX(i,k)=GRID(ng)%pmon_p(i,j)*                      &
     &                       0.5_r8*(Khx(i-1)+Khx(i))*                  &
     &                       (tl_Awrk(i,k)-tl_Awrk(i-1,k))
#  ifdef MASKING
!^                FX(i,k)=FX(i,k)*GRID(ng)%pmask(i,j)
!^
                  tl_FX(i,k)=tl_FX(i,k)*GRID(ng)%pmask(i,j)
#  endif
              END DO
            END DO
          END IF
        END IF
!
      END SELECT
!
!  Compute lateral boundary edge K-Laplacian.
!
      SELECT CASE (ctype)
!
        CASE (r3dvar)                           ! 3D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
!^                Awrk(j,k)=Awrk(j,k)-                                  &
!^   &                      Hfac(j)*(FE(j+1,k)-FE(j,k))
!^
                  tl_Awrk(j,k)=tl_Awrk(j,k)-                            &
     &                         Hfac(j)*(tl_FE(j+1,k)-tl_FE(j,k))
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
!^                Awrk(i,k)=Awrk(i,k)-                                  &
!^   &                      Hfac(i)*(FX(i+1,k)-FX(i,k))
!^
                  tl_Awrk(i,k)=tl_Awrk(i,k)-                            &
     &                         Hfac(i)*(tl_FX(i+1,k)-tl_FX(i,k))
                END DO
              END DO
            END IF
          END IF
!
        CASE (u3dvar)                           ! 3D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
!^                Awrk(j,k)=Awrk(j,k)-                                  &
!^   &                      Hfac(j)*(FE(j+1,k)-FE(j,k))
!^
                  tl_Awrk(j,k)=tl_Awrk(j,k)-                            &
     &                         Hfac(j)*(tl_FE(j+1,k)-tl_FE(j,k))
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO k=1,N(ng)
                DO i=IstrU,Iend
!^                Awrk(i,k)=Awrk(i,k)-                                  &
!^   &                      Hfac(i)*(FX(i,k)-FX(i-1,k))
!^
                  tl_Awrk(i,k)=tl_Awrk(i,k)-                            &
     &                         Hfac(i)*(tl_FX(i,k)-tl_FX(i-1,k))
                END DO
              END DO
            END IF
          END IF
!
        CASE (v3dvar)                           ! 3D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO k=1,N(ng)
                DO j=JstrV,Jend
!^                Awrk(j,k)=Awrk(j,k)-                                  &
!^   &                      Hfac(j)*(FE(j,k)-FE(j-1,k))
!^
                  tl_Awrk(j,k)=tl_Awrk(j,k)-                            &
     &                         Hfac(j)*(tl_FE(j,k)-tl_FE(j-1,k))
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
!^                Awrk(i,k)=Awrk(i,k)-                                  &
!^   &                      Hfac(i)*(FX(i+1,k)-FX(i,k))
!^
                  tl_Awrk(i,k)=tl_Awrk(i,k)-                            &
     &                         Hfac(i)*(tl_FX(i+1,k)-tl_FX(i,k))
                END DO
              END DO
            END IF
          END IF
!
      END SELECT
!
!  Load K-Laplacian solution.
!
      IF (Lboundary(ibry)) THEN
        IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jmin,Jmax
!^            A(j,k)=Awrk(j,k)
!^
              tl_A(j,k)=tl_Awrk(j,k)
            END DO
          END DO
        ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Imin,Imax
!^            A(i,k)=Awrk(i,k)
!^
              tl_A(i,k)=tl_Awrk(i,k)
            END DO
          END DO
        END IF
      END IF

#  ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        Awrk)
!^
      CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        tl_Awrk)
#  endif
!
      RETURN
      END SUBROUTINE multiscale_Klap_b2d_tl
!
!-----------------------------------------------------------------------
!  It computes the adjoint K-Laplacian, [1 + Del(K*Del)], of lateral
!  boundary conditions for a 2D control variable.
!
      SUBROUTINE multiscale_Klap_b2d_ad (self, ng, tile, model,         &
     &                                   ifield, ibry, ctype, ms,       &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Mlap, i, j, k
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
!
      real (r8)                         :: adfac, cffx, cffy
!
      real (r8), dimension(LBij:UBij,N(ng))   :: ad_Awrk
!
      real (r8), dimension(LBij:UBij)         :: Hfac
      real (r8), dimension(IminS:ImaxS)       :: Khx
      real (r8), dimension(JminS:JmaxS)       :: Khy
      real (r8), dimension(JminS:JmaxS,N(ng)) :: ad_FE
      real (r8), dimension(IminS:ImaxS,N(ng)) :: ad_FX
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational indices
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
!  Initialize private adjoint variables.
!
      ad_Awrk=0.0_r8
      ad_FE=0.0_r8
      ad_FX=0.0_r8
      Hfac=0.0_r8
      Khx=0.0_r8
      Khy=0.0_r8     
!
!  Assign contol variable isotropic or anisotropic correlation length
!  scales.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward', 'v', 'v_northward')
          Mlap=self%Mlap(ifield,ms)
        CASE ('temp', 'salt')
          Mlap=self%Mlap(ifield,ms)
      END SELECT
!
!  Compute metrics factor.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               DO j=Jstr-1,Jend+1
                 i=BOUNDS(ng)%edge(ibry,r2dvar)
                 Hfac(j)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               DO i=Istr-1,Iend+1
                 j=BOUNDS(ng)%edge(ibry,r2dvar)
                 Hfac(i)=GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
               END DO
             END IF
           CASE (u3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               DO j=Jstr-1,Jend+1
                 i=BOUNDS(ng)%edge(ibry,u2dvar)
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               DO i=IstrU-1,Iend+1
                 j=BOUNDS(ng)%edge(ibry,u2dvar)
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i-1,j)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i-1,j)+GRID(ng)%pn(i,j))
               END DO
             END IF
           CASE (v3dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               DO j=JstrV-1,Jend+1
                 i=BOUNDS(ng)%edge(ibry,v2dvar)
                 Hfac(j)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               DO i=Istr-1,Iend+1
                 j=BOUNDS(ng)%edge(ibry,v2dvar)
                 Hfac(i)=0.25_r8*(GRID(ng)%pm(i,j-1)+GRID(ng)%pm(i,j))* &
     &                           (GRID(ng)%pn(i,j-1)+GRID(ng)%pn(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Set horizontal diffusion coefficients (Khx, Khy) with units of
!  correlation length squared (Equation 44, Weaver and Mirouze, 2013).
!  For d=1, kappa=D*D/(2*M-d-2), where D is the Daley length scale
!  and d is the space dimension.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO j=Jmin-1,Jmax+1
            cffy=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khy(j)=cffy/REAL(2*Mlap-3,r8)
          END DO
        CASE (isouth, inorth)
          DO i=Imin-1,Imax+1
            cffx=HdecayB(ifield,ms,ibry,ng)*HdecayB(ifield,ms,ibry,ng)
            Khx(i)=cffx/REAL(2*Mlap-3,r8)
          END DO
      END SELECT
!
!  Adjoint of load K-Laplacian solution.
!
#  ifdef DISTRIBUTE
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_Awrk)
!^
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
!
#  endif
      IF (Lboundary(ibry)) THEN
        IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jmin,Jmax
!^            tl_A(j,k)=tl_Awrk(j,k)
!^
              ad_Awrk(j,k)=ad_Awrk(j,k)+ad_A(j,k)
              ad_A(j,k)=0.0_r8
            END DO
          END DO
        ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Imin,Imax
!^            tl_A(i,k)=tl_Awrk(i,k)
!^
              ad_Awrk(i,k)=ad_Awrk(i,k)+ad_A(i,k)
              ad_A(i,k)=0.0_r8
            END DO
          END DO
        END IF
      END IF
!
!  Adjoint of compute lateral boundary edge K-Laplacian.
!
      SELECT CASE (ctype)
!
        CASE (r3dvar)                           ! 3D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
!^                tl_Awrk(j,k)=tl_Awrk(j,k)-                            &
!^   &                         Hfac(j)*(tl_FE(j+1,k)-tl_FE(j,k))
!^
                  adfac=-Hfac(j)*ad_Awrk(j,k)
                  ad_FE(j  ,k)=ad_FE(j  ,k)-adfac
                  ad_FE(j+1,k)=ad_FE(j+1,k)+adfac
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
!^                tl_Awrk(i,k)=tl_Awrk(i,k)-                            &
!^   &                         Hfac(i)*(tl_FX(i+1,k)-tl_FX(i,k))
!^
                  adfac=-Hfac(i)*ad_Awrk(i,k)
                  ad_FX(i  ,k)=ad_FX(i  ,k)-adfac
                  ad_FX(i+1,k)=ad_FX(i+1,k)+adfac
                END DO
              END DO
            END IF
          END IF
!
        CASE (u3dvar)                           ! 3D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO k=1,N(ng)
                DO j=Jstr,Jend
!^                tl_Awrk(j,k)=tl_Awrk(j,k)-                            &
!^   &                         Hfac(j)*(tl_FE(j+1,k)-tl_FE(j,k))
!^
                  adfac=-Hfac(j)*ad_Awrk(j,k)
                  ad_FE(j  ,k)=ad_FE(j  ,k)-adfac
                  ad_FE(j+1,k)=ad_FE(j+1,k)+adfac
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO k=1,N(ng)
                DO i=IstrU,Iend
!^                tl_Awrk(i,k)=tl_Awrk(i,k)-                            &
!^   &                         Hfac(i)*(tl_FX(i,k)-tl_FX(i-1,k))
!^
                  adfac=-Hfac(i)*ad_Awrk(i,k)
                  ad_FX(i-1,k)=ad_FX(i-1,k)-adfac
                  ad_FX(i  ,k)=ad_FX(i  ,k)+adfac
                END DO
              END DO
            END IF
          END IF
!
        CASE (v3dvar)                           ! 3D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              DO k=1,N(ng)
                DO j=JstrV,Jend
!^                tl_Awrk(j,k)=tl_Awrk(j,k)-                            &
!^   &                         Hfac(j)*(tl_FE(j,k)-tl_FE(j-1,k))
!^
                  adfac=-Hfac(j)*ad_Awrk(j,k)
                  ad_FE(j-1,k)=ad_FE(j-1,k)-adfac
                  ad_FE(j  ,k)=ad_FE(j  ,k)+adfac
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              DO k=1,N(ng)
                DO i=Istr,Iend
!^                tl_Awrk(i,k)=tl_Awrk(i,k)-                            &
!^   &                         Hfac(i)*(tl_FX(i+1,k)-tl_FX(i,k))
!^
                  adfac=-Hfac(i)*ad_Awrk(i,k)
                  ad_FX(i  ,k)=ad_FX(i  ,k)-adfac
                  ad_FX(i+1,k)=ad_FX(i+1,k)+adfac
                END DO
              END DO
            END IF
          END IF
!
      END SELECT
!
!  Adjoint of compute XI- or ETA-components of diffusive flux.
!
      SELECT CASE (ctype)
!
        CASE (r3dvar)                           ! 3D RHO-grid variables
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,r2dvar)
              DO k=1,N(ng)
                DO j=Jstr,Jend+1
#  ifdef MASKING
!^                tl_FE(j,k)=tl_FE(j,k)*GRID(ng)%vmask(i,j)
!^
                  ad_FE(j,k)=ad_FE(j,k)*GRID(ng)%vmask(i,j)
#  endif
!^                tl_FE(j,k)=GRID(ng)%pnom_v(i,j)*                      &
!^   &                       0.5_r8*(Khy(j-1)+Khy(j))*                  &
!^   &                       (tl_Awrk(j,k)-tl_Awrk(j-1,k))
!^
                  adfac=GRID(ng)%pnom_v(i,j)*                           &
     &                  0.5_r8*(Khy(j-1)+Khy(j))*ad_FE(j,k)
                  ad_Awrk(j-1,k)=ad_Awrk(j-1,k)-adfac
                  ad_Awrk(j  ,k)=ad_Awrk(j  ,k)+adfac
                  ad_FE(j,k)=0.0_r8
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,r2dvar)
              DO k=1,N(ng)
                DO i=Istr,Iend+1
#  ifdef MASKING
!^                tl_FX(i,k)=tl_FX(i,k)*GRID(ng)%umask(i,j)
!^
                  ad_FX(i,k)=ad_FX(i,k)*GRID(ng)%umask(i,j)
#  endif
!^                tl_FX(i,k)=GRID(ng)%pmon_u(i,j)*                      &
!^   &                       0.5_r8*(Khx(i-1)+Khx(i))*                  &
!^   &                       (tl_Awrk(i,k)-tl_Awrk(i-1,k))
!^
                  adfac=GRID(ng)%pmon_u(i,j)*                           &
     &                  0.5_r8*(Khx(i-1)+Khx(i))*ad_FX(i,k)
                  ad_Awrk(i-1,k)=ad_Awrk(i-1,k)-adfac
                  ad_Awrk(i  ,k)=ad_Awrk(i  ,k)+adfac
                  ad_FX(i,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
!
        CASE (u3dvar)                           ! 3D U-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,u2dvar)
              DO k=1,N(ng)
                DO j=Jstr,Jend+1
#  ifdef MASKING
!^                tl_FE(j,k)=tl_FE(j,k)*GRID(ng)%pmask(i,j)
!^
                  ad_FE(j,k)=ad_FE(j,k)*GRID(ng)%pmask(i,j)
#  endif
!^                tl_FE(j,k)=GRID(ng)%pnom_p(i,j)*                      &
!^   &                       0.5_r8*(Khy(j)+Khy(j-1))*                  &
!^   &                       (tl_Awrk(j,k)-tl_Awrk(j-1,k))
!^
                  adfac=GRID(ng)%pnom_p(i,j)*                           &
     &                  0.5_r8*(Khy(j)+Khy(j-1))*ad_FE(j,k)
                  ad_Awrk(j-1,k)=ad_Awrk(j-1,k)-adfac
                  ad_Awrk(j  ,k)=ad_Awrk(j  ,k)+adfac
                  ad_FE(j,k)=0.0_r8
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,u2dvar)
              DO k=1,N(ng)
                DO i=IstrU-1,Iend
!^                tl_FX(i,k)=GRID(ng)%pmon_r(i,j)*                      &
!^   &                       Khx(i)*(tl_Awrk(i+1,k)-tl_Awrk(i,k))
!^
                  adfac=GRID(ng)%pmon_r(i,j)*Khx(i)*ad_FX(i,k)
                  ad_Awrk(i  ,k)=ad_Awrk(i  ,k)-adfac
                  ad_Awrk(i+1,k)=ad_Awrk(i+1,k)+adfac
                  ad_FX(i,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
!
        CASE (v3dvar)                           ! 3D V-grid variables
!
          IF (Lboundary(ibry)) THEN
            IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
              i=BOUNDS(ng)%edge(ibry,v2dvar)
              DO k=1,N(ng)
                DO j=JstrV-1,Jend
!^                tl_FE(j,k)=GRID(ng)%pnom_r(i,j)*                      &
!^   &                       Khy(j)*(tl_Awrk(j+1,k)-tl_Awrk(j,k))
!^
                  adfac=GRID(ng)%pnom_r(i,j)*Khy(j)*ad_FE(j,k)
                  ad_Awrk(j  ,k)=ad_Awrk(j  ,k)-adfac
                  ad_Awrk(j+1,k)=ad_Awrk(j+1,k)+adfac
                  ad_FE(j,k)=0.0_r8
                END DO
              END DO
            ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
              j=BOUNDS(ng)%edge(ibry,v2dvar)
              DO k=1,N(ng)
                DO i=Istr,Iend+1
#  ifdef MASKING
!^                tl_FX(i,k)=tl_FX(i,k)*GRID(ng)%pmask(i,j)
!^
                  ad_FX(i,k)=ad_FX(i,k)*GRID(ng)%pmask(i,j)
#  endif
!^                tl_FX(i,k)=GRID(ng)%pmon_p(i,j)*                      &
!^   &                       0.5_r8*(Khx(i-1)+Khx(i))*                  &
!^   &                       (tl_Awrk(i,k)-tl_Awrk(i-1,k))
!^
                  adfac=GRID(ng)%pmon_p(i,j)*                           &
     &                  0.5_r8*(Khx(i-1)+Khx(i))*ad_FX(i,k)
                  ad_Awrk(i-1,k)=ad_Awrk(i-1,k)-adfac
                  ad_Awrk(i  ,k)=ad_Awrk(i  ,k)+adfac
                  ad_FX(i,k)=0.0_r8
                END DO
              END DO
            END IF
          END IF
!
      END SELECT
!
!  Adjoint of set operator initial conditions.
!
      IF (Lboundary(ibry)) THEN
        IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
          DO k=1,N(ng)
            DO j=Jmin-1,Jmax+1
!^            Awrk(j,k)=A(j,k)
!^
              ad_Awrk(j,k)=ad_A(j,k)
            END DO
          END DO
        ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
          DO k=1,N(ng)
            DO i=Imin-1,Imax+1
!^            Awrk(i,k)=A(i,k)
!^
              ad_Awrk(i,k)=ad_A(i,k)
            END DO
          END DO
        END IF
      END IF

#  ifdef DISTRIBUTE
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
#  endif
!
      SELECT CASE (ctype)
        CASE (r3dvar)
!^        CALL bc_r3d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij, 1, N(ng),                   &
!^   &                          tl_A)
!^
          CALL ad_bc_r3d_bry_tile (ng, tile, ibry,                      &
     &                             LBij, UBij, 1, N(ng),                &
     &                             ad_A)
        CASE (u3dvar)
!^        CALL bc_u3d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij, 1, N(ng),                   &
!^   &                          tl_A)
!^
          CALL ad_bc_u3d_bry_tile (ng, tile, ibry,                      &
     &                             LBij, UBij, 1, N(ng),                &
     &                             ad_A)
        CASE (v3dvar)
!^        CALL bc_v3d_bry_tile (ng, tile, ibry,                         &
!^   &                          LBij, UBij, 1, N(ng),                   &
!^   &                          tl_A)
!^
          CALL ad_bc_v3d_bry_tile (ng, tile, ibry,                      &
     &                             LBij, UBij, 1, N(ng),                &
     &                             ad_A)
      END SELECT
!
      RETURN
      END SUBROUTINE multiscale_Klap_b2d_ad

# endif /* SOLVE3D */
#endif /* ADJUST_BOUNDARY */
