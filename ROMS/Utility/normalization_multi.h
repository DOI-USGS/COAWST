#include "cppdefs.h"
      MODULE normalization_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Normalization Coefficients Multiscale Approach.                     !
!                                                                      !
!  This module computes normalization factors to model the spreading   !
!  of the background-error symmetric covariance matrix (B). These      !
!  coefficients ensure that the diagonal elements of B are equal to    !
!  unity. In applications involving land or sea masking, the           !
!  normalization factors can significantly modify the covariance       !
!  structure near boundaries.                                          !
!                                                                      !
!  Normalization coefficients may be computed using either the exact   !
!  method or the approximated randomization method:                    !
!                                                                      !
!  The exact method is computationally intensive. In this approach,    !
!  normalization factors are determined by perturbing each model grid  !
!  cell with a delta function scaled by area (for 2D factors) or       !
!  volume (for 3D factors), then convolving with the square-root       !
!  adjoint and tangent diffusion operators.                            !
!                                                                      !
!  The approximated method is less computationally demanding. In this  !
!  approach, normalization factors are computed using the randomization!
!  technique described by Fisher and Courtier (1995). Factors are      !
!  initialized with random numbers drawn from a normal distribution    !
!  with zero mean and unit variance. These factors are then scaled by  !
!  the inverse-square-root of cell area (for 2D factors) or volume     !
!  (for 3D factors), and convolved with the square-root adjoint and    !
!  tangent diffusion operators over a specified number of iterations,  !
!  denoted as Nrandom.                                                 !
!                                                                      !
!  Multiscale Background-error Covariance (B) Modeling/Spreading:      !
!                                                                      !
!  Spreading is modeled using pseudo-diffusion operators in spatial    !
!  correlation space, where diffusion coefficients (K) are proportional!
!  to the square of the correlation length scale (Daley, 1992).        !
!  Correlation scales may be constant or spatially varying for each    !
!  variable in the control vector. The background-error covariance     !
!  matrix (B) can be represented using either a multiscale aproach.    !
!  approach.                                                           !
!                                                                      !
!  In the multiscale formulation, B is expressed as a linear           !
!  combination of distinct spatial scales, ranging from large to small.!
!  This approach enables the representation of both broad structures   !
!  and fine-scale features, while reducing scale aliasing in the data  !
!  assimilation cost function (Weaver et al., 2016). In contrast, the  !
!  monoscale formulation employs a single scale.                       !
!                                                                      !
!  The multiscale approach uses an implicit horizontal pseudo-         !
!  diffusion operator implemented via Conjugate Gradient (CG) and      !
!  Chebyshev Iterations (CI). Conversely, the monoscale default method !
!  employs an explicit horizontal algorithm. Error correlations are    !
!  considered separable in the horizontal and vertical directions.     !
!  The vertical diffusion operator is also implicit.                   !
!                                                                      !
!  The multiscale concept, expressed as B = SUM(Wi*Bi), represents a   !
!  weighted sum of B values corresponding to various scales. Typically,!
!  two to four different scales are combined in practical applications.!
!  The weight coefficients Wi are constrained to sum to unity. The     !
!  resulting correlation functions belong to the Matern-class family,  !
!  which accommodates complex functional shapes.                       !
!                                                                      !
!  Horizontally spatially varying correlation length scales for the    !
!  K-diffusion coefficient are permitted. A horizontal map of isotropic!
!  or anisotropic correlations can be specified and is read from an    !
!  input NetCDF file. For computational efficiency, these correlations !
!  require an implicit diffusion CG/CI solver.                         !
!                                                                      !
!  References:                                                         !
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
      USE mod_kinds
      USE mod_param
      USE mod_parallel
#ifdef ADJUST_BOUNDARY
      USE mod_boundary
#endif
      USE mod_grid
#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
#endif
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_mixing
      USE mod_ocean
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
#if defined SEDIMENT && defined SED_MORPH && defined SOLVE3D
      USE mod_sedbed
#endif
!
      USE ad_conv_2d_mod
#ifdef SOLVE3D
      USE ad_conv_3d_mod
#endif
#ifdef ADJUST_BOUNDARY
      USE ad_conv_bry2d_mod
# ifdef SOLVE3D
      USE ad_conv_bry3d_mod
# endif
#endif
      USE bc_2d_mod
      USE ad_bc_2d_mod
#ifdef SOLVE3D
      USE bc_3d_mod
      USE ad_bc_3d_mod
#endif
#ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod
# ifdef SOLVE3D
      USE bc_bry3d_mod
# endif
#endif
#ifdef DISTRIBUTE
      USE mp_exchange_mod
#endif
      USE set_depth_mod
      USE tl_conv_2d_mod
#ifdef SOLVE3D
      USE tl_conv_3d_mod
#endif
#ifdef ADJUST_BOUNDARY
      USE tl_conv_bry2d_mod
# ifdef SOLVE3D
      USE tl_conv_bry3d_mod
# endif
#endif
      USE white_noise_mod
!
#ifdef DISTRIBUTE
# ifdef ADJUST_BOUNDARY
      USE distribute_mod,       ONLY : mp_collect
# endif
      USE distribute_mod,       ONLY : mp_reduce
#endif
      USE nf_fwrite2d_mod,      ONLY : nf_fwrite2d
      USE nf_fwrite3d_mod,      ONLY : nf_fwrite3d
      USE strings_mod,          ONLY : FoundError
!
      USE multiscale_eigen_mod, ONLY : multiscale_eigen  ! Eigenvalues
      USE roms_multiscale_mod,  ONLY : multiscale        ! CLASS object
!
      implicit none
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  Module public and private routines:
!
      PUBLIC  :: normalization
!
      PRIVATE :: normalization_tile
      PRIVATE :: randomization_tile
!
      PRIVATE :: dot_prod2d
#ifdef SOLVE3D
      PRIVATE :: dot_prod3d
#endif
!
      PRIVATE :: wrt_norm2d_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_norm2d_pio
#endif
      PRIVATE :: wrt_norm3d_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_norm3d_pio
#endif
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE normalization (ng, tile, ifac)
!***********************************************************************
!
      USE mod_stepping,        ONLY : nnew, nstp
      USE roms_multiscale_mod, ONLY : MSB
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng     ! nested grid
      integer, intent(in) :: tile   ! domain partition
      integer, intent(in) :: ifac   ! iteraction factor
                                    ! (squared-root filter, ifac=2)
!
!  Local variable declarations.
!
      integer :: my_ifac = 1        ! K-Laplacian iteraction factor
!
#include "tile.h"
!
!  Compute background error covariance normalization factors using
!  the very expensive exact method.
!
      IF (Nmethod(ng).eq.0) THEN
        CALL normalization_tile (MSB(ng), ng, tile,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp(ng), nnew(ng), ifac,              &
     &                           MIXING(ng) % Kh,                       &
#ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                           BOUNDARY(ng) % b_t_obc,                &
     &                           BOUNDARY(ng) % b_u_obc,                &
     &                           BOUNDARY(ng) % b_v_obc,                &
# endif
     &                           BOUNDARY(ng) % b_ubar_obc,             &
     &                           BOUNDARY(ng) % b_vbar_obc,             &
     &                           BOUNDARY(ng) % b_zeta_obc,             &
#endif
#ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                           FORCES(ng) % b_stflx,                  &
#endif
#ifdef SOLVE3D
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
#endif
     &                           OCEAN(ng) % b_zeta,                    &
     &                           OCEAN(ng) % b_ubar,                    &
     &                           OCEAN(ng) % b_vbar)
!
!  Compute background error covariance normalization factors using
!  the approximated randomization method.
!
      ELSE IF (Nmethod(ng).eq.1) THEN
        CALL randomization_tile (MSB(ng), ng, tile,                     &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           LBij, UBij,                            &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           nstp(ng), nnew(ng), ifac,              &
     &                           MIXING(ng) % Kh,                       &
#ifdef SOLVE3D
     &                           MIXING(ng) % Kv,                       &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                           BOUNDARY(ng) % b_t_obc,                &
     &                           BOUNDARY(ng) % b_u_obc,                &
     &                           BOUNDARY(ng) % b_v_obc,                &
# endif
     &                           BOUNDARY(ng) % b_ubar_obc,             &
     &                           BOUNDARY(ng) % b_vbar_obc,             &
     &                           BOUNDARY(ng) % b_zeta_obc,             &
#endif
#ifdef ADJUST_WSTRESS
     &                           FORCES(ng) % b_sustr,                  &
     &                           FORCES(ng) % b_svstr,                  &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                           FORCES(ng) % b_stflx,                  &
#endif
#ifdef SOLVE3D
     &                           OCEAN(ng) % b_t,                       &
     &                           OCEAN(ng) % b_u,                       &
     &                           OCEAN(ng) % b_v,                       &
#endif
     &                           OCEAN(ng) % b_zeta,                    &
     &                           OCEAN(ng) % b_ubar,                    &
     &                           OCEAN(ng) % b_vbar)
      END IF
!
      RETURN
      END SUBROUTINE normalization
!
!***********************************************************************
      SUBROUTINE normalization_tile (self, ng, tile,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nstp, nnew, ifac,                  &
     &                               Kh,                                &
#ifdef SOLVE3D
     &                               Kv,                                &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                               VnormRobc, VnormUobc, VnormVobc,   &
# endif
     &                               HnormRobc, HnormUobc, HnormVobc,   &
#endif
#ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                               HnormSTF,                          &
#endif
#ifdef SOLVE3D
     &                               VnormR, VnormU, VnormV,            &
#endif
     &                               HnormR, HnormU, HnormV)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, ifac
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Kh(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent(out) :: VnormUobc(LBij:,:,:)
      real(r8), intent(out) :: VnormVobc(LBij:,:,:)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:,:)
      real(r8), intent(out) :: HnormUobc(LBij:,:)
      real(r8), intent(out) :: HnormVobc(LBij:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:,:)
# endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:,:)

#else

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(out) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent(out) :: VnormVobc(LBij:UBij,N(ng),4)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormUobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormVobc(LBij:UBij,4)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng),NSA)
# endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj,NSA)
#endif
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      logical :: Ldiffer, Lsame
#endif
#ifdef ADJUST_BOUNDARY
      logical :: bounded
      logical :: Lconvolve(4)
#endif
      logical :: Lweak
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, ic, ifile, is, j, jc, rec
      integer :: ifield, nfield
      integer :: ns
#ifdef SOLVE3D
      integer :: UBt, itrc, k, kc, ntrc
#endif
#ifdef ADJUST_BOUNDARY
      integer :: Bmin, Bmax, IDmeta, IJlen, IJKlen
      integer :: ib, ibry, kb
!
      real(r8), parameter :: Aspv = 0.0_r8
#endif
      real(dp) :: my_time
      real(r8) :: cff, compute
      real(r8) :: my_dot, Gdotp
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
#endif
#ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
      real(r8), dimension(LBij:UBij) :: HscaleB
# ifdef SOLVE3D
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
      real(r8), dimension(LBij:UBij,1:N(ng)) :: VscaleB
#  ifdef DISTRIBUTE
      real(r8), dimension((UBij-LBij+1)*N(ng)) :: Bwrk
#  endif
# endif
#endif
!
#ifdef DISTRIBUTE
      character (len=3  ) :: op_handle
#endif
      character (len=40 ) :: Text
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", normalization_tile"

#if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
#endif

#include "set_bounds.h"
!
      SourceFile=MyFile

      my_time=tdays(ng)*day2sec

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8                   ! free surface
        END DO
      END DO

      CALL set_depth_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     SEDBED(ng) % bed_thick,                      &
# endif
     &                     A2d,                                         &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  The prior (initial conditions) and model normalization factors are
!  computed for modeling the spreading of the multiscale background-
!  error covariance matrix (B) using the EXACT METHOD. This method
!  involves convolving an implicit pseudo-diffusion operator within
!  the correlation-length function in K-space for each spatial grid
!  point independently.
!
!  Specifically, each point is perturbed with a delta function that is
!  scaled by the inverse square root of the area in two dimensions (2D)
!  or by the inverse square root of the volume in three dimensions (3D),
!  followed by spatial convolution.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      IF (Master) WRITE (stdout,10)

      FILE_LOOP : DO ifile=1,NSA

        LWRITE_NRM : IF (LwrtNRM(ifile,ng)) THEN
          Lweak=.FALSE.
          IF (ifile.eq.1) THEN
            Text='initial conditions'
          ELSE IF (ifile.eq.2) THEN
            Text='model'
            Lweak=.TRUE.
          END IF
!
!  Set time record index to write in normalization NetCDF file.
!
          ncname=NRM(ifile,ng)%name
          NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
          NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,idtime), my_time,           &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

#if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,idtime), my_time,       &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
#endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  2D normalization at RHO-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_R2D : IF (Cnorm(ifile,isFsur)) THEN
            MS_R2D_LOOP : DO ns=1,Nscale(ng)
              Imin=1
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                '2D normalization factors at RHO-points'
                FLUSH (stdout)
              END IF
!
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_r(i,j)*           &
     &                                    GRID(ng)%on_r(i,j))
                END DO
              END DO
!
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (GRID(ng)%rmask(ic,jc).gt.0) compute=1.0_r8
                  END IF
# ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
# endif
#else
                  compute=1.0_r8
#endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=1.0_r8
                    END IF
!
!  Apply lateral boundary conditions.
!
                    CALL ad_dabc_r2d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_2d (ng, tile, iADM, isFsur, r2dvar, &
     &                                  ns, NiterCI(ns,ng), ifac,       &
     &                                  Lweak,                          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  A2d)
!
                    DO j=JstrT,JendT
                      DO i=IstrT,IendT
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
!
!  Set normalization factors.
!
                    Gdotp=dot_prod2d (ng, tile, iADM, r2dvar,           &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d, A2d)
                    cff=1.0_r8/SQRT(Gdotp)
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    HnormR(ic,jc,ifile)=HnormR(ic,jc,ifile)+            &
     &                                  self%Bwgt(isFsur,ns)*cff
                  END IF
                END DO
              END DO
            END DO MS_R2D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormR(:,:,ifile))
#ifdef DISTRIBUTE
!
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormR(:,:,ifile))
#endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idFsur,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idFsur),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                GRID(ng)%rmask,                   &
#endif
     &                                HnormR(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idFsur)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idFsur,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idFsur),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               GRID(ng)%rmask,                    &
# endif
     &                               HnormR(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_R2D
!
!-----------------------------------------------------------------------
!  2D normalization at U-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_U2D : IF (Cnorm(ifile,isUbar)) THEN
            MS_U2D_LOOP : DO ns=1,Nscale(ng)
              IF (EWperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=2
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              END IF
!
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                '2D normalization factors at   U-points'
                FLUSH (stdout)
              END IF
!
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_u(i,j)*           &
     &                                    GRID(ng)%on_u(i,j))
                END DO
              END DO
!
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (GRID(ng)%umask(ic,jc).gt.0) compute=1.0_r8
                  END IF
# ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
# endif
#else
                  compute=1.0_r8
#endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=1.0_r8
                    END IF
!
!  Apply lateral boundary conditions.
!
                    CALL ad_dabc_u2d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_2d (ng, tile, iADM, isUbar, u2dvar, &
     &                                  ns, NiterCI(ns,ng), ifac,       &
     &                                  Lweak,                          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  A2d)
!
                    DO j=JstrT,JendT
                      DO i=IstrP,IendT
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
!
!  Set normalization factors.
!
                    Gdotp=dot_prod2d (ng, tile, iADM, u2dvar,           &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d, A2d)
                    cff=1.0_r8/SQRT(Gdotp)
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    HnormU(ic,jc,ifile)=HnormU(ic,jc,ifile)+            &
     &                                  self%Bwgt(isUbar,ns)*cff
                  END IF
                END DO
              END DO
            END DO MS_U2D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormU(:,:,ifile))
#ifdef DISTRIBUTE
!
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormU(:,:,ifile))
#endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idUbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idUbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                GRID(ng)%umask,                   &
#endif
     &                                HnormU(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idUbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idUbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               GRID(ng)%umask,                    &
# endif
     &                               HnormU(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_U2D
!
!-----------------------------------------------------------------------
!  2D normalization at V-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_V2D : IF (Cnorm(ifile,isVbar)) THEN
            MS_V2D_LOOP : DO ns=1,Nscale(ng)
              IF (NSperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=1
                Imax=Lm(ng)
                Jmin=2
                Jmax=Mm(ng)
              END IF
!
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                '2D normalization factors at   V-points'
                FLUSH (stdout)
              END IF
!
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_v(i,j)*           &
     &                                    GRID(ng)%on_v(i,j))
                END DO
              END DO
!
              DO jc=Jmin,Jmax
                DO ic=Imin,Imax
#ifdef MASKING
                  compute=0.0_r8
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    IF (GRID(ng)%vmask(ic,jc).gt.0) compute=1.0_r8
                  END IF
# ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
# endif
#else
                  compute=1.0_r8
#endif
                  IF (compute.gt.0.0_r8) THEN
                    DO j=LBj,UBj
                      DO i=LBi,UBi
                        A2d(i,j)=0.0_r8
                      END DO
                    END DO
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      A2d(ic,jc)=1.0_r8
                    END IF
!
!  Apply lateral boundary conditions.
!
                    CALL ad_dabc_v2d_tile (ng, tile,                    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_2d (ng, tile, iADM, isVbar, v2dvar, &
     &                                  ns, NiterCI(ns,ng), ifac,       &
     &                                  Lweak,                          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  A2d)
!
                    DO j=JstrP,JendT
                      DO i=IstrT,IendT
                        A2d(i,j)=A2d(i,j)*Hscale(i,j)
                      END DO
                    END DO
!
!  Set normalization factors.
!
                    Gdotp=dot_prod2d (ng, tile, iADM, v2dvar,           &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d, A2d)
                    cff=1.0_r8/SQRT(Gdotp)
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    HnormV(ic,jc,ifile)=HnormV(ic,jc,ifile)+            &
     &                                  self%Bwgt(isVbar,ns)*cff
                  END IF
                END DO
              END DO
            END DO MS_V2D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormV(:,:,ifile))
#ifdef DISTRIBUTE
!
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormV(:,:,ifile))
#endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idVbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idVbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                GRID(ng)%vmask,                   &
#endif
     &                                HnormV(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idVbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idVbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               GRID(ng)%vmask,                    &
# endif
     &                               HnormV(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_V2D

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D normalization U-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_U3D : IF (Cnorm(ifile,isUvel)) THEN
            MS_U3D_LOOP : DO ns=1,Nscale(ng)
              IF (EWperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=2
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              END IF
!
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                '3D normalization factors at   U-points'
                FLUSH (stdout)
              END IF
!
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  cff=GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+       &
     &                                      GRID(ng)%Hz(i  ,j,k)))
                  END DO
                END DO
              END DO
!
              DO kc=1,N(ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
# ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (GRID(ng)%umask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#  ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                    compute=1.0_r8
# endif
                    IF (compute.gt.0.0_r8) THEN
                      DO k=1,N(ng)
                        DO j=LBj,UBj
                          DO i=LBi,UBi
                            A3d(i,j,k)=0.0_r8
                          END DO
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A3d(ic,jc,kc)=1.0_r8
                      END IF
!
!  Implicit vertical convolution.
!
                      CALL self%ad_Vdiff_u3d (ng, tile, iADM, isUvel,   &
     &                                      NVsteps(ifile,isUvel)/ifac, &
     &                                        LBi, UBi, LBj, UBj,       &
     &                                        IminS, ImaxS,             &
     &                                        JminS, JmaxS,             &
     &                                        DTsizeV(ifile,isUvel),    &
     &                                        Kv, A3d)
!
!  Apply lateral boundary conditions.
!
                      CALL ad_dabc_u3d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                      CALL self%ad_CI_3d (ng, tile, iADM, isUvel,       &
     &                                    u3dvar,                       &
     &                                    ns, NiterCI(ns,ng), ifac,     &
     &                                    Lweak,                        &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    A3d)
!
                      DO k=1,N(ng)
                        DO j=JstrT,JendT
                          DO i=IstrP,IendT
                            A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                          END DO
                        END DO
                      END DO
!
!  Set normalization factors.
!
                      Gdotp=dot_prod3d (ng, tile, iADM, u3dvar,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  A3d, A3d)
                      cff=1.0_r8/SQRT(Gdotp)
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      VnormU(ic,jc,kc,ifile)=VnormU(ic,jc,kc,ifile)+    &
     &                                       self%Bwgt(isUvel,ns)*cff
                    END IF
                  END DO
                END DO
              END DO
            END DO MS_U3D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormU(:,:,:,ifile))
# ifdef DISTRIBUTE
!
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormU(:,:,:,ifile))
# endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idUvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idUvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                GRID(ng)%umask,                   &
# endif
     &                                VnormU(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idUvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idUvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               GRID(ng)%umask,                    &
#  endif
     &                               VnormU(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_U3D
!
!-----------------------------------------------------------------------
!  3D normalization at V-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_V3D : IF (Cnorm(ifile,isVvel)) THEN
            MS_V3D_LOOP : DO ns=1,Nscale(ng)
              IF (NSperiodic(ng)) THEN
                Imin=1
                Imax=Lm(ng)
                Jmin=1
                Jmax=Mm(ng)
              ELSE
                Imin=1
                Imax=Lm(ng)
                  Jmin=2
                Jmax=Mm(ng)
              END IF
!
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Text),                           &
     &                '3D normalization factors at   V-points'
                FLUSH (stdout)
              END IF
!
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  cff=GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+       &
     &                                      GRID(ng)%Hz(i,j  ,k)))
                  END DO
                END DO
              END DO
!
              DO kc=1,N(ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
# ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (GRID(ng)%vmask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#  ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                    compute=1.0_r8
# endif
                    IF (compute.gt.0.0_r8) THEN
                      DO k=1,N(ng)
                        DO j=LBj,UBj
                          DO i=LBi,UBi
                            A3d(i,j,k)=0.0_r8
                          END DO
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A3d(ic,jc,kc)=1.0_r8
                      END IF
!
!  Implicit vertical convolution.
!
                      CALL self%ad_Vdiff_v3d (ng, tile, iADM, isVvel,   &
     &                                      NVsteps(ifile,isVvel)/ifac, &
     &                                        LBi, UBi, LBj, UBj,       &
     &                                        IminS, ImaxS,             &
     &                                        JminS, JmaxS,             &
     &                                        DTsizeV(ifile,isVvel),    &
     &                                        Kv, A3d)
!
!  Apply lateral boundary conditions.
!
                      CALL ad_dabc_v3d_tile (ng, tile,                  &
     &                                       LBi, UBi, LBj, UBj,        &
     &                                       1, N(ng),                  &
     &                                       A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                      CALL self%ad_CI_3d (ng, tile, iADM, isVvel,       &
     &                                    v3dvar,                       &
     &                                    ns, NiterCI(ns,ng), ifac,     &
     &                                    Lweak,                        &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    A3d)
!
                      DO k=1,N(ng)
                        DO j=JstrP,JendT
                          DO i=IstrT,IendT
                            A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                          END DO
                        END DO
                      END DO
!
!  Set normalization factors.
!
                      Gdotp=dot_prod3d (ng, tile, iADM, v3dvar,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  A3d, A3d)
                      cff=1.0_r8/SQRT(Gdotp)
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      VnormV(ic,jc,kc,ifile)=VnormV(ic,jc,kc,ifile)+    &
     &                                       self%Bwgt(isVvel,ns)*cff
                    END IF
                  END DO
                END DO
              END DO
            END DO MS_V3D_LOOP
!
!  Exchange bounndary data.
!
            CALL dabc_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormV(:,:,:,ifile))
# ifdef DISTRIBUTE
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormV(:,:,:,ifile))
# endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idVvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idVvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                GRID(ng)%vmask,                   &
# endif
     &                                VnormV(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idVvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idVvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               GRID(ng)%vmask,                    &
#  endif
     &                               VnormV(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_V3D
!
!-----------------------------------------------------------------------
!  3D normalization at RHO-points.
!-----------------------------------------------------------------------
!
          IF (Master) THEN
            Lsame=.FALSE.
            DO itrc=1,NT(ng)
              ifield=isTvar(itrc)
              IF (Cnorm(ifile,ifield)) Lsame=.TRUE.
            END DO
            IF (Lsame) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '3D normalization factors at RHO-points'
              FLUSH (stdout)
            END IF
          END IF
!
!  Check if the decorrelation scales for all the tracers are different.
!  If not, just compute the normalization factors for the first tracer
!  and assign the same value to the rest.  Recall that this computation
!  is very expensive.
!
          Ldiffer=.FALSE.
          DO ns=1,Nscale(ng)
            Imin=1
            Imax=Lm(ng)
            Jmin=1
            Jmax=Mm(ng)

# ifdef NONUNIFORM_SCALES
            Ldiffer=.TRUE.
# else
            DO itrc=2,NT(ng)
              IF ((HdecayX(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayX(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (HdecayY(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayY(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (Vdecay(ifile,isTvar(itrc  ),ng).ne.                  &
     &             Vdecay(ifile,isTvar(itrc-1),ng))) THEN
                Ldiffer=.TRUE.
              END IF
            END DO
# endif
            IF (.not.Ldiffer) THEN
              Lsame=.TRUE.
              UBt=1
            ELSE
              Lsame=.FALSE.
              UBt=NT(ng)
            END IF
          END DO
!
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              cff=GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j)
              DO k=1,N(ng)
                Vscale(i,j,k)=1.0_r8/SQRT(cff*GRID(ng)%Hz(i,j,k))
              END DO
            END DO
          END DO
!
          TRACER_LOOP : DO itrc=1,UBt
            ifield=isTvar(itrc)
            MS_COMPUTE_R3D : IF (Cnorm(ifile,ifield)) THEN
              MS_R3D_LOOP : DO ns=1,Nscale(ng)
!
                DO kc=1,N(ng)
                  DO jc=Jmin,Jmax
                    DO ic=Imin,Imax
# ifdef MASKING
                      compute=0.0_r8
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        IF (GRID(ng)%rmask(ic,jc).gt.0) compute=1.0_r8
                      END IF
#  ifdef DISTRIBUTE
                      CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#  endif
# else
                      compute=1.0_r8
# endif
                      IF (compute.gt.0.0_r8) THEN
                        DO k=1,N(ng)
                          DO j=LBj,UBj
                            DO i=LBi,UBi
                              A3d(i,j,k)=0.0_r8
                            END DO
                          END DO
                        END DO
                        IF (((Jstr.le.jc).and.(jc.le.Jend)).and.        &
     &                      ((Istr.le.ic).and.(ic.le.Iend))) THEN
                          A3d(ic,jc,kc)=1.0_r8
                        END IF
!
!  Implicit vertical convolution.
!
                        CALL self%ad_Vdiff_r3d (ng, tile, iADM, ifield, &
     &                                      NVsteps(ifile,ifield)/ifac, &
     &                                          LBi, UBi, LBj, UBj,     &
     &                                          IminS, ImaxS,           &
     &                                          JminS, JmaxS,           &
     &                                          DTsizeV(ifile,ifield),  &
     &                                          Kv, A3d)
!
!  Apply lateral boundary conditions.
!
                        CALL ad_dabc_r3d_tile (ng, tile,                &
     &                                         LBi, UBi, LBj, UBj,      &
     &                                         1, N(ng),                &
     &                                         A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                        CALL self%ad_CI_3d (ng, tile, iADM, ifield,     &
     &                                      r3dvar,                     &
     &                                      ns, NiterCI(ns,ng), ifac,   &
     &                                      Lweak,                      &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      IminS, ImaxS, JminS, JmaxS, &
     &                                      A3d)
!
                        DO k=1,N(ng)
                          DO j=JstrT,JendT
                            DO i=IstrT,IendT
                              A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                            END DO
                          END DO
                        END DO
!
!  Set normalization factors.
!
                        Gdotp=dot_prod3d (ng, tile, iADM, r3dvar,       &
     &                                    LBi, UBi, LBj, UBj, 1, N(ng), &
     &                                    A3d, A3d)
                        cff=1.0_r8/SQRT(Gdotp)
                      ELSE
                        cff=0.0_r8
                      END IF
!
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        IF (Lsame) THEN
                          DO ntrc=1,NT(ng)
                            VnormR(ic,jc,kc,ifile,ntrc)=                &
     &                                   VnormR(ic,jc,kc,ifile,ntrc)+   &
     &                                   self%Bwgt(isTvar(ntrc),ns)*cff
                          END DO
                        ELSE
                          VnormR(ic,jc,kc,ifile,itrc)=                  &
     &                                   VnormR(ic,jc,kc,ifile,itrc)+   &
     &                                   self%Bwgt(ifield,ns)*cff
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END DO MS_R3D_LOOP
            END IF MS_COMPUTE_R3D
          END DO TRACER_LOOP
!
!  Exchange boundary data.
!
          DO itrc=1,NT(ng)
            ifield=isTvar(itrc)
            IF (Cnorm(ifile,ifield)) THEN
              CALL dabc_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            VnormR(:,:,:,ifile,itrc))
# ifdef DISTRIBUTE
!
              CALL mp_exchange3d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            VnormR(:,:,:,ifile,itrc))
# endif
!
!  Write out into output NetCDF file.
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  idTvar(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTvar(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
# ifdef MASKING
     &                                  GRID(ng)%rmask,                 &
# endif
     &                                  VnormR(:,:,:,ifile,itrc))

# if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioTrc(itrc)%dkind.eq.              &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r3dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r3dvar(ng)
                  END IF
                  CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 idTvar(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioTrc(itrc),         &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                  ioDesc,                         &
#  ifdef MASKING
     &                                  GRID(ng)%rmask,                 &
#  endif
     &                                  VnormR(:,:,:,ifile,itrc))
# endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#endif
        END IF LWRITE_NRM
      END DO FILE_LOOP

#ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute open boundaries error covariance, B, normalization factors
!  using the exact method.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      ifile=3
      BRY_LWRITE_NRM : IF (LwrtNRM(ifile,ng)) THEN
        Text='boundary conditions'
        IJlen=UBij-LBij+1
# ifdef SOLVE3D
        IJKlen=IJlen*N(ng)
# endif
        Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
        Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
        Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
        Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        my_time=tdays(ng)*day2sec

        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  2D boundary normalization at RHO-points.
!-----------------------------------------------------------------------
!
        HnormRobc=Aspv

        IF (Master.and.ANY(CnormB(isFsur,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '2D normalization factors at RHO-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_R2D : DO ibry=1,4
          HscaleB=0.0_r8
          BRY_COMPUTE_R2D : IF (CnormB(isFsur,ibry)) THEN
            BRY_MS_LOOP_R2D : DO ns=1,Nscale(ng)
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,r2dvar)
                  Bmin=1
                  Bmax=Mm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      HscaleB(j)=1.0_r8/SQRT(GRID(ng)%on_r(i,j))
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,r2dvar)
                  Bmin=1
                  Bmax=Lm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      HscaleB(i)=1.0_r8/SQRT(GRID(ng)%om_r(i,j))
                    END DO
                  END IF
              END SELECT
!
              DO ib=Bmin,Bmax
!
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Jstr.le.ib).and.(ib.le.Jend))
                    i=BOUNDS(ng)%edge(ibry,r2dvar)
                    j=ib
                  CASE (isouth, inorth)
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Istr.le.ib).and.(ib.le.Iend))
                    i=ib
                    j=BOUNDS(ng)%edge(ibry,r2dvar)
                END SELECT
#  ifdef MASKING
                IF (bounded) THEN
                  compute=GRID(ng)%rmask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#   ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                compute=1.0_r8
#  endif
                IF (compute.gt.0.0_r8) THEN
                  B2d=0.0_r8
                  IF (bounded) THEN
                    B2d(ib)=1.0_r8
                  END IF
!
!  Implicit adjoint convolution, CG/CI solver.
!
                  CALL self%ad_CI_b1d (ng, tile, iADM, isFsur, ibry,    &
     &                                 r2dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B2d)
!
!  HscaleB must be applied twice.
!
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO j=JstrT,JendT
                        B2d(j)=B2d(j)*HscaleB(j)
                      END DO
!
                      DO j=JstrT,JendT
                        B2d(j)=B2d(j)*HscaleB(j)
                      END DO
                    CASE (isouth, inorth)
                      DO i=IstrT,IendT
                        B2d(i)=B2d(i)*HscaleB(i)
                      END DO
!
                      DO i=IstrT,IendT
                        B2d(i)=B2d(i)*HscaleB(i)
                      END DO
                  END SELECT
!
!  Implicit tangent linear convolution, CG/CI solver.
!
                  CALL self%tl_CI_b1d (ng, tile, iTLM, isFsur, ibry,    &
     &                                 r2dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B2d)
!
!  Set normalization factors.
!
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B2d(ib))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  HnormRobc(ib,ibry)=HnormRobc(ib,ibry)+                &
     &                               self%Bwgt(isFsur,ns)*cff
                END IF
              END DO
            END DO BRY_MS_LOOP_R2D
!
!  Exchange boundary data.
!
            CALL bc_r2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormRobc(:,ibry))
# ifdef DISTRIBUTE
!
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormRobc(LBij:,ibry))
# endif
          END IF BRY_COMPUTE_R2D
        END DO BRY_LOOP_R2D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isFsur,:))) THEN
          IDmeta=idSbry(isFsur)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              HnormRobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  HnormRobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)

# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D boundary normalization at U-points.
!-----------------------------------------------------------------------
!
        HnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '2D normalization factors at   U-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_U2D : DO ibry=1,4
          HscaleB=0.0_r8
          BRY_COMPUTE_U2D : IF (CnormB(isUbar,ibry)) THEN
            BRY_MS_LOOP_U2D : DO ns=1,Nscale(ng)
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,u2dvar)
                  Bmin=1
                  Bmax=Mm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      HscaleB(j)=1.0_r8/SQRT(GRID(ng)%on_u(i,j))
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,u2dvar)
                  IF (EWperiodic(ng)) THEN
                    Bmin=1
                    Bmax=Lm(ng)
                  ELSE
                    Bmin=2
                    Bmax=Lm(ng)
                  END IF
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrP,IendT
                      HscaleB(i)=1.0_r8/SQRT(GRID(ng)%om_u(i,j))
                    END DO
                  END IF
              END SELECT
!
              DO ib=Bmin,Bmax
!
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Jstr.le.ib).and.(ib.le.Jend))
                    i=BOUNDS(ng)%edge(ibry,u2dvar)
                    j=ib
                  CASE (isouth, inorth)
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Istr.le.ib).and.(ib.le.Iend))
                    i=ib
                    j=BOUNDS(ng)%edge(ibry,u2dvar)
                END SELECT
#  ifdef MASKING
                IF (bounded) THEN
                  compute=GRID(ng)%umask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#   ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                compute=1.0_r8
#  endif
                IF (compute.gt.0.0_r8) THEN
                  B2d=0.0_r8
                  IF (bounded) THEN
                    B2d(ib)=1.0_r8
                  END IF
!
!  Implicit adjoint convolution, CG/CI solver.
!
                  CALL self%ad_CI_b1d (ng, tile, iADM, isUbar, ibry,    &
     &                                 u2dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B2d)
!
!  HscaleB must be applied twice.
!
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO j=JstrT,JendT
                        B2d(j)=B2d(j)*HscaleB(j)
                      END DO
!
                      DO j=JstrT,JendT
                        B2d(j)=B2d(j)*HscaleB(j)
                      END DO
                    CASE (isouth, inorth)
                      DO i=IstrP,IendT
                        B2d(i)=B2d(i)*HscaleB(i)
                      END DO
!
                      DO i=IstrP,IendT
                        B2d(i)=B2d(i)*HscaleB(i)
                      END DO
                  END SELECT
!
!  Implicit tangent linear convolution, CG/CI solver.
!
                  CALL self%tl_CI_b1d (ng, tile, iTLM, isUbar, ibry,    &
     &                                 u2dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B2d)
!
!  Set normalization factors.
!
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B2d(ib))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  HnormUobc(ib,ibry)=HnormUobc(ib,ibry)+                &
     &                               self%Bwgt(isUbar,ns)*cff
                END IF
              END DO
            END DO BRY_MS_LOOP_U2D
!
!  Exchange boundary data.
!
            CALL bc_u2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormUobc(:,ibry))
# ifdef DISTRIBUTE
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormUobc(LBij:,ibry))
# endif
          END IF BRY_COMPUTE_U2D
        END DO BRY_LOOP_U2D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isUbar,:))) THEN
          IDmeta=idSbry(isUbar)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              HnormUobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  HnormUobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D boundary normalization at V-points.
!-----------------------------------------------------------------------
!
        HnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '2D normalization factors at   V-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_V2D : DO ibry=1,4
          HscaleB=0.0_r8
          BRY_COMPUTE_V2D : IF (CnormB(isVbar,ibry)) THEN
            BRY_MS_LOOP_V2D : DO ns=1,Nscale(ng)
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,v2dvar)
                  IF (NSperiodic(ng)) THEN
                    Bmin=1
                    Bmax=Mm(ng)
                  ELSE
                    Bmin=2
                    Bmax=Mm(ng)
                  END IF
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      HscaleB(j)=1.0_r8/SQRT(GRID(ng)%on_v(i,j))
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,v2dvar)
                  Bmin=1
                  Bmax=Lm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      HscaleB(i)=1.0_r8/SQRT(GRID(ng)%om_v(i,j))
                    END DO
                  END IF
              END SELECT
!
              DO ib=Bmin,Bmax
!
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Jstr.le.ib).and.(ib.le.Jend))
                    i=BOUNDS(ng)%edge(ibry,v2dvar)
                    j=ib
                  CASE (isouth, inorth)
                    bounded=Lconvolve(ibry).and.                        &
     &                      ((Istr.le.ib).and.(ib.le.Iend))
                    i=ib
                    j=BOUNDS(ng)%edge(ibry,v2dvar)
                END SELECT
#  ifdef MASKING
                IF (bounded) THEN
                  compute=GRID(ng)%vmask(i,j)
                ELSE
                  compute=0.0_r8
                END IF
#   ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#   endif
#  else
                compute=1.0_r8
#  endif
                IF (compute.gt.0.0_r8) THEN
                  B2d=0.0_r8
                  IF (bounded) THEN
                    B2d(ib)=1.0_r8
                  END IF
!
!  Implicit adjoint convolution, CG/CI solver.
!
                  CALL self%ad_CI_b1d (ng, tile, iADM, isVbar, ibry,    &
     &                                 v2dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B2d)
!
!  HscaleB must be applied twice.
!
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO j=JstrP,JendT
                        B2d(j)=B2d(j)*HscaleB(j)
                      END DO
!
                      DO j=JstrP,JendT
                        B2d(j)=B2d(j)*HscaleB(j)
                      END DO
                    CASE (isouth, inorth)
                      DO i=IstrT,IendT
                        B2d(i)=B2d(i)*HscaleB(i)
                      END DO
!
                      DO i=IstrT,IendT
                        B2d(i)=B2d(i)*HscaleB(i)
                      END DO
                  END SELECT
!
!  Implicit tangent linear convolution, CG/CI solver.
!
                  CALL self%tl_CI_b1d (ng, tile, iTLM, isVbar, ibry,    &
     &                                 v2dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B2d)
!
!  Set normalization factors.
!
                  IF (bounded) THEN
                    cff=1.0_r8/SQRT(B2d(ib))
                  END IF
                ELSE
                  cff=0.0_r8
                END IF
                IF (bounded) THEN
                  HnormVobc(ib,ibry)=HnormVobc(ib,ibry)+                &
     &                               self%Bwgt(isVbar,ns)*cff
                END IF
              END DO
            END DO BRY_MS_LOOP_V2D
!
!  Exchange boundary data.
!
            CALL bc_v2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormVobc(:,ibry))
# ifdef DISTRIBUTE
!
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormVobc(LBij:,ibry))
# endif
          END IF BRY_COMPUTE_V2D
        END DO BRY_LOOP_V2D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isVbar,:))) THEN
          IDmeta=idSbry(isVbar)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              HnormVobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  HnormVobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D boundary normalization at U-points.
!-----------------------------------------------------------------------
!
        VnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '3D normalization factors at   U-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_U3D : DO ibry=1,4
          VscaleB=0.0_r8
          BRY_COMPUTE_U3D : IF (CnormB(isUvel,ibry)) THEN
            BRY_MS_LOOP_U3D : DO ns=1,Nscale(ng)
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,u2dvar)
                  Bmin=1
                  Bmax=Mm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      cff=GRID(ng)%on_u(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(j,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+    &
     &                                         GRID(ng)%Hz(i  ,j,k)))
                      END DO
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,u2dvar)
                  IF (EWperiodic(ng)) THEN
                    Bmin=1
                    Bmax=Lm(ng)
                  ELSE
                    Bmin=2
                    Bmax=Lm(ng)
                  END IF
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrP,IendT
                      cff=GRID(ng)%om_u(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(i,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+    &
     &                                         GRID(ng)%Hz(i  ,j,k)))
                      END DO
                    END DO
                  END IF
              END SELECT
!
              DO kb=1,N(ng)
                DO ib=Bmin,Bmax
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      bounded=Lconvolve(ibry).and.                      &
     &                        ((Jstr.le.ib).and.(ib.le.Jend))
                      i=BOUNDS(ng)%edge(ibry,u2dvar)
                      j=ib
                    CASE (isouth, inorth)
                      bounded=Lconvolve(ibry).and.                      &
     &                        ((Istr.le.ib).and.(ib.le.Iend))
                      i=ib
                      j=BOUNDS(ng)%edge(ibry,u2dvar)
                  END SELECT
#   ifdef MASKING
                  IF (bounded) THEN
                    compute=GRID(ng)%umask(i,j)
                  ELSE
                    compute=0.0_r8
                  END IF
#    ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                  compute=1.0_r8
#   endif
                  IF (compute.gt.0.0_r8) THEN
                    B3d=0.0_r8
                    IF (bounded) THEN
                      B3d(ib,kb)=1.0_r8
                    END IF
!
!  First, implicit adjoint vertical convolution.
!
                    CALL self%ad_bry_Vdiff (ng, tile, iADM,             &
     &                                      isUvel, ibry, u3dvar,       &
     &                                      NVstepsB(ibry,isUvel)/ifac, &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      LBij, UBij,                 &
     &                                      DTsizeVB(ibry,isUvel), Kv,  &
     &                                      B3d)
!
!  Now, implicit adjoint horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_b2d (ng, tile, iADM,                &
     &                                   isUvel, ibry, u3dvar,          &
     &                                   ns, NiterCI(ns,ng), ifac,      &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   B3d)
!
!  VscaleB must be applied twice.
!
                    SELECT CASE (ibry)
                      CASE (iwest, ieast)
                        DO k=1,N(ng)
                          DO j=JstrT,JendT
                            B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                          END DO
                        END DO
!
                        DO k=1,N(ng)
                          DO j=JstrT,JendT
                            B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                          END DO
                        END DO
                      CASE (isouth, inorth)
                        DO k=1,N(ng)
                          DO i=IstrP,IendT
                            B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                          END DO
                        END DO
!
                        DO k=1,N(ng)
                          DO i=IstrP,IendT
                            B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                          END DO
                        END DO
                    END SELECT
!
!  First, implicit tangent linear horizontal convolution, CG/CI solver.
!
                    CALL self%tl_CI_b2d (ng, tile, iTLM,                &
     &                                   isUvel, ibry, u3dvar,          &
     &                                   ns, NiterCI(ns,ng), ifac,      &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   B3d)
!
!  Now, implicit tangent linear vertical convolution.
!
                    CALL self%tl_bry_Vdiff (ng, tile, iTLM,             &
     &                                      isUvel, ibry, u3dvar,       &
     &                                      NVstepsB(ibry,isUvel)/ifac, &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      LBij, UBij,                 &
     &                                      DTsizeVB(ibry,isUvel), Kv,  &
     &                                      B3d)
!
!  Set normalization factors.
!
                    IF (bounded) THEN
                      cff=1.0_r8/SQRT(B3d(ib,kb))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (bounded) THEN
                    VnormUobc(ib,kb,ibry)=VnormUobc(ib,kb,ibry)+        &
     &                                    self%Bwgt(isUvel,ns)*cff
                  END IF
                END DO
              END DO
            END DO BRY_MS_LOOP_U3D
!
!  Exchange boundary data.
!
            CALL bc_u3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormUobc(:,:,ibry))

#  ifdef DISTRIBUTE
!
            Bwrk=RESHAPE(VnormUobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormUobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF BRY_COMPUTE_U3D
        END DO BRY_LOOP_U3D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isUvel,:))) THEN
          IDmeta=idSbry(isUvel)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              VnormUobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  VnormUobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D boundary normalization at V-points.
!-----------------------------------------------------------------------
!
        VnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '3D normalization factors at   V-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_V3D : DO ibry=1,4
          VscaleB=0.0_r8
          BRY_COMPUTE_V3D : IF (CnormB(isVvel,ibry)) THEN
            BRY_MS_LOOP_V3D : DO ns=1,Nscale(ng)
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,v2dvar)
                  IF (NSperiodic(ng)) THEN
                    Bmin=1
                    Bmax=Mm(ng)
                  ELSE
                    Bmin=2
                    Bmax=Mm(ng)
                  END IF
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrP,JendT
                      cff=GRID(ng)%on_v(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(j,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+    &
     &                                         GRID(ng)%Hz(i,j  ,k)))
                      END DO
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,v2dvar)
                  Bmin=1
                  Bmax=Lm(ng)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      cff=GRID(ng)%om_v(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(i,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+    &
     &                                         GRID(ng)%Hz(i,j,  k)))
                      END DO
                    END DO
                  END IF
              END SELECT
!
              DO kb=1,N(ng)
                DO ib=Bmin,Bmax
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      bounded=Lconvolve(ibry).and.                      &
     &                        ((Jstr.le.ib).and.(ib.le.Jend))
                      i=BOUNDS(ng)%edge(ibry,v2dvar)
                      j=ib
                    CASE (isouth, inorth)
                      bounded=Lconvolve(ibry).and.                      &
     &                        ((Istr.le.ib).and.(ib.le.Iend))
                      i=ib
                      j=BOUNDS(ng)%edge(ibry,v2dvar)
                  END SELECT
#   ifdef MASKING
                  IF (bounded) THEN
                    compute=GRID(ng)%vmask(i,j)
                  ELSE
                    compute=0.0_r8
                  END IF
#    ifdef DISTRIBUTE
                  CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                  compute=1.0_r8
#   endif
                  IF (compute.gt.0.0_r8) THEN
                    B3d=0.0_r8
                    IF (bounded) THEN
                      B3d(ib,kb)=1.0_r8
                    END IF
!
!  First, implicit adjoint vertical convolution.
!
                    CALL self%ad_bry_Vdiff (ng, tile, iADM,             &
     &                                      isVvel, ibry, v3dvar,       &
     &                                      NVstepsB(ibry,isUvel)/ifac, &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      LBij, UBij,                 &
     &                                      DTsizeVB(ibry,isVvel), Kv,  &
     &                                      B3d)
!
!  Now, implicit adjoint horizontal convolution, CG/CI solver.
!
                    CALL self%ad_CI_b2d (ng, tile, iADM,                &
     &                                   isVvel, ibry, v3dvar,          &
     &                                   ns, NiterCI(ns,ng), ifac,      &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   B3d)
!
!  VscaleB must be applied twice.
!
                    SELECT CASE (ibry)
                      CASE (iwest, ieast)
                        DO k=1,N(ng)
                          DO j=JstrP,JendT
                            B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                          END DO
                        END DO
!
                        DO k=1,N(ng)
                          DO j=JstrP,JendT
                            B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                          END DO
                        END DO
                      CASE (isouth, inorth)
                        DO k=1,N(ng)
                          DO i=IstrT,IendT
                            B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                          END DO
                        END DO
!
                        DO k=1,N(ng)
                          DO i=IstrT,IendT
                            B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                          END DO
                        END DO
                    END SELECT
!
!  First, implicit tangent linear horizontal convolution, CG/CI solver.
!
                    CALL self%tl_CI_b2d (ng, tile, iTLM,                &
     &                                   isVvel, ibry, v3dvar,          &
     &                                   ns, NiterCI(ns,ng), ifac,      &
     &                                   LBij, UBij,                    &
     &                                   IminS, ImaxS, JminS, JmaxS,    &
     &                                   B3d)
!
!  Now, implicit tangent linear vertical convolution.
!
                    CALL self%tl_bry_Vdiff (ng, tile, iTLM,             &
     &                                      isVvel, ibry, v3dvar,       &
     &                                      NVstepsB(ibry,isUvel)/ifac, &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      LBij, UBij,                 &
     &                                      DTsizeVB(ibry,isVvel), Kv,  &
     &                                      B3d)
!
!  Set normalization factors.
!
                    IF (bounded) THEN
                      cff=1.0_r8/SQRT(B3d(ib,kb))
                    END IF
                  ELSE
                    cff=0.0_r8
                  END IF
                  IF (bounded) THEN
                    VnormVobc(ib,kb,ibry)=VnormVobc(ib,kb,ibry)+        &
     &                                    self%Bwgt(isVvel,ns)*cff
                  END IF
                END DO
              END DO
            END DO BRY_MS_LOOP_V3D
!
!  Exchange boundary data.
!
            CALL bc_v3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormVobc(:,:,ibry))

#  ifdef DISTRIBUTE
!
            Bwrk=RESHAPE(VnormVobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormVobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF BRY_COMPUTE_V3D
        END DO BRY_LOOP_V3D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isVvel,:))) THEN
          IDmeta=idSbry(isVvel)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              VnormVobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  VnormVobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D boundary normalization at RHO-points.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          DO itrc=1,NT(ng)
            ifield=isTvar(itrc)
            IF (ANY(CnormB(ifield,:))) THEN
              Lsame=.TRUE.
              EXIT
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &            '3D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF
!
        BRY_TRACER_LOOP : DO itrc=1,NT(ng)
          VnormRobc=Aspv
          ifield=isTvar(itrc)
          BRY_LOOP_R3D : DO ibry=1,4
            VscaleB=0.0_r8
            BRY_COMPUTE_R3D : IF (CnormB(ifield,ibry)) THEN
              BRY_MS_LOOP_R3D : DO ns=1,Nscale(ng)
!
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,r2dvar)
                    Bmin=1
                    Bmax=Mm(ng)
                    IF (Lconvolve(ibry)) THEN
                      DO j=JstrT,JendT
                        cff=GRID(ng)%on_r(i,j)
                        DO k=1,N(ng)
                          VscaleB(j,k)=1.0_r8/                          &
     &                                 SQRT(cff*GRID(ng)%Hz(i,j,k))
                        END DO
                      END DO
                    END IF
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,r2dvar)
                    Bmin=1
                    Bmax=Lm(ng)
                    IF (Lconvolve(ibry)) THEN
                      DO i=IstrT,IendT
                        cff=GRID(ng)%om_r(i,j)
                        DO k=1,N(ng)
                          VscaleB(i,k)=1.0_r8/                          &
     &                                 SQRT(cff*GRID(ng)%Hz(i,j,k))
                        END DO
                      END DO
                    END IF
                END SELECT
!
                DO kb=1,N(ng)
                  DO ib=Bmin,Bmax
                    SELECT CASE (ibry)
                      CASE (iwest, ieast)
                        bounded=Lconvolve(ibry).and.                    &
     &                          ((Jstr.le.ib).and.(ib.le.Jend))
                        i=BOUNDS(ng)%edge(ibry,r2dvar)
                        j=ib
                      CASE (isouth, inorth)
                        bounded=Lconvolve(ibry).and.                    &
     &                          ((Istr.le.ib).and.(ib.le.Iend))
                        i=ib
                        j=BOUNDS(ng)%edge(ibry,r2dvar)
                    END SELECT
#   ifdef MASKING
                    IF (bounded) THEN
                      compute=GRID(ng)%rmask(i,j)
                    ELSE
                      compute=0.0_r8
                    END IF
#    ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                    compute=1.0_r8
#   endif
                    IF (compute.gt.0.0_r8) THEN
                      B3d=0.0_r8
                      IF (bounded) THEN
                        B3d(ib,kb)=1.0_r8
                      END IF
!
!  First, implicit adjoint vertical convolution.
!
                      CALL self%ad_bry_Vdiff (ng, tile, iADM,           &
     &                                        ifield, ibry, r3dvar,     &
     &                                      NVstepsB(ibry,ifield)/ifac, &
     &                                        LBi, UBi, LBj, UBj,       &
     &                                        LBij, UBij,               &
     &                                        DTsizeVB(ibry,ifield), Kv,&
     &                                        B3d)
!
!  Now, implicit adjoint horizontal convolution, CG/CI solver.
!
                      CALL self%ad_CI_b2d (ng, tile, iADM,              &
     &                                     ifield, ibry, r3dvar,        &
     &                                     ns, NiterCI(ns,ng), ifac,    &
     &                                     LBij, UBij,                  &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     B3d)
!
!  VscaleB must be applied twice.
!
                      SELECT CASE (ibry)
                        CASE (iwest, ieast)
                          DO k=1,N(ng)
                            DO j=JstrT,JendT
                              B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                            END DO
                          END DO
!
                          DO k=1,N(ng)
                            DO j=JstrT,JendT
                              B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                            END DO
                          END DO
                        CASE (isouth, inorth)
                          DO k=1,N(ng)
                            DO i=IstrT,IendT
                              B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                            END DO
                          END DO
!
                          DO k=1,N(ng)
                            DO i=IstrT,IendT
                              B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                            END DO
                          END DO
                      END SELECT
!
!  First, implicit tangent linear horizontal convolution, CG/CI solver.
!
                      CALL self%tl_CI_b2d (ng, tile, iTLM,              &
     &                                     ifield, ibry, r3dvar,        &
     &                                     ns, NiterCI(ns,ng), ifac,    &
     &                                     LBij, UBij,                  &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     B3d)
!
!  Now, implicit tangent linear vertical convolution.
!
                      CALL self%tl_bry_Vdiff (ng, tile, iTLM,           &
     &                                        ifield, ibry, r3dvar,     &
     &                                      NVstepsB(ibry,ifield)/ifac, &
     &                                        LBi, UBi, LBj, UBj,       &
     &                                        LBij, UBij,               &
     &                                        DTsizeVB(ibry,ifield), Kv,&
     &                                        B3d)
!
!  Set normalization factors.
!
                      IF (bounded) THEN
                        cff=1.0_r8/SQRT(B3d(ib,kb))
                      END IF
                    ELSE
                      cff=0.0_r8
                    END IF
                    IF (bounded) THEN
                      VnormRobc(ib,kb,ibry,itrc)=                       &
     &                                     VnormRobc(ib,kb,ibry,itrc)+  &
     &                                     self%Bwgt(ifield,ns)*cff
                    END IF
                  END DO
                END DO
              END DO BRY_MS_LOOP_R3D
!
!  Exchange boundary data.
!
              CALL bc_r3d_bry_tile (ng, tile, ibry,                     &
     &                              LBij, UBij, 1, N(ng),               &
     &                              VnormRobc(:,:,ibry,itrc))

#  ifdef DISTRIBUTE
!
              Bwrk=RESHAPE(VnormRobc(:,:,ibry,itrc), (/IJKlen/))
              CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
              ic=0
              DO k=1,N(ng)
                DO ib=LBij,UBij
                  ic=ic+1
                  VnormRobc(ib,k,ibry,itrc)=Bwrk(ic)
                END DO
              END DO
#  endif
            END IF BRY_COMPUTE_R3D
          END DO BRY_LOOP_R3D
!
!  Write out into output NetCDF file.
!
          IF (ANY(CnormB(ifield,:))) THEN
            IDmeta=idSbry(isTvar(itrc))
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,IDmeta),                  &
     &                                VnormRobc(LBij:,:,:,itrc),        &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

#  if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_put_fvar (ng, iTLM, ncname,             &
     &                                    Vname(1,IDmeta),              &
     &                                    VnormRobc(LBij:,:,:,itrc),    &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
#  endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

        END DO BRY_TRACER_LOOP

# endif /* SOLVE3D */
!
!  Synchronize open boundaries normalization NetCDF file to disk to
!  allow other processes to access data immediately after it is
!  written.
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_sync (ng, iTLM, ncname,                         &
     &                        NRM(ifile,ng)%ncid)
# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_sync (ng, iTLM, ncname,                     &
     &                            NRM(ifile,ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF BRY_LWRITE_NRM

#endif /* ADJUST_BOUNDARY */

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute surface forcing error covariance, B, normalization factors
!  using the exact method.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      ifile=4
      FRC_LWRITE_NRM : IF (LwrtNRM(ifile,ng)) THEN
        rec=1
        Text='surface forcing'
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        my_time=tdays(ng)*day2sec

        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  2D normalization at U-stress points.
!-----------------------------------------------------------------------
!
        MS_COMPUTE_SUS : IF (Cnorm(rec,isUstr)) THEN
          MS_SUS_LOOP : DO ns=1,Nscale(ng)
            IF (EWperiodic(ng)) THEN
              Imin=1
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
            ELSE
              Imin=2
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
            END IF
!
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '2D normalization factors at U-stress points'
              FLUSH (stdout)
            END IF
!
            DO j=JstrT,JendT
              DO i=IstrP,IendT
                Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_u(i,j)*             &
     &                                  GRID(ng)%on_u(i,j))
              END DO
            END DO
!
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
#   ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (GRID(ng)%umask(ic,jc).gt.0) compute=1.0_r8
                END IF
#    ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                compute=1.0_r8
#   endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=1.0_r8
                  END IF
!
!  Implicit horizontal convolution, CG/CI solver.
!
                  CALL self%ad_CI_2d (ng, tile, iADM, isUstr, u2dvar,   &
     &                                ns, NiterCI(ns,ng), ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                A2d)
!
                  DO j=JstrT,JendT
                    DO i=IstrP,IendT
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
!
!  Set normalization factors.
!
                  Gdotp=dot_prod2d (ng, tile, iADM, u2dvar,             &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A2d, A2d)
                  cff=1.0_r8/SQRT(Gdotp)
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormSUS(ic,jc)=HnormSUS(ic,jc)+                      &
     &                            self%Bwgt(isUstr,ns)*cff
                END IF
              END DO
            END DO
          END DO MS_SUS_LOOP
!
!  Exchange boundary data.
!
          CALL dabc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSUS)
#  ifdef DISTRIBUTE
!
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSUS)
#  endif
!
!  Write out into output NetCDF file.
!
         SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idUsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idUsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              GRID(ng)%umask,                     &
#  endif
     &                              HnormSUS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idUsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idUsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idUsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             GRID(ng)%umask,                      &
#   endif
     &                             HnormSUS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF MS_COMPUTE_SUS
!
!-----------------------------------------------------------------------
!  2D normalization at V-stress points.
!-----------------------------------------------------------------------
!
        MS_COMPUTE_SVS : IF (Cnorm(rec,isVstr)) THEN
          MS_SVS_LOOP : DO ns=1,Nscale(ng)
            IF (NSperiodic(ng)) THEN
              Imin=1
              Imax=Lm(ng)
              Jmin=1
              Jmax=Mm(ng)
            ELSE
              Imin=1
              Imax=Lm(ng)
              Jmin=2
              Jmax=Mm(ng)
            END IF
!
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '2D normalization factors at V-stress points'
              FLUSH (stdout)
            END IF
!
            DO j=JstrP,JendT
              DO i=IstrT,IendT
                Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_v(i,j)*             &
     &                                  GRID(ng)%on_v(i,j))
              END DO
            END DO
!
            DO jc=Jmin,Jmax
              DO ic=Imin,Imax
#   ifdef MASKING
                compute=0.0_r8
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  IF (GRID(ng)%vmask(ic,jc).gt.0) compute=1.0_r8
                END IF
#    ifdef DISTRIBUTE
                CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                compute=1.0_r8
#   endif
                IF (compute.gt.0.0_r8) THEN
                  DO j=LBj,UBj
                    DO i=LBi,UBi
                      A2d(i,j)=0.0_r8
                    END DO
                  END DO
                  IF (((Jstr.le.jc).and.(jc.le.Jend)).and.              &
     &                ((Istr.le.ic).and.(ic.le.Iend))) THEN
                    A2d(ic,jc)=1.0_r8
                  END IF
!
!  Implicit horizontal convolution, CG/CI solver.
!
                  CALL self%ad_CI_2d (ng, tile, iADM, isVstr, v2dvar,   &
     &                                ns, NiterCI(ns,ng), ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                A2d)
!
                  DO j=JstrP,JendT
                    DO i=IstrT,IendT
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
!
!  Set normalization factors.
!
                  Gdotp=dot_prod2d (ng, tile, iADM, v2dvar,             &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A2d, A2d)
                  cff=1.0_r8/SQRT(Gdotp)
                ELSE
                  cff=0.0_r8
                END IF
                IF (((Jstr.le.jc).and.(jc.le.Jend)).and.                &
     &              ((Istr.le.ic).and.(ic.le.Iend))) THEN
                  HnormSVS(ic,jc)=HnormSVS(ic,jc)+                      &
     &                            self%Bwgt(isVstr,ns)*cff
                END IF
              END DO
            END DO
          END DO MS_SVS_LOOP
!
!  Exchange boundary data.
!
          CALL dabc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSVS)
#  ifdef DISTRIBUTE
!
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSVS)
#  endif
!
!  Write out into output NetCDF file.
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idVsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idVsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              GRID(ng)%vmask,                     &
#  endif
     &                              HnormSVS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idVsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idVsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idVsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             GRID(ng)%vmask,                      &
#   endif
     &                             HnormSVS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF MS_COMPUTE_SVS
# endif /* ADJUST_WSTRESS */

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  2D normalization at surface tracer fluxes points.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          Lsame=.FALSE.
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              ifield=isTsur(itrc)
              IF (Cnorm(rec,ifield)) Lsame=.TRUE.
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
                  '2D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF
!
!  Check if the decorrelation scales for all the surface tracer fluxes
!  are different. If not, just compute the normalization factors for the
!  first tracer and assign the same value to the rest.  Recall that this
!  computation is very expensive.
!
        Ldiffer=.FALSE.
        DO ns=1,Nscale(ng)
          Imin=1
          Imax=Lm(ng)
          Jmin=1
          Jmax=Mm(ng)

# ifdef NONUNIFORM_SCALES
          Ldiffer=.TRUE.
# else
          DO itrc=2,NT(ng)
            IF ((HdecayX(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayX(rec,isTsur(itrc-1),ns,ng)).or.                 &
     &          (HdecayY(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayY(rec,isTsur(itrc-1),ns,ng))) THEN
              Ldiffer=.TRUE.
            END IF
          END DO
# endif
          IF (.not.Ldiffer) THEN
            Lsame=.TRUE.
            UBt=1
          ELSE
            Lsame=.FALSE.
            UBt=NT(ng)
          END IF
        END DO
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_r(i,j)*                 &
     &                              GRID(ng)%on_r(i,j))
          END DO
        END DO
!
        TRACER_FLUX_LOOP : DO itrc=1,UBt
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            MS_COMPUTE_STF : IF (Cnorm(rec,ifield)) THEN
              MS_STF_LOOP : DO ns=1,Nscale(ng)
                DO jc=Jmin,Jmax
                  DO ic=Imin,Imax
#   ifdef MASKING
                    compute=0.0_r8
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (GRID(ng)%rmask(ic,jc).gt.0) compute=1.0_r8
                    END IF
#    ifdef DISTRIBUTE
                    CALL mp_reduce (ng, iTLM, 1, compute, 'SUM')
#    endif
#   else
                    compute=1.0_r8
#   endif
                    IF (compute.gt.0.0_r8) THEN
                      DO j=LBj,UBj
                        DO i=LBi,UBi
                          A2d(i,j)=0.0_r8
                        END DO
                      END DO
                      IF (((Jstr.le.jc).and.(jc.le.Jend)).and.          &
     &                    ((Istr.le.ic).and.(ic.le.Iend))) THEN
                        A2d(ic,jc)=1.0_r8
                      END IF
!
!  Implicit horizontal convolution, CG/CI solver.
!
                      CALL self%ad_CI_2d (ng, tile, iADM, ifield,       &
     &                                    r2dvar,                       &
     &                                    ns, NiterCI(ns,ng), ifac,     &
     &                                    Lweak,                        &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    A2d)
!
                      DO j=JstrT,JendT
                        DO i=IstrT,IendT
                          A2d(i,j)=A2d(i,j)*Hscale(i,j)
                        END DO
                      END DO
!
!  Set normalization factors.
!
                      Gdotp=dot_prod2d (ng, tile, iADM, r2dvar,         &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  A2d, A2d)
                      cff=1.0_r8/SQRT(Gdotp)
                    ELSE
                      cff=0.0_r8
                    END IF
!
                    IF (((Jstr.le.jc).and.(jc.le.Jend)).and.            &
     &                  ((Istr.le.ic).and.(ic.le.Iend))) THEN
                      IF (Lsame) THEN
                        DO ntrc=1,NT(ng)
                          nfield=isTsur(ntrc)
                          IF (Lstflux(ntrc,ng)) THEN
                            HnormSTF(ic,jc,ntrc)=HnormSTF(ic,jc,ntrc)+  &
     &                                           self%Bwgt(nfield,ns)*  &
     &                                           cff
                          END IF
                        END DO
                      ELSE
                        HnormSTF(ic,jc,itrc)=HnormSTF(ic,jc,itrc)+      &
     &                                       self%Bwgt(ifield,ns)*cff
                      END IF
                    END IF
                  END DO
                END DO
              END DO MS_STF_LOOP
            END IF MS_COMPUTE_STF
          END IF
        END DO TRACER_FLUX_LOOP
!
!  Exchange boundary data.
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            IF (Cnorm(rec,ifield)) THEN
              CALL dabc_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            HnormSTF(:,:,itrc))
#  ifdef DISTRIBUTE
!
              CALL mp_exchange2d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            HnormSTF(:,:,itrc))
#  endif
!
!  Write out into output NetCDF file.
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  idTsur(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTsur(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                                  GRID(ng)%rmask,                 &
#  endif
     &                                  HnormSTF(:,:,itrc))

#  if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioVar(idTsur(itrc))%dkind.eq.      &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r2dvar(ng)
                  END IF
                  CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 idTsur(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioVar(idTsur(itrc)), &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                 ioDesc,                          &
#   ifdef MASKING
     &                                 GRID(ng)%rmask,                  &
#   endif
     &                                 HnormSTF(:,:,itrc))
#  endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
# endif /* ADJUST_STFLUX */
      END IF FRC_LWRITE_NRM
#endif
!
      IF (Master) THEN
        WRITE (stdout,30)
      END IF
!
 10   FORMAT (/,' Error Covariance Normalization Factors: ',            &
     &        'Exact Method',/)
 20   FORMAT (4x,'Computing',1x,a,1x,a)
 30   FORMAT (/)
!
      RETURN
      END SUBROUTINE normalization_tile

!
!***********************************************************************
      SUBROUTINE randomization_tile (self, ng, tile,                    &
     &                               LBi, UBi, LBj, UBj,                &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               nstp, nnew, ifac,                  &
     &                               Kh,                                &
#ifdef SOLVE3D
     &                               Kv,                                &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                               VnormRobc, VnormUobc, VnormVobc,   &
# endif
     &                               HnormRobc, HnormUobc, HnormVobc,   &
#endif
#ifdef ADJUST_WSTRESS
     &                               HnormSUS, HnormSVS,                &
#endif
#if defined ADJUST_STFLUX && defined SOLVE3D
     &                               HnormSTF,                          &
#endif
#ifdef SOLVE3D
     &                               VnormR, VnormU, VnormV,            &
#endif
     &                               HnormR, HnormU, HnormV)
!***********************************************************************
!
!  Imported variable declarations.
!
      CLASS (multiscale), intent(inout) :: self
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew, ifac
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Kh(LBi:,LBj:)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:,LBj:,0:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:,:,:,:)
      real(r8), intent(out) :: VnormUobc(LBij:,:,:)
      real(r8), intent(out) :: VnormVobc(LBij:,:,:)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:,:)
      real(r8), intent(out) :: HnormUobc(LBij:,:)
      real(r8), intent(out) :: HnormVobc(LBij:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:,LBj:)
      real(r8), intent(out) :: HnormSVS(LBi:,LBj:)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:,LBj:,:)
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: VnormU(LBi:,LBj:,:,:)
      real(r8), intent(out) :: VnormV(LBi:,LBj:,:,:)
# endif
      real(r8), intent(out) :: HnormR(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormU(LBi:,LBj:,:)
      real(r8), intent(out) :: HnormV(LBi:,LBj:,:)

#else

      real(r8), intent(in) :: Kh(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Kv(LBi:UBi,LBj:UBj,0:N(ng))
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(out) :: VnormRobc(LBij:UBij,N(ng),4,NT(ng))
      real(r8), intent(out) :: VnormUobc(LBij:UBij,N(ng),4)
      real(r8), intent(out) :: VnormVobc(LBij:UBij,N(ng),4)
#  endif
      real(r8), intent(out) :: HnormRobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormUobc(LBij:UBij,4)
      real(r8), intent(out) :: HnormVobc(LBij:UBij,4)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(out) :: HnormSUS(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: HnormSVS(LBi:UBi,LBj:UBj)
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(out) :: HnormSTF(LBi:UBi,LBj:UBj,NT(ng))
# endif
# ifdef SOLVE3D
      real(r8), intent(out) :: VnormR(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(out) :: VnormU(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(out) :: VnormV(LBi:UBi,LBj:UBj,N(ng),NSA)
# endif
      real(r8), intent(out) :: HnormR(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormU(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(out) :: HnormV(LBi:UBi,LBj:UBj,NSA)
# ifdef SOLVE3D
      real(r8), intent(out) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(out) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
# endif
#endif
!
!  Local variable declarations.
!
#ifdef SOLVE3D
      logical :: Ldiffer, Lsame
#endif
#ifdef ADJUST_BOUNDARY
      logical :: Lconvolve(4)
#endif
      logical :: Lweak
!
      integer :: i, ifile, is, iter, j, rec
      integer :: ifield, nfield
      integer :: ns
#ifdef SOLVE3D
      integer :: UBt, itrc, k
#endif
#ifdef ADJUST_BOUNDARY
      integer :: IDmeta, IJlen, IJKlen, ib, ibry, ic
#endif
      integer :: start(4), total(4)
!
      real(dp) :: my_time
      real(r8) :: Aavg, Amax, Amin, Asqr, FacAvg, FacSqr
      real(r8) :: cff

#ifdef ADJUST_BOUNDARY
      real(r8) :: Bavg, Bmin, Bmax, Bsqr

      real(r8), parameter :: Aspv = 0.0_r8
#endif
!
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2d
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2davg
      real(r8), dimension(LBi:UBi,LBj:UBj) :: A2dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj) :: Hscale
#ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij) :: B2d
      real(r8), dimension(LBij:UBij) :: B2davg
      real(r8), dimension(LBij:UBij) :: B2dsqr
      real(r8), dimension(LBij:UBij) :: HscaleB
#endif
#ifdef SOLVE3D
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3d
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3davg
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: A3dsqr
      real(r8), dimension(LBi:UBi,LBj:UBj,1:N(ng)) :: Vscale
# ifdef ADJUST_BOUNDARY
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3d
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3davg
      real(r8), dimension(LBij:UBij,1:N(ng)) :: B3dsqr
      real(r8), dimension(LBij:UBij,1:N(ng)) :: VscaleB
#  ifdef DISTRIBUTE
      real(r8), dimension((UBij-LBij+1)*N(ng)) :: Bwrk
#  endif
# endif
#endif
!
      character (len=40 ) :: Text
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", randomization_tile"

#if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
#endif

#include "set_bounds.h"
!
      SourceFile=MyFile

      my_time=tdays(ng)*day2sec

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Compute time invariant depths (use zero free-surface).
!-----------------------------------------------------------------------
!
      DO i=LBi,UBi
        DO j=LBj,UBj
          A2d(i,j)=0.0_r8                   ! free surface
        END DO
      END DO

      CALL set_depth_tile (ng, tile, iNLM,                              &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     GRID(ng) % h,                                &
# ifdef ICESHELF
     &                     GRID(ng) % zice,                             &
# endif
# if defined SEDIMENT && defined SED_MORPH
     &                     SEDBED(ng) % bed_thick,                      &
# endif
     &                     A2d,                                         &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
#endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  The prior (initial conditions) and model normalization factors are
!  computed for modeling the spreading of the multiscale background-
!  error covariance matrix (B) using the computationally efficient
!  RANDOMIZATION approach described by Fisher and Courtier (1995).
!
!  These factors ensure that the diagonal elements of B are unity. In
!  applications with land/sea masking, boundary conditions can
!  substantially modify covariance structures near boundaries.
!
!  The factors are initialized with random numbers ("white noise")
!  sampled from a uniform distribution with zero mean and unit variance.
!  These values are scaled by the inverse square root of the area (2D)
!  or volume (3D), then convolved using implicit horizontal and vertical
!  diffusion operators. This process is repeated for a specified number
!  of ensemble members, denoted as Nrandom.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      IF (Master) WRITE (stdout,10)

      FILE_LOOP : DO ifile=1,NSA

        LWRITE_NRM : IF (LwrtNRM(ifile,ng)) THEN
          Lweak=.FALSE.
          IF (ifile.eq.1) THEN
            Text='initial conditions'
          ELSE IF (ifile.eq.2) THEN
            Text='model'
            Lweak=.TRUE.
          END IF
!
!  Set randomization summation factors.
!
          FacAvg=1.0_r8/REAL(Nrandom,r8)
          FacSqr=SQRT(REAL(Nrandom,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
          ncname=NRM(ifile,ng)%name
          NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
          NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,idtime), my_time,           &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))
#if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,idtime), my_time,       &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
#endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  2D normalization at RHO-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_R2D : IF (Cnorm(ifile,isFsur)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '2D normalization factors at RHO-points'
              FLUSH (stdout)
            END IF
!
            MS_R2D_LOOP : DO ns=1,Nscale(ng)
              DO j=JstrT,JendT
                DO i=IstrT,IendT
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                  Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_r(i,j)*           &
     &                                    GRID(ng)%on_r(i,j))
                END DO
              END DO
!
              RANDOM_R2D : DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),      &
     &                              IstrR, IendR, JstrR, JendR,         &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Amin, Amax, A2d)
!
                DO j=JstrT,JendT
                  DO i=IstrT,IendT
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
!
!  Apply lateral boundary conditions.
!
                CALL dabc_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_2d (ng, tile, iTLM, isFsur, r2dvar,     &
     &                              ns, NiterCI(ns,ng), ifac, Lweak,    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              A2d)
!
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO RANDOM_R2D
!
!  Set normalization factors.
!
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  Aavg=FacAvg*A2davg(i,j)
                  Asqr=FacAvg*A2dsqr(i,j)
# ifdef MASKING
                  IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                    HnormR(i,j,ifile)=HnormR(i,j,ifile)+                &
     &                                self%Bwgt(isFsur,ns)/SQRT(Asqr)
                  ELSE
                    HnormR(i,j,ifile)=0.0_r8
                  END IF
# else
                  HnormR(i,j,ifile)=HnormR(i,j,ifile)+                  &
     &                              self%Bwgt(isFsur,ns)/SQRT(Asqr)
# endif
                END DO
              END DO
            END DO MS_R2D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormR(:,:,ifile))
#ifdef DISTRIBUTE
!
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormR(:,:,ifile))
#endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idFsur,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idFsur),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                GRID(ng)%rmask,                   &
#endif
     &                                HnormR(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idFsur)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idFsur,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idFsur),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               GRID(ng)%rmask,                    &
# endif
     &                               HnormR(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_R2D
!
!-----------------------------------------------------------------------
!  2D normalization at U-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_U2D : IF (Cnorm(ifile,isUbar)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '2D normalization factors at   U-points'
              FLUSH (stdout)
            END IF
!
            MS_U2D_LOOP : DO ns=1,Nscale(ng)
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                  Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_u(i,j)*           &
     &                                    GRID(ng)%on_u(i,j))
                END DO
              END DO
!
              RANDOM_U2D : DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),      &
     &                              Istr, IendR, JstrR, JendR,          &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Amin, Amax, A2d)
!
                DO j=JstrT,JendT
                  DO i=IstrP,IendT
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
!
!  Apply lateral boundary conditions.
!
                CALL dabc_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_2d (ng, tile, iTLM, isUbar, u2dvar,     &
     &                              ns, NiterCI(ns,ng), ifac, Lweak,    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              A2d)
!
                DO j=Jstr,Jend
                  DO i=IstrU,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO RANDOM_U2D
!
!  Set normalization factors.
!
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  Aavg=FacAvg*A2davg(i,j)
                  Asqr=FacAvg*A2dsqr(i,j)
# ifdef MASKING
                  IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                    HnormU(i,j,ifile)=HnormU(i,j,ifile)+                &
     &                                self%Bwgt(isUbar,ns)/SQRT(Asqr)
                  ELSE
                    HnormU(i,j,ifile)=0.0_r8
                  END IF
# else
                  HnormU(i,j,ifile)=HnormU(i,j,ifile)+                  &
     &                              self%Bwgt(isUbar,ns)/SQRT(Asqr)
# endif
                END DO
              END DO
            END DO MS_U2D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormU(:,:,ifile))
#ifdef DISTRIBUTE
!
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormU(:,:,ifile))
#endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idUbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idUbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                GRID(ng)%umask,                   &
#endif
     &                                HnormU(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idUbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idUbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               GRID(ng)%umask,                    &
# endif
     &                               HnormU(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_U2D
!
!-----------------------------------------------------------------------
!  2D normalization at V-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_V2D : IF (Cnorm(ifile,isVbar)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '2D normalization factors at   V-points'
              FLUSH (stdout)
            END IF
!
            MS_V2D_LOOP : DO ns=1,Nscale(ng)
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  A2davg(i,j)=0.0_r8
                  A2dsqr(i,j)=0.0_r8
                  Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_v(i,j)*           &
     &                                    GRID(ng)%on_v(i,j))
                END DO
              END DO
!
              RANDOM_V2D : DO iter=1,Nrandom
                CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),      &
     &                              IstrR, IendR, Jstr, JendR,          &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              Amin, Amax, A2d)
!
                DO j=JstrP,JendT
                  DO i=IstrT,IendT
                    A2d(i,j)=A2d(i,j)*Hscale(i,j)
                  END DO
                END DO
!
!  Apply lateral boundary conditions.
!
                CALL dabc_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A2d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_2d (ng, tile, iTLM, isVbar, v2dvar,     &
     &                              ns, NiterCI(ns,ng), ifac, Lweak,    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              A2d)
!
                DO j=JstrV,Jend
                  DO i=Istr,Iend
                    A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                    A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                  END DO
                END DO
              END DO RANDOM_V2D
!
!  Set normalization factors.
!
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  Aavg=FacAvg*A2davg(i,j)
                  Asqr=FacAvg*A2dsqr(i,j)
# ifdef MASKING
                  IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                    HnormV(i,j,ifile)=HnormV(i,j,ifile)+                &
     &                                self%Bwgt(isVbar,ns)/SQRT(Asqr)
                  ELSE
                    HnormV(i,j,ifile)=0.0_r8
                  END IF
# else
                  HnormV(i,j,ifile)=HnormV(i,j,ifile)+                  &
     &                              self%Bwgt(isVbar,ns)/SQRT(Asqr)
# endif
                END DO
              END DO
            END DO MS_V2D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          HnormV(:,:,ifile))
#ifdef DISTRIBUTE
!
            CALL mp_exchange2d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          HnormV(:,:,ifile))
#endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, idVbar,       &
     &                                NRM(ifile,ng)%ncid,               &
     &                                NRM(ifile,ng)%Vid(idVbar),        &
     &                                NRM(ifile,ng)%Rindex,             &
#ifdef MASKING
     &                                GRID(ng)%vmask,                   &
#endif
     &                                HnormV(:,:,ifile))

#if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVbar)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
                CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, idVbar,        &
     &                               NRM(ifile,ng)%pioFile,             &
     &                               NRM(ifile,ng)%pioVar(idVbar),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
# ifdef MASKING
     &                               GRID(ng)%vmask,                    &
# endif
     &                               HnormV(:,:,ifile))
#endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_V2D

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D normalization U-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_U3D : IF (Cnorm(ifile,isUvel)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '3D normalization factors at   U-points'
              FLUSH (stdout)
            END IF
!
            MS_U3D_LOOP : DO ns=1,Nscale(ng)
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  cff=GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j)*0.5_r8
                  DO k=1,N(ng)
                    A3davg(i,j,k)=0.0_r8
                    A3dsqr(i,j,k)=0.0_r8
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+       &
     &                                      GRID(ng)%Hz(i  ,j,k)))
                  END DO
                END DO
              END DO
!
              RANDOM_U3D : DO iter=1,Nrandom
                CALL white_noise3d (ng, iTLM, u3dvar, Rscheme(ng),      &
     &                              Istr, IendR, JstrR, JendR,          &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              Amin, Amax, A3d)
!
                DO k=1,N(ng)
                  DO j=JstrT,JendT
                    DO i=IstrP,IendT
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
!
!  Apply lateeral boundary conditions.
!
                CALL dabc_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_3d (ng, tile, iTLM, isUvel, u3dvar,     &
     &                              ns, NiterCI(ns,ng), ifac, Lweak,    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              A3d)
!
!  Implicit vertical convolution.
!
                CALL self%tl_Vdiff_u3d (ng, tile, iTLM, isUvel,         &
     &                                  NVsteps(ifile,isUvel)/ifac,     &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  DTsizeV(ifile,isUvel),          &
     &                                  Kv, A3d)
!
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    DO i=IstrU,Iend
                      A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                      A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                    END DO
                  END DO
                END DO
              END DO RANDOM_U3D
!
!  Set normalization factors.
!
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=IstrU,Iend
                    Aavg=FacAvg*A3davg(i,j,k)
                    Asqr=FacAvg*A3dsqr(i,j,k)
#  ifdef MASKING
                    IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                      VnormU(i,j,k,ifile)=VnormU(i,j,k,ifile)+          &
     &                                    self%Bwgt(isUvel,ns)/         &
     &                                    SQRT(Asqr)
                    ELSE
                      VnormU(i,j,k,ifile)=0.0_r8
                    END IF
#  else
                    VnormU(i,j,k,ifile)=VnormU(i,j,k,ifile)+            &
     &                                  self%Bwgt(isUvel,ns)/SQRT(Asqr)
#  endif
                  END DO
                END DO
              END DO
            END DO MS_U3D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormU(:,:,:,ifile))
# ifdef DISTRIBUTE
!
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormU(:,:,:,ifile))
# endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idUvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idUvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                GRID(ng)%umask,                   &
# endif
     &                                VnormU(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idUvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_u3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idUvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idUvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               GRID(ng)%umask,                    &
#  endif
     &                               VnormU(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_U3D
!
!-----------------------------------------------------------------------
!  3D normalization at V-points.
!-----------------------------------------------------------------------
!
          MS_COMPUTE_V3D : IF (Cnorm(ifile,isVvel)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '3D normalization factors at   V-points'
              FLUSH (stdout)
            END IF
!
            MS_V3D_LOOP : DO ns=1,Nscale(ng)
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  cff=GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j)*0.5_r8
                  DO k=1,N(ng)
                    A3davg(i,j,k)=0.0_r8
                    A3dsqr(i,j,k)=0.0_r8
                    Vscale(i,j,k)=1.0_r8/                               &
     &                            SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+       &
     &                                      GRID(ng)%Hz(i,j  ,k)))
                  END DO
                END DO
              END DO
!
              RANDOM_V3D : DO iter=1,Nrandom
                CALL white_noise3d (ng, iTLM, v3dvar, Rscheme(ng),      &
     &                              IstrR, IendR, Jstr, JendR,          &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              Amin, Amax, A3d)
!
                DO k=1,N(ng)
                  DO j=JstrP,JendT
                    DO i=IstrT,IendT
                      A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                    END DO
                  END DO
                END DO
!
!  Apply lateeral boundary conditions.
!
                CALL dabc_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_3d (ng, tile, iTLM, isVvel, v3dvar,     &
     &                              ns, NiterCI(ns,ng), ifac, Lweak,    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              A3d)
!
!  Implicit vertical convolution.
!
                CALL self%tl_Vdiff_v3d (ng, tile, iTLM, isVvel,         &
     &                                  NVsteps(ifile,isUvel)/ifac,     &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  IminS, ImaxS, JminS, JmaxS,     &
     &                                  DTsizeV(ifile,isUvel),          &
     &                                  Kv, A3d)
!
                DO k=1,N(ng)
                  DO j=JstrV,Jend
                    DO i=Istr,Iend
                      A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                      A3dsqr(i,j,k)=A3dsqr(i,j,k)+A3d(i,j,k)*A3d(i,j,k)
                    END DO
                  END DO
                END DO
              END DO RANDOM_V3D
!
!  Set normalization factors.
!
              DO k=1,N(ng)
                DO j=JstrV,Jend
                  DO i=Istr,Iend
                    Aavg=FacAvg*A3davg(i,j,k)
                    Asqr=FacAvg*A3dsqr(i,j,k)
#  ifdef MASKING
                    IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                      VnormV(i,j,k,ifile)=VnormV(i,j,k,ifile)+          &
     &                                    self%Bwgt(isVvel,ns)/         &
     &                                    SQRT(Asqr)
                    ELSE
                      VnormV(i,j,k,ifile)=0.0_r8
                    END IF
#  else
                    VnormV(i,j,k,ifile)=VnormV(i,j,k,ifile)+            &
     &                                  self%Bwgt(isVvel,ns)/SQRT(Asqr)
#  endif
                  END DO
                END DO
              END DO
            END DO MS_V3D_LOOP
!
!  Exchange boundary data.
!
            CALL dabc_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          VnormV(:,:,:,ifile))
# ifdef DISTRIBUTE
!
            CALL mp_exchange3d (ng, tile, iTLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          VnormV(:,:,:,ifile))
# endif
!
!  Write out into output NetCDF file.
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,           &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                idVvel, NRM(ifile,ng)%ncid,       &
     &                                NRM(ifile,ng)%Vid(idVvel),        &
     &                                NRM(ifile,ng)%Rindex,             &
# ifdef MASKING
     &                                GRID(ng)%vmask,                   &
# endif
     &                                VnormV(:,:,:,ifile))

# if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                IF (NRM(ifile,ng)%pioVar(idVvel)%dkind.eq.              &
     &              PIO_double) THEN
                  ioDesc => ioDesc_dp_v3dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v3dvar(ng)
                END IF
                CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,            &
     &                               LBi, UBi, LBj, UBj, 1, N(ng),      &
     &                               idVvel, NRM(ifile,ng)%pioFile,     &
     &                               NRM(ifile,ng)%pioVar(idVvel),      &
     &                               NRM(ifile,ng)%Rindex,              &
     &                               ioDesc,                            &
#  ifdef MASKING
     &                               GRID(ng)%vmask,                    &
#  endif
     &                               VnormV(:,:,:,ifile))
# endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF MS_COMPUTE_V3D
!
!-----------------------------------------------------------------------
!  3D normalization at RHO-points.
!-----------------------------------------------------------------------
!
          IF (Master) THEN
            Lsame=.FALSE.
            DO itrc=1,NT(ng)
              ifield=isTvar(itrc)
              IF (Cnorm(ifile,ifield)) Lsame=.TRUE.
            END DO
            IF (Lsame) THEN
              WRITE (stdout,20) TRIM(Text),                             &
     &              '3D normalization factors at RHO-points'
              FLUSH (stdout)
            END IF
          END IF
!
!  Check if the decorrelation scales for all the tracers are different.
!  If not, just compute the normalization factors for the first tracer
!  and assign the same value to the rest.  Recall that this computation
!  is very expensive.
!
          Ldiffer=.FALSE.
          DO ns=1,Nscale(ng)

# ifdef NONUNIFORM_SCALES
            Ldiffer=.TRUE.
# else
            DO itrc=2,NT(ng)
              IF ((HdecayX(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayX(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (HdecayY(ifile,isTvar(itrc  ),ns,ng).ne.              &
     &             HdecayY(ifile,isTvar(itrc-1),ns,ng)).or.             &
     &            (Vdecay(ifile,isTvar(itrc  ),ng).ne.                  &
     &             Vdecay(ifile,isTvar(itrc-1),ng))) THEN
                Ldiffer=.TRUE.
              END IF
            END DO
# endif
            IF (.not.Ldiffer) THEN
              Lsame=.TRUE.
              UBt=1
            ELSE
              Lsame=.FALSE.
              UBt=NT(ng)
            END IF
          END DO
!
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              cff=GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j)
              DO k=1,N(ng)
                Vscale(i,j,k)=1.0_r8/SQRT(cff*GRID(ng)%Hz(i,j,k))
              END DO
            END DO
          END DO
!
          TRACER_LOOP : DO itrc=1,UBt
            ifield=isTvar(itrc)
            MS_COMPUTE_R3D : IF (Cnorm(ifile,ifield)) THEN
              MS_R3D_LOOP : DO ns=1,Nscale(ng)
!
                DO k=1,N(ng)
                  DO j=JstrT,JendT
                    DO i=IstrT,IendT
                      A3davg(i,j,k)=0.0_r8
                      A3dsqr(i,j,k)=0.0_r8
                    END DO
                  END DO
                END DO
!
                RANDOM_R3D : DO iter=1,Nrandom
                  CALL white_noise3d (ng, iTLM, r3dvar, Rscheme(ng),    &
     &                                IstrR, IendR, JstrR, JendR,       &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                Amin, Amax, A3d)
!
                  DO k=1,N(ng)
                    DO j=JstrT,JendT
                      DO i=IstrT,IendT
                        A3d(i,j,k)=A3d(i,j,k)*Vscale(i,j,k)
                      END DO
                    END DO
                  END DO
!
!  Apply lateral boundary conditions.
!
                  CALL dabc_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                A3d)
!
!  Implicit horizontal convolution, CG/CI solver.
!
                  CALL self%tl_CI_3d (ng, tile, iTLM, ifield, r3dvar,   &
     &                                ns, NiterCI(ns,ng), ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                A3d)
!
!  Implicit vertical convolution.
!
                  CALL self%tl_Vdiff_r3d (ng, tile, iTLM, ifield,       &
     &                                    NVsteps(ifile,ifield)/ifac,   &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    IminS, ImaxS, JminS, JmaxS,   &
     &                                    DTsizeV(ifile,ifield),        &
     &                                    Kv, A3d)
!
                  DO k=1,N(ng)
                    DO j=Jstr,Jend
                      DO i=Istr,Iend
                        A3davg(i,j,k)=A3davg(i,j,k)+A3d(i,j,k)
                        A3dsqr(i,j,k)=A3dsqr(i,j,k)+                    &
     &                                A3d(i,j,k)*A3d(i,j,k)
                      END DO
                    END DO
                  END DO
                END DO RANDOM_R3D
!
!  Set normalization factors.
!
                DO k=1,N(ng)
                  DO j=Jstr,Jend
                    DO i=Istr,Iend
                      Aavg=FacAvg*A3davg(i,j,k)
                      Asqr=FacAvg*A3dsqr(i,j,k)
#  ifdef MASKING
                      IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                        VnormR(i,j,k,ifile,itrc)=                       &
     &                                    VnormR(i,j,k,ifile,itrc)+     &
     &                                    self%Bwgt(ifield,ns)/         &
     &                                    SQRT(Asqr)
                      ELSE
                        VnormR(i,j,k,ifile,itrc)=0.0_r8
                      END IF
#  else
                      VnormR(i,j,k,ifile,itrc)=                         &
     &                                  VnormR(i,j,k,ifile,itrc)+       &
     &                                  self%Bwgt(ifield,ns)/           &
     &                                  SQRT(Asqr)
#  endif
                    END DO
                  END DO
                END DO
              END DO MS_R3D_LOOP
            END IF MS_COMPUTE_R3D
          END DO TRACER_LOOP
!
!  If same correlation parameters, replicate normalization factors
!  for the remaining tracer fields.
!
          IF (Lsame) THEN
            DO itrc=2,NT(ng)
              DO k=1,N(ng)
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    VnormR(i,j,k,ifile,itrc)=VnormR(i,j,k,ifile,1)
                  END DO
                END DO
              END DO
            END DO
          END IF
!
!  Exchange boundary data.
!
          DO itrc=1,NT(ng)
            ifield=isTvar(itrc)
            IF (Cnorm(ifile,ifield)) THEN
              CALL dabc_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            VnormR(:,:,:,ifile,itrc))
# ifdef DISTRIBUTE
!
              CALL mp_exchange3d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            VnormR(:,:,:,ifile,itrc))
# endif
!
!  Write out into output NetCDF file.
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm3d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                                  idTvar(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTvar(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
# ifdef MASKING
     &                                  GRID(ng)%rmask,                 &
# endif
     &                                  VnormR(:,:,:,ifile,itrc))

# if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioTrc(itrc)%dkind.eq.              &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r3dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r3dvar(ng)
                  END IF
                  CALL wrt_norm3d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj, 1, N(ng),    &
     &                                 idTvar(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioTrc(itrc),         &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                 ioDesc,                          &
#  ifdef MASKING
     &                                 GRID(ng)%rmask,                  &
#  endif
     &                                 VnormR(:,:,:,ifile,itrc))
# endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END DO
#endif
        END IF LWRITE_NRM
      END DO FILE_LOOP

#ifdef ADJUST_BOUNDARY
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute open boundaries multiscale background-error covariance, B,
!  normalization factors using the randomization approach of Fisher and
!  Courtier (1995).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      ifile=3
      BRY_LWRITE_NRM : IF (LwrtNRM(ifile,ng)) THEN
        Text='boundary conditions'
        IJlen=UBij-LBij+1
# ifdef SOLVE3D
        IJKlen=IJlen*N(ng)
# endif
        Lconvolve(iwest )=DOMAIN(ng)%Western_Edge (tile)
        Lconvolve(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
        Lconvolve(isouth)=DOMAIN(ng)%Southern_Edge(tile)
        Lconvolve(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
!  Set randomization summation factors.
!
        FacAvg=1.0_r8/REAL(Nrandom,r8)
        FacSqr=SQRT(REAL(Nrandom,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  2D boundary normalization at RHO-points.
!-----------------------------------------------------------------------
!
        HnormRobc=Aspv

        IF (Master.and.ANY(CnormB(isFsur,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '2D normalization factors at RHO-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_R2D : DO ibry=1,4
          BRY_COMPUTE_R2D : IF (CnormB(isFsur,ibry)) THEN
            BRY_MS_LOOP_R2D : DO ns=1,Nscale(ng)
              HscaleB=0.0_r8
              B2davg=0.0_r8
              B2dsqr=0.0_r8
              B2d=0.0_r8
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,r2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      HscaleB(j)=1.0_r8/SQRT(GRID(ng)%on_r(i,j))
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,r2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      HscaleB(i)=1.0_r8/SQRT(GRID(ng)%om_r(i,j))
                    END DO
                  END IF
              END SELECT
!
              BRY_RANDOM_R2D : DO iter=1,Nrandom
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    CALL white_noise2d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      JstrR, JendR,               &
     &                                      LBij, UBij,                 &
     &                                      Bmin, Bmax, B2d)
                    DO j=JstrT,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  CASE (isouth, inorth)
                    CALL white_noise2d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      IstrR, IendR,               &
     &                                      LBij, UBij,                 &
     &                                      Bmin, Bmax, B2d)
                    DO i=IstrT,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                END SELECT
!
!  Implicit tangent linear convolution, CG/CI solver.
!
                CALL self%tl_CI_b1d (ng, tile, iTLM, isFsur, ibry,      &
     &                               r2dvar,                            &
     &                               ns, NiterCI(ns,ng), ifac,          &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               B2d)
!
                IF (Lconvolve(ibry)) THEN
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO j=Jstr,Jend
                        B2davg(j)=B2davg(j)+B2d(j)
                        B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                      END DO
                    CASE (isouth, inorth)
                      DO i=Istr,Iend
                        B2davg(i)=B2davg(i)+B2d(i)
                        B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                      END DO
                  END SELECT
                END IF
              END DO BRY_RANDOM_R2D
!
!  Set normalization factors.
!
              IF (Lconvolve(ibry)) THEN
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,r2dvar)
                    DO j=Jstr,Jend
                      Bavg=FacAvg*B2davg(j)
                      Bsqr=FacAvg*B2dsqr(j)
#  ifdef MASKING
                      IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                        HnormRobc(j,ibry)=HnormRobc(j,ibry)+            &
     &                                    self%Bwgt(isFsur,ns)/         &
     &                                    SQRT(Bsqr)
                      ELSE
                        HnormRobc(j,ibry)=0.0_r8
                      END IF
#  else
                      HnormRobc(j,ibry)=HnormRobc(j,ibry)+              &
     &                                  self%Bwgt(isFsur,ns)/SQRT(Bsqr)
#  endif
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,r2dvar)
                    DO i=Istr,Iend
                      Bavg=FacAvg*B2davg(i)
                      Bsqr=FacAvg*B2dsqr(i)
#  ifdef MASKING
                      IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                        HnormRobc(i,ibry)=HnormRobc(i,ibry)+            &
     &                                    self%Bwgt(isFsur,ns)/         &
     &                                    SQRT(Bsqr)
                      ELSE
                        HnormRobc(i,ibry)=0.0_r8
                      END IF
#  else
                      HnormRobc(i,ibry)=HnormRobc(i,ibry)+              &
     &                                  self%Bwgt(isFsur,ns)/SQRT(Bsqr)
#  endif
                    END DO
                END SELECT
              END IF
            END DO BRY_MS_LOOP_R2D
!
!  Exchange boundary data.
!
            CALL bc_r2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormRobc(:,ibry))
# ifdef DISTRIBUTE
!
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormRobc(LBij:,ibry))
# endif
          END IF BRY_COMPUTE_R2D
        END DO BRY_LOOP_R2D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isFsur,:))) THEN
          IDmeta=idSbry(isFsur)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              HnormRobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  HnormRobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)

# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D boundary normalization at U-points.
!-----------------------------------------------------------------------
!
        HnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '2D normalization factors at   U-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_U2D : DO ibry=1,4
          BRY_COMPUTE_U2D : IF (CnormB(isUbar,ibry)) THEN
            BRY_MS_LOOP_U2D : DO ns=1,Nscale(ng)
              HscaleB=0.0_r8
              B2davg=0.0_r8
              B2dsqr=0.0_r8
              B2d=0.0_r8
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,u2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      HscaleB(j)=1.0_r8/SQRT(GRID(ng)%on_u(i,j))
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,u2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrP,IendT
                      HscaleB(i)=1.0_r8/SQRT(GRID(ng)%om_u(i,j))
                    END DO
                  END IF
              END SELECT
!
              BRY_RANDOM_U2D : DO iter=1,Nrandom
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    CALL white_noise2d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      JstrR, JendR,               &
     &                                      LBij, UBij,                 &
     &                                      Bmin, Bmax, B2d)
                    DO j=JstrT,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  CASE (isouth, inorth)
                    CALL white_noise2d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      Istr, IendR,                &
     &                                      LBij, UBij,                 &
     &                                      Bmin, Bmax, B2d)
                    DO i=IstrP,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                END SELECT
!
!  Implicit tangent linear convolution, CG/CI solver.
!
                CALL self%tl_CI_b1d (ng, tile, iTLM, isUbar, ibry,      &
     &                               u2dvar,                            &
     &                               ns, NiterCI(ns,ng), ifac,          &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               B2d)
!
                IF (Lconvolve(ibry)) THEN
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO j=Jstr,Jend
                        B2davg(j)=B2davg(j)+B2d(j)
                        B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                      END DO
                    CASE (isouth, inorth)
                      DO i=IstrU,Iend
                        B2davg(i)=B2davg(i)+B2d(i)
                        B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                      END DO
                  END SELECT
                END IF
              END DO BRY_RANDOM_U2D
!
!  Set normalization factors.
!
              IF (Lconvolve(ibry)) THEN
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,u2dvar)
                    DO j=Jstr,Jend
                      Bavg=FacAvg*B2davg(j)
                      Bsqr=FacAvg*B2dsqr(j)
#  ifdef MASKING
                      IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                        HnormUobc(j,ibry)=HnormUobc(j,ibry)+            &
     &                                    self%Bwgt(isUbar,ns)/         &
     &                                    SQRT(Bsqr)
                      ELSE
                        HnormUobc(j,ibry)=0.0_r8
                      END IF
#  else
                      HnormUobc(j,ibry)=HnormUobc(j,ibry)+              &
     &                                  self%Bwgt(isUbar,ns)/SQRT(Bsqr)
# endif
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,u2dvar)
                    DO i=IstrU,Iend
                      Bavg=FacAvg*B2davg(i)
                      Bsqr=FacAvg*B2dsqr(i)
#  ifdef MASKING
                      IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                        HnormUobc(i,ibry)=HnormUobc(i,ibry)+            &
     &                                    self%Bwgt(isUbar,ns)/         &
     &                                    SQRT(Bsqr)
                      ELSE
                        HnormUobc(i,ibry)=0.0_r8
                      END IF
#  else
                      HnormUobc(i,ibry)=HnormUobc(i,ibry)+              &
     &                                  self%Bwgt(isUbar,ns)/SQRT(Bsqr)
# endif
                    END DO
                END SELECT
              END IF
            END DO BRY_MS_LOOP_U2D
!
!  Exchange boundary data.
!
            CALL bc_u2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormUobc(:,ibry))
# ifdef DISTRIBUTE
!
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormUobc(LBij:,ibry))
# endif
          END IF BRY_COMPUTE_U2D
        END DO BRY_LOOP_U2D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isUbar,:))) THEN
          IDmeta=idSbry(isUbar)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              HnormUobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  HnormUobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D boundary normalization at V-points.
!-----------------------------------------------------------------------
!
        HnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVbar,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '2D normalization factors at   V-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_V2D : DO ibry=1,4
          BRY_COMPUTE_V2D : IF (CnormB(isVbar,ibry)) THEN
            BRY_MS_LOOP_V2D : DO ns=1,Nscale(ng)
              HscaleB=0.0_r8
              B2davg=0.0_r8
              B2dsqr=0.0_r8
              B2d=0.0_r8
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,v2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrP,JendT
                      HscaleB(j)=1.0_r8/SQRT(GRID(ng)%on_v(i,j))
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,v2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      HscaleB(i)=1.0_r8/SQRT(GRID(ng)%om_v(i,j))
                    END DO
                  END IF
              END SELECT
!
              BRY_RANDOM_V2D : DO iter=1,Nrandom
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    CALL white_noise2d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      Jstr, JendR,                &
     &                                      LBij, UBij,                 &
     &                                      Bmin, Bmax, B2d)
                    DO j=JstrP,JendT
                      B2d(j)=B2d(j)*HscaleB(j)
                    END DO
                  CASE (isouth, inorth)
                    CALL white_noise2d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      IstrR, IendR,               &
     &                                      LBij, UBij,                 &
     &                                      Bmin, Bmax, B2d)
                    DO i=IstrT,IendT
                      B2d(i)=B2d(i)*HscaleB(i)
                    END DO
                END SELECT
!
!  Implicit tangent linear convolution, CG/CI solver.
!
                CALL self%tl_CI_b1d (ng, tile, iTLM, isVbar, ibry,      &
     &                               v2dvar,                            &
     &                               ns, NiterCI(ns,ng), ifac,          &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               B2d)
!
                IF (Lconvolve(ibry)) THEN
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO j=JstrV,Jend
                        B2davg(j)=B2davg(j)+B2d(j)
                        B2dsqr(j)=B2dsqr(j)+B2d(j)*B2d(j)
                      END DO
                    CASE (isouth, inorth)
                      DO i=Istr,Iend
                        B2davg(i)=B2davg(i)+B2d(i)
                        B2dsqr(i)=B2dsqr(i)+B2d(i)*B2d(i)
                      END DO
                  END SELECT
                END IF
              END DO BRY_RANDOM_V2D
!
!  Set normalization factors.
!
              IF (Lconvolve(ibry)) THEN
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,v2dvar)
                    DO j=JstrV,Jend
                      Bavg=FacAvg*B2davg(j)
                      Bsqr=FacAvg*B2dsqr(j)
#  ifdef MASKING
                      IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                        HnormVobc(j,ibry)=HnormVobc(j,ibry)+            &
     &                                    self%Bwgt(isVbar,ns)/         &
     &                                    SQRT(Bsqr)
                      ELSE
                        HnormVobc(j,ibry)=0.0_r8
                      END IF
#  else
                      HnormVobc(j,ibry)=HnormVobc(j,ibry)+              &
     &                                  self%Bwgt(isVbar,ns)/SQRT(Bsqr)
#  endif
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,v2dvar)
                    DO i=Istr,Iend
                      Bavg=FacAvg*B2davg(i)
                      Bsqr=FacAvg*B2dsqr(i)
#  ifdef MASKING
                      IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                        HnormVobc(i,ibry)=HnormVobc(i,ibry)+            &
     &                                    self%Bwgt(isVbar,ns)/         &
     &                                    SQRT(Bsqr)
                      ELSE
                        HnormVobc(i,ibry)=0.0_r8
                      END IF
#  else
                      HnormVobc(i,ibry)=HnormVobc(i,ibry)+              &
     &                                  self%Bwgt(isVbar,ns)/           &
     &                                  SQRT(Bsqr)
#  endif
                    END DO
                END SELECT
              END IF
            END DO BRY_MS_LOOP_V2D
!
!  Exchange boundary data.
!
            CALL bc_v2d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij,                           &
     &                            HnormVobc(:,ibry))
# ifdef DISTRIBUTE
!
            CALL mp_collect (ng, iTLM, IJlen, Aspv,                     &
     &                       HnormVobc(LBij:,ibry))
# endif
          END IF BRY_COMPUTE_V2D
        END DO BRY_LOOP_V2D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isVbar,:))) THEN
          IDmeta=idSbry(isVbar)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              HnormVobc(LBij:,:),                 &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  HnormVobc(LBij:,:),             &
     &                         start = (/1,1,NRM(ifile,ng)%Rindex/),    &
     &                         total = (/IJlen,4,1/),                   &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
# endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D boundary norm at U-points.
!-----------------------------------------------------------------------
!
        VnormUobc=Aspv

        IF (Master.and.ANY(CnormB(isUvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '3D normalization factors at   U-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_U3D : DO ibry=1,4
          BRY_COMPUTE_U3D : IF (CnormB(isUvel,ibry)) THEN
            BRY_MS_LOOP_U3D : DO ns=1,Nscale(ng)
              VscaleB=0.0_r8
              B3davg=0.0_r8
              B3dsqr=0.0_r8
              B3d=0.0_r8
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,u2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrT,JendT
                      cff=GRID(ng)%on_u(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(j,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+    &
     &                                         GRID(ng)%Hz(i  ,j,k)))
                      END DO
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,u2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrP,IendT
                      cff=GRID(ng)%om_u(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(i,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i-1,j,k)+    &
     &                                         GRID(ng)%Hz(i  ,j,k)))
                      END DO
                    END DO
                  END IF
              END SELECT
!
              BRY_RANDOM_U3D : DO iter=1,Nrandom
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    CALL white_noise3d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      JstrR, JendR,               &
     &                                      LBij, UBij, 1, N(ng),       &
     &                                      Bmin, Bmax, B3d)
                    DO k=1,N(ng)
                      DO j=JstrT,JendT
                        B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    CALL white_noise3d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      Istr, IendR,                &
     &                                      LBij, UBij, 1, N(ng),       &
     &                                      Bmin, Bmax, B3d)
                    DO k=1,N(ng)
                      DO i=IstrP,IendT
                        B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                      END DO
                    END DO
                END SELECT
!
!  Implicit tangent linear horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_b2d (ng, tile, iTLM, isUvel, ibry,      &
     &                               u3dvar,                            &
     &                               ns, NiterCI(ns,ng), ifac,          &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               B3d)
!
!  Implicit tangent linear vertical convolution.
!
                CALL self%tl_bry_Vdiff (ng, tile, iTLM, isUvel, ibry,   &
     &                                  u3dvar,                         &
     &                                  NVstepsB(ibry,isUvel)/ifac,     &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  LBij, UBij,                     &
     &                                  DTsizeVB(ibry,isUvel), Kv,      &
     &                                  B3d)
!
                IF (Lconvolve(ibry)) THEN
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO k=1,N(ng)
                        DO j=Jstr,Jend
                          B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                          B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                        END DO
                      END DO
                    CASE (isouth, inorth)
                      DO k=1,N(ng)
                        DO i=IstrU,Iend
                          B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                          B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                        END DO
                      END DO
                  END SELECT
                END IF
              END DO BRY_RANDOM_U3D
!
!  Set normalization factors.
!
              IF (Lconvolve(ibry)) THEN
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,u2dvar)
                    DO k=1,N(ng)
                      DO j=Jstr,Jend
                        Bavg=FacAvg*B3davg(j,k)
                        Bsqr=FacAvg*B3dsqr(j,k)
#  ifdef MASKING
                        IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                          VnormUobc(j,k,ibry)=VnormUobc(j,k,ibry)+      &
     &                                        self%Bwgt(isUvel,ns)/     &
     &                                        SQRT(Bsqr)
                        ELSE
                          VnormUobc(j,k,ibry)=0.0_r8
                        END IF
#  else
                        VnormUobc(j,k,ibry)=VnormUobc(j,k,ibry)+        &
     &                                      self%Bwgt(isUvel,ns)/       &
     &                                      SQRT(Bsqr)
#  endif
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,u2dvar)
                    DO k=1,N(ng)
                      DO i=IstrU,Iend
                        Bavg=FacAvg*B3davg(i,k)
                        Bsqr=FacAvg*B3dsqr(i,k)
#  ifdef MASKING
                        IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                          VnormUobc(i,k,ibry)=VnormUobc(i,k,ibry)+      &
     &                                        self%Bwgt(isUvel,ns)/     &
     &                                        SQRT(Bsqr)
                        ELSE
                          VnormUobc(i,k,ibry)=0.0_r8
                        END IF
#  else
                        VnormUobc(i,k,ibry)=VnormUobc(i,k,ibry)+        &
     &                                      self%Bwgt(isUvel,ns)/       &
     &                                      SQRT(Bsqr)
#  endif
                      END DO
                    END DO
                END SELECT
              END IF
            END DO BRY_MS_LOOP_U3D
!
!  Exchange boundary data.
!
            CALL bc_u3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormUobc(:,:,ibry))
#  ifdef DISTRIBUTE
!
            Bwrk=RESHAPE(VnormUobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormUobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF BRY_COMPUTE_U3D
        END DO BRY_LOOP_U3D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isUvel,:))) THEN
          IDmeta=idSbry(isUvel)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              VnormUobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  VnormUobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D boundary normalization at V-points.
!-----------------------------------------------------------------------
!
        VnormVobc=Aspv

        IF (Master.and.ANY(CnormB(isVvel,:))) THEN
          WRITE (stdout,20) TRIM(Text),                                 &
     &          '3D normalization factors at   V-points'
          FLUSH (stdout)
        END IF
!
        BRY_LOOP_V3D : DO ibry=1,4
          BRY_COMPUTE_V3D : IF (CnormB(isVvel,ibry)) THEN
            BRY_MS_LOOP_V3D : DO ns=1,Nscale(ng)
              VscaleB=0.0_r8
              B3davg=0.0_r8
              B3dsqr=0.0_r8
              B3d=0.0_r8
!
              SELECT CASE (ibry)
                CASE (iwest, ieast)
                  i=BOUNDS(ng)%edge(ibry,v2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO j=JstrP,JendT
                      cff=GRID(ng)%on_v(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(j,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+    &
     &                                         GRID(ng)%Hz(i,j  ,k)))
                      END DO
                    END DO
                  END IF
                CASE (isouth, inorth)
                  j=BOUNDS(ng)%edge(ibry,v2dvar)
                  IF (Lconvolve(ibry)) THEN
                    DO i=IstrT,IendT
                      cff=GRID(ng)%om_v(i,j)*0.5_r8
                      DO k=1,N(ng)
                        VscaleB(i,k)=1.0_r8/                            &
     &                               SQRT(cff*(GRID(ng)%Hz(i,j-1,k)+    &
     &                                         GRID(ng)%Hz(i,j,  k)))
                      END DO
                    END DO
                  END IF
              END SELECT
!
              BRY_RANDOM_V3D : DO iter=1,Nrandom
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    CALL white_noise3d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      Jstr, JendR,                &
     &                                      LBij, UBij, 1, N(ng),       &
     &                                      Bmin, Bmax, B3d)
                    DO k=1,N(ng)
                      DO j=JstrP,JendT
                        B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    CALL white_noise3d_bry (ng, tile, iTLM, ibry,       &
     &                                      Rscheme(ng),                &
     &                                      IstrR, IendR,               &
     &                                      LBij, UBij, 1, N(ng),       &
     &                                      Bmin, Bmax, B3d)
                    DO k=1,N(ng)
                      DO i=IstrT,IendT
                        B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                      END DO
                    END DO
                END SELECT
!
!  Implicit tangent linear horizontal convolution, CG/CI solver.
!
                CALL self%tl_CI_b2d (ng, tile, iTLM, isVvel, ibry,      &
     &                               v3dvar,                            &
     &                               ns, NiterCI(ns,ng), ifac,          &
     &                               LBij, UBij,                        &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               B3d)
!
!  Implicit tangent linnear vertical convolution.
!
                CALL self%tl_bry_Vdiff (ng, tile, iTLM, isVvel, ibry,   &
     &                                  v3dvar,                         &
     &                                  NVstepsB(ibry,isVvel)/ifac,     &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  LBij, UBij,                     &
     &                                  DTsizeVB(ibry,isVvel), Kv,      &
     &                                  B3d)
!
                IF (Lconvolve(ibry)) THEN
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      DO k=1,N(ng)
                        DO j=JstrV,Jend
                          B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                          B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                        END DO
                      END DO
                    CASE (isouth, inorth)
                      DO k=1,N(ng)
                        DO i=Istr,Iend
                          B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                          B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                        END DO
                      END DO
                  END SELECT
                END IF
              END DO BRY_RANDOM_V3D
!
!  Set normalization factors.
!
              IF (Lconvolve(ibry)) THEN
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,v2dvar)
                    DO k=1,N(ng)
                      DO j=JstrV,Jend
                        Bavg=FacAvg*B3davg(j,k)
                        Bsqr=FacAvg*B3dsqr(j,k)
#  ifdef MASKING
                        IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                          VnormVobc(j,k,ibry)=VnormVobc(j,k,ibry)+      &
     &                                        self%Bwgt(isVvel,ns)/     &
     &                                        SQRT(Bsqr)
                        ELSE
                          VnormVobc(j,k,ibry)=0.0_r8
                        END IF
#  else
                        VnormVobc(j,k,ibry)=VnormVobc(j,k,ibry)+        &
     &                                      self%Bwgt(isVvel,ns)/       &
     &                                      SQRT(Bsqr)
#  endif
                      END DO
                    END DO
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,v2dvar)
                    DO k=1,N(ng)
                      DO i=Istr,Iend
                        Bavg=FacAvg*B3davg(i,k)
                        Bsqr=FacAvg*B3dsqr(i,k)
#  ifdef MASKING
                        IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                          VnormVobc(i,k,ibry)=VnormVobc(i,k,ibry)+      &
     &                                        self%Bwgt(isVvel,ns)/     &
     &                                        SQRT(Bsqr)
                        ELSE
                          VnormVobc(i,k,ibry)=0.0_r8
                        END IF
#  else
                        VnormVobc(i,k,ibry)=VnormVobc(i,k,ibry)+        &
     &                                      self%Bwgt(isVvel,ns)/       &
     &                                      SQRT(Bsqr)
#  endif
                      END DO
                    END DO
                END SELECT
              END IF
            END DO BRY_MS_LOOP_V3D
!
!  Exchange boundary data.
!
            CALL bc_v3d_bry_tile (ng, tile, ibry,                       &
     &                            LBij, UBij, 1, N(ng),                 &
     &                            VnormVobc(:,:,ibry))
#  ifdef DISTRIBUTE
!
            Bwrk=RESHAPE(VnormVobc(:,:,ibry), (/IJKlen/))
            CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
            ic=0
            DO k=1,N(ng)
              DO ib=LBij,UBij
                ic=ic+1
                VnormVobc(ib,k,ibry)=Bwrk(ic)
              END DO
            END DO
#  endif
          END IF BRY_COMPUTE_V3D
        END DO BRY_LOOP_V3D
!
!  Write out into output NetCDF file.
!
        IF (ANY(CnormB(isVvel,:))) THEN
          IDmeta=idSbry(isVvel)
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_put_fvar (ng, iTLM, ncname,                   &
     &                              Vname(1,IDmeta),                    &
     &                              VnormVobc(LBij:,:,:),               &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_put_fvar (ng, iTLM, ncname,               &
     &                                  Vname(1,IDmeta),                &
     &                                  VnormVobc(LBij:,:,:),           &
     &                         start = (/1,1,1,NRM(ifile,ng)%Rindex/),  &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D boundary normalization at RHO-points.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          DO itrc=1,NT(ng)
            ifield=isTvar(itrc)
            IF (ANY(CnormB(ifield,:))) THEN
              Lsame=.TRUE.
              EXIT
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &            '3D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF
!
        BRY_TRACER_LOOP : DO itrc=1,NT(ng)
          VnormRobc=Aspv
          ifield=isTvar(itrc)
          BRY_LOOP_R3D : DO ibry=1,4
            BRY_COMPUTE_R3D : IF (CnormB(ifield,ibry)) THEN
              BRY_MS_LOOP_R3D : DO ns=1,Nscale(ng)
                VscaleB=0.0_r8
                B3davg=0.0_r8
                B3dsqr=0.0_r8
                B3d=0.0_r8
!
                SELECT CASE (ibry)
                  CASE (iwest, ieast)
                    i=BOUNDS(ng)%edge(ibry,r2dvar)
                    IF (Lconvolve(ibry)) THEN
                      DO j=JstrT,JendT
                        cff=GRID(ng)%on_r(i,j)
                        DO k=1,N(ng)
                          VscaleB(j,k)=1.0_r8/                          &
     &                                 SQRT(cff*GRID(ng)%Hz(i,j,k))
                        END DO
                      END DO
                    END IF
                  CASE (isouth, inorth)
                    j=BOUNDS(ng)%edge(ibry,r2dvar)
                    IF (Lconvolve(ibry)) THEN
                      DO i=IstrT,IendT
                        cff=GRID(ng)%om_r(i,j)
                        DO k=1,N(ng)
                          VscaleB(i,k)=1.0_r8/                          &
     &                                 SQRT(cff*GRID(ng)%Hz(i,j,k))
                        END DO
                      END DO
                    END IF
                END SELECT
!
                BRY_RANDOM_R3D : DO iter=1,Nrandom
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      CALL white_noise3d_bry (ng, tile, iTLM, ibry,     &
     &                                        Rscheme(ng),              &
     &                                        JstrR, JendR,             &
     &                                        LBij, UBij, 1, N(ng),     &
     &                                        Bmin, Bmax, B3d)
                      DO k=1,N(ng)
                        DO j=JstrT,JendT
                          B3d(j,k)=B3d(j,k)*VscaleB(j,k)
                        END DO
                      END DO
                    CASE (isouth, inorth)
                      CALL white_noise3d_bry (ng, tile, iTLM, ibry,     &
     &                                        Rscheme(ng),              &
     &                                        IstrR, IendR,             &
     &                                        LBij, UBij, 1, N(ng),     &
     &                                        Bmin, Bmax, B3d)
                      DO k=1,N(ng)
                        DO i=IstrT,IendT
                          B3d(i,k)=B3d(i,k)*VscaleB(i,k)
                        END DO
                      END DO
                  END SELECT
!
!  Implicit tangent linear horizontal convolution, CG/CI solver.
!
                  CALL self%tl_CI_b2d (ng, tile, iTLM, ifield, ibry,    &
     &                                 r3dvar,                          &
     &                                 ns, NiterCI(ns,ng), ifac,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 B3d)
!
!  Implicit tangent linear vertical convolution.
!
                  CALL self%tl_bry_Vdiff (ng, tile, iTLM, ifield, ibry, &
     &                                    r3dvar,                       &
     &                                    NVstepsB(ibry,ifield)/ifac,   &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    LBij, UBij,                   &
     &                                    DTsizeVB(ibry,ifield), Kv,    &
     &                                    B3d)
!
                  IF (Lconvolve(ibry)) THEN
                    SELECT CASE (ibry)
                      CASE (iwest, ieast)
                        DO k=1,N(ng)
                          DO j=Jstr,Jend
                            B3davg(j,k)=B3davg(j,k)+B3d(j,k)
                            B3dsqr(j,k)=B3dsqr(j,k)+B3d(j,k)*B3d(j,k)
                          END DO
                        END DO
                      CASE (isouth, inorth)
                        DO k=1,N(ng)
                          DO i=Istr,Iend
                            B3davg(i,k)=B3davg(i,k)+B3d(i,k)
                            B3dsqr(i,k)=B3dsqr(i,k)+B3d(i,k)*B3d(i,k)
                          END DO
                        END DO
                    END SELECT
                  END IF
                END DO BRY_RANDOM_R3D
!
!  Set normalization factors.
!
                IF (Lconvolve(ibry)) THEN
                  SELECT CASE (ibry)
                    CASE (iwest, ieast)
                      i=BOUNDS(ng)%edge(ibry,r2dvar)
                      DO k=1,N(ng)
                        DO j=Jstr,Jend
                          Bavg=FacAvg*B3davg(j,k)
                          Bsqr=FacAvg*B3dsqr(j,k)
#  ifdef MASKING
                          IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                            VnormRobc(j,k,ibry,itrc)=                   &
     &                                       VnormRobc(j,k,ibry,itrc)+  &
     &                                       self%Bwgt(ifield,ns)/      &
     &                                       SQRT(Bsqr)

                          ELSE
                            VnormRobc(j,k,ibry,itrc)=0.0_r8
                          END IF
#  else
                          VnormRobc(j,k,ibry,itrc)=                     &
     &                                       VnormRobc(j,k,ibry,itrc)+  &
     &                                       self%Bwgt(ifield,ns)/      &
     &                                       SQRT(Bsqr)
#  endif
                        END DO
                      END DO
                    CASE (isouth, inorth)
                      j=BOUNDS(ng)%edge(ibry,r2dvar)
                      DO k=1,N(ng)
                        DO i=Istr,Iend
                          Bavg=FacAvg*B3davg(i,k)
                          Bsqr=FacAvg*B3dsqr(i,k)
#  ifdef MASKING
                          IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                            VnormRobc(i,k,ibry,itrc)=                   &
       &                                     VnormRobc(i,k,ibry,itrc)+  &
       &                                     self%Bwgt(ifield,ns)/      &
       &                                     SQRT(Bsqr)
                          ELSE
                            VnormRobc(i,k,ibry,itrc)=0.0_r8
                          END IF
#  else
                          VnormRobc(i,k,ibry,itrc)=                     &
       &                                     VnormRobc(i,k,ibry,itrc)+  &
       &                                     self%Bwgt(ifield,ns)/      &
       &                                     SQRT(Bsqr)
#  endif
                        END DO
                      END DO
                  END SELECT
                END IF
              END DO BRY_MS_LOOP_R3D
!
!  Exchange boundary data.
!
              CALL bc_r3d_bry_tile (ng, tile, ibry,                     &
     &                              LBij, UBij, 1, N(ng),               &
     &                              VnormRobc(:,:,ibry,itrc))
#  ifdef DISTRIBUTE
!
              Bwrk=RESHAPE(VnormRobc(:,:,ibry,itrc), (/IJKlen/))
              CALL mp_collect (ng, iTLM, IJKlen, Aspv, Bwrk)
              ic=0
              DO k=1,N(ng)
                DO ib=LBij,UBij
                  ic=ic+1
                  VnormRobc(ib,k,ibry,itrc)=Bwrk(ic)
                END DO
              END DO
#  endif
            END IF BRY_COMPUTE_R3D
          END DO BRY_LOOP_R3D
!
!  Write out into output NetCDF file.
!
          IF (ANY(CnormB(ifield,:))) THEN
            IDmeta=idSbry(isTvar(itrc))
!
            SELECT CASE (NRM(ifile,ng)%IOtype)
              CASE (io_nf90)
                CALL netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,IDmeta),                  &
     &                                VnormRobc(LBij:,:,:,itrc),        &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(IDmeta))

#  if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)
                CALL pio_netcdf_put_fvar (ng, iTLM, ncname,             &
     &                                    Vname(1,IDmeta),              &
     &                                    VnormRobc(LBij:,:,:,itrc),    &
     &                         start =(/1,1,1,NRM(ifile,ng)%Rindex/),   &
     &                         total = (/IJlen,N(ng),4,1/),             &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(IDmeta)%vd)
#  endif
            END SELECT
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF

        END DO BRY_TRACER_LOOP

# endif /* SOLVE3D */
!
!  Synchronize open boundaries normalization NetCDF file to disk to
!  allow other processes to access data immediately after it is
!  written.
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_sync (ng, iTLM, ncname,                         &
     &                        NRM(ifile,ng)%ncid)
# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_sync (ng, iTLM, ncname,                     &
     &                            NRM(ifile,ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF BRY_LWRITE_NRM
#endif /* ADJUST_BOUNDARY */

#if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Compute surface forcing multiscale background-error covariance, B,
!  normalization factors using the randomization approach of Fisher and
!  Courtier (1995).
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      ifile=4
      FRC_LWRITE_NRM : IF (LwrtNRM(ifile,ng)) THEN
        rec=1
        Text='surface forcing'
!
!  Set randomization summation factors.
!
        FacAvg=1.0_r8/REAL(Nrandom,r8)
        FacSqr=SQRT(REAL(Nrandom,r8))
!
!  Set time record index to write in normalization NetCDF file.
!
        ncname=NRM(ifile,ng)%name
        NRM(ifile,ng)%Rindex=NRM(ifile,ng)%Rindex+1
        NRM(ifile,ng)%Nrec=NRM(ifile,ng)%Nrec+1
!
!  Write out model time (s).
!
        SELECT CASE (NRM(ifile,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_put_fvar (ng, iTLM, ncname,                     &
     &                            Vname(1,idtime), my_time,             &
     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         ncid = NRM(ifile,ng)%ncid,               &
     &                         varid = NRM(ifile,ng)%Vid(idtime))

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_put_fvar (ng, iTLM, ncname,                 &
     &                                Vname(1,idtime), my_time,         &

     &                         start = (/NRM(ifile,ng)%Rindex/),        &
     &                         total = (/1/),                           &
     &                         pioFile = NRM(ifile,ng)%pioFile,         &
     &                         pioVar = NRM(ifile,ng)%pioVar(idtime)%vd)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_WSTRESS
!
!-----------------------------------------------------------------------
!  2D normalization at U-stress points.
!-----------------------------------------------------------------------
!
        MS_COMPUTE_SUS : IF (Cnorm(rec,isUstr)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &            '2D normalization factors at U-points'
            FLUSH (stdout)
          END IF
!
          MS_SUS_LOOP : DO ns=1,Nscale(ng)
            DO j=JstrT,JendT
              DO i=IstrP,IendT
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_u(i,j)*             &
     &                                  GRID(ng)%on_u(i,j))
              END DO
            END DO
!
            RANDOM_SUS : DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, u2dvar, Rscheme(ng),        &
     &                            Istr, IendR, JstrR, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
              DO j=JstrT,JendT
                DO i=IstrP,IendT
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
!
!  Apply lateral boundary conditions.
!
              CALL dabc_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            A2d)
!
!  Implicit tangent linear horizontal convolution, CG/CI solver.
!
              CALL self%tl_CI_2d (ng, tile, iTLM, isUstr, u2dvar,       &
     &                            ns, NiterCI(ns,ng), ifac, Lweak,      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            A2d)
!
              DO j=Jstr,Jend
                DO i=IstrU,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO RANDOM_SUS
!
!  Set normalization factors.
!
            DO j=Jstr,Jend
              DO i=IstrU,Iend
                Aavg=FacAvg*A2davg(i,j)
                Asqr=FacAvg*A2dsqr(i,j)
#   ifdef MASKING
                IF (GRID(ng)%umask(i,j).gt.0.0_r8) THEN
                  HnormSUS(i,j)=HnormSUS(i,j)+                          &
     &                          self%Bwgt(isUstr,ns)/SQRT(Asqr)
                ELSE
                  HnormSUS(i,j)=0.0_r8
                END IF
#   else
                HnormSUS(i,j)=HnormSUS(i,j)+                            &
     &                        self%Bwgt(isUstr,ns)/SQRT(Asqr)
#   endif
              END DO
            END DO
          END DO MS_SUS_LOOP
!
!  Exchange boundary data.
!
          CALL dabc_u2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSUS)
#  ifdef DISTRIBUTE
!
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSUS)
#  endif
!
!  Write out into output NetCDF file.
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idUsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idUsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              GRID(ng)%umask,                     &
#  endif
     &                              HnormSUS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idUsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_u2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_u2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idUsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idUsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             GRID(ng)%umask,                      &
#   endif
     &                             HnormSUS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF MS_COMPUTE_SUS
!
!-----------------------------------------------------------------------
!  2D normalization at V-stress points.
!-----------------------------------------------------------------------
!
        MS_COMPUTE_SVS : IF (Cnorm(rec,isVstr)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &            '2D normalization factors at V-points'
            FLUSH (stdout)
          END IF
!
          MS_SVS_LOOP : DO ns=1,Nscale(ng)
            DO j=JstrP,JendT
              DO i=IstrT,IendT
                A2davg(i,j)=0.0_r8
                A2dsqr(i,j)=0.0_r8
                Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_v(i,j)*             &
     &                                  GRID(ng)%on_v(i,j))
              END DO
            END DO
!
            RANDOM_SVS : DO iter=1,Nrandom
              CALL white_noise2d (ng, iTLM, v2dvar, Rscheme(ng),        &
     &                            IstrR, IendR, Jstr, JendR,            &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Amin, Amax, A2d)
!
              DO j=JstrP,JendT
                DO i=IstrT,IendT
                  A2d(i,j)=A2d(i,j)*Hscale(i,j)
                END DO
              END DO
!
!  Apply lateral boundary conditions.
!
              CALL dabc_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            A2d)
!
!  Implicit tangent linear horizontal convolution, CG/CI solver.
!
              CALL self%tl_CI_2d (ng, tile, iTLM, isVstr, v2dvar,       &
     &                            ns, NiterCI(ns,ng), ifac, Lweak,      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            A2d)
!
              DO j=JstrV,Jend
                DO i=Istr,Iend
                  A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                  A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                END DO
              END DO
            END DO RANDOM_SVS
!
!  Set normalization factors.
!
            DO j=JstrV,Jend
              DO i=Istr,Iend
                Aavg=FacAvg*A2davg(i,j)
                Asqr=FacAvg*A2dsqr(i,j)
#   ifdef MASKING
                IF (GRID(ng)%vmask(i,j).gt.0.0_r8) THEN
                  HnormSVS(i,j)=HnormSVS(i,j)+                          &
     &                          self%Bwgt(isVstr,ns)/SQRT(Asqr)
                ELSE
                  HnormSVS(i,j)=0.0_r8
                END IF
#   else
                HnormSVS(i,j)=HnormSVS(i,j)+                            &
     &                        self%Bwgt(isVstr,ns)/SQRT(Asqr)
#   endif
              END DO
            END DO
          END DO MS_SVS_LOOP
!
!  Exchange boundary data.
!
          CALL dabc_v2d_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        HnormSVS)
#  ifdef DISTRIBUTE
!
          CALL mp_exchange2d (ng, tile, iTLM, 1,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        HnormSVS)
#  endif
!
!  Write out into output NetCDF file.
!
          SELECT CASE (NRM(ifile,ng)%IOtype)
            CASE (io_nf90)
              CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,             &
     &                              LBi, UBi, LBj, UBj, idVsms,         &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idVsms),          &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                              GRID(ng)%vmask,                     &
#  endif
     &                              HnormSVS)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              IF (NRM(ifile,ng)%pioVar(idVsms)%dkind.eq.                &
     &            PIO_double) THEN
                ioDesc => ioDesc_dp_v2dvar(ng)
              ELSE
                ioDesc => ioDesc_sp_v2dvar(ng)
              END IF
              CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,              &
     &                             LBi, UBi, LBj, UBj, idVsms,          &
     &                             NRM(ifile,ng)%pioFile,               &
     &                             NRM(ifile,ng)%pioVar(idVsms),        &
     &                             NRM(ifile,ng)%Rindex,                &
     &                             ioDesc,                              &
#   ifdef MASKING
     &                             GRID(ng)%vmask,                      &
#   endif
     &                             HnormSVS)
#  endif
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF MS_COMPUTE_SVS
# endif /* ADJUST_WSTRESS */

# if defined ADJUST_STFLUX && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  2D normalization for surface tracer flux.
!-----------------------------------------------------------------------
!
        IF (Master) THEN
          Lsame=.FALSE.
          DO itrc=1,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              ifield=isTsur(itrc)
              IF (Cnorm(rec,ifield)) Lsame=.TRUE.
            END IF
          END DO
          IF (Lsame) THEN
            WRITE (stdout,20) TRIM(Text),                               &
     &            '2D normalization factors at RHO-points'
            FLUSH (stdout)
          END IF
        END IF
!
!  Check if the decorrelation scales for all the surface tracer fluxes
!  are different. If not, just compute the normalization factors for the
!  first tracer and assign the same value to the rest.  Recall that this
!  computation is very expensive.
!
        Ldiffer=.FALSE.
        DO ns=1,Nscale(ng)

# ifdef NONUNIFORM_SCALES
          Ldiffer=.TRUE.
# else
          DO itrc=2,NT(ng)
            IF ((HdecayX(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayX(rec,isTsur(itrc-1),ns,ng)).or.                 &
     &          (HdecayY(rec,isTsur(itrc  ),ns,ng).ne.                  &
     &           HdecayY(rec,isTsur(itrc-1),ns,ng))) THEN
              Ldiffer=.TRUE.
            END IF
          END DO
# endif
          IF (.not.Ldiffer) THEN
            Lsame=.TRUE.
            UBt=1
          ELSE
            Lsame=.FALSE.
            UBt=NT(ng)
          END IF
        END DO
!
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            Hscale(i,j)=1.0_r8/SQRT(GRID(ng)%om_r(i,j)*                 &
     &                              GRID(ng)%on_r(i,j))
          END DO
        END DO
!
        TRACER_FLUX_LOOP : DO itrc=1,UBt
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            MS_COMPUTE_STF : IF (Cnorm(rec,ifield)) THEN
              MS_STF_LOOP : DO ns=1,Nscale(ng)
                DO j=JstrT,JendT
                  DO i=IstrT,IendT
                    A2davg(i,j)=0.0_r8
                    A2dsqr(i,j)=0.0_r8
                  END DO
                END DO
!
                RANDOM_STF : DO iter=1,Nrandom
                  CALL white_noise2d (ng, iTLM, r2dvar, Rscheme(ng),    &
     &                                IstrR, IendR, JstrR, JendR,       &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Amin, Amax, A2d)
!
                  DO j=JstrT,JendT
                    DO i=IstrT,IendT
                      A2d(i,j)=A2d(i,j)*Hscale(i,j)
                    END DO
                  END DO
!
!  Apply lateral boundary conditions.
!
                  CALL dabc_r2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                A2d)
!
!  Implicit tangent linear horizontal convolution, CG/CI solver.
!
                  CALL self%tl_CI_2d (ng, tile, iTLM, ifield, r2dvar,   &
     &                                ns, NiterCI(ns,ng), ifac, Lweak,  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                A2d)
!
                  DO j=Jstr,Jend
                    DO i=Istr,Iend
                      A2davg(i,j)=A2davg(i,j)+A2d(i,j)
                      A2dsqr(i,j)=A2dsqr(i,j)+A2d(i,j)*A2d(i,j)
                    END DO
                  END DO
                END DO RANDOM_STF
!
!  Set normalization factors.
!
                DO j=Jstr,Jend
                  DO i=Istr,Iend
                    Aavg=FacAvg*A2davg(i,j)
                    Asqr=FacAvg*A2dsqr(i,j)
#   ifdef MASKING
                    IF (GRID(ng)%rmask(i,j).gt.0.0_r8) THEN
                      HnormSTF(i,j,itrc)=HnormSTF(i,j,itrc)+            &
     &                                   self%Bwgt(ifield,ns)/          &
     &                                   SQRT(Asqr)
                    ELSE
                      HnormSTF(i,j,itrc)=0.0_r8
                    END IF
#   else
                    HnormSTF(i,j,itrc)=HnormSTF(i,j,itrc)+              &
     &                                 self%Bwgt(ifield,ns)/            &
     &                                 SQRT(Asqr)
#   endif
                  END DO
                END DO
              END DO MS_STF_LOOP
            END IF MS_COMPUTE_STF
          END IF
        END DO TRACER_FLUX_LOOP
!
!  If same correlation parameters, replicate normalization factors
!  for the remaining tracer fields.
!
        IF (Lsame) THEN
          DO itrc=2,NT(ng)
            IF (Lstflux(itrc,ng)) THEN
              DO j=Jstr,Jend
                DO i=Istr,Iend
                  HnormSTF(i,j,itrc)=HnormSTF(i,j,1)
                END DO
              END DO
            END IF
          END DO
        END IF
!
!  Exchange boundary data.
!
        DO itrc=1,NT(ng)
          IF (Lstflux(itrc,ng)) THEN
            ifield=isTsur(itrc)
            IF (Cnorm(rec,ifield)) THEN
              CALL dabc_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            HnormSTF(:,:,itrc))
#  ifdef DISTRIBUTE
!
              CALL mp_exchange2d (ng, tile, iTLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            HnormSTF(:,:,itrc))
#  endif
!
!  Write out into output NetCDF file.
!
              SELECT CASE (NRM(ifile,ng)%IOtype)
                CASE (io_nf90)
                  CALL wrt_norm2d_nf90 (ng, tile, iTLM, ncname,         &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  idTsur(itrc),                   &
     &                              NRM(ifile,ng)%ncid,                 &
     &                              NRM(ifile,ng)%Vid(idTsur(itrc)),    &
     &                              NRM(ifile,ng)%Rindex,               &
#  ifdef MASKING
     &                                  GRID(ng)%rmask,                 &
#  endif
     &                                  HnormSTF(:,:,itrc))

#  if defined PIO_LIB && defined DISTRIBUTE
                CASE (io_pio)
                  IF (NRM(ifile,ng)%pioVar(idTsur(itrc))%dkind.eq.      &
     &                PIO_double) THEN
                    ioDesc => ioDesc_dp_r2dvar(ng)
                  ELSE
                    ioDesc => ioDesc_sp_r2dvar(ng)
                  END IF
                  CALL wrt_norm2d_pio (ng, tile, iTLM, ncname,          &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 idTsur(itrc),                    &
     &                              NRM(ifile,ng)%pioFile,              &
     &                              NRM(ifile,ng)%pioVar(idTsur(itrc)), &
     &                              NRM(ifile,ng)%Rindex,               &
     &                                 ioDesc,                          &
#   ifdef MASKING
     &                                 GRID(ng)%rmask,                  &
#   endif
     &                                 HnormSTF(:,:,itrc))
#  endif
              END SELECT
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
        END DO
# endif /* ADJUST_STFLUX */
      END IF FRC_LWRITE_NRM
#endif
!
      IF (Master) THEN
        WRITE (stdout,30)
      END IF
!
 10   FORMAT (/,' Error Covariance Factors: Randomization Method',/)
 20   FORMAT (4x,'Computing',1x,a,1x,a)
 30   FORMAT (/)
!
      RETURN
      END SUBROUTINE randomization_tile
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Support PRIVATE functions.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
!***********************************************************************
      FUNCTION dot_prod2d (ng, tile, model, ctype,                      &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     A1, A2) RESULT (DotProd)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,   intent(in   ) :: ng                  ! nested grid
      integer,   intent(in   ) :: tile                ! domain partition
      integer,   intent(in   ) :: model               ! kernel ID
      integer,   intent(in   ) :: ctype               ! C-grid type
      integer,   intent(in   ) :: LBi, UBi, LBj, UBj
      real (r8), intent(in   ) :: A1(LBi:,LBj:)
      real (r8), intent(in   ) :: A2(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer   :: IstrT, IstrP, IendT, Imin, Imax
      integer   :: JstrT, JstrP, JendT, Jmin, Jmax
      integer   :: NSUB, i, j
      real (r8) :: DotProd
      real (r8) :: cff, my_DotProd
#ifdef DISTRIBUTE
      character (len=3) :: op_handle
#endif
!
!-----------------------------------------------------------------------
!  It computes the dot product between two 2D tiled arrays.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IstrT=BOUNDS(ng)%IstrT(tile)   ! tile computational range
      IstrP=BOUNDS(ng)%IstrP(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrP=BOUNDS(ng)%JstrP(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
      Imin=IstrT
      Imax=IendT
      Jmin=JstrT
      Jmax=JendT
      SELECT CASE (ctype)
        CASE (u2dvar, u3dvar)
          Imin=IstrP
        CASE (v2dvar, v3dvar)
          Jmin=JstrP
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          cff=A1(i,j)*A2(i,j)
          my_DotProd=my_DotProd+cff
        END DO
      END DO
!
!  Perform parallel global reduction operation.
!
#ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
#endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END FUNCTION dot_prod2d

#ifdef SOLVE3D
!
!***********************************************************************
      FUNCTION dot_prod3d (ng, tile, model, ctype,                      &
     &                     LBi, UBi, LBj, UBj, LBk, UBk,                &
     &                     A1, A2) RESULT (DotProd)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,   intent(in   ) :: ng                  ! nested grid
      integer,   intent(in   ) :: tile                ! domain partition
      integer,   intent(in   ) :: model               ! kernel ID
      integer,   intent(in   ) :: ctype               ! C-grid type
      integer,   intent(in   ) :: LBi, UBi, LBj, UBj, LBk, UBk
      real (r8), intent(in   ) :: A1(LBi:,LBj:,LBk:)
      real (r8), intent(in   ) :: A2(LBi:,LBj:,LBk:)
!
!  Local variable declarations.
!
      integer   :: IstrT, IstrP, IendT, Imin, Imax
      integer   :: JstrT, JstrP, JendT, Jmin, Jmax
      integer   :: NSUB, i, j, k
      real (r8) :: DotProd
      real (r8) :: cff, my_DotProd
# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif
!
!-----------------------------------------------------------------------
!  It computes the dot product between two 2D tiled arrays.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IstrT=BOUNDS(ng)%IstrT(tile)   ! tile computational range
      IstrP=BOUNDS(ng)%IstrP(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrP=BOUNDS(ng)%JstrP(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
      Imin=IstrT
      Imax=IendT
      Jmin=JstrT
      Jmax=JendT
      SELECT CASE (ctype)
        CASE (u2dvar, u3dvar)
          Imin=IstrP
        CASE (v2dvar, v3dvar)
          Jmin=JstrP
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
      DO k=LBk,UBk
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            cff=A1(i,j,k)*A2(i,j,k)
            my_DotProd=my_DotProd+cff
          END DO
        END DO
      END DO
!
!  Perform parallel global reduction operation.
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
# endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END FUNCTION dot_prod3d
#endif

!
!***********************************************************************
      SUBROUTINE wrt_norm2d_nf90 (ng, tile, model, ncname,              &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            ifield, ncid, ncvarid, tindex,        &
#ifdef MASKING
     &                            Amask,                                &
#endif
     &                            A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: A(LBi:,LBj:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
#endif
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm2d_nf90"

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 2D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
        gfactor=-1
#else
        gfactor=1
#endif
!
!  Write out 2D normalization field.
!
      gtype=gfactor*Iinfo(1,ifield,ng)
      scale=1.0_dp
      status=nf_fwrite2d(ng, model, ncid, ifield,                       &
     &                   ncvarid, tindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#ifdef MASKING
     &                   Amask,                                         &
#endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex,      &
     &                              TRIM(ncname)
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM2D_NF90 - error while writing variable: ',a,  &
     &        /,19x 'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0,', file: ',a)
!
      RETURN
      END SUBROUTINE wrt_norm2d_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_norm2d_pio (ng, tile, model, ncname,               &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           ifield, pioFile, pioVar, tindex,       &
     &                           pioDesc,                               &
# ifdef MASKING
     &                           Amask,                                 &
# endif
     &                           A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ifield, tindex

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj)
# endif
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (IO_Desc_t),   intent(inout) :: pioDesc
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm2d_pio"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 2D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out 2D normalization field.
!
      scale=1.0_dp
      status=nf_fwrite2d(ng, model, pioFile, ifield,                    &
     &                   pioVar, tindex, pioDesc,                       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   Amask,                                         &
# endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL pio_netcdf_sync (ng, model, ncname, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex,      &
     &                              TRIM(ncname)
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM2D_PIO - error while writing variable: ',a,   &
     &        /,18x 'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0,', file: ',a)
!
      RETURN
      END SUBROUTINE wrt_norm2d_pio
#endif

!
!***********************************************************************
      SUBROUTINE wrt_norm3d_nf90 (ng, tile, model, ncname,              &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            ifield, ncid, ncvarid, tindex,        &
#ifdef MASKING
     &                            Amask,                                &
#endif
     &                            A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
#else
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
#endif
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm3d_nf90"

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 3D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
        gfactor=-1
#else
        gfactor=1
#endif
!
!  Write out 3D normalization field.
!
      gtype=gfactor*Iinfo(1,ifield,ng)
      scale=1.0_dp
      status=nf_fwrite3d(ng, model, ncid, ifield,                       &
     &                   ncvarid, tindex, gtype,                        &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, scale,           &
#ifdef MASKING
     &                   Amask,                                         &
#endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL netcdf_sync (ng, model, ncname, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex,      &
     &                              TRIM(ncname)
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM3D_NF90 - error while writing variable: ',a,  &
     &        /,19x,'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0,', file: ',a)
!
      RETURN
      END SUBROUTINE wrt_norm3d_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_norm3d_pio (ng, tile, model, ncname,               &
     &                           LBi, UBi, LBj, UBj, LBk, UBk,          &
     &                           ifield, pioFile, pioVar, tindex,       &
     &                           pioDesc,                               &
# ifdef MASKING
     &                           Amask,                                 &
# endif
     &                           A)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: ifield, tindex

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: A(LBi:,LBj:,LBk:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)
# endif
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (IO_Desc_t),   intent(inout) :: pioDesc
      TYPE (My_VarDesc),  intent(inout) :: pioVar
!
!  Local variable declarations.
!
      integer :: status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_norm3d_pio"

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Write out requested 3D field into normalization NetCDF file. Since
!  the computation of normalization coefficients is a very expensive
!  computation, synchronize file to disk.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out 3D normalization field.
!
      scale=1.0_dp
      status=nf_fwrite3d(ng, model, pioFile, ifield,                    &
     &                   pioVar, tindex, pioDesc,                       &
     &                   LBi, UBi, LBj, UBj, LBk, UBk, scale,           &
# ifdef MASKING
     &                   Amask,                                         &
# endif
     &                   A)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,ifield)), tindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Synchronize normalization NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!
      CALL pio_netcdf_sync (ng, model, ncname, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)), tindex,      &
     &                              TRIM(ncname)
      FLUSH (stdout)
!
  10  FORMAT (/,' WRT_NORM3D_PIO - error while writing variable: ',a,   &
     &        /,18x,'into normalization NetCDF file for time record: ', &
     &        i0)
  20  FORMAT (7x,'wrote  ',a, t21,'normalization factors into record ', &
     &        i0,', file: ',a)
!
      RETURN
      END SUBROUTINE wrt_norm3d_pio
#endif
      END MODULE normalization_mod
