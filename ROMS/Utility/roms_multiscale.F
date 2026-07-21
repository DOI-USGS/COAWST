#include "cppdefs.h"
      MODULE roms_multiscale_mod
#ifdef MULTI_SCALE_B
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
! The module implements routines and objects for modeling multiscale   !
! background-error covariance (B) through implicit diffusion-based     !
! operators, as outlined by Weaver et al. (2016, 2018) and Mirouze     !
! et al. (2016).                                                       !
!                                                                      !
! The implicit two-dimensional diffusion solver, applied level by      !
! level, utilizes Chebyshev iterations (CI), which require estimates   !
! of the maximum and minimum eigenvalues of the K-Laplacian operator.  !
! These eigenvalues are precomputed for a fixed grid and specified     !
! k-values (diffusion coefficient with units of length squared) using  !
! the Lanczos variant of the conjugate gradients (CG) iterative        !
! algorithm, initialized with random vectors, to obtain the necessary  !
! eigenvalue estimates.                                                !
!                                                                      !
! This combined CG and CI methodology is employed for modeling the     !
! implicit multiscale backgrond-error covariance.                      !
!                                                                      !
! Horizontal spatially varying isotropic or anisotropic correlation    !
! scales (HscaleX, HscaleY) can be specified for each variable in the  !
! control vector (zeta, u, v, T, S). When provided, these correlation  !
! scales are read from an input NetCDF file. Error correlations are    !
! assumed to be separable in the horizontal and vertical directions.   !
!                                                                      !
! References:                                                          !
!                                                                      !
! Weaver, A.T. and I. Mirouze, 2013: On the diffusion equation and     !
!   its application to isotropic and anisotropic correlation modeling  !
!   in variational assimilation, Q. J. R. Meteorol. Soc., 139, 242-260,!
!   doi:10.1002/qj.1955.                                               !
!                                                                      !
! Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation     !
!   operators based on an implicitly formulated diffusion equation     !
!   solved with the Chebyshev iteration, Q. J. R. Meteorol. Soc., 142, !
!   455-471, doi:10.1002/qj.2664.                                      !
!                                                                      !
! Weaver, A.T., S. Gurol, J. Tshimanga, M. Chrust, and A. Piacentini,  !
!   2018: "Time"-parallel diffusion-based correlation operators. Q. J. !
!   R. Meteorol. Soc., 144, 2067-2088, doi:10.1002/qj.3302.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netCDF
      USE mod_scalars
!
      USE ad_bc_2d_mod,        ONLY : ad_dabc_r2d_tile,                 &
     &                                ad_dabc_u2d_tile,                 &
     &                                ad_dabc_v2d_tile
# ifdef ADJUST_BOUNDARY
      USE ad_bc_bry2d_mod,     ONLY : ad_bc_r2d_bry_tile,               &
     &                                ad_bc_u2d_bry_tile,               &
     &                                ad_bc_v2d_bry_tile
# endif
# ifdef SOLVE3D
      USE ad_bc_3d_mod,        ONLY : ad_dabc_r3d_tile,                 &
     &                                ad_dabc_u3d_tile,                 &
     &                                ad_dabc_v3d_tile
#  ifdef ADJUST_BOUNDARY
      USE ad_bc_bry3d_mod,     ONLY : ad_bc_r3d_bry_tile,               &
     &                                ad_bc_u3d_bry_tile,               &
     &                                ad_bc_v3d_bry_tile
#  endif
# endif
      USE bc_2d_mod,           ONLY : dabc_r2d_tile,                    &
     &                                dabc_u2d_tile,                    &
     &                                dabc_v2d_tile
# ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod,        ONLY : bc_r2d_bry_tile,                  &
     &                                bc_u2d_bry_tile,                  &
     &                                bc_v2d_bry_tile
# endif
# ifdef SOLVE3D
      USE bc_3d_mod,           ONLY : dabc_r3d_tile,                    &
     &                                dabc_u3d_tile,                    &
     &                                dabc_v3d_tile
#  ifdef ADJUST_BOUNDARY
      USE bc_bry3d_mod,        ONLY : bc_r3d_bry_tile,                  &
     &                                bc_u3d_bry_tile,                  &
     &                                bc_v3d_bry_tile
#  endif
# endif
      USE lapack_mod,          ONLY : DSTEQR
# ifdef DISTRIBUTE
      USE mp_exchange_mod,     ONLY : ad_mp_exchange2d
#  ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,     ONLY : ad_mp_exchange2d_bry
#  endif
#  ifdef SOLVE3D
      USE mp_exchange_mod,     ONLY : ad_mp_exchange3d
#   ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,     ONLY : ad_mp_exchange3d_bry
#   endif

#  endif
      USE mp_exchange_mod,     ONLY : mp_exchange2d
#  ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,     ONLY : mp_exchange2d_bry
#  endif
#  ifdef SOLVE3D
      USE mp_exchange_mod,     ONLY : mp_exchange3d
#   ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,     ONLY : mp_exchange3d_bry
#   endif
#  endif
      USE distribute_mod,      ONLY : mp_reduce
# endif
      USE nf_fread3d_mod,      ONLY : nf90_fread3d
      USE strings_mod,         ONLY : FoundError,                       &
     &                                find_string
      USE tracer_metadata_mod, ONLY : tracer_index
!
      implicit none
!
!-----------------------------------------------------------------------
!  ROMS multiscale background-error covariance CLASS object.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: multiscale
!
!  Kernel variables lower and upper bound array dimensions, and parallel
!  tile partition indices.
!
        integer :: tile                           ! domain partition
        integer :: LBi,  UBi,  LBj,  UBj          ! array bounds
# ifdef ADJUST_BOUNDARY
        integer :: LBij, UBij                     ! boundary bounds
# endif
!
!  Multiscale background covariance matrix parameters:
!
!    Nscale   number of correlation scales to be combined
!    Mlap     number of K-Laplacian smoother applications: implicit
!             M iterations representing Matern d-dimensional correlation
!             functions approaching a Gaussian in the asymptotic limit
!    NiterCG  number of CG iterations to obtain each scale eigenvalue
!    NiterCI  number of implicit Chebyshev solver iterations
!    Bwgt     weights assigned to each correlation function
!
        integer              :: Nscale
        integer, allocatable :: Mlap(:,:)         ! [MstateVar, Nscale]
        integer, allocatable :: NiterCG(:)        ! [Nscale]
        integer, allocatable :: NiterCI(:)        ! [Nscale]
!
        real(r8), allocatable :: Bwgt(:,:)        ! [MstateVar, Nscale]
!
!  The extrema eigenvalues of the K-Laplacian for each 2D variable in
!  the control vector as functions of both the number of iterations
!  and specified correlation scales, [Mlap,Nscale,range]. The range
!  dimension has a value of 2, where 1 corresponds to the minimum and
!  2 to the maximum eigenvalues.
!
        real(r8), pointer :: zeta_eigen(:,:,:)    ! free surface
        real(r8), pointer :: ubar_eigen(:,:,:)    ! 2D U-momentum
        real(r8), pointer :: vbar_eigen(:,:,:)    ! 2D V-momentum

# ifdef ADJUST_WSTRESS
!
        real(r8), pointer :: sustr_eigen(:,:,:)   ! surface U-stress
        real(r8), pointer :: svstr_eigen(:,:,:)   ! surface V-stress
# endif

# ifdef SOLVE3D
!
!  The extrema eigenvalues of the K-Laplacian for each 3D  variable in
!  the control vector as functions of both the number of iterations,
!  number of vertical level, and the specified correlation scales,
!  [Mlap,N,Nscale,range] and [Mlap,N,Nscale,range,NT]. The range
!  dimension has a value of 2, where 1 corresponds to the minimum and
!  2 to the maximum eigenvalues.
!
        real(r8), pointer :: u_eigen(:,:,:,:)      ! 3D U-momentum
        real(r8), pointer :: v_eigen(:,:,:,:)      ! 3D V-momentum
        real(r8), pointer :: t_eigen(:,:,:,:,:)    ! tracer variables

#  ifdef ADJUST_STFLUX
!
        real(r8), pointer :: stflux_eigen(:,:,:,:) ! surf. Tracer flux
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!  The extrema eigenvalues of the K-Laplacian for each 2D variable
!  lateral boundary adjustment in the control vector as functions of
!  number of iterations, boundary index, and the specified correlation
!  scales, [Mlap,4,Nscale,range]. The range dimension has a value
!  of 2, where 1 corresponds to the minimum and 2 to the maximum
!  eigenvalues.
!
        real(r8), pointer :: zeta_obc_eigen(:,:,:,:)  ! free surface
        real(r8), pointer :: ubar_obc_eigen(:,:,:,:)  ! 2D U-momentum
        real(r8), pointer :: vbar_obc_eigen(:,:,:,:)  ! 2D V-momentum

#  ifdef SOLVE3D
!
!  The eextrema eigenvalues of the K-Laplacian for each 3D  variable
!  lateral boundary adjustment in the control vector as functions of
!  number of iterations, vertical levels, boundary index, and the
!  specified correlation scales, [Mlap,4,Nscale,range] and
!  [Mlap,N,4,Nscale,range,NT]. The range dimension has a value
!  of 2, where 1 corresponds to the minimum and 2 to the maximum
!  eigenvalues.
!
        real(r8), pointer :: u_obc_eigen(:,:,:,:,:)   ! 3D U-momentum
        real(r8), pointer :: v_obc_eigen(:,:,:,:,:)   ! 3D V-momentum
        real(r8), pointer :: t_obc_eigen(:,:,:,:,:,:) ! tracers
#  endif
# endif
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 2D variables in control vector [LBi:UBi,LBj:UBj].
!
        real(r8), allocatable :: cg2d_p(:,:)         ! direction, p(k)
        real(r8), allocatable :: cg2d_q(:,:)         ! auxiliary, A*p(k)
        real(r8), allocatable :: cg2d_r(:,:)         ! residual,  r(k)
        real(r8), allocatable :: cg2d_x(:,:)         ! estimate,  x(k)

# ifdef ADJUST_BOUNDARY
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 1D lateral boundary adjustment [LBij:UBij].
!
        real(r8), allocatable :: cgB1d_p(:)          ! direction, p(k)
        real(r8), allocatable :: cgB1d_q(:)          ! auxiliary, A*p(k)
        real(r8), allocatable :: cgB1d_r(:)          ! residual,  r(k)
        real(r8), allocatable :: cgB1d_x(:)          ! estimate,  x(k)

#  ifdef SOLVE3D
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 2D lateral boundary adjustment [LBij:UBij,N(ng)].
!
        real(r8), allocatable :: cgB2d_p(:,:)        ! direction, p(k)
        real(r8), allocatable :: cgB2d_q(:,:)        ! auxiliary, A*p(k)
        real(r8), allocatable :: cgB2d_r(:,:)        ! residual,  r(k)
        real(r8), allocatable :: cgB2d_x(:,:)        ! estimate,  x(k)
#  endif
# endif

# ifdef SOLVE3D
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 3D variables in control vector [LBi:UBi,LBj:UBj,N(ng)].
!
        real(r8), allocatable :: cg3d_p(:,:,:)       ! direction, p(k)
        real(r8), allocatable :: cg3d_q(:,:,:)       ! auxiliary, A*p(k)
        real(r8), allocatable :: cg3d_r(:,:,:)       ! residual,  r(k)
        real(r8), allocatable :: cg3d_x(:,:,:)       ! estimate,  x(k)
# endif
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  2D variables in control vector, [LBi:UBi,LBj:UBj].
!
        real(r8), allocatable :: ci2d_p(:,:)         ! direction, p(k)
        real(r8), allocatable :: ci2d_q(:,:)         ! auxiliary, A*p(k)
        real(r8), allocatable :: ci2d_r(:,:)         ! residual,  r(k)
        real(r8), allocatable :: ci2d_x(:,:)         ! estimate,  x(k)

# ifdef ADJUST_BOUNDARY
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  lateral 1D boundary adjustment in control vector, [LBij:UBij].
!
        real(r8), allocatable :: ciB1d_p(:)          ! direction, p(k)
        real(r8), allocatable :: ciB1d_q(:)          ! auxiliary, A*p(k)
        real(r8), allocatable :: ciB1d_r(:)          ! residual,  r(k)
        real(r8), allocatable :: ciB1d_x(:)          ! estimate,  x(k)

#  ifdef SOLVE3D
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  lateral 2D boundary adjustment in control vector, [LBij:UBij].
!
        real(r8), allocatable :: ciB2d_p(:,:)        ! direction, p(k)
        real(r8), allocatable :: ciB2d_q(:,:)        ! auxiliary, A*p(k)
        real(r8), allocatable :: ciB2d_r(:,:)        ! residual,  r(k)
        real(r8), allocatable :: ciB2d_x(:,:)        ! estimate,  x(k)
#  endif
# endif

# ifdef SOLVE3D
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  3D variable in control vector, [LBi:UBi,LBj:UBj,N(ng)].
!
        real(r8), allocatable :: ci3d_p(:,:,:)       ! direction, p(k)
        real(r8), allocatable :: ci3d_q(:,:,:)       ! auxiliary, A*p(k)
        real(r8), allocatable :: ci3d_r(:,:,:)       ! residual,  r(k)
        real(r8), allocatable :: ci3d_x(:,:,:)       ! estimate,  x(k)
# endif

# ifdef NONUNIFORM_SCALES
!
!  Spatially varying isotropic or anisotropic horizontal Daley Length
!  Scales (in meters) are specified for each correlation function and
!  are combined to represent multiscale background-error covariance
!  matrix effects. These length scales are sourced from a designated
!  input NetCDF file. The data assimilation algorithm requires
!  correlation scale values for each variable in the control vector.
!
!  Currently, horizontal scale values are applied uniformly across
!  all vertical levels in the 3D control variables because determining
!  unique values for each level is complex. This strategy also reduces
!  memory requirements in large-scale applications.
!
!  Users are responsible for selecting appropriate scale values, such
!  as the first Rossby radius of deformation, or alternative
!  methodologies.
!
!  At present, these correlation scales are structured as
!    (lon, lat, axis, Nscale)      for regular variables and
!    (lon, lat, axis, Nscale, NT)  for tracer-type variables.
!  The "axis" dimension has a value of 2, where 1 corresponds to the
!  x-axis (abscissa) and 2 to the y-axis (ordinate) of the anisotropic
!  correlation.
!
        real(r8), pointer :: zeta_Bcorr(:,:,:,:)      ! free surface
        real(r8), pointer :: ubar_Bcorr(:,:,:,:)      ! 2D U-momentum
        real(r8), pointer :: vbar_Bcorr(:,:,:,:)      ! 2D V-momentum

#  ifdef ADJUST_WSTRESS
        real(r8), pointer :: sustr_Bcorr(:,:,:,:)     ! surface U-stress
        real(r8), pointer :: svstr_Bcorr(:,:,:,:)     ! surface V-stress
#  endif

#  ifdef SOLVE3D
        real(r8), pointer :: u_Bcorr(:,:,:,:)         ! 3D U-momentum
        real(r8), pointer :: v_Bcorr(:,:,:,:)         ! 3D V-momentum
        real(r8), pointer :: t_Bcorr(:,:,:,:,:)       ! tracer variables
#   ifdef ADJUST_STFLUX
        real(r8), pointer :: stflux_Bcorr(:,:,:,:,:)  ! surf tracer flux
#   endif
#  endif
# endif
!
        CONTAINS
!
!  Implicit Chebyshev Iterations (CI) solver.
!
        PROCEDURE :: tl_CI_2d      => multiscale_CI_2d_tl
        PROCEDURE :: ad_CI_2d      => multiscale_CI_2d_ad
# ifdef SOLVE3D
        PROCEDURE :: tl_CI_3d      => multiscale_CI_3d_tl
        PROCEDURE :: ad_CI_3d      => multiscale_CI_3d_ad
# endif
# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_CI_b1d     => multiscale_CI_b1d_tl
        PROCEDURE :: ad_CI_b1d     => multiscale_CI_b1d_ad
#  ifdef SOLVE3D
        PROCEDURE :: tl_CI_b2d     => multiscale_CI_b2d_tl
        PROCEDURE :: ad_CI_b2d     => multiscale_CI_b2d_ad
#  endif
# endif
!
!  Implicit Conjugate Gradient (CG) solver.
!
        PROCEDURE :: tl_CG_2d      => multiscale_CG_2d_tl
# ifdef SOLVE3D
        PROCEDURE :: tl_CG_3d      => multiscale_CG_3d_tl
# endif
# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_CG_b1d     => multiscale_CG_b1d_tl
#  ifdef SOLVE3D
        PROCEDURE :: tl_CG_b2d     => multiscale_CG_b2d_tl
#  endif
# endif

# ifdef SOLVE3D
!
!  Implicit Vertical Diffusion solver.
!
        PROCEDURE :: tl_Vdiff_r3d  => multiscale_Vdiff_r3d_tl
        PROCEDURE :: ad_Vdiff_r3d  => multiscale_Vdiff_r3d_ad
        PROCEDURE :: tl_Vdiff_u3d  => multiscale_Vdiff_u3d_tl
        PROCEDURE :: ad_Vdiff_u3d  => multiscale_Vdiff_u3d_ad
        PROCEDURE :: tl_Vdiff_v3d  => multiscale_Vdiff_v3d_tl
        PROCEDURE :: ad_Vdiff_v3d  => multiscale_Vdiff_v3d_ad
#  ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_bry_Vdiff  => multiscale_bry_Vdiff_tl
        PROCEDURE :: ad_bry_Vdiff  => multiscale_bry_Vdiff_ad
#  endif
# endif
!
!  K-Laplacian Operator, [1 + Del(K*Del)].
!
        PROCEDURE :: tl_Klap_r2d   => multiscale_Klap_r2d_tl
        PROCEDURE :: ad_Klap_r2d   => multiscale_Klap_r2d_ad
        PROCEDURE :: tl_Klap_u2d   => multiscale_Klap_u2d_tl
        PROCEDURE :: ad_Klap_u2d   => multiscale_Klap_u2d_ad
        PROCEDURE :: tl_Klap_v2d   => multiscale_Klap_v2d_tl
        PROCEDURE :: ad_Klap_v2d   => multiscale_Klap_v2d_ad

# ifdef SOLVE3D
!
        PROCEDURE :: tl_Klap_r3d   => multiscale_Klap_r3d_tl
        PROCEDURE :: ad_Klap_r3d   => multiscale_Klap_r3d_ad
        PROCEDURE :: tl_Klap_u3d   => multiscale_Klap_u3d_tl
        PROCEDURE :: ad_Klap_u3d   => multiscale_Klap_u3d_ad
        PROCEDURE :: tl_Klap_v3d   => multiscale_Klap_v3d_tl
        PROCEDURE :: ad_Klap_v3d   => multiscale_Klap_v3d_ad
# endif
# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_Klap_b1d   => multiscale_Klap_b1d_tl
        PROCEDURE :: ad_Klap_b1d   => multiscale_Klap_b1d_ad
#  ifdef SOLVE3D
!
        PROCEDURE :: tl_Klap_b2d   => multiscale_Klap_b2d_tl
        PROCEDURE :: ad_Klap_b2d   => multiscale_Klap_b2d_ad
#  endif
# endif

# ifdef NONUNIFORM_SCALES
!
!  Process spatially-varying horizontal decorrelation length scales.
!
        PROCEDURE :: get_scales    => multiscale_get_scales
# endif
!
      END TYPE multiscale
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      TYPE (multiscale), allocatable :: MSB(:)
!
      PUBLIC  :: deallocate_multiscale
      PUBLIC  :: initialize_multiscale
!
# ifdef MULTI_SCALE_DEBUG
      PRIVATE :: multiscale_print
# endif
# ifdef ADJUST_BOUNDARY
      PRIVATE :: dot_prod_1d
# endif
      PRIVATE :: dot_prod_2d
!
      CONTAINS
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  It allocates and initializes the background-error covariance (B)
!  multiscale object.
!
      SUBROUTINE initialize_multiscale (ng, tile, model)
!
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
!
!  Allocate multiscale background-error covariance object.
!
      IF (ng.eq.1) THEN
        IF (.not.allocated(MSB)) allocate ( MSB(Ngrids) )
      END IF
!
!  Create mutiscale background-error convariance object.
!
      CALL multiscale_create (MSB(ng), ng, tile)
!
      RETURN
      END SUBROUTINE initialize_multiscale
!
!-----------------------------------------------------------------------
!  It deallocates the multiscale background-error covariance object.
!
      SUBROUTINE deallocate_multiscale (ng)
!
      integer,            intent(in   ) :: ng       ! nested grid
!
!  Deallocate background-error ocovariance derived-type structure.
!  Note that it also deallocates its internal arrays.
!
      IF (ng.eq.1) THEN
        IF (allocated(MSB)) deallocate (MSB)
      END IF
!
      RETURN
      END SUBROUTINE deallocate_multiscale
!
!-----------------------------------------------------------------------
!  It allocates and initializes the multiscale background-error
!  covariance object arrays used by the implicit diffusion operators.
!
      SUBROUTINE multiscale_create (self, ng, tile)
!
      CLASS (multiscale), intent(inout) :: self     ! multiscale object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
!
      integer                           :: Mlap, is, ms
      integer                           :: LBi, UBi, LBj, UBj
# ifdef ADJUST_BOUNDARY
      integer                           :: LBij, UBij
# endif
!
      real (r8), parameter              :: IniVal = 0.0_r8
      real (r8)                         :: size2d
!
!  Set tiled arrays bounds.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
# ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
# endif
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Initiziale CLASS parameters.
!
      self%tile=tile
!
      self%LBi=LBi
      self%UBi=UBi
      self%LBj=LBj
      self%UBj=UBj
# ifdef ADJUST_BOUNDARY
      self%LBij=LBij
      self%UBij=UBij
# endif
!
      self%Nscale=Nscale(ng)
!
!  Allocate CLASS parameters.
!
      IF (.not.allocated(self % Mlap)) THEN
        allocate ( self % Mlap(MstateVar,Nscale(ng)) )
        DO ms=1,Nscale(ng)
          DO is=1,MstateVar
            self % Mlap(is,ms) = Nlap(is,ms,ng)
          END DO
        END DO
      END IF
!
      IF (.not.allocated(self % NiterCG)) THEN
        allocate ( self % NiterCG(Nscale(ng)) )
        DO ms=1,Nscale(ng)
          self % NiterCG(ms) = NiterCG(ms,ng)
        END DO
      END IF
!
      IF (.not.allocated(self % NiterCI)) THEN
        allocate ( self % NiterCI(Nscale(ng)) )
        DO ms=1,Nscale(ng)
          self % NiterCI(ms) = NiterCI(ms,ng)
        END DO
      END IF
!
      IF (.not.allocated(self % Bwgt)) THEN
        allocate ( self % Bwgt(MstateVar,Nscale(ng)) )
        DO ms=1,Nscale(ng)
          DO is=1,MstateVar
            self % Bwgt(is,ms) = Bwgt(is,ms,ng)
          END DO
        END DO
      END IF
!
      Mlap=MAXVAL(self % Mlap)
!
!  Allocate the maximum and minimum eigenvalues of the K-Laplacian.
!
      nullify (self % zeta_eigen)
      IF (.not.associated(self % zeta_eigen)) THEN
        allocate ( self % zeta_eigen(Mlap,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*Nscale(ng)*2,r8)
        self % zeta_eigen = IniVal
      END IF
!
      nullify (self % ubar_eigen)
      IF (.not.associated(self % ubar_eigen)) THEN
        allocate ( self % ubar_eigen(Mlap,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*Nscale(ng)*2,r8)
        self % ubar_eigen = IniVal
      END IF
!
      nullify (self % vbar_eigen)
      IF (.not.associated(self % vbar_eigen)) THEN
        allocate ( self % vbar_eigen(Mlap,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*Nscale(ng)*2,r8)
        self % vbar_eigen = IniVal
      END IF

# ifdef ADJUST_WSTRESS
!
      nullify (self % sustr_eigen)
      IF (.not.associated(self % sustr_eigen)) THEN
        allocate ( self % sustr_eigen(Mlap,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*Nscale(ng)*2,r8)
        self % sustr_eigen = IniVal
      END IF
!
      nullify (self % svstr_eigen)
      IF (.not.associated(self % svstr_eigen)) THEN
        allocate ( self % svstr_eigen(Mlap,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*Nscale(ng)*2,r8)
        self % svstr_eigen = IniVal
      END IF
# endif

# ifdef SOLVE3D
!
      nullify (self % u_eigen)
      IF (.not.associated(self % u_eigen)) THEN
        allocate ( self % u_eigen(Mlap,N(ng),Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*N(ng)*Nscale(ng)*2,r8)
        self % u_eigen = IniVal
      END IF
!
      nullify (self % v_eigen)
      IF (.not.associated(self % v_eigen)) THEN
        allocate ( self % v_eigen(Mlap,N(ng),Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*N(ng)*Nscale(ng)*2,r8)
        self % v_eigen = IniVal
      END IF
!
      nullify (self % t_eigen)
      IF (.not.associated(self % t_eigen)) THEN
        allocate ( self % t_eigen(Mlap,N(ng),Nscale(ng),2,NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*N(ng)*Nscale(ng)*2*NT(ng),r8)
        self % t_eigen = IniVal
      END IF

#  ifdef ADJUST_STFLUX
!
      nullify (self % stflux_eigen)
      IF (.not.associated(self % stflux_eigen)) THEN
        allocate ( self % stflux_eigen(Mlap,Nscale(ng),2,NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*Nscale(ng)*2*NT(ng),r8)
        self % stflux_eigen = IniVal
      END IF
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
      nullify (self % zeta_obc_eigen)
      IF (.not.associated(self % zeta_obc_eigen)) THEN
        allocate ( self % zeta_obc_eigen(Mlap,4,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*4*Nscale(ng)*2,r8)
        self % zeta_obc_eigen = IniVal
      END IF
!
      nullify (self % ubar_obc_eigen)
      IF (.not.associated(self % ubar_obc_eigen)) THEN
        allocate ( self % ubar_obc_eigen(Mlap,4,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*4*Nscale(ng)*2,r8)
        self % ubar_obc_eigen = IniVal
      END IF
!
      nullify (self % vbar_obc_eigen)
      IF (.not.associated(self % vbar_obc_eigen)) THEN
        allocate ( self % vbar_obc_eigen(Mlap,4,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*4*Nscale(ng)*2,r8)
        self % vbar_obc_eigen = IniVal
      END IF

#  ifdef SOLVE3D
!
      nullify (self % u_obc_eigen)
      IF (.not.associated(self % u_obc_eigen)) THEN
        allocate ( self % u_obc_eigen(Mlap,N(ng),4,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*N(ng)*4*Nscale(ng)*2,r8)
        self % u_obc_eigen = IniVal
      END IF
!
      nullify (self % v_obc_eigen)
      IF (.not.associated(self % v_obc_eigen)) THEN
        allocate ( self % v_obc_eigen(Mlap,N(ng),4,Nscale(ng),2) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*N(ng)*4*Nscale(ng)*2,r8)
        self % v_obc_eigen = IniVal
      END IF
!
      nullify (self % t_obc_eigen)
      IF (.not.associated(self % t_obc_eigen)) THEN
        allocate ( self % t_obc_eigen(Mlap,N(ng),4,Nscale(ng),2,        &
     &                                NT(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Mlap*N(ng)*4*Nscale(ng)*2*NT(ng),r8)
        self % t_obc_eigen = IniVal
      END IF
#  endif
# endif
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors.
!
      IF (.not.allocated(self % cg2d_p)) THEN
        allocate ( self % cg2d_p(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_p = IniVal
      END IF
      IF (.not.allocated(self % cg2d_q)) THEN
        allocate ( self % cg2d_q(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_q = IniVal
      END IF
      IF (.not.allocated(self % cg2d_r)) THEN
        allocate ( self % cg2d_r(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_r = IniVal
      END IF
      IF (.not.allocated(self % cg2d_x)) THEN
        allocate ( self % cg2d_x(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_x = IniVal
      END IF

# ifdef ADJUST_BOUNDARY
!
      IF (.not.allocated(self % cgB1d_p)) THEN
        allocate ( self % cgB1d_p(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_p = IniVal
      END IF
      IF (.not.allocated(self % cgB1d_q)) THEN
        allocate ( self % cgB1d_q(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_q = IniVal
      END IF
      IF (.not.allocated(self % cgB1d_r)) THEN
        allocate ( self % cgB1d_r(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_r = IniVal
      END IF
      IF (.not.allocated(self % cgB1d_x)) THEN
        allocate ( self % cgB1d_x(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_x = IniVal
      END IF

#  ifdef SOLVE3D
!
      IF (.not.allocated(self % cgB2d_p)) THEN
        allocate ( self % cgB2d_p(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_p = IniVal
      END IF
      IF (.not.allocated(self % cgB2d_q)) THEN
        allocate ( self % cgB2d_q(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_q = IniVal
      END IF
      IF (.not.allocated(self % cgB2d_r)) THEN
        allocate ( self % cgB2d_r(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_r = IniVal
      END IF
      IF (.not.allocated(self % cgB2d_x)) THEN
        allocate ( self % cgB2d_x(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_x = IniVal
      END IF
#  endif
# endif

# ifdef SOLVE3D
!
      IF (.not.allocated(self % cg3d_p)) THEN
        allocate ( self % cg3d_p(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_p = IniVal
      END IF
      IF (.not.allocated(self % cg3d_q)) THEN
        allocate ( self % cg3d_q(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_q = IniVal
      END IF
      IF (.not.allocated(self % cg3d_r)) THEN
        allocate ( self % cg3d_r(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_r = IniVal
      END IF
      IF (.not.allocated(self % cg3d_x)) THEN
        allocate ( self % cg3d_x(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_x = IniVal
      END IF
# endif
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors.
!
      IF (.not.allocated(self % ci2d_p)) THEN
        allocate ( self % ci2d_p(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_p = IniVal
      END IF
      IF (.not.allocated(self % ci2d_q)) THEN
        allocate ( self % ci2d_q(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_q = IniVal
      END IF
      IF (.not.allocated(self % ci2d_r)) THEN
        allocate ( self % ci2d_r(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_r = IniVal
      END IF
      IF (.not.allocated(self % ci2d_x)) THEN
        allocate ( self % ci2d_x(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_x = IniVal
      END IF

# ifdef ADJUST_BOUNDARY
!
      IF (.not.allocated(self % ciB1d_p)) THEN
        allocate ( self % ciB1d_p(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_p = IniVal
      END IF
      IF (.not.allocated(self % ciB1d_q)) THEN
        allocate ( self % ciB1d_q(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_q = IniVal
      END IF
      IF (.not.allocated(self % ciB1d_r)) THEN
        allocate ( self % ciB1d_r(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_r = IniVal
      END IF
      IF (.not.allocated(self % ciB1d_x)) THEN
        allocate ( self % ciB1d_x(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_x = IniVal
      END IF

#  ifdef SOLVE3D
!
      IF (.not.allocated(self % ciB2d_p)) THEN
        allocate ( self % ciB2d_p(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_p = IniVal
      END IF
      IF (.not.allocated(self % ciB2d_q)) THEN
        allocate ( self % ciB2d_q(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_q = IniVal
      END IF
      IF (.not.allocated(self % ciB2d_r)) THEN
        allocate ( self % ciB2d_r(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_r = IniVal
      END IF
      IF (.not.allocated(self % ciB2d_x)) THEN
        allocate ( self % ciB2d_x(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_x = IniVal
      END IF
#  endif
# endif

# ifdef SOLVE3D
!
      IF (.not.allocated(self % ci3d_p)) THEN
        allocate ( self % ci3d_p(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_p = IniVal
      END IF
      IF (.not.allocated(self % ci3d_q)) THEN
        allocate ( self % ci3d_q(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_q = IniVal
      END IF
      IF (.not.allocated(self % ci3d_r)) THEN
        allocate ( self % ci3d_r(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_r = IniVal
      END IF
      IF (.not.allocated(self % ci3d_x)) THEN
        allocate ( self % ci3d_x(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_x = IniVal
      END IF
# endif

# ifdef NONUNIFORM_SCALES
!
!  Spatially varying isotropic or anisotropic horizontal Daley length
!  scales.
!
      nullify (self % zeta_Bcorr)
      IF (.not.associated(self % zeta_Bcorr)) THEN
        allocate ( self % zeta_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % zeta_Bcorr = IniVal
      END IF
!
      nullify (self % ubar_Bcorr)
      IF (.not.associated(self % ubar_Bcorr)) THEN
        allocate ( self % ubar_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % ubar_Bcorr = IniVal
      END IF
!
      nullify (self % vbar_Bcorr)
      IF (.not.associated(self % vbar_Bcorr)) THEN
        allocate ( self % vbar_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % vbar_Bcorr = IniVal
      END IF

#  ifdef ADJUST_WSTRESS
!
      nullify (self % sustr_Bcorr)
      IF (.not.associated(self % sustr_Bcorr)) THEN
        allocate ( self % sustr_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % sustr_Bcorr = IniVal
      END IF
!
      nullify (self % svstr_Bcorr)
      IF (.not.associated(self % svstr_Bcorr)) THEN
        allocate ( self % svstr_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % svstr_Bcorr = IniVal
      END IF
#  endif

#  ifdef SOLVE3D
!
      nullify (self % u_Bcorr)
      IF (.not.associated(self % u_Bcorr)) THEN
        allocate ( self % u_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % u_Bcorr = IniVal
      END IF
!
      nullify (self % v_Bcorr)
      IF (.not.associated(self % v_Bcorr)) THEN
        allocate ( self % v_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng),r8)
        self % v_Bcorr = IniVal
      END IF
!
      nullify (self % t_Bcorr)
      IF (.not.associated(self % t_Bcorr)) THEN
        allocate ( self % t_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng),         &
     &                            NT(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng)*NT(ng),r8)
        self % t_Bcorr = IniVal
      END IF

#   ifdef ADJUST_STFLUX
!
      nullify (self % stflux_Bcorr)
      IF (.not.associated(self % stflux_Bcorr)) THEN
        allocate ( self % stflux_Bcorr(LBi:UBi,LBj:UBj,2,Nscale(ng),    &
     &                                 NT(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(2*Nscale(ng)*NT(ng),r8)
        self % stflux_Bcorr = IniVal
      END IF
#   endif
#  endif
# endif
!
      RETURN
      END SUBROUTINE multiscale_create
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! <><><> CLASS MULTISCALE:  Chebyshev Iterations Solver <><><><><><><><>
!
# include "multiscale_CIsolver.h"

!
! <><><> CLASS MULTISCALE:  Conjugate Gradient Solver <><><><><><><><><>
!
# include "multiscale_CGsolver.h"

# ifdef SOLVE3D
!
! <><><> CLASS MULTISCALE:  Vertical Diffusion Solver <><><><><><><><><>
!
#  include "multiscale_Vdiff.h"
# endif

!
! <><><> CLASS MULTISCALE:  K-Laplacian Operators <><><><><><><><><><><>
!
# include "multiscale_Klaplacian.h"

# ifdef NONUNIFORM_SCALES
!
! <><><> CLASS MULTISCALE:  Process Correlation Scales  <><><><><><><><>
!
!  Read and load spatially-varying, background-error correlation 
!  horizontal length scales (Daley, 1992) from input NetCDF file.
!
!  Currently, there is not depth dependency for the horizontal length
!  scales, so their values are replicated for each veritcal level
!  in the implicit pseudo-diffusion solvers.
!
      SUBROUTINE multiscale_get_scales (self, ng, tile, model)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
!
      logical                           :: foundit
!                            
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: is, ie, js, je
      integer                           :: IDmeta, Nscale, Vsize(4)
      integer                           :: gtype, i, itrc, ns
      integer                           :: ncid, status, varid
!
      real (r8), pointer                :: Bscale(:,:,:) => NULL()
#  ifdef MASKING
      real (r8), pointer                :: mask(:,:) => NULL()
#  endif
      real (r8)                         :: Fmin, Fmax, Fscl
      real (r8), dimension(4)           :: Fvalue
!
#  ifdef DISTRIBUTE
      character (len=3), dimension(4)   :: op_handle
#  endif
      character (len=40 )               :: Cname
      character (len=256)               :: ncname
      character (len=*), parameter      :: MyFile =                     &
     &  __FILE__//", multiscale_get_scales"
!
      SourceFile=MyFile
!
!  Initialize.
!
      LBi=self%LBi
      UBi=self%UBi
      LBj=self%LBj
      UBj=self%UBj
!
      is=LBi
      ie=UBi
      js=LBj
      je=UBj
!
      Nscale=self%Nscale
      Vsize=0
      Fscl=1.0_r8

#  ifdef DISTRIBUTE
!
      op_handle(1)='MIN'
      op_handle(2)='MAX'
      op_handle(3)='MIN'
      op_handle(4)='MAX'
#  endif
!
!  Open input NetCDF file.
!
      ncname=TRIM(SVC(ng)%name)
!
      CALL netcdf_open (ng, model, ncname, 1, SVC(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = SVC(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = SVC(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in horizontal correlation scales (lon, lat, axis, Nscale).
!  (Note: by default the pointer sections remap bounds begins at
!         index 1, so we need to specify the all bounds)
!
      IF (Master) THEN
        WRITE (stdout,20) 'Reading horizontal correlation length '//    &
                          'scales file: ', ng, TRIM(ncname)
      END IF
!
      CONTROL_VAR : DO i=1,MstateVar
        IF (icSvar(i).gt.0) THEN
          IDmeta=icSvar(i)
          Cname=TRIM(BscaleVarName(i))
          foundit=find_string(var_name, n_var, Cname, varid)
!
          MS_LOOP : DO ns=1,Nscale
            LREAD_CORRELATION : IF (foundit) THEN
              SELECT CASE (TRIM(Cname))
                CASE ('zeta_Bcorr')
                  Bscale(is:,js:,1:) => self%zeta_Bcorr(is:ie,js:je,    &
     &                                                  1:2,ns)
#  ifdef MASKING
                  mask => GRID(ng)%rmask
#  endif
                CASE ('ubar_Bcorr')
                  Bscale(is:,js:,1:) => self%ubar_Bcorr(is:ie,js:je,    &
     &                                                  1:2,ns)
#  ifdef MASKING
                  mask => GRID(ng)%umask
#  endif
                CASE ('vbar_Bcorr')
                  Bscale(is:,js:,1:) => self%vbar_Bcorr(is:ie,js:je,    &
     &                                                  1:2,ns)
#  ifdef MASKING
                  mask => GRID(ng)%vmask
#  endif
#  ifdef ADJUST_WSTRESS
                CASE ('sustr_Bcorr')
                  Bscale(is:,js:,1:) => self%sustr_Bcorr(is:ie,js:je,   &
     &                                                   1:2,ns)
#   ifdef MASKING
                  mask => GRID(ng)%umask
#   endif
                CASE ('svstr_Bcorr')
                  Bscale(is:,js:,1:) => self%svstr_Bcorr(is:ie,js:je,   &
     &                                                   1:2,ns)
#   ifdef MASKING
                  mask => GRID(ng)%vmask
#   endif
#  endif
#  ifdef SOLVE3D
                CASE ('u_Bcorr')
                  Bscale(is:,js:,1:) => self%u_Bcorr(is:ie,js:je,       &
     &                                               1:2,ns)
#   ifdef MASKING
                  mask => GRID(ng)%umask
#   endif
                CASE ('v_Bcorr')
                  Bscale(is:,js:,1:) => self%v_Bcorr(is:ie,js:je,       &
     &                                               1:2,ns)
#   ifdef MASKING
                  mask => GRID(ng)%vmask
#   endif
                CASE ('temp_Bcorr', 'salt_Bcorr')
                  itrc = tracer_index(TRIM(Cname))
                  Bscale(is:,js:,1:) => self%t_Bcorr(is:ie,js:je,       &
     &                                               1:2,ns,itrc)
#   ifdef MASKING
                  mask => GRID(ng)%rmask
#   endif
#   ifdef ADJUST_STFLUX
                CASE ('shflux_Bcorr', 'ssflux_Bcorr')
                  itrc = tracer_index(TRIM(Cname))
                  Bscale(is:,js:,1:) => self%stflux_Bcorr(is:ie,js:je,  &
     &                                                    1:2,ns,itrc)
#    ifdef MASKING
                  mask => GRID(ng)%rmask
#    endif
#   endif
#  endif
              END SELECT
!
              gtype=Iinfo(1,IDmeta,ng)
              status=nf90_fread3d(ng, model, ncname, SVC(ng)%ncid,      &
     &                            Cname, varid,                         &
     &                            1, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj, 1, 2,             &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            mask,                                 &
#  endif
     &                            Bscale)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,30) TRIM(Vname(1,IDmeta)), TRIM(ncname)
                END IF
                exit_flag=2
                ioerror=status
                RETURN              
              END IF
!
              Fvalue(1)=MINVAL(Bscale(:,:,1),                           &
     &                         MASK = Bscale(:,:,1) > 0.0_r8)
              Fvalue(2)=MAXVAL(Bscale(:,:,1))
              Fvalue(3)=MINVAL(Bscale(:,:,2),                           &
     &                         MASK = Bscale(:,:,2) > 0.0_r8)
              Fvalue(4)=MAXVAL(Bscale(:,:,2))
#  ifdef DISTRIBUTE
              CALL mp_reduce (ng, model, 4, Fvalue, op_handle)
#  endif
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(2,IDmeta)), ns,            &
     &                            Fvalue(1), Fvalue(2),                 &
     &                            Fvalue(3), Fvalue(4)
              END IF

#  ifdef DISTRIBUTE
!
!  Exchange boundary data.

              CALL mp_exchange3d (ng, tile, model, 1,                   &
     &                            LBi, UBi, LBj, UBj, 1, 2,             &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            Bscale)
#  endif
              nullify (Bscale)
            END IF LREAD_CORRELATION
          END DO MS_LOOP
        END IF
      END DO CONTROL_VAR        
!
  10  FORMAT (/,2x,'MULTISCALE_GET_SCALES - Unable to open input',      &
     &        ' NetCDF file: ',a)
  20  FORMAT (/,1x,'MULTISCALE_GET_SCALES - SVC: ',a,/,26x,             &
     &        '(Grid ',i2.2,', File: ',a,')')
  30  FORMAT (/,2x,'MULTISCALE_GET_SCALES - Error while reading',       &
     &        ' variable: ',2x,a,/,26x,'in input NetCDF file: 'a)
  40  FORMAT (23x,'- ',a,/,26x,'(Bscale = ', i2.2,                      &
     &            ', Xmin = ',1p,e15.8, ', Xmax = ',1p,e15.8,/, 38x,    &
     &            '  Ymin = ',1p,e15.8, ', Ymax = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE multiscale_get_scales
# endif

!
!::::::: PRIVATE FUNCTIONS :::::::::::::::::::::::::::::::::::::::::::::

# ifdef MULTI_SCALE_DEBUG
!
!-----------------------------------------------------------------------
!  It prints valuable diagnostic information about the convergence of
!  the CG and CI algorithms for tuning the implicit multiscale diffusion
!  operator.
!  
      FUNCTION multiscale_print (ng, ifield, ms, iterDiff, Ndiff,       &
     &                           iter, Niter, label, deps,              &
     &                           level, boundary) RESULT (status)
!
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: iterDiff  ! diffusion step
      integer,            intent(in   ) :: Ndiff     ! diff applications
      integer,            intent(in   ) :: iter      ! solver iteration
      integer,            intent(in   ) :: Niter     ! CG/CI iterations
      integer, optional,  intent(in   ) :: level     ! vertical level
      integer, optional,  intent(in   ) :: boundary  ! boundary edge
      real (r8),          intent(in   ) :: deps      ! convergence norm
      character (len=*),  intent(in   ) :: label     ! solver label
!
      integer, parameter                :: iunit = 45
      integer                           :: status
!
!  Report CG or CI algorithm convergence, normalized ratio of initial
!  and iterative inner products.
!
      status=0
!
      IF (Master) THEN
        IF (PRESENT(level)) THEN
          IF ((iterDiff.eq.1).and.(level.eq.1).and.(iter.eq.0)) THEN
            IF (PRESENT(boundary)) THEN
              WRITE (iunit,10) label, TRIM(Vname(1,idSvar(ifield))),    &
     &                         ms, boundary
            ELSE
              WRITE (iunit,20) label, TRIM(Vname(1,idSvar(ifield))), ms
            END IF
          END IF
        ELSE
          IF ((iterDiff.eq.1).and.(iter.eq.0)) THEN
            IF (PRESENT(boundary)) THEN
              WRITE (iunit,10) label, TRIM(Vname(1,idSvar(ifield))),    &
     &                         ms, boundary
            ELSE
              WRITE (iunit,20) label, TRIM(Vname(1,idSvar(ifield))), ms
            END IF
          END IF
        END IF
        IF (iter.eq.Niter) THEN
          IF (PRESENT(level)) THEN
            WRITE (iunit,30) iterDiff, iter, deps, level
          ELSE
            WRITE (iunit,40) iterDiff, iter, deps
          END IF
        END IF
        IF (iterDiff.eq.Ndiff) THEN
          IF (PRESENT(level)) THEN
            WRITE (iunit,50) iterDiff, iter, deps, level
          ELSE
            WRITE (iunit,60) iterDiff, iter, deps
          END IF
          FLUSH (iunit)
        END IF
      END IF
!
  10  FORMAT (a,': Convergence for ''',a,'''', ', B-MultiScale = ',i0,  &
     &        ', boundary = ',i0,/)
  20  FORMAT (a,': Convergence for ''',a,'''', ', B-MultiScale = ',i0,/)
  30  FORMAT (7x,'K-LAP  => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5,2x,'level = ',i0) 
  40  FORMAT (7x,'K-LAP  => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5) 
  50  FORMAT (7x,'K-DIFF => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5,2x,'level = ',i0) 
  60  FORMAT (7x,'K-DIFF => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5) 
!
      RETURN
      END FUNCTION multiscale_print
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  It computes the dot product between two 1D tiled boundary vectors.
!  
      FUNCTION dot_prod_1d (ng, tile, model, ctype, ibry,               &
     &                      LBij, UBij, B1, B2) RESULT (DotProd)
!
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: LBij, UBij
      real (r8),          intent(in   ) :: B1(LBij:)
      real (r8),          intent(in   ) :: B2(LBij:)
!
      logical, dimension(4)             :: Lboundary
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: NSUB, i, j
      real (r8)                         :: DotProd
      real (r8)                         :: cff, my_DotProd
#  ifdef DISTRIBUTE
      character (len=3)                 :: op_handle
#  endif
!
!  Initialize.
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
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
        CASE (u2dvar, u3dvar)
          Imin=IstrU
        CASE (v2dvar, v3dvar)
          Jmin=JstrV
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              cff=B1(j)*B2(j)
              my_DotProd=my_DotProd+cff
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              cff=B1(i)*B2(i)
              my_DotProd=my_DotProd+cff
            END DO
        END SELECT
      END IF
!
!  Perform parallel global reduction operation.
!
#  ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#  else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#  endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
#  endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END FUNCTION dot_prod_1d
# endif
!
!-----------------------------------------------------------------------
!  It compute the dot product between two 2D tiled arrays.
!  
      FUNCTION dot_prod_2d (ng, tile, model, ctype,                     &
     &                      LBi, UBi, LBj, UBj, A1, A2) RESULT (DotProd)
!
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      real (r8),          intent(in   ) :: A1(LBi:,LBj:)
      real (r8),          intent(in   ) :: A2(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: NSUB, i, j
      real (r8)                         :: DotProd
      real (r8)                         :: cff, my_DotProd
# ifdef DISTRIBUTE
      character (len=3)                 :: op_handle
# endif
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
        CASE (u2dvar, u3dvar)
          Imin=IstrU
        CASE (v2dvar, v3dvar)
          Jmin=JstrV
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
      END FUNCTION dot_prod_2d
#endif /* MULTI_SCALE_B */
      END MODULE roms_multiscale_mod
