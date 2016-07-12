#include "cppdefs.h"
      MODULE inner2state_mod
#if defined HESSIAN_SV || defined HESSIAN_SO || defined HESSIAN_FSV
!
!svn $Id: inner2state.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the tangent linear model  initial conditions  !
!  as the weighted sum of all Lanczos vectors computed from the first  !
!  outer loop of the I4DVAR Lanczos algorithm. It converts from 4DVAR  !
!  inner loop space span by assimilation to physical state space.      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
!
      PUBLIC  :: ad_inner2state
      PUBLIC  :: tl_inner2state
      PUBLIC  :: ini_C_norm
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_inner2state (ng, tile, Lini, state)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lini

      real(r8), intent(in) :: state(Ninner)
!
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
# endif
      CALL tl_inner2state_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lini, state,                            &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % tl_t_obc,                &
     &                          BOUNDARY(ng) % tl_u_obc,                &
     &                          BOUNDARY(ng) % tl_v_obc,                &
#   ifdef LCZ_FINAL
     &                          BOUNDARY(ng) % ad_t_obc,                &
     &                          BOUNDARY(ng) % ad_u_obc,                &
     &                          BOUNDARY(ng) % ad_v_obc,                &
#   endif
#  endif
     &                          BOUNDARY(ng) % tl_ubar_obc,             &
     &                          BOUNDARY(ng) % tl_vbar_obc,             &
     &                          BOUNDARY(ng) % tl_zeta_obc,             &
#  ifdef LCZ_FINAL
     &                          BOUNDARY(ng) % ad_ubar_obc,             &
     &                          BOUNDARY(ng) % ad_vbar_obc,             &
     &                          BOUNDARY(ng) % ad_zeta_obc,             &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % tl_ustr,                   &
     &                          FORCES(ng) % tl_vstr,                   &
#  ifdef LCZ_FINAL
     &                          FORCES(ng) % ad_ustr,                   &
     &                          FORCES(ng) % ad_vstr,                   &
#  endif
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % tl_tflux,                  &
#  ifdef LCZ_FINAL
     &                          FORCES(ng) % ad_tflux,                  &
#  endif
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
#  endif
#  ifdef LCZ_FINAL
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#   endif
#  endif
# else
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
#  ifdef LCZ_FINAL
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#  endif
# endif
     &                          OCEAN(ng) % tl_zeta                     &
# ifdef LCZ_FINAL
     &                          ,OCEAN(ng) % ad_zeta                    &
# endif
# ifdef HESSIAN_FSV
#  ifdef SOLVE3D
     &                          ,GRID(ng) % Hz,                         &
     &                          OCEAN(ng) % f_t,                        &
     &                          OCEAN(ng) % f_u,                        &
     &                          OCEAN(ng) % f_v,                        &
     &                          OCEAN(ng) % f_ubar,                     &
     &                          OCEAN(ng) % f_vbar,                     &
#  else
     &                          OCEAN(ng) % f_ubar,                     &
     &                          OCEAN(ng) % f_vbar,                     &
#  endif
     &                          OCEAN(ng) % f_zeta                      &
# endif
     &                                             )
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
# endif

      RETURN
      END SUBROUTINE tl_inner2state
!
!***********************************************************************
      SUBROUTINE tl_inner2state_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, LBij, UBij,   &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Lini, state,                      &
# ifdef MASKING
     &                                rmask, umask, vmask,              &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                tl_t_obc, tl_u_obc, tl_v_obc,     &
#   ifdef LCZ_FINAL
     &                                ad_t_obc, ad_u_obc, ad_v_obc,     &
#   endif
#  endif
     &                                tl_ubar_obc, tl_vbar_obc,         &
     &                                tl_zeta_obc,                      &
#  ifdef LCZ_FINAL
     &                                ad_ubar_obc, ad_vbar_obc,         &
     &                                ad_zeta_obc,                      &
#  endif
# endif
# ifdef ADJUST_WSTRESS
     &                                tl_ustr, tl_vstr,                 &
#  ifdef LCZ_FINAL
     &                                ad_ustr, ad_vstr,                 &
#  endif
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                tl_tflux,                         &
#  ifdef LCZ_FINAL
     &                                ad_tflux,                         &
#  endif
# endif
# ifdef SOLVE3D
     &                                tl_t, tl_u, tl_v,                 &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                                tl_ubar, tl_vbar,                 &
#   endif
#  ifdef LCZ_FINAL
     &                                ad_t, ad_u, ad_v,                 &
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                                ad_ubar, ad_vbar,                 &
#   endif
#  endif
# else
     &                                tl_ubar, tl_vbar,                 &
#  ifdef LCZ_FINAL
     &                                ad_ubar, ad_vbar,                 &
#  endif
# endif
     &                                tl_zeta                           &
# ifdef LCZ_FINAL
     &                                ,ad_zeta                          &
# endif
# ifdef HESSIAN_FSV
#  ifdef SOLVE3D
     &                                , Hz, f_t, f_u, f_v,              &
     &                                f_ubar, f_vbar  ,                 &
#  else
     &                                f_ubar, f_vbar  ,                 &
#  endif
     &                                f_zeta                            &
# endif
     &                                       )
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
      USE state_copy_mod, ONLY : state_copy
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcastf
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lini
!
      real(r8), intent(in) :: state(Ninner)
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#    ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#    endif
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#   endif
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#    endif
#   endif
#   ifdef HESSIAN_FSV
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#   endif
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#   endif
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#   ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
#   endif
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#   endif
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#    ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#    endif
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#   endif
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#   ifdef HESSIAN_FSV
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#    if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#    endif
#   endif
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#   endif
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
#  endif
#   ifdef LCZ_FINAL
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   endif
# endif
!
!  Local variable declarations.
!
      integer :: Lwrk, i, j, lstr, outLoop, rec, info

      integer :: ndefLCZ = 1
# ifdef LCZ_FINAL
      integer :: ndefLZE = 1
      integer :: ncidsav
# endif
# ifdef SOLVE3D
      integer :: itrc, k
#  ifdef HESSIAN_FSV
      real(r8) :: cff1, cff2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
#  endif
# endif
      real(r8) :: fac, fac1, fac2
      real(r8) :: zbeta

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Ninner) :: DotProd
      real(r8), dimension(Ninner) :: bvector
      real(r8), dimension(Ninner) :: work
# ifdef LCZ_FINAL
      real(r8), dimension(Ninner,Ninner) :: GStemp
      real(r8), dimension(Ninner,Ninner) :: GSsub
      real(r8), dimension(Ninner) :: work1
# endif

# ifdef LCZ_FINAL
      real :: sum
      logical, save :: first = .TRUE.
      logical, save :: first1 = .TRUE.
# endif

      character (len=256) :: ncname

# include "set_bounds.h"
!
      SourceFile='tl_inner2state.F, tl_inner2state_tile'
!
!  Compute a Cholesky factorization of the tridiagonal matrix
!  of inner-loop Lanczos vector coefficients of the form
!  L*D*L' where L is a lower unit lower bidiagonal matrix,
!  and D is a diagonal matrix.
!
      outLoop=Nouter
      IF (Master) THEN
        DO i=1,Ninner
          zLanczos_diag(i)=cg_delta(i,outLoop)
        END DO
        DO i=2,Ninner
          zLanczos_offdiag(i-1)=cg_beta(i,outLoop)
        END DO
!
        CALL dpttrf (Ninner, zLanczos_diag, zLanczos_offdiag, info)
!
!  Overwrite zLanczos_diag with its square root.
!
        DO i=1,Ninner
          zLanczos_diag(i)=SQRT(zLanczos_diag(i))
        END DO
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, info)
      CALL mp_bcastf (ng, iTLM, zLanczos_diag)
      CALL mp_bcastf (ng, iTLM, zLanczos_offdiag)
# endif
      IF (info.ne.0) THEN
        IF (Master) PRINT *,' Error in DPTTRF: info = ', info
        exit_flag=8
        RETURN
      END IF
!
!  Compute inv(L')*inv(sqrt(D))*state.
!  Since L is lower unit bidiagonal, then L'=U, a unit
!  upper bidiagonal. The following loops solve for the inverse of U.
!
      DO i=1,Ninner
        work(i)=state(i)/zLanczos_diag(i)
      END DO
      DO i=Ninner-1,1,-1
        work(i)=work(i)-zLanczos_offdiag(i)*work(i+1)
      END DO
# ifdef LCZ_FINAL
!
!  If using the evolved Lanczos vectors, then we need to perform a
!  Gramm-Schmidt orthogonalization to compute a new set of orthonormal
!  basis functions since the evolved Lanczos vectors are no longer
!  orthonormal. The evolved Lanczos vectors from I4D-Var are computed
!  using the EVOLVED_LCZ cpp option in 4D-Var.
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
      CALL netcdf_get_ivar (ng, iTLM, TRIM(LCZ(ng)%name), 'ndefADJ',    &
     &                      ndefLCZ)
      IF (exit_flag.ne.NoError) RETURN
!
      IF (first) THEN
      first=.FALSE.
!
!  Determine Lanczos vector file to read.  The Lanczos vectors are
!  written into the adjoint NetCDF in the I4D-Var Lanczos algorithm.
!  The Lanczos vector for each inner loop is accumulated in the
!  unlimited dimension. The name of this file is provided here in
!  the LCZ(ng)%name variable since the ADM(ng)%name value will be
!  use in the adjoint sensitivity part.
!
      IF (ndefLCZ.gt.0) THEN
        lstr=LEN_TRIM(LCZ(ng)%name)
        WRITE (ncname,10) LCZ(ng)%name(1:lstr-8), Nouter
      ELSE
        ncname=LCZ(ng)%name
      END IF
      Lwrk=2
!
!   Initialize the Gramm-Schmidt matrix.
!
        DO j=1,Ninner
          DO i=1,Ninner
            GSmatrix(i,j)=0.0_r8
          END DO
        END DO
        DO j=1,Ninner
            GSmatrix(j,j)=1.0_r8
        END DO
!
      DO inner=1,Ninner
!
!   Initialize the Gramm-Schmidt sub-matrices.
!
        DO j=1,Ninner
          DO i=1,Ninner
            GStemp(i,j)=0.0_r8
          END DO
        END DO
        DO j=1,Ninner
            GStemp(j,j)=1.0_r8
        END DO
!
        ncname=LCZ(ng)%name
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, inner,                                   &
     &                   ndefLCZ, LCZ(ng)%ncid, TRIM(ncname),           &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   tl_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   tl_t, tl_u, tl_v,                              &
# else
     &                   tl_ubar, tl_vbar,                              &
# endif
     &                   tl_zeta)
!
!  Copy tl_var(Lwrk) into tl_var(Lini).
!
!
!  NOTE: In the case of WEAK_CONSTRAINT and TIME_CONV, tl_ubar, tl_vbar
!        ad_ubar and ad_vbar are only passed as required but are not
!        used in subsequent calculations.
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, Lini,                                    &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_t_obc,                            &
     &                   tl_u_obc, tl_u_obc,                            &
     &                   tl_v_obc, tl_v_obc,                            &
#  endif
     &                   tl_ubar_obc, tl_ubar_obc,                      &
     &                   tl_vbar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc, tl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_ustr,                              &
     &                   tl_vstr, tl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   tl_tflux, tl_tflux,                            &
#  endif
     &                   tl_t, tl_t,                                    &
     &                   tl_u, tl_u,                                    &
     &                   tl_v, tl_v,                                    &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                   tl_ubar, tl_ubar,                              &
     &                   tl_vbar, tl_vbar,                              &
#  endif
# else
     &                   tl_ubar, tl_ubar,                              &
     &                   tl_vbar, tl_vbar,                              &
# endif
     &                   tl_zeta, tl_zeta)
!
! Orthonormalize Lanczos vectors.
!
        DO rec=1,inner-1
!
!  Read in gradient just computed Hessian eigenvectors into tangent
!  linear state array index Lwrk.
!
        ncname=LZE(ng)%name
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefLZE, LZE(ng)%ncid, TRIM(ncname),           &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   ad_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product.
!
          CALL state_dotprod (ng, tile, iTLM,                           &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        NstateVar(ng), dot(0:),                   &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        tl_t_obc(:,:,:,:,Lini,:),                 &
     &                        ad_t_obc(:,:,:,:,Lwrk,:),                 &
     &                        tl_u_obc(:,:,:,:,Lini),                   &
     &                        ad_u_obc(:,:,:,:,Lwrk),                   &
     &                        tl_v_obc(:,:,:,:,Lini),                   &
     &                        ad_v_obc(:,:,:,:,Lwrk),                   &
#  endif
     &                        tl_ubar_obc(:,:,:,Lini),                  &
     &                        ad_ubar_obc(:,:,:,Lwrk),                  &
     &                        tl_vbar_obc(:,:,:,Lini),                  &
     &                        ad_vbar_obc(:,:,:,Lwrk),                  &
     &                        tl_zeta_obc(:,:,:,Lini),                  &
     &                        ad_zeta_obc(:,:,:,Lwrk),                  &
# endif
# ifdef ADJUST_WSTRESS
     &                        tl_ustr(:,:,:,Lini), ad_ustr(:,:,:,Lwrk), &
     &                        tl_vstr(:,:,:,Lini), ad_vstr(:,:,:,Lwrk), &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        tl_tflux(:,:,:,Lini,:),                   &
     &                        ad_tflux(:,:,:,Lwrk,:),                   &
#  endif
     &                        tl_t(:,:,:,Lini,:), ad_t(:,:,:,Lwrk,:),   &
     &                        tl_u(:,:,:,Lini), ad_u(:,:,:,Lwrk),       &
     &                        tl_v(:,:,:,Lini), ad_v(:,:,:,Lwrk),       &
# else
     &                        tl_ubar(:,:,Lini), ad_ubar(:,:,Lwrk),     &
     &                        tl_vbar(:,:,Lini), ad_vbar(:,:,Lwrk),     &
# endif
     &                        tl_zeta(:,:,Lini), ad_zeta(:,:,Lwrk))
!
          GStemp(rec,inner)=-dot(0)
!
!    tl_var(Lini) = fac1 * tl_var(Lini) + fac2 * ad_var(Lwrk)
!
          fac1=1.0_r8
          fac2=-dot(0)

!
!  NOTE: In the case of WEAK_CONSTRAINT and TIME_CONV, tl_ubar, tl_vbar
!        ad_ubar and ad_vbar are only passed as required but are not
!        used in subsequent calculations.
!
          CALL state_addition (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lini, Lwrk, Lini, fac1, fac2,            &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         tl_t_obc, ad_t_obc,                      &
     &                         tl_u_obc, ad_u_obc,                      &
     &                         tl_v_obc, ad_v_obc,                      &
#  endif
     &                         tl_ubar_obc, ad_ubar_obc,                &
     &                         tl_vbar_obc, ad_vbar_obc,                &
     &                         tl_zeta_obc, ad_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                         tl_ustr, ad_ustr,                        &
     &                         tl_vstr, ad_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         tl_tflux, ad_tflux,                      &
#  endif
     &                         tl_t, ad_t,                              &
     &                         tl_u, ad_u,                              &
     &                         tl_v, ad_v,                              &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                         tl_ubar, ad_ubar,                        &
     &                         tl_vbar, ad_vbar,                        &
#  endif
# else
     &                         tl_ubar, ad_ubar,                        &
     &                         tl_vbar, ad_vbar,                        &
# endif
     &                       tl_zeta, ad_zeta)
        END DO
!
!  Compute normalization factor.
!
        CALL state_dotprod (ng, tile, iTLM,                             &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      tl_t_obc(:,:,:,:,Lini,:),                   &
     &                      tl_t_obc(:,:,:,:,Lini,:),                   &
     &                      tl_u_obc(:,:,:,:,Lini),                     &
     &                      tl_u_obc(:,:,:,:,Lini),                     &
     &                      tl_v_obc(:,:,:,:,Lini),                     &
     &                      tl_v_obc(:,:,:,:,Lini),                     &
#  endif
     &                      tl_ubar_obc(:,:,:,Lini),                    &
     &                      tl_ubar_obc(:,:,:,Lini),                    &
     &                      tl_vbar_obc(:,:,:,Lini),                    &
     &                      tl_vbar_obc(:,:,:,Lini),                    &
     &                      tl_zeta_obc(:,:,:,Lini),                    &
     &                      tl_zeta_obc(:,:,:,Lini),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      tl_ustr(:,:,:,Lini), tl_ustr(:,:,:,Lini),   &
     &                      tl_vstr(:,:,:,Lini), tl_vstr(:,:,:,Lini),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      tl_tflux(:,:,:,Lini,:),                     &
     &                      tl_tflux(:,:,:,Lini,:),                     &
#  endif
     &                      tl_t(:,:,:,Lini,:), tl_t(:,:,:,Lini,:),     &
     &                      tl_u(:,:,:,Lini), tl_u(:,:,:,Lini),         &
     &                      tl_v(:,:,:,Lini), tl_v(:,:,:,Lini),         &
# else
     &                      tl_ubar(:,:,Lini), tl_ubar(:,:,Lini),       &
     &                      tl_vbar(:,:,Lini), tl_vbar(:,:,Lini),       &
# endif
     &                      tl_zeta(:,:,Lini), tl_zeta(:,:,Lini))
!
        DO i=1,inner
          GStemp(i,inner)=GStemp(i,inner)/SQRT(dot(0))
        END DO
!
!  Normalize the evolved Lanczos vectors:
!
!    tl_var(Lini) = fac * tl_var(Lini)
!
        fac=1.0_r8/SQRT(dot(0))

        CALL state_scale (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    Lini, Lini, fac,                              &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc, tl_u_obc, tl_v_obc,                 &
#  endif
     &                    tl_ubar_obc, tl_vbar_obc,                     &
     &                    tl_zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr, tl_vstr,                             &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                    tl_tflux,                                     &
#  endif
     &                    tl_t, tl_u, tl_v,                             &
# else
     &                    tl_ubar, tl_vbar,                             &
# endif
     &                    tl_zeta)
!
!  Copy tl_var(Lini) into  ad_var(Lini) and write into the Hessian
!  netcdf file.
!
!
!  NOTE: In the case of WEAK_CONSTRAINT and TIME_CONV, tl_ubar, tl_vbar
!        ad_ubar and ad_vbar are only passed as required but are not
!        used in subsequent calculations.
!
        CALL state_copy (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lini, Lini,                                    &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, tl_t_obc,                            &
     &                   ad_u_obc, tl_u_obc,                            &
     &                   ad_v_obc, tl_v_obc,                            &
#  endif
     &                   ad_ubar_obc, tl_ubar_obc,                      &
     &                   ad_vbar_obc, tl_vbar_obc,                      &
     &                   ad_zeta_obc, tl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, tl_ustr,                              &
     &                   ad_vstr, tl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux, tl_tflux,                            &
#  endif
     &                   ad_t, tl_t,                                    &
     &                   ad_u, tl_u,                                    &
     &                   ad_v, tl_v,                                    &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                   ad_ubar, tl_ubar,                              &
     &                   ad_vbar, tl_vbar,                              &
#  endif
# else
     &                   ad_ubar, tl_ubar,                              &
     &                   ad_vbar, tl_vbar,                              &
# endif
     &                   ad_zeta, tl_zeta)
!
! AMM: Had to add this temporary save of LZEncid otherwise the write
!      of ocean_time to the LZE netcdf file fails.
!
        IF (inner.gt.1) THEN
             LZE(ng)%ncid=ncidsav
        ELSE
             ncidsav=LZE(ng)%ncid
        END IF
        LZE(ng)%Rindex=inner-1
        CALL wrt_evolved (ng, Lini, Lini)
        IF (exit_flag.ne.NoERRor) RETURN
!
!   Build the full Gramm-Schmidt matrix from as the product of each
!   sub-matrix.
!
        DO j=1,Ninner
         DO i=1,Ninner
           sum=0.0_r8
           DO k=1,Ninner
            sum=sum+GSmatrix(i,k)*GStemp(k,j)
           END DO
           GSsub(i,j)=sum
         END DO
        END DO
        DO j=1,Ninner
         DO i=1,Ninner
          GSmatrix(i,j)=GSsub(i,j)
         END DO
        END DO

      END DO
!
!  Test the orthonormalization of the evolved Lanczos vectors.
!
!
      IF (Master) THEN
        print *,'Test of orthonormalization'
      END IF
!
      DO inner=1,Ninner
!
        ncname=LZE(ng)%name
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, inner,                                   &
     &                   ndefLZE, LZE(ng)%ncid, TRIM(ncname),           &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   tl_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   tl_t, tl_u, tl_v,                              &
# else
     &                   tl_ubar, tl_vbar,                              &
# endif
     &                   tl_zeta)
        DO rec=1,Ninner
!
        ncname=LZE(ng)%name
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, rec,                                     &
     &                   ndefLZE, LZE(ng)%ncid, TRIM(ncname),           &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   ad_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
          IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product.
!
          CALL state_dotprod (ng, tile, iTLM,                           &
     &                        LBi, UBi, LBj, UBj, LBij, UBij,           &
     &                        NstateVar(ng), dot(0:),                   &
# ifdef MASKING
     &                        rmask, umask, vmask,                      &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                        tl_t_obc(:,:,:,:,Lwrk,:),                 &
     &                        ad_t_obc(:,:,:,:,Lwrk,:),                 &
     &                        tl_u_obc(:,:,:,:,Lwrk),                   &
     &                        ad_u_obc(:,:,:,:,Lwrk),                   &
     &                        tl_v_obc(:,:,:,:,Lwrk),                   &
     &                        ad_v_obc(:,:,:,:,Lwrk),                   &
#  endif
     &                        tl_ubar_obc(:,:,:,Lwrk),                  &
     &                        ad_ubar_obc(:,:,:,Lwrk),                  &
     &                        tl_vbar_obc(:,:,:,Lwrk),                  &
     &                        ad_vbar_obc(:,:,:,Lwrk),                  &
     &                        tl_zeta_obc(:,:,:,Lwrk),                  &
     &                        ad_zeta_obc(:,:,:,Lwrk),                  &
# endif
# ifdef ADJUST_WSTRESS
     &                        tl_ustr(:,:,:,Lwrk), ad_ustr(:,:,:,Lwrk), &
     &                        tl_vstr(:,:,:,Lwrk), ad_vstr(:,:,:,Lwrk), &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                        tl_tflux(:,:,:,Lwrk,:),                   &
     &                        ad_tflux(:,:,:,Lwrk,:),                   &
#  endif
     &                        tl_t(:,:,:,Lwrk,:), ad_t(:,:,:,Lwrk,:),   &
     &                        tl_u(:,:,:,Lwrk), ad_u(:,:,:,Lwrk),       &
     &                        tl_v(:,:,:,Lwrk), ad_v(:,:,:,Lwrk),       &
# else
     &                        tl_ubar(:,:,Lwrk), ad_ubar(:,:,Lwrk),     &
     &                        tl_vbar(:,:,Lwrk), ad_vbar(:,:,Lwrk),     &
# endif
     &                        tl_zeta(:,:,Lwrk), ad_zeta(:,:,Lwrk))
!
         IF (Master) THEN
           print *,'inner=',inner,' rec=',rec,' dot-product=',dot(0)
         END IF
!
        END DO
      END DO
!
!  Compute the inverse of GSmatrix.
!
      IF (Master) THEN
        DO j=1,Ninner
         DO i=1,Ninner
          GSmatinv(i,j)=GSmatrix(i,j)
         END DO
        END DO
        CALL dtrtri('U','N',Ninner,GSmatinv,Ninner,info)
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, info)
      CALL mp_bcastf (ng, iTLM, GSmatrix)
      CALL mp_bcastf (ng, iTLM, GSmatinv)
# endif
      IF (info.ne.0) THEN
        IF (Master) PRINT *,' Error in DPTTRF: info = ', info
        exit_flag=8
        RETURN
      END IF
!
      END IF
!
!  Now multiply by GSmatinv.
!
      DO i=1,Ninner
        sum=0.0_r8
        DO j=1,Ninner
         sum=sum+GSmatinv(i,j)*work(j)
        END DO
        work1(i)=sum
      END DO
      DO i=1,Ninner
        work(i)=work1(i)
      END DO

# endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear model initial conditions from the weighted
!  sum of the Lanczos vectors.
!-----------------------------------------------------------------------
!
# ifdef LCZ_FINAL
!  Determine if single or multiple Lanczos vector NetCDF files.
!
      CALL netcdf_get_ivar (ng, iTLM, TRIM(LZE(ng)%name), 'ndefADJ',    &
     &                      ndefLZE)
      IF (exit_flag.ne.NoError) RETURN
# else
!  Determine if single or multiple Lanczos vector NetCDF files.
!
      CALL netcdf_get_ivar (ng, iTLM, TRIM(LCZ(ng)%name), 'ndefADJ',    &
     &                      ndefLCZ)
      IF (exit_flag.ne.NoError) RETURN
# endif
!
# ifdef LCZ_FINAL
!  Determine Lanczos vector file to read. The reorthonormalized
!  Lanczos vectors are written in the Hessian netcdf file.
!
        IF (ndefLZE.gt.0) THEN
          lstr=LEN_TRIM(LZE(ng)%name)
          WRITE (ncname,10) LZE(ng)%name(1:lstr-8), Nouter
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=LZE(ng)%name
        END IF
# else
!  Determine Lanczos vector file to read.  The Lanczos vectors are
!  written into the adjoint NetCDF in the I4D-Var Lanczos algorithm.
!  The Lanczos vector for each inner loop is accumulated in the
!  unlimited dimension. The name of this file is provided here in
!  the LCZ(ng)%name variable since the ADM(ng)%name value will be
!  use in the adjoint sensitivity part.
!
        IF (ndefLCZ.gt.0) THEN
          lstr=LEN_TRIM(LCZ(ng)%name)
          WRITE (ncname,10) LCZ(ng)%name(1:lstr-8), Nouter
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=LCZ(ng)%name
        END IF
# endif
!
      Lwrk=2
# ifdef LCZ_FINAL
!
! AMM: Had to add this temporary save of LZEncid otherwise the write
!      of ocean_time to the LZE netcdf file fails.
!
      IF (first1) THEN
        first1=.FALSE.
        LZE(ng)%ncid=ncidsav
      END IF
!
# endif
      DO inner=1,Ninner
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the I4D-Var algorithm first outer loop. Load
!  Lanczos vectors into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, inner,                                   &
# ifdef LCZ_FINAL
     &                   ndefLZE, LZE(ng)%ncid, TRIM(ncname),           &
# else
     &                   ndefLCZ, LCZ(ng)%ncid, TRIM(ncname),           &
# endif
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   tl_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   tl_t, tl_u, tl_v,                              &
# else
     &                   tl_ubar, tl_vbar,                              &
# endif
     &                   tl_zeta)
        IF (exit_flag.ne.NoError) RETURN

        IF (inner.eq.1) THEN
          fac=0.0_r8
          CALL state_initialize (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           Lini, fac,                             &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           tl_t_obc, tl_u_obc, tl_v_obc,          &
#  endif
     &                           tl_ubar_obc, tl_vbar_obc,              &
     &                           tl_zeta_obc,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                           tl_tflux,                              &
# endif
# ifdef SOLVE3D
     &                           tl_t, tl_u, tl_v,                      &
# else
     &                           tl_ubar, tl_vbar,                      &
# endif
     &                           tl_zeta)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!
!  NOTE: In the case of WEAK_CONSTRAINT and TIME_CONV, tl_ubar, tl_vbar
!        ad_ubar and ad_vbar are only passed as required but are not
!        used in subsequent calculations.
!
        fac1=1.0_r8
        fac2=work(inner)
        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lini, Lwrk, Lini,                          &
     &                       fac1, fac2,                                &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       tl_t_obc, tl_t_obc,                        &
     &                       tl_u_obc, tl_u_obc,                        &
    &                        tl_v_obc, tl_v_obc,                        &
#  endif
     &                       tl_ubar_obc, tl_ubar_obc,                  &
     &                       tl_vbar_obc, tl_vbar_obc,                  &
     &                       tl_zeta_obc, tl_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, tl_ustr,                          &
    &                        tl_vstr, tl_vstr,                          &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       tl_tflux, tl_tflux,                        &
# endif
# ifdef SOLVE3D
     &                       tl_t, tl_t,                                &
     &                       tl_u, tl_u,                                &
     &                       tl_v, tl_v,                                &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       tl_ubar, tl_ubar,                          &
     &                       tl_vbar, tl_vbar,                          &
#  endif
# else
     &                       tl_ubar, tl_ubar,                          &
     &                       tl_vbar, tl_vbar,                          &
# endif
     &                       tl_zeta, tl_zeta)
        IF (exit_flag.ne.NoError) RETURN

      END DO
# ifdef HESSIAN_FSV
!
! Copy tl_var into f_var arrays then clear the tl_var arrays.
!
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f_zeta(i,j)=tl_zeta(i,j,Lini)
        END DO
      END DO
#  ifndef SOLVE3D
!
!  Tangent linear 2D momentum.

      DO j=JstrR,JendR
        DO i=Istr,IendR
          f_ubar(i,j)=tl_ubar(i,j,Lini)
       END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          f_vbar(i,j)=tl_vbar(i,j,Lini)
        END DO
      END DO

#  else
!
!  Tangent linear 3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            f_u(i,j,k)=tl_u(i,j,k,Lini)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            f_v(i,j,k)=tl_v(i,j,k,Lini)
          END DO
        END DO
      END DO
!
!   Compute the forcing for tl_ubar based on f_u.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*f_u(i,j,k)
          END DO
        END DO
        DO i=IstrU,Iend
          cff1=1.0_r8/DC(i,0)
          cff2=CF(i,0)*cff1
#  ifdef MASKING
          cff2=cff2*umask(i,j)
#  endif
          f_ubar(i,j)=cff2
        END DO
      END DO
!
!   Compute the forcing for tl_vbar based on f_v.
!
      DO j=JstrV,Jend
        IF (j.ge.JstrM) THEN
          DO i=Istr,Iend
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*f_v(i,j,k)
            END DO
          END DO
          DO i=Istr,Iend
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
#  ifdef MASKING
            cff2=cff2*vmask(i,j)
#  endif
            f_vbar(i,j)=cff2
          END DO
        END IF
      END DO
!
!  Tangent linear tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              f_t(i,j,k,itrc)=tl_t(i,j,k,Lini,itrc)
            END DO
          END DO
        END DO
      END DO
#  endif
      fac=0.0_r8
       CALL state_initialize (ng, tile,                                 &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           Lini, fac,                             &
#  ifdef MASKING
     &                           rmask, umask, vmask,                   &
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
     &                           tl_t_obc, tl_u_obc, tl_v_obc,          &
#   endif
     &                           tl_ubar_obc, tl_vbar_obc,              &
     &                           tl_zeta_obc,                           &
#  endif
#  ifdef ADJUST_WSTRESS
     &                           tl_ustr, tl_vstr,                      &
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
     &                           tl_tflux,                              &
#  endif
#  ifdef SOLVE3D
     &                           tl_t, tl_u, tl_v,                      &
#  else
     &                           tl_ubar, tl_vbar,                      &
#  endif
     &                           tl_zeta)

# endif

      RETURN
      END SUBROUTINE tl_inner2state_tile
!
!***********************************************************************
      SUBROUTINE ad_inner2state (ng, tile, Lini, ad_state)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lini
!
      real(r8), intent(inout) :: ad_state(Ninner)
!
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
# endif
      CALL ad_inner2state_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          Lini, ad_state,                         &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % ad_t_obc,                &
     &                          BOUNDARY(ng) % ad_u_obc,                &
     &                          BOUNDARY(ng) % ad_v_obc,                &
#  endif
     &                          BOUNDARY(ng) % ad_ubar_obc,             &
     &                          BOUNDARY(ng) % ad_vbar_obc,             &
     &                          BOUNDARY(ng) % ad_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % ad_ustr,                   &
     &                          FORCES(ng) % ad_vstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % ad_tflux,                  &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#  endif
# else
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
# endif
     &                          OCEAN(ng) % ad_zeta,                    &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                          BOUNDARY(ng) % tl_t_obc,                &
     &                          BOUNDARY(ng) % tl_u_obc,                &
     &                          BOUNDARY(ng) % tl_v_obc,                &
#  endif
     &                          BOUNDARY(ng) % tl_ubar_obc,             &
     &                          BOUNDARY(ng) % tl_vbar_obc,             &
     &                          BOUNDARY(ng) % tl_zeta_obc,             &
# endif
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % tl_ustr,                   &
     &                          FORCES(ng) % tl_vstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % tl_tflux,                  &
# endif
# ifdef SOLVE3D
     &                          OCEAN(ng) % tl_t,                       &
     &                          OCEAN(ng) % tl_u,                       &
     &                          OCEAN(ng) % tl_v,                       &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
#  endif
# else
     &                          OCEAN(ng) % tl_ubar,                    &
     &                          OCEAN(ng) % tl_vbar,                    &
# endif
     &                          OCEAN(ng) % tl_zeta                     &
# ifdef HESSIAN_FSV
#  ifdef SOLVE3D
     &                          ,GRID(ng) % Hz,                         &
     &                          OCEAN(ng) % f_t,                        &
     &                          OCEAN(ng) % f_u,                        &
     &                          OCEAN(ng) % f_v,                        &
     &                          OCEAN(ng) % f_ubar,                     &
     &                          OCEAN(ng) % f_vbar,                     &
#  else
     &                          OCEAN(ng) % f_ubar,                     &
     &                          OCEAN(ng) % f_vbar,                     &
#  endif
     &                          OCEAN(ng) % f_zeta                      &
# endif
     &                                  )
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
# endif

      RETURN
      END SUBROUTINE ad_inner2state
!
!***********************************************************************
      SUBROUTINE ad_inner2state_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, LBij, UBij,   &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                Lini, ad_state,                   &
# ifdef MASKING
     &                                rmask, umask, vmask,              &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                ad_t_obc, ad_u_obc, ad_v_obc,     &
#  endif
     &                                ad_ubar_obc, ad_vbar_obc,         &
     &                                ad_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                                ad_ustr, ad_vstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                ad_tflux,                         &
# endif
# ifdef SOLVE3D
     &                                ad_t, ad_u, ad_v,                 &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                                ad_ubar, ad_vbar,                 &
#  endif
# else
     &                                ad_ubar, ad_vbar,                 &
# endif
     &                                ad_zeta,                          &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                                tl_t_obc, tl_u_obc, tl_v_obc,     &
#  endif
     &                                tl_ubar_obc, tl_vbar_obc,         &
     &                                tl_zeta_obc,                      &
# endif
# ifdef ADJUST_WSTRESS
     &                                tl_ustr, tl_vstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                tl_tflux,                         &
# endif
# ifdef SOLVE3D
     &                                tl_t, tl_u, tl_v,                 &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                                tl_ubar, tl_vbar,                 &
#  endif
# else
     &                                tl_ubar, tl_vbar,                 &
# endif
     &                                tl_zeta                           &
# ifdef HESSIAN_FSV
#  ifdef SOLVE3D
     &                                ,Hz, f_t, f_u, f_v,               &
     &                                f_ubar, f_vbar  ,                 &
#  else
     &                                f_ubar, f_vbar  ,                 &
#  endif
     &                                f_zeta                            &
# endif
     &                                       )
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lini
!
      real(r8), intent(inout) :: ad_state(Ninner)
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
#   ifdef HESSIAN_FSV
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#   endif
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#   ifdef HESSIAN_FSV
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#   endif
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef HESSIAN_FSV
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   if defined WEAK_CONSTRAINT && defined TIME_CONV
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
# endif
!
      real(r8) :: work(Ninner)
      real(r8) :: dot_sav(Ninner)
# ifdef LCZ_FINAL
      real(r8) :: work1(Ninner)
      real(r8) :: sum
# endif
!
!  Local variable declarations.
!
      integer :: Lwrk, i, j, lstr, outLoop, rec, info

      integer :: ndefLCZ = 1
# ifdef LCZ_FINAL
      integer :: ndefLZE = 1
# endif
      integer :: kin
# ifdef SOLVE3D
      integer :: itrc, k, nin
#  ifdef HESSIAN_FSV
      real(r8) :: cff1, cff2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
#  endif
# endif
      real(r8) :: fac, fac1, fac2
      real(r8) :: zbeta

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Ninner) :: DotProd
      real(r8), dimension(Ninner) :: bvector

      character (len=256) :: ncname

# include "set_bounds.h"
!
      SourceFile='ad_inner2state.F, ad_inner2state_tile'
!
!  Clear ad_state array.
!
      DO i=1,Ninner
        ad_state(i)=0.0_r8
      END DO
!
# ifdef SOLVE3D
      nin=nstp(ng)
# endif
# ifdef HESSIAN_SO
      kin=knew(ng)
# else
      kin=kstp(ng)
# endif
# ifdef HESSIAN_FSV
!
! Copy f_var into ad_var arrays then clear the f_var arrays.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          ad_zeta(i,j,kin)=f_zeta(i,j)
          f_zeta(i,j)=0.0_r8
        END DO
      END DO
#  ifndef SOLVE3D
!
!  Tangent linear 2D momentum.

      DO j=JstrR,JendR
        DO i=Istr,IendR
          ad_ubar(i,j,kin)=f_ubar(i,j)
          f_ubar(i,j)=0.0_r8
       END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          ad_vbar(i,j,kin)=f_vbar(i,j)
          f_vbar(i,j)=0.0_r8
        END DO
      END DO

#  else
!
!   Compute the contribution of f_ubar to f_u.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          DC(i,0)=0.0_r8
          CF(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
            DC(i,0)=DC(i,0)+DC(i,k)
          END DO
        END DO
        DO i=IstrU,Iend
          cff2=f_ubar(i,j)
          f_ubar(i,j)=0.0_r8
#  ifdef MASKING
          cff2=cff2*umask(i,j)
#  endif
          cff1=1.0_r8/DC(i,0)
          CF(i,0)=cff2*cff1
          cff2=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            f_u(i,j,k)=f_u(i,j,k)+DC(i,k)*CF(i,0)
          END DO
        END DO
        DO i=IstrU,Iend
          CF(i,0)=0.0_r8
        END DO
      END DO
!
!   Compute the contribution of f_vbar to f_v.
!
      DO j=JstrV,Jend
        IF (j.ge.JstrM) THEN
          DO i=Istr,Iend
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
              DC(i,0)=DC(i,0)+DC(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            cff2=f_vbar(i,j)
            f_vbar(i,j)=0.0_r8
#  ifdef MASKING
            cff2=cff2*vmask(i,j)
#  endif
            cff1=1.0_r8/DC(i,0)
            CF(i,0)=cff2*cff1
            cff2=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              f_v(i,j,k)=f_v(i,j,k)+DC(i,k)*CF(i,0)
            END DO
          END DO
          DO i=Istr,Iend
            CF(i,0)=0.0_r8
          END DO
        END IF
      END DO
!
!  Tangent linear 3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ad_u(i,j,k,nin)=f_u(i,j,k)
            f_u(i,j,k)=0.0_r8
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            ad_v(i,j,k,nin)=f_v(i,j,k)
            f_v(i,j,k)=0.0_r8
          END DO
        END DO
      END DO
!
!  Tangent linear tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              ad_t(i,j,k,nin,itrc)=f_t(i,j,k,itrc)
              f_t(i,j,k,itrc)=0.0_r8
            END DO
          END DO
        END DO
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Compute tangent linear model initial conditions from the weighted
!  sum of the Lanczos vectors.
!-----------------------------------------------------------------------
!
!  Determine if single or multiple Lanczos vector NetCDF files.
# ifdef LCZ_FINAL
!
      CALL netcdf_get_ivar (ng, iADM, TRIM(LZE(ng)%name), 'ndefADJ',    &
     &                      ndefLZE)
      IF (exit_flag.ne.NoError) RETURN
# else
!
      CALL netcdf_get_ivar (ng, iADM, TRIM(LCZ(ng)%name), 'ndefADJ',    &
     &                      ndefLCZ)
      IF (exit_flag.ne.NoError) RETURN
# endif
!
      Lwrk=2
      DO inner=1,Ninner
# ifdef LCZ_FINAL
!
!  Determine Lanczos vector file to read.  The orthonormalized evolved
!  Lanczos vectors are written into the Hessian.
!
        IF (ndefLZE.gt.0) THEN
          lstr=LEN_TRIM(LZE(ng)%name)
          WRITE (ncname,10) LZE(ng)%name(1:lstr-8), inner
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=LZE(ng)%name
        END IF
# else
!
!  Determine Lanczos vector file to read.  The Lanczos vectors are
!  written into the adjoint NetCDF in the I4D-Var Lanczos algorithm.
!  The Lanczos vector for each inner loop is accumulated in the
!  unlimited dimension. The name of this file is provided here in
!  the LCZ(ng)%name variable since the ADM(ng)%name value will be
!  use in the adjoint sensitivity part.
!
        IF (ndefLCZ.gt.0) THEN
          lstr=LEN_TRIM(LCZ(ng)%name)
          WRITE (ncname,10) LCZ(ng)%name(1:lstr-8), inner
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=LCZ(ng)%name
        END IF
# endif
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the I4D-Var algorithm first outer loop. Load
!  Lanczos vectors into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, iADM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, inner,                                   &
# ifdef LCZ_FINAL
     &                   ndefLZE, LZE(ng)%ncid, TRIM(ncname),           &
# else
     &                   ndefLCZ, LCZ(ng)%ncid, TRIM(ncname),           &
# endif
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   tl_t_obc, tl_u_obc, tl_v_obc,                  &
#  endif
     &                   tl_ubar_obc, tl_vbar_obc,                      &
     &                   tl_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   tl_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   tl_t, tl_u, tl_v,                              &
# else
     &                   tl_ubar, tl_vbar,                              &
# endif
     &                   tl_zeta)
        IF (exit_flag.ne.NoError) RETURN

        CALL state_dotprod (ng, tile, iADM,                             &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      ad_t_obc(:,:,:,:,Lnew,:),                   &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      ad_u_obc(:,:,:,:,Lnew),                     &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      ad_v_obc(:,:,:,:,Lnew),                     &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#  endif
     &                      ad_ubar_obc(:,:,:,Lnew),                    &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      ad_vbar_obc(:,:,:,Lnew),                    &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      ad_zeta_obc(:,:,:,Lnew),                    &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Lnew), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Lnew), tl_vstr(:,:,:,Lwrk),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Lnew,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
#  endif
     &                      ad_t(:,:,:,nin,:), tl_t(:,:,:,Lwrk,:),      &
     &                      ad_u(:,:,:,nin), tl_u(:,:,:,Lwrk),          &
     &                      ad_v(:,:,:,nin), tl_v(:,:,:,Lwrk),          &
# else
     &                      ad_ubar(:,:,kin), tl_ubar(:,:,Lwrk),        &
     &                      ad_vbar(:,:,kin), tl_vbar(:,:,Lwrk),        &
# endif
     &                      ad_zeta(:,:,kin), tl_zeta(:,:,Lwrk))

        dot_sav(inner)=dot(0)
      END DO
      DO i=1,Ninner
        work(i)=dot_sav(i)
      END DO
# ifdef LCZ_FINAL
!
!  Now multiply by Tranpose(GSmatinv).
!
      DO i=1,Ninner
        sum=0.0_r8
        DO j=1,Ninner
         sum=sum+GSmatinv(j,i)*work(j)
        END DO
        work1(i)=sum
      END DO
      DO i=1,Ninner
        work(i)=work1(i)
      END DO
# endif
!
!  Compute inv(sqrt(D))*inv(L).
!
      DO i=1,Ninner-1
       work(i+1)=work(i+1)-zLanczos_offdiag(i)*work(i)
      END DO
      DO i=1,Ninner
       ad_state(i)=work(i)/zLanczos_diag(i)
       work(i)=0.0_r8
      END DO
!
      RETURN
      END SUBROUTINE ad_inner2state_tile
!
!***********************************************************************
      SUBROUTINE ini_C_norm (ng, tile, Kinp, Ninp, StateNorm)
!***********************************************************************
!
      USE mod_param
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kinp, Ninp

      real(r8), intent(out) :: StateNorm
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
# endif
      CALL ini_C_norm_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Kinp, Ninp,                                 &
     &                      StateNorm,                                  &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      BOUNDARY(ng) % ad_t_obc,                    &
     &                      BOUNDARY(ng) % ad_u_obc,                    &
     &                      BOUNDARY(ng) % ad_v_obc,                    &
#  endif
     &                      BOUNDARY(ng) % ad_ubar_obc,                 &
     &                      BOUNDARY(ng) % ad_vbar_obc,                 &
     &                      BOUNDARY(ng) % ad_zeta_obc,                 &
# endif
# ifdef ADJUST_WSTRESS
     &                      FORCES(ng) % ad_ustr,                       &
     &                      FORCES(ng) % ad_vstr,                       &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      FORCES(ng) % ad_tflux,                      &
# endif
# ifdef SOLVE3D
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
# else
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
# endif
     &                      OCEAN(ng) % ad_zeta,                        &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      BOUNDARY(ng) % tl_t_obc,                    &
     &                      BOUNDARY(ng) % tl_u_obc,                    &
     &                      BOUNDARY(ng) % tl_v_obc,                    &
#  endif
     &                      BOUNDARY(ng) % tl_ubar_obc,                 &
     &                      BOUNDARY(ng) % tl_vbar_obc,                 &
     &                      BOUNDARY(ng) % tl_zeta_obc,                 &
# endif
# ifdef ADJUST_WSTRESS
     &                      FORCES(ng) % tl_ustr,                       &
     &                      FORCES(ng) % tl_vstr,                       &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      FORCES(ng) % tl_tflux,                      &
# endif
# ifdef SOLVE3D
     &                      OCEAN(ng) % tl_t,                           &
     &                      OCEAN(ng) % tl_u,                           &
     &                      OCEAN(ng) % tl_v,                           &
# else
     &                      OCEAN(ng) % tl_ubar,                        &
     &                      OCEAN(ng) % tl_vbar,                        &
# endif
     &                      OCEAN(ng) % tl_zeta,                        &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      BOUNDARY(ng) % t_obc,                       &
     &                      BOUNDARY(ng) % u_obc,                       &
     &                      BOUNDARY(ng) % v_obc,                       &
#  endif
     &                      BOUNDARY(ng) % ubar_obc,                    &
     &                      BOUNDARY(ng) % vbar_obc,                    &
     &                      BOUNDARY(ng) % zeta_obc,                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      FORCES(ng) % ustr,                          &
     &                      FORCES(ng) % vstr,                          &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      FORCES(ng) % tflux,                         &
# endif
# ifdef SOLVE3D
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % u,                              &
     &                      OCEAN(ng) % v,                              &
# else
     &                      OCEAN(ng) % ubar,                           &
     &                      OCEAN(ng) % vbar,                           &
# endif
     &                      OCEAN(ng) % zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
# endif

      RETURN
      END SUBROUTINE ini_C_norm
!
!***********************************************************************
      SUBROUTINE ini_C_norm_tile (ng, tile ,                            &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Kinp, Ninp,                           &
     &                            StateNorm,                            &
# ifdef MASKING
     &                            rmask, umask, vmask,                  &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                            ad_t_obc, ad_u_obc, ad_v_obc,         &
#  endif
     &                            ad_ubar_obc, ad_vbar_obc,             &
     &                            ad_zeta_obc,                          &
# endif
# ifdef ADJUST_WSTRESS
     &                            ad_ustr, ad_vstr,                     &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                            ad_tflux,                             &
# endif
# ifdef SOLVE3D
     &                            ad_t, ad_u, ad_v,                     &
# else
     &                            ad_ubar, ad_vbar,                     &
# endif
     &                            ad_zeta,                              &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                            tl_t_obc, tl_u_obc, tl_v_obc,         &
#  endif
     &                            tl_ubar_obc, tl_vbar_obc,             &
     &                            tl_zeta_obc,                          &
# endif
# ifdef ADJUST_WSTRESS
     &                            tl_ustr, tl_vstr,                     &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                            tl_tflux,                             &
# endif
# ifdef SOLVE3D
     &                            tl_t, tl_u, tl_v,                     &
# else
     &                            tl_ubar, tl_vbar,                     &
# endif
     &                            tl_zeta,                              &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                            t_obc, u_obc, v_obc,                  &
#  endif
     &                            ubar_obc,_vbar_obc,                   &
     &                            zeta_obc,                             &
# endif
# ifdef ADJUST_WSTRESS
     &                            ustr, vstr,                           &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                            tflux,                                &
# endif
# ifdef SOLVE3D
     &                            t, u, v,                              &
# else
     &                            ubar, vbar,                           &
# endif
     &                            zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE state_addition_mod, ONLY : state_addition
      USE state_dotprod_mod, ONLY : state_dotprod
      USE state_initialize_mod, ONLY : state_initialize
      USE state_scale_mod, ONLY : state_scale
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Ninp

      real(r8), intent(out) :: StateNorm
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: ad_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: ad_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: tl_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: tl_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: ad_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ad_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: ad_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: tl_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: tl_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: tl_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: t_obc(LBij:UBij,N(ng),4,            &
     &                                    Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Lwrk1, Lwrk2, i, j, lstr, outLoop, rec

      integer :: ndefLCZ = 1
# ifdef LCZ_FINAL
      integer :: ndefLZE = 1
# endif
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: fac, fac1, fac2
      real(r8) :: zbeta

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Ninner) :: DotProd
      real(r8), dimension(Ninner) :: bvector
# ifdef LCZ_FINAL
      real(r8) :: sum
      real(r8), dimension(Ninner) :: work1
# endif

      character (len=256) :: ncname

# include "set_bounds.h"
!
!
!-----------------------------------------------------------------------
!  Compute tangent linear model initial conditions from the weighted
!  sum of the Lanczos vectors.
!-----------------------------------------------------------------------
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
# ifdef LCZ_FINAL
      CALL netcdf_get_ivar (ng, iADM, TRIM(LZE(ng)%name), 'ndefADJ',    &
     &                      ndefLZE)
      IF (exit_flag.ne.NoError) RETURN
# else
      CALL netcdf_get_ivar (ng, iADM, TRIM(LCZ(ng)%name), 'ndefADJ',    &
     &                      ndefLCZ)
      IF (exit_flag.ne.NoError) RETURN
# endif
!
      Lwrk1=1
      Lwrk2=2
      DO inner=1,Ninner
# ifdef LCZ_FINAL
!
!  Determine Lanczos vector file to read.  The orthonormalized
!  evolved Lanczos vectors are written into the hessian NetCDF.
!
        IF (ndefLZE.gt.0) THEN
          lstr=LEN_TRIM(LZE(ng)%name)
          WRITE (ncname,10) LZE(ng)%name(1:lstr-8), outLoop
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=LZE(ng)%name
        END IF
# else
!
!  Determine Lanczos vector file to read.  The Lanczos vectors are
!  written into the adjoint NetCDF in the I4D-Var Lanczos algorithm.
!  The Lanczos vector for each inner loop is accumulated in the
!  unlimited dimension. The name of this file is provided here in
!  the LCZ(ng)%name variable since the ADM(ng)%name value will be
!  use in the adjoint sensitivity part.
!
        IF (ndefLCZ.gt.0) THEN
          lstr=LEN_TRIM(LCZ(ng)%name)
          WRITE (ncname,10) LCZ(ng)%name(1:lstr-8), outLoop
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=LCZ(ng)%name
        END IF
# endif
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the I4D-Var algorithm first outer loop. Load
!  Lanczos vectors into ADJOINT LINEAR STATE ARRAYS at index Lwrk1.
!
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk1, inner,                                  &
# ifdef LCZ_FINAL
     &                   ndefLZE, LZE(ng)%ncid, TRIM(ncname),           &
# else
     &                   ndefLCZ, LCZ(ng)%ncid, TRIM(ncname),           &
# endif
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   ad_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product between the tangent linear solution, x(0),
!  and Lanczos vectors, q_i.
!
!    DotProd(inner) = a_i = < x(0), q_i) >
!
        CALL state_dotprod (ng, tile, iTLM,                             &
     &                      LBi, UBi, LBj, UBj, LBij, UBij,             &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                      tl_t_obc(:,:,:,:,Ladj,:),                   &
     &                      ad_t_obc(:,:,:,:,Lwrk1,:),                  &
     &                      tl_u_obc(:,:,:,:,Ladj),                     &
     &                      ad_u_obc(:,:,:,:,Lwrk1),                    &
     &                      tl_v_obc(:,:,:,:,Ladj),                     &
     &                      ad_v_obc(:,:,:,:,Lwrk1),                    &
#  endif
     &                      tl_ubar_obc(:,:,:,Ladj),                    &
     &                      ad_ubar_obc(:,:,:,Lwrk1),                   &
     &                      tl_vbar_obc(:,:,:,Ladj),                    &
     &                      ad_vbar_obc(:,:,:,Lwrk1),                   &
     &                      tl_zeta_obc(:,:,:,Ladj),                    &
     &                      ad_zeta_obc(:,:,:,Lwrk1),                   &
# endif
# ifdef ADJUST_WSTRESS
     &                      tl_ustr(:,:,:,Ladj), ad_ustr(:,:,:,Lwrk1),  &
     &                      tl_vstr(:,:,:,Ladj), ad_vstr(:,:,:,Lwrk1),  &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      tl_tflux(:,:,:,Ladj,:),                     &
     &                      ad_tflux(:,:,:,Lwrk1,:),                    &
# endif
# ifdef SOLVE3D
     &                      tl_t(:,:,:,Ninp,:), ad_t(:,:,:,Lwrk1,:),    &
     &                      tl_u(:,:,:,Ninp), ad_u(:,:,:,Lwrk1),        &
     &                      tl_v(:,:,:,Ninp), ad_v(:,:,:,Lwrk1),        &
# else
     &                      tl_ubar(:,:,Kinp), ad_ubar(:,:,Lwrk1),      &
     &                      tl_vbar(:,:,Kinp), ad_vbar(:,:,Lwrk1),      &
# endif
     &                      tl_zeta(:,:,Kinp), ad_zeta(:,:,Lwrk1))
!
!  Store dot product.
!
        DotProd(inner)=dot(0)
      END DO
# ifdef LCZ_FINAL
!
!  Now multiply by GSmatrix.
!
      DO i=1,Ninner
        sum=0.0_r8
        DO j=1,Ninner
         sum=sum+GSmatrix(i,j)*DotProd(j)
        END DO
        work1(i)=sum
      END DO
      DO i=1,Ninner
        DotProd(i)=work1(i)
      END DO
# endif
!
!  Multiply by the tridiagonal matrix associated with the
!  Lanczos recursion relation.
!
!  For now, we can only use the first outer loop. A different scaling
!  is required for additional outer loops.
!
      outLoop=1
!
      bvector(1)=cg_delta(1,outLoop)*DotProd(1)+                        &
     &           cg_beta (2,outLoop)*DotProd(2)
      DO i=2,Ninner-1
        bvector(i)=cg_delta(i,outLoop)*DotProd(i)+                      &
     &             cg_beta(i+1,outLoop)*DotProd(i+1)+                   &
     &             cg_beta(i,outLoop)*DotProd(i-1)
      END DO
      bvector(Ninner)=cg_delta(Ninner,outLoop)*DotProd(Ninner)+         &
     &                cg_beta (Ninner,outLoop)*DotProd(Ninner-1)
# ifdef LCZ_FINAL
!
!  Now multiply by Transpose(GSmatrix).
!
      DO i=1,Ninner
        sum=0.0_r8
        DO j=1,Ninner
         sum=sum+GSmatrix(j,i)*bvector(j)
        END DO
        work1(i)=sum
      END DO
      DO i=1,Ninner
        bvector(i)=work1(i)
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Compute Lanczos vectors weigthed sum.
!-----------------------------------------------------------------------
!
!  Initialize non-linear state arrays: var(Lwrk1) = fac
!
      fac=0.0_r8

      CALL state_initialize (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk1, fac,                                &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       t_obc, u_obc, v_obc,                       &
#  endif
     &                       ubar_obc, vbar_obc,                        &
     &                       zeta_obc,                                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       ustr, vstr,                                &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       tflux,                                     &
# endif
# ifdef SOLVE3D
     &                       t, u, v,                                   &
# else
     &                       ubar, vbar,                                &
# endif
     &                       zeta)
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the IS4DVAR algorithm first outer loop. Load
!  Lanczos vectors into ADJOINT STATE ARRAYS at index Lwrk1.
!
      DO inner=1,Ninner
# ifdef LCZ_FINAL
        IF (ndefLZE.gt.0) THEN
          lstr=LEN_TRIM(LZE(ng)%name)
          WRITE (ncname,10) LZE(ng)%name(1:lstr-8), inner
        ELSE
          ncname=LZE(ng)%name
        END IF
# else
        IF (ndefLCZ.gt.0) THEN
          lstr=LEN_TRIM(LCZ(ng)%name)
          WRITE (ncname,10) LCZ(ng)%name(1:lstr-8), inner
        ELSE
          ncname=LCZ(ng)%name
        END IF
# endif
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk1, inner,                                  &
# ifdef LCZ_FINAL
     &                   ndefLZE, LZE(ng)%ncid, ncname,                 &
# else
     &                   ndefLCZ, LCZ(ng)%ncid, ncname,                 &
# endif
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                   ad_t_obc, ad_u_obc, ad_v_obc,                  &
#  endif
     &                   ad_ubar_obc, ad_vbar_obc,                      &
     &                   ad_zeta_obc,                                   &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                   ad_tflux,                                      &
# endif
# ifdef SOLVE3D
     &                   ad_t, ad_u, ad_v,                              &
# else
     &                   ad_ubar, ad_vbar,                              &
# endif
     &                   ad_zeta)
        IF (exit_flag.ne.NoError) RETURN
!
!  Sum over all Lanczos vectors:
!
!    var(Lwrk2) = fac1 * var(Lwrk2) + fac2 * ad_var(Lwrk1)
!
!  This will become the tangent linear model initial conditions at
!  time index Lnew.
!
        fac1=1.0_r8
        fac2=bvector(inner)

!
!  NOTE: In the case of WEAK_CONSTRAINT and TIME_CONV, tl_ubar, tl_vbar
!        ad_ubar and ad_vbar are only passed as required but are not
!        used in subsequent calculations.
!
        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk2, Lwrk1, Lwrk2, fac1, fac2,           &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       t_obc, ad_t_obc,                           &
     &                       u_obc, ad_u_obc,                           &
     &                       v_obc, ad_v_obc,                           &
#  endif
     &                       ubar_obc, ad_ubar_obc,                     &
     &                       vbar_obc, ad_vbar_obc,                     &
     &                       zeta_obc, ad_zeta_obc,                     &
# endif
# ifdef ADJUST_WSTRESS
     &                       ustr, ad_ustr,                             &
     &                       vstr, ad_vstr,                             &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       tflux, ad_tflux,                           &
# endif
# ifdef SOLVE3D
     &                       t, ad_t,                                   &
     &                       u, ad_u,                                   &
     &                       v, ad_v,                                   &
#  if defined WEAK_CONSTRAINT && defined TIME_CONV
     &                       ubar, ad_ubar,                             &
     &                       vbar, ad_vbar,                             &
#  endif
# else
     &                       ubar, ad_ubar,                             &
     &                       vbar, ad_vbar,                             &
# endif
     &                       zeta, ad_zeta)
      END DO
!
! Finally compute the dot-product with the input tl vector.
!
      CALL state_dotprod (ng, tile, iTLM,                               &
     &                    LBi, UBi, LBj, UBj, LBij, UBij,               &
     &                    NstateVar(ng), dot(0:),                       &
# ifdef MASKING
     &                    rmask, umask, vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                    tl_t_obc(:,:,:,:,Ladj,:),                     &
     &                    t_obc(:,:,:,:,Lwrk2,:),                       &
     &                    tl_u_obc(:,:,:,:,Ladj),                       &
     &                    u_obc(:,:,:,:,Lwrk2),                         &
     &                    tl_v_obc(:,:,:,:,Ladj),                       &
     &                    v_obc(:,:,:,:,Lwrk2),                         &
#  endif
     &                    tl_ubar_obc(:,:,:,Ladj),                      &
     &                    ubar_obc(:,:,:,Lwrk2),                        &
     &                    tl_vbar_obc(:,:,:,Ladj),                      &
     &                    vbar_obc(:,:,:,Lwrk2),                        &
     &                    tl_zeta_obc(:,:,:,Ladj),                      &
     &                    zeta_obc(:,:,:,Lwrk2),                        &
# endif
# ifdef ADJUST_WSTRESS
     &                    tl_ustr(:,:,:,Ladj), ustr(:,:,:,Lwrk2),       &
     &                    tl_vstr(:,:,:,Ladj), vstr(:,:,:,Lwrk2),       &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                    tl_tflux(:,:,:,Ladj,:),                       &
     &                    tflux(:,:,:,Lwrk2,:),                         &
# endif
# ifdef SOLVE3D
     &                    tl_t(:,:,:,Ninp,:), t(:,:,:,Lwrk2,:),         &
     &                    tl_u(:,:,:,Ninp), u(:,:,:,Lwrk2),             &
     &                    tl_v(:,:,:,Ninp), v(:,:,:,Lwrk2),             &
# else
     &                    tl_ubar(:,:,Kinp), ubar(:,:,Lwrk2),           &
     &                    tl_vbar(:,:,Kinp), vbar(:,:,Lwrk2),           &
# endif
     &                    tl_zeta(:,:,Kinp), zeta(:,:,Lwrk2))
!
      StateNorm=dot(0)
!
      RETURN
      END SUBROUTINE ini_C_norm_tile
!
!
!***********************************************************************
      SUBROUTINE read_state (ng, tile, model,                           &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       s_t_obc, s_u_obc, s_v_obc,                 &
#  endif
     &                       s_ubar_obc, s_vbar_obc,                    &
     &                       s_zeta_obc,                                &
# endif
# ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       s_tflux,                                   &
# endif
# ifdef SOLVE3D
     &                       s_t, s_u, s_v,                             &
# else
     &                       s_ubar, s_vbar,                            &
# endif
     &                       s_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
#ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
# ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
# endif
#endif
      USE nf_fread2d_mod, ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod, ONLY : nf_fread3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lwrk, rec, ndef

      integer, intent(inout) :: ncfileid

      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: i, ifield, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif
      integer :: gtype, ncid, status, varid

      integer, dimension(4) :: Vsize
      integer, dimension(NV) :: vid

      real(r8) :: Fmin, Fmax, scale

# include "set_bounds.h"
!
      SourceFile='inner2state, read_state'
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF ((ndef.gt.0).or.(ncfileid.eq.-1)) THEN
        CALL netcdf_open (ng, model, ncname, 0, ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
        ncfileid=ncid
      ELSE
        ncid=ncfileid
      END IF

      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idFsur), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  rmask,                                          &
# endif
     &                  s_zeta(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        ifield=idSbry(isFsur)
        gtype=r2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_zeta_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifndef SOLVE3D
!
!  Read in 2D U-momentum component.
!
      gtype=u2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ubar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in 2D V-momentum component.
!
      gtype=v2dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVbar),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread2d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vbar(:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        ifield=idSbry(isUbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_ubar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Read in 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        ifield=idSbry(isVbar)
        gtype=u2dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread2d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_vbar_obc(:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ustr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVsms),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vstr(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
!
!  Read in 3D U-momentum component.
!
      gtype=u3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idUvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_u(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        ifield=idSbry(isUvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_u_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Read in 3D V-momentum component.
!
      gtype=v3dvar
      scale=1.0_r8
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),        &
     &                       ncid, varid)
      IF (exit_flag.ne.NoError) RETURN

      status=nf_fread3d(ng, model, ncname, ncid,                        &
     &                  Vname(1,idVvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_v(:,:,:,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        ifield=idSbry(isVvel)
        gtype=u3dvar
        scale=1.0_r8
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d_bry (ng, model, ncname, ncid,                 &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_v_obc(:,:,:,:,Lwrk))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Read in tracers.
!
      gtype=r3dvar
      scale=1.0_r8
      DO itrc=1,NT(ng)
        CALL netcdf_inq_varid (ng, model, ncname,                       &
     &                         Vname(1,idTvar(itrc)), ncid, varid)
        IF (exit_flag.ne.NoError) RETURN

        status=nf_fread3d(ng, model, ncname, ncid,                      &
     &                    Vname(1,idTvar(itrc)), varid,                 &
     &                    rec, gtype, Vsize,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask,                                        &
#  endif
     &                    s_t(:,:,:,Lwrk,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Read in tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          ifield=idSbry(isTvar(itrc))
          gtype=r3dvar
          scale=1.0_r8
          CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),    &
     &                           ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d_bry (ng, model, ncname, ncid,               &
     &                           Vname(1,ifield), varid,                &
     &                           rec, gtype,                            &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale, Fmin, Fmax,                     &
     &                           s_t_obc(:,:,:,:,Lwrk,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      scale=1.0_r8
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          CALL netcdf_inq_varid (ng, model, ncname,                     &
     &                           Vname(1,idTsur(itrc)), ncid, varid)
          IF (exit_flag.ne.NoError) RETURN

          status=nf_fread3d(ng, model, ncname, ncid,                    &
     &                      Vname(1,idTsur(itrc)), varid,               &
     &                      rec, gtype, Vsize,                          &
     &                      LBi, UBi, LBj, UBj, 1,Nfrec(ng),            &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      rmask,                                      &
#   endif
     &                      s_tflux(:,:,:,Lwrk,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))), rec,       &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
!
!  If multiple files, close current gfile.
!
      IF (ndef.gt.0) THEN
        CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state
#endif
      END MODULE inner2state_mod

