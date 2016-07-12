#include "cppdefs.h"
      MODULE ini_lanczos_mod

#ifdef IS4DVAR_SENSITIVITY
!
!svn $Id: ini_lanczos.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the tangent linear model  initial conditions  !
!  as the weighted sum of all Lanczos vectors computed from the first  !
!  outer loop of the  IS4DVAR Lanczos  algorithm. It is used to study  !
!  the spatial/temporal impact of observations on the circulation.     !
!                                                                      !
!  Notice:                                                             !
!                                                                      !
!    (1) Additional outer loops will required different scaling and    !
!        saving of the Lanczos vectors. Currently, IS4DVAR destroys    !
!        the Lanczos vectors in each outer loop.                       !
!                                                                      !
!    (2) The IS4DVAR algorithm computes Ninner+1 Lanczos vectors in    !
!        the inner loop (0:Ninner).  The input NetCDF file contains    !
!        Ninner+1 records.  We will ignore the last record since it    !
!        has gradient information that is only relevant to the next    !
!        inner loop. The coefficients "cg_beta" and "cg_delta" take    !
!        this inner loop design into consideration.                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE

      PUBLIC :: ini_lanczos

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_lanczos (ng, tile, Ladj, Lini)
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
      integer, intent(in) :: ng, tile, Ladj, Lini
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
# endif
      CALL ini_lanczos_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Ladj, Lini,                                &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % ad_t_obc,                   &
     &                       BOUNDARY(ng) % ad_u_obc,                   &
     &                       BOUNDARY(ng) % ad_v_obc,                   &
#  endif
     &                       BOUNDARY(ng) % ad_ubar_obc,                &
     &                       BOUNDARY(ng) % ad_vbar_obc,                &
     &                       BOUNDARY(ng) % ad_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % ad_ustr,                      &
     &                       FORCES(ng) % ad_vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       FORCES(ng) % ad_tflux,                     &
# endif
# ifdef SOLVE3D
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
# else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
# endif
     &                       OCEAN(ng) % ad_zeta,                       &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       BOUNDARY(ng) % tl_t_obc,                   &
     &                       BOUNDARY(ng) % tl_u_obc,                   &
     &                       BOUNDARY(ng) % tl_v_obc,                   &
#  endif
     &                       BOUNDARY(ng) % tl_ubar_obc,                &
     &                       BOUNDARY(ng) % tl_vbar_obc,                &
     &                       BOUNDARY(ng) % tl_zeta_obc,                &
# endif
# ifdef ADJUST_WSTRESS
     &                       FORCES(ng) % tl_ustr,                      &
     &                       FORCES(ng) % tl_vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       FORCES(ng) % tl_tflux,                     &
# endif
# ifdef SOLVE3D
     &                       OCEAN(ng) % tl_t,                          &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
# else
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
# endif
     &                       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
# endif

      RETURN
      END SUBROUTINE ini_lanczos
!
!***********************************************************************
      SUBROUTINE ini_lanczos_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj, LBij, UBij,      &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Ladj, Lini,                          &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                             ad_t_obc, ad_u_obc, ad_v_obc,        &
#  endif
     &                             ad_ubar_obc, ad_vbar_obc,            &
     &                             ad_zeta_obc,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                             ad_ustr, ad_vstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                             ad_tflux,                            &
# endif
# ifdef SOLVE3D
     &                             ad_t, ad_u, ad_v,                    &
# else
     &                             ad_ubar, ad_vbar,                    &
# endif
     &                             ad_zeta,                             &
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                             tl_t_obc, tl_u_obc, tl_v_obc,        &
#  endif
     &                             tl_ubar_obc, tl_vbar_obc,            &
     &                             tl_zeta_obc,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                             tl_ustr, tl_vstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                             tl_tflux,                            &
# endif
# ifdef SOLVE3D
     &                             tl_t, tl_u, tl_v,                    &
# else
     &                             tl_ubar, tl_vbar,                    &
# endif
     &                             tl_zeta)
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
      integer, intent(in) :: Ladj, Lini
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
# endif
!
!  Local variable declarations.
!
      integer :: Lwrk, i, j, lstr, ndefLCZ, outLoop, rec
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: fac, fac1, fac2
      real(r8) :: zbeta

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Ninner) :: DotProd
      real(r8), dimension(Ninner) :: bvector
      real(r8), dimension(Ninner) :: zgamma

      character (len=256) :: ncname

# include "set_bounds.h"
!
      SourceFile='ini_lanczos.F, ini_lanczos_tile'
!
!-----------------------------------------------------------------------
!  Compute tangent linear model initial conditions from the weighted
!  sum of the Lanczos vectors.
!-----------------------------------------------------------------------
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
      CALL netcdf_get_ivar (ng, iADM, TRIM(LCZ(ng)%name), 'ndefADJ',    &
     &                      ndefLCZ)
      IF (exit_flag.ne.NoError) RETURN
!
      Lwrk=1
      DO inner=1,Ninner                 ! last record ignored
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
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the I4D-Var algorithm first outer loop. Load
!  Lanczos vectors into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
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
        IF (exit_flag.ne.NoError) RETURN
!
!  Compute dot product between the adjoint sensitivity solution, x(0),
!  and Lanczos vectors, q_i. The x(0) solution is assumed to be in
!  ADJOINT STATE ARRAYS at index Ladj.
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
     &                      ad_t_obc(:,:,:,:,Ladj,:),                   &
     &                      tl_t_obc(:,:,:,:,Lwrk,:),                   &
     &                      ad_u_obc(:,:,:,:,Ladj),                     &
     &                      tl_u_obc(:,:,:,:,Lwrk),                     &
     &                      ad_v_obc(:,:,:,:,Ladj),                     &
     &                      tl_v_obc(:,:,:,:,Lwrk),                     &
#  endif
     &                      ad_ubar_obc(:,:,:,Ladj),                    &
     &                      tl_ubar_obc(:,:,:,Lwrk),                    &
     &                      ad_vbar_obc(:,:,:,Ladj),                    &
     &                      tl_vbar_obc(:,:,:,Lwrk),                    &
     &                      ad_zeta_obc(:,:,:,Ladj),                    &
     &                      tl_zeta_obc(:,:,:,Lwrk),                    &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Ladj), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Ladj), tl_vstr(:,:,:,Lwrk),   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                      ad_tflux(:,:,:,Ladj,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
# endif
# ifdef SOLVE3D
     &                      ad_t(:,:,:,Ladj,:), tl_t(:,:,:,Lwrk,:),     &
     &                      ad_u(:,:,:,Ladj), tl_u(:,:,:,Lwrk),         &
     &                      ad_v(:,:,:,Ladj), tl_v(:,:,:,Lwrk),         &
# else
     &                      ad_ubar(:,:,Ladj), tl_ubar(:,:,Lwrk),       &
     &                      ad_vbar(:,:,Ladj), tl_vbar(:,:,Lwrk),       &
# endif
     &                      ad_zeta(:,:,Ladj), tl_zeta(:,:,Lwrk))
!
!  Store dot product.
!
        DotProd(inner)=dot(0)
      END DO
!
!-----------------------------------------------------------------------
!  Invert tri-diagonal matrix, T, associated with the Lanczos vectors.
!
!    T * b_i = a_i,       where i = 1,2,...k
!
!  Here T is (k,k) matrix computed from the IS4DVAR Lanczos algorithm
!  and b_i is the solution to the tri-diagonal system.  The Lanczos
!  algorithms coefficients (cg_beta, cg_gamma) used to build the
!  tri-diagonal system are assumed to be read elsewhere.
!-----------------------------------------------------------------------
!
!  For now, we can only use the first outer loop. A different scaling
!  is required for additional outer loops.
!
      outLoop=1
!
!  Decomposition and forward substitution.
!
      zbeta=cg_delta(1,outLoop)
      bvector(1)=DotProd(1)/zbeta
      DO i=2,Ninner
        zgamma(i)=cg_beta(i,outLoop)/zbeta
        zbeta=cg_delta(i,outLoop)-cg_beta(i,outLoop)*zgamma(i)
        bvector(i)=(DotProd(i)-cg_beta(i,outLoop)*bvector(i-1))/zbeta
      END DO
!
!  Back substitution.
!
      DO i=Ninner-1,1,-1
        bvector(i)=bvector(i)-zgamma(i+1)*bvector(i+1)
      END DO
!
!-----------------------------------------------------------------------
!  Compute Lanczos vectors weigthed sum.
!-----------------------------------------------------------------------
!
!  Initialize tangent linear state arrays: tl_var(Lini) = fac
!
      fac=0.0_r8

      CALL state_initialize (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lini, fac,                                 &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       tl_t_obc, tl_u_obc, tl_v_obc,              &
#  endif
     &                       tl_ubar_obc, tl_vbar_obc,                  &
     &                       tl_zeta_obc,                               &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, tl_vstr,                          &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       tl_tflux,                                  &
# endif
# ifdef SOLVE3D
     &                       tl_t, tl_u, tl_v,                          &
# else
     &                       tl_ubar, tl_vbar,                          &
# endif
     &                       tl_zeta)
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the IS4DVAR algorithm first outer loop. Load
!  Lanczos vectors into ADJOINT STATE ARRAYS at index Lwrk.
!
      IF (Ladj.eq.3) THEN
        Lwrk=1
      ELSE
        Lwrk=3-Ladj
      END IF
      DO inner=1,Ninner                 ! last record ignored
        IF (ndefLCZ.gt.0) THEN
          lstr=LEN_TRIM(LCZ(ng)%name)
          WRITE (ncname,10) LCZ(ng)%name(1:lstr-8), inner
        ELSE
          ncname=LCZ(ng)%name
        END IF
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj, LBij, UBij,                &
     &                   Lwrk, inner,                                   &
     &                   ndefLCZ, LCZ(ng)%ncid, ncname,                 &
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
!    tl_var(Lini) = fac1 * tl_var(Lini) + fac2 * ad_var(Lwrk)
!
!  This will become the tangent linear model initial conditions at
!  time index Lnew.
!
        fac1=1.0_r8
        fac2=bvector(inner)

        CALL state_addition (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lini, Lwrk, Lini, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                       tl_t_obc, ad_t_obc,                        &
     &                       tl_u_obc, ad_u_obc,                        &
     &                       tl_v_obc, ad_v_obc,                        &
#  endif
     &                       tl_ubar_obc, ad_ubar_obc,                  &
     &                       tl_vbar_obc, ad_vbar_obc,                  &
     &                       tl_zeta_obc, ad_zeta_obc,                  &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, ad_ustr,                          &
     &                       tl_vstr, ad_vstr,                          &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                       tl_tflux, ad_tflux,                        &
# endif
# ifdef SOLVE3D
     &                       tl_t, ad_t,                                &
     &                       tl_u, ad_u,                                &
     &                       tl_v, ad_v,                                &
# else
     &                       tl_ubar, ad_ubar,                          &
     &                       tl_vbar, ad_vbar,                          &
# endif
     &                       tl_zeta, ad_zeta)
      END DO

      RETURN
      END SUBROUTINE ini_lanczos_tile
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
      SourceFile='ini_lanczos.F, read_state'
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
!  If multiple files, close current file.
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
      END MODULE ini_lanczos_mod

