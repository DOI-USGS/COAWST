#include "cppdefs.h"
      MODULE ini_lanczos_mod

#ifdef OBS_SENSITIVITY
!
!svn $Id: ini_lanczos.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group       Andrew M. Moore   !
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
      CALL wclock_on (ng, iTLM, 1)
# endif
      CALL ini_lanczos_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Ladj, Lini,                                &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef ADJUST_WSTRESS
     &                       OCEAN(ng) % ad_ustr,                       &
     &                       OCEAN(ng) % ad_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       OCEAN(ng) % ad_tflux,                      &
#  endif
     &                       OCEAN(ng) % ad_t,                          &
     &                       OCEAN(ng) % ad_u,                          &
     &                       OCEAN(ng) % ad_v,                          &
# else
     &                       OCEAN(ng) % ad_ubar,                       &
     &                       OCEAN(ng) % ad_vbar,                       &
# endif
     &                       OCEAN(ng) % ad_zeta,                       &
# ifdef ADJUST_WSTRESS
     &                       OCEAN(ng) % tl_ustr,                       &
     &                       OCEAN(ng) % tl_vstr,                       &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       OCEAN(ng) % tl_tflux,                      &
#  endif
     &                       OCEAN(ng) % tl_t,                          &
     &                       OCEAN(ng) % tl_u,                          &
     &                       OCEAN(ng) % tl_v,                          &
# else
     &                       OCEAN(ng) % tl_ubar,                       &
     &                       OCEAN(ng) % tl_vbar,                       &
# endif
     &                       OCEAN(ng) % tl_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 1)
# endif

      RETURN
      END SUBROUTINE ini_lanczos
!
!***********************************************************************
      SUBROUTINE ini_lanczos_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             Ladj, Lini,                          &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef ADJUST_WSTRESS
     &                             ad_ustr, ad_vstr,                    &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                             ad_tflux,                            &
#  endif
     &                             ad_t, ad_u, ad_v,                    &
# else
     &                             ad_ubar, ad_vbar,                    &
# endif
     &                             ad_zeta,                             &
# ifdef ADJUST_WSTRESS
     &                             tl_ustr, tl_vstr,                    &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                             tl_tflux,                            &
#  endif
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
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Ladj, Lini
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#   endif
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
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#   endif
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
      integer :: Lwrk, i, j, lstr, ndefHSS, outLoop, rec
# ifdef SOLVE3D
      integer :: itrc, k
# endif
      real(r8) :: fac, fac1, fac2
      real(r8) :: zbeta

      real(r8), dimension(0:NstateVar(ng)) :: dot
      real(r8), dimension(Ninner) :: DotProd
      real(r8), dimension(Ninner) :: bvector
      real(r8), dimension(Ninner) :: zgamma

      character (len=80) :: ncname

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Compute tangent linear model initial conditions from the weighted
!  sum of the Lanczos vectors.
!------------------------------------------------------------------------
!
!  Determine if single or multiple Lanczos vector NetCDF files.
!
      CALL netcdf_get_ivar (ng, iADM, TRIM(HSSname(ng)), 'ndefADJ',      &
     &                      ndefHSS)
      IF (exit_flag.ne.NoError) RETURN
! 
      Lwrk=1
      DO inner=1,Ninner                 ! last record ignored
!
!  Determine Lanczos vector file to read.  The Lanczos vectors are
!  written into the adjoint NetCDF in the IS4DVAR Lanczos algorithm.
!  The Lanczos vector for each inner loop is accumulated in the 
!  unlimited dimension. The name of this file is provided here in
!  the HSSname variable since the ADJname value will be use in the
!  adjoint sensitivity part.
!
        IF (ndefHSS.gt.0) THEN
          lstr=LEN_TRIM(HSSname(ng))
          WRITE (ncname,10) HSSname(ng)(1:lstr-8), inner
 10       FORMAT (a,'_',i4.4,'.nc')
        ELSE
          ncname=HSSname(ng)
        END IF
!
!  Read in the Lanczos vectors (q_i, where i=1,2,...k) computed from
!  k inner-loops of the IS4DVAR algorithm first outer loop. Load
!  Lanczos vectors into TANGENT LINEAR STATE ARRAYS at index Lwrk.
!
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, inner,                                   &
     &                   ndefHSS, ncHSSid(ng), TRIM(ncname),            &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                   tl_ustr, tl_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   tl_tflux,                                      &
#  endif
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
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NstateVar(ng), dot(0:),                     &
# ifdef MASKING
     &                      rmask, umask, vmask,                        &
# endif
# ifdef ADJUST_WSTRESS
     &                      ad_ustr(:,:,:,Ladj), tl_ustr(:,:,:,Lwrk),   &
     &                      ad_vstr(:,:,:,Ladj), tl_vstr(:,:,:,Lwrk),   &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                      ad_tflux(:,:,:,Ladj,:),                     &
     &                      tl_tflux(:,:,:,Lwrk,:),                     &
#  endif
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
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lini, fac,                                 &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, tl_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       tl_tflux,                                  &
#  endif
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
        IF (ndefHSS.gt.0) THEN
          lstr=LEN_TRIM(HSSname(ng))
          WRITE (ncname,10) HSSname(ng)(1:lstr-8), inner
        ELSE
          ncname=HSSname(ng)
        END IF
        CALL read_state (ng, tile, iTLM,                                &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   Lwrk, inner,                                   &
     &                   ndefHSS, ncHSSid(ng), ncname,                  &
# ifdef MASKING
     &                   rmask, umask, vmask,                           &
# endif
# ifdef ADJUST_WSTRESS
     &                   ad_ustr, ad_vstr,                              &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                   ad_tflux,                                      &
#  endif
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
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lini, Lwrk, Lini, fac1, fac2,              &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_WSTRESS
     &                       tl_ustr, ad_ustr,                          &
     &                       tl_vstr, ad_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       tl_tflux, ad_tflux,                        &
#  endif
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
     &                       LBi, UBi, LBj, UBj,                        &
     &                       Lwrk, rec,                                 &
     &                       ndef, ncfileid, ncname,                    &
# ifdef MASKING
     &                       rmask, umask, vmask,                       &
# endif
# ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
#  endif
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

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
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
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#   endif
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
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#   endif
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
      integer :: i, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif
      integer :: gtype, ncid, status
      integer, dimension(NV) :: vid
      integer, dimension(4) :: Vsize

      integer :: nf_fread2d
# ifdef SOLVE3D
      integer :: nf_fread3d
# endif

      real(r8) :: Fmin, Fmax, scale

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lwrk.
!-----------------------------------------------------------------------
!
!  Determine file and variables ids.
!
      IF ((ndef.gt.0).or.(ncfileid.eq.-1)) THEN
        IF (InpThread) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF            
          ncfileid=ncid
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, ncfileid, 1)
        CALL mp_bcasti (ng, model, exit_flag, 1)
        IF (exit_flag.ne.NoError) RETURN
# endif
      ELSE
        ncid=ncfileid
      END IF
      IF (InpThread) THEN
# ifndef SOLVE3D
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUbar)), vid(idUbar))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVbar)), vid(idVbar))
# endif
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idFsur)), vid(idFsur))
# ifdef ADJUST_WSTRESS
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUsms)), vid(idUsms))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVsms)), vid(idVsms))
# endif
# ifdef SOLVE3D
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idUvel)), vid(idUvel))
        status=nf90_inq_varid(ncid, TRIM(Vname(1,idVvel)), vid(idVvel))
        DO itrc=1,NT(ng)
          status=nf90_inq_varid(ncid, TRIM(Vname(1,idTvar(itrc))),      &
     &                          vid(idTvar(itrc)))
#  ifdef ADJUST_STFLUX
          status=nf90_inq_varid(ncid, TRIM(Vname(1,idTsur(itrc))),      &
     &                          vid(idTsur(itrc)))
#  endif
        END DO
# endif
      END IF
      DO i=1,4
        Vsize(i)=0
      END DO
      scale=1.0_r8
!
!  Read in free-surface.
!
      gtype=r2dvar
      status=nf_fread2d(ng, iTLM, ncid, vid(idFsur), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  rmask(LBi,LBj),                                 &
# endif
     &                  s_zeta(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifndef SOLVE3D
!
!  Read in 2D momentum.
!
      gtype=u2dvar
      status=nf_fread2d(ng, iTLM, ncid, vid(idUbar), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask(LBi,LBj),                                 &
#  endif
     &                  s_ubar(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v2dvar
      status=nf_fread2d(ng, iTLM, ncid, vid(idVbar), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj,                      &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
#  endif
     &                  s_vbar(LBi,LBj,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idUsms), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask(LBi,LBj),                                 &
#  endif
     &                  s_ustr(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idVsms), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, Nfrec(ng),        &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
#  endif
     &                  s_vstr(LBi,LBj,1,Lwrk))
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
!  Read in 3D momentum.
!
      gtype=u3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idUvel), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask(LBi,LBj),                                 &
#  endif
     &                  s_u(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

      gtype=v3dvar
      status=nf_fread3d(ng, iTLM, ncid, vid(idVvel), rec, gtype,        &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask(LBi,LBj),                                 &
#  endif
     &                  s_v(LBi,LBj,1,Lwrk))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Read in tracers.
!
      gtype=r3dvar
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTvar(itrc)), rec,       &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj, 1, N(ng),   &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask(LBi,LBj),                               &
#  endif
     &                    s_t(LBi,LBj,1,Lwrk,itrc))
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

#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      DO itrc=1,NT(ng)
        status=nf_fread3d(ng, iTLM, ncid, vid(idTsur(itrc)), rec,       &
     &                    gtype, Vsize, LBi, UBi, LBj, UBj, 1,Nfrec(ng),&
     &                    scale, Fmin, Fmax,                            &
#   ifdef MASKING
     &                    rmask(LBi,LBj),                               &
#   endif
     &                    s_tflux(LBi,LBj,1,Lwrk,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
# endif
!
!  If multiple files, close current file.
!
      IF (InpThread.and.(ndef.gt.0)) THEN
        status=nf90_close(ncid)
      END IF
!
 10   FORMAT (' READ_STATE - unable to open NetCDF file: ',a)
 20   FORMAT (' READ_STATE - error while reading variable: ',a,2x,      &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)

      RETURN
      END SUBROUTINE read_state
#endif
      END MODULE ini_lanczos_mod

