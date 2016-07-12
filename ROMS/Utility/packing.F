#include "cppdefs.h"

#if defined OPT_PERTURBATION || defined FORCING_SV
# define ENERGYNORM_SCALE
#endif
#if defined STOCHASTIC_OPT
# ifndef HESSIAN_SO
#  define ENERGYNORM_SCALE
# endif
#endif

#ifdef FULL_GRID
# define IR_RANGE IstrT,IendT
# define IU_RANGE IstrP,IendT
# define JR_RANGE JstrT,JendT
# define JV_RANGE JstrP,JendT
#else
# define IR_RANGE Istr,Iend
# define IU_RANGE IstrU,Iend
# define JR_RANGE Jstr,Jend
# define JV_RANGE JstrV,Jend
#endif

      MODULE packing_mod

#ifdef PROPAGATOR
!
!svn $Id: packing.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines pack and unpack model state varaibles into/from a    !
!  single vector to interface with  ARPACKs Arnoldi Method  for the    !
!  computation Ritz eigenfunctions.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: c_norm2
      PUBLIC :: r_norm2
# ifdef ADJOINT
#  ifdef STOCHASTIC_OPT
#   ifdef STOCH_OPT_WHITE
      PUBLIC :: ad_so_pack
#   else
      PUBLIC :: ad_so_pack_red
#   endif
#  else
      PUBLIC :: ad_pack
#  endif
      PUBLIC :: ad_unpack
# endif
# ifdef TANGENT
      PUBLIC :: tl_pack
      PUBLIC :: tl_unpack
# endif
# ifdef SO_SEMI
#  ifdef SO_SEMI_WHITE
      PUBLIC :: so_semi_white
#  else
      PUBLIC :: so_semi_red
#  endif
# endif
!
      CONTAINS
!
      SUBROUTINE c_norm2 (ng, model, Mstr, Mend,                        &
     &                    EvalueR, EvalueI, EvectorR, EvectorI,         &
     &                    state, norm2)
!
!=======================================================================
!                                                                      !
!  This function computes the Euclidean norm between the  propagator   !
!  real/imaginary Ritz eigenvalue (EvalueR, EvalueI) and eigenvector   !
!  (EvectorR, EvectorI) with state vector (state):                     !
!                                                                      !
!     norm2 = Euclidean NORM (state(:) + EvalueR * EvectorR(:) +       !
!                                        EvalueI * EvectorI(:))        !
!                                                                      !
!  WARNING: This function is only intended for serial or distributed   !
!           memory applications.  There is not tiled partitions. All   !
!           quantities are vectors. It replaces the calls to "daxpy"   !
!           and "dnrm2" from the BLAS library. This "legacy" library   !
!           gives  different  results when called inside modules and   !
!           the input arguments are pointers (specially using ifort).  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend

      real(r8), intent(in) :: EvalueR
      real(r8), intent(in) :: EvalueI

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: EvectorR(Mstr:)
      real(r8), intent(in) :: EvectorI(Mstr:)
      real(r8), intent(in) :: state(Mstr:)
# else
      real(r8), intent(in) :: EvectorR(Mstr:Mend)
      real(r8), intent(in) :: EvectorI(Mstr:Mend)
      real(r8), intent(in) :: state(Mstr:Mend)
# endif
      real(r8), intent(out) :: norm2
!
!  Local variable declarations.
!
      integer :: NSUB, is

      real(r8) :: cff, my_norm2

# ifdef DISTRIBUTE
      character (len=3) :: op_handle
# endif
!
!-----------------------------------------------------------------------
!  Compute the Euclidean norm of: state(:) + Rvalue * Rvector(:)
!-----------------------------------------------------------------------
!
!  Accumulate squared sum.
!
      my_norm2=0.0_r8
      DO is=Mstr,Mend
        cff=state(is)+EvalueR*EvectorR(is)+                             &
     &                EvalueI*EvectorI(is)
        my_norm2=my_norm2+cff*cff
      END DO
!
!  Take sum squared-root: perform global reduction.
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      NSUB=NtileX(ng)*NtileE(ng)         ! tiled application
# endif
!$OMP CRITICAL (C_NORM)
      IF (tile_count.eq.0) THEN
        norm2=my_norm2
      ELSE
        norm2=norm2+my_norm2
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, norm2, op_handle)
# endif
      END IF
!$OMP END CRITICAL (C_NORM)
      norm2=SQRT(norm2)

      RETURN
      END SUBROUTINE c_norm2
!
# if defined HESSIAN_FSV || defined HESSIAN_SO || defined HESSIAN_SV

      SUBROUTINE r_norm2 (ng, model, Mstr, Mend,                        &
     &                    Evalue, Evector, state, norm2)
!
!=======================================================================
!                                                                      !
!  This function computes the Euclidean norm between the propagator    !
!  real Ritz eigenvalue (Evalue) and eigenvector (Evector) with the    !
!  state vector (state):                                               !
!                                                                      !
!     norm2 = Euclidean NORM (state(:) + Evalue * Evector(:))          !
!                                                                      !
!  WARNING: The norm is computed by the master thread and broadcasted  !
!           to all the nodes in the group.  It is used when the state  !
!           vector is not partitioned between all nodes.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend

      real(r8), intent(in) :: Evalue

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Evector(Mstr:)
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: Evector(Mstr:Mend)
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
      real(r8), intent(out) :: norm2
!
!  Local variable declarations.
!
      integer :: NSUB, is

      real(r8) :: cff, my_norm2
!
!-----------------------------------------------------------------------
!  Compute the Euclidean norm of: state(:) + Rvalue * Rvector(:)
!-----------------------------------------------------------------------
!
!  Accumulate squared sum.
!
      IF (Master) THEN
        my_norm2=0.0_r8
        DO is=Mstr,Mend
          cff=state(is)+Evalue*Evector(is)
          my_norm2=my_norm2+cff*cff
        END DO
        norm2=SQRT(my_norm2)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, norm2)
#  endif

      RETURN
      END SUBROUTINE r_norm2

# else

      SUBROUTINE r_norm2 (ng, model, Mstr, Mend,                        &
     &                    Evalue, Evector, state, norm2)
!
!=======================================================================
!                                                                      !
!  This function computes the Euclidean norm between the propagator    !
!  real Ritz eigenvalue (Evalue) and eigenvector (Evector) with the    !
!  state vector (state):                                               !
!                                                                      !
!     norm2 = Euclidean NORM (state(:) + Evalue * Evector(:))          !
!                                                                      !
!  WARNING: This function is only intended for serial or distributed   !
!           memory applications.  There is not tiled partitions. All   !
!           quantities are vectors. It replaces the calls to "daxpy"   !
!           and "dnrm2" from the BLAS library. This "legacy" library   !
!           gives  different  results when called inside modules and   !
!           the input arguments are pointers (specially using ifort).  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend

      real(r8), intent(in) :: Evalue

#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Evector(Mstr:)
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: Evector(Mstr:Mend)
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
      real(r8), intent(out) :: norm2
!
!  Local variable declarations.
!
      integer :: NSUB, is

      real(r8) :: cff, my_norm2

#  ifdef DISTRIBUTE
      character (len=3) :: op_handle
#  endif
!
!-----------------------------------------------------------------------
!  Compute the Euclidean norm of: state(:) + Rvalue * Rvector(:)
!-----------------------------------------------------------------------
!
!  Accumulate squared sum.
!
      my_norm2=0.0_r8
      DO is=Mstr,Mend
        cff=state(is)+Evalue*Evector(is)
        my_norm2=my_norm2+cff*cff
      END DO
!
!  Take sum squared-root: perform global reduction.
!
#  ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#  else
      NSUB=NtileX(ng)*NtileE(ng)         ! tiled application
#  endif
!$OMP CRITICAL (R_NORM)
      IF (tile_count.eq.0) THEN
        norm2=my_norm2
      ELSE
        norm2=norm2+my_norm2
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, norm2, op_handle)
#  endif
      END IF
!$OMP END CRITICAL (R_NORM)
      norm2=SQRT(norm2)

      RETURN
      END SUBROUTINE r_norm2
# endif

# if defined ADJOINT && (defined FORCING_SV || defined SO_SEMI)
!
      SUBROUTINE ad_pack (ng, tile, Mstr, Mend, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp(ng),                                      &
#  ifdef SOLVE3D
     &                   nstp(ng),                                      &
#  endif
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork,                          &
#  else
     &                   Mstr, Mend, ad_state,                          &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % h,                                  &
#  ifdef SOLVE3D
     &                   GRID(ng) % Hz,                                 &
     &                   OCEAN(ng) % f_t,                               &
     &                   OCEAN(ng) % f_u,                               &
     &                   OCEAN(ng) % f_v,                               &
     &                   FORCES(ng) % ad_stflx,                         &
#  endif
     &                   OCEAN(ng) % f_ubar,                            &
     &                   OCEAN(ng) % f_vbar,                            &
     &                   OCEAN(ng) % f_zeta,                            &
     &                   FORCES(ng) % ad_sustr,                         &
     &                   FORCES(ng) % ad_svstr)

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

      RETURN
      END SUBROUTINE ad_pack
!
!***********************************************************************
      SUBROUTINE ad_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp,                                    &
#  ifdef SOLVE3D
     &                         nstp,                                    &
#  endif
     &                         Mstr, Mend, ad_state,                    &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         h,                                       &
#  ifdef SOLVE3D
     &                         Hz,                                      &
     &                         f_t, f_u, f_v, ad_stflx,                 &
#  endif
     &                         f_ubar, f_vbar,                          &
     &                         f_zeta, ad_sustr, ad_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_scalars
      USE mod_ocean
!
#  ifdef FORCING_SV
      USE ad_exchange_2d_mod
#   ifdef SOLVE3D
      USE ad_exchange_3d_mod
#   endif
#   ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : ad_mp_exchange2d
#    ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : ad_mp_exchange3d, ad_mp_exchange4d
#    endif
#   endif
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, icount, is, itrc, j, k

#  ifdef SALINITY
      integer, dimension(7+2*NT(ng)) :: offset
#  else
      integer, dimension(7+2*(NT(ng)+1)) :: offset
#  endif

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

#  ifdef SOLVE3D
      real(r8) :: cff1, cff2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
#  endif

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        ad_state(is)=Aspv
      END DO
#  endif

#  ifdef FORCING_SV
!
!  Impose adjoint periodic boundary conditions as appropriate.
!
#   ifdef DISTRIBUTE
      CALL ad_mp_exchange2d (ng, tile, iADM, 1,                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f_zeta)
#    ifndef SOLVE3D
      CALL ad_mp_exchange2d (ng, tile, iADM, 2,                         &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f_ubar, f_vbar)
#    else
      CALL ad_mp_exchange3d (ng, tile, iADM, 2,                         &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng), f_u, f_v)
      CALL ad_mp_exchange4d (ng, tile, iADM, 1,                         &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng), f_t)
#    endif
#   endif
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL ad_exchange_r2d_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj, f_zeta)
#   ifndef SOLVE3D
        CALL ad_exchange_u2d_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj, f_ubar)
        CALL ad_exchange_v2d_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj, f_vbar)
#   else
        CALL ad_exchange_u3d_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), f_u)
        CALL ad_exchange_v3d_tile (ng, tile,                            &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), f_v)
      DO itrc=1,NT(ng)
          CALL ad_exchange_r3d_tile (ng, tile,                          &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            f_t(:,:,:,itrc))
      END DO
#   endif
      END IF

#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint STATE variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
!  First clear the "offset" array.
!
      offset=0
!
#  ifdef SOLVE3D
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      END IF
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
         offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
         iadd=NwaterR(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+NwaterV(ng)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+NwaterR(ng)
          END IF
        END IF
      END DO
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      END IF
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+                        &
     &                           (Lm(ng)+2)*(Mm(ng)+1)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+                &
     &                           (Lm(ng)+2)*(Mm(ng)+2)
          END IF
        END IF
      END DO
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      END IF
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+Lm(ng)*(Mm(ng)-Voff)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+Lm(ng)*Mm(ng)
          END IF
        END IF
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=offset(isUbar)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
#    endif
#   endif
#  endif
!
!  Load adjoint of free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
        scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
        scale=1.0_r8
#  endif
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
              is=IJwaterR(i,j)+offset(isFsur)
              ad_state(is)=scale*f_zeta(i,j)
              f_zeta(i,j)=0.0_r8
            ELSE
              f_zeta(i,j)=0.0_r8
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
            ad_state(is)=scale*f_zeta(i,j)
            f_zeta(i,j)=0.0_r8
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Load adjoint of 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              ad_state(is)=scale*f_ubar(i,j)
              f_ubar(i,j)=0.0_r8
            ELSE
              f_ubar(i,j)=0.0_r8
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            ad_state(is)=scale*f_ubar(i,j)
            f_ubar(i,j)=0.0_r8
#   endif
          END DO
        END DO
      END IF
!
!  Load adjoint of 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              ad_state(is)=scale*f_vbar(i,j)
              f_vbar(i,j)=0.0_r8
            ELSE
              f_vbar(i,j)=0.0_r8
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
            ad_state(is)=scale*f_vbar(i,j)
            f_vbar(i,j)=0.0_r8
#   endif
          END DO
        END DO
      END IF

#  else
!
!  Load adjoint of 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
!
!   Compute the adjoint forcing for tl_ubar based on f_u.
!
        DO j=JR_RANGE
          DO i=IU_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IU_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
            END DO
          END DO
          DO i=IU_RANGE
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
            DO i=IU_RANGE
              f_u(i,j,k)=f_u(i,j,k)+DC(i,k)*CF(i,0)
            END DO
          END DO
          DO i=IU_RANGE
            CF(i,0)=0.0_r8
          END DO
        END DO
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Jmax=Mm(ng)
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#   endif
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+iadd
                ad_state(is)=scale*f_u(i,j,k)
                f_u(i,j,k)=0.0_r8
              ELSE
                f_u(i,j,k)=0.0_r8
              END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*f_u(i,j,k)
              f_u(i,j,k)=0.0_r8
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Load adjoint of 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
!
!   Compute the adjoint forcing for tl_vbar based on f_v.
!
        DO j=JV_RANGE
          IF (j.ge.JstrM) THEN
            DO i=IR_RANGE
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
              DO i=IR_RANGE
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
              END DO
            END DO
            DO i=IR_RANGE
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
              DO i=IR_RANGE
                f_v(i,j,k)=f_v(i,j,k)+DC(i,k)*CF(i,0)
              END DO
            END DO
            DO i=IR_RANGE
              CF(i,0)=0.0_r8
            END DO
          END IF
        END DO
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-Voff
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#   endif
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+iadd
                ad_state(is)=scale*f_v(i,j,k)
                f_v(i,j,k)=0.0_r8
              ELSE
                f_v(i,j,k)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*f_v(i,j,k)
              f_v(i,j,k)=0.0_r8
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Load adjoint of tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
          IF (itrc.eq.itemp) THEN
            cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
          ELSE IF (itrc.eq.isalt) THEN
            cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
          ELSE
            cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
          END IF
#   else
          scale=1.0_r8
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#   endif
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
                  is=IJwaterR(i,j)+iadd
                  ad_state(is)=scale*f_t(i,j,k,itrc)
                  f_t(i,j,k,itrc)=0.0_r8
                ELSE
                  f_t(i,j,k,itrc)=0.0_r8
                END IF
#   else
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                ad_state(is)=scale*f_t(i,j,k,itrc)
                f_t(i,j,k,itrc)=0.0_r8
#   endif
              END DO
            END DO
          END DO
        END IF
     END DO
#  endif
!
!  Load adjoint of surface U-stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              ad_state(is)=scale*ad_sustr(i,j)
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            ad_state(is)=scale*ad_sustr(i,j)
#  endif
          END DO
        END DO
      END IF
!
!  Load adjoint of surface V-stress.
!
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              ad_state(is)=scale*ad_svstr(i,j)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            ad_state(is)=scale*ad_svstr(i,j)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Load adjoint of tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTsur(itrc))
                ad_state(is)=scale*ad_stflx(i,j,itrc)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTsur(itrc))
              ad_state(is)=scale*ad_stflx(i,j,itrc)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_pack_tile

# elif defined STOCHASTIC_OPT
#  ifdef STOCH_OPT_WHITE
!
      SUBROUTINE ad_so_pack (ng, tile, Mstr, Mend, IntTrap, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: IntTrap
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_so_pack_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      kstp(ng),                                   &
#  ifdef SOLVE3D
     &                      nstp(ng),                                   &
#  endif
     &                      IntTrap,                                    &
#  ifdef DISTRIBUTE
     &                      1, Mstate(ng), Swork,                       &
#  else
     &                      Mstr, Mend, ad_state,                       &
#  endif
#  ifdef MASKING
     &                      GRID(ng) % IJwaterR,                        &
     &                      GRID(ng) % IJwaterU,                        &
     &                      GRID(ng) % IJwaterV,                        &
     &                      GRID(ng) % rmask,                           &
     &                      GRID(ng) % umask,                           &
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      GRID(ng) % h,                               &
#  ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
     &                      FORCES(ng) % ad_stflx,                      &
#  else
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
#  endif
     &                      OCEAN(ng) % ad_zeta,                        &
     &                      FORCES(ng) % ad_sustr,                      &
     &                      FORCES(ng) % ad_svstr)

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

      RETURN
      END SUBROUTINE ad_so_pack
!
!***********************************************************************
      SUBROUTINE ad_so_pack_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            kstp,                                 &
#  ifdef SOLVE3D
     &                            nstp,                                 &
#  endif
     &                            IntTrap,                              &
     &                            Mstr, Mend, ad_state,                 &
#  ifdef MASKING
     &                            IJwaterR, IJwaterU, IJwaterV,         &
     &                            rmask, umask, vmask,                  &
#  endif
     &                            h,                                    &
#  ifdef SOLVE3D
     &                            Hz,                                   &
     &                            ad_t, ad_u, ad_v, ad_stflx,           &
#  else
     &                            ad_ubar, ad_vbar,                     &
#  endif
     &                            ad_zeta, ad_sustr, ad_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_scalars
      USE mod_ocean
      USE mod_storage
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
      integer, intent(in) :: IntTrap
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, icount, is, itrc, j, k

#  ifdef SALINITY
      integer, dimension(7+2*NT(ng)) :: offset
#  else
      integer, dimension(7+2*(NT(ng)+1)) :: offset
#  endif

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, cff1, scale

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        ad_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint STATE variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Initialize local summations.
!
      IF (IntTrap.eq.1) THEN
        DO is=Mstr,Mend
          STORAGE(ng)%ad_Work(is)=0.0_r8
        END DO
      END IF
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
!  First clear the "offset" array.
!
      offset=0
!
#  ifdef SOLVE3D
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      END IF
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
         offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
         iadd=NwaterR(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+NwaterV(ng)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+NwaterR(ng)
          END IF
        END IF
      END DO
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      END IF
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+                        &
     &                           (Lm(ng)+2)*(Mm(ng)+1)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+                &
     &                           (Lm(ng)+2)*(Mm(ng)+2)
          END IF
        END IF
      END DO
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      END IF
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+Lm(ng)*(Mm(ng)-Voff)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+Lm(ng)*Mm(ng)
          END IF
        END IF
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=offset(isUbar)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
#    endif
#   endif
#  endif
!
      cff1=dt(ng)*REAL(ntimes(ng)/Nintervals,r8)
      IF ((IntTrap.eq.1).or.(IntTrap.eq.Nintervals+1)) THEN
        cff1=0.5_r8*cff1
      ENDIF
!
!  Load adjoint of free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
        scale=cff1
#  ifdef ENERGYNORM_SCALE
        scale=scale/SQRT(0.5_r8*g*rho0)
#  endif
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
              is=IJwaterR(i,j)+offset(isFsur)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_zeta(i,j,kstp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
            ad_state(is)=STORAGE(ng)%ad_Work(is)+                       &
     &                              scale*ad_zeta(i,j,kstp)
            STORAGE(ng)%ad_Work(is)=ad_state(is)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Load adjoint of 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=cff1/SQRT(cff*(h(i-1,j)+h(i,j)))
#   else
            scale=cff1
#   endif
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_ubar(i,j,kstp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            ad_state(is)=STORAGE(ng)%ad_Work(is)+                       &
     &                   scale*ad_ubar(i,j,kstp)
            STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
          END DO
        END DO
      END IF
!
!  Load adjoint of 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=cff1/SQRT(cff*(h(i,j-1)+h(i,j)))
#   else
            scale=cff1
#   endif
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_vbar(i,j,kstp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
            END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_vbar(i,j,kstp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
          END DO
        END DO
      END IF

#  else
!
!  Load adjoint of 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Jmax=Mm(ng)
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=cff1/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    else
                scale=cff1
#    endif
                is=IJwaterU(i,j)+iadd
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scale*ad_u(i,j,k,nstp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
                scale=cff1/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    else
                scale=cff1
#    endif
                is=(i-Ioff)+(j-Joff)*Imax+iadd
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scale*ad_u(i,j,k,nstp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Load adjoint of 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-Voff
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=cff1/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    else
                scale=cff1
#    endif
                is=IJwaterV(i,j)+iadd
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scale*ad_v(i,j,k,nstp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=cff1/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    else
              scale=cff1
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_v(i,j,k,nstp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Load adjoint of tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
          IF (itrc.eq.itemp) THEN
            cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
          ELSE IF (itrc.eq.isalt) THEN
            cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
          ELSE
            cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
          END IF
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                  scale=cff1/SQRT(cff*Hz(i,j,k))
#    else
                  scale=cff1
#    endif

                  is=IJwaterR(i,j)+iadd
                  ad_state(is)=STORAGE(ng)%ad_Work(is)+                 &
     &                         scale*ad_t(i,j,k,nstp,itrc)
                  STORAGE(ng)%ad_Work(is)=ad_state(is)
                END IF
#   else
#    ifdef ENERGYNORM_SCALE
                scale=cff1/SQRT(cff*Hz(i,j,k))
#    else
                scale=cff1
#    endif
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scale*ad_t(i,j,k,nstp,itrc)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
     END DO
#  endif
!
!  Load adjoint of surface U-stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#   endif
#  endif
        scale=cff1
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_sustr(i,j)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
            END IF
#  else
              is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_sustr(i,j)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Load adjoint of surface V-stress.
!
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#   endif
#  endif
        scale=cff1
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_svstr(i,j)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            ad_state(is)=STORAGE(ng)%ad_Work(is)+                       &
     &                   scale*ad_svstr(i,j)
            STORAGE(ng)%ad_Work(is)=ad_state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Load adjoint of surface tracer flux variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=cff1
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTsur(itrc))
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scale*ad_stflx(i,j,itrc)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTsur(itrc))
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scale*ad_stflx(i,j,itrc)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_so_pack_tile

#  else
!
      SUBROUTINE ad_so_pack_red (ng, tile, Mstr, Mend, IntTrap,         &
     &                           ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: IntTrap
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_so_pack_red_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          kstp(ng),                               &
#  ifdef SOLVE3D
     &                          nstp(ng),                               &
#  endif
     &                          IntTrap,                                &
#  ifdef DISTRIBUTE
     &                          1, Mstate(ng), Swork,                   &
#  else
     &                          Mstr, Mend, ad_state,                   &
#  endif
#  ifdef MASKING
     &                          GRID(ng) % IJwaterR,                    &
     &                          GRID(ng) % IJwaterU,                    &
     &                          GRID(ng) % IJwaterV,                    &
     &                          GRID(ng) % rmask,                       &
     &                          GRID(ng) % umask,                       &
     &                          GRID(ng) % vmask,                       &
#  endif
     &                          GRID(ng) % h,                           &
#  ifdef SOLVE3D
     &                          GRID(ng) % Hz,                          &
     &                          OCEAN(ng) % ad_t,                       &
     &                          OCEAN(ng) % ad_u,                       &
     &                          OCEAN(ng) % ad_v,                       &
     &                          FORCES(ng) % ad_stflx,                  &
#  else
     &                          OCEAN(ng) % ad_ubar,                    &
     &                          OCEAN(ng) % ad_vbar,                    &
#  endif
     &                          OCEAN(ng) % ad_zeta,                    &
     &                          FORCES(ng) % ad_sustr,                  &
     &                          FORCES(ng) % ad_svstr)

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

      RETURN
      END SUBROUTINE ad_so_pack_red
!
!***********************************************************************
      SUBROUTINE ad_so_pack_red_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                kstp,                             &
#  ifdef SOLVE3D
     &                                nstp,                             &
#  endif
     &                                IntTrap,                          &
     &                                Mstr, Mend, ad_state,             &
#  ifdef MASKING
     &                                IJwaterR, IJwaterU, IJwaterV,     &
     &                                rmask, umask, vmask,              &
#  endif
     &                                h,                                &
#  ifdef SOLVE3D
     &                                Hz,                               &
     &                                ad_t, ad_u, ad_v, ad_stflx,       &
#  else
     &                                ad_ubar, ad_vbar,                 &
#  endif
     &                                ad_zeta, ad_sustr, ad_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_storage
!
      USE nf_fread2d_mod, ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod, ONLY : nf_fread3d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
      integer, intent(in) :: IntTrap
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, icount, is, itrc, j, k
      integer :: Fcount, Irec, Nrec
      integer :: gtype, status
      integer :: Vsize(4), ntAD, ntTL, Iinp

#  ifdef SALINITY
      integer, dimension(7+2*NT(ng)) :: offset
#  else
      integer, dimension(7+2*(NT(ng)+1)) :: offset
#  endif

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: Fmin, Fmax, scale

      real(r8) :: cff, cff1, afac, scalev

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        ad_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint STATE variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Initialize local summations.
!
      IF (IntTrap.eq.1) THEN
        DO is=Mstr,Mend
          STORAGE(ng)%ad_Work(is)=0.0_r8
        END DO
      END IF
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
!  First clear the "offset" array.
!
      offset=0
!
#  ifdef SOLVE3D
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      END IF
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
         offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
         iadd=NwaterR(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+NwaterV(ng)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+NwaterR(ng)
          END IF
        END IF
      END DO
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      END IF
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+(Lm(ng)+2)*(Mm(ng)+1)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+                &
     &                           (Lm(ng)+2)*(Mm(ng)+2)
          END IF
        END IF
      END DO
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      END IF
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+Lm(ng)*(Mm(ng)-Voff)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+Lm(ng)*Mm(ng)
          END IF
        END IF
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=offset(isUbar)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
#    endif
#   endif
#  endif
!
!  Read in adjoint state variables records.
!
      cff1=dt(ng)*REAL(ntimes(ng)/Nintervals,r8)
      cff1=cff1*cff1
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Loop over the number of number of ADJ netcdf file records.
!
      Iinp=1
      scale=1.0_r8
      ntTL=(IntTrap-1)*nADJ(ng)+1
      Fcount=ADM(ng)%Fcount
      Nrec=ADM(ng)%Nrec(Fcount)
!
      ADREC_LOOP : DO Irec=1,Nrec
!
!  Autocorrelation factor.
!
        ntAD=(Nrec-Irec)*nADJ(ng)+1
        CALL sp_bcoef (ng, ntAD, ntTL, afac)
!AMM    afac=afac*cff1
!
!  Process free-surface.
!
        IF (SCALARS(ng)%Fstate(isFsur)) THEN
          gtype=r2dvar
          status=nf_fread2d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,       &
     &                      Vname(1,idFsur), ADM(ng)%Vid(idFsur),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      ad_zeta(:,:,Iinp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
          Imax=Lm(ng)+2
          Ioff=1
          Joff=0
#   else
          Imax=Lm(ng)
          Ioff=0
          Joff=1
#   endif
#  endif
          scalev=afac
#  if defined ENERGYNORM_SCALE
          scalev=scalev/SQRT(0.5_r8*g*rho0)
#  endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#  ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isFsur)
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scalev*ad_zeta(i,j,Iinp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#  else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scalev*ad_zeta(i,j,Iinp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#  endif
            END DO
          END DO
        END IF

#  ifndef SOLVE3D
!
!  Process 2D U-momentum.
!
        IF (SCALARS(ng)%Fstate(isUbar)) THEN
          gtype=u2dvar
          status=nf_fread2d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,       &
     &                      Vname(1,idUbar), ADM(ng)%Vid(idUbar),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      ad_ubar(:,:,Iinp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+1
          Ioff=0
          Joff=0
#    else
          Imax=Lm(ng)-Uoff
          Ioff=Uoff
          Joff=1
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
          cff=0.25_r8*rho0
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   if defined ENERGYNORM_SCALE
              scalev=afac/SQRT(cff*(h(i-1,j)+h(i,j)))
#   else
              scalev=afac
#   endif
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+offset(isUbar)
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scalev*ad_ubar(i,j,Iinp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scalev*ad_ubar(i,j,Iinp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
            END DO
          END DO
        END IF
!
!  Process 2D V-momentum.
!
        IF (SCALARS(ng)%Fstate(isVbar)) THEN
          gtype=v2dvar
          status=nf_fread2d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,       &
     &                      Vname(1,idVbar), ADM(ng)%Vid(idVbar),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      ad_vbar(:,:,Iinp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Ioff=1
          Joff=1
#    else
          Imax=Lm(ng)
          Ioff=0
          Joff=1+Voff
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
          cff=0.25_r8*rho0
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   if defined ENERGYNORM_SCALE
              scalev=afac/SQRT(cff*(h(i,j-1)+h(i,j)))
#   else
              scalev=afac
#   endif
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+offset(isVbar)
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scalev*ad_vbar(i,j,Iinp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     scalev*ad_vbar(i,j,Iinp)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
            END DO
          END DO
        END IF

#  else
!
!  Process 3D U-momentum.
!
        IF (SCALARS(ng)%Fstate(isUvel)) THEN
          gtype=u3dvar
          status=nf_fread3d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,       &
     &                      Vname(1,idUvel), ADM(ng)%Vid(idUvel),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      ad_u(:,:,:,Iinp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+1
          Jmax=Mm(ng)+2
          Ioff=0
          Joff=0
#    else
          Imax=Lm(ng)-Uoff
          Jmax=Mm(ng)
          Ioff=Uoff
          Joff=1
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
          cff=0.25_r8*rho0
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
            iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
            DO j=JR_RANGE
              DO i=IU_RANGE
#   ifdef MASKING
                IF (umask(i,j).gt.0.0_r8) THEN
#    if defined ENERGYNORM_SCALE
                  scalev=afac/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    else
                  scalev=afac
#    endif
                  is=IJwaterU(i,j)+iadd
                  ad_state(is)=STORAGE(ng)%ad_Work(is)+                 &
     &                         scalev*ad_u(i,j,k,Iinp)
                  STORAGE(ng)%ad_Work(is)=ad_state(is)
                END IF
#   else
#    if defined ENERGYNORM_SCALE
                scalev=afac/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    else
                scalev=afac
#    endif
                is=(i-Ioff)+(j-Joff)*Imax+iadd
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scalev*ad_u(i,j,k,Iinp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
!
!  Process 3D V-momentum.
!
        IF (SCALARS(ng)%Fstate(isVvel)) THEN
          gtype=v3dvar
          status=nf_fread3d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,       &
     &                      Vname(1,idVvel), ADM(ng)%Vid(idVvel),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      ad_v(:,:,:,Iinp))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+1
          Ioff=1
          Joff=1
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)-Voff
          Ioff=0
          Joff=1+Voff
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
          cff=0.25_r8*rho0
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
            iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
            DO j=JV_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (vmask(i,j).gt.0.0_r8) THEN
#    if defined ENERGYNORM_SCALE
                  scalev=afac/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    else
                  scalev=afac
#    endif
                  is=IJwaterV(i,j)+iadd
                  ad_state(is)=STORAGE(ng)%ad_Work(is)+                 &
     &                         scalev*ad_v(i,j,k,Iinp)
                  STORAGE(ng)%ad_Work(is)=ad_state(is)
                END IF
#   else
#    if defined ENERGYNORM_SCALE
                scalev=afac/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    else
                scalev=afac
#    endif
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       scalev*ad_v(i,j,k,Iinp)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
!
!  Process tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
            gtype=r3dvar
            status=nf_fread3d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,     &
     &                        Vname(1,idTvar(itrc)), ADM(ng)%Tid(itrc), &
     &                        Irec, gtype, Vsize,                       &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        ad_t(:,:,:,Iinp,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), Irec,    &
     &                            TRIM(ADM(ng)%name)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
!  Load adjoint of tracers variables.
!
#   ifndef MASKING
#    ifdef FULL_GRID
            Imax=Lm(ng)+2
            Jmax=Mm(ng)+2
            Ioff=1
            Joff=0
#    else
            Imax=Lm(ng)
            Jmax=Mm(ng)
            Ioff=0
            Joff=1
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
            IF (itrc.eq.itemp) THEN
              cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
            ELSE IF (itrc.eq.isalt) THEN
              cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
            ELSE
              cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
            END IF
#   endif
            DO k=1,N(ng)
#   ifdef MASKING
              iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
              iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
              DO j=JR_RANGE
                DO i=IR_RANGE
#   ifdef MASKING
                  IF (rmask(i,j).gt.0.0_r8) THEN
#    if defined ENERGYNORM_SCALE
                    scalev=afac/SQRT(cff*Hz(i,j,k))
#    else
                    scalev=afac
#    endif
                    is=IJwaterR(i,j)+iadd
                    ad_state(is)=STORAGE(ng)%ad_Work(is)+               &
     &                           scalev*ad_t(i,j,k,Iinp,itrc)
                    STORAGE(ng)%ad_Work(is)=ad_state(is)
                  END IF
#   else
#    if defined ENERGYNORM_SCALE
                  scalev=afac/SQRT(cff*Hz(i,j,k))
#    else
                  scalev=afac
#    endif
                  is=(i+Ioff)+(j-Joff)*Imax+iadd
                  ad_state(is)=STORAGE(ng)%ad_Work(is)+                 &
     &                         scalev*ad_t(i,j,k,Iinp,itrc)
                  STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
!
!  Process 2D U-momentum stress.
!
        IF (SCALARS(ng)%Fstate(isUstr)) THEN
          gtype=u2dvar
          status=nf_fread2d(ng, iADM, ADM(ng)%name,  ADM(ng)%ncid,      &
     &                      Vname(1,idUsms), ADM(ng)%Vid(idUsms),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask,                           &
#  endif
     &                      ad_sustr(:,:))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUsms)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of surface U-stress.
!
#  ifndef MASKING
#   ifdef FULL_GRID
          Imax=Lm(ng)+1
          Ioff=0
          Joff=0
#   else
          Imax=Lm(ng)-Uoff
          Ioff=Uoff
          Joff=1
#   endif
#  endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#  ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+offset(isUstr)
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       afac*ad_sustr(i,j)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#  else
              is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     afac*ad_sustr(i,j)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#  endif
            END DO
          END DO
        END IF
!
!  Process 2D V-momentum stress.
!
        IF (SCALARS(ng)%Fstate(isVstr)) THEN
          gtype=v2dvar
          status=nf_fread2d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,       &
     &                      Vname(1,idVsms), ADM(ng)%Vid(idVsms),       &
     &                      Irec, gtype, Vsize,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask,                           &
#  endif
     &                      ad_svstr(:,:))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVsms)), Irec,            &
     &                          TRIM(ADM(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Load adjoint of surface V-stress.
!
#  ifndef MASKING
#   ifdef FULL_GRID
          Imax=Lm(ng)+2
          Ioff=1
          Joff=1
#   else
          Imax=Lm(ng)
          Ioff=0
          Joff=1+Voff
#   endif
#  endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#  ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+offset(isVstr)
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       afac*ad_svstr(i,j)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
              END IF
#  else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
              ad_state(is)=STORAGE(ng)%ad_Work(is)+                     &
     &                     afac*ad_svstr(i,j)
              STORAGE(ng)%ad_Work(is)=ad_state(is)
#  endif
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
!
!  Process surface tracer flux.
!
        DO itrc=1,NT(ng)
          IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
            gtype=r2dvar
            status=nf_fread2d(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,     &
     &                        Vname(1,idTsur(itrc)),                    &
     &                        ADM(ng)%Vid(idTsur(itrc)),                &
     &                        Irec, gtype, Vsize,                       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, Fmin, Fmax,                        &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        ad_stflx(:,:,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))), Irec,    &
     &                            TRIM(ADM(ng)%name)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
!  Load surface tracer flux.
!
#   ifndef MASKING
#    ifdef FULL_GRID
            Imax=Lm(ng)+2
            Jmax=Mm(ng)+2
            Ioff=1
            Joff=0
#    else
            Imax=Lm(ng)
            Jmax=Mm(ng)
            Ioff=0
            Joff=1
#    endif
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
                  is=IJwaterR(i,j)+offset(isTsur(itrc))
                  ad_state(is)=STORAGE(ng)%ad_Work(is)+                 &
     &                         afac*ad_stflx(i,j,itrc)
                  STORAGE(ng)%ad_Work(is)=ad_state(is)
                END IF
#   else
                is=(i+Ioff)+(j-Joff)*Imax+offset(isTsur(itrc))
                ad_state(is)=STORAGE(ng)%ad_Work(is)+                   &
     &                       afac*ad_stflx(i,j,itrc)
                STORAGE(ng)%ad_Work(is)=ad_state(is)
#   endif
              END DO
            END DO
          END IF
        END DO
#  endif
      END DO ADREC_LOOP
!
  10  FORMAT (/,' AD_SO_PACK_RED - error while reading variable: ',a,2x,&
     &        'at time record = ',i3,/,17x,'in input NetCDF file: ',a)

      RETURN
      END SUBROUTINE ad_so_pack_red_tile

#  endif

# elif defined ADJOINT
!
      SUBROUTINE ad_pack (ng, tile, Mstr, Mend, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp(ng),                                      &
#  ifdef SOLVE3D
     &                   nstp(ng),                                      &
#  endif
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork,                          &
#  else
     &                   Mstr, Mend, ad_state,                          &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % h,                                  &
#  ifdef SOLVE3D
     &                   GRID(ng) % Hz,                                 &
     &                   OCEAN(ng) % ad_t,                              &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
#  else
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % ad_vbar,                           &
#  endif
     &                   OCEAN(ng) % ad_zeta)

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

      RETURN
      END SUBROUTINE ad_pack
!
!***********************************************************************
      SUBROUTINE ad_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp,                                    &
#  ifdef SOLVE3D
     &                         nstp,                                    &
#  endif
     &                         Mstr, Mend, ad_state,                    &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         h,                                       &
#  ifdef SOLVE3D
     &                         Hz,                                      &
     &                         ad_t, ad_u, ad_v,                        &
#  else
     &                         ad_ubar, ad_vbar,                        &
#  endif
     &                         ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        ad_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint STATE variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#    ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Load adjoint free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            ad_state(is)=scale*ad_zeta(i,j,kstp)
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          ad_state(is)=scale*ad_zeta(i,j,kstp)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  Load adjoint 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            ad_state(is)=scale*ad_ubar(i,j,kstp)
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          ad_state(is)=scale*ad_ubar(i,j,kstp)
#   endif
        END DO
      END DO
!
!  Load adjoint 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            ad_state(is)=scale*ad_vbar(i,j,kstp)
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          ad_state(is)=scale*ad_vbar(i,j,kstp)
#   endif
        END DO
      END DO
#  else
!
!  Load adjoint 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Jmax=Mm(ng)
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              ad_state(is)=scale*ad_u(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_u(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load adjoint 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-Voff
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              ad_state(is)=scale*ad_v(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_v(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load adjoint tracers variables. For now, use salinity scale for
!  passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif
      RETURN
      END SUBROUTINE ad_pack_tile

# endif

# if defined ADJOINT && (defined SO_SEMI || defined STOCHASTIC_OPT)
!
      SUBROUTINE ad_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the requested adjoint state and/or  surface    !
!  forcing variables used in stochastic optimals.  The state vector    !
!  contains only interior water points.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) adjoint state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iNLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif

#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
#  ifdef STOCHASTIC_OPT
     &                     knew(ng),                                    &
#  else
     &                     kstp(ng),                                    &
#  endif
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
#  ifdef ENERGYNORM_SCALE
     &                     GRID(ng) % h,                                &
#   ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
#   endif
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork)
#  else
     &                     Mstr, Mend, state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

      RETURN
      END SUBROUTINE ad_unpack
!
!***********************************************************************
      SUBROUTINE ad_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kout,                                  &
#  ifdef SOLVE3D
     &                           nout,                                  &
#  endif
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
#  ifdef ENERGYNORM_SCALE
     &                           h,                                     &
#   ifdef SOLVE3D
     &                           Hz,                                    &
#   endif
#  endif
     &                           Mstr, Mend, state)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kout
#  ifdef SOLVE3D
      integer, intent(in) :: nout
#  endif
      integer, intent(in) :: Mstr, Mend
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef ENERGYNORM_SCALE
      real(r8), intent(in) :: h(LBi:,LBj:)
#    ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#    endif
#   endif
      real(r8), intent(in) :: state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ENERGYNORM_SCALE
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#    ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#    endif
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, icount, is, itrc, j, k

#  ifdef SOLVE3D
#    ifdef SALINITY
      integer, dimension(7+2*NT(ng)) :: offset
#    else
      integer, dimension(7+2*(NT(ng)+1)) :: offset
#    endif
#  else
      integer, dimension(5) :: offset
#  endif

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract adjoint FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
!  First clear the "offset" array.
!
      offset=0
!
#  ifdef SOLVE3D
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      END IF
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
         offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
         iadd=NwaterR(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+NwaterV(ng)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+NwaterR(ng)
          END IF
        END IF
      END DO
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      END IF
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+                        &
     &                           (Lm(ng)+2)*(Mm(ng)+1)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+                &
     &                           (Lm(ng)+2)*(Mm(ng)+2)
          END IF
        END IF
      END DO
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      END IF
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+Lm(ng)*(Mm(ng)-Voff)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+Lm(ng)*Mm(ng)
          END IF
        END IF
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=offset(isUbar)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
#    endif
#   endif
#  endif
!
!  Extract adjoint free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
        scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
        scale=1.0_r8
#  endif
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
              is=IJwaterR(i,j)+offset(isFsur)
              OCEAN(ng)%ad_zeta(i,j,kout)=scale*state(is)
            ELSE
              OCEAN(ng)%ad_zeta(i,j,kout)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isFsur)
            OCEAN(ng)%ad_zeta(i,j,kout)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Extract adjoint 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              OCEAN(ng)%ad_ubar(i,j,kout)=scale*state(is)
            ELSE
              OCEAN(ng)%ubar(i,j,kout)=0.0_r8
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            OCEAN(ng)%ad_ubar(i,j,kout)=scale*state(is)
#   endif
          END DO
        END DO
      END IF
!
!  Extract adjoint 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              OCEAN(ng)%ad_vbar(i,j,kout)=scale*state(is)
            ELSE
              OCEAN(ng)%ad_vbar(i,j,kout)=0.0_r8
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            OCEAN(ng)%ad_vbar(i,j,kout)=scale*state(is)
#   endif
          END DO
        END DO
      END IF

#  else
!
!  Extract adjoint 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Jmax=Mm(ng)
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
                is=IJwaterU(i,j)+iadd
                OCEAN(ng)%ad_u(i,j,k,nout)=scale*state(is)
              ELSE
                OCEAN(ng)%ad_u(i,j,k,nout)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              OCEAN(ng)%ad_u(i,j,k,nout)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract adjoint 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-Voff
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
                is=IJwaterV(i,j)+iadd
                OCEAN(ng)%ad_v(i,j,k,nout)=scale*state(is)
              ELSE
                OCEAN(ng)%ad_v(i,j,k,nout)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              OCEAN(ng)%ad_v(i,j,k,nout)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract adjoint tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
          IF (itrc.eq.itemp) THEN
            cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
          ELSE IF (itrc.eq.isalt) THEN
            cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
          ELSE
            cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
          END IF
#   else
          scale=1.0_r8
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                  scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                  is=IJwaterR(i,j)+iadd
                  OCEAN(ng)%ad_t(i,j,k,nout,itrc)=scale*state(is)
                ELSE
                  OCEAN(ng)%ad_t(i,j,k,nout,itrc)=0.0_r8
                END IF
#   else
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                OCEAN(ng)%ad_t(i,j,k,nout,itrc)=scale*state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
!
!  Extract adjoint surface U-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              FORCES(ng)%ad_sustr(i,j)=scale*state(is)
            ELSE
              FORCES(ng)%ad_sustr(i,j)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            FORCES(ng)%ad_sustr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Extract adjoint surface V-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              FORCES(ng)%ad_svstr(i,j)=scale*state(is)
            ELSE
              FORCES(ng)%ad_svstr(i,j)=0.0_r8
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            FORCES(ng)%ad_svstr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Extract adjoint surface tracer flux variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTvar(itrc))
                FORCES(ng)%ad_stflx(i,j,itrc)=scale*state(is)
              ELSE
                FORCES(ng)%ad_stflx(i,j,itrc)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTvar(itrc))
              FORCES(ng)%ad_stflx(i,j,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_unpack_tile

# elif defined ADJOINT
!
      SUBROUTINE ad_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the adjoint model variables from the state     !
!  vector. If applicable,  the state vector includes only unmasked     !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) adjoint state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL ad_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     knew(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork,                        &
#  else
     &                     Mstr, Mend, state,                           &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % h,                                &
#  ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
#  else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
#  endif
     &                     OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif
      RETURN
      END SUBROUTINE ad_unpack
!
!***********************************************************************
      SUBROUTINE ad_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
     &                           Mstr, Mend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  ifdef SOLVE3D
     &                           Hz,                                    &
     &                           ad_t, ad_u, ad_v,                      &
#  else
     &                           ad_ubar, ad_vbar,                      &
#  endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: knew
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract adjoint state variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Extract adjoint free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  if defined ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            ad_zeta(i,j,knew)=scale*state(is)
          ELSE
            ad_zeta(i,j,knew)=0.0_r8
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          ad_zeta(i,j,knew)=scale*state(is)
#  endif
        END DO
      END DO
#  ifndef SOLVE3D
!
!  Extract adjoint 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   if define ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            ad_ubar(i,j,knew)=scale*state(is)
          ELSE
            ad_ubar(i,j,knew)=0.0_r8
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          ad_ubar(i,j,knew)=scale*state(is)
#   endif
        END DO
      END DO
!
!  Extract adjoint 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   if defined ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            ad_vbar(i,j,knew)=scale*state(is)
          ELSE
            ad_vbar(i,j,knew)=0.0_r8
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          ad_vbar(i,j,knew)=scale*state(is)
#   endif
        END DO
      END DO
#  else
!
!  Extract adjoint 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Jmax=Mm(ng)
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    if defined ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              ad_u(i,j,k,nstp)=scale*state(is)
            ELSE
              ad_u(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    if defined ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            ad_u(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract adjoint 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-Voff
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   if defined ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    if defined ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              ad_v(i,j,k,nstp)=scale*state(is)
            ELSE
              ad_v(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    if defined ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            ad_v(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   if defined ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    if defined ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                ad_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                ad_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#   else
#    if defined ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_unpack_tile
# endif

# ifdef TANGENT
!
      SUBROUTINE tl_pack (ng, tile, Mstr, Mend, tl_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the tangent linear variables into the state      !
!  vetor. The state vector contains only interior water points.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tl_state(Mstr:)
#  else
      real(r8), intent(out) :: tl_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif
      CALL tl_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs(ng), kstp(ng), knew(ng),                  &
#  ifdef SOLVE3D
     &                   nstp(ng),                                      &
#  endif
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork,                          &
#  else
     &                   Mstr, Mend, tl_state,                          &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % h,                                  &
#  ifdef SOLVE3D
     &                   GRID(ng) % Hz,                                 &
     &                   OCEAN(ng) % tl_t,                              &
     &                   OCEAN(ng) % tl_u,                              &
     &                   OCEAN(ng) % tl_v,                              &
#  else
     &                   OCEAN(ng) % tl_ubar,                           &
     &                   OCEAN(ng) % tl_vbar,                           &
#  endif
     &                   OCEAN(ng) % tl_zeta)

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) tangent linear state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iTLM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, tl_state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif
      RETURN
      END SUBROUTINE tl_pack
!
!***********************************************************************
      SUBROUTINE tl_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs, kstp, knew,                        &
#  ifdef SOLVE3D
     &                         nstp,                                    &
#  endif
     &                         Mstr, Mend, tl_state,                    &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         h,                                       &
#  ifdef SOLVE3D
     &                         Hz,                                      &
     &                         tl_t, tl_u, tl_v,                        &
#  else
     &                         tl_ubar, tl_vbar,                        &
#  endif
     &                         tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: krhs, kstp, knew
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

      real(r8), intent(out) :: tl_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(out) :: tl_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize tangent linear state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        tl_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load tangent linear state variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)

      offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Load tangent linear free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            tl_state(is)=scale*tl_zeta(i,j,knew)
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          tl_state(is)=scale*tl_zeta(i,j,knew)
#  endif
        END DO
      END DO
#  ifndef SOLVE3D
!
!  Load tangent linear 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            tl_state(is)=scale*tl_ubar(i,j,knew)
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          tl_state(is)=scale*tl_ubar(i,j,knew)
#   endif
        END DO
      END DO
!
!  Load tangent linear 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            tl_state(is)=scale*tl_vbar(i,j,knew)
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          tl_state(is)=scale*tl_vbar(i,j,knew)
#   endif
        END DO
      END DO
#  else
!
!  Load tangent linear 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Jmax=Mm(ng)
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              tl_state(is)=scale*tl_u(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            tl_state(is)=scale*tl_u(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load tangent linear 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-Voff
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              tl_state(is)=scale*tl_v(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            tl_state(is)=scale*tl_v(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                tl_state(is)=scale*tl_t(i,j,k,nstp,itrc)
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_state(is)=scale*tl_t(i,j,k,nstp,itrc)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_pack_tile
# endif

# if defined TANGENT && (defined STOCHASTIC_OPT || \
     defined HESSIAN_SV )
!
      SUBROUTINE tl_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the tangent linear variables from the state    !
!  vector.  If applicable,  the state vector includes only unmasked    !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_forces
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) tangent linear state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL tl_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork,                        &
#  else
     &                     Mstr, Mend, state,                           &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % h,                                &
#  ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#  else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#  endif
     &                     OCEAN(ng) % tl_zeta,                         &
#  ifdef SOLVE3D
     &                     FORCES(ng) % tl_stflx,                       &
#  endif
     &                     FORCES(ng) % tl_sustr,                       &
     &                     FORCES(ng) % tl_svstr)

#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif

      RETURN
      END SUBROUTINE tl_unpack
!
!***********************************************************************
      SUBROUTINE tl_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
     &                           Mstr, Mend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  ifdef SOLVE3D
     &                           Hz,                                    &
     &                           tl_t, tl_u, tl_v,                      &
#  else
     &                           tl_ubar, tl_vbar,                      &
#  endif
     &                           tl_zeta,                               &
#  ifdef SOLVE3D
     &                           tl_stflx,                              &
#  endif
     &                           tl_sustr, tl_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_svstr(LBi:,LBj:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_svstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, icount, is, itrc, j, k

#  ifdef SALINITY
      integer, dimension(7+2*NT(ng)) :: offset
#  else
      integer, dimension(7+2*(NT(ng)+1)) :: offset
#  endif

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract tangent linear FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
!  First clear the "offset" array.
!
      offset=0
!
#  ifdef SOLVE3D
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      END IF
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
         offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
         iadd=NwaterR(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+NwaterV(ng)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+NwaterR(ng)
          END IF
        END IF
      END DO
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      END IF
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+                        &
     &                           (Lm(ng)+2)*(Mm(ng)+1)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+                &
     &                           (Lm(ng)+2)*(Mm(ng)+2)
          END IF
        END IF
      END DO
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      END IF
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+Lm(ng)*(Mm(ng)-Voff)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+Lm(ng)*Mm(ng)
          END IF
        END IF
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=offset(isUbar)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
#    endif
#   endif
#  endif
!
!  Extract tangent linear free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
              is=IJwaterR(i,j)+offset(isFsur)
              tl_zeta(i,j,kstp)=scale*state(is)
            ELSE
              tl_zeta(i,j,kstp)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isFsur)
            tl_zeta(i,j,kstp)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Extract tangent linear 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              tl_ubar(i,j,kstp)=scale*state(is)
            ELSE
              tl_ubar(i,j,kstp)=0.0_r8
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            tl_ubar(i,j,kstp)=scale*state(is)
#   endif
          END DO
        END DO
      END IF
!
!  Extract tangent linear 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              tl_vbar(i,j,kstp)=scale*state(is)
            ELSE
              tl_vbar(i,j,kstp)=0.0_r8
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            tl_vbar(i,j,kstp)=scale*state(is)
#   endif
          END DO
        END DO
      END IF

#  else
!
!  Extract tangent linear 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Jmax=Mm(ng)
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
                is=IJwaterU(i,j)+iadd
                tl_u(i,j,k,nstp)=scale*state(is)
              ELSE
                tl_u(i,j,k,nstp)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              tl_u(i,j,k,nstp)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract tangent linear 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-Voff
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
                is=IJwaterV(i,j)+iadd
                tl_v(i,j,k,nstp)=scale*state(is)
              ELSE
                tl_v(i,j,k,nstp)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_v(i,j,k,nstp)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract tangent linear tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                  scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                  is=IJwaterR(i,j)+iadd
                  tl_t(i,j,k,nstp,itrc)=scale*state(is)
                ELSE
                  tl_t(i,j,k,nstp,itrc)=0.0_r8
                END IF
#   else
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                tl_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
!
!  Extract tangent linear surface U-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              tl_sustr(i,j)=scale*state(is)
            ELSE
              tl_sustr(i,j)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            tl_sustr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Extract tangent linear surface V-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              tl_svstr(i,j)=scale*state(is)
            ELSE
              tl_svstr(i,j)=0.0_r8
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            tl_svstr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Extract tangent linear surface tracer flux variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTsur(itrc))
                tl_stflx(i,j,itrc)=scale*state(is)
              ELSE
                tl_stflx(i,j,itrc)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTsur(itrc))
              tl_stflx(i,j,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_unpack_tile

# elif defined TANGENT && defined FORCING_SV
!
      SUBROUTINE tl_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the tangent linear variables from the state    !
!  vector.  If applicable,  the state vector includes only unmasked    !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_forces
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) tangent linear state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL tl_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork,                        &
#  else
     &                     Mstr, Mend, state,                           &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % h,                                &
#  ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % f_t,                             &
     &                     OCEAN(ng) % f_u,                             &
     &                     OCEAN(ng) % f_v,                             &
#  endif
     &                     OCEAN(ng) % f_ubar,                          &
     &                     OCEAN(ng) % f_vbar,                          &
     &                     OCEAN(ng) % f_zeta,                          &
#  ifdef SOLVE3D
     &                     FORCES(ng) % tl_stflx,                       &
#  endif
     &                     FORCES(ng) % tl_sustr,                       &
     &                     FORCES(ng) % tl_svstr)

#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif

      RETURN
      END SUBROUTINE tl_unpack
!
!***********************************************************************
      SUBROUTINE tl_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
     &                           Mstr, Mend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  ifdef SOLVE3D
     &                           Hz,                                    &
     &                           f_t, f_u, f_v,                         &
#  endif
     &                           f_ubar, f_vbar,                        &
     &                           f_zeta,                                &
#  ifdef SOLVE3D
     &                           tl_stflx,                              &
#  endif
     &                           tl_sustr, tl_svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
#  ifdef FORCING_SV
      USE exchange_2d_mod
#   ifdef SOLVE3D
      USE exchange_3d_mod
#   endif
#   ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#    ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
#    endif
#   endif
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
      real(r8), intent(inout) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_svstr(LBi:,LBj:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_svstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, icount, is, itrc, j, k

#  ifdef SALINITY
      integer, dimension(7+2*NT(ng)) :: offset
#  else
      integer, dimension(7+2*(NT(ng)+1)) :: offset
#  endif

      real(r8) :: cff, scale

# ifdef SOLVE3D
      real(r8) :: cff1, cff2
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
# endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract tangent linear FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
!  First clear the "offset" array.
!
      offset=0
!
#  ifdef SOLVE3D
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      END IF
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
         offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
         iadd=NwaterR(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+NwaterV(ng)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+NwaterR(ng)
          END IF
        END IF
      END DO
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      END IF
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+                        &
     &                           (Lm(ng)+2)*(Mm(ng)+1)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+                &
     &                           (Lm(ng)+2)*(Mm(ng)+2)
          END IF
        END IF
      END DO
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
     END IF.
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)*N(ng)
        END IF
      END DO
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          IF (itrc.eq.1) THEN
            offset(isTsur(itrc))=offset(isVstr)+Lm(ng)*(Mm(ng)-Voff)
          ELSE
            offset(isTsur(itrc))=offset(isTsur(itrc-1))+Lm(ng)*Mm(ng)
          END IF
        END IF
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+NwaterR(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=offset(isUbar)+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+NwaterU(ng)
      END IF
#   else
#    ifdef FULL_GRID
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
      END IF
#    else
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar) THEN
        offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=0
      END IF
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isVstr)=offset(isUstr)+(Lm(ng)-Uoff)*Mm(ng)
      END IF
#    endif
#   endif
#  endif
!
!  Extract tangent linear free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
        scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
        scale=1.0_r8
#  endif
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
              is=IJwaterR(i,j)+offset(isFsur)
              f_zeta(i,j)=scale*state(is)
            ELSE
              f_zeta(i,j)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isFsur)
            f_zeta(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Extract tangent linear 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              f_ubar(i,j)=scale*state(is)
            ELSE
              f_ubar(i,j)=0.0_r8
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            f_ubar(i,j)=scale*state(is)
#   endif
          END DO
        END DO
      END IF
!
!  Extract tangent linear 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              f_vbar(i,j)=scale*state(is)
            ELSE
              f_vbar(i,j)=0.0_r8
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            f_vbar(i,j)=scale*state(is)
#   endif
          END DO
        END DO
      END IF

#  else
!
!  Extract tangent linear 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-Uoff
        Jmax=Mm(ng)
        Ioff=Uoff
        Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#   endif
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+iadd
                f_u(i,j,k)=scale*state(is)
              ELSE
                f_u(i,j,k)=0.0_r8
              END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              f_u(i,j,k)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
!
!   Compute the forcing forcing for tl_ubar based on f_u.
!
        DO j=JR_RANGE
          DO i=IU_RANGE
            DC(i,0)=0.0_r8
            CF(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=IU_RANGE
              DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*f_u(i,j,k)
            END DO
          END DO
          DO i=IU_RANGE
            cff1=1.0_r8/DC(i,0)
            cff2=CF(i,0)*cff1
#  ifdef MASKING
            cff2=cff2*umask(i,j)
#  endif
            f_ubar(i,j)=cff2
          END DO
        END DO
      END IF
!
!  Extract tangent linear 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-Voff
        Ioff=0
        Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        cff=0.25_r8*rho0
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#   endif
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+iadd
                f_v(i,j,k)=scale*state(is)
              ELSE
                f_v(i,j,k)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              f_v(i,j,k)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
!
!   Compute the forcing forcing for tl_vbar based on f_v.
!
        DO j=JV_RANGE
          IF (j.ge.JstrM) THEN
            DO i=IR_RANGE
              DC(i,0)=0.0_r8
              CF(i,0)=0.0_r8
            END DO
            DO k=1,N(ng)
              DO i=IR_RANGE
                DC(i,k)=0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))
                DC(i,0)=DC(i,0)+DC(i,k)
                CF(i,0)=CF(i,0)+DC(i,k)*f_v(i,j,k)
              END DO
            END DO
            DO i=IR_RANGE
              cff1=1.0_r8/DC(i,0)
              cff2=CF(i,0)*cff1
#  ifdef MASKING
              cff2=cff2*vmask(i,j)
#  endif
              f_vbar(i,j)=cff2
            END DO
          END IF
        END DO
      END IF
!
!  Extract tangent linear tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
          scale=1.0_r8
#   endif
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#   endif
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
                  is=IJwaterR(i,j)+iadd
                  f_t(i,j,k,itrc)=scale*state(is)
                ELSE
                  f_t(i,j,k,itrc)=0.0_r8
                END IF
#   else
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                tl_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif

#  ifdef FORCING_SV
!
!  Impose periodic boundary conditions as appropriate.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, f_zeta)
#   ifndef SOLVE3D
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, f_ubar)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, f_vbar)
#   else
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), f_u)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng), f_v)
      DO itrc=1,NT(ng)
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            f_t(:,:,:,itrc))
      END DO
#   endif
      END IF

#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f_zeta)
#    ifndef SOLVE3D
      CALL mp_exchange2d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f_ubar, f_vbar)
#    else
      CALL mp_exchange3d (ng, tile, iTLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng), f_u, f_v)
      CALL mp_exchange4d (ng, tile, iTLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 1, NT(ng),      &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng), f_t)
#    endif
#   endif
#  endif
!
!  Extract tangent linear surface U-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-Uoff
        Ioff=Uoff
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              tl_sustr(i,j)=scale*state(is)
            ELSE
              tl_sustr(i,j)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            tl_sustr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Extract tangent linear surface V-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1+Voff
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              tl_svstr(i,j)=scale*state(is)
            ELSE
              tl_svstr(i,j)=0.0_r8
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            tl_svstr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Extract tangent linear surface tracer flux variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTsur(itrc))
                tl_stflx(i,j,itrc)=scale*state(is)
              ELSE
                tl_stflx(i,j,itrc)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTsur(itrc))
              tl_stflx(i,j,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_unpack_tile

# elif defined TANGENT
!
      SUBROUTINE tl_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the tangent linear variables from the state    !
!  vector.  If applicable,  the state vector includes only unmasked    !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) tangent linear state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL tl_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork,                        &
#  else
     &                     Mstr, Mend, state,                           &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % h,                                &
#  ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#  else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#  endif
     &                     OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif
      RETURN
      END SUBROUTINE tl_unpack
!
!***********************************************************************
      SUBROUTINE tl_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
     &                           Mstr, Mend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  ifdef SOLVE3D
     &                           Hz,                                    &
     &                           tl_t, tl_u, tl_v,                      &
#  else
     &                           tl_ubar, tl_vbar,                      &
#  endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: Uoff, Voff
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract tangent linear STATE variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Set offsets for momentum variables due to periodic boundary
!  conditions. Recall that in East-West periodic boundary conditions
!  IstrU=1. Otherwise, IstrU=2. Similarly, in North-South periodic
!  applications IstrV=1 or else IstrV=2.
!
      IF (EWperiodic(ng)) THEN
        Uoff=0
      ELSE
        Uoff=1
      END IF
!
      IF (NSperiodic(ng)) THEN
        Voff=0
      ELSE
        Voff=1
      END IF
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-Uoff)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-Voff)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-Uoff)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Extract tangent linear free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            tl_zeta(i,j,kstp)=scale*state(is)
          ELSE
            tl_zeta(i,j,kstp)=0.0_r8
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          tl_zeta(i,j,kstp)=scale*state(is)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  Extract tangent linear 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            tl_ubar(i,j,kstp)=scale*state(is)
          ELSE
            tl_ubar(i,j,kstp)=0.0_r8
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          tl_ubar(i,j,kstp)=scale*state(is)
#   endif
        END DO
      END DO
!
!  Extract tangent linear 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            tl_vbar(i,j,kstp)=scale*state(is)
          ELSE
            tl_vbar(i,j,kstp)=0.0_r8
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          tl_vbar(i,j,kstp)=scale*state(is)
#   endif
        END DO
      END DO

#  else
!
!  Extract tangent linear 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-Uoff
      Jmax=Mm(ng)
      Ioff=Uoff
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              tl_u(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_u(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            tl_u(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-Voff
      Ioff=0
      Joff=1+Voff
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              tl_v(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_v(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            tl_v(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                tl_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                tl_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_unpack_tile
# endif

# ifdef SO_SEMI
#  ifdef SO_SEMI_WHITE
!
      SUBROUTINE so_semi_white (ng, tile, Mstr, Mend, state, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine computes a new stochastic optimals perturbation vector !
!  (seminorm estimation) assuming white noise forcing using ARPACK.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_storage
#   ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#   endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#   ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(out) :: ad_state(Mstr:)
#   else
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#   endif
!
!  Local variable declarations.
!
      integer :: NSUB, is, rec

      real(r8) :: SOnorm, my_SOnorm, my_TRnorm
      real(r8) :: SOnorm1, my_SOnorm1
      real(r8) :: cff, cff1, cff2

#   ifdef DISTRIBUTE
      real(r8), dimension(3) :: buffer

      character (len=3), dimension(3) :: op_handle
#   endif

#   include "tile.h"
!
!-----------------------------------------------------------------------
!  Compute seminorm, stochastic optimals adjoint perturbation vector.
!-----------------------------------------------------------------------
!
!  Initialize adjoint state vector.
!
      DO is=Mstr,Mend
        ad_state(is)=0.0_r8
      END DO
!
!  Sum over all adjoint surface forcing records available in "so_state'.
!
      IF (Master) THEN
        WRITE (stdout,'(/)')
      END IF
      my_TRnorm=0.0_r8
!
      DO rec=1,Nsemi(ng)
!
!  Compute normalization factor.
!
        cff=REAL((nADJ(ng)-1)*(2*nADJ(ng)-1),r8)/REAL(6*nADJ(ng),r8)
        cff1=1.0_r8+cff
        cff2=0.5_r8*REAL((nADJ(ng)-1))-cff
!
        my_SOnorm=0.0_r8
        my_SOnorm1=0.0_r8
        DO is=Mstr,Mend
          my_SOnorm=my_SOnorm+                                          &
     &              STORAGE(ng)%so_state(is,rec)*state(is)
        END DO
!
        IF (rec.ne.Nsemi(ng)) THEN
          DO is=Mstr,Mend
            my_SOnorm1=my_SOnorm1+                                      &
     &                 STORAGE(ng)%so_state(is,rec+1)*state(is)
            my_TRnorm=my_TRnorm+                                        &
     &                cff1*STORAGE(ng)%so_state(is,rec)*                &
     &                     STORAGE(ng)%so_state(is,rec)+                &
     &                2.0_r8*cff2*STORAGE(ng)%so_state(is,rec  )*       &
     &                            STORAGE(ng)%so_state(is,rec+1)+       &
     &                cff*STORAGE(ng)%so_state(is,rec+1)*               &
     &                    STORAGE(ng)%so_state(is,rec+1)
          END DO
        ELSE
          DO is=Mstr,Mend
            my_TRnorm=my_TRnorm+                                        &
     &                STORAGE(ng)%so_state(is,rec)*                     &
     &                STORAGE(ng)%so_state(is,rec)
          END DO
        END IF
!
!  Global reduction of normalization factor.
!
#   ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
#   else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
#   endif
!$OMP CRITICAL (SO_NORM)
        IF (tile_count.eq.0) THEN
          SOnorm=0.0_r8
          SOnorm1=0.0_r8
          IF (rec.eq.1) THEN
            TRnorm(ng)=0.0_r8
          END IF
        END IF
        SOnorm=SOnorm+my_SOnorm
        SOnorm1=SOnorm1+my_SOnorm1
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#   ifdef DISTRIBUTE
          buffer(1)=SOnorm
          buffer(2)=SOnorm1
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          CALL mp_reduce (ng, iADM, 2, buffer, op_handle)
          SOnorm=buffer(1)
          SOnorm1=buffer(2)
#   endif
        END IF
!$OMP END CRITICAL (SO_NORM)
!
!  Report normalization factors.
!
        IF (Master) THEN
          WRITE (stdout,10) rec, SOnorm, SOnorm1
 10       FORMAT (3x,'Rec = ',i2.2,2x,'SOnorm = ',1p,e15.8,0p,          &
     &            2x,'SOnorm1 = ',1p,e15.8)
        END IF
!
!  Compute new perturbation vector.
!
        IF (rec.ne.Nsemi(ng)) THEN
          DO is=Mstr,Mend
            ad_state(is)=ad_state(is)+                                  &
     &                   cff1*SOnorm *STORAGE(ng)%so_state(is,rec  )+   &
     &                   cff2*SOnorm1*STORAGE(ng)%so_state(is,rec  )+   &
     &                   cff2*SOnorm *STORAGE(ng)%so_state(is,rec+1)+   &
     &                   cff *SOnorm1*STORAGE(ng)%so_state(is,rec+1)
          END DO
        ELSE
          DO is=Mstr,Mend
            ad_state(is)=ad_state(is)+                                  &
     &                   SOnorm*STORAGE(ng)%so_state(is,rec)
          END DO
        END IF
      END DO
!
!  Global reduction of normalization factor, TRnorm.
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
!$OMP CRITICAL (TR_NORM)
      IF (tile_count.eq.0) THEN
        TRnorm(ng)=0.0_r8
      END IF
      TRnorm(ng)=TRnorm(ng)+my_TRnorm
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        op_handle(1)='SUM'
        CALL mp_reduce (ng, iADM, 1, TRnorm(ng), op_handle(1))
#   endif
      END IF
!$OMP END CRITICAL (TR_NORM)

      RETURN
      END SUBROUTINE so_semi_white

#  else
!
      SUBROUTINE so_semi_red (ng, tile, Mstr, Mend, state, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine computes a new stochastic optimals perturbation vector !
!  (seminorm estimation) assuming red noise forcing using ARPACK.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_storage
      USE mod_scalars
#   ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#   endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#   ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(out) :: ad_state(Mstr:)
#   else
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#   endif
!
!  Local variable declarations.
!
      integer :: NSUB, is, ntAD, ntTL, rec, rec1

      real(r8) :: SOnorm, my_TRnorm

      real(r8), dimension(Nsemi(ng)) :: Bcoef
      real(r8), dimension(Nsemi(ng)) :: SOdotprod
      real(r8), dimension(Nsemi(ng)) :: my_dotprod

#   ifdef DISTRIBUTE
      character (len=3), dimension(Nsemi(ng)) :: op_handle
#   endif

#   include "tile.h"
!
!-----------------------------------------------------------------------
!  Compute seminorm, stochastic optimals adjoint perturbation vector.
!-----------------------------------------------------------------------
!
!  Initialize adjoint state vector.
!
      DO is=Mstr,Mend
        ad_state(is)=0.0_r8
      END DO
!
!  First compute the dot-products.
!
      DO rec=1,Nsemi(ng)
        my_dotprod(rec)=0.0_r8
        DO is=Mstr,Mend
          my_dotprod(rec)=my_dotprod(rec)+                              &
     &                    STORAGE(ng)%so_state(is,rec)*state(is)
        END DO
      END DO
!
!  Global reduction of dot products.
!
#   ifdef DISTRIBUTE
      NSUB=1                           ! distributed-memory
#   else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                         ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
      END IF
#   endif
!$OMP CRITICAL (SO_DOT)
      IF (tile_count.eq.0) THEN
        DO rec=1,Nsemi(ng)
          SOdotprod(rec)=0.0_r8
        END DO
      END IF
      DO rec=1,Nsemi(ng)
        SOdotprod(rec)=SOdotprod(rec)+my_dotprod(rec)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        DO rec=1,Nsemi(ng)
          op_handle(rec)='SUM'
        END DO
        CALL mp_reduce (ng, iADM, Nsemi(ng), SOdotprod, op_handle)
#   endif
      END IF
!$OMP END CRITICAL (SO_DOT)
!
! Second, loop over time twice allowing for the decorrelation due to the
! red noise AR(1) process with assumed decorrelation time SOdecay.
!
      IF (Master) THEN
        WRITE (stdout,'(/)')
      END IF
      my_TRnorm=0.0_r8
!
      DO rec=1,Nsemi(ng)
        ntAD=(rec-1)*nADJ(ng)+1
        SOnorm=0.0_r8
        DO rec1=1,Nsemi(ng)
          ntTL=(rec1-1)*nADJ(ng)+1
          CALL sp_bcoef (ng, ntAD, ntTL, Bcoef(rec1))
          SOnorm=SOnorm+Bcoef(rec1)*SOdotprod(rec1)
          DO is=Mstr,Mend
            my_TRnorm=my_TRnorm+                                        &
     &                STORAGE(ng)%so_state(is,rec )*Bcoef(rec1)*        &
     &                STORAGE(ng)%so_state(is,rec1)
          END DO
        END DO
!
!  Report normalization factors.
!
        IF (Master) THEN
          WRITE (stdout,10) rec, SOdotprod(rec), Bcoef(rec), SOnorm
 10       FORMAT (1x,'Rec = ',i2.2,1x,'SOdotprod = ',1p,e13.6,0p,       &
     &            1x,'Bcoef = ',1p,e13.6,0p,1x,'SOnorm = ',1p,e13.6)
        END IF
!
!  Compute new perturbation vector.
!
        DO is=Mstr,Mend
          ad_state(is)=ad_state(is)+                                    &
     &                 SOnorm*STORAGE(ng)%so_state(is,rec)
        END DO
      END DO
!
!  Global reduction of normalization factor, TRnorm.
!
#   ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#   else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#   endif
!$OMP CRITICAL (TR_NORM)
      IF (tile_count.eq.0) THEN
        TRnorm(ng)=0.0_r8
      END IF
      TRnorm(ng)=TRnorm(ng)+my_TRnorm
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        op_handle(1)='SUM'
        CALL mp_reduce (ng, iADM, 1, TRnorm(ng), op_handle(1))
#   endif
      END IF
!$OMP END CRITICAL (TR_NORM)

      RETURN
      END SUBROUTINE so_semi_red
#  endif
# endif

# if defined SO_SEMI || !defined STOCH_OPT_WHITE
!
      SUBROUTINE sp_bcoef (ng, ntAD, ntTL, Bcoef)
!
!=======================================================================
!                                                                      !
!  This routine is used to compute  red noise stochastic processes     !
!  time-lagged  coefficient,  Bcoef,  used  to  evaluate  discrete     !
!  double-time integrals.  Currently, a discrete-time Markov chain     !
!  model is assumed with autoregressive order-one processes, AR(1).    !
!  Notice that the routine SP_ACOEF is called to compute the inner     !
!  integral.                                                           !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ntAD, ntTL

      real(r8), intent(out):: Bcoef
!
!  Local variable declarations.
!
      integer :: i, it1, it2

      real(r8) :: Acoef, Acoef1, Acoef2, df1, rov

!
!-----------------------------------------------------------------------
!  Compute red noise stochastic process time-lagged coefficient to
!  evaluate discrete double time-integrals. Currently, only Markov
!  processes, AR(1), are considered.
!-----------------------------------------------------------------------
!
!  Here, ntAD is the current model timestep and ntTL is the timestep
!  associated with forcing.
!
      rov=1.0_r8/REAL(nADJ(ng),r8)
!
      IF ((ntAD.gt.1).and.(ntAD.lt.ntimes(ng)+1)) THEN
        it1=ntAD
        it2=ntAD-nADJ(ng)
        CALL sp_acoef (ng, it1, ntTL, Acoef)
        Bcoef=Acoef
        DO i=1,nADJ(ng)-1
          CALL sp_acoef (ng, it1+i, ntTL, Acoef1)
          CALL sp_acoef (ng, it2+i, ntTL, Acoef2)
          df1=REAL(i,r8)*rov
          Bcoef=Bcoef+(1.0_r8-df1)*Acoef1+df1*Acoef2
        END DO
      ELSE IF (ntAD.eq.1) THEN
        CALL sp_acoef (ng, 1, ntTL, Acoef)
        Bcoef=Acoef
        DO i=1,nADJ(ng)-1
          CALL sp_acoef (ng, 1+i, ntTL, Acoef1)
          df1=REAL(i,r8)*rov
          Bcoef=Bcoef+(1.0_r8-df1)*Acoef1
        END DO
      ELSE IF (ntAD.eq.ntimes(ng)+1) THEN
        CALL sp_acoef (ng, ntimes(ng)+1, ntTL, Acoef)
        Bcoef=Acoef
        DO i=1,nADJ(ng)-1
          CALL sp_acoef (ng, ntimes(ng)+1-nADJ(ng)+i, ntTL, Acoef1)
          df1=REAL(i,r8)*rov
          Bcoef=Bcoef+df1*Acoef1
        END DO
      END IF

      RETURN
      END SUBROUTINE sp_bcoef
!
      SUBROUTINE sp_acoef (ng, ntAD, ntTL, Acoef)
!
!=======================================================================
!                                                                      !
!  This routine is used to compute red noise stochastic processes      !
!  time-lagged coefficient, Acoef,  used  to  evaluate inner time      !
!  integral.  Currently,  a  discrete-time Markov chain model  is      !
!  assumed with autoregressive order-one processes, AR(1). Notice      !
!  that function SP_AUTOC is used to set autocorrelation model.        !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ntAD, ntTL

      real(r8), intent(out):: Acoef
!
!  Local variable declarations.
!
      integer :: i, idf1, idf2, idf4

      real(r8) :: df3, rov
!
!-----------------------------------------------------------------------
!  Compute red noise stochastic process time-lagged coefficients to
!  evaluate discrete inner time-integral. Currently, only Markov
!  processes, AR(1), are considered.
!-----------------------------------------------------------------------
!
!  Here, ntAD is the current timestep corresponding to time when
!  solution is saved.
!
      rov=1.0_r8/REAL(nADJ(ng),r8)
      IF ((ntTL.gt.1).and.(ntTL.lt.ntimes(ng)+1)) THEN
        Acoef=0.0_r8
        DO i=1,nADJ(ng)-1
          idf1=IABS(ntAD-ntTL-i)+1
          idf2=IABS(ntAD-(ntTL-nADJ(ng))-i)+1
          df3=REAL(i,r8)*rov
          Acoef=Acoef+sp_autoc(ng,idf1)*(1.0_r8-df3)+                   &
     &                sp_autoc(ng,idf2)*df3
        END DO
        idf4=IABS(ntAD-ntTL)+1
        Acoef=Acoef+sp_autoc(ng,idf4)
      ELSE IF (ntTL.eq.1) THEN
        Acoef=0.0_r8
        DO i=1,nADJ(ng)-1
          idf1=IABS(ntAD-1-i)+1
          df3=REAL(i,r8)*rov
          Acoef=Acoef+sp_autoc(ng,idf1)*(1.0_r8-df3)
        END DO
        idf4=IABS(ntAD-1)+1
        Acoef=Acoef+sp_autoc(ng,idf4)
      ELSE IF (ntTL.eq.ntimes(ng)+1) THEN
        Acoef=0.0_r8
        DO i=1,nADJ(ng)-1
          idf2=IABS(ntAD-ntimes(ng)-1+nADJ(ng)-i)+1
          df3=REAL(i,r8)*rov
          Acoef=Acoef+sp_autoc(ng,idf2)*df3
        END DO
        idf4=IABS(ntAD-ntimes(ng)-1)+1
        Acoef=Acoef+sp_autoc(ng,idf4)
      END IF

      RETURN
      END SUBROUTINE sp_acoef
!
      FUNCTION sp_autoc (ng, idf)
!
!=======================================================================
!                                                                      !
!  This routine is used to compute red noise stochastic processes      !
!  autocorrelation model.  Notice that only  AR(1)  processes are      !
!  considered. However, other models can be easily implemented in      !
!  terms of look tables.                                               !
!
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, idf
!
!  Function result.
!
      real(r8) :: sp_autoc
!
!-----------------------------------------------------------------------
!  Set autocorrelation model.
!-----------------------------------------------------------------------
#  ifdef SO_NON_AR1
!
!  Use a user-defined temporal decorrelation function such as in the
!  form of a look-up table computed from actual data.
!
      sp_autoc=0.0_r8
#  else
!
!  Assume an AR(1) process with decorrelation time SO_decay.
!
      sp_autoc=EXP(-ABS(REAL(idf-1,r8))*dt(ng)/SO_decay(ng))
#  endif
!
      RETURN
      END FUNCTION sp_autoc
# endif
#endif

#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE

      END MODULE packing_mod
