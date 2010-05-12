#include "cppdefs.h"

#if defined OPT_PERTURBATION
# define ENERGYNORM_SCALE
#endif

#ifdef FULL_GRID
# define IR_RANGE IstrR,IendR
# define IU_RANGE Istr,IendR
# define JR_RANGE JstrR,JendR
# define JV_RANGE Jstr,JendR
#else
# define IR_RANGE Istr,Iend
# define IU_RANGE IstrU,Iend
# define JR_RANGE Jstr,Jend
# define JV_RANGE JstrV,Jend
#endif

      MODULE packing_mod

#ifdef PROPAGATOR
!
!svn $Id: packing.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group       Andrew M. Moore   !
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

      PRIVATE
# ifdef ADJOINT
      PUBLIC :: ad_pack
      PUBLIC :: ad_unpack
# endif
# ifdef TANGENT
      PUBLIC :: tl_pack
      PUBLIC :: tl_unpack
# endif
# ifdef SO_SEMI
      PUBLIC :: so_pack
      PUBLIC :: so_unpack
#  ifdef SO_SEMI_WHITE
      PUBLIC :: so_semi_white
#  else
      PUBLIC :: so_semi_red
#  endif
# endif
      CONTAINS

# ifdef ADJOINT
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
#  if defined FORCING_SV
      USE mod_forces
#  else
      USE mod_ocean
#  endif
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
      CALL wclock_on (ng, iADM, 1)
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
#  if defined FORCING_SV
     &                   FORCES(ng) % ad_sustr,                         &
     &                   FORCES(ng) % ad_svstr)
#  else
#   ifdef SOLVE3D
     &                   GRID(ng) % Hz,                                 &
     &                   OCEAN(ng) % ad_t,                              &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
#   else
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % ad_vbar,                           &
#   endif
     &                   OCEAN(ng) % ad_zeta)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 1)
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
     &                          Nstr, Nend, ad_state,                   &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         h,                                       &
#  if defined FORCING_SV
     &                         ad_sustr, ad_svstr)
#  else
#   ifdef SOLVE3D
     &                         Hz,                                      &
     &                         ad_t, ad_u, ad_v,                        &
#   else
     &                         ad_ubar, ad_vbar,                        &
#   endif
     &                         ad_zeta)
#  endif
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
      integer, intent(in) :: Nstr, Nend
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
#   if defined FORCING_SV
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
#   else
#    ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#    endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: ad_state(Nstr:)
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
#   if defined FORCING_SV
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
#   else
#    ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#    else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#    endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(out) :: ad_state(Nstr:Nend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

#  if defined FORCING_SV
      integer, dimension(2) :: offset
#  else
      integer, dimension(5+NT(ng)) :: offset
#  endif

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
      DO is=Nstr,Nend
        ad_state(is)=Aspv
      END DO
#  endif

#  if defined FORCING_SV
!
!-----------------------------------------------------------------------
!  Load adjoint FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#   ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#   endif
!
#   ifdef MASKING
      offset(isUstr)=0
      offset(isVstr)=offset(isUstr)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isUstr)=0
      offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isUstr)=0
      offset(isVstr)=offset(isUstr)+(Lm(ng)-1)*Mm(ng)
#    endif
#   endif

!
!  Load adjoint surface U-momentum stress.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#    endif
#   endif
      scale=1.0_r8
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUstr)
            ad_state(is)=scale*ad_sustr(i,j)
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isVstr)
          ad_state(is)=scale*ad_sustr(i,j)
#   endif
        END DO
      END DO
!
!  Load adjoint surface V-momentum stress.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#    endif
#   endif
      scale=1.0_r8
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVstr)
            ad_state(is)=scale*ad_svstr(i,j)
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
          ad_state(is)=scale*ad_svstr(i,j)
#   endif
        END DO
      END DO

#  else
!
!-----------------------------------------------------------------------
!  Load adjoint STATE variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#   ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#   endif
!
#   ifdef SOLVE3D
#    ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#    else
#     ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#     else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#     endif
#    endif
#   else
#    ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#    else
#     ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#     else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#     endif
#    endif
#   endif
!
!  Load adjoint free-surface.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#   ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            ad_state(is)=scale*ad_zeta(i,j,kstp)
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          ad_state(is)=scale*ad_zeta(i,j,kstp)
#   endif
        END DO
      END DO
#   ifndef SOLVE3D
!
!  Load adjoint 2D U-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#     else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#    ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#    endif
#    ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            ad_state(is)=scale*ad_ubar(i,j,kstp)
          END IF
#    else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          ad_state(is)=scale*ad_ubar(i,j,kstp)
#    endif
        END DO
      END DO
!
!  Load adjoint 2D V-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#     else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#    ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#    endif
#    ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            ad_state(is)=scale*ad_vbar(i,j,kstp)
          END IF
#    else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          ad_state(is)=scale*ad_vbar(i,j,kstp)
#    endif
        END DO
      END DO
#   else
!
!  Load adjoint 3D U-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#     else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO k=1,N(ng)
#    ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#    else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#    endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#    ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#     ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#     endif
              is=IJwaterU(i,j)+iadd
              ad_state(is)=scale*ad_u(i,j,k,nstp)
            END IF
#    else
#     ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#     endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_u(i,j,k,nstp)
#    endif
          END DO
        END DO
      END DO
!
!  Load adjoint 3D V-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#     else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO k=1,N(ng)
#    ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#    else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#    endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#    ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#     ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#     endif
              is=IJwaterV(i,j)+iadd
              ad_state(is)=scale*ad_v(i,j,k,nstp)
            END IF
#    else
#     ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#     endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_v(i,j,k,nstp)
#    endif
          END DO
        END DO
      END DO
!
!  Load adjoint tracers variables. For now, use salinity scale for
!  passive tracers.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#     else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#     endif
#    endif
      DO itrc=1,NT(ng)
#    ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#    else
        scale=1.0_r8
#    endif
        DO k=1,N(ng)
#    ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#    else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#    endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#    ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#     ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#     endif
                is=IJwaterR(i,j)+iadd
                ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
              END IF
#    else
#     ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#     endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
#    endif
            END DO
          END DO
        END DO
      END DO
#   endif
#  endif
      RETURN
      END SUBROUTINE ad_pack_tile

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
      CALL wclock_on (ng, iADM, 1)
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
      CALL wclock_off (ng, iADM, 1)
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
     &                           Nstr, Nend, state,                     &
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
      integer, intent(in) :: Nstr, Nend
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
      real(r8), intent(in) :: state(Nstr:)
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
      real(r8), intent(in) :: state(Nstr:Nend)
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
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract adjoint state variables from full 1D state vector.
!-----------------------------------------------------------------------
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
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
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
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
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
      Imax=Lm(ng)-1
      Ioff=1
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
      Joff=2
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
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
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
              ad_u(i,j,k,nstp)=scale*state(is)
            ELSE
              ad_u(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
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
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
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
              ad_v(i,j,k,nstp)=scale*state(is)
            ELSE
              ad_v(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
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
                ad_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                ad_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
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
      CALL wclock_on (ng, iTLM, 1)
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
      CALL wclock_off (ng, iTLM, 1)
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
     &                         Nstr, Nend, tl_state,                    &
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
      integer, intent(in) :: Nstr, Nend
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

      real(r8), intent(out) :: tl_state(Nstr:)
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

      real(r8), intent(out) :: tl_state(Nstr:Nend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
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
      DO is=Nstr,Nend
        tl_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load tangent linear state variables into full 1D state vector.
!-----------------------------------------------------------------------
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
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
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
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
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
      Imax=Lm(ng)-1
      Ioff=1
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
      Joff=2
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
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
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
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
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
#  if defined FORCING_SV
      USE mod_forces
#  else
      USE mod_ocean
#  endif
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
      CALL wclock_on (ng, iTLM, 1)
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
#  if defined FORCING_SV
     &                     FORCES(ng) % tl_sustr,                       &
     &                     FORCES(ng) % tl_svstr)
#  else
#   ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#   else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#   endif
     &                     OCEAN(ng) % tl_zeta)
#  endif
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 1)
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
     &                           Nstr, Nend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  if defined FORCING_SV
     &                           tl_sustr, tl_svstr)
#  else
#   ifdef SOLVE3D
     &                           Hz,                                    &
     &                           tl_t, tl_u, tl_v,                      &
#   else
     &                           tl_ubar, tl_vbar,                      &
#   endif
     &                           tl_zeta)
#  endif
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
      integer, intent(in) :: Nstr, Nend
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
      real(r8), intent(in) :: state(Nstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   if defined FORCING_SV
      real(r8), intent(inout) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_svstr(LBi:,LBj:)
#   else
#    ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#    endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#   endif
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Nstr:Nend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   if defined FORCING_SV
      real(r8), intent(inout) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_svstr(LBi:UBi,LBj:UBj)
#   else
#    ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#    else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#    endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   endif
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

#  if defined FORCING_SV
      integer, dimension(2) :: offset
#  else
      integer, dimension(5+NT(ng)) :: offset
#  endif
      real(r8) :: cff, scale

#  include "set_bounds.h"

#  if defined FORCING_SV
!
!-----------------------------------------------------------------------
!  Extract tangent linear FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#   ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#   endif
!
#   ifdef MASKING
      offset(isUstr)=0
      offset(isVstr)=offset(isUstr)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isUstr)=0
      offset(isVstr)=offset(isUstr)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isUstr)=0
      offset(isVstr)=offset(isUstr)+(Lm(ng)-1)*Mm(ng)
#    endif
#   endif
!
!  Extract tangent linear surface U-momentum stress.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#    endif
#   endif
      scale=1.0_r8
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUstr)
            tl_sustr(i,j)=scale*state(is)
          ELSE
            tl_sustr(i,j)=0.0_r8
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
          tl_sustr(i,j)=scale*state(is)
#   endif
        END DO
      END DO
!
!  Extract tangent linear surface V-momentum stress.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#    endif
#   endif
      scale=1.0_r8
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVstr)
            tl_svstr(i,j)=scale*state(is)
          ELSE
            tl_svstr(i,j)=0.0_r8
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
          tl_svstr(i,j)=scale*state(is)
#   endif
        END DO
      END DO

#  else
!
!-----------------------------------------------------------------------
!  Extract tangent linear STATE variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#   ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#   endif
!
#   ifdef SOLVE3D
#    ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#    else
#     ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#     else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#     endif
#    endif
#   else
#    ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#    else
#     ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#     else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#     endif
#    endif
#   endif
!
!  Extract tangent linear free-surface.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#   ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            tl_zeta(i,j,kstp)=scale*state(is)
          ELSE
            tl_zeta(i,j,kstp)=0.0_r8
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          tl_zeta(i,j,kstp)=scale*state(is)
#   endif
        END DO
      END DO
#   ifndef SOLVE3D
!
!  Extract tangent linear 2D U-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#     else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#    ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#    endif
#    ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            tl_ubar(i,j,kstp)=scale*state(is)
          ELSE
            tl_ubar(i,j,kstp)=0.0_r8
          END IF
#    else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          tl_ubar(i,j,kstp)=scale*state(is)
#    endif
        END DO
      END DO
!
!  Extract tangent linear 2D V-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#     else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#    ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#    endif
#    ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            tl_vbar(i,j,kstp)=scale*state(is)
          ELSE
            tl_vbar(i,j,kstp)=0.0_r8
          END IF
#    else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          tl_vbar(i,j,kstp)=scale*state(is)
#    endif
        END DO
      END DO
#   else
!
!  Extract tangent linear 3D U-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#     else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO k=1,N(ng)
#    ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#    else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#    endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#    ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#     ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#     endif
              is=IJwaterU(i,j)+iadd
              tl_u(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_u(i,j,k,nstp)=0.0_r8
            END IF
#    else
#     ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#     endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            tl_u(i,j,k,nstp)=scale*state(is)
#    endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear 3D V-velocity.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#     else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#     endif
#    endif
#    ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#    else
      scale=1.0_r8
#    endif
      DO k=1,N(ng)
#    ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#    else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#    endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#    ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#     ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#     endif
              is=IJwaterV(i,j)+iadd
              tl_v(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_v(i,j,k,nstp)=0.0_r8
            END IF
#    else
#     ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#     endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            tl_v(i,j,k,nstp)=scale*state(is)
#    endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#    ifndef MASKING
#     ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#     else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#     endif
#    endif
      DO itrc=1,NT(ng)
#    ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#    else
        scale=1.0_r8
#    endif
        DO k=1,N(ng)
#    ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#    else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#    endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#    ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#     ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#     endif
                is=IJwaterR(i,j)+iadd
                tl_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                tl_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#    else
#     ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#     endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_t(i,j,k,nstp,itrc)=scale*state(is)
#    endif
            END DO
          END DO
        END DO
      END DO
#   endif
#  endif
      RETURN
      END SUBROUTINE tl_unpack_tile
# endif

# if defined SO_SEMI
      SUBROUTINE so_pack (ng, tile, Mstr, Mend, Irec)
!
!=======================================================================
!                                                                      !
!  This routine packs the  adjoint surface forcing  variables used     !
!  in stochastic optimals. If land/sea masking, the forcing vector     !
!  contains only water points.  The  stochastic optimal  vector is     !
!  assumed to be  nondimensional to facilitate its interpretation.     !
!  They are  non-dimensionalized  by a representative value of the     !
!  surface forcing standard deviation, SO_sigma.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_storage
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
#  endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend, Irec
!
!  Local variable declarations.
!
#  ifdef DISTRIBUTE
      integer :: i
      real(r8), parameter :: Aspv = 0.0_r8
#  endif
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 1)
#  endif
      CALL so_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork,                          &
#  else
     &                   Mstr, Mend, so_state(Mstr:,Irec),              &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
#  ifdef SOLVE3D
     &                   FORCES(ng) % ad_stflx,                         &
#  endif
     &                   FORCES(ng) % ad_sustr,                         &
     &                   FORCES(ng) % ad_svstr)

#  ifdef DISTRIBUTE
!
!  Collect data from all nodes. Then, load to threaded stochastic
!  optimals state array.
!
      CALL mp_collect (ng, iADM, Mstate(ng), Aspv, Swork)
      DO i=Mstr,Mend
        so_state(i,Irec)=Swork(i)
      END DO
#  endif
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 1)
#  endif
      RETURN
      END SUBROUTINE so_pack
!
!***********************************************************************
      SUBROUTINE so_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Nstr, Nend,                              &
     &                         so_state,                                &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
#  ifdef SOLVE3D
     &                         ad_stflx,                                &
#  endif
     &                         ad_sustr, ad_svstr)

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
      integer, intent(in) :: Nstr, Nend
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
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
      real(r8), intent(out) :: so_state(Nstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: so_state(Nstr:Nend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

#  ifdef SOLVE3D
      integer, dimension(2+NT(ng)) :: offset
#  else
      integer, dimension(2) :: offset
#  endif

      real(r8) :: cff

      real(r8), parameter :: Aspv = 0.0_r8

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Nstr,Nend
        so_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
      iadd=0
      offset=0
#  ifdef MASKING
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        offset(isUstr)=iadd
        iadd=NwaterU(ng)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        offset(isVstr)=offset(isVstr-1)+iadd
        iadd=NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=offset(isTsur(itrc)-1)+iadd
          iadd=NwaterR(ng)
        END IF
      END DO
#   endif
#  else
#   ifdef FULL_GRID
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        offset(isUstr)=iadd
        iadd=(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        offset(isVstr)=offset(isVstr-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+1)
      END IF
#    ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)
        END IF
      END DO
#    endif
#   else
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        offset(isUstr)=iadd
        iadd=(Lm(ng)-1)*Mm(ng)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        offset(isVstr)=offset(isVstr-1)+iadd
        iadd=(Lm(ng))*(Mm(ng)-1)
      END IF
#    ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=offset(isTsur(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)
        END IF
      END DO
#    endif
#   endif
#  endif
!
!  Load scaled adjoint surface U-momentum stress.
!
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#   endif
#  endif
        cff=SO_sdev(isUstr,ng)/rho0
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              so_state(is)=cff*ad_sustr(i,j)
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            so_state(is)=cff*ad_sustr(i,j)
#  endif
          END DO
        END DO
      END IF
!
!  Load scaled adjoint surface V-momentum stress.
!
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#   endif
#  endif
        cff=SO_sdev(isVstr,ng)/rho0
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              so_state(is)=cff*ad_svstr(i,j)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            so_state(is)=cff*ad_svstr(i,j)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Load scaled adjoint surface tracer fluxes.
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
        IF (itrc.eq.1) THEN
          cff=SO_sdev(isTsur(itrc),ng)/(rho0*Cp)
        ELSE
          cff=SO_sdev(isTsur(itrc),ng)
        END IF
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          iadd=offset(isTsur(itrc))
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+iadd
                so_state(is)=cff*ad_stflx(i,j,itrc)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              so_state(is)=cff*ad_stflx(i,j,itrc)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE so_pack_tile

      SUBROUTINE so_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the adjoint surface forcing variables  used    !
!  stochastic optimals.  The forcing vector  contains only interior    !
!  water points. For output purposes, the surface forcing variables    !
!  are saved in the nonlinear state variables.                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_forces
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect
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
      integer :: i
      real(r8), parameter :: Aspv = 0.0_r8

#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 1)
#  endif
#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      DO i=1,Mstate(ng)
        Swork(i)=Aspv
      END DO
      DO i=Mstr,Mend
        Swork(i)=state(i)
      END DO
      CALL mp_collect (ng, iADM, Mstate(ng), Aspv, Swork)
#  endif
      CALL so_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
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
#  ifdef SOLVE3D
     &                     FORCES(ng) % stflx,                          &
#  endif
     &                     FORCES(ng) % sustr,                          &
     &                     FORCES(ng) % svstr)

#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 1)
#  endif
      RETURN
      END SUBROUTINE so_unpack
!
!***********************************************************************
      SUBROUTINE so_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           Nstr, Nend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
#  ifdef SOLVE3D
     &                           stflx,                                 &
#  endif
     &                           sustr, svstr)
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
      integer, intent(in) :: Nstr, Nend
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
      real(r8), intent(in) :: state(Nstr:)
#   ifdef SOLVE3D
      real(r8), intent(out) :: stflx(LBi:,LBj:,:)
#   endif
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Nstr:Nend)

#   ifdef SOLVE3D
      real(r8), intent(out) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#   endif
      real(r8), intent(out) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: svstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

#  ifdef SOLVE3D
      integer, dimension(2+NT(ng)) :: offset
#  else
      integer, dimension(2) :: offset
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract adjoint FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
      iadd=0
      offset=0
#  ifdef MASKING
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        offset(isUstr)=iadd
        iadd=NwaterU(ng)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        offset(isVstr)=offset(isVstr-1)+iadd
        iadd=NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=offset(isTsur(itrc)-1)+iadd
          iadd=NwaterR(ng)
        END IF
      END DO
#   endif
#  else
#   ifdef FULL_GRID
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        offset(isUstr)=iadd
        iadd=(Lm(ng)+1)*(Mm(ng)+2)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        offset(isVstr)=offset(isVstr-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+1)
      END IF
#    ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
          iadd=(Lm(ng)+2)*(Mm(ng)+2)
        END IF
      END DO
#    endif
#   else
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
        offset(isUstr)=iadd
        iadd=(Lm(ng)-1)*Mm(ng)
      END IF
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
        offset(isVstr)=offset(isVstr-1)+iadd
        iadd=(Lm(ng))*(Mm(ng)-1)
      END IF
#    ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=offset(isTsur(itrc)-1)+iadd
          iadd=Lm(ng)*Mm(ng)
        END IF
      END DO
#    endif
#   endif
#  endif
!
!  Extract adjoint surface U-momentum stress.
!
      IF (SCALARS(ng)%SOstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#   endif
#  endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              sustr(i,j)=state(is)
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            sustr(i,j)=state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Extract adjoint surface V-momentum stress.
!
      IF (SCALARS(ng)%SOstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#   endif
#  endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              svstr(i,j)=state(is)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            svstr(i,j)=state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Extract adjoint surface tracer fluxes.
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
        IF (SCALARS(ng)%SOstate(isTsur(itrc))) THEN
          iadd=offset(isTsur(itrc))
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+iadd
                stflx(i,j,itrc)=state(is)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              stflx(i,j,itrc)=state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE so_unpack_tile

#  ifdef SO_SEMI_WHITE
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
      DO rec=1,Nsemi
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
          my_SOnorm=my_SOnorm+so_state(is,rec)*state(is)
        END DO
!
        IF (rec.ne.Nsemi) THEN
          DO is=Mstr,Mend
            my_SOnorm1=my_SOnorm1+so_state(is,rec+1)*state(is)
            my_TRnorm=my_TRnorm+                                        &
     &                cff1*so_state(is,rec)*so_state(is,rec)+           &
     &                2.0_r8*cff2*so_state(is,rec)*so_state(is,rec+1)+  &
     &                cff*so_state(is,rec+1)*so_state(is,rec+1)
          END DO
        ELSE
          DO is=Mstr,Mend
            my_TRnorm=my_TRnorm+so_state(is,rec)*so_state(is,rec)
          END DO
        END IF
!
!  Global reduction of normalization factor.
!
        IF (SOUTH_WEST_CORNER.and.                                      &
     &      NORTH_EAST_CORNER) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
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
        IF (rec.ne.Nsemi) THEN
          DO is=Mstr,Mend
            ad_state(is)=ad_state(is)+                                  &
     &                   cff1*SOnorm *so_state(is,rec  )+               &
     &                   cff2*SOnorm1*so_state(is,rec  )+               &
     &                   cff2*SOnorm *so_state(is,rec+1)+               &
     &                   cff *SOnorm1*so_state(is,rec+1)
          END DO
        ELSE
          DO is=Mstr,Mend
            ad_state(is)=ad_state(is)+SOnorm*so_state(is,rec)
          END DO
        END IF
      END DO
!
!  Global reduction of normalization factor, TRnorm.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
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

      real(r8), dimension(Nsemi) :: Bcoef
      real(r8), dimension(Nsemi) :: SOdotprod
      real(r8), dimension(Nsemi) :: my_dotprod

#   ifdef DISTRIBUTE
      character (len=3), dimension(Nsemi) :: op_handle
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
      DO rec=1,Nsemi
        my_dotprod(rec)=0.0_r8
        DO is=Mstr,Mend
          my_dotprod(rec)=my_dotprod(rec)+so_state(is,rec)*state(is)
        END DO
      END DO
!
!  Global reduction of dot products.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                         ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
      END IF
!$OMP CRITICAL (SO_DOT)
      IF (tile_count.eq.0) THEN
        DO rec=1,Nsemi
          SOdotprod(rec)=0.0_r8
        END DO
      END IF
      DO rec=1,Nsemi
        SOdotprod(rec)=SOdotprod(rec)+my_dotprod(rec)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        DO rec=1,Nsemi
          op_handle(rec)='SUM'
        END DO
        CALL mp_reduce (ng, iADM, Nsemi, SOdotprod, op_handle)
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
      DO rec=1,Nsemi
        ntAD=(rec-1)*nADJ(ng)+1
        SOnorm=0.0_r8
        DO rec1=1,Nsemi
          ntTL=(rec1-1)*nADJ(ng)+1
          CALL sp_bcoef (ng, ntAD, ntTL, Bcoef(rec1))
          SOnorm=SOnorm+Bcoef(rec1)*SOdotprod(rec1)
          DO is=Mstr,Mend
            my_TRnorm=my_TRnorm+                                        &
     &                so_state(is,rec)*Bcoef(rec1)*so_state(is,rec1)
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
          ad_state(is)=ad_state(is)+SOnorm*so_state(is,rec)
        END DO
      END DO
!
!  Global reduction of normalization factor, TRnorm.
!
      IF (SOUTH_WEST_CORNER.and.                                        &
     &    NORTH_EAST_CORNER) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
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
# ifdef SO_NON_AR1
!
!  Use a user-defined temporal decorrelation function such as in the
!  form of a look-up table computed from actual data.
!
      sp_autoc=0.0_r8
# else
!
!  Assume an AR(1) process with decorrelation time SO_decay.
!
      sp_autoc=EXP(-ABS(REAL(idf-1,r8))*dt(ng)/SO_decay(ng))
# endif
!
      RETURN
      END FUNCTION sp_autoc

# undef IR_RANGE
# undef IU_RANGE
# undef JR_RANGE
# undef JV_RANGE
#endif

      END MODULE packing_mod
