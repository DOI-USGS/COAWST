#include "cppdefs.h"
      MODULE dotproduct_mod

#if defined TLM_CHECK || defined PROPAGATOR
!
!svn $Id: dotproduct.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines compute various dot products between nonlinear,      !
!  tangent linear, and adjoint state variables.                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef TLM_CHECK
      PUBLIC :: ad_dotproduct
# endif
# ifdef TLM_CHECK
      PUBLIC :: nl_dotproduct, tl_dotproduct
# endif
# ifdef PROPAGATOR
#  ifdef ADJOINT
      PUBLIC :: ad_statenorm
#  endif
#  ifdef TANGENT
      PUBLIC :: tl_statenorm
#  endif
# endif

      CONTAINS

# ifdef TLM_CHECK
      SUBROUTINE nl_dotproduct (ng, tile, Linp)
!
!=======================================================================
!                                                                      !
!  This routine computes dot product function (g1) associated with     !
!  the nonlinear solution.                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
#  endif
      CALL nl_dotproduct_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         KOUT, Linp,                              &
#  ifdef SOLVE3D
     &                         NOUT,                                    &
#  endif
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
#  endif
#  ifdef SOLVE3D
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % t,                           &
     &                         OCEAN(ng) % tG,                          &
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % u,                           &
     &                         OCEAN(ng) % uG,                          &
     &                         OCEAN(ng) % ad_v,                        &
     &                         OCEAN(ng) % v,                           &
     &                         OCEAN(ng) % vG,                          &
#  endif
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ubar,                        &
     &                         OCEAN(ng) % ubarG,                       &
     &                         OCEAN(ng) % ad_vbar,                     &
     &                         OCEAN(ng) % vbar,                        &
     &                         OCEAN(ng) % vbarG,                       &
     &                         OCEAN(ng) % ad_zeta,                     &
     &                         OCEAN(ng) % zeta,                        &
     &                         OCEAN(ng) % zetaG)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
#  endif
      RETURN
      END SUBROUTINE nl_dotproduct
!
!***********************************************************************
      SUBROUTINE nl_dotproduct_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Kinp, Linp,                        &
#  ifdef SOLVE3D
     &                               Ninp,                              &
#  endif
#  ifdef MASKING
     &                               rmask, umask, vmask,               &
#  endif
#  ifdef SOLVE3D
     &                               ad_t, t, tG,                       &
     &                               ad_u, u, uG,                       &
     &                               ad_v, v, vG,                       &
#  endif
     &                               ad_ubar, ubar, ubarG,              &
     &                               ad_vbar, vbar, vbarG,              &
     &                               ad_zeta, zeta, zetaG)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_ncparam
      USE mod_iounits
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Linp
#  ifdef SOLVE3D
      integer, intent(in) :: Ninp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: t (LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tG(LBi:,LBj:,:,:,:)

      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: u (LBi:,LBj:,:,:)
      real(r8), intent(in) :: uG(LBi:,LBj:,:,:)

      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v (LBi:,LBj:,:,:)
      real(r8), intent(in) :: vG(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ubar (LBi:,LBj:,:)
      real(r8), intent(in) :: ubarG(LBi:,LBj:,:)

      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar (LBi:,LBj:,:)
      real(r8), intent(in) :: vbarG(LBi:,LBj:,:)

      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: zeta (LBi:,LBj:,:)
      real(r8), intent(in) :: zetaG(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: t (LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng))

      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: u (LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: uG(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v (LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: vG(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ubar (LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ubarG(LBi:UBi,LBj:UBj,2)

      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar (LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbarG(LBi:UBi,LBj:UBj,2)

      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zeta (LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: zetaG(LBi:UBi,LBj:UBj,2)
#  endif
!
!  Local variable declarations.
!
      integer :: NSUB, Tindex, i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: my_dot, p, scale

#  ifdef DISTRIBUTE
      character(len=3 ) :: op_handle
#  endif
      character(len=12) :: text

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variables: steepest descent
!  direction times the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
      IF ((outer.ge.1).and.(MOD(iic(ng)-1,nHIS(ng)).eq.0).and.          &
     &    ((nrrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
!
!  Set perturbation amplitude as a function of the inner loop.
!
        p=10.0_r8**REAL(-inner,r8)
        scale=1.0_r8/SQRT(adDotProduct)
        my_dot=0.0_r8
        IF (time(ng).eq.Tintrp(1,idFsur,ng)) THEN
          Tindex=1
        ELSE IF (time(ng).eq.Tintrp(2,idFsur,ng)) THEN
          Tindex=2
        END IF
        IF (iic(ng).eq.ntstart(ng)) THEN
          text='------------'
        ELSE
          text='            '
        END IF
!
!  Free-surface.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            my_dot=my_dot+                                              &
     &             p*scale*ad_zeta(i,j,Linp)*                           &
     &             (zeta(i,j,Kinp)-zetaG(i,j,Tindex))
#  ifdef MASKING
            my_dot=my_dot*rmask(i,j)
#  endif
          END DO
        END DO
!
!  2D u-momentum component.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            my_dot=my_dot+                                              &
     &             p*scale*ad_ubar(i,j,Linp)*                           &
     &             (ubar(i,j,Kinp)-ubarG(i,j,Tindex))
#  ifdef MASKING
            my_dot=my_dot*umask(i,j)
#  endif
          END DO
        END DO
!
!  2D v-momentum component.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            my_dot=my_dot+                                              &
     &             p*scale*ad_vbar(i,j,Linp)*                           &
     &             (vbar(i,j,Kinp)-vbarG(i,j,Tindex))
#  ifdef MASKING
            my_dot=my_dot*vmask(i,j)
#  endif
          END DO
        END DO
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              my_dot=my_dot+                                            &
     &               p*scale*ad_u(i,j,k,Linp)*                          &
     &               (u(i,j,k,Ninp)-uG(i,j,k,Tindex))
#   ifdef MASKING
              my_dot=my_dot*umask(i,j)
#   endif
            END DO
          END DO
        END DO
!
!  3D v-momentum component.
!
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              my_dot=my_dot+                                            &
     &               p*scale*ad_v(i,j,k,Linp)*                          &
     &               (v(i,j,k,Ninp)-vG(i,j,k,Tindex))
#   ifdef MASKING
              my_dot=my_dot*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
!
!  Tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                my_dot=my_dot+                                          &
     &                 p*scale*ad_t(i,j,k,Linp,itrc)*                   &
     &                 (t(i,j,k,Ninp,itrc)-tG(i,j,k,Tindex,itrc))
#   ifdef MASKING
                my_dot=my_dot*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END DO
#  endif
!
!  Perform parallel global reduction operation: dot product.
!
#  ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
#  else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
#  endif
!$OMP CRITICAL (NL_DOT)
        IF (tile_count.eq.0) THEN
          DotProduct=my_dot
        ELSE
          DotProduct=DotProduct+my_dot
        END IF
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#  ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, iNLM, 1, DotProduct, op_handle)
#  endif
          ig1count=ig1count+1
          g1(ig1count)=DotProduct
          IF (Master) THEN
            WRITE (stdout,10) outer, inner, p, text, Tindex, DotProduct
 10         FORMAT (5x,'(Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,     &
     &              'TLM Check, p  = ',1p,e21.14,0p,/,5x,a,4x,          &
     &              '(Index=',i1,')',5x,'TLM Check, g1 = ',1p,e21.14,0p)
          END IF
        END IF
!$OMP END CRITICAL (NL_DOT)
      END IF

      RETURN
      END SUBROUTINE nl_dotproduct_tile
# endif

# if defined ADJOINT && defined PROPAGATOR
      SUBROUTINE ad_statenorm (ng, tile, Kinp, Ninp, StateNorm)
!
!=======================================================================
!                                                                      !
!  This routine computes the dot product norm of requested adjoint     !
!  state.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
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
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      CALL ad_statenorm_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        Kinp, Ninp,                               &
     &                        StateNorm,                                &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
#  endif
#  ifdef SOLVE3D
     &                        GRID(ng) % Hz,                            &
#  else
     &                        GRID(ng) % h,                             &
#  endif
#  ifdef SOLVE3D
     &                        OCEAN(ng) % ad_t,                         &
     &                        OCEAN(ng) % ad_u,                         &
     &                        OCEAN(ng) % ad_v,                         &
#  else
     &                        OCEAN(ng) % ad_ubar,                      &
     &                        OCEAN(ng) % ad_vbar,                      &
#  endif
     &                        OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif
      RETURN
      END SUBROUTINE ad_statenorm
!
!***********************************************************************
      SUBROUTINE ad_statenorm_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Kinp, Ninp,                         &
     &                              StateNorm,                          &
#  ifdef MASKING
     &                              rmask, umask, vmask,                &
#  endif
#  ifdef SOLVE3D
     &                              Hz,                                 &
#  else
     &                              h,                                  &
#  endif
#  ifdef SOLVE3D
     &                              ad_t, ad_u, ad_v,                   &
#  else
     &                              ad_ubar, ad_vbar,                   &
#  endif
     &                              ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Ninp

      real(r8), intent(out) :: StateNorm
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   else
      real(r8), intent(in) :: h(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:Ubi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: cff, cff1, my_norm, scale

#  ifdef DISTRIBUTE
      character(len=3) :: op_handle
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute dot product norm of current adjoint solution.
!-----------------------------------------------------------------------

#  ifdef FULL_GRID
#   define IR_RANGE IstrT,IendT
#   define IU_RANGE IstrP,IendT
#   define JR_RANGE JstrT,JendT
#   define JV_RANGE JstrP,JendT
#  else
#   define IR_RANGE Istr,Iend
#   define IU_RANGE IstrU,Iend
#   define JR_RANGE Jstr,Jend
#   define JV_RANGE JstrV,Jend
#  endif
!
!  Initialize private norm.
!
      my_norm=0.0_r8
!
!  Free-surface.
!
      scale=0.5_r8*g*rho0
      DO j=JR_RANGE
        DO i=IR_RANGE
          cff1=scale*ad_zeta(i,j,Kinp)*ad_zeta(i,j,Kinp)
#  ifdef MASKING
          cff1=cff1*rmask(i,j)
#  endif
          my_norm=my_norm+cff1
        END DO
      END DO

#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      cff=0.25_r8*rho0
      DO j=JR_RANGE
        DO i=IU_RANGE
          scale=cff*(h(i-1,j)+h(i,j))
          cff1=scale*ad_ubar(i,j,Kinp)*ad_ubar(i,j,Kinp)
#   ifdef MASKING
          cff1=cff1*umask(i,j)
#   endif
          my_norm=my_norm+cff1
        END DO
      END DO
!
!  2D v-momentum component.
!
      cff=0.25_r8*rho0
      DO j=JV_RANGE
        DO i=IR_RANGE
          scale=cff*(h(i,j-1)+h(i,j))
          cff1=scale*ad_vbar(i,j,Kinp)*ad_vbar(i,j,Kinp)
#   ifdef MASKING
          cff1=cff1*vmask(i,j)
#   endif
          my_norm=my_norm+cff1
        END DO
      END DO
#  else
!
!  3D u-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=JR_RANGE
          DO i=IU_RANGE
            scale=cff*(Hz(i-1,j,k)+Hz(i,j,k))
            cff1=scale*ad_u(i,j,k,Ninp)*ad_u(i,j,k,Ninp)
#   ifdef MASKING
            cff1=cff1*umask(i,j)
#   endif
            my_norm=my_norm+cff1
          END DO
        END DO
      END DO
!
!  3D v-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=JV_RANGE
          DO i=IR_RANGE
            scale=cff*(Hz(i,j-1,k)+Hz(i,j,k))
            cff1=scale*ad_v(i,j,k,Ninp)*ad_v(i,j,k,Ninp)
#   ifdef MASKING
            cff1=cff1*vmask(i,j)
#   endif
            my_norm=my_norm+cff1
          END DO
        END DO
      END DO
!
!  Tracers. For now, use salinity scale for passive tracers.
!
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
        DO k=1,N(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
              scale=cff*Hz(i,j,k)
              cff1=scale*ad_t(i,j,k,Ninp,itrc)*                         &
     &                   ad_t(i,j,k,Ninp,itrc)
#   ifdef MASKING
              cff1=cff1*rmask(i,j)
#   endif
              my_norm=my_norm+cff1
            END DO
          END DO
        END DO
      END DO
#  endif
!
!  Perform parallel global reduction operation: dot product.
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
!$OMP CRITICAL (TL_NORM)
      IF (tile_count.eq.0) THEN
        StateNorm=my_norm
      ELSE
        StateNorm=StateNorm+my_norm
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, iTLM, 1, StateNorm, op_handle)
#  endif
      END IF
!$OMP END CRITICAL (TL_NORM)

#  undef IR_RANGE
#  undef IU_RANGE
#  undef JR_RANGE
#  undef JV_RANGE

      RETURN
      END SUBROUTINE ad_statenorm_tile
# endif

# if defined TANGENT && defined PROPAGATOR
      SUBROUTINE tl_statenorm (ng, tile, Kinp, Ninp, StateNorm)
!
!=======================================================================
!                                                                      !
!  This routine computes the dot product norm of requested tangent     !
!  linear state.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
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
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL tl_statenorm_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        Kinp, Ninp,                               &
     &                        StateNorm,                                &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
#  endif
#  ifdef BNORM
#    ifdef SOLVE3D
     &                        OCEAN(ng) % e_t,                          &
     &                        OCEAN(ng) % e_u,                          &
     &                        OCEAN(ng) % e_v,                          &
#    else
     &                        OCEAN(ng) % e_ubar,                       &
     &                        OCEAN(ng) % e_vbar,                       &
#    endif
     &                        OCEAN(ng) % e_zeta,                       &
#  endif
#  ifdef SOLVE3D
     &                        GRID(ng) % Hz,                            &
#  else
     &                        GRID(ng) % h,                             &
#  endif
#  ifdef SOLVE3D
     &                        OCEAN(ng) % tl_t,                         &
     &                        OCEAN(ng) % tl_u,                         &
     &                        OCEAN(ng) % tl_v,                         &
#  else
     &                        OCEAN(ng) % tl_ubar,                      &
     &                        OCEAN(ng) % tl_vbar,                      &
#  endif
     &                        OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif
      RETURN
      END SUBROUTINE tl_statenorm

!
!***********************************************************************
      SUBROUTINE tl_statenorm_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              Kinp, Ninp,                         &
     &                              StateNorm,                          &
#  ifdef MASKING
     &                              rmask, umask, vmask,                &
#  endif
#  ifdef BNORM
#    ifdef SOLVE3D
     &                              e_t, e_u, e_v,                      &
#    else
     &                              e_ubar, e_vbar,                     &
#    endif
     &                              e_zeta,                             &
#  endif
#  ifdef SOLVE3D
     &                              Hz,                                 &
#  else
     &                              h,                                  &
#  endif
#  ifdef SOLVE3D
     &                              tl_t, tl_u, tl_v,                   &
#  else
     &                              tl_ubar, tl_vbar,                   &
#  endif
     &                              tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Ninp

      real(r8), intent(out) :: StateNorm
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   else
      real(r8), intent(in) :: h(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#   ifdef BNORM
      real(r8), intent(in) :: e_zeta(LBi:,LBj:,:)
#    ifdef SOLVE3D
      real(r8), intent(in) :: e_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: e_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: e_v(LBi:,LBj:,:,:)
#    else
      real(r8), intent(in) :: e_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: e_vbar(LBi:,LBj:,:)
#    endif
#   endif
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:Ubi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   else
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#   ifdef BNORM
      real(r8), intent(in) :: e_zeta(LBi:UBi,LBj:UBj,NSA)
#    ifdef SOLVE3D
      real(r8), intent(in) :: e_t(LBi:UBi,LBj:UBj,N(ng),NSA,NT(ng))
      real(r8), intent(in) :: e_u(LBi:UBi,LBj:UBj,N(ng),NSA)
      real(r8), intent(in) :: e_v(LBi:UBi,LBj:UBj,N(ng),NSA)
#    else
      real(r8), intent(in) :: e_ubar(LBi:UBi,LBj:UBj,NSA)
      real(r8), intent(in) :: e_vbar(LBi:UBi,LBj:UBj,NSA)
#    endif
#   endif
#  endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: cff, cff1, my_norm, scale

#  ifdef DISTRIBUTE
      character(len=3) :: op_handle
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute dot product norm of current tangent linear solution.
!-----------------------------------------------------------------------

#  ifdef FULL_GRID
#   define IR_RANGE IstrT,IendT
#   define IU_RANGE IstrP,IendT
#   define JR_RANGE JstrT,JendT
#   define JV_RANGE JstrP,JendT
#  else
#   define IR_RANGE Istr,Iend
#   define IU_RANGE IstrU,Iend
#   define JR_RANGE Jstr,Jend
#   define JV_RANGE JstrV,Jend
#  endif
!
!  Initialize private norm.
!
      my_norm=0.0_r8
!
!  Free-surface.
!
      scale=0.5_r8*g*rho0
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef BNORM
          IF (e_zeta(i,j,1).gt.0.0_r8) THEN
            cff1=1.0_r8/(e_zeta(i,j,1)*e_zeta(i,j,1))
          ELSE
            cff1=1.0_r8
          END IF
          cff1=cff1*tl_zeta(i,j,Kinp)*tl_zeta(i,j,Kinp)
#  else
          cff1=scale*tl_zeta(i,j,Kinp)*tl_zeta(i,j,Kinp)
#  endif
#  ifdef MASKING
          cff1=cff1*rmask(i,j)
#  endif
          my_norm=my_norm+cff1
        END DO
      END DO

#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      cff=0.25_r8*rho0
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef BNORM
          IF (e_ubar(i,j,1).gt.0.0_r8) THEN
            cff1=1.0_r8/(e_ubar(i,j,1)*e_ubar(i,j,1))
          ELSE
            cff1=1.0_r8
          END IF
          cff1=cff1*tl_ubar(i,j,Kinp)*tl_ubar(i,j,Kinp)
#   else
          scale=cff*(h(i-1,j)+h(i,j))
          cff1=scale*tl_ubar(i,j,Kinp)*tl_ubar(i,j,Kinp)
#   endif
#   ifdef MASKING
          cff1=cff1*umask(i,j)
#   endif
          my_norm=my_norm+cff1
        END DO
      END DO
!
!  2D v-momentum component.
!
      cff=0.25_r8*rho0
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef BNORM
          IF (e_vbar(i,j,1).gt.0.0_r8) THEN
            cff1=1.0_r8/(e_vbar(i,j,1)*e_vbar(i,j,1))
          ELSE
            cff1=1.0_r8
          END IF
          cff1=cff1*tl_vbar(i,j,Kinp)*tl_vbar(i,j,Kinp)
#   else
          scale=cff*(h(i,j-1)+h(i,j))
          cff1=scale*tl_vbar(i,j,Kinp)*tl_vbar(i,j,Kinp)
#   endif
#   ifdef MASKING
          cff1=cff1*vmask(i,j)
#   endif
          my_norm=my_norm+cff1
        END DO
      END DO
#  else
!
!  3D u-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef BNORM
            IF (e_u(i,j,k,1).gt.0.0_r8) THEN
              cff1=1.0_r8/(e_u(i,j,k,1)*e_u(i,j,k,1))
            ELSE
              cff1=1.0_r8
            END IF
            cff1=cff1*tl_u(i,j,k,Ninp)*tl_u(i,j,k,Ninp)
#   else
            scale=cff*(Hz(i-1,j,k)+Hz(i,j,k))
            cff1=scale*tl_u(i,j,k,Ninp)*tl_u(i,j,k,Ninp)
#   endif
#   ifdef MASKING
            cff1=cff1*umask(i,j)
#   endif
            my_norm=my_norm+cff1
          END DO
        END DO
      END DO
!
!  3D v-momentum component.
!
      cff=0.25_r8*rho0
      DO k=1,N(ng)
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef BNORM
            IF (e_v(i,j,k,1).gt.0.0_r8) THEN
              cff1=1.0_r8/(e_v(i,j,k,1)*e_v(i,j,k,1))
            ELSE
              cff1=1.0_r8
            END IF
            cff1=cff1*tl_v(i,j,k,Ninp)*tl_v(i,j,k,Ninp)
#   else
            scale=cff*(Hz(i,j-1,k)+Hz(i,j,k))
            cff1=scale*tl_v(i,j,k,Ninp)*tl_v(i,j,k,Ninp)
#   endif
#   ifdef MASKING
            cff1=cff1*vmask(i,j)
#   endif
            my_norm=my_norm+cff1
          END DO
        END DO
      END DO
!
!  Tracers. For now, use salinity scale for passive tracers.
!
      DO itrc=1,NT(ng)
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
        DO k=1,N(ng)
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef BNORM
              IF (e_t(i,j,k,1,itrc).gt.0.0_r8) THEN
                cff1=1.0_r8/(e_t(i,j,k,1,itrc)*e_t(i,j,k,1,itrc))
              ELSE
                cff1=1.0_r8
              END IF
              cff1=cff1*tl_t(i,j,k,Ninp,itrc)*                          &
     &                  tl_t(i,j,k,Ninp,itrc)
#   else
              scale=cff*Hz(i,j,k)
              cff1=scale*tl_t(i,j,k,Ninp,itrc)*                         &
     &                   tl_t(i,j,k,Ninp,itrc)
#   endif
#   ifdef MASKING
              cff1=cff1*rmask(i,j)
#   endif
              my_norm=my_norm+cff1
            END DO
          END DO
        END DO
      END DO
#  endif
!
!  Perform parallel global reduction operation: dot product.
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
!$OMP CRITICAL (TL_NORM)
      IF (tile_count.eq.0) THEN
        StateNorm=my_norm
      ELSE
        StateNorm=StateNorm+my_norm
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, iTLM, 1, StateNorm, op_handle)
#  endif
      END IF
!$OMP END CRITICAL (TL_NORM)

#  undef IR_RANGE
#  undef IU_RANGE
#  undef JR_RANGE
#  undef JV_RANGE

      RETURN
      END SUBROUTINE tl_statenorm_tile
# endif

# ifdef TLM_CHECK
      SUBROUTINE tl_dotproduct (ng, tile, Linp)
!
!=======================================================================
!                                                                      !
!  This routine computes dot product function (g1) associated with     !
!  the tangent linear solution.                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 7)
#  endif
      CALL tl_dotproduct_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         KOUT, Linp,                              &
#  ifdef SOLVE3D
     &                         NOUT,                                    &
#  endif
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
#  endif
#  ifdef SOLVE3D
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % tl_t,                        &
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % tl_u,                        &
     &                         OCEAN(ng) % ad_v,                        &
     &                         OCEAN(ng) % tl_v,                        &
#  endif
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % tl_ubar,                     &
     &                         OCEAN(ng) % ad_vbar,                     &
     &                         OCEAN(ng) % tl_vbar,                     &
     &                         OCEAN(ng) % ad_zeta,                     &
     &                         OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
#  endif
      RETURN
      END SUBROUTINE tl_dotproduct

!
!***********************************************************************
      SUBROUTINE tl_dotproduct_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Kinp, Linp,                        &
#  ifdef SOLVE3D
     &                               Ninp,                              &
#  endif
#  ifdef MASKING
     &                               rmask, umask, vmask,               &
#  endif
#  ifdef SOLVE3D
     &                               ad_t, tl_t,                        &
     &                               ad_u, tl_u,                        &
     &                               ad_v, tl_v,                        &
#  endif
     &                               ad_ubar, tl_ubar,                  &
     &                               ad_vbar, tl_vbar,                  &
     &                               ad_zeta, tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kinp, Linp
#  ifdef SOLVE3D
      integer, intent(in) :: Ninp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: tl_t(LBi:,LBj:,:,:,:)

      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_u(LBi:,LBj:,:,:)

      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_ubar(LBi:,LBj:,:)

      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_vbar(LBi:,LBj:,:)

      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))

      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)

      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_ubar(LBi:UBi,LBj:UBj,3)

      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_vbar(LBi:UBi,LBj:UBj,3)

      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
      integer :: NSUB, Tindex, i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: my_dot, p, scale

#  ifdef DISTRIBUTE
      character(len=3 ) :: op_handle
#  endif
      character(len=12) :: text

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add a perturbation to nonlinear state variables: steepest descent
!  direction times the perturbation amplitude "p".
!-----------------------------------------------------------------------
!
      IF ((outer.ge.1).and.(MOD(iic(ng)-1,nTLM(ng)).eq.0).and.          &
     &    ((nrrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
!
!  Set perturbation amplitude as a function of the inner loop.
!
        p=10.0_r8**REAL(-inner,r8)
        my_dot=0.0_r8
        scale=1.0_r8/SQRT(adDotProduct)
        IF (iic(ng).eq.ntstart(ng)) THEN
          text='------------'
        ELSE
          text='            '
        END IF
!
!  Free-surface.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            my_dot=my_dot+                                              &
     &             p*scale*ad_zeta(i,j,Linp)*tl_zeta(i,j,Kinp)
#  ifdef MASKING
            my_dot=my_dot*rmask(i,j)
#  endif
          END DO
        END DO
!
!  2D u-momentum component.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            my_dot=my_dot+                                              &
     &             p*scale*ad_ubar(i,j,Linp)*tl_ubar(i,j,Kinp)
#  ifdef MASKING
            my_dot=my_dot*umask(i,j)
#  endif
          END DO
        END DO
!
!  2D v-momentum component.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            my_dot=my_dot+                                              &
     &             p*scale*ad_vbar(i,j,Linp)*tl_vbar(i,j,Kinp)
#  ifdef MASKING
            my_dot=my_dot*vmask(i,j)
#  endif
          END DO
        END DO
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              my_dot=my_dot+                                            &
     &               p*scale*ad_u(i,j,k,Linp)*tl_u(i,j,k,Ninp)
#   ifdef MASKING
              my_dot=my_dot*umask(i,j)
#   endif
            END DO
          END DO
        END DO
!
!  3D v-momentum component.
!
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              my_dot=my_dot+                                            &
     &               p*scale*ad_v(i,j,k,Linp)*tl_v(i,j,k,Ninp)
#   ifdef MASKING
              my_dot=my_dot*vmask(i,j)
#   endif
            END DO
          END DO
        END DO
!
!  Tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                my_dot=my_dot+                                          &
     &                 p*scale*ad_t(i,j,k,Linp,itrc)*                   &
     &                         tl_t(i,j,k,Ninp,itrc)
#   ifdef MASKING
                my_dot=my_dot*rmask(i,j)
#   endif
              END DO
            END DO
          END DO
        END DO
#  endif
!
!  Perform parallel global reduction operation: dot product.
!
#  ifdef DISTRIBUTE
        NSUB=1                           ! distributed-memory
#  else
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
#  endif
!$OMP CRITICAL (TL_DOT)
        IF (tile_count.eq.0) THEN
          DotProduct=my_dot
        ELSE
          DotProduct=DotProduct+my_dot
        END IF
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#  ifdef DISTRIBUTE
          op_handle='SUM'
          CALL mp_reduce (ng, iTLM, 1, DotProduct, op_handle)
#  endif
          ig2count=ig2count+1
          g2(ig2count)=DotProduct
          IF (Master) THEN
            WRITE (stdout,10) outer, inner, p, text, DotProduct
 10         FORMAT (5x,'(Outer,Inner) = ','(',i4.4,',',i4.4,')',3x,     &
     &              'TLM Check, p  = ',1p,e21.14,0p,/,5x,a,18x,         &
     &              'TLM Check, g2 = ',1p,e21.14,0p)
          END IF
        END IF
!$OMP END CRITICAL (TL_DOT)
      END IF

      RETURN
      END SUBROUTINE tl_dotproduct_tile
# endif

# ifdef TLM_CHECK
      SUBROUTINE ad_dotproduct (ng, tile, Linp)
!
!=======================================================================
!                                                                      !
!  This routine computes the adjoint solution dot product to scale     !
!  the adjoint solution.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
#  endif
      CALL ad_dotproduct_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Linp,                                    &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
     &                         GRID(ng) % umask,                        &
     &                         GRID(ng) % vmask,                        &
#  endif
#  ifdef SOLVE3D
     &                         OCEAN(ng) % ad_t,                        &
     &                         OCEAN(ng) % ad_u,                        &
     &                         OCEAN(ng) % ad_v,                        &
#  endif
     &                         OCEAN(ng) % ad_ubar,                     &
     &                         OCEAN(ng) % ad_vbar,                     &
     &                         OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
#  endif
      RETURN
      END SUBROUTINE ad_dotproduct
!
!***********************************************************************
      SUBROUTINE ad_dotproduct_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
     &                               Linp,                              &
#  ifdef MASKING
     &                               rmask, umask, vmask,               &
#  endif
#  ifdef SOLVE3D
     &                               ad_t, ad_u, ad_v,                  &
#  endif
     &                               ad_ubar, ad_vbar, ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_scalars

#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Linp
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(in) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ad_v(LBi:,LBj:,:,:)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
#   ifdef SOLVE3D
      real(r8), intent(in) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(in) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
      real(r8), intent(in) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
      integer :: NSUB, i, j
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
      real(r8) :: my_dot

#  ifdef DISTRIBUTE
      character(len=3) :: op_handle
#  endif

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute adjoint solution dot product.
!-----------------------------------------------------------------------
!
      my_dot=0.0_r8
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          my_dot=my_dot+                                                &
     &           ad_zeta(i,j,Linp)*ad_zeta(i,j,Linp)
#  ifdef MASKING
          my_dot=my_dot*rmask(i,j)
#  endif
        END DO
      END DO
#  ifndef SOLVE3D
!
!  2D u-momentum component.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          my_dot=my_dot+                                                &
     &           ad_ubar(i,j,Linp)*ad_ubar(i,j,Linp)
#   ifdef MASKING
          my_dot=my_dot*umask(i,j)
#   endif
        END DO
      END DO
!
!  2D v-momentum component.
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          my_dot=my_dot+                                                &
     &           ad_vbar(i,j,Linp)*ad_vbar(i,j,Linp)
#   ifdef MASKING
          my_dot=my_dot*vmask(i,j)
#   endif
        END DO
      END DO
#  endif
#  ifdef SOLVE3D
!
!  3D u-momentum component.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            my_dot=my_dot+                                              &
     &             ad_u(i,j,k,Linp)*ad_u(i,j,k,Linp)
#   ifdef MASKING
            my_dot=my_dot*umask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  3D v-momentum component.
!
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            my_dot=my_dot+                                              &
     &             ad_v(i,j,k,Linp)*ad_v(i,j,k,Linp)
#   ifdef MASKING
            my_dot=my_dot*vmask(i,j)
#   endif
          END DO
        END DO
      END DO
!
!  Tracers.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              my_dot=my_dot+                                            &
     &               ad_t(i,j,k,Linp,itrc)*ad_t(i,j,k,Linp,itrc)
#   ifdef MASKING
              my_dot=my_dot*rmask(i,j)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif
!
!  Perform parallel global reduction operation: dot product.
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
!$OMP CRITICAL (AD_DOT)
      IF (tile_count.eq.0) THEN
        adDotProduct=my_dot
      ELSE
        adDotProduct=adDotProduct+my_dot
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, iADM, 1, adDotProduct, op_handle)
#  endif
        IF (Master) THEN
          WRITE (stdout,10) adDotProduct
 10       FORMAT (/,' Adjoint Solution Dot Product: ',1p,e21.14)
        END IF
      END IF
!$OMP END CRITICAL (AD_DOT)

      RETURN
      END SUBROUTINE ad_dotproduct_tile
# endif

#endif
      END MODULE dotproduct_mod
