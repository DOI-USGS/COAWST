#include "cppdefs.h"
      MODULE frc_iau_mod

#if defined RBL4DVAR && defined RPCG
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group           Andrew M. Moore    !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module is used to compute the Incremental Analysis Update      !
!  (IAU) forcing for the NLM state vector during the "analysis phase". !
!                                                                      !
!  After 4D-Var increment phase, a fraction of the variable increment  !
!  (dt/timeIAU) is added as force term to the state variables to       !
!  reduce the "initialization shocks" from large amplitude inertia-    !
!  gravity waves generated due to unbalanced increment contribution.   !
!                                                                      !
!  In weak constraint 4D-Var, the IAU forcing term is added to the     !
!  convolved impulse term already time intepolated into the f_var      !
!  arrays in routine "frc_iau" below.                                  !
!                                                                      !
!  The method is described in:                                         !
!                                                                      !
!  Bloom et al., 1996: Data Assimilation Using Incremental Analysis    !
!                      Update. Mon. Wea. Rev., 124, 1256-1271.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_fourdvar
      USE mod_grid
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
      PUBLIC  :: frc_iau
      PUBLIC  :: frc_iau_ini
      PRIVATE
!
      CONTAINS
!
      SUBROUTINE frc_iau_ini (ng, tile, model, Linp)
!
!=======================================================================
!                                                                      !
!  This subroutine computes the incremental analysis update (IAU)      !
!  that is applied as a forcing term in the non-linear model during    !
!  the first timeIAU seconds of the assimilation window.               !
!                                                                      !
!  The tangent linear state arrays tl_var(Linp) are used as temporary  !
!  storage for the IAU.                                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     model     Calling model identifier.                              !
!     Linp      TLM variable time index containing the increment.      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 7, __LINE__, __FILE__)
# endif
      CALL frc_iau_ini_tile (ng, tile, model, Linp,                     &
     &                       LBi, UBi, LBj, UBj,                        &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
# ifdef SOLVE3D
     &                       OCEAN(ng) % d_u,                           &
     &                       OCEAN(ng) % d_v,                           &
     &                       OCEAN(ng) % d_t,                           &
# else
     &                       OCEAN(ng) % d_ubar,                        &
     &                       OCEAN(ng) % d_vbar,                        &
# endif
     &                       OCEAN(ng) % d_zeta,                        &
# ifdef SOLVE3D
     &                       OCEAN(ng) % f_u,                           &
     &                       OCEAN(ng) % f_v,                           &
     &                       OCEAN(ng) % f_t,                           &
# else
     &                       OCEAN(ng) % f_ubar,                        &
     &                       OCEAN(ng) % f_vbar,                        &
# endif
     &                       OCEAN(ng) % f_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, model, 7, __LINE__, __FILE__)
# endif
!
      RETURN
      END SUBROUTINE frc_iau_ini
!
!***********************************************************************
      SUBROUTINE frc_iau_ini_tile (ng, tile, model, Linp,               &
     &                             LBi, UBi, LBj, UBj,                  &
# ifdef MASKING
     &                             rmask, umask, vmask,                 &
# endif
# ifdef SOLVE3D
     &                             d_u, d_v, d_t,                       &
# else
     &                             d_ubar, d_vbar,                      &
# endif
     &                             d_zeta,                              &
# ifdef SOLVE3D
     &                             f_u, f_v, f_t,                       &
# else
     &                             f_ubar, f_vbar,                      &
# endif
     &                             f_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Linp
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:,LBj:)
      real(r8), intent(in   ) :: umask(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: d_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: d_t(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: d_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)

# else

#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  else
      real(r8), intent(inout) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: d_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  else
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations
      integer :: i, it1, it2, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif
      real(r8) :: fac1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize Incremental Analysis Update (IAU) forcing terms used in
!  the nonlinear model. The cummulative increment from 4D-Var outer
!  loop is stored in the work "d_*" arrays, and it is gradually added
!  over the specified time interval. 
!
!  The "d_*" arrays are overwritten by the fractional increment to be
!  added during time-stepping until the end filtering period.
!-----------------------------------------------------------------------
!
!  Set uniform weights in time.
!
      fac1=REAL(timeIAU(ng)/dt(ng),r8)
      fac1=1.0_r8/fac1
!
!  Compute free-surface IAU forcing term.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          f_zeta(i,j)=fac1*d_zeta(i,j)
# ifdef MASKING
          f_zeta(i,j)=f_zeta(i,j)*rmask(i,j)
# endif
          d_zeta(i,j)=f_zeta(i,j)
        END DO
      END DO

# ifndef SOLVE3D
!
!  Compute 2D-momentum IAU forcing terms.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          f_ubar(i,j)=fac1*d_ubar(i,j)
#  ifdef MASKING
          f_ubar(i,j)=f_ubar(i,j)*umask(i,j)
#  endif
          d_ubar(i,j)=f_ubar(i,j)
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          f_vbar(i,j)=fac1*d_vbar(i,j)
#  ifdef MASKING
          f_vbar(i,j)=f_vbar(i,j)*vmask(i,j)
#  endif
          d_vbar(i,j)=f_vbar(i,j)
        END DO
      END DO
# endif
# ifdef SOLVE3D
!
!  Compute 3D-momentum IAU forcing terms.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            f_u(i,j,k)=fac1*d_u(i,j,k)
#  ifdef MASKING
            f_u(i,j,k)=f_u(i,j,k)*umask(i,j)
#  endif
            d_u(i,j,k)=f_u(i,j,k)
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            f_v(i,j,k)=fac1*d_v(i,j,k)
#  ifdef MASKING
            f_v(i,j,k)=f_v(i,j,k)*vmask(i,j)
#  endif
            d_v(i,j,k)=f_v(i,j,k)
          END DO
        END DO
      END DO
!
!  Compute tracer IAU forcing terms.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              f_t(i,j,k,itrc)=fac1*d_t(i,j,k,itrc)
#  ifdef MASKING
              f_t(i,j,k,itrc)=f_t(i,j,k,itrc)*rmask(i,j)
#  endif
              d_t(i,j,k,itrc)=f_t(i,j,k,itrc)
            END DO
          END DO
        END DO
      END DO
# endif
!
      RETURN
      END SUBROUTINE frc_iau_ini_tile
!
      SUBROUTINE frc_iau (ng, tile, model, Linp)
!
!=======================================================================
!                                                                      !
!  This subroutine computes the combined weak constraint forcing       !
!  and the Incremental Analysis Update (IAU) when both weak constraint !
!  data assimilation and the IAU are activated.                        !
!                                                                      !
!  The input tangent linear state arrays "tl_var(Linp)" hold the IAU   !
!  forcing.                                                            !
!                                                                      !
!  The "f_var" arrays hold the time interpolated weak constraint       !
!  impulses processed from the TLF 4D-Var file.                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     model     Calling model identifier.                              !
!     Linp      TLM variable time index containing the increment.      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Linp
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 7, __LINE__, __FILE__)
# endif
      CALL frc_iau_tile (ng, tile, model, Linp,                         &
     &                   LBi, UBi, LBj, UBj,                            &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
# endif
# ifdef SOLVE3D
     &                   OCEAN(ng) % d_u,                               &
     &                   OCEAN(ng) % d_v,                               &
     &                   OCEAN(ng) % d_t,                               &
# else
     &                   OCEAN(ng) % d_ubar,                            &
     &                   OCEAN(ng) % d_vbar,                            &
# endif
     &                   OCEAN(ng) % d_zeta,                            &
# ifdef SOLVE3D
     &                   OCEAN(ng) % f_u,                               &
     &                   OCEAN(ng) % f_v,                               &
     &                   OCEAN(ng) % f_t,                               &
# else
     &                   OCEAN(ng) % f_ubar,                            &
     &                   OCEAN(ng) % f_vbar,                            &
# endif
     &                   OCEAN(ng) % f_zeta)
# ifdef PROFILE
      CALL wclock_off (ng, model, 7, __LINE__, __FILE__)
# endif
!
      RETURN
      END SUBROUTINE frc_iau
!
!***********************************************************************
      SUBROUTINE frc_iau_tile (ng, tile, model, Linp,                   &
     &                         LBi, UBi, LBj, UBj,                      &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef SOLVE3D
     &                         d_u, d_v, d_t,                           &
# else
     &                         d_ubar, d_vbar,                          &
# endif
     &                         d_zeta,                                  &
# ifdef SOLVE3D
     &                         f_u, f_v, f_t,                           &
# else
     &                         f_ubar, f_vbar,                          &
# endif
     &                         f_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Linp
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:,LBj:)
      real(r8), intent(in   ) :: umask(LBi:,LBj:)
      real(r8), intent(in   ) :: vmask(LBi:,LBj:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in   ) :: d_u(LBi:,LBj:,:)
      real(r8), intent(in   ) :: d_v(LBi:,LBj:,:)
      real(r8), intent(in   ) :: d_t(LBi:,LBj:,:,:)
#  else
      real(r8), intent(in   ) :: d_ubar(LBi:,LBj:)
      real(r8), intent(in   ) :: d_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(in   ) :: d_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)

# else

#  ifdef MASKING
      real(r8), intent(in   ) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(in   ) :: d_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in   ) :: d_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in   ) :: d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  else
      real(r8), intent(in   ) :: d_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(in   ) :: d_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in   ) :: d_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
#  else
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j, k
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Update the weak constraint forcing arrays with the IAU frational
!  forcing contribution, which is stored in the work "d_*" arrays.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          f_zeta(i,j)=f_zeta(i,j)+d_zeta(i,j)
# ifdef MASKING
          f_zeta(i,j)=f_zeta(i,j)*rmask(i,j)
# endif
        END DO
      END DO

# ifndef SOLVE3D
!
!  Add 2D-momentum IAU forcing terms.
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          f_ubar(i,j)=f_ubar(i,j)+d_ubar(i,j)
#  ifdef MASKING
          f_ubar(i,j)=f_ubar(i,j)*umask(i,j)
#  endif
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          f_vbar(i,j)=f_vbar(i,j)+d_vbar(i,j)
#  ifdef MASKING
          f_vbar(i,j)=f_vbar(i,j)*vmask(i,j)
#  endif
        END DO
      END DO
# endif
# ifdef SOLVE3D
!
!  Add 3D-momentum IAU forcing terms.
!
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            f_u(i,j,k)=f_u(i,j,k)+d_u(i,j,k)
#  ifdef MASKING
            f_u(i,j,k)=f_u(i,j,k)*umask(i,j)
#  endif
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            f_v(i,j,k)=f_v(i,j,k)+d_v(i,j,k)
#  ifdef MASKING
            f_v(i,j,k)=f_v(i,j,k)*vmask(i,j)
#  endif
          END DO
        END DO
      END DO
!
!  Add tracer iau forcing terms.
!
      DO itrc=1,NT(ng)
        DO k=1,N(ng)
          DO j=JstrT,JendT
            DO i=IstrT,IendT
              f_t(i,j,k,itrc)=f_t(i,j,k,itrc)+d_t(i,j,k,itrc)
#  ifdef MASKING
              f_t(i,j,k,itrc)=f_t(i,j,k,itrc)*rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
!
      RETURN
      END SUBROUTINE frc_iau_tile
#endif
      END MODULE frc_iau_mod
