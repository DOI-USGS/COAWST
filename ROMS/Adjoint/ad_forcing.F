#include "cppdefs.h"
      MODULE ad_forcing_mod

#if defined ADJOINT && defined FORCING_SV
!
!svn $Id: ad_forcing.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Andrew M. Moore    !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is used to force the adjoint state equations.          !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ad_forcing

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_forcing (ng, tile, Kfrc, Nfrc)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
# ifdef SOLVE3D
      USE mod_coupling
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Kfrc, Nfrc
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ad_forcing_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      Kfrc, Nfrc,                                 &
# ifdef SOLVE3D
     &                      OCEAN(ng) % f_t,                            &
     &                      OCEAN(ng) % f_u,                            &
     &                      OCEAN(ng) % f_v,                            &
#  ifdef FORCING_SV
     &                      OCEAN(ng) % f_ubar,                         &
     &                      OCEAN(ng) % f_vbar,                         &
#  endif
# else
     &                      OCEAN(ng) % f_ubar,                         &
     &                      OCEAN(ng) % f_vbar,                         &
# endif
     &                      OCEAN(ng) % f_zeta,                         &
# ifdef SOLVE3D
     &                      OCEAN(ng) % ad_t,                           &
     &                      OCEAN(ng) % ad_u,                           &
     &                      OCEAN(ng) % ad_v,                           &
#  ifdef FORCING_SV
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
#  endif
# else
     &                      OCEAN(ng) % ad_ubar,                        &
     &                      OCEAN(ng) % ad_vbar,                        &
# endif
# ifdef SOLVE3D
     &                      COUPLING(ng) % ad_Zt_avg1,                  &
# endif
     &                      OCEAN(ng) % ad_zeta,                        &
# ifdef SOLVE3D
     &                      OCEAN(ng) % ad_t_sol,                       &
     &                      OCEAN(ng) % ad_u_sol,                       &
     &                      OCEAN(ng) % ad_v_sol,                       &
# else
     &                      OCEAN(ng) % ad_ubar_sol,                    &
     &                      OCEAN(ng) % ad_vbar_sol,                    &
# endif
     &                      OCEAN(ng) % ad_zeta_sol)

      RETURN
      END SUBROUTINE ad_forcing
!
!***********************************************************************
      SUBROUTINE ad_forcing_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Kfrc, Nfrc,                           &
# ifdef SOLVE3D
     &                            f_t, f_u, f_v,                        &
#  ifdef FORCING_SV
     &                            f_ubar, f_vbar,                       &
#  endif
# else
     &                            f_ubar, f_vbar,                       &
# endif
     &                            f_zeta,                               &
# ifdef SOLVE3D
     &                            ad_t, ad_u, ad_v,                     &
#  ifdef FORCING_SV
     &                            ad_ubar, ad_vbar,                     &
#  endif
# else
     &                            ad_ubar, ad_vbar,                     &
# endif
# ifdef SOLVE3D
     &                            ad_Zt_avg1,                           &
# endif
     &                            ad_zeta,                              &
# ifdef SOLVE3D
     &                            ad_t_sol, ad_u_sol, ad_v_sol,         &
# else
     &                            ad_ubar_sol, ad_vbar_sol,             &
# endif
     &                            ad_zeta_sol)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kfrc
      integer, intent(in) :: Nfrc
!
# ifdef ASSUMED_SHAPE
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_t(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: f_u(LBi:,LBj:,:)
      real(r8), intent(inout) :: f_v(LBi:,LBj:,:)
#   ifdef FORCING_SV
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#   endif
#  else
      real(r8), intent(inout) :: f_ubar(LBi:,LBj:)
      real(r8), intent(inout) :: f_vbar(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   ifdef FORCING_SV
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_sol(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_u_sol(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_v_sol(LBi:,LBj:,:)
#  else
      real(r8), intent(inout) :: ad_ubar_sol(LBi:,LBj:)
      real(r8), intent(inout) :: ad_vbar(_solLBi:,LBj:)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:,LBj:)
# else
#  ifdef SOLVE3D
      real(r8), intent(inout) :: f_t(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: f_u(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: f_v(LBi:UBi,LBj:UBj,N(ng))
#   ifdef FORCING_SV
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#   endif
#  else
      real(r8), intent(inout) :: f_ubar(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: f_vbar(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: f_zeta(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   ifdef FORCING_SV
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
#  else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#  endif
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_Zt_avg1(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  ifdef SOLVE3D
      real(r8), intent(inout) :: ad_t_sol(LBi:UBi,LBj:UBj,N(ng),NT(ng))
      real(r8), intent(inout) :: ad_u_sol(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: ad_v_sol(LBi:UBi,LBj:UBj,N(ng))
#  else
      real(r8), intent(inout) :: ad_ubar_sol(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_vbar_sol(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: ad_zeta_sol(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc, k
# endif

# include "set_bounds.h"
!
!
!  Adjoint linear free-surface. The two different cases in the
!  case of SOLVE3D are due to the fact that tl_ini_fields is
!  also called on the first timestep. tl_forcing MUST be called
!  before tl_ini_fields.
!
#  ifdef SOLVE3D
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f_zeta(i,j)=f_zeta(i,j)+ad_Zt_avg1(i,j)
        END DO
      END DO
#  else
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f_zeta(i,j)=f_zeta(i,j)+ad_zeta(i,j,Kfrc)
        END DO
      END DO
#  endif

# ifndef SOLVE3D
!
!  Adjoint linear 2D momentum.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            f_ubar(i,j)=f_ubar(i,j)+ad_ubar(i,j,Kfrc)
          END DO
        END DO
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            f_vbar(i,j)=f_vbar(i,j)+ad_vbar(i,j,Kfrc)
          END DO
        END DO

# else
!
!  Adjoint linear 2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          f_ubar(i,j)=f_ubar(i,j)+ad_ubar(i,j,1)+ad_ubar(i,j,2)
        END DO
      END DO
!
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          f_vbar(i,j)=f_vbar(i,j)+ad_vbar(i,j,1)+ad_vbar(i,j,2)
        END DO
      END DO
!
!  Adjoint linear 3D momentum.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              f_u(i,j,k)=f_u(i,j,k)+ad_u(i,j,k,Nfrc)
            END DO
          END DO
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              f_v(i,j,k)=f_v(i,j,k)+ad_v(i,j,k,Nfrc)
            END DO
          END DO
        END DO
!
!  Adjoint linear tracers.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                f_t(i,j,k,itrc)=f_t(i,j,k,itrc)+                        &
     &                            ad_t(i,j,k,Nfrc,itrc)
            END DO
            END DO
          END DO
        END DO
# endif

      RETURN
      END SUBROUTINE ad_forcing_tile
#endif
      END MODULE ad_forcing_mod
