#include "cppdefs.h"
      MODULE state_addition_mod
!
!svn $Id: state_addition.F 584 2008-03-18 20:01:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the following model state addition:           !
!                                                                      !
!      s1_var(...,Lout) = fac1 * s1_var(...,Linp1) +                   !
!                         fac2 * s2_var(...,Linp2)                     !
!                                                                      !
!  where alpha and beta are scalars.                                   !
!                                                                      !
!-----------------------------------------------------------------------
!
      implicit none

      PUBLIC  :: state_addition

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_addition (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Linp1, Linp2, Lout,                    &
     &                           fac1, fac2,                            &
#ifdef MASKING
     &                           rmask, umask, vmask,                   &
#endif
#ifdef ADJUST_WSTRESS
     &                           s1_sustr, s2_sustr,                    &
     &                           s1_svstr, s2_svstr,                    &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                           s1_tflux, s2_tflux,                    &
# endif
     &                           s1_t, s2_t,                            &
     &                           s1_u, s2_u,                            &
     &                           s1_v, s2_v,                            &
#else
     &                           s1_ubar, s2_ubar,                      &
     &                           s1_vbar, s2_vbar,                      &
#endif
     &                           s1_zeta, s2_zeta)
!***********************************************************************
!
      USE mod_param
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_scalars
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp1, Linp2, Lout
!
      real(r8), intent(in) :: fac1, fac2
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s2_sustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s2_svstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s2_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s2_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s2_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s2_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s2_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s2_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s2_zeta(LBi:,LBj:,:)

# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_svstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s1_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s1_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s1_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s1_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s1_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s1_zeta(LBi:,LBj:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s2_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s2_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s2_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s2_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s2_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s2_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s2_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s2_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s2_zeta(LBi:UBi,LBj:UBj,3)

# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s1_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s1_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s1_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s1_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s1_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s1_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s1_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s1_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s1_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
#ifdef SOLVE3D
      integer :: itrc
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Compute the following operation between S1 and S2 model state
!  trajectories:
!                 S1(Lout) = fac1 * S1(Linp1) + fac2 * S2(Linp2)
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          s1_zeta(i,j,Lout)=fac1*s1_zeta(i,j,Linp1)+                    &
     &                      fac2*s2_zeta(i,j,Linp2)
#ifdef MASKING
          s1_zeta(i,j,Lout)=s1_zeta(i,j,Lout)*rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          s1_ubar(i,j,Lout)=fac1*s1_ubar(i,j,Linp1)+                    &
     &                      fac2*s2_ubar(i,j,Linp2)
# ifdef MASKING
          s1_ubar(i,j,Lout)=s1_ubar(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          s1_vbar(i,j,Lout)=fac1*s1_vbar(i,j,Linp1)+                    &
     &                      fac2*s2_vbar(i,j,Linp2)
# ifdef MASKING
          s1_vbar(i,j,Lout)=s1_vbar(i,j,Lout)*vmask(i,j)
# endif
        END DO
      END DO
#endif

#ifdef ADJUST_WSTRESS
!
!  Surface momentum stress.
!
      DO k=1,Nfrec(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s1_sustr(i,j,k,Lout)=fac1*s1_sustr(i,j,k,Linp1)+            &
     &                           fac2*s2_sustr(i,j,k,Linp2)
# ifdef MASKING
            s1_sustr(i,j,k,Lout)=s1_sustr(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s1_svstr(i,j,k,Lout)=fac1*s1_svstr(i,j,k,Linp1)+            &
     &                           fac2*s2_svstr(i,j,k,Linp2)
# ifdef MASKING
            s1_svstr(i,j,k,Lout)=s1_svstr(i,j,k,Lout)*vmask(i,j)
# endif
          END DO
        END DO
      END DO
#endif

#ifdef SOLVE3D
!
!  3D momentum.
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=Istr,IendR
            s1_u(i,j,k,Lout)=fac1*s1_u(i,j,k,Linp1)+                    &
     &                       fac2*s2_u(i,j,k,Linp2)
# ifdef MASKING
            s1_u(i,j,k,Lout)=s1_u(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s1_v(i,j,k,Lout)=fac1*s1_v(i,j,k,Linp1)+                    &
     &                       fac2*s2_v(i,j,k,Linp2)
# ifdef MASKING
            s1_v(i,j,k,Lout)=s1_v(i,j,k,Lout)*vmask(i,j)
# endif
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
              s1_t(i,j,k,Lout,itrc)=fac1*s1_t(i,j,k,Linp1,itrc)+        &
     &                              fac2*s2_t(i,j,k,Linp2,itrc)
# ifdef MASKING
              s1_t(i,j,k,Lout,itrc)=s1_t(i,j,k,Lout,itrc)*rmask(i,j)
# endif
            END DO
          END DO
        END DO
      END DO

# ifdef ADJUST_STFLUX
!
!  Surface tracers flux.
!
      DO itrc=1,NT(ng)
        DO k=1,Nfrec(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              s1_tflux(i,j,k,Lout,itrc)=fac1*s1_tflux(i,j,k,Linp1,itrc)+&
     &                                  fac2*s2_tflux(i,j,k,Linp2,itrc)
#  ifdef MASKING
              s1_tflux(i,j,k,Lout,itrc)=s1_tflux(i,j,k,Lout,itrc)*      &
     &                                  rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE state_addition

      END MODULE state_addition_mod
