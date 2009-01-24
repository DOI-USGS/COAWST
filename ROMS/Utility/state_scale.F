#include "cppdefs.h"
      MODULE state_scale_mod
!
!svn $Id: state_scale.F 584 2008-03-18 20:01:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine scales a state variable by a constant:                 !
!                                                                      !
!     s_var(...,Lout) = fac * s_var(...,Linp)                          !
!                                                                      !
!-----------------------------------------------------------------------
!
      implicit none

      PUBLIC  :: state_scale

      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_scale (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Linp, Lout, fac,                          &
#ifdef MASKING
     &                        rmask, umask, vmask,                      &
#endif
#ifdef ADJUST_WSTRESS
     &                        s_sustr, s_svstr,                         &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                        s_tflux,                                  &
# endif
     &                        s_t, s_u, s_v,                            &
#else
     &                        s_ubar, s_vbar,                           &
#endif
     &                        s_zeta)
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
      integer, intent(in) :: Linp, Lout
!
      real(r8), intent(in) :: fac
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_sustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_svstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_sustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_svstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
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
!  Scale model state variable by a constant.
!-----------------------------------------------------------------------
!
!  Free-surface.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          s_zeta(i,j,Lout)=fac*s_zeta(i,j,Linp)
#ifdef MASKING
          s_zeta(i,j,Lout)=s_zeta(i,j,Lout)*rmask(i,j)
#endif
        END DO
      END DO

#ifndef SOLVE3D
!
!  2D momentum.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          s_ubar(i,j,Lout)=fac*s_ubar(i,j,Linp)
# ifdef MASKING
          s_ubar(i,j,Lout)=s_ubar(i,j,Lout)*umask(i,j)
# endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          s_vbar(i,j,Lout)=fac*s_vbar(i,j,Linp)
# ifdef MASKING
          s_vbar(i,j,Lout)=s_vbar(i,j,Lout)*vmask(i,j)
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
            s_sustr(i,j,k,Lout)=fac*s_sustr(i,j,k,Linp)
# ifdef MASKING
            s_sustr(i,j,k,Lout)=s_sustr(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s_svstr(i,j,k,Lout)=fac*s_svstr(i,j,k,Linp)
# ifdef MASKING
            s_svstr(i,j,k,Lout)=s_svstr(i,j,k,Lout)*vmask(i,j)
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
            s_u(i,j,k,Lout)=fac*s_u(i,j,k,Linp)
# ifdef MASKING
            s_u(i,j,k,Lout)=s_u(i,j,k,Lout)*umask(i,j)
# endif
          END DO
        END DO
      END DO
      DO k=1,N(ng)
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            s_v(i,j,k,Lout)=fac*s_v(i,j,k,Linp)
# ifdef MASKING
            s_v(i,j,k,Lout)=s_v(i,j,k,Lout)*vmask(i,j)
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
              s_t(i,j,k,Lout,itrc)=fac*s_t(i,j,k,Linp,itrc)
# ifdef MASKING
              s_t(i,j,k,Lout,itrc)=s_t(i,j,k,Lout,itrc)*rmask(i,j)
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
              s_tflux(i,j,k,Lout,itrc)=fac*s_tflux(i,j,k,Linp,itrc)
#  ifdef MASKING
              s_tflux(i,j,k,Lout,itrc)=s_tflux(i,j,k,Lout,itrc)*        &
     &                                 rmask(i,j)
#  endif
            END DO
          END DO
        END DO
      END DO
# endif
#endif

      RETURN
      END SUBROUTINE state_scale

      END MODULE state_scale_mod
