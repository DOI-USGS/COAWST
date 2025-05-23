#include "cppdefs.h"
      MODULE rp_frc_adjust_mod

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!git $Id$
!svn $Id: rp_frc_adjust.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-interpolates  4DVar representer model surface     !
!  forcing increments. The increments can be constant (Nfrec=1) or     !
!  time interpolated between snapshots (Nfrec>1).                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      Representer model state time index to process.         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: rp_frc_adjust

      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_frc_adjust (ng, tile, Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iRPM, 7, __LINE__, MyFile)
# endif
      CALL rp_frc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef ADJUST_WSTRESS
     &                         FORCES(ng) % tl_ustr,                    &
     &                         FORCES(ng) % tl_vstr,                    &
     &                         FORCES(ng) % tl_sustr,                   &
     &                         FORCES(ng) % tl_svstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                         FORCES(ng) % tl_tflux,                   &
     &                         FORCES(ng) % tl_stflx,                   &
# endif
     &                         Linp)
# ifdef PROFILE
      CALL wclock_off (ng, iRPM, 7, __LINE__, MyFile)
# endif
!
      RETURN
      END SUBROUTINE rp_frc_adjust
!
!***********************************************************************
      SUBROUTINE rp_frc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
# ifdef ADJUST_WSTRESS
     &                               tl_ustr,  tl_vstr,                 &
     &                               tl_sustr, tl_svstr,                &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                               tl_tflux, tl_stflx,                &
# endif
     &                               Linp)
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
      integer, intent(in) :: Linp
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: tl_svstr(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_stflx(LBi:,LBj:,:)
#  endif
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tl_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: tl_svstr(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, it1, it2, j
# ifdef SOLVE3D
      integer :: itrc
# endif
      real(r8) :: fac, fac1, fac2

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Adjust representer model surface forcing fields with 4DVar
!  increments.
!-----------------------------------------------------------------------
!
!  Set time records and interpolation factor, if any.
!
      IF (Nfrec(ng).eq.1) THEN
        it1=1
        it2=1
        fac1=1.0_r8
        fac2=0.0_r8
      ELSE
# ifdef GENERIC_DSTART
        it1=MAX(0,(iic(ng)-ntstart(ng))/nSFF(ng))+1
# else
        it1=MAX(0,(iic(ng)-1)/nSFF(ng))+1
# endif
        it2=MIN(it1+1,Nfrec(ng))
        fac1=SF_time(it2,ng)-(time(ng)+dt(ng))
        fac2=(time(ng)+dt(ng))-SF_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF

# ifdef ADJUST_WSTRESS
!
!  Adjust surface wind stress. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
           tl_sustr(i,j)=tl_sustr(i,j)+                                 &
     &                   fac1*tl_ustr(i,j,it1,Linp)+                    &
     &                   fac2*tl_ustr(i,j,it2,Linp)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
           tl_svstr(i,j)=tl_svstr(i,j)+                                 &
     &                   fac1*tl_vstr(i,j,it1,Linp)+                    &
     &                   fac2*tl_vstr(i,j,it2,Linp)
        END DO
      END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjust surface tracer fluxes. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
               tl_stflx(i,j,itrc)=tl_stflx(i,j,itrc)+                   &
     &                            fac1*tl_tflux(i,j,it1,Linp,itrc)+     &
     &                            fac2*tl_tflux(i,j,it2,Linp,itrc)
            END DO
          END DO
        END IF
      END DO
# endif
!
      RETURN
      END SUBROUTINE rp_frc_adjust_tile
#endif
      END MODULE rp_frc_adjust_mod
