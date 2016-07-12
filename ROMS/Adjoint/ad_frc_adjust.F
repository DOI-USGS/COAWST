#include "cppdefs.h"
      MODULE ad_frc_adjust_mod

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!svn $Id: ad_frc_adjust.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine is the adjoint of the time-interpolation of 4DVar      !
!  surface forcing increments.  The  increments  can be  constant      !
!  (Nfrec=1) or time interpolated between snapshots (Nfrec>1).         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      4DVar increment time index to process.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: ad_frc_adjust

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_frc_adjust (ng, tile, Linp)
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
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
# endif
      CALL ad_frc_adjust_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef ADJUST_WSTRESS
     &                         FORCES(ng) % ad_ustr,                    &
     &                         FORCES(ng) % ad_vstr,                    &
     &                         FORCES(ng) % ad_sustr,                   &
     &                         FORCES(ng) % ad_svstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                         FORCES(ng) % ad_tflux,                   &
     &                         FORCES(ng) % ad_stflx,                   &
# endif
     &                         Linp)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
# endif

      RETURN
      END SUBROUTINE ad_frc_adjust
!
!***********************************************************************
      SUBROUTINE ad_frc_adjust_tile (ng, tile,                          &
     &                               LBi, UBi, LBj, UBj,                &
     &                               IminS, ImaxS, JminS, JmaxS,        &
# ifdef ADJUST_WSTRESS
     &                               ad_ustr,  ad_vstr,                 &
     &                               ad_sustr, ad_svstr,                &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                               ad_tflux, ad_stflx,                &
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
      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_stflx(LBi:,LBj:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:,LBj:,:,:,:)
#  endif
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ad_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: ad_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: ad_tflux(LBi:UBi,LBj:UBj,              &
     &                                 Nfrec(ng),2,NT(ng))
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
!  Adjust nonlinear surface forcing fields with 4DVAR increments.
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
        it1=MAX(0,(iic(ng)-1)/nSFF(ng))+1
        it2=MIN(it1+1,Nfrec(ng))
        fac1=SF_time(it2,ng)-(time(ng)+dt(ng))
        fac2=(time(ng)+dt(ng))-SF_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF

# ifdef ADJUST_WSTRESS
!
!  Adjoint of surface wind stress adjustment.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
!>        tl_sustr(i,j)=fac1*tl_ustr(i,j,it1,Linp)+                     &
!>   &                  fac2*tl_ustr(i,j,it2,Linp)
!>
          ad_ustr(i,j,it1,Linp)=ad_ustr(i,j,it1,Linp)+                  &
     &                          fac1*ad_sustr(i,j)
          ad_ustr(i,j,it2,Linp)=ad_ustr(i,j,it2,Linp)+                  &
     &                           fac2*ad_sustr(i,j)
          ad_sustr(i,j)=0.0_r8
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
!>        tl_svstr(i,j)=fac1*tl_vstr(i,j,it1,Linp)+                     &
!>   &                  fac2*tl_vstr(i,j,it2,Linp)
!>
          ad_vstr(i,j,it1,Linp)=ad_vstr(i,j,it1,Linp)+                  &
     &                          fac1*ad_svstr(i,j)
          ad_vstr(i,j,it2,Linp)=ad_vstr(i,j,it2,Linp)+                  &
     &                          fac2*ad_svstr(i,j)
          ad_svstr(i,j)=0.0_r8
        END DO
      END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjoint of surface tracer fluxes adjustmenst.
!
      IF (Lstflux(itemp,ng)) THEN
#  ifdef QCORRECTION
        DO j=JstrR,JendR
          DO i=IstrR,IendR
!>          tl_stflx(i,j,itemp)=tl_stflx(i,j,itemp)+                    &
!>   &                          fac1*tl_tflux(i,j,it1,Linp,itemp)+      &
!>   &                          fac2*tl_tflux(i,j,it2,Linp,itemp)
!>
            ad_tflux(i,j,it1,Linp,itemp)=ad_tflux(i,j,it1,Linp,itemp)+  &
     &                                   fac1*ad_stflx(i,j,itemp)
            ad_tflux(i,j,it2,Linp,itemp)=ad_tflux(i,j,it2,Linp,itemp)+  &
     &                                   fac2*ad_stflx(i,j,itemp)
          END DO
        END DO
#  else
        DO j=JstrR,JendR
          DO i=IstrR,IendR
!>          tl_stflx(i,j,itemp)=fac1*tl_tflux(i,j,it1,Linp,itemp)+      &
!>   &                          fac2*tl_tflux(i,j,it2,Linp,itemp)
!>
            ad_tflux(i,j,it1,Linp,itemp)=ad_tflux(i,j,it1,Linp,itemp)+  &
     &                                   fac1*ad_stflx(i,j,itemp)
            ad_tflux(i,j,it2,Linp,itemp)=ad_tflux(i,j,it2,Linp,itemp)+  &
     &                                   fac2*ad_stflx(i,j,itemp)
            ad_stflx(i,j,itemp)=0.0_r8
          END DO
        END DO
#  endif
      END IF

      IF (Lstflux(isalt,ng)) THEN
        DO j=JstrR,JendR
          DO i=IstrR,IendR
!>          tl_stflx(i,j,isalt)=tl_stflx(i,j,isalt)+                    &
!>   &                          fac1*tl_tflux(i,j,it1,Linp,isalt)+      &
!>   &                          fac2*tl_tflux(i,j,it2,Linp,isalt)
!>
            ad_tflux(i,j,it1,Linp,isalt)=ad_tflux(i,j,it1,Linp,isalt)+  &
     &                                   fac1*ad_stflx(i,j,isalt)
            ad_tflux(i,j,it2,Linp,isalt)=ad_tflux(i,j,it2,Linp,isalt)+  &
     &                                   fac2*ad_stflx(i,j,isalt)
          END DO
        END DO
      END IF

      DO itrc=NAT+1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
!>            tl_stflx(i,j,itrc)=fac1*tl_tflux(i,j,it1,Linp,itrc)+      &
!>   &                           fac2*tl_tflux(i,j,it2,Linp,itrc)
!>
              ad_tflux(i,j,it1,Linp,itrc)=ad_tflux(i,j,it1,Linp,itrc)+  &
     &                                    fac1*ad_stflx(i,j,itrc)
              ad_tflux(i,j,it2,Linp,itrc)=ad_tflux(i,j,it2,Linp,itrc)+  &
     &                                  fac2*ad_stflx(i,j,itrc)
              ad_stflx(i,j,itrc)=0.0_r8
            END DO
          END DO
        END IF
      END DO
# endif

      RETURN
      END SUBROUTINE ad_frc_adjust_tile
#endif
      END MODULE ad_frc_adjust_mod
