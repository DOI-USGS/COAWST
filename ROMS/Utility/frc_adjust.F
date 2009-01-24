#include "cppdefs.h"
      MODULE frc_adjust_mod

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS

# ifdef EW_PERIODIC
#  define IU_RANGE Istr,Iend
#  define IV_RANGE Istr,Iend
# else
#  define IU_RANGE Istr,IendR
#  define IV_RANGE IstrR,IendR
# endif
# ifdef NS_PERIODIC
#  define JU_RANGE Jstr,Jend
#  define JV_RANGE Jstr,Jend
# else
#  define JU_RANGE JstrR,JendR
#  define JV_RANGE Jstr,JendR
# endif
!
!svn $Id: frc_adjust.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are using during the adjustment of surface forcing.  !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: frc_ADadjust
      PUBLIC :: frc_NLadjust
      PUBLIC :: frc_TLadjust
      PUBLIC :: load_NLforcing
      PUBLIC :: load_TLforcing

      CONTAINS

      SUBROUTINE frc_NLadjust (ng, tile, Linp, Lout)
!
!=======================================================================
!                                                                      !
!  This routine adds 4DVAR surface forcing increments to nonlinear     !
!  model fields.  The increments can be constant (Nfrec=1) or time     !
!  interpolated between snapshots (Nfrec>1).                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      Tangent linear state time index to add.                !
!     Lout      Background state time index to save.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Linp, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
# endif
      CALL frc_NLadjust_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
     &                        GRID(ng) % umask,                         &
     &                        GRID(ng) % vmask,                         &
# endif
# ifdef ADJUST_WSTRESS
     &                        FORCES(ng) % tl_ustr,                     &
     &                        FORCES(ng) % tl_vstr,                     &
     &                        FORCES(ng) % ustr,                        &
     &                        FORCES(ng) % vstr,                        &
     &                        FORCES(ng) % sustr,                       &
     &                        FORCES(ng) % svstr,                       &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                        FORCES(ng) % tl_tflux,                    &
     &                        FORCES(ng) % tflux,                       &
     &                        FORCES(ng) % stflx,                       &
# endif
     &                        Linp, Lout)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
# endif
      RETURN
      END SUBROUTINE frc_NLadjust
!
!***********************************************************************
      SUBROUTINE frc_NLadjust_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
# ifdef MASKING
     &                              rmask, umask, vmask,                &
# endif
# ifdef ADJUST_WSTRESS
     &                              tl_ustr,  tl_vstr,                  &
     &                              ustr, vstr,                         &
     &                              sustr, svstr,                       &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                              tl_tflux, tflux, stflx,             &
# endif
     &                              Linp, Lout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Linp, Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(in) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: vstr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
#  endif
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(in) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tflux(LBi:UBi,LBj:UBj,                    &
     &                              Nfrec(ng),2,NT(ng))
      real(r8), intent(in) :: tl_tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
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
        it1=MAX(1,iic(ng)/nSFF(ng))
        it2=MIN(it1+1,Nfrec(ng))
        fac1=SF_time(it2,ng)-time(ng)
        fac2=time(ng)-SF_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF

# ifdef ADJUST_WSTRESS
!
!  If first pass, load background surface momentum stress to store
!  arrays so it can be written into INIname NetCDF file for latter
!  use.  Notice that "sustr" and "svstr" contains the initial
!  interpolated forcing from snapshots.
!
      IF ((time(ng).eq.SF_time(it1,ng))) THEN
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ustr(i,j,it1,Lout)=sustr(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vstr(i,j,it1,Lout)=svstr(i,j)
          END DO
        END DO
      END IF
!
!  Adjust surface momentum stress. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=sustr(i,j)+                                        &
     &               fac1*tl_ustr(i,j,it1,Linp)+                        &
     &               fac2*tl_ustr(i,j,it2,Linp)
#  ifdef MASKING
          sustr(i,j)=sustr(i,j)*umask(i,j)
#  endif
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
          svstr(i,j)=svstr(i,j)+                                        &
     &               fac1*tl_vstr(i,j,it1,Linp)+                        &
     &               fac2*tl_vstr(i,j,it2,Linp)
#  ifdef MASKING
          svstr(i,j)=svstr(i,j)*vmask(i,j)
#  endif
        END DO
      END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  If first pass, load background surface tracers flux to store
!  arrays so it can be written into INIname NetCDF file for latter
!  use.  Notice that "stflx" contain the initial interpolated
!  forcing from snapshots.
!
      IF ((time(ng).eq.SF_time(it1,ng))) THEN
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tflux(i,j,it1,Lout,itrc)=stflx(i,j,itrc)
            END DO
          END DO
        END DO
      END IF
!
!  Adjust surface tracer fluxes. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            stflx(i,j,itrc)=stflx(i,j,itrc)+                            &
     &                      fac1*tl_tflux(i,j,it1,Linp,itrc)+           &
     &                      fac2*tl_tflux(i,j,it2,Linp,itrc)
#  ifdef MASKING
            stflx(i,j,itrc)=stflx(i,j,itrc)*rmask(i,j)
#  endif
          END DO
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE frc_NLadjust_tile

      SUBROUTINE frc_ADadjust (ng, tile, Lout)
!
!=======================================================================
!                                                                      !
!  This routine adds 4DVAR surface forcing adjoint outputs to adjoint  !
!  adjusted fields.  The adjustments can be constant (Nfrec=1) or time !
!  interpolated between snapshots (Nfrec>1).                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Lout      Adjoint state time index were to add.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iADM, 7)
# endif
      CALL frc_ADadjust_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
# ifdef ADJUST_WSTRESS
     &                        FORCES(ng) % ad_ustr,                     &
     &                        FORCES(ng) % ad_vstr,                     &
     &                        FORCES(ng) % ad_sustr,                    &
     &                        FORCES(ng) % ad_svstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                        FORCES(ng) % ad_tflux,                    &
     &                        FORCES(ng) % ad_stflx,                    &
# endif
     &                        Lout)
# ifdef PROFILE
      CALL wclock_off (ng, iADM, 7)
# endif
      RETURN
      END SUBROUTINE frc_ADadjust
!
!***********************************************************************
      SUBROUTINE frc_ADadjust_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
# ifdef ADJUST_WSTRESS
     &                              ad_ustr,  ad_vstr,                  &
     &                              ad_sustr, ad_svstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                              ad_tflux, ad_stflx,                 &
# endif
     &                              Lout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(in) :: ad_svstr(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: ad_stflx(LBi:,LBj:,:)
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
      real(r8), intent(in) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ad_svstr(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: ad_stflx(LBi:UBi,LBj:UBj,NT(ng))
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
        it1=MAX(1,iic(ng)/nSFF(ng))
        it2=MIN(it1+1,Nfrec(ng))
        fac1=SF_time(it2,ng)-time(ng)
        fac2=time(ng)-SF_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF

# ifdef ADJUST_WSTRESS
!
!  Adjoint of : Interpolate between surface forcing
!  increments, if appropriate.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
           ad_ustr(i,j,it1,Lout)=ad_ustr(i,j,it1,Lout)+                 &
     &                           fac1*ad_sustr(i,j)
           ad_ustr(i,j,it2,Lout)=ad_ustr(i,j,it2,Lout)+                 &
     &                           fac2*ad_sustr(i,j)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
           ad_vstr(i,j,it1,Lout)=ad_vstr(i,j,it1,Lout)+                 &
     &                           fac1*ad_svstr(i,j)
           ad_vstr(i,j,it2,Lout)=ad_vstr(i,j,it2,Lout)+                 &
     &                           fac2*ad_svstr(i,j)
        END DO
      END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Adjust surface tracer fluxes. Interpolate between surface forcing
!  increments, if appropriate.
!
      DO itrc=1,NT(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
           ad_tflux(i,j,it1,Lout,itrc)=ad_tflux(i,j,it1,Lout,itrc)+     &
     &                                 fac1*ad_stflx(i,j,itrc)
           ad_tflux(i,j,it2,Lout,itrc)=ad_tflux(i,j,it2,Lout,itrc)+     &
     &                                 fac2*ad_stflx(i,j,itrc)
          END DO
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE frc_ADadjust_tile

      SUBROUTINE frc_TLadjust (ng, tile, Linp)
!
!=======================================================================
!                                                                      !
!  This routine addsinterpolates 4DVAR surface forcing increments.     !
!  The increments can be constant (Nfrec=1) or time                    !
!  interpolated between snapshots (Nfrec>1).                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Linp      Tangent linear state time index to add                 !
!                                                                      !
!=======================================================================
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
      CALL wclock_on (ng, iTLM, 7)
# endif
      CALL frc_TLadjust_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
# ifdef ADJUST_WSTRESS
     &                        FORCES(ng) % tl_ustr,                     &
     &                        FORCES(ng) % tl_vstr,                     &
     &                        FORCES(ng) % tl_sustr,                    &
     &                        FORCES(ng) % tl_svstr,                    &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                        FORCES(ng) % tl_tflux,                    &
     &                        FORCES(ng) % tl_stflx,                    &
# endif
     &                        Linp)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 7)
# endif
      RETURN
      END SUBROUTINE frc_TLadjust
!
!***********************************************************************
      SUBROUTINE frc_TLadjust_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
# ifdef ADJUST_WSTRESS
     &                              tl_ustr,  tl_vstr,                  &
     &                              tl_sustr, tl_svstr,                 &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                              tl_tflux, tl_stflx,                 &
# endif
     &                              Linp)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
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
!  Adjust tangent linear surface forcing fields with 4DVAR increments.
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
        it1=MAX(1,iic(ng)/nSFF(ng))
        it2=MIN(it1+1,Nfrec(ng))
        fac1=SF_time(it2,ng)-time(ng)
        fac2=time(ng)-SF_time(it1,ng)
        fac=1.0_r8/(fac1+fac2)
        fac1=fac*fac1
        fac2=fac*fac2
      END IF

# ifdef ADJUST_WSTRESS
!
!  Adjoint of : Interpolate between surface forcing
!  increments, if appropriate.
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
           tl_sustr(i,j)=fac1*tl_ustr(i,j,it1,Linp)+                    &
     &                   fac2*tl_ustr(i,j,it2,Linp)
        END DO
      END DO
      DO j=Jstr,JendR
        DO i=IstrR,IendR
           tl_svstr(i,j)=fac1*tl_vstr(i,j,it1,Linp)+                    &
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
        DO j=JstrR,JendR
          DO i=IstrR,IendR
             tl_stflx(i,j,itrc)=fac1*tl_tflux(i,j,it1,Linp,itrc)+       &
     &                          fac2*tl_tflux(i,j,it2,Linp,itrc)
          END DO
        END DO
      END DO
# endif
      RETURN
      END SUBROUTINE frc_TLadjust_tile

      SUBROUTINE load_NLforcing (ng, tile, Lout)
!
!=======================================================================
!                                                                      !
!  This routine loads surface forcing into nonlinear storage arrays.   !
!  In  4DVAR  surface forcing adjustment,  the fluxes are stored in    !
!  arrays with extra dimensions to facilitate minimization at other    !
!  times in addition to initialization time.                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Lout      Time index to process in storage arrays.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 8)
# endif
      CALL load_NLforcing_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % sustr,                     &
     &                          FORCES(ng) % svstr,                     &
     &                          FORCES(ng) % ustr,                      &
     &                          FORCES(ng) % vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % stflx,                     &
     &                          FORCES(ng) % tflux,                     &
# endif
     &                          Lout)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 8)
# endif

      RETURN
      END SUBROUTINE load_NLforcing
!
!***********************************************************************
      SUBROUTINE load_NLforcing_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
# ifdef ADJUST_WSTRESS
     &                                sustr, svstr,                     &
     &                                ustr,  vstr,                      &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                stflx, tflux,                     &
# endif
     &                                Lout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx(LBi:,LBj:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:,LBj:,:,:,:)
#  endif
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tflux(LBi:UBi,LBj:UBj,                 &
     &                                 Nfrec(ng),2,NT(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Load nonlinear surface forcing fields into storage arrays.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng)-1,nSFF(ng)).eq.0) THEN
        SFcount(ng)=SFcount(ng)+1

# ifdef ADJUST_WSTRESS
!
!  Load surface momentum stress.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            ustr(i,j,SFcount(ng),Lout)=sustr(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            vstr(i,j,SFcount(ng),Lout)=svstr(i,j)
          END DO
        END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Load surface tracer fluxes.
!
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tflux(i,j,SFcount(ng),Lout,itrc)=stflx(i,j,itrc)
            END DO
          END DO
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE load_NLforcing_tile

      SUBROUTINE load_TLforcing (ng, tile, Lout)
!
!=======================================================================
!                                                                      !
!  This routine loads surface forcing into tangent storage arrays. In  !
!  4DVAR surface forcing adjustment,  the fluxes are stored in arrays  !
!  with extra dimensions to facilitate minimization at other times in  !
!  addition to initialization time.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     tile      Domain partition.                                      !
!     Lout      Time index to process in storage arrays.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lout
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iTLM, 8)
# endif
      CALL load_TLforcing_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
# ifdef ADJUST_WSTRESS
     &                          FORCES(ng) % tl_sustr,                  &
     &                          FORCES(ng) % tl_svstr,                  &
     &                          FORCES(ng) % tl_ustr,                   &
     &                          FORCES(ng) % tl_vstr,                   &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                          FORCES(ng) % tl_stflx,                  &
     &                          FORCES(ng) % tl_tflux,                  &
# endif
     &                          Lout)
# ifdef PROFILE
      CALL wclock_off (ng, iTLM, 8)
# endif

      RETURN
      END SUBROUTINE load_TLforcing
!
!***********************************************************************
      SUBROUTINE load_TLforcing_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
# ifdef ADJUST_WSTRESS
     &                                tl_sustr, tl_svstr,               &
     &                                tl_ustr,  tl_vstr,                &
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
     &                                tl_stflx, tl_tflux,               &
# endif
     &                                Lout)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Lout
!
# ifdef ASSUMED_SHAPE
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(in) :: tl_svstr(LBi:,LBj:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tl_stflx(LBi:,LBj:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_vstr(LBi:,LBj:,:,:)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:,LBj:,:,:,:)
#  endif
# else
#  ifdef ADJUST_WSTRESS
      real(r8), intent(in) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_svstr(LBi:UBi,LBj:UBj)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(in) :: tl_stflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: tl_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: tl_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      real(r8), intent(inout) :: tl_tflux(LBi:UBi,LBj:UBj,              &
     &                                    Nfrec(ng),2,NT(ng))
#  endif
# endif
!
!  Local variable declarations.
!
      integer :: i, j
# ifdef SOLVE3D
      integer :: itrc
# endif

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Load tangent linear surface forcing fields into storage arrays.
!-----------------------------------------------------------------------
!
      IF (MOD(iic(ng)-1,nSFF(ng)).eq.0) THEN
        SFcount(ng)=SFcount(ng)+1

# ifdef ADJUST_WSTRESS
!
!  Load surface momentum stress.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            tl_ustr(i,j,SFcount(ng),Lout)=tl_sustr(i,j)
          END DO
        END DO
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            tl_vstr(i,j,SFcount(ng),Lout)=tl_svstr(i,j)
          END DO
        END DO
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Load surface tracer fluxes.
!
        DO itrc=1,NT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              tl_tflux(i,j,SFcount(ng),Lout,itrc)=tl_stflx(i,j,itrc)
            END DO
          END DO
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE load_TLforcing_tile

# undef IU_RANGE
# undef IV_RANGE
# undef JU_RANGE
# undef JV_RANGE
#endif
      END MODULE frc_adjust_mod
