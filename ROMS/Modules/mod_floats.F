#include "cppdefs.h"
      MODULE mod_floats
#ifdef FLOATS
!
!svn $Id: mod_floats.F 737 2008-09-07 02:06:44Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Findex     Indices of spherical coordinates entries in initial      !
!               location arrays, if any.                               !
!  Flon       Initial longitude locations, if any.                     !
!  Flat       Initial latitude locations, if any.                      !
!  Ftype      Float trajectory type:                                   !
!               Ftype(:) = 1,  neutral density 3D Lagrangian           !
!               Ftype(:) = 2,  isobaric (constant depth) float.        !
!  Tinfo      Float trajectory initial information.                    !
!  bounded    Float bounded status switch.                             !
!  track      Multivariate float trajectory data at several time       !
!               time levels.                                           !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

        TYPE T_FLT

          logical, pointer  :: bounded(:)

          integer, pointer :: Findex(:)
          integer, pointer :: Ftype(:)

          real(r8), pointer :: Flon(:)
          real(r8), pointer :: Flat(:)
          real(r8), pointer :: Tinfo(:,:)
          real(r8), pointer :: track(:,:,:)
          real(r8), pointer :: Fz0(:)

        END TYPE T_FLT

        TYPE (T_FLT), allocatable :: FLT(:)

      CONTAINS

      SUBROUTINE allocate_floats (ng)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initialize all variables in the module   !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, iflt

      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FLT(Ngrids) )
!
      allocate ( FLT(ng) % bounded(Nfloats(ng)) )

      allocate ( FLT(ng) % Findex(0:Nfloats(ng)) )

      allocate ( FLT(ng) % Ftype(Nfloats(ng)) )

      allocate ( FLT(ng) % Flon(Nfloats(ng)) )

      allocate ( FLT(ng) % Flat(Nfloats(ng)) )

      allocate ( FLT(ng) % Tinfo(0:izrhs,Nfloats(ng)) )

      allocate ( FLT(ng) % track(NFV(ng),0:NFT,Nfloats(ng)) )

      allocate ( FLT(ng) % Fz0(Nfloats(ng)) )
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      FLT(ng) % Findex(0) = 0
      DO iflt=1,Nfloats(ng)
        FLT(ng) % bounded(iflt) = .FALSE.
        FLT(ng) % Findex(iflt) = 0
        FLT(ng) % Ftype(iflt) = 0
        FLT(ng) % Flon(iflt) = IniVal
        FLT(ng) % Flat(iflt) = IniVal
        FLT(ng) % Fz0(iflt) = 0
        DO i=0,izrhs
          FLT(ng) % Tinfo(i,iflt) = IniVal
        END DO
        DO i=1,NFV(ng)
          FLT(ng) % track(i,0,iflt) = IniVal
          FLT(ng) % track(i,1,iflt) = IniVal
          FLT(ng) % track(i,2,iflt) = IniVal
          FLT(ng) % track(i,3,iflt) = IniVal
          FLT(ng) % track(i,4,iflt) = IniVal
        END DO
      END DO

      RETURN
      END SUBROUTINE allocate_floats
#endif
      END MODULE mod_floats
