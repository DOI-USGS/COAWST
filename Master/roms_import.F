#include "cppdefs.h"
      MODULE ROMS_import_mod

#ifdef ROMS_COUPLING
!
!svn $Id: roms_import.F 594 2008-04-01 18:11:31Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to interpolate imported       !
!  to particular model grid during coupling.                           !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      PUBLIC :: ROMS_import2d

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ROMS_import2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
# if defined MCT_LIB
     &                          Npts, InpField,                         &
# elif defined ESMF_LIB
     &                          InpField,                               &
# endif
     &                          Imin, Imax, Jmin, Jmax,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OutFmin, OutFmax,                       &
     &                          OutField,                               &
     &                          status)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE distribute_mod,  ONLY : mp_reduce
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, id, gtype
      integer, intent(in) :: Imin, Imax, Jmin, Jmax
      integer, intent(in) :: LBi, UBi, LBj, UBj
# ifdef MCT_LIB
      integer, intent(in) :: Npts
# endif
      integer, intent(out) :: status

      real(r8), intent(in) :: scale, add_offset
      real(r8), intent(out) :: OutFmin, OutFmax

# ifdef ASSUMED_SHAPE
#  if defined MCT_LIB
      real(r8), intent(in) ::  InpField(:)
#  elif defined ESMF_LIB
      real(r8), intent(in) ::  InpField(:,:)
#  endif
      real(r8), intent(out) :: OutField(LBi:,LBj:)
# else
#  if defined MCT_LIB
      real(r8), intent(in) ::  InpField(Npts)
#  elif defined ESMF_LIB
      real(r8), intent(in) ::  InpField(:,:)
#  endif
      real(r8), intent(out) :: OutField(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, ij, j

      real(r8), parameter :: Large = 1.0E+20_r8

      real(r8), dimension(2) :: range
# ifdef DISTRIBUTE
      character (len=3), dimension(2) :: op_handle
# endif
!
!-----------------------------------------------------------------------
!  Import 2D field.
!-----------------------------------------------------------------------
!
      status=0
      range(1)= Large
      range(2)=-Large

# if defined MCT_LIB
!
!  For now couple models grids are identical so no interpolation is
!  necessary. Interpolation logic will be provided in the future.
!
      ij=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          ij=ij+1
          OutField(i,j)=scale*InpField(ij)+add_offset
          range(1)=MIN(range(1),OutField(i,j))
          range(2)=MAX(range(2),OutField(i,j))
        END DO
      END DO
# elif defined ESMF_LIB
# endif
# ifdef DISTRIBUTE
!
!  Global reduction for imported field range values.
!
      op_handle(1)='MIN'
      op_handle(2)='MAX'
      CALL mp_reduce (ng, iNLM, 2, range, op_handle)
      OutFmin=range(1)
      OutFmax=range(2)
# endif
!
!-----------------------------------------------------------------------
!  Exchange boundary information.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        IF (gtype.eq.r2dvar) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OutField)
        ELSE IF (gtype.eq.u2dvar) THEN
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OutField)
        ELSE IF (gtype.eq.v2dvar) THEN
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            OutField)
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    OutField)
# endif

      END SUBROUTINE ROMS_import2d
#endif
      END MODULE ROMS_import_mod
