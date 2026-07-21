#include "cppdefs.h"
      MODULE ct3dbc_dir_mod
#ifdef SOLVE3D
!
!svn $Id: t3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: ct3dbc_dir_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE ct3dbc_dir (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_inwave_vars
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL ct3dbc_dir_tile (ng, tile,                                   &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 WAVEP(ng)% ct)
      RETURN
      END SUBROUTINE ct3dbc_dir

!
!***********************************************************************
      SUBROUTINE ct3dbc_dir_tile (ng, tile,                             &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       ct)
!***********************************************************************
!
      USE mod_param
      USE mod_inwave_params
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
# ifdef REFINED_GRID
      USE mod_stepping
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj

!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ct(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ct(LBi:UBi,LBj:UBj,ND+1)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d


# include "set_bounds.h"

      DO j=Jstr,Jend
        DO i=Istr,Iend
#  if defined UP_CT_GRADIENT
          ct(i,j,ND+1)=ct(i,j,ND)
#  elif defined UP_CT_WALL
          ct(i,j,ND+1)=0.0_r8
#  endif
#  if defined DOWN_CT_GRADIENT
          ct(i,j,1)=ct(i,j,2)
#  elif defined DOWN_CT_WALL
          ct(i,j,1)=0.0_r8
#  endif
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ct3dbc_dir_tile
#endif
      END MODULE ct3dbc_dir_mod
