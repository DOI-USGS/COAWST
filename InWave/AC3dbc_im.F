#include "cppdefs.h"
      MODULE AC3dbc_mod
#ifdef INWAVE_COUPLING
!
!svn $Id: t3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!                                                                      !
!  This subroutine sets lateral boundary conditions for the ITRC-th    !
!  tracer field.                                                       !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: AC3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE AC3dbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
      USE mod_inwave_vars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL AC3dbc_tile (ng, tile,                                       &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 nout,                                            &
     &                 WAVEP(ng)% AC)
      RETURN
      END SUBROUTINE AC3dbc
!
!***********************************************************************
      SUBROUTINE AC3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       nout,                                      &
     &                       AC)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      USE mod_inwave_params
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: nout
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: AC(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: AC(LBi:UBi,LBj:UBj,ND(ng),3)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d

      real(r8), parameter :: eps =1.0E-20_r8

      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau

      real(r8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: grad

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN

# if defined WEST_ACGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO d=1,ND(ng)
          DO j=Jstr,Jend
            AC(Istr-1,j,d,nout)=AC(Istr,j,d,nout)
#  ifdef MASKING
            AC(Istr-1,j,d,nout)=AC(Istr-1,j,d,nout)*                    &
     &                              GRID(ng)%rmask(Istr-1,j)
#  endif
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (EASTERN_EDGE) THEN
# if defined EAST_ACGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO d=1,ND(ng)
          DO j=Jstr,Jend
            AC(Iend+1,j,d,nout)=AC(Iend,j,d,nout)
#  ifdef MASKING
            AC(Iend+1,j,d,nout)=AC(Iend+1,j,d,nout)*                     &
     &                              GRID(ng)%rmask(Iend+1,j)
#  endif
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (SOUTHERN_EDGE) THEN
# if defined SOUTH_ACGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO d=1,ND(ng)
          DO i=Istr,Iend
            AC(i,Jstr-1,d,nout)=AC(i,Jstr,d,nout)
#  ifdef MASKING
            AC(i,Jstr-1,d,nout)=AC(i,Jstr-1,d,nout)*                     &
     &                              GRID(ng)%rmask(i,Jstr-1)
#  endif
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (NORTHERN_EDGE) THEN
# if defined NORTH_ACGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO d=1,ND(ng)
          DO i=Istr,Iend
            AC(i,Jend+1,d,nout)=AC(i,Jend,d,nout)
#  ifdef MASKING
            AC(i,Jend+1,d,nout)=AC(i,Jend+1,d,nout)*                    &
     &                              GRID(ng)%rmask(i,Jend+1)
#  endif
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF ((SOUTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND(ng)
          AC(Istr-1,Jstr-1,d,nout)=0.5_r8*                              &
     &                                 (AC(Istr  ,Jstr-1,d,nout)+       &
     &                                  AC(Istr-1,Jstr  ,d,nout))
        END DO
      END IF
      IF ((SOUTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND(ng)
          AC(Iend+1,Jstr-1,d,nout)=0.5_r8*                              &
     &                                 (AC(Iend  ,Jstr-1,d,nout)+       &
     &                                  AC(Iend+1,Jstr  ,d,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(WESTERN_EDGE)) THEN
        DO d=1,ND(ng)
          AC(Istr-1,Jend+1,d,nout)=0.5_r8*                              &
     &                                 (AC(Istr-1,Jend  ,d,nout)+       &
     &                                  AC(Istr  ,Jend+1,d,nout))
        END DO
      END IF
      IF ((NORTHERN_EDGE).and.(EASTERN_EDGE)) THEN
        DO d=1,ND(ng)
          AC(Iend+1,Jend+1,d,nout)=0.5_r8*                              &
     &                                 (AC(Iend+1,Jend  ,d,nout)+       &
     &                                  AC(Iend  ,Jend+1,d,nout))
        END DO
      END IF
      RETURN
      END SUBROUTINE AC3dbc_tile
#endif
      END MODULE AC3dbc_mod
