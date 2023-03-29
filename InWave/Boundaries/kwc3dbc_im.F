#include "cppdefs.h"
      MODULE kwc3dbc_mod
#ifdef SOLVE3D
!
!svn $Id: kwc3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!                                                                      !
!  This subroutine sets lateral boundary conditions for the            !
!  wave number kwc field.                                              !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: kwc3dbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE kwc3dbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         kwc)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
      USE mod_inwave_params
      USE mod_inwave_swan
      USE mod_inwave_vars
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
      real(r8), intent(inout) :: kwc(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: kwc(LBi:UBi,LBj:UBj,ND)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, d

      real(r8) :: twopi, otwopi
      real(r8) :: error, Tr_min
      real(r8) :: L0, k0, k1, kh, wr
      real(r8) :: F, FDER, tanhkh
      real(r8) :: cff, cff1, cff2
      real(r8), parameter :: maxErr = 0.1_r8
!     real(r8), parameter :: kwc_max = 10.0_r8
      real(r8), parameter :: kwc_max = 1.0_r8
      real(r8), parameter :: kwc_min = 0.015_r8

# include "set_bounds.h"

      twopi=2.0_r8*pi
      otwopi=1.0_r8/twopi
      Tr_min=1.0_r8

!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng)) THEN
        IF (DOMAIN(ng)%Western_Edge(tile)) THEN
          IF (LBC(iwest,isAC3d,ng)%acquire) THEN
            DO d=1,ND
              DO j=Jstr,Jend
                L0=g*otwopi*MAX(Tr_min,WAVEP(ng)%Tr(Istr-1,j,d))**2.0_r8
                k0=twopi/L0
                error=100.0_r8
                wr=twopi/MAX(Tr_min,WAVEP(ng)%Tr(Istr-1,j,d))
                DO WHILE(error.gt.maxErr)
                  kh=k0*WAVEP(ng)%h_tot(Istr-1,j)
                  tanhkh=TANH(kh)
                  cff1=wr**2.0_r8
                  cff2=-g*k0*tanhkh
                  F=cff1+cff2
                  cff1=-g*tanhkh
                  cff2=-g*kh/COSH(kh)**2.0_r8
                  FDER=cff1+cff2
                  k1=k0-F/FDER
                  error=100.0_r8*ABS((k1-k0)/k0)
                  k0=k1
                END DO
                kwc(Istr-1,j,d)=MAX(kwc_min,MIN(k0,kwc_max))
!#   ifdef MASKING
!                kwc(Istr-1,j,d)=kwc(Istr-1,j,d)*                        &
!     &                          GRID(ng)%rmask(Istr-1,j)
!#   endif
              END DO
            END DO
          ELSE
            DO d=1,ND
              DO j=Jstr,Jend
                kwc(Istr-1,j,d)=kwc(Istr,j,d)
!#   ifdef MASKING
!                kwc(Istr-1,j,d)=kwc(Istr-1,j,d)*                        &
!     &                          GRID(ng)%rmask(Istr-1,j)
!#   endif
              END DO
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
          IF (LBC(ieast,isAC3d,ng)%acquire) THEN
            DO d=1,ND
              DO j=Jstr,Jend
                L0=g*otwopi*MAX(Tr_min,WAVEP(ng)%Tr(Iend+1,j,d))**2.0_r8
                k0=twopi/L0
                error=100.0_r8
                wr=twopi/MAX(Tr_min,WAVEP(ng)%Tr(Iend+1,j,d))
                DO WHILE(error.gt.maxErr)
                  kh=k0*WAVEP(ng)%h_tot(Iend+1,j)
                  tanhkh=TANH(kh)
                  cff1=wr**2.0_r8
                  cff2=-g*k0*tanhkh
                  F=cff1+cff2
                  cff1=-g*tanhkh
                  cff2=-g*kh/COSH(kh)**2.0_r8
                  FDER=cff1+cff2
                  k1=k0-F/FDER
                  error=100.0_r8*ABS((k1-k0)/k0)
                  k0=k1
                END DO
                kwc(Iend+1,j,d)=MAX(kwc_min,MIN(k0,kwc_max))
!#   ifdef MASKING
!                kwc(Iend+1,j,d)=kwc(Iend+1,j,d)*                        &
!     &                          GRID(ng)%rmask(Iend+1,j)
!#   endif
              END DO
            END DO
          ELSE
            DO d=1,ND
              DO j=Jstr,Jend
                kwc(Iend+1,j,d)=kwc(Iend,j,d)
!#   ifdef MASKING
!                kwc(Iend+1,j,d)=kwc(Iend+1,j,d)*                        &
!     &                         GRID(ng)%rmask(Iend+1,j)
!#   endif
              END DO
            END DO
          END IF
        END IF
      END IF
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (.not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
          IF (LBC(isouth,isAC3d,ng)%acquire) THEN
            DO d=1,ND
              DO i=Istr,Iend
!               L0=g*otwopi*MAX(Tr_min,WAVEP(ng)%Tr(i,Jstr-1,d))**2.0_r8
                L0=g*otwopi*MAX(Tr_min,WAVEG(ng)%Trep)**2.0_r8
                k0=twopi/L0
                error=100.0_r8
!               wr=twopi/MAX(Tr_min,WAVEP(ng)%Tr(i,Jstr-1,d))
                wr=twopi/MAX(Tr_min,WAVEG(ng)%Trep)
                DO WHILE(error.gt.maxErr)
                  kh=k0*WAVEP(ng)%h_tot(i,Jstr-1)
                  tanhkh=TANH(kh)
                  cff1=wr**2.0_r8
                  cff2=-g*k0*tanhkh
                  F=cff1+cff2
                  cff1=-g*tanhkh
                  cff2=-g*kh/COSH(kh)**2.0_r8
                  FDER=cff1+cff2
                  k1=k0-F/FDER
                  error=100.0_r8*ABS((k1-k0)/k0)
                  k0=k1
                END DO
                kwc(i,Jstr-1,d)=MAX(kwc_min,MIN(k0,kwc_max))
!#   ifdef MASKING
!                kwc(i,Jstr-1,d)=kwc(i,Jstr-1,d)*                        &
!     &                          GRID(ng)%rmask(i,Jstr-1)
!#   endif
              END DO
            END DO
          ELSE
            DO d=1,ND
              DO i=Istr,Iend
                kwc(i,Jstr-1,d)=kwc(i,Jstr,d)
!#   ifdef MASKING
!                kwc(i,Jstr-1,d)=kwc(i,Jstr-1,d)*                        &
!     &                         GRID(ng)%rmask(i,Jstr-1)
!#   endif
              END DO
            END DO
          END IF
        END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
        IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
          IF (LBC(inorth,isAC3d,ng)%acquire) THEN
            DO d=1,ND
              DO i=Istr,Iend
                L0=g*otwopi*MAX(Tr_min,WAVEP(ng)%Tr(i,Jend+1,d))**2.0_r8
                k0=twopi/L0
                error=100.0_r8
                wr=twopi/MAX(Tr_min,WAVEP(ng)%Tr(i,Jend+1,d))
                DO WHILE(error.gt.maxErr)
                  kh=k0*WAVEP(ng)%h_tot(i,Jend+1)
                  tanhkh=TANH(kh)
                  cff1=wr**2.0_r8
                  cff2=-g*k0*tanhkh
                  F=cff1+cff2
                  cff1=-g*tanhkh
                  cff2=-g*kh/COSH(kh)**2.0_r8
                  FDER=cff1+cff2
                  k1=k0-F/FDER
                  error=100.0_r8*ABS((k1-k0)/k0)
                  k0=k1
                END DO
                kwc(i,Jend+1,d)=MAX(kwc_min,MIN(k0,kwc_max))
!#   ifdef MASKING
!                kwc(i,Jend+1,d)=kwc(i,Jend+1,d)*                        &
!     &                          GRID(ng)%rmask(i,Jend+1)
!#   endif
              END DO
            END DO
          ELSE
            DO d=1,ND
              DO i=Istr,Iend
                kwc(i,Jend+1,d)=kwc(i,Jend,d)
!#   ifdef MASKING
!                kwc(i,Jend+1,d)=kwc(i,Jend+1,d)*                          &
!     &                          GRID(ng)%rmask(i,Jend+1)
!#   endif
              END DO
            END DO
          END IF
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO d=1,ND
              kwc(Istr-1,Jstr-1,d)=0.5_r8*                              &
     &                            (kwc(Istr,Jstr-1,d)+                  &
     &                             kwc(Istr-1,Jstr,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO d=1,ND
              kwc(Iend+1,Jstr-1,d)=0.5_r8*                              &
     &                            (kwc(Iend  ,Jstr-1,d)+                &
     &                             kwc(Iend+1,Jstr  ,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO d=1,ND
              kwc(Istr-1,Jend+1,d)=0.5_r8*                              &
     &                            (kwc(Istr-1,Jend  ,d)+                &
     &                             kwc(Istr  ,Jend+1,d))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO d=1,ND
              kwc(Iend+1,Jend+1,d)=0.5_r8*                              &
     &                            (kwc(Iend+1,Jend  ,d)+                &
     &                             kwc(Iend  ,Jend+1,d))
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE kwc3dbc_tile
#endif
      END MODULE kwc3dbc_mod
