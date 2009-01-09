      SUBROUTINE ana_passive (ng, tile, model)
!
!! svn $Id: ana_passive.h 38 2007-04-28 01:11:25Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for passive inert tracers      !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_ocean
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_passive_tile (ng, model, Istr, Iend, Jstr, Jend,         &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       OCEAN(ng) % t)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME(18),'(a,a)') TRIM(Adir), '/ana_passive.h'
      END IF

      RETURN
      END SUBROUTINE ana_passive
!
!***********************************************************************
      SUBROUTINE ana_passive_tile (ng, model, Istr, Iend, Jstr, Jend,   &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(out) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
      integer :: i, itrc, j, k

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set analytical initial conditions for passive inert tracers.
!-----------------------------------------------------------------------
!
#if defined MY_APPLICATION
      DO ip=1,NPT
        itrc=inert(ip)        
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              t(i,j,k,1,itrc)=???
              t(i,j,k,2,itrc)=t(i,j,k,1,itrc)
            END DO
          END DO
        END DO
      END DO
#else
      ana_passive.h: No values provided for passive tracers.
#endif

      RETURN
      END SUBROUTINE ana_passive_tile
