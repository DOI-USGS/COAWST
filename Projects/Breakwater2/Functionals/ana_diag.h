      SUBROUTINE ana_diag (ng, tile, model)
!
!! svn $Id: ana_diag.h 34 2007-04-27 04:40:21Z arango $
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine is provided so the USER can compute any specialized    !
!  diagnostics.  If activated, this routine is call at end of every    !
!  3D-equations timestep.                                              !
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
      CALL ana_diag_tile (ng, model, Istr, Iend, Jstr, Jend,            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    OCEAN(ng) % ubar,                             &
     &                    OCEAN(ng) % vbar,                             &
     &                    OCEAN(ng) % u,                                &
     &                    OCEAN(ng) % v)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        WRITE (ANANAME( 5),'(a,a)') TRIM(Adir), '/ana_diag.h'
      END IF

      RETURN
      END SUBROUTINE ana_diag
!
!***********************************************************************
      SUBROUTINE ana_diag_tile (ng, model, Istr, Iend, Jstr, Jend,      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          ubar, vbar, u, v)
!***********************************************************************
!
      USE mod_param
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Iend, Istr, Jend, Jstr
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#else
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#endif
!
!  Local variable declarations.
!
      integer :: i, j, k
!
!-----------------------------------------------------------------------
!  Report any user diagnostics
!-----------------------------------------------------------------------
!
!  Open USER file.
!
      IF (iic(ng).eq.ntstart(ng)) THEN
        OPEN (usrout,file=USRname,form='formatted',status='unknown',    &
     &        err=40)
        GO TO 60
  40    WRITE (stdout,50) USRname
  50    FORMAT (' ANA_DIAG - unable to open output file: ',a)
        exit_flag=2
  60    CONTINUE
      END IF
!
!  Write out maximum values of velocity.
!
      WRITE (usrout,70) 'No user diagnostics computed.'
  70  FORMAT (a)

      RETURN
      END SUBROUTINE ana_diag_tile
