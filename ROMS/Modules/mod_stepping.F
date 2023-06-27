#include "cppdefs.h"
      MODULE mod_stepping
!
!git $Id$
!svn $Id: mod_stepping.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains time stepping indices.                         !
!                                                                      !
#ifdef ADJUST_BOUNDARY
!  Lbinp     Open boundary adjustment input  fields index.             !
!  Lbout     Open boundary adjustment output fields index.             !
#endif
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
!  Lfinp     Surface forcing adjustment input  fields index.           !
!  Lfout     Surface forcing adjustment output fields index.           !
#endif
!  Lnew      New descent algorithm state solution index.               !
!  Lold      Previous descent algorithm state solution index.          !
!                                                                      !
!  knew      Barotropic (fast) time-step index corresponding to the    !
!              newest values for 2D primitive equation variables.      !
!  krhs      Barotropic (fast) time-step index used to compute the     !
!              right-hand-terms of 2D primitive equation variables.    !
!  kstp      Barotropic (fast) time-step index to which the current    !
!              changes are added to compute new 2D primitive equation  !
!              variables.                                              !
!                                                                      !
!  nfm3      Float index for time level "n-3".                         !
!  nfm2      Float index for time level "n-2".                         !
!  nfm1      Float index for time level "n-1".                         !
!  nf        Float index for time level "n".                           !
!  nfp1      Float index for time level "n+1".                         !
!                                                                      !
!  nnew      Baroclinic (slow) time-step index corresponding to the    !
!              newest values for 3D primitive equation variables.      !
!  nrhs      Baroclinic (slow) time-step index used to compute the     !
!              right-hand-terms of 3D primitive equation variables.    !
!  nstp      Baroclinic (slow) time-step index to which the current    !
!              changes are added to compute new 3D primitive equation  !
!              variables.                                              !
#if defined SSH_TIDES || defined UV_TIDES
!                                                                      !
!  NTC       Number of tidal components to consider.                   !
#endif
!                                                                      !
!=======================================================================
!
!       USE mod_param
!
        implicit none
!
        PUBLIC :: allocate_stepping
        PUBLIC :: deallocate_stepping
!
!-----------------------------------------------------------------------
!  Define module variables.
!-----------------------------------------------------------------------
!
        integer, allocatable :: knew(:)
        integer, allocatable :: krhs(:)
        integer, allocatable :: kstp(:)
!$OMP THREADPRIVATE (knew, krhs, kstp)
!
        integer, allocatable :: nnew(:)
        integer, allocatable :: nrhs(:)
        integer, allocatable :: nstp(:)
!$OMP THREADPRIVATE (nnew, nrhs, nstp)

#ifdef FLOATS
!
        integer, allocatable :: nf(:)
        integer, allocatable :: nfp1(:)
        integer, allocatable :: nfm3(:)
        integer, allocatable :: nfm2(:)
        integer, allocatable :: nfm1(:)
!$OMP THREADPRIVATE (nf, nfp1, nfm3, nfm2, nfm1)
#endif
!
#ifdef ADJUST_BOUNDARY
        integer, allocatable :: Lbinp(:)
        integer, allocatable :: Lbout(:)
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
        integer, allocatable :: Lfinp(:)
        integer, allocatable :: Lfout(:)
#endif

        integer, allocatable :: Lnew(:)
        integer, allocatable :: Lold(:)

#if defined SSH_TIDES || defined UV_TIDES
        integer, allocatable :: NTC(:)
#endif
!
      CONTAINS
!
      SUBROUTINE allocate_stepping (Ngrids)
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
       integer, intent(in) :: Ngrids
!
!-----------------------------------------------------------------------
!  Allocate and intialize time indices.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
      IF (.not.allocated(knew)) THEN
        allocate ( knew(Ngrids) )
      END IF
      knew(1:Ngrids)=1

      IF (.not.allocated(krhs)) THEN
        allocate ( krhs(Ngrids) )
      END IF
      krhs(1:Ngrids)=1

      IF (.not.allocated(kstp)) THEN
        allocate ( kstp(Ngrids) )
      END IF
      kstp(1:Ngrids)=1

      IF (.not.allocated(nnew)) THEN
        allocate ( nnew(Ngrids) )
      END IF
      nnew(1:Ngrids)=1

      IF (.not.allocated(nrhs)) THEN
        allocate ( nrhs(Ngrids) )
      END IF
      nrhs(1:Ngrids)=1

      IF (.not.allocated(nstp)) THEN
        allocate ( nstp(Ngrids) )
      END IF
      nstp(1:Ngrids)=1

#ifdef FLOATS
      IF (.not.allocated(nf)) THEN
        allocate ( nf(Ngrids) )
      END IF
      nf(1:Ngrids)=0

      IF (.not.allocated(nfp1)) THEN
        allocate ( nfp1(Ngrids) )
      END IF
      nfp1(1:Ngrids)=1

      IF (.not.allocated(nfm3)) THEN
        allocate ( nfm3(Ngrids) )
      END IF
      nfm3(1:Ngrids)=2

      IF (.not.allocated(nfm2)) THEN
        allocate ( nfm2(Ngrids) )
      END IF
      nfm2(1:Ngrids)=3

      IF (.not.allocated(nfm1)) THEN
        allocate ( nfm1(Ngrids) )
      END IF
      nfm1(1:Ngrids)=4
#endif
!$OMP END PARALLEL

#ifdef ADJUST_BOUNDARY
      IF (.not.allocated(Lbinp)) THEN
        allocate ( Lbinp(Ngrids) )
      END IF
      Lbinp(1:Ngrids)=1

      IF (.not.allocated(Lbout)) THEN
        allocate ( Lbout(Ngrids) )
      END IF
      Lbout(1:Ngrids)=1
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
      IF (.not.allocated(Lfinp)) THEN
        allocate ( Lfinp(Ngrids) )
      END IF
      Lfinp(1:Ngrids)=1

      IF (.not.allocated(Lfout)) THEN
        allocate ( Lfout(Ngrids) )
      END IF
      Lfout(1:Ngrids)=1
#endif

      IF (.not.allocated(Lnew)) THEN
        allocate ( Lnew(Ngrids) )
      END IF
      Lnew(1:Ngrids)=1

      IF (.not.allocated(Lold)) THEN
        allocate ( Lold(Ngrids) )
      END IF
      Lold(1:Ngrids)=1

#if defined SSH_TIDES || defined UV_TIDES
      IF (.not.allocated(NTC)) THEN
        allocate ( NTC(Ngrids) )
      END IF
#endif
!
      RETURN
      END SUBROUTINE allocate_stepping
!
      SUBROUTINE deallocate_stepping
!
!=======================================================================
!                                                                      !
!  This routine deallocates several variables in the module that       !
!  depend on the number of nested grids.                               !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Deallocate variables in module.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
      IF (allocated(knew))  deallocate ( knew )

      IF (allocated(krhs))  deallocate ( krhs )

      IF (allocated(kstp))  deallocate ( kstp )

      IF (allocated(nnew))  deallocate ( nnew )

      IF (allocated(nrhs))  deallocate ( nrhs )

      IF (allocated(nstp))  deallocate ( nstp )

#ifdef FLOATS
      IF (allocated(nf))    deallocate ( nf )

      IF (allocated(nfp1))  deallocate ( nfp1 )

      IF (allocated(nfm3))  deallocate ( nfm3 )

      IF (allocated(nfm2))  deallocate ( nfm2 )

      IF (allocated(nfm1))  deallocate ( nfm1 )
#endif
!$OMP END PARALLEL

#ifdef ADJUST_BOUNDARY
      IF (allocated(Lbinp)) deallocate ( Lbinp )

      IF (allocated(Lbout)) deallocate ( Lbout )
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
      IF (allocated(Lfinp)) deallocate ( Lfinp )

      IF (allocated(Lfout)) deallocate ( Lfout )
#endif

      IF (allocated(Lnew))  deallocate ( Lnew )

      IF (allocated(Lold))  deallocate ( Lold )

#if defined SSH_TIDES || defined UV_TIDES
      IF (allocated(NTC))   deallocate ( NTC )
#endif
!
      RETURN
      END SUBROUTINE deallocate_stepping
!
      END MODULE mod_stepping

