#include "cppdefs.h"
      MODULE mod_stepping
!
!svn $Id: mod_stepping.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
        USE mod_param
!
        implicit none
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

#ifdef ICE_MODEL
        integer, allocatable :: linew(:)
        integer, allocatable :: liold(:)
        integer, allocatable :: liunw(:)
        integer, allocatable :: liuol(:)
        integer, allocatable :: lienw(:)
        integer, allocatable :: lieol(:)
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

#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
        integer, allocatable :: NTC(:)
#endif
!
      CONTAINS
!
      SUBROUTINE allocate_stepping
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate and intialize time indices.
!-----------------------------------------------------------------------
!
!$OMP PARALLEL
      allocate ( knew(Ngrids) )
      knew(1:Ngrids)=1

      allocate ( krhs(Ngrids) )
      krhs(1:Ngrids)=1

      allocate ( kstp(Ngrids) )
      kstp(1:Ngrids)=1

      allocate ( nnew(Ngrids) )
      nnew(1:Ngrids)=1

      allocate ( nrhs(Ngrids) )
      nrhs(1:Ngrids)=1

      allocate ( nstp(Ngrids) )
      nstp(1:Ngrids)=1

#ifdef FLOATS
      allocate ( nf(Ngrids) )
      nf(1:Ngrids)=0

      allocate ( nfp1(Ngrids) )
      nfp1(1:Ngrids)=1

      allocate ( nfm3(Ngrids) )
      nfm3(1:Ngrids)=2

      allocate ( nfm2(Ngrids) )
      nfm2(1:Ngrids)=3

      allocate ( nfm1(Ngrids) )
      nfm1(1:Ngrids)=4
#endif
!$OMP END PARALLEL

#ifdef ICE_MODEL
      allocate ( linew(Ngrids) )
      linew(1:Ngrids) = 1
      allocate ( liold(Ngrids) )
      liold(1:Ngrids) = 1
      allocate ( liunw(Ngrids) )
      liunw(1:Ngrids) = 1
      allocate ( liuol(Ngrids) )
      liuol(1:Ngrids) = 1
      allocate ( lienw(Ngrids) )
      lienw(1:Ngrids) = 1
      allocate ( lieol(Ngrids) )
      lieol(1:Ngrids) = 1
#endif

#ifdef ADJUST_BOUNDARY
      allocate ( Lbinp(Ngrids) )
      Lbinp(1:Ngrids)=1

      allocate ( Lbout(Ngrids) )
      Lbout(1:Ngrids)=1
#endif

#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
      allocate ( Lfinp(Ngrids) )
      Lfinp(1:Ngrids)=1

      allocate ( Lfout(Ngrids) )
      Lfout(1:Ngrids)=1
#endif

      allocate ( Lnew(Ngrids) )
      Lnew(1:Ngrids)=1

      allocate ( Lold(Ngrids) )
      Lold(1:Ngrids)=1

#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
      allocate ( NTC(Ngrids) )
#endif

      RETURN
      END SUBROUTINE allocate_stepping

      END MODULE mod_stepping

