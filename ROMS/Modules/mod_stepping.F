#include "cppdefs.h"
      MODULE mod_stepping
!
!svn $Id: mod_stepping.F 1404 2008-02-20 01:43:42Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains time stepping indices.                         !
!                                                                      !
#ifdef BOUNDARY
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

        implicit none

        integer, private :: ig

        integer, dimension(Ngrids) :: knew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: krhs = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: kstp = (/ (1, ig=1,Ngrids) /)

        integer, dimension(Ngrids) :: nnew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nrhs = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nstp = (/ (1, ig=1,Ngrids) /)

#ifdef FLOATS
        integer, dimension(Ngrids) :: nf   = (/ (0, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfp1 = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfm3 = (/ (2, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfm2 = (/ (3, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: nfm1 = (/ (4, ig=1,Ngrids) /)
#endif

#ifdef ICE_MODEL
        integer, dimension(Ngrids) :: linew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: liold = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: liunw = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: liuol = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: lienw = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: lieol = (/ (1, ig=1,Ngrids) /)
#endif

#if defined SSH_TIDES || defined UV_TIDES || defined POT_TIDES
        integer, dimension(Ngrids) :: NTC
#endif
#ifdef ADJUST_BOUNDARY
        integer, dimension(Ngrids) :: Lbinp = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: Lbout = (/ (1, ig=1,Ngrids) /)
#endif
#if defined ADJUST_BOUNDARY || defined ADJUST_STFLUX || \
    defined ADJUST_WSTRESS
        integer, dimension(Ngrids) :: Lfinp = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: Lfout = (/ (1, ig=1,Ngrids) /)
#endif
        integer, dimension(Ngrids) :: Lnew = (/ (1, ig=1,Ngrids) /)
        integer, dimension(Ngrids) :: Lold = (/ (1, ig=1,Ngrids) /)

      END MODULE mod_stepping

