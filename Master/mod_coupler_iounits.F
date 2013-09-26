#include "cppdefs.h"
      MODULE mod_coupler_iounits
!
!svn $Id: mod_iounits.F 755 2008-09-14 19:07:08Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  Aname       Atmosphere model stadard input file name.               !
!  IWBNDname   Input boundary data file name for InWave model          !
!  IWSWNname   Input spectral SWAN data file name for InWave model     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
#ifdef SWAN_MODEL
      USE swan_iounits
#endif
#ifdef ROMS_MODEL
      USE mod_param
#endif
      implicit none

#ifndef ROMS_MODEL
      integer, parameter :: stdinp = 5
      integer, parameter :: stdout = 6
#endif

#ifdef ROMS_COUPLING
      character (len=80), dimension(Ngrids) :: O2Aname
      character (len=80), dimension(Ngrids) :: O2Wname
      character (len=80), dimension(Ngrids) :: W2Oname
      character (len=80), dimension(Ngrids) :: A2Oname
#endif
#ifdef SWAN_COUPLING
      character (len=80), dimension(Ngridss) :: W2Aname
      character (len=80), dimension(Ngridss) :: A2Wname
#endif
#ifdef WRF_MODEL
      character (len=80) :: Aname
#endif
#ifdef INWAVE_MODEL
      character (len=80), dimension(Ngridss) :: IWINIname
      character (len=80), dimension(Ngridss) :: IWBNDname
      character (len=80), dimension(Ngridss) :: IWINIname
      character (len=80), dimension(Ngridss) :: IWBNDname
      character (len=80), dimension(Ngridss) :: IWSWNname
#endif

      END MODULE mod_coupler_iounits
