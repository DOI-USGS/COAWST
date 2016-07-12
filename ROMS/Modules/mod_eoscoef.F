#include "cppdefs.h"
      MODULE mod_eoscoef
#ifdef NONLIN_EOS
!
!svn $Id: mod_eoscoef.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Polynomial  expansion  coefficients for the computation of          !
!  "in situ" density  and other associated quantities via the          !
!  nonlinear equation of state for seawater  as a function of          !
!  potential temperature, salinity, and pressure (Jackett and          !
!  McDougall, 1992).                                                   !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        real(r8), parameter :: A00 = +1.909256e+04_r8
        real(r8), parameter :: A01 = +2.098925e+02_r8
        real(r8), parameter :: A02 = -3.041638e+00_r8
        real(r8), parameter :: A03 = -1.852732e-03_r8
        real(r8), parameter :: A04 = -1.361629e-05_r8
        real(r8), parameter :: B00 = +1.044077e+02_r8
        real(r8), parameter :: B01 = -6.500517e+00_r8
        real(r8), parameter :: B02 = +1.553190e-01_r8
        real(r8), parameter :: B03 = +2.326469e-04_r8
        real(r8), parameter :: D00 = -5.587545e+00_r8
        real(r8), parameter :: D01 = +7.390729e-01_r8
        real(r8), parameter :: D02 = -1.909078e-02_r8
        real(r8), parameter :: E00 = +4.721788e-01_r8
        real(r8), parameter :: E01 = +1.028859e-02_r8
        real(r8), parameter :: E02 = -2.512549e-04_r8
        real(r8), parameter :: E03 = -5.939910e-07_r8
        real(r8), parameter :: F00 = -1.571896e-02_r8
        real(r8), parameter :: F01 = -2.598241e-04_r8
        real(r8), parameter :: F02 = +7.267926e-06_r8
        real(r8), parameter :: G00 = +2.042967e-03_r8
        real(r8), parameter :: G01 = +1.045941e-05_r8
        real(r8), parameter :: G02 = -5.782165e-10_r8
        real(r8), parameter :: G03 = +1.296821e-07_r8
        real(r8), parameter :: H00 = -2.595994e-07_r8
        real(r8), parameter :: H01 = -1.248266e-09_r8
        real(r8), parameter :: H02 = -3.508914e-09_r8
        real(r8), parameter :: Q00 = +9.99842594e+02_r8
        real(r8), parameter :: Q01 = +6.793952e-02_r8
        real(r8), parameter :: Q02 = -9.095290e-03_r8
        real(r8), parameter :: Q03 = +1.001685e-04_r8
        real(r8), parameter :: Q04 = -1.120083e-06_r8
        real(r8), parameter :: Q05 = +6.536332e-09_r8
        real(r8), parameter :: U00 = +8.24493e-01_r8
        real(r8), parameter :: U01 = -4.08990e-03_r8
        real(r8), parameter :: U02 = +7.64380e-05_r8
        real(r8), parameter :: U03 = -8.24670e-07_r8
        real(r8), parameter :: U04 = +5.38750e-09_r8
        real(r8), parameter :: V00 = -5.72466e-03_r8
        real(r8), parameter :: V01 = +1.02270e-04_r8
        real(r8), parameter :: V02 = -1.65460e-06_r8
        real(r8), parameter :: W00 = +4.8314e-04_r8

#endif

      END MODULE mod_eoscoef
