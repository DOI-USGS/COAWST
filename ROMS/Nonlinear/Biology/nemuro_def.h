/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines Nemuro ecosystem model input parameters in output NetCDF  **
**  files. It is included in routine "def_info.F".                    **
**                                                                    **
************************************************************************
*/

!
!  Define Nemuro ecosystem model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation by seawater'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttPS'
      Vinfo( 2)='light attenuation due to small phytoplankton, '//      &
     &          'self-shading coefficient'
      Vinfo( 3)='meter2 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttPL'
      Vinfo( 2)='light attenuation due to large phytoplankton, '//      &
     &          'self-shading coefficient'
      Vinfo( 3)='meter2 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PARfrac'
      Vinfo( 2)='photosynthetically available radiation fraction'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AlphaPS'
      Vinfo( 2)='small phytoplankton initial slope of P-I curve'
      Vinfo( 3)='meter2 watt-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AlphaPL'
      Vinfo( 2)='Large phytoplankton initial slope of P-I curve'
      Vinfo( 3)='meter2 watt-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BetaPS'
      Vinfo( 2)='small phytoplankton photoinhibition coefficient'
      Vinfo( 3)='meter2 watt-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BetaPL'
      Vinfo( 2)='large phytoplankton photoinhibition coefficient'
      Vinfo( 3)='meter2 watt-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='VmaxS'
      Vinfo( 2)='small phytoplankton maximum photosynthetic rate '//    &
     &          'at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='VmaxL'
      Vinfo( 2)='large phytoplankton maximum photosynthetic rate '//    &
     &          'at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NO3S'
      Vinfo( 2)='small phytoplankton half-saturation constant '//       &
     &          'for nitrate'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NO3L'
      Vinfo( 2)='large phytoplankton half saturation constant '//       &
     &          'for nitrate'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NH4S'
      Vinfo( 2)=                                                        &
          'small phytoplankton half saturation constant for ammonium'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NH4L'
      Vinfo( 2)='large phytoplankton half saturation constant '//       &
     &          'for ammonium'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_SiL'
      Vinfo( 2)='large phytoplankton half saturation constant '//       &
     &          'for silicate'
      Vinfo( 3)='millimole_Si meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PusaiS'
      Vinfo( 2)='small phytoplankton ammonium inhibition coefficient'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PusaiL'
      Vinfo( 2)='small phytoplankton ammonium inhibition coefficient'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KGppS'
      Vinfo( 2)='small phytoplankton temperature coefficient for '//    &
     &          'photosynthetic rate'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KGppL'
      Vinfo( 2)='large phytoplankton temperature coefficient for '//    &
     &          'photosynthetic rate'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ResPS0'
      Vinfo( 2)='small phytoplankton respiration rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ResPL0'
      Vinfo( 2)='large phytoplankton respiration rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KResPS'
      Vinfo( 2)='small phytoplankton temperature coefficient '//        &
     &          'for respiration'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KResPL'
      Vinfo( 2)='large phytoplankton temperature coefficient '//        &
     &          'for respiration'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GammaS'
      Vinfo( 2)='small phytoplankton ratio of extracellular '//         &
     &          'excretion to photosynthesis'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GammaL'
      Vinfo( 2)='large phytoplankton ratio of extracellular '//         &
     &          'excretion to photosynthesis'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MorPS0'
      Vinfo( 2)='small phytoplankton mortality rate at 0 Celsius'
      Vinfo( 3)='meter3 millimole_N-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MorPL0'
      Vinfo( 2)='large phytoplankton mortality rate at 0 Celsius'
      Vinfo( 3)='meter3 millimole_N-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KMorPS'
      Vinfo( 2)='small phytoplankton temperature coefficient for '//    &
     &           'mortality'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KMorPL'
      Vinfo( 2)='large phytoplankton temperature coefficient for '//    &
     &           'mortality'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxSps'
      Vinfo( 2)='small zooplankton maximum grazing rate on '//          &
     &          'small phytoplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxSpl'
      Vinfo( 2)='small zooplankton maximum grazing rate on '//          &
     &          'large phytoplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxLps'
      Vinfo( 2)='large zooplankton maximum grazing rate on '//          &
     &          'small phytoplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxLpl'
      Vinfo( 2)='large zooplankton maximum grazing rate on '//          &
     &          'large phytoplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxLzs'
      Vinfo( 2)='large zooplankton maximum grazing rate on '//          &
     &          'small zooplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxPpl'
      Vinfo( 2)='predator zooplankton maximum grazing rate on '//       &
     &          'large phytoplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxPzs'
      Vinfo( 2)='predator zooplankton maximum grazing rate on '//       &
     &          'small zooplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GRmaxPzl'
      Vinfo( 2)='predator zooplankton maximum grazing rate on '//       &
     &          'large zooplankton at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KGraS'
      Vinfo( 2)='small zooplankton temperature coefficient '//          &
     &          'for grazing'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KGraL'
      Vinfo( 2)='large zooplankton temperature coefficient '//          &
     &          'for grazing'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KGraP'
      Vinfo( 2)='predator zooplankton temperature coefficient '//       &
     &          'for grazing'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LamS'
      Vinfo( 2)='small zooplankton Ivlev constant'
      Vinfo( 3)='meter3 millimole_N'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LamL'
      Vinfo( 2)='large zooplankton Ivlev constant'
      Vinfo( 3)='meter3 millimole_N'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LamP'
      Vinfo( 2)='predator zooplankton Ivlev constant'
      Vinfo( 3)='meter3 millimole_N'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef HOLLING_GRAZING
      Vinfo( 1)='K_PS2ZS'
      Vinfo( 2)='small zooplankton squared half-saturation '//          &
     &          'coefficient for ingestion on small phytoplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_PL2ZS'
      Vinfo( 2)='small zooplankton squared half-saturation '//          &
     &          'coefficient for ingestion on large phytoplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_PS2ZL'
      Vinfo( 2)='large zooplankton squared half-saturation '//          &
     &          'coefficient for ingestion on small phytoplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_PL2ZL'
      Vinfo( 2)='large zooplankton squared half-saturation '//          &
     &          'coefficient for ingestion on large phytoplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_ZS2ZL'
      Vinfo( 2)='large zooplankton squared half-saturation '//          &
     &          'coefficient for ingestion on small zooplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_PL2ZP'
      Vinfo( 2)='predator zooplankton squared half-saturation '//       &
     &          'coefficient for ingestion on large phytoplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_ZS2ZP'
      Vinfo( 2)='predator zooplankton squared half-saturation '//       &
     &          'coefficient for ingestion on small zooplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_ZL2ZP'
      Vinfo( 2)='predator zooplankton squared half-saturation '//       &
     &          'coefficient for ingestion on large zooplankton'
      Vinfo( 3)='millimole_N2 meter-6'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#else

      Vinfo( 1)='PS2ZSstar'
      Vinfo( 2)='small zooplankton threshold value for grazing '//      &
     &          'on small phytoplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PL2ZSstar'
      Vinfo( 2)='small zooplankton threshold value for grazing '//      &
     &          'on large phytoplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PS2ZLstar'
      Vinfo( 2)='large zooplankton threshold value for grazing '//      &
     &          'on small phytoplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PL2ZLstar'
      Vinfo( 2)='large zooplankton threshold value for grazing '//      &
     &          'on large phytoplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZS2ZLstar'
      Vinfo( 2)='large zooplankton threshold value for grazing '//      &
     &          'on small zooplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PL2ZPstar'
      Vinfo( 2)='predator zooplankton threshold value for grazing '//   &
     &          'on large phytoplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZS2ZPstar'
      Vinfo( 2)='predator zooplankton threshold value for grazing '//   &
     &          'on small zooplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZL2ZPstar'
      Vinfo( 2)='predator zooplankton threshold value for grazing '//   &
     &          'on large zooplankton'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

      Vinfo( 1)='PusaiPL'
      Vinfo( 2)='predator zooplankton grazing inhibition coefficient'// &
     &          'on large phytoplankton'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PusaiZS'
      Vinfo( 2)='predator zooplankton grazing inhibition coefficient'// &
     &          'on small zooplankton'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MorZS0'
      Vinfo( 2)='small zooplankton mortality rate at 0 Celsius'
      Vinfo( 3)='meter3 millimole_N-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MorZL0'
      Vinfo( 2)='large zooplankton mortality rate at 0 Celsius'
      Vinfo( 3)='meter3 millimole_N-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MorZP0'
      Vinfo( 2)='predator zooplankton mortality rate at 0 Celsius'
      Vinfo( 3)='meter3 millimole_N-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AlphaZS'
      Vinfo( 2)='small zooplankton assimilation efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AlphaZL'
      Vinfo( 2)='large zooplankton assimilation efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AlphaZP'
      Vinfo( 2)='predator zooplankton assimilation efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BetaZS'
      Vinfo( 2)='small zooplankton growth efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BetaZL'
      Vinfo( 2)='large zooplankton growth efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BetaZP'
      Vinfo( 2)='predator zooplankton growth efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Nit0'
      Vinfo( 2)='NH4 nitrification rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='VP2N0'
      Vinfo( 2)='PON decomposition to NH4 rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='VP2D0'
      Vinfo( 2)='PON decomposition to DON rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='VD2N0'
      Vinfo( 2)='DON decomposition to NH4 rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='VO2S0'
      Vinfo( 2)='opal decomposition to silicate rate at 0 Celsius'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KNit'
      Vinfo( 2)='temperature coefficient for NH4 nitrification'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KP2D'
      Vinfo( 2)='temperature coefficient for PON decomposition to DON'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KP2N'
      Vinfo( 2)='temperature coefficient for PON decomposition to NH4'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KD2N'
      Vinfo( 2)='temperature coefficient for DON decomposition to NH4'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KO2S'
      Vinfo( 2)='temperature coefficient for opal decomposition to '//  &
     &          'silicate'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RSiN'
      Vinfo( 2)='Si:N ratio'
      Vinfo( 3)='millimole_Si millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='setVPON'
      Vinfo( 2)='PON setting velocity'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='setVOpal'
      Vinfo( 2)='opal setting velocity'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef IRON_LIMIT
      Vinfo( 1)='T_Fe'
      Vinfo( 2)='Iron updake timescale'
      Vinfo( 3)='day'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='A_Fe'
      Vinfo( 2)='Empirical FE:C power'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='B_Fe'
      Vinfo( 2)='Empirical FE:C coefficient'
      Vinfo( 3)='meter-1 C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SK_FeC'
      Vinfo( 2)='Small phytoplankton Fe:C at F=0.5'
      Vinfo( 3)='muM-Fe/M-C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LK_FeC'
      Vinfo( 2)='Large phytoplankton Fe:C at F=0.5'
      Vinfo( 3)='muM-Fe/M-C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='FeRR'
      Vinfo( 2)='Fe remineralization rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
