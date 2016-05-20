/*
** svn $Id: estuarybgc_def.h 429 2009-12-20 17:30:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2010 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines Fennel et al. (2006) ecosystem model input parameters in  **
**  output NetCDF files. It is included in routine "def_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Define Fennel et al. (2006) ecosystem model parameters.
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

      Vinfo( 1)='AttChl'
      Vinfo( 2)='light attenuation by chlorophyll'
      Vinfo( 3)='meter-2 milligram_Chl-1'
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

      Vinfo( 1)='Vp0'
      Vinfo( 2)='Eppley temperature-limited growth parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='I_thNH4'
      Vinfo( 2)='radiation threshold for nitrification'
      Vinfo( 3)='watt meter-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='D_p5NH4'
      Vinfo( 2)='half-saturation radiation for nitrification'
      Vinfo( 3)='watt meter-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='NitriR'
      Vinfo( 2)='nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NH4'
      Vinfo( 2)='inverse half-saturation for phytoplankton NH4 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_Phy'
      Vinfo( 2)='zooplankton half-saturation constant for ingestion'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Chl2C_m'
      Vinfo( 2)='maximum chlorophyll to carbon ratio'
      Vinfo( 3)='milligram_chl milligram_carbon-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ChlMin'
      Vinfo( 2)='minimum chlorophyll threshold'
      Vinfo( 3)='milligram_chl meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyCN'
      Vinfo( 2)='phytoplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyIP'
      Vinfo( 2)='phytoplankton, NH4 inhibition parameter'
      Vinfo( 3)='millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton, initial slope of P-I curve'
      Vinfo( 3)='milligram_C milligram_Chl-1 watt-1 meter2 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyMin'
      Vinfo( 2)='minimum phytoplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyMR'
      Vinfo( 2)='phytoplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooAE_N'
      Vinfo( 2)='zooplankton mitrogen assimilation efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooBM'
      Vinfo( 2)='rate for zooplankton basal metabolism'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooCN'
      Vinfo( 2)='zooplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooER'
      Vinfo( 2)='zooplankton specific excretion rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton maximum growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooMin'
      Vinfo( 2)='minimum zooplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooMR'
      Vinfo( 2)='zooplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LDeRRN'
      Vinfo( 2)='rate of large detritus nitrogen re-mineralization'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LDeRRC'
      Vinfo( 2)='rate of large detritus carbon re-mineralization'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='CoagR'
      Vinfo( 2)='coagulation rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SDeRRN'
      Vinfo( 2)='remineralization rate for small detritus nitrogen'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SDeRRC'
      Vinfo( 2)='remineralization rate for small detritus carbon'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wPhy'
      Vinfo( 2)='vertical sinking velocity for phytoplankton'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wLDet'
      Vinfo( 2)='vertical sinking velocity for large detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wSDet'
      Vinfo( 2)='vertical sinking velocity for small detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='pCO2air'
      Vinfo( 2)='partial pressure of CO2 in the air'
      Vinfo( 3)='parts per million by volume'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SPECTRAL_LIGHT
      Vinfo( 1)='SIGATRB'
      Vinfo( 2)='NAP absorb. cross section at 440 nm'
      Vinfo( 3)='meter2 gram-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='STRB'
      Vinfo( 2)='Spectral slope of NAP absorption'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BLTRB'
      Vinfo( 2)='Baseline NAP absorption'
      Vinfo( 3)='meter2 gram-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SIGBTRB'
      Vinfo( 2)='Scattering cross section of turbidity'
      Vinfo( 3)='meter2 gram-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ETASPEC'
      Vinfo( 2)='Scattering spectral exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BB2B'
      Vinfo( 2)='Particulate backscattering ratio'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SAV_BIOMASS
      Vinfo( 1)='GMODopt'
      Vinfo( 2)='SAV growth rate options'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KNSED'
      Vinfo( 2)='Half saturation coeff. for sediment N uptake'
      Vinfo( 3)='mu-M'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KNWC'
      Vinfo( 2)='Half saturation coeff. for water-column N uptake'
      Vinfo( 3)='mu-M'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='TOPT'
      Vinfo( 2)='Optimum temperature for SAV growth'
      Vinfo( 3)='deg C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='THTA'
      Vinfo( 2)='Temperature growth theta for growth model#1'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='THTA2'
      Vinfo( 2)='Temperature growth theta for growth model#2'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SCL'
      Vinfo( 2)='SAV growth fraction for growth model#1'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SCL2'
      Vinfo( 2)='SAV growth fraction for growth model#2'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KI'
      Vinfo( 2)='Half saturation for light'
      Vinfo( 3)='mu-E m-2 s-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SR'
      Vinfo( 2)='Surface reflectance of light'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LMBAMX'
      Vinfo( 2)='Maximum AGB in self-shading formulation'
      Vinfo( 3)='g C m-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KMAG'
      Vinfo( 2)='Above biomass ground mortality rate'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ARSC'
      Vinfo( 2)='Maximum fraction of PP respired'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ARC'
      Vinfo( 2)='Active respiration coefficient'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BSRC'
      Vinfo( 2)='Maximum fraction of biomass required'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)

      Vinfo( 1)='RC'
      Vinfo( 2)='Basal respiration coefficient'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)

      Vinfo( 1)='RtStTL'
      Vinfo( 2)='Seasonal root storage coefficient'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='DOWNt'
      Vinfo( 2)='Downward translocation coefficient'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)

      Vinfo( 1)='TRNS'
      Vinfo( 2)='Upward translocation coefficient'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)

      Vinfo( 1)='TCRIT'
      Vinfo( 2)='Critical temp. for development of abovegrd. biomass'
      Vinfo( 3)='deg C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)

      Vinfo( 1)='KM'
      Vinfo( 2)='Below ground biomass mortality'
      Vinfo( 3)='unitless'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
#endif
      
