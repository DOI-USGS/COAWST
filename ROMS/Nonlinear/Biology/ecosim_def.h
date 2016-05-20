/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines EcoSim bio-optical model input parameters in output       **
**  NetCDF files.  It is included in routine "def_info.F".            **
**                                                                    **
************************************************************************
*/

!
!  Define EcoSim bio-optical model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='maximum number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RtUVR_flag'
      Vinfo( 2)='switch to calculate CDOC UV photolysis.'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='NFIX_flag'
      Vinfo( 2)='switch to calculate temperature based N fixation'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Regen_flag'
      Vinfo( 2)='switch to calculate fecal matter regeneration'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsNO3'
      Vinfo( 2)='half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='micromole_NO3 liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsNH4'
      Vinfo( 2)='half-saturation for phytoplankton NH4 uptake'
      Vinfo( 3)='micromole_NH4 liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsSiO'
      Vinfo( 2)='half-saturation for phytoplankton SiO uptake'
      Vinfo( 3)='micromole_SiO liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsPO4'
      Vinfo( 2)='half-saturation for phytoplankton PO4 uptake'
      Vinfo( 3)='micromole_PO4 liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsFe'
      Vinfo( 2)='half-saturation for phytoplankton Fe uptake'
      Vinfo( 3)='micromole_Fe liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GtALG_max'
      Vinfo( 2)='maximum phytoplankton 24 hour growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyTbase'
      Vinfo( 2)=                                                        &
     &         'phytoplankton temperature base for exponential response'
      Vinfo( 3)='Celsius'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyTfac'
      Vinfo( 2)='phytoplankton exponential temperature factor'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BET_'
      Vinfo( 2)='nitrate uptake inhibition for NH4'
      Vinfo( 3)='micromole-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='maxC2nALG'
      Vinfo( 2)='maximum phytoplankton C:N ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='minC2nALG'
      Vinfo( 2)='balanced phytoplankton C:N ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2nALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:N ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='maxC2SiALG'
      Vinfo( 2)='maximum phytoplankton C:Si ratio'
      Vinfo( 3)='micromole_C micromole_Si-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='minC2SiALG'
      Vinfo( 2)='balanced phytoplankton C:Si ratio'
      Vinfo( 3)='micromole_C micromole_Si-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2SiALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:Si ratio'
      Vinfo( 3)='micromole_C micromole_Si-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='maxC2pALG'
      Vinfo( 2)='maximum phytoplankton C:P ratio'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='minC2pALG'
      Vinfo( 2)='balanced phytoplankton C:P ratio'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2pALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:P ratio'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='maxC2FeALG'
      Vinfo( 2)='maximum phytoplankton C:Fe ratio'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='minC2FeALG'
      Vinfo( 2)='balanced phytoplankton C:Fe ratio'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2FeALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:Fe ratio'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='qu_yld'
      Vinfo( 2)='maximum quantum yield'
      Vinfo( 3)='micromole_C micromole_quanta-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='E0_comp'
      Vinfo( 2)='compensation light level'
      Vinfo( 3)='micromole_quanta'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='E0_inhib'
      Vinfo( 2)='light level for onset of photoinhibition'
      Vinfo( 3)='micromole_quanta'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='inhib_fac'
      Vinfo( 2)='exponential decay factor for light limited growth'
      Vinfo( 3)='micromole_quanta-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2CHL_max'
      Vinfo( 2)='maximum lighted limited C:Chl ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxC2Cl'
      Vinfo( 2)='rate of change in light limited C:CHL ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1 micromole_quanta-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_C2Cl'
      Vinfo( 2)='minimum lighted limited C:Chl ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxC2Cn'
      Vinfo( 2)='rate of change in nutient limited C:Chl ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1 micromole_N micromole_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_C2Cn'
      Vinfo( 2)='minimum nutrient limited C:CHL ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxPacEff'
      Vinfo( 2)='rate of change in package effect'
      Vinfo( 3)='microgram_Chl microgram_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_PacEff'
      Vinfo( 2)='maximum package effect'
      Vinfo( 3)='microgram_Chl microgram_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxChlB'
      Vinfo( 2)='rate of change in the Chl_b:CHL_a ratio'
      Vinfo( 3)=                                                        &
     & 'microgram_Chl_b microgram_Chl_a-1 microgram_Chl_a microgrma_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_ChlB'
      Vinfo( 2)='maximum Chl_b:Chl_a ratio'
      Vinfo( 3)='microgram_Chl_b microgram_Chl_a-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxChlC'
      Vinfo( 2)='rate of change in the Chl_c:Chl_a ratio'
      Vinfo( 3)=                                                        &
     & 'microgram_Chl_c microgram_Chl_a-1 microgram_Chl_a microgrma_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_ChlC'
      Vinfo( 2)='maximum Chl_c:Chl_a ratio'
      Vinfo( 3)='microgram_Chl_c microgram_Chl_a-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxPSC'
      Vinfo( 2)='rate of change in the PSC:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &   'microgram_PSC microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_PSC'
      Vinfo( 2)='maximum PSC:Chl_a ratio'
      Vinfo( 3)='microgram_PSC microgram_Chl_a-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxPPC'
      Vinfo( 2)='rate of change in the PPC:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &   'microgram_PPC microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_PPC'
      Vinfo( 2)='maximum PPC:Chl_a ratio'
      Vinfo( 3)='microgram_Chl_c microgram_Chl_a-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxLPUb'
      Vinfo( 2)='rate of change in the LPUb:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &  'microgram_LPUb microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_LPUb'
      Vinfo( 2)='Maximum LPUb:Chl_a ratio'
      Vinfo( 3)='microgram_HPUb microgram_Chl_a-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='mxHPUb'
      Vinfo( 2)='rate of change in the HPUb:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &  'microgram_HPUb microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='b_HPUb'
      Vinfo( 2)='maximum HPUb:Chl_a ratio'
      Vinfo( 3)='microgram_HPUb microgram_Chl_a-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='FecDOC'
      Vinfo( 2)='proportion of grazing stress apportioned to DOM'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='FecPEL'
      Vinfo( 2)='proportion of grazing stress apportioned to fecal'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='FecCYC'
      Vinfo( 2)='proportion of grazing stress that is recycled'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ExALG'
      Vinfo( 2)='proportion of daily production lost to excretion'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='WS'
      Vinfo( 2)='phytoplankton sinking speed'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsGRZ'
      Vinfo( 2)='phytoplankton grazing parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MinRefuge'
      Vinfo( 2)='refuge phytoplankton population'
      Vinfo( 3)='micromole_C liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RefugeDep'
      Vinfo( 2)='maximum refuge phytoplankton depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Norm_Vol'
      Vinfo( 2)='normalized volume factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Norm_Surf'
      Vinfo( 2)='normalized surface area factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsDOP'
      Vinfo( 2)='half-saturation constant for DOP uptake'
      Vinfo( 3)='micromole_DOP liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2pALKPHOS'
      Vinfo( 2)='C:P ratio where DOP uptake begins'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsDON'
      Vinfo( 2)='Half Saturation Constant for DON uptake'
      Vinfo( 3)='micromole_DON liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2nNupDON'
      Vinfo( 2)='C:N ratio where DON uptake begins'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/phydim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsDOC_ba'
      Vinfo( 2)='half-saturation constant for bacteria DOC uptake'
      Vinfo( 3)='micromole_DOC liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/bacdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='GtBAC_max'
      Vinfo( 2)='maximum 24 hour bacterial growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/bacdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BacTbase'
      Vinfo( 2)='bacteria temperature base for exponential response'
      Vinfo( 3)='Celsius'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/bacdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BacTfac'
      Vinfo( 2)='bacteria exponential temperature factor'
      Vinfo( 3)='Celsius'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/bacdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2nBAC'
      Vinfo( 2)='carbon to nitrogen ratio of bacteria'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2pBAC'
      Vinfo( 2)='carbon to phosphorus ratio of bacteria'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='C2FeBAC'
      Vinfo( 2)='carbon to iron ratio of bacteria'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BacDOC'
      Vinfo( 2)=                                                        &
     &        'proportion of bacteria grazing stress apportioned to DOM'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BacPEL'
      Vinfo( 2)=                                                        &
     &      'proportion of bacteria grazing stress apportioned to fecal'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='BacCYC'
      Vinfo( 2)='proportion of bacteria grazing stress recycled'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ExBAC_c'
      Vinfo( 2)=                                                        &
     &      'bacterial recalcitrant C excretion as proportion of uptake'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ExBacC2N'
      Vinfo( 2)='bacterial recalcitrant excretion carbon:nitrogen ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Bac_Ceff'
      Vinfo( 2)='bacterial gross growth carbon efficiency'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RtNIT'
      Vinfo( 2)='maximum bacterial nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HsNIT'
      Vinfo( 2)='half-saturation constant for bacterial nitrification'
      Vinfo( 3)='micromole_NH4 liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='cDOCfrac_c'
      Vinfo( 2)=                                                        &
     &  'colored fraction of DOC from phytoplankton and bacteria losses'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/domdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RtUVR_DIC'
      Vinfo( 2)='UV degradation of DOC into DIC at 410 nanometer'
      Vinfo( 3)='micromole meter-1 liter-1 hour-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RtUVR_DOC'
      Vinfo( 2)=                                                        &
     &  'UV degradation of DOC to colorless labile DOC at 410 nanometer'
      Vinfo( 3)='micromole meter-1 liter-1 hour-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='WF'
      Vinfo( 2)='fecal sinking flux'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegTbase'
      Vinfo( 2)=                                                        &
     &    'fecal regeneration temperature base for exponential response'
      Vinfo( 3)='Celsius'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegTfac'
      Vinfo( 2)='fecal regeneration exponential temperature factor'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegCR'
      Vinfo( 2)='fecal carbon regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegNR'
      Vinfo( 2)='fecal nitrogen regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegSR'
      Vinfo( 2)='fecal silica regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegPR'
      Vinfo( 2)='fecal phosphorus regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='RegFR'
      Vinfo( 2)='fecal iron regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/fecdim/), Aval, Vinfo, ncname)
      IF (exit_flag.ne.NoError) RETURN
