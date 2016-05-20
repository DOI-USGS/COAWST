/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines Powell et al. (2006) ecosystem model input parameters in  **
**  output NetCDF file.  It is included in routine "def_info.F".      **
**                                                                    **
************************************************************************
*/

!
!  Define NPZD (Powell et al., 2006) biological model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation due to sea water'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='AttPhy'
      Vinfo( 2)='light attenuation due to phytoplankton'
      Vinfo( 3)='meter2 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton initial slope of the P-I curve'
      Vinfo( 3)='meter2 watt-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Vm_NO3'
      Vinfo( 2)='nitrate uptake rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyMRD'
      Vinfo( 2)='phytoplankton mortality rate to the detritus pool'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PhyMRN'
      Vinfo( 2)='phytoplankton mortality rate to the Nitrogen pool'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ivlev'
      Vinfo( 2)='Ivlev constant for zooplankton grazing'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton grazing rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooEED'
      Vinfo( 2)='zooplankton excretion efficiency to detritus pool'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooEEN'
      Vinfo( 2)='zooplankton excretion efficiency to Nitrogen pool'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooMRD'
      Vinfo( 2)='zooplankton mortality rate to the detritus pool'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ZooMRN'
      Vinfo( 2)='zootoplankton mortality rate to the Nitrogen pool'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='DetRR'
      Vinfo( 2)='detritus remineralization rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wPhy'
      Vinfo( 2)='phytoplankton sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wDet'
      Vinfo( 2)='detrital sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
