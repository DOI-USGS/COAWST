/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines oyster floats model input parameters in output NetCDF     **
**  files.  It is included in routine "def_info.F".                   **
**                                                                    **
************************************************************************
*/

!
!  Define oyster float model (Dekshenieks et al. 1993, 1996, 1997;
!  Narvaez et al. 2012a,b).
!
      Vinfo( 1)='Larvae_GR0'
      Vinfo( 2)='oyster larvae growth rate'
      Vinfo( 3)='micrometer day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Larvae_size0'
      Vinfo( 2)='oyster larvae size in terms of length'
      Vinfo( 3)='micrometer'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='food_supply'
      Vinfo( 2)='constant oyster larvae food supply'
      Vinfo( 3)='gram Carbon liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='settle_size'
      Vinfo( 2)='oyster larvae settlement size'
      Vinfo( 3)='micrometer'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Turbidity effect parameters on planktonic oyster larvae growth.
!
      Vinfo( 1)='turb_ambi'
      Vinfo( 2)='ambient turbidity level effect on larvae growth'
      Vinfo( 3)='gram liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_axis'
      Vinfo( 2)='turbidity linear axis crossing in growth curve'
      Vinfo( 3)='gram liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_base'
      Vinfo( 2)='turbidity exponential base factor'
      Vinfo( 3)='gram liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_crit'
      Vinfo( 2)='critical turbidity value'
      Vinfo( 3)='gram liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_mean'
      Vinfo( 2)='turbidity exponential mean'
      Vinfo( 3)='gram liter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_rate'
      Vinfo( 2)='turbidity exponential rate'
      Vinfo( 3)='liter gram-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_size'
      Vinfo( 2)='minimum oyster larvae size affected by turbidity'
      Vinfo( 3)='micrometer'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='turb_slop'
      Vinfo( 2)='turbidity linear slope on oyster larvae growth curve'
      Vinfo( 3)='liter gram-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Planktonic oyster larvae vertical migration (swimming) parameters.
!
      Vinfo( 1)='slope_Sdec'
      Vinfo( 2)='oyster larvae swimming coefficient due to '//          &
     &          'decreasing salinity'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='slope_Sinc'
      Vinfo( 2)='oyster larvae swimming coefficient due to '//          &
     &          'increasing salinity'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='swim_Sdec'
      Vinfo( 2)='oyster larvae active swimming fraction due to '//      &
     &          'decreasing salinity'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='swim_Sinc'
      Vinfo( 2)='oyster larvae active swimming fraction due to '//      &
     &          'increasing salinity'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='swim_Tmin'
      Vinfo( 2)='oyster larvae minimum swimming time fraction'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='swim_Tmax'
      Vinfo( 2)='oyster larvae maximum swimming time fraction'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Planktonic larvae sinking parameters.
!
      Vinfo( 1)='sink_base'
      Vinfo( 2)='oyster larvae exponential base factor in sinking curve'
      Vinfo( 3)='millimeter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='sink_rate'
      Vinfo( 2)='oyster larvae exponential rate factor in sinking curve'
      Vinfo( 3)='micrometer-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='sink_size'
      Vinfo( 2)='oyster larvae exponential mean size in sinking curve'
      Vinfo( 3)='micrometer'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
