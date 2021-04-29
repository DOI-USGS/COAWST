/*
** svn $Id: sediment_def.h 1054 2021-03-06 19:47:12Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines sediment model input parameters in output NetCDF files.   **
**  It is included in routine "def_info.F".                           **
**                                                                    **
************************************************************************
*/

!
!  Define sediment model parameters.
!
      Vinfo( 1)='minlayer_thick'
      Vinfo( 2)='depositional bed layer minimum thickness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='newlayer_thick'
      Vinfo( 2)='depositional bed thickness criteria to crate a '//     &
     &          'new layer'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef BEDLOAD
      Vinfo( 1)='bedload_coeff'
      Vinfo( 2)='bedload transport rate coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
!
!# ifdef BEDLOAD_VANDERA
      Vinfo( 1)='sg_zwbl'
      Vinfo( 2)='input elevation to get near-bottom current vel.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      Vinfo( 1)='sedslope_crit_wet'
      Vinfo( 2)='critical wet bed slope for slumping.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      Vinfo( 1)='sedslope_crit_dry'
      Vinfo( 2)='critical dry bed slope for slumping.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      Vinfo( 1)='slopefac_wet'
      Vinfo( 2)='bedload wet bed slumping factor.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      Vinfo( 1)='slopefac_dry'
      Vinfo( 2)='bedload dry bed slumping factor.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      Vinfo( 1)='bedload_vandera_alphac'
      Vinfo( 2)='bedload scale factor for currents contribution.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      Vinfo( 1)='bedload_vandera_alphaw'
      Vinfo( 2)='bedload scale factor for waves contribution.'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!# endif
!
#endif
!
!#ifdef ANA_SEDIMENT
      Vinfo( 1)='Sd50'
      Vinfo( 2)='median sediment grain diameter used in '//             &
     &          'uniform initial conditions'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='Srho'
      Vinfo( 2)='sediment grain density used in '//                     &
     &          'uniform initial conditions'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='Csed'
      Vinfo( 2)='sediment concentration used in '//                     &
     &          'uniform initial conditions'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!#endif

      Vinfo( 1)='Wsed'
      Vinfo( 2)='sediment particle settling velocity'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='Erate'
      Vinfo( 2)='sediment surface layer erosing rate'
      Vinfo( 3)='kilogram meter-2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='tau_ce'
      Vinfo( 2)='sediment critical shear for erosion'
      Vinfo( 3)='meter-2 second-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='tau_cd'
      Vinfo( 2)='sediment critical shear for deposition'
      Vinfo( 3)='meter-2 second-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='poros'
      Vinfo( 2)='sediment porosity'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/seddim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
