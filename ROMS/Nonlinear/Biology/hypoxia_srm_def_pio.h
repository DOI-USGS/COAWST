/*
** git $Id$
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines Hypoxia Simple Respiration Model input parameters in      **
**  output NetCDF files. It is included in routine  "def_info.F".     **
**                                                                    **
************************************************************************
*/

!
!  Define Hypoxia Simple Respiration Model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='ResRate'
      Vinfo( 2)='total biological respiration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
