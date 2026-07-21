/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Hypoxia Simple Respiration       **
**  Model biological variables that are used in input and output      **
**  NetCDF files. The metadata information is read from               **
**  "varinfo.yaml".                                                   **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

          CASE ('idResR')
            idResR=varid
          CASE ('idTvar(iOxyg)')
            idTvar(iOxyg)=varid

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

          CASE ('idTads(iOxyg)')
            idTads(iOxyg)=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

          CASE ('idTbry(iwest,iOxyg)')
            idTbry(iwest,iOxyg)=varid
          CASE ('idTbry(ieast,iOxyg)')
            idTbry(ieast,iOxyg)=varid
          CASE ('idTbry(isouth,iOxyg)')
            idTbry(isouth,iOxyg)=varid
          CASE ('idTbry(inorth,iOxyg)')
            idTbry(inorth,iOxyg)=varid

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

          CASE ('idRtrc(iOxyg)')
            idRtrc(iOxyg)=varid
