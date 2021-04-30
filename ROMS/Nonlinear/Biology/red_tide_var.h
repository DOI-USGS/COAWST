/*
** svn $Id: red_tide_var.h 1054 2021-03-06 19:47:12Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for red tide (Stock et al., 2005;        **
**  He et al., 2008) biological model variables that are used         **
**  in input and output NetCDF files. The metadata information        **
**  is read from "varinfo.dat".                                       **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

            CASE ('idAsrf')
              idAsrf=varid
            CASE ('idCyst')
              idCyst=varid
            CASE ('idODIN')
              idODIN=varid
            CASE ('idTvar(iDino)')
              idTvar(iDino)=varid

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

            CASE ('idTads(iDino)')
              idTads(iDino)=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

            CASE ('idTbry(iwest,iDino)')
              idTbry(iwest,iDino)=varid
            CASE ('idTbry(ieast,iDino)')
              idTbry(ieast,iDino)=varid
            CASE ('idTbry(isouth,iDino)')
              idTbry(isouth,iDino)=varid
            CASE ('idTbry(inorth,iDino)')
              idTbry(inorth,iDino)=varid

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

            CASE ('idRtrc(iDino)')
              idRtrc(iDino)=varid
