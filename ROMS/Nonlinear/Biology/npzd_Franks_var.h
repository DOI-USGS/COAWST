/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Franks et al. (1986) ecosystem   **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iSDet)')
                idTvar(iSDet)=varid

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iSDet)')
                idTads(iSDet)=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid

              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid

              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid

              CASE ('idTbry(iwest,iSDet)')
                idTbry(iwest,iSDet)=varid
              CASE ('idTbry(ieast,iSDet)')
                idTbry(ieast,iSDet)=varid
              CASE ('idTbry(isouth,iSDet)')
                idTbry(isouth,iSDet)=varid
              CASE ('idTbry(inorth,iSDet)')
                idTbry(inorth,iSDet)=varid


/*
**  Biological tracers point Source/Sinks (river runoff).
*/


              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iSDet)')
                idRtrc(iSDet)=varid
