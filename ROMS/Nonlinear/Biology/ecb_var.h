/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2026 The ROMS Group                             **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Estuarine Carbon Biogeochemistry **
**  (Feng et al., 2015, and follow-ups) model variables that are used **
**  in input and output NetCDF files.  The metadata information is    **
**  read from "varinfo.yaml".                                         **
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
          CASE ('idTvar(iNH4_)')
            idTvar(iNH4_)=varid
          CASE ('idTvar(iPO4_)')
            idTvar(iPO4_)=varid
          CASE ('idTvar(iPhyt)')
            idTvar(iPhyt)=varid
          CASE ('idTvar(iZoop)')
            idTvar(iZoop)=varid
          CASE ('idTvar(iLDeN)')
            idTvar(iLDeN)=varid
          CASE ('idTvar(iSDeN)')
            idTvar(iSDeN)=varid
          CASE ('idTvar(iChlo)')
            idTvar(iChlo)=varid
          CASE ('idTvar(iTIC_)')
            idTvar(iTIC_)=varid
          CASE ('idTvar(iTAlk)')
            idTvar(iTAlk)=varid
          CASE ('idTvar(iLDeC)')
            idTvar(iLDeC)=varid
          CASE ('idTvar(iSDeC)')
            idTvar(iSDeC)=varid
          CASE ('idTvar(iOxyg)')
            idTvar(iOxyg)=varid
          CASE ('idTvar(iDON_)')
            idTvar(iDON_)=varid
          CASE ('idTvar(iDOC_)')
            idTvar(iDOC_)=varid
          CASE ('idTvar(irDON)')
            idTvar(irDON)=varid
          CASE ('idTvar(irDOC)')
            idTvar(irDOC)=varid
          CASE ('idTvar(iCalc)')
            idTvar(iCalc)=varid

/*
**  Adjoint sensitivity state biological tracers.
*/

#if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR      || \
    defined SO_SEMI
          CASE ('idTads(iNO3_)')
            idTads(iNO3_)=varid
          CASE ('idTads(iNH4_)')
            idTads(iNH4_)=varid
          CASE ('idTads(iPO4_)')
            idTads(iPO4_)=varid
          CASE ('idTads(iPhyt)')
            idTads(iPhyt)=varid
          CASE ('idTads(iZoop)')
            idTads(iZoop)=varid
          CASE ('idTads(iLDeN)')
            idTads(iLDeN)=varid
          CASE ('idTads(iSDeN)')
            idTads(iSDeN)=varid
          CASE ('idTads(iChlo)')
            idTads(iChlo)=varid
          CASE ('idTads(iTIC_)')
            idTads(iTIC_)=varid
          CASE ('idTads(iTAlk)')
            idTads(iTAlk)=varid
          CASE ('idTads(iLDeC)')
            idTads(iLDeC)=varid
          CASE ('idTads(iSDeC)')
            idTads(iSDeC)=varid
          CASE ('idTads(iOxyg)')
            idTads(iOxyg)=varid
          CASE ('idTads(iDON_)')
            idTads(iDON_)=varid
          CASE ('idTads(iDOC_)')
            idTads(iDOC_)=varid
          CASE ('idTads(irDON)')
            idTads(irDON)=varid
          CASE ('idTads(irDOC)')
            idTads(irDOC)=varid
          CASE ('idTads(iCalc)')
            idTads(iCalc)=varid
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

          CASE ('idTbry(iwest,iNH4_)')
            idTbry(iwest,iNH4_)=varid
          CASE ('idTbry(ieast,iNH4_)')
            idTbry(ieast,iNH4_)=varid
          CASE ('idTbry(isouth,iNH4_)')
            idTbry(isouth,iNH4_)=varid
          CASE ('idTbry(inorth,iNH4_)')
            idTbry(inorth,iNH4_)=varid

          CASE ('idTbry(iwest,iPO4_)')
            idTbry(iwest,iPO4_)=varid
          CASE ('idTbry(ieast,iPO4_)')
            idTbry(ieast,iPO4_)=varid
          CASE ('idTbry(isouth,iPO4_)')
            idTbry(isouth,iPO4_)=varid
          CASE ('idTbry(inorth,iPO4_)')
            idTbry(inorth,iPO4_)=varid
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

          CASE ('idTbry(iwest,iSDeN)')
            idTbry(iwest,iSDeN)=varid
          CASE ('idTbry(ieast,iSDeN)')
            idTbry(ieast,iSDeN)=varid
          CASE ('idTbry(isouth,iSDeN)')
            idTbry(isouth,iSDeN)=varid
          CASE ('idTbry(inorth,iSDeN)')
            idTbry(inorth,iSDeN)=varid

          CASE ('idTbry(iwest,iLDeN)')
            idTbry(iwest,iLDeN)=varid
          CASE ('idTbry(ieast,iLDeN)')
            idTbry(ieast,iLDeN)=varid
          CASE ('idTbry(isouth,iLDeN)')
            idTbry(isouth,iLDeN)=varid
          CASE ('idTbry(inorth,iLDeN)')
            idTbry(inorth,iLDeN)=varid

          CASE ('idTbry(iwest,iChlo)')
            idTbry(iwest,iChlo)=varid
          CASE ('idTbry(ieast,iChlo)')
            idTbry(ieast,iChlo)=varid
          CASE ('idTbry(isouth,iChlo)')
            idTbry(isouth,iChlo)=varid
          CASE ('idTbry(inorth,iChlo)')
            idTbry(inorth,iChlo)=varid

          CASE ('idTbry(iwest,iSDeC)')
            idTbry(iwest,iSDeC)=varid
          CASE ('idTbry(ieast,iSDeC)')
            idTbry(ieast,iSDeC)=varid
          CASE ('idTbry(isouth,iSDeC)')
            idTbry(isouth,iSDeC)=varid
          CASE ('idTbry(inorth,iSDeC)')
            idTbry(inorth,iSDeC)=varid

          CASE ('idTbry(iwest,iLDeC)')
            idTbry(iwest,iLDeC)=varid
          CASE ('idTbry(ieast,iLDeC)')
            idTbry(ieast,iLDeC)=varid
          CASE ('idTbry(isouth,iLDeC)')
            idTbry(isouth,iLDeC)=varid
          CASE ('idTbry(inorth,iLDeC)')
            idTbry(inorth,iLDeC)=varid

          CASE ('idTbry(iwest,iTIC_)')
            idTbry(iwest,iTIC_)=varid
          CASE ('idTbry(ieast,iTIC_)')
            idTbry(ieast,iTIC_)=varid
          CASE ('idTbry(isouth,iTIC_)')
            idTbry(isouth,iTIC_)=varid
          CASE ('idTbry(inorth,iTIC_)')
            idTbry(inorth,iTIC_)=varid

          CASE ('idTbry(iwest,iTAlk)')
            idTbry(iwest,iTAlk)=varid
          CASE ('idTbry(ieast,iTAlk)')
            idTbry(ieast,iTAlk)=varid
          CASE ('idTbry(isouth,iTAlk)')
            idTbry(isouth,iTAlk)=varid
          CASE ('idTbry(inorth,iTAlk)')
            idTbry(inorth,iTAlk)=varid
          CASE ('idTbry(iwest,iOxyg)')
            idTbry(iwest,iOxyg)=varid
          CASE ('idTbry(ieast,iOxyg)')
            idTbry(ieast,iOxyg)=varid
          CASE ('idTbry(isouth,iOxyg)')
            idTbry(isouth,iOxyg)=varid
          CASE ('idTbry(inorth,iOxyg)')
            idTbry(inorth,iOxyg)=varid
          CASE ('idTbry(iwest,iDON_)')
            idTbry(iwest,iDON_)=varid
          CASE ('idTbry(ieast,iDON_)')
            idTbry(ieast,iDON_)=varid
          CASE ('idTbry(isouth,iDON_)')
            idTbry(isouth,iDON_)=varid
          CASE ('idTbry(inorth,iDON_)')
            idTbry(inorth,iDON_)=varid
          CASE ('idTbry(iwest,iDOC_)')
            idTbry(iwest,iDOC_)=varid
          CASE ('idTbry(ieast,iDOC_)')
            idTbry(ieast,iDOC_)=varid
          CASE ('idTbry(isouth,iDOC_)')
            idTbry(isouth,iDOC_)=varid
          CASE ('idTbry(inorth,iDOC_)')
            idTbry(inorth,iDOC_)=varid
          CASE ('idTbry(iwest,irDON)')
            idTbry(iwest,irDON)=varid
          CASE ('idTbry(ieast,irDON)')
            idTbry(ieast,irDON)=varid
          CASE ('idTbry(isouth,irDON)')
            idTbry(isouth,irDON)=varid
          CASE ('idTbry(inorth,irDON)')
            idTbry(inorth,irDON)=varid
          CASE ('idTbry(iwest,irDOC)')
            idTbry(iwest,irDOC)=varid
          CASE ('idTbry(ieast,irDOC)')
            idTbry(ieast,irDOC)=varid
          CASE ('idTbry(isouth,irDOC)')
            idTbry(isouth,irDOC)=varid
          CASE ('idTbry(inorth,irDOC)')
            idTbry(inorth,irDOC)=varid
          CASE ('idTbry(iwest,iCalc)')
            idTbry(iwest,iCalc)=varid
          CASE ('idTbry(ieast,iCalc)')
            idTbry(ieast,iCalc)=varid
          CASE ('idTbry(isouth,iCalc)')
            idTbry(isouth,iCalc)=varid
          CASE ('idTbry(inorth,iCalc)')
            idTbry(inorth,iCalc)=varid

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

          CASE ('idRtrc(iNO3_)')
            idRtrc(iNO3_)=varid
          CASE ('idRtrc(iNH4_)')
            idRtrc(iNH4_)=varid
          CASE ('idRtrc(iPO4_)')
            idRtrc(iPO4_)=varid
          CASE ('idRtrc(iPhyt)')
            idRtrc(iPhyt)=varid
          CASE ('idRtrc(iZoop)')
            idRtrc(iZoop)=varid
          CASE ('idRtrc(iLDeN)')
            idRtrc(iLDeN)=varid
          CASE ('idRtrc(iSDeN)')
            idRtrc(iSDeN)=varid
          CASE ('idRtrc(iChlo)')
            idRtrc(iChlo)=varid
          CASE ('idRtrc(iTIC_)')
            idRtrc(iTIC_)=varid
          CASE ('idRtrc(iTAlk)')
            idRtrc(iTAlk)=varid
          CASE ('idRtrc(iLDeC)')
            idRtrc(iLDeC)=varid
          CASE ('idRtrc(iSDeC)')
            idRtrc(iSDeC)=varid
          CASE ('idRtrc(iOxyg)')
            idRtrc(iOxyg)=varid
          CASE ('idRtrc(iDON_)')
            idRtrc(iDON_)=varid
          CASE ('idRtrc(iDOC_)')
            idRtrc(iDOC_)=varid
          CASE ('idRtrc(irDON)')
            idRtrc(irDON)=varid
          CASE ('idRtrc(irDOC)')
            idRtrc(irDOC)=varid
          CASE ('idRtrc(iCalc)')
            idRtrc(iCalc)=varid

#ifdef DIAGNOSTICS_BIO

/*
**  Biological tracers term diagnostics.
*/

          CASE ('iDbio2(iDNIT)')
            iDbio2(iDNIT)=varid
          CASE ('iDbio2(iCOfx)')
            iDbio2(iCOfx)=varid
          CASE ('iDbio2(ipCO2)')
            iDbio2(ipCO2)=varid
          CASE ('iDbio2(iO2fx)')
            iDbio2(iO2fx)=varid
          CASE ('iDbio3(iPPro)')
            iDbio3(iPPro)=varid
          CASE ('iDbio3(iNO3u)')
            iDbio3(iNO3u)=varid
          CASE ('iDbio2(iNbur)')
            iDbio2(iNbur)=varid
          CASE ('iDbio2(iCbot)')
            iDbio2(iCbot)=varid
          CASE ('iDbio2(iCbur)')
            iDbio2(iCbur)=varid
          CASE ('iDbio3(iPPCe)')
            iDbio3(iPPCe)=varid
          CASE ('iDbio3(iwdeN)')
            iDbio3(iwdeN)=varid
          CASE ('iDbio2(iSoxy)')
            iDbio2(iSoxy)=varid
          CASE ('iDbio3(iNifx)')
            iDbio3(iNifx)=varid
#endif
