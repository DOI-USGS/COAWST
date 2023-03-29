/*
** svn $Id: fennel_var.h 1054 2021-03-06 19:47:12Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Fennel et al. (2006) ecosystem   **
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
            CASE ('idTvar(iNH4_)')
              idTvar(iNH4_)=varid
#ifdef PO4
            CASE ('idTvar(iPO4_)')
              idTvar(iPO4_)=varid
#endif
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
            CASE ('idTvar(iODU_)')
              idTvar(iODU_)=varid
#ifdef CARBON
            CASE ('idTvar(iTIC_)')
              idTvar(iTIC_)=varid
            CASE ('idTvar(iTAlk)')
              idTvar(iTAlk)=varid
            CASE ('idTvar(iLDeC)')
              idTvar(iLDeC)=varid
            CASE ('idTvar(iSDeC)')
              idTvar(iSDeC)=varid
#endif
#ifdef OXYGEN
            CASE ('idTvar(iOxyg)')
              idTvar(iOxyg)=varid
#endif
#ifdef RIVER_DON
            CASE ('idTvar(iRDeN)')
              idTvar(iRDeN)=varid
# ifdef CARBON
            CASE ('idTvar(iRDeC)')
              idTvar(iRDeC)=varid
            CASE ('idTvar(iTIC_)')
             idTvar(iTIC_)=varid
            CASE ('idTvar(iTAlk)')
              idTvar(iTAlk)=varid
            CASE ('idTvar(iLDeC)')
              idTvar(iLDeC)=varid
            CASE ('idTvar(iSDeC)')
              idTvar(iSDeC)=varid
# endif
# ifdef OXYGEN
            CASE ('idTvar(iOxyg)')
              idTvar(iOxyg)=varid
# endif
# ifdef ODU
            CASE ('idTvar(iODU_)')
              idTvar(iODU_)=varid
# endif
#endif

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
# ifdef PO4
            CASE ('idTads(iPO4_)')
              idTads(iPO4_)=varid
# endif
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
# ifdef CARBON
            CASE ('idTads(iTIC_)')
              idTads(iTIC_)=varid
            CASE ('idTads(iTAlk)')
              idTads(iTAlk)=varid
            CASE ('idTads(iLDeC)')
              idTads(iLDeC)=varid
            CASE ('idTads(iSDeC)')
              idTads(iSDeC)=varid
# endif
# ifdef OXYGEN
            CASE ('idTads(iOxyg)')
              idTads(iOxyg)=varid
# endif
# ifdef ODU
              CASE ('idTads(iODU_)')
                idTads(iODU_)=varid
# endif
# ifdef RIVER_DON
            CASE ('idTads(iRDeN)')
              idTads(iRDeN)=varid
#  ifdef CARBON
            CASE ('idTads(iRDeC)')
              idTads(iRDeC)=varid
#  endif
# endif
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

#ifdef PO4
            CASE ('idTbry(iwest,iPO4_)')
              idTbry(iwest,iPO4_)=varid
            CASE ('idTbry(ieast,iPO4_)')
              idTbry(ieast,iPO4_)=varid
            CASE ('idTbry(isouth,iPO4_)')
              idTbry(isouth,iPO4_)=varid
            CASE ('idTbry(inorth,iPO4_)')
              idTbry(inorth,iPO4_)=varid
#endif
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

#ifdef CARBON
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
#endif
#ifdef OXYGEN
            CASE ('idTbry(iwest,iOxyg)')
              idTbry(iwest,iOxyg)=varid
            CASE ('idTbry(ieast,iOxyg)')
              idTbry(ieast,iOxyg)=varid
            CASE ('idTbry(isouth,iOxyg)')
              idTbry(isouth,iOxyg)=varid
            CASE ('idTbry(inorth,iOxyg)')
              idTbry(inorth,iOxyg)=varid
#endif
#ifdef RIVER_DON
            CASE ('idTbry(iwest,iRDeN)')
              idTbry(iwest,iRDeN)=varid
            CASE ('idTbry(ieast,iRDeN)')
              idTbry(ieast,iRDeN)=varid
            CASE ('idTbry(isouth,iRDeN)')
              idTbry(isouth,iRDeN)=varid
            CASE ('idTbry(inorth,iRDeN)')
              idTbry(inorth,iRDeN)=varid
# ifdef CARBON
            CASE ('idTbry(iwest,iRDeC)')
              idTbry(iwest,iRDeC)=varid
            CASE ('idTbry(ieast,iRDeC)')
              idTbry(ieast,iRDeC)=varid
            CASE ('idTbry(isouth,iRDeC)')
              idTbry(isouth,iRDeC)=varid
            CASE ('idTbry(inorth,iRDeC)')
              idTbry(inorth,iRDeC)=varid
# endif
#endif
#ifdef ODU
              CASE ('idTbry(iwest,iODU_)')
                idTbry(iwest,iODU_)=varid
              CASE ('idTbry(ieast,iODU_)')
                idTbry(ieast,iODU_)=varid
              CASE ('idTbry(isouth,iODU_)')
                idTbry(isouth,iODU_)=varid
              CASE ('idTbry(inorth,iODU_)')
                idTbry(inorth,iODU_)=varid
#endif


/*
**  Biological tracers point Source/Sinks (river runoff).
*/

            CASE ('idRtrc(iNO3_)')
              idRtrc(iNO3_)=varid
            CASE ('idRtrc(iNH4_)')
              idRtrc(iNH4_)=varid
#ifdef PO4
            CASE ('idRtrc(iPO4_)')
              idRtrc(iPO4_)=varid
#endif
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
#ifdef CARBON
            CASE ('idRtrc(iTIC_)')
              idRtrc(iTIC_)=varid
            CASE ('idRtrc(iTAlk)')
              idRtrc(iTAlk)=varid
            CASE ('idRtrc(iLDeC)')
              idRtrc(iLDeC)=varid
            CASE ('idRtrc(iSDeC)')
              idRtrc(iSDeC)=varid
#endif
#ifdef OXYGEN
            CASE ('idRtrc(iOxyg)')
              idRtrc(iOxyg)=varid
#endif
#ifdef RIVER_DON
            CASE ('idRtrc(iRDeN)')
              idRtrc(iRDeN)=varid
# ifdef CARBON
            CASE ('idRtrc(iRDeC)')
              idRtrc(iRDeC)=varid
# endif
#endif
#ifdef ODU
              CASE ('idRtrc(iODU_)')
                idRtrc(iODU_)=varid
#endif


#ifdef DIAGNOSTICS_BIO

/*
**  Biological tracers term diagnostics.
*/

# ifdef DENITRIFICATION
            CASE ('iDbio2(iDNIT)')
              iDbio2(iDNIT)=varid
# endif
# ifdef CARBON
            CASE ('iDbio2(iCOfx)')
              iDbio2(iCOfx)=varid
            CASE ('iDbio2(ipCO2)')
              iDbio2(ipCO2)=varid
# endif
# ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
# endif
# ifdef SEDBIO_COUP
              CASE ('iDbio2(isdO2)')
                iDbio2(isdO2)=varid
              CASE ('iDbio2(iseO2)')
                iDbio2(iseO2)=varid
              CASE ('iDbio2(isdNO)')
                iDbio2(isdNO)=varid
              CASE ('iDbio2(iseNO)')
                iDbio2(iseNO)=varid
              CASE ('iDbio2(isdNH)')
                iDbio2(isdNH)=varid
              CASE ('iDbio2(iseNH)')
                iDbio2(iseNH)=varid
              CASE ('iDbio2(isdOD)')
                iDbio2(isdOD)=varid
              CASE ('iDbio2(iseOD)')
                iDbio2(iseOD)=varid
# endif
            CASE ('iDbio3(iPPro)')
              iDbio3(iPPro)=varid
            CASE ('iDbio3(iNO3u)')
              iDbio3(iNO3u)=varid
            CASE ('iDbio3(iNifx)')
              iDbio3(iNifx)=varid
#endif
