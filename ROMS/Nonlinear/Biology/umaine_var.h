/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices to the Nemuro ecosystem model            **
**  variables that are used in input and output NetCDF files.         **
**  The metadata nformation is read from "varinfo.dat".               **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

!
!  Model state biological tracers.
!

              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iSiOH)')
                idTvar(iSiOH)=varid
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
              CASE ('idTvar(iS1_N)')
                idTvar(iS1_N)=varid
              CASE ('idTvar(iS1_C)')
                idTvar(iS1_C)=varid
              CASE ('idTvar(iS1CH)')
                idTvar(iS1CH)=varid
              CASE ('idTvar(iS2_N)')
                idTvar(iS2_N)=varid
              CASE ('idTvar(iS2_C)')
                idTvar(iS2_C)=varid
              CASE ('idTvar(iS2CH)')
                idTvar(iS2CH)=varid
              CASE ('idTvar(iS3_N)')
                idTvar(iS3_N)=varid
              CASE ('idTvar(iS3_C)')
                idTvar(iS3_C)=varid
              CASE ('idTvar(iS3CH)')
                idTvar(iS3CH)=varid
              CASE ('idTvar(iZ1_N)')
                idTvar(iZ1_N)=varid
              CASE ('idTvar(iZ1_C)')
                idTvar(iZ1_C)=varid
              CASE ('idTvar(iZ2_N)')
                idTvar(iZ2_N)=varid
              CASE ('idTvar(iZ2_C)')
                idTvar(iZ2_C)=varid
              CASE ('idTvar(iBAC_)')
                idTvar(iBAC_)=varid
              CASE ('idTvar(iDD_N)')
                idTvar(iDD_N)=varid
              CASE ('idTvar(iDD_C)')
                idTvar(iDD_C)=varid
              CASE ('idTvar(iDDSi)')
                idTvar(iDDSi)=varid
              CASE ('idTvar(iLDON)')
                idTvar(iLDON)=varid
              CASE ('idTvar(iLDOC)')
                idTvar(iLDOC)=varid
              CASE ('idTvar(iSDON)')
                idTvar(iSDON)=varid
              CASE ('idTvar(iSDOC)')
                idTvar(iSDOC)=varid
              CASE ('idTvar(iCLDC)')
                idTvar(iCLDC)=varid
              CASE ('idTvar(iCSDC)')
                idTvar(iCSDC)=varid
              CASE ('idTvar(iDDCA)')
                idTvar(iDDCA)=varid

#ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
#endif
#ifdef IRON_LIMIT
              CASE ('idTvar(iS1_Fe)')
                idTvar(iS1_Fe)=varid
              CASE ('idTvar(iS2_Fe)')
                idTvar(iS2_Fe)=varid
              CASE ('idTvar(iS3_Fe)')
                idTvar(iS3_Fe)=varid
              CASE ('idTvar(iFeD_)')
                idTvar(iFeD_)=varid
#endif
#ifdef PRIMARY_PROD
              CASE ('idNPP')
                idNPP=varid
#endif

!
!  Biological tracers open boundary conditions.
!

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
              CASE ('idTbry(iwest,iSiOH)')
                idTbry(iwest,iSiOH)=varid
              CASE ('idTbry(ieast,iSiOH)')
                idTbry(ieast,iSiOH)=varid
              CASE ('idTbry(isouth,iSiOH)')
                idTbry(isouth,iSiOH)=varid
              CASE ('idTbry(inorth,iSiOH)')
                idTbry(inorth,iSiOH)=varid
              CASE ('idTbry(iwest,iPO4_)')
                idTbry(iwest,iPO4_)=varid
              CASE ('idTbry(ieast,iPO4_)')
                idTbry(ieast,iPO4_)=varid
              CASE ('idTbry(isouth,iPO4_)')
                idTbry(isouth,iPO4_)=varid
              CASE ('idTbry(inorth,iPO4_)')
                idTbry(inorth,iPO4_)=varid
              CASE ('idTbry(iwest,iS1_N)')
                idTbry(iwest,iS1_N)=varid
              CASE ('idTbry(ieast,iS1_N)')
                idTbry(ieast,iS1_N)=varid
              CASE ('idTbry(isouth,iS1_N)')
                idTbry(isouth,iS1_N)=varid
              CASE ('idTbry(inorth,iS1_N)')
                idTbry(inorth,iS1_N)=varid
              CASE ('idTbry(iwest,iS1_C)')
                idTbry(iwest,iS1_C)=varid
              CASE ('idTbry(ieast,iS1_C)')
                idTbry(ieast,iS1_C)=varid
              CASE ('idTbry(isouth,iS1_C)')
                idTbry(isouth,iS1_C)=varid
              CASE ('idTbry(inorth,iS1_C)')
                idTbry(inorth,iS1_C)=varid
              CASE ('idTbry(iwest,iS1CH)')
                idTbry(iwest,iS1CH)=varid
              CASE ('idTbry(ieast,iS1CH)')
                idTbry(ieast,iS1CH)=varid
              CASE ('idTbry(isouth,iS1CH)')
                idTbry(isouth,iS1CH)=varid
              CASE ('idTbry(inorth,iS1CH)')
                idTbry(inorth,iS1CH)=varid
              CASE ('idTbry(iwest,iS2_N)')
                idTbry(iwest,iS2_N)=varid
              CASE ('idTbry(ieast,iS2_N)')
                idTbry(ieast,iS2_N)=varid
              CASE ('idTbry(isouth,iS2_N)')
                idTbry(isouth,iS2_N)=varid
              CASE ('idTbry(inorth,iS2_N)')
                idTbry(inorth,iS2_N)=varid
              CASE ('idTbry(iwest,iS2_C)')
                idTbry(iwest,iS2_C)=varid
              CASE ('idTbry(ieast,iS2_C)')
                idTbry(ieast,iS2_C)=varid
              CASE ('idTbry(isouth,iS2_C)')
                idTbry(isouth,iS2_C)=varid
              CASE ('idTbry(inorth,iS2_C)')
                idTbry(inorth,iS2_C)=varid
              CASE ('idTbry(iwest,iS2CH)')
                idTbry(iwest,iS2CH)=varid
              CASE ('idTbry(ieast,iS2CH)')
                idTbry(ieast,iS2CH)=varid
              CASE ('idTbry(isouth,iS2CH)')
                idTbry(isouth,iS2CH)=varid
              CASE ('idTbry(inorth,iS2CH)')
                idTbry(inorth,iS2CH)=varid
              CASE ('idTbry(iwest,iS3_N)')
                idTbry(iwest,iS3_N)=varid
              CASE ('idTbry(ieast,iS3_N)')
                idTbry(ieast,iS3_N)=varid
              CASE ('idTbry(isouth,iS3_N)')
                idTbry(isouth,iS3_N)=varid
              CASE ('idTbry(inorth,iS3_N)')
                idTbry(inorth,iS3_N)=varid
              CASE ('idTbry(iwest,iS3_C)')
                idTbry(iwest,iS3_C)=varid
              CASE ('idTbry(ieast,iS3_C)')
                idTbry(ieast,iS3_C)=varid
              CASE ('idTbry(isouth,iS3_C)')
                idTbry(isouth,iS3_C)=varid
              CASE ('idTbry(inorth,iS3_C)')
                idTbry(inorth,iS3_C)=varid
              CASE ('idTbry(iwest,iS3CH)')
                idTbry(iwest,iS3CH)=varid
              CASE ('idTbry(ieast,iS3CH)')
                idTbry(ieast,iS3CH)=varid
              CASE ('idTbry(isouth,iS3CH)')
                idTbry(isouth,iS3CH)=varid
              CASE ('idTbry(inorth,iS3CH)')
                idTbry(inorth,iS3CH)=varid
              CASE ('idTbry(iwest,iZ1_N)')
                idTbry(iwest,iZ1_N)=varid
              CASE ('idTbry(ieast,iZ1_N)')
                idTbry(ieast,iZ1_N)=varid
              CASE ('idTbry(isouth,iZ1_N)')
                idTbry(isouth,iZ1_N)=varid
              CASE ('idTbry(inorth,iZ1_N)')
                idTbry(inorth,iZ1_N)=varid
              CASE ('idTbry(iwest,iZ1_C)')
                idTbry(iwest,iZ1_C)=varid
              CASE ('idTbry(ieast,iZ1_C)')
                idTbry(ieast,iZ1_C)=varid
              CASE ('idTbry(isouth,iZ1_C)')
                idTbry(isouth,iZ1_C)=varid
              CASE ('idTbry(inorth,iZ1_C)')
                idTbry(inorth,iZ1_C)=varid
              CASE ('idTbry(iwest,iZ2_N)')
                idTbry(iwest,iZ2_N)=varid
              CASE ('idTbry(ieast,iZ2_N)')
                idTbry(ieast,iZ2_N)=varid
              CASE ('idTbry(isouth,iZ2_N)')
                idTbry(isouth,iZ2_N)=varid
              CASE ('idTbry(inorth,iZ2_N)')
                idTbry(inorth,iZ2_N)=varid
              CASE ('idTbry(iwest,iZ2_C)')
                idTbry(iwest,iZ2_C)=varid
              CASE ('idTbry(ieast,iZ2_C)')
                idTbry(ieast,iZ2_C)=varid
              CASE ('idTbry(isouth,iZ2_C)')
                idTbry(isouth,iZ2_C)=varid
              CASE ('idTbry(inorth,iZ2_C)')
                idTbry(inorth,iZ2_C)=varid
              CASE ('idTbry(iwest,iBAC_)')
                idTbry(iwest,iBAC_)=varid
              CASE ('idTbry(ieast,iBAC_)')
                idTbry(ieast,iBAC_)=varid
              CASE ('idTbry(isouth,iBAC_)')
                idTbry(isouth,iBAC_)=varid
              CASE ('idTbry(inorth,iBAC_)')
                idTbry(inorth,iBAC_)=varid
              CASE ('idTbry(iwest,iDD_N)')
                idTbry(iwest,iDD_N)=varid
              CASE ('idTbry(ieast,iDD_N)')
                idTbry(ieast,iDD_N)=varid
              CASE ('idTbry(isouth,iDD_N)')
                idTbry(isouth,iDD_N)=varid
              CASE ('idTbry(inorth,iDD_N)')
                idTbry(inorth,iDD_N)=varid
              CASE ('idTbry(iwest,iDD_C)')
                idTbry(iwest,iDD_C)=varid
              CASE ('idTbry(ieast,iDD_C)')
                idTbry(ieast,iDD_C)=varid
              CASE ('idTbry(isouth,iDD_C)')
                idTbry(isouth,iDD_C)=varid
              CASE ('idTbry(inorth,iDD_C)')
                idTbry(inorth,iDD_C)=varid
              CASE ('idTbry(iwest,iDDSi)')
                idTbry(iwest,iDDSi)=varid
              CASE ('idTbry(ieast,iDDSi)')
                idTbry(ieast,iDDSi)=varid
              CASE ('idTbry(isouth,iDDSi)')
                idTbry(isouth,iDDSi)=varid
              CASE ('idTbry(inorth,iDDSi)')
                idTbry(inorth,iDDSi)=varid
              CASE ('idTbry(iwest,iLDON)')
                idTbry(iwest,iLDON)=varid
              CASE ('idTbry(ieast,iLDON)')
                idTbry(ieast,iLDON)=varid
              CASE ('idTbry(isouth,iLDON)')
                idTbry(isouth,iLDON)=varid
              CASE ('idTbry(inorth,iLDON)')
                idTbry(inorth,iLDON)=varid
              CASE ('idTbry(iwest,iLDOC)')
                idTbry(iwest,iLDOC)=varid
              CASE ('idTbry(ieast,iLDOC)')
                idTbry(ieast,iLDOC)=varid
              CASE ('idTbry(isouth,iLDOC)')
                idTbry(isouth,iLDOC)=varid
              CASE ('idTbry(inorth,iLDOC)')
                idTbry(inorth,iLDOC)=varid
              CASE ('idTbry(iwest,iSDON)')
                idTbry(iwest,iSDON)=varid
              CASE ('idTbry(ieast,iSDON)')
                idTbry(ieast,iSDON)=varid
              CASE ('idTbry(isouth,iSDON)')
                idTbry(isouth,iSDON)=varid
              CASE ('idTbry(inorth,iSDON)')
                idTbry(inorth,iSDON)=varid
              CASE ('idTbry(iwest,iSDOC)')
                idTbry(iwest,iSDOC)=varid
              CASE ('idTbry(ieast,iSDOC)')
                idTbry(ieast,iSDOC)=varid
              CASE ('idTbry(isouth,iSDOC)')
                idTbry(isouth,iSDOC)=varid
              CASE ('idTbry(inorth,iSDOC)')
                idTbry(inorth,iSDOC)=varid
              CASE ('idTbry(iwest,iCLDC)')
                idTbry(iwest,iCLDC)=varid
              CASE ('idTbry(ieast,iCLDC)')
                idTbry(ieast,iCLDC)=varid
              CASE ('idTbry(isouth,iCLDC)')
                idTbry(isouth,iCLDC)=varid
              CASE ('idTbry(inorth,iCLDC)')
                idTbry(inorth,iCLDC)=varid
              CASE ('idTbry(iwest,iCSDC)')
                idTbry(iwest,iCSDC)=varid
              CASE ('idTbry(ieast,iCSDC)')
                idTbry(ieast,iCSDC)=varid
              CASE ('idTbry(isouth,iCSDC)')
                idTbry(isouth,iCSDC)=varid
              CASE ('idTbry(inorth,iCSDC)')
                idTbry(inorth,iCSDC)=varid
              CASE ('idTbry(iwest,iDDCA)')
                idTbry(iwest,iDDCA)=varid
              CASE ('idTbry(ieast,iDDCA)')
                idTbry(ieast,iDDCA)=varid
              CASE ('idTbry(isouth,iDDCA)')
                idTbry(isouth,iDDCA)=varid
              CASE ('idTbry(inorth,iDDCA)')
                idTbry(inorth,iDDCA)=varid

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
#ifdef CARBON
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
# ifdef IRON_LIMIT
              CASE ('idTbry(iwest,iS1_Fe)')
                idTbry(iwest,iS1_Fe)=varid
              CASE ('idTbry(ieast,iS1_Fe)')
                idTbry(ieast,iS1_Fe)=varid
              CASE ('idTbry(isouth,iS1_Fe)')
                idTbry(isouth,iS1_Fe)=varid
              CASE ('idTbry(inorth,iS1_Fe)')
                idTbry(inorth,iS1_Fe)=varid
              CASE ('idTbry(iwest,iS2_Fe)')
                idTbry(iwest,iS2_Fe)=varid
              CASE ('idTbry(ieast,iS2_Fe)')
                idTbry(ieast,iS2_Fe)=varid
              CASE ('idTbry(isouth,iS2_Fe)')
                idTbry(isouth,iS2_Fe)=varid
              CASE ('idTbry(inorth,iS2_Fe)')
                idTbry(inorth,iS2_Fe)=varid
              CASE ('idTbry(iwest,iS3_Fe)')
                idTbry(iwest,iS3_Fe)=varid
              CASE ('idTbry(ieast,iS3_Fe)')
                idTbry(ieast,iS3_Fe)=varid
              CASE ('idTbry(isouth,iS3_Fe)')
                idTbry(isouth,iS3_Fe)=varid
              CASE ('idTbry(inorth,iS3_Fe)')
                idTbry(inorth,iS3_Fe)=varid
              CASE ('idTbry(iwest,iFeD_)')
                idTbry(iwest,iFeD_)=varid
              CASE ('idTbry(ieast,iFeD_)')
                idTbry(ieast,iFeD_)=varid
              CASE ('idTbry(isouth,iFeD_)')
                idTbry(isouth,iFeD_)=varid
              CASE ('idTbry(inorth,iFeD_)')
                idTbry(inorth,iFeD_)=varid
# endif

#ifdef TS_PSOURCE

!
!  Biological tracers point Source/Sinks (river runoff).
!

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iSiOH)')
                idRtrc(iSiOH)=varid
              CASE ('idRtrc(iPO4_)')
                idRtrc(iPO4_)=varid
              CASE ('idRtrc(iS1_N)')
                idRtrc(iS1_N)=varid
              CASE ('idRtrc(iS1_C)')
                idRtrc(iS1_C)=varid
              CASE ('idRtrc(iS1CH)')
                idRtrc(iS1CH)=varid
              CASE ('idRtrc(iS2_N)')
                idRtrc(iS2_N)=varid
              CASE ('idRtrc(iS2_C)')
                idRtrc(iS2_C)=varid
              CASE ('idRtrc(iS2CH)')
                idRtrc(iS2CH)=varid
              CASE ('idRtrc(iS3_N)')
                idRtrc(iS3_N)=varid
              CASE ('idRtrc(iS3_C)')
                idRtrc(iS3_C)=varid
              CASE ('idRtrc(iS3CH)')
                idRtrc(iS3CH)=varid
              CASE ('idRtrc(iZ1_N)')
                idRtrc(iZ1_N)=varid
              CASE ('idRtrc(iZ1_C)')
                idRtrc(iZ1_C)=varid
              CASE ('idRtrc(iZ2_N)')
                idRtrc(iZ2_N)=varid
              CASE ('idRtrc(iZ2_C)')
                idRtrc(iZ2_C)=varid
              CASE ('idRtrc(iBAC_)')
                idRtrc(iBAC_)=varid
              CASE ('idRtrc(iDD_N)')
                idRtrc(iDD_N)=varid
              CASE ('idRtrc(iDD_C)')
                idRtrc(iDD_C)=varid
              CASE ('idRtrc(iDDSi)')
                idRtrc(iDDSi)=varid
              CASE ('idRtrc(iLDON)')
                idRtrc(iLDON)=varid
              CASE ('idRtrc(iLDOC)')
                idRtrc(iLDOC)=varid
              CASE ('idRtrc(iSDON)')
                idRtrc(iSDON)=varid
              CASE ('idRtrc(iSDOC)')
                idRtrc(iSDOC)=varid
              CASE ('idRtrc(iCLDC)')
                idRtrc(iCLDC)=varid
              CASE ('idRtrc(iCSDC)')
                idRtrc(iCSDC)=varid
              CASE ('idRtrc(iDDCA)')
                idRtrc(iDDCA)=varid

# ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
# endif
# ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
# endif
#endif
#ifdef DIAGNOSTICS_BIO
/*
**  Biological tracers term diagnostics.
*/
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
              CASE ('iDbio3(iPPro1)')
                iDbio3(iPPro1)=varid
              CASE ('iDbio3(iPPro2)')
                iDbio3(iPPro2)=varid
              CASE ('iDbio3(iPPro3)')
                iDbio3(iPPro3)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
#endif
