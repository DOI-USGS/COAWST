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

              CASE ('idTvar(iNO3)')
                idTvar(iNO3)=varid
              CASE ('idTvar(iNH4)')
                idTvar(iNH4)=varid
              CASE ('idTvar(iPhS)')
                idTvar(iPhS)=varid
              CASE ('idTvar(iPhL)')
                idTvar(iPhL)=varid
              CASE ('idTvar(iMZS)')
                idTvar(iMZS)=varid
              CASE ('idTvar(iMZL)')
                idTvar(iMZL)=varid
              CASE ('idTvar(iCop)')
                idTvar(iCop)=varid
              CASE ('idTvar(iNCaS)')
                idTvar(iNCaS)=varid
              CASE ('idTvar(iEupS)')
                idTvar(iEupS)=varid
              CASE ('idTvar(iNCaO)')
                idTvar(iNCaO)=varid
              CASE ('idTvar(iEupO)')
                idTvar(iEupO)=varid
              CASE ('idTvar(iDet)')
                idTvar(iDet)=varid
              CASE ('idTvar(iDetF)')
                idTvar(iDetF)=varid
#ifdef JELLY
              CASE ('idTvar(iJel)')
                idTvar(iJel)=varid
#endif
#ifdef IRON_LIMIT
              CASE ('idTvar(iFe)')
                idTvar(iFe)=varid
#endif

#ifdef BENTHIC
              CASE ('idBeTvar(iBen)')
               idBeTvar(iBen)=varid
              CASE ('idBeTvar(iBenDet)')
                idBeTvar(iBenDet)=varid
#endif
#ifdef ICE_BIO
# ifdef CLIM_ICE_1D
              CASE ('idIceBvar(iIcePhL)')
                idIceBvar(iIcePhL)=varid
              CASE ('idIceBvar(iIceNO3)')
                idIceBvar(iIceNO3)=varid
              CASE ('idIceBvar(iIceNH4)')
                idIceBvar(iIceNH4)=varid
              CASE ('idIceBvar(iIceZ)')
                idIceBvar(iIceZ)=varid
              CASE ('idIceBvar(iIceLog)')
                idIceBvar(iIceLog)=varid
              CASE ('idIceLogvar(iIceLog)')
                idIceLogvar(iIceLog)=varid
# elif defined BERING_10K 
              CASE ('idIcePhL')
                idIcePhL=varid

              CASE ('idIcePhLbc(iwest)')
                idIcePhLbc(iwest)=varid
              CASE ('idIcePhLbc(ieast)')
                idIcePhLbc(ieast)=varid
              CASE ('idIcePhLbc(isouth)')
                idIcePhLbc(isouth)=varid
              CASE ('idIcePhLbc(inorth)')
                idIcePhLbc(inorth)=varid
               
              CASE ('idIceNO3')
                idIceNO3=varid
               
              CASE ('idIceNO3bc(iwest)')
                idIceNO3bc(iwest)=varid
              CASE ('idIceNO3bc(ieast)')
                idIceNO3bc(ieast)=varid
              CASE ('idIceNO3bc(isouth)')
                idIceNO3bc(isouth)=varid
              CASE ('idIceNO3bc(inorth)')
                idIceNO3bc(inorth)=varid
               
              CASE ('idIceNH4')
                idIceNH4=varid
               
              CASE ('idIceNH4bc(iwest)')
                idIceNH4bc(iwest)=varid
              CASE ('idIceNH4bc(ieast)')
                idIceNH4bc(ieast)=varid
              CASE ('idIceNH4bc(isouth)')
                idIceNH4bc(isouth)=varid
              CASE ('idIceNH4bc(inorth)')
                idIceNH4bc(inorth)=varid

              CASE ('idIceLog')
                idIceLog=varid
# endif
#endif
#ifdef STATIONARY
             !-----------------------------
             !3D Stationary Array
             !-----------------------------
              CASE ('idTSvar(i3Stat1)')
                idTSvar(i3Stat1) = varid
              CASE ('idTSvar(i3Stat2)')
                idTSvar(i3Stat2) = varid
              CASE ('idTSvar(i3Stat3)')
                idTSvar(i3Stat3) = varid
              CASE ('idTSvar(i3Stat4)')
                idTSvar(i3Stat4) = varid
              CASE ('idTSvar(i3Stat5)')
                idTSvar(i3Stat5) = varid
              CASE ('idTSvar(i3Stat6)')
                idTSvar(i3Stat6) = varid
              CASE ('idTSvar(i3Stat7)')
                idTSvar(i3Stat7) = varid

              CASE ('idTSvar(i3Stat8)')
                idTSvar(i3Stat8) = varid
              CASE ('idTSvar(i3Stat9)')
                idTSvar(i3Stat9) = varid
              CASE ('idTSvar(i3Stat10)')
                idTSvar(i3Stat10) = varid
              CASE ('idTSvar(i3Stat11)')
                idTSvar(i3Stat11) = varid
              CASE ('idTSvar(i3Stat12)')
                idTSvar(i3Stat12) = varid
              CASE ('idTSvar(i3Stat13)')
                idTSvar(i3Stat13) = varid
              CASE ('idTSvar(i3Stat14)')
                idTSvar(i3Stat14) = varid
              CASE ('idTSvar(i3Stat15)')
                idTSvar(i3Stat15) = varid
              CASE ('idTSvar(i3Stat16)')
                idTSvar(i3Stat16) = varid
#endif
#ifdef STATIONARY2
             !-----------------------------
             !2D Stationary Array
             !-----------------------------
              CASE ('idTS2var(i2Stat1)')
                idTS2var(i2Stat1) = varid
              CASE ('idTS2var(i2Stat2)')
                idTS2var(i2Stat2) = varid
              CASE ('idTS2var(i2Stat3)')
                idTS2var(i2Stat3) = varid
              CASE ('idTS2var(i2Stat4)')
                idTS2var(i2Stat4) = varid
              CASE ('idTS2var(i2Stat5)')
                idTS2var(i2Stat5) = varid
              CASE ('idTS2var(i2Stat6)')
                idTS2var(i2Stat6) = varid
              CASE ('idTS2var(i2Stat7)')
                idTS2var(i2Stat7) = varid
              CASE ('idTS2var(i2Stat8)')
                idTS2var(i2Stat8) = varid

#endif
#ifdef PROD3
             !-----------------------------
             !3D Production Array
             !-----------------------------
              CASE ('idPT3var(iPhSprd)')
                idPT3var(iPhSprd) = varid
              CASE ('idPT3var(iPhLprd)')
                idPT3var(iPhLprd) = varid
              CASE ('idPT3var(iMZSprd)')
                idPT3var(iMZSprd) = varid
              CASE ('idPT3var(iMZLprd)')
                idPT3var(iMZLprd) = varid
              CASE ('idPT3var(iCopPrd)')
                idPT3var(iCopPrd) = varid
              CASE ('idPT3var(iNCaPrd)')
                idPT3var(iNCaPrd) = varid
              CASE ('idPT3var(iEupPrd)')
                idPT3var(iEupPrd) = varid

# ifdef JELLY
              CASE ('idPT3var(iJelPrd)')
                idPT3var(iJelPrd) = varid
# endif
#endif

#ifdef PROD2
             !-----------------------------
             !2D Production Array
             !-----------------------------
              CASE ('idPT2var(iBenPrd)')
                idPT2var(iBenPrd) = varid
              CASE ('idPT2var(iIAPrd)')
                idPT2var(iIAPrd) = varid
              CASE ('idPT2var(iXPrd)')
                idPT2var(iXPrd) = varid
#endif
#ifdef BIOFLUX
             !-----------------------------
             !bio flux Array
             !-----------------------------
              CASE ('idTBFvar(iBF)')
                idTBFvar(iBF) = varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iNO3)')
                idTbry(iwest,iNO3)=varid
              CASE ('idTbry(ieast,iNO3)')
                idTbry(ieast,iNO3)=varid
              CASE ('idTbry(isouth,iNO3)')
                idTbry(isouth,iNO3)=varid
              CASE ('idTbry(inorth,iNO3)')
                idTbry(inorth,iNO3)=varid

              CASE ('idTbry(iwest,iNH4)')
                idTbry(iwest,iNH4)=varid
              CASE ('idTbry(ieast,iNH4)')
                idTbry(ieast,iNH4)=varid
              CASE ('idTbry(isouth,iNH4)')
                idTbry(isouth,iNH4)=varid
              CASE ('idTbry(inorth,iNH4)')
                idTbry(inorth,iNH4)=varid

              CASE ('idTbry(iwest,iPhS)')
                idTbry(iwest,iPhS)=varid
              CASE ('idTbry(ieast,iPhS)')
                idTbry(ieast,iPhS)=varid
              CASE ('idTbry(isouth,iPhS)')
                idTbry(isouth,iPhS)=varid
              CASE ('idTbry(inorth,iPhS)')
                idTbry(inorth,iPhS)=varid

              CASE ('idTbry(iwest,iPhL)')
                idTbry(iwest,iPhL)=varid
              CASE ('idTbry(ieast,iPhL)')
                idTbry(ieast,iPhL)=varid
              CASE ('idTbry(isouth,iPhL)')
                idTbry(isouth,iPhL)=varid
              CASE ('idTbry(inorth,iPhL)')
                idTbry(inorth,iPhL)=varid

              CASE ('idTbry(iwest,iMZS)')
                idTbry(iwest,iMZS)=varid
              CASE ('idTbry(ieast,iMZS)')
                idTbry(ieast,iMZS)=varid
              CASE ('idTbry(isouth,iMZS)')
                idTbry(isouth,iMZS)=varid
              CASE ('idTbry(inorth,iMZS)')
                idTbry(inorth,iMZS)=varid

              CASE ('idTbry(iwest,iMZL)')
                idTbry(iwest,iMZL)=varid
              CASE ('idTbry(ieast,iMZL)')
                idTbry(ieast,iMZL)=varid
              CASE ('idTbry(isouth,iMZL)')
                idTbry(isouth,iMZL)=varid
              CASE ('idTbry(inorth,iMZL)')
                idTbry(inorth,iMZL)=varid

              CASE ('idTbry(iwest,iCop)')
                idTbry(iwest,iCop)=varid
              CASE ('idTbry(ieast,iCop)')
                idTbry(ieast,iCop)=varid
              CASE ('idTbry(isouth,iCop)')
                idTbry(isouth,iCop)=varid
              CASE ('idTbry(inorth,iCop)')
                idTbry(inorth,iCop)=varid

              CASE ('idTbry(iwest,iNCaS)')
                idTbry(iwest,iNCaS)=varid
              CASE ('idTbry(ieast,iNCaS)')
                idTbry(ieast,iNCaS)=varid
              CASE ('idTbry(isouth,iNCaS)')
                idTbry(isouth,iNCaS)=varid
              CASE ('idTbry(inorth,iNCaS)')
                idTbry(inorth,iNCaS)=varid

              CASE ('idTbry(iwest,iEupS)')
                idTbry(iwest,iEupS)=varid
              CASE ('idTbry(ieast,iEupS)')
                idTbry(ieast,iEupS)=varid
              CASE ('idTbry(isouth,iEupS)')
                idTbry(isouth,iEupS)=varid
              CASE ('idTbry(inorth,iEupS)')
                idTbry(inorth,iEupS)=varid

              CASE ('idTbry(iwest,iNCaO)')
                idTbry(iwest,iNCaO)=varid
              CASE ('idTbry(ieast,iNCaO)')
                idTbry(ieast,iNCaO)=varid
              CASE ('idTbry(isouth,iNCaO)')
                idTbry(isouth,iNCaO)=varid
              CASE ('idTbry(inorth,iNCaO)')
                idTbry(inorth,iNCaO)=varid

              CASE ('idTbry(iwest,iEupO)')
                idTbry(iwest,iEupO)=varid
              CASE ('idTbry(ieast,iEupO)')
                idTbry(ieast,iEupO)=varid
              CASE ('idTbry(isouth,iEupO)')
                idTbry(isouth,iEupO)=varid
              CASE ('idTbry(inorth,iEupO)')
                idTbry(inorth,iEupO)=varid

              CASE ('idTbry(iwest,iDet)')
                idTbry(iwest,iDet)=varid
              CASE ('idTbry(ieast,iDet)')
                idTbry(ieast,iDet)=varid
              CASE ('idTbry(isouth,iDet)')
                idTbry(isouth,iDet)=varid
              CASE ('idTbry(inorth,iDet)')
                idTbry(inorth,iDet)=varid

              CASE ('idTbry(iwest,iDetF)')
                idTbry(iwest,iDetF)=varid
              CASE ('idTbry(ieast,iDetF)')
                idTbry(ieast,iDetF)=varid
              CASE ('idTbry(isouth,iDetF)')
                idTbry(isouth,iDetF)=varid
              CASE ('idTbry(inorth,iDetF)')
                idTbry(inorth,iDetF)=varid

#ifdef JELLY
              CASE ('idTbry(iwest,iJel)')
                idTbry(iwest,iJel)=varid
              CASE ('idTbry(ieast,iJel)')
                idTbry(ieast,iJel)=varid
              CASE ('idTbry(isouth,iJel)')
                idTbry(isouth,iJel)=varid
              CASE ('idTbry(inorth,iJel)')
                idTbry(inorth,iJel)=varid
#endif

#ifdef IRON_LIMIT
              CASE ('idTbry(iwest,iFe)')
                idTbry(iwest,iFe)=varid
              CASE ('idTbry(ieast,iFe)')
                idTbry(ieast,iFe)=varid
              CASE ('idTbry(isouth,iFe)')
                idTbry(isouth,iFe)=varid
              CASE ('idTbry(inorth,iFe)')
                idTbry(inorth,iFe)=varid
#endif


#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iNO3)')
                idRtrc(iNO3)=varid
#endif
