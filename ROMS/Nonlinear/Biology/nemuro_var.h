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

/*
**  Model state biological tracers.
*/

              CASE ('idTvar(iLphy)')
                idTvar(iLphy)=varid
              CASE ('idTvar(iSphy)')
                idTvar(iSphy)=varid
              CASE ('idTvar(iLzoo)')
                idTvar(iLzoo)=varid
              CASE ('idTvar(iSzoo)')
                idTvar(iSzoo)=varid
              CASE ('idTvar(iPzoo)')
                idTvar(iPzoo)=varid
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPON_)')
                idTvar(iPON_)=varid
              CASE ('idTvar(iDON_)')
                idTvar(iDON_)=varid
              CASE ('idTvar(iSiOH)')
                idTvar(iSiOH)=varid
              CASE ('idTvar(iopal)')
                idTvar(iopal)=varid
# ifdef IRON_LIMIT
              CASE ('idTvar(iFeSp)')
                idTvar(iFeSp)=varid
              CASE ('idTvar(iFeLp)')
                idTvar(iFeLp)=varid
              CASE ('idTvar(iFeD_)')
                idTvar(iFeD_)=varid
# endif
# ifdef NEMURO_SED1
              CASE ('idPONsed')
                idPONsed=varid
              CASE ('idOPALsed')
                idOPALsed=varid
              CASE ('idDENITsed')
                idDENITsed=varid
              CASE ('idPONbur')
                idPONbur=varid
              CASE ('idOPALbur')
                idOPALbur=varid
# endif
# ifdef PRIMARY_PROD
              CASE ('idNPP')
                idNPP=varid
# endif
!
!  Biological tracers open boundary conditions.
!

              CASE ('idTbry(iwest,iLphy)')
                idTbry(iwest,iLphy)=varid
              CASE ('idTbry(ieast,iLphy)')
                idTbry(ieast,iLphy)=varid
              CASE ('idTbry(isouth,iLphy)')
                idTbry(isouth,iLphy)=varid
              CASE ('idTbry(inorth,iLphy)')
                idTbry(inorth,iLphy)=varid

              CASE ('idTbry(iwest,iSphy)')
                idTbry(iwest,iSphy)=varid
              CASE ('idTbry(ieast,iSphy)')
                idTbry(ieast,iSphy)=varid
              CASE ('idTbry(isouth,iSphy)')
                idTbry(isouth,iSphy)=varid
              CASE ('idTbry(inorth,iSphy)')
                idTbry(inorth,iSphy)=varid

              CASE ('idTbry(iwest,iLzoo)')
                idTbry(iwest,iLzoo)=varid
              CASE ('idTbry(ieast,iLzoo)')
                idTbry(ieast,iLzoo)=varid
              CASE ('idTbry(isouth,iLzoo)')
                idTbry(isouth,iLzoo)=varid
              CASE ('idTbry(inorth,iLzoo)')
                idTbry(inorth,iLzoo)=varid

              CASE ('idTbry(iwest,iSzoo)')
                idTbry(iwest,iSzoo)=varid
              CASE ('idTbry(ieast,iSzoo)')
                idTbry(ieast,iSzoo)=varid
              CASE ('idTbry(isouth,iSzoo)')
                idTbry(isouth,iSzoo)=varid
              CASE ('idTbry(inorth,iSzoo)')
                idTbry(inorth,iSzoo)=varid

              CASE ('idTbry(iwest,iPzoo)')
                idTbry(iwest,iPzoo)=varid
              CASE ('idTbry(ieast,iPzoo)')
                idTbry(ieast,iPzoo)=varid
              CASE ('idTbry(isouth,iPzoo)')
                idTbry(isouth,iPzoo)=varid
              CASE ('idTbry(inorth,iPzoo)')
                idTbry(inorth,iPzoo)=varid

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

              CASE ('idTbry(iwest,iPON_)')
                idTbry(iwest,iPON_)=varid
              CASE ('idTbry(ieast,iPON_)')
                idTbry(ieast,iPON_)=varid
              CASE ('idTbry(isouth,iPON_)')
                idTbry(isouth,iPON_)=varid
              CASE ('idTbry(inorth,iPON_)')
                idTbry(inorth,iPON_)=varid

              CASE ('idTbry(iwest,iDON_)')
                idTbry(iwest,iDON_)=varid
              CASE ('idTbry(ieast,iDON_)')
                idTbry(ieast,iDON_)=varid
              CASE ('idTbry(isouth,iDON_)')
                idTbry(isouth,iDON_)=varid
              CASE ('idTbry(inorth,iDON_)')
                idTbry(inorth,iDON_)=varid

              CASE ('idTbry(iwest,iSiOH)')
                idTbry(iwest,iSiOH)=varid
              CASE ('idTbry(ieast,iSiOH)')
                idTbry(ieast,iSiOH)=varid
              CASE ('idTbry(isouth,iSiOH)')
                idTbry(isouth,iSiOH)=varid
              CASE ('idTbry(inorth,iSiOH)')
                idTbry(inorth,iSiOH)=varid

              CASE ('idTbry(iwest,iopal)')
                idTbry(iwest,iopal)=varid
              CASE ('idTbry(ieast,iopal)')
                idTbry(ieast,iopal)=varid
              CASE ('idTbry(isouth,iopal)')
                idTbry(isouth,iopal)=varid
              CASE ('idTbry(inorth,iopal)')
                idTbry(inorth,iopal)=varid
# ifdef IRON_LIMIT
              CASE ('idTbry(iwest,iFeSp)')
                idTbry(iwest,iFeSp)=varid
              CASE ('idTbry(ieast,iFeSp)')
                idTbry(ieast,iFeSp)=varid
              CASE ('idTbry(isouth,iFeSp)')
                idTbry(isouth,iFeSp)=varid
              CASE ('idTbry(inorth,iFeSp)')
                idTbry(inorth,iFeSp)=varid
              CASE ('idTbry(iwest,iFeLp)')
                idTbry(iwest,iFeLp)=varid
              CASE ('idTbry(ieast,iFeLp)')
                idTbry(ieast,iFeLp)=varid
              CASE ('idTbry(isouth,iFeLp)')
                idTbry(isouth,iFeLp)=varid
              CASE ('idTbry(inorth,iFeLp)')
                idTbry(inorth,iFeLp)=varid
              CASE ('idTbry(iwest,iFeD_)')
                idTbry(iwest,iFeD_)=varid
              CASE ('idTbry(ieast,iFeD_)')
                idTbry(ieast,iFeD_)=varid
              CASE ('idTbry(isouth,iFeD_)')
                idTbry(isouth,iFeD_)=varid
              CASE ('idTbry(inorth,iFeD_)')
                idTbry(inorth,iFeD_)=varid
# endif

!
! Biological tracers point Source/Sinks (river runoff).
!
              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iDON_)')
                idRtrc(iDON_)=varid
              CASE ('idRtrc(iPON_)')
                idRtrc(iPON_)=varid
