/*
** svn $Id: sediment_var.h 439 2010-01-11 19:29:29Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the sediment model variables that    **
**  are used in input and output NetCDF files.  The metadata          **
**  information is read from "varinfo.dat".                           **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state sediment tracers.
*/

              CASE ('idSbed(ithck)')
                idSbed(ithck)=varid
              CASE ('idSbed(iaged)')
                idSbed(iaged)=varid
              CASE ('idSbed(iporo)')
                idSbed(iporo)=varid
              CASE ('idSbed(idiff)')
                idSbed(idiff)=varid
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
              CASE ('idSbed(ibtcr)')
                idSbed(ibtcr)=varid
#endif
              CASE ('idBott(isd50)')
                idBott(isd50)=varid
              CASE ('idBott(idens)')
                idBott(idens)=varid
              CASE ('idBott(iwsed)')
                idBott(iwsed)=varid
              CASE ('idBott(itauc)')
                idBott(itauc)=varid
              CASE ('idBott(irlen)')
                idBott(irlen)=varid
              CASE ('idBott(irhgt)')
                idBott(irhgt)=varid
              CASE ('idBott(ibwav)')
                idBott(ibwav)=varid
              CASE ('idBott(izdef)')
                idBott(izdef)=varid
              CASE ('idBott(izapp)')
                idBott(izapp)=varid
              CASE ('idBott(izNik)')
                idBott(izNik)=varid
              CASE ('idBott(izbio)')
                idBott(izbio)=varid
              CASE ('idBott(izbfm)')
                idBott(izbfm)=varid
              CASE ('idBott(izbld)')
                idBott(izbld)=varid
              CASE ('idBott(izwbl)')
                idBott(izwbl)=varid
              CASE ('idBott(iactv)')
                idBott(iactv)=varid
              CASE ('idBott(ishgt)')
                idBott(ishgt)=varid
              CASE ('idBott(imaxD)')
                idBott(imaxD)=varid
              CASE ('idBott(idnet)')
                idBott(idnet)=varid
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
              CASE ('idBott(idoff)')
                idBott(idoff)=varid
              CASE ('idBott(idslp)')
                idBott(idslp)=varid
              CASE ('idBott(idtim)')
                idBott(idtim)=varid
              CASE ('idBott(idbmx)')
                idBott(idbmx)=varid
              CASE ('idBott(idbmm)')
                idBott(idbmm)=varid
              CASE ('idBott(idbzs)')
                idBott(idbzs)=varid
              CASE ('idBott(idbzm)')
                idBott(idbzm)=varid
              CASE ('idBott(idbzp)')
                idBott(idbzp)=varid
#endif
#if defined MIXED_BED
              CASE ('idBott(idprp)')
                idBott(idprp)=varid
#endif
              CASE ('idTvar(idmud(i))')
                load=.FALSE.
                IF (NCS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idTvar(idmud(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTvar(isand(i))')
                load=.FALSE.
                IF (NNS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idTvar(isand(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idfrac')
                load=.FALSE.
                IF ((NCS.gt.0).and.                                     &
     &              (Vinfo(1)(1:8).eq.'mudfrac_')) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idfrac(i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
                IF ((NNS.gt.0).and.                                     &
     &              (Vinfo(1)(1:9).eq.'sandfrac_')) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idfrac(NCS+i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idBmas')
                load=.FALSE.
                IF ((NCS.gt.0).and.                                     &
     &              (Vinfo(1)(1:8).eq.'mudmass_')) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idBmas(i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
                IF ((NNS.gt.0).and.                                     &
     &              (Vinfo(1)(1:9).eq.'sandmass_')) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idBmas(NCS+i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
#ifdef BEDLOAD
              CASE ('idUbld')
                load=.FALSE.
                IF ((NCS.gt.0).and.                                     &
     &              (Vinfo(1)(1:13).eq.'bedload_Umud_')) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idUbld(i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
                IF ((NNS.gt.0).and.                                     &
     &              (Vinfo(1)(1:14).eq.'bedload_Usand_')) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idUbld(NCS+i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idVbld')
                load=.FALSE.
                IF ((NCS.gt.0).and.                                     &
     &              (Vinfo(1)(1:13).eq.'bedload_Vmud_')) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idVbld(i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
                IF ((NNS.gt.0).and.                                     &
     &              (Vinfo(1)(1:14).eq.'bedload_Vsand_')) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idVbld(NCS+i)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
#endif

/*
**  Sediment tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,idmud(i))')
                load=.FALSE.
                IF (NCS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idTbry(iwest,idmud(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(ieast,idmud(i))')
                load=.FALSE.
                IF (NCS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idTbry(ieast,idmud(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(isouth,idmud(i))')
                load=.FALSE.
                IF (NCS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idTbry(isouth,idmud(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(inorth,idmud(i))')
                load=.FALSE.
                IF (NCS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idTbry(inorth,idmud(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(iwest,isand(i))')
                load=.FALSE.
                IF (NNS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idTbry(iwest,isand(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(ieast,isand(i))')
                load=.FALSE.
                IF (NNS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idTbry(ieast,isand(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(isouth,isand(i))')
                load=.FALSE.
                IF (NNS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idTbry(isouth,isand(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idTbry(inorth,isand(i))')
                load=.FALSE.
                IF (NNS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idTbry(inorth,isand(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF


/*
**  Sediment tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(idmud)')
                load=.FALSE.
                IF (NCS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NCS
                    varid=varid+1
                    idRtrc(idmud(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
              CASE ('idRtrc(isand)')
                load=.FALSE.
                IF (NNS.gt.0) THEN
                  varid=varid-1
                  DO i=1,NNS
                    varid=varid+1
                    idRtrc(isand(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i2.2)')                   &
     &                    TRIM(ADJUSTL(Vinfo(1))), i
                    WRITE (Vname(2,varid),'(a,a,i2.2)')                 &
     &                    TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                  varid=varid+1
                END IF
