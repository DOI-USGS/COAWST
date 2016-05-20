      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in Fennel et al. (2006) ecosystem model input    !
!  parameters. They are specified in input script "estuarybgc.in".     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status

      integer :: decode_line, load_i, load_l, load_lbc, load_r

      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(200) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in Fennel et al. (2006) biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lbiology')
              Npts=load_l(Nval, Cval, Ngrids, Lbiology)
            CASE ('BioIter')
              Npts=load_i(Nval, Rval, Ngrids, BioIter)
            CASE ('AttSW')
              Npts=load_r(Nval, Rval, Ngrids, AttSW)
            CASE ('AttChl')
              Npts=load_r(Nval, Rval, Ngrids, AttChl)
            CASE ('PARfrac')
              Npts=load_r(Nval, Rval, Ngrids, PARfrac)
            CASE ('Vp0')
              Npts=load_r(Nval, Rval, Ngrids, Vp0)
            CASE ('I_thNH4')
              Npts=load_r(Nval, Rval, Ngrids, I_thNH4)
            CASE ('D_p5NH4')
              Npts=load_r(Nval, Rval, Ngrids, D_p5NH4)
            CASE ('NitriR')
              Npts=load_r(Nval, Rval, Ngrids, NitriR)
            CASE ('K_NO3')
              Npts=load_r(Nval, Rval, Ngrids, K_NO3)
            CASE ('K_NH4')
              Npts=load_r(Nval, Rval, Ngrids, K_NH4)
            CASE ('K_Phy')
              Npts=load_r(Nval, Rval, Ngrids, K_Phy)
            CASE ('Chl2C_m')
              Npts=load_r(Nval, Rval, Ngrids, Chl2C_m)
            CASE ('ChlMin')
              Npts=load_r(Nval, Rval, Ngrids, ChlMin)
            CASE ('PhyCN')
              Npts=load_r(Nval, Rval, Ngrids, PhyCN)
            CASE ('PhyIP')
              Npts=load_r(Nval, Rval, Ngrids, PhyIP)
            CASE ('PhyIS')
              Npts=load_r(Nval, Rval, Ngrids, PhyIS)
            CASE ('PhyMin')
              Npts=load_r(Nval, Rval, Ngrids, PhyMin)
            CASE ('PhyMR')
              Npts=load_r(Nval, Rval, Ngrids, PhyMR)
            CASE ('ZooAE_N')
              Npts=load_r(Nval, Rval, Ngrids, ZooAE_N)
            CASE ('ZooBM')
              Npts=load_r(Nval, Rval, Ngrids, ZooBM)
            CASE ('ZooCN')
              Npts=load_r(Nval, Rval, Ngrids, ZooCN)
            CASE ('ZooER')
              Npts=load_r(Nval, Rval, Ngrids, ZooER)
            CASE ('ZooGR')
              Npts=load_r(Nval, Rval, Ngrids, ZooGR)
            CASE ('ZooMin')
              Npts=load_r(Nval, Rval, Ngrids, ZooMin)
            CASE ('ZooMR')
              Npts=load_r(Nval, Rval, Ngrids, ZooMR)
            CASE ('LDeRRN')
              Npts=load_r(Nval, Rval, Ngrids, LDeRRN)
            CASE ('LDeRRC')
              Npts=load_r(Nval, Rval, Ngrids, LDeRRC)
            CASE ('CoagR')
              Npts=load_r(Nval, Rval, Ngrids, CoagR)
            CASE ('SDeRRN')
              Npts=load_r(Nval, Rval, Ngrids, SDeRRN)
            CASE ('SDeRRC')
              Npts=load_r(Nval, Rval, Ngrids, SDeRRC)
            CASE ('wPhy')
              Npts=load_r(Nval, Rval, Ngrids, wPhy)
            CASE ('wLDet')
              Npts=load_r(Nval, Rval, Ngrids, wLDet)
            CASE ('wSDet')
              Npts=load_r(Nval, Rval, Ngrids, wSDet)
            CASE ('pCO2air')
              Npts=load_r(Nval, Rval, Ngrids, pCO2air)
#if defined SPECTRAL_LIGHT
            CASE ('SIGATRB')
              Npts=load_r(Nval, Rval, Ngrids, SIGATRB)
            CASE ('STRB')
              Npts=load_r(Nval, Rval, Ngrids, STRB)
            CASE ('BLTRB')
              Npts=load_r(Nval, Rval, Ngrids, BLTRB)
            CASE ('SIGBTRB')
              Npts=load_r(Nval, Rval, Ngrids, SIGBTRB)
            CASE ('ETASPEC')
              Npts=load_r(Nval, Rval, Ngrids, ETASPEC)
            CASE ('BB2B')
              Npts=load_r(Nval, Rval, Ngrids, BB2B)
#endif
#if defined SAV_BIOMASS
            CASE ('GMODopt')
              Npts=load_i(Nval, Rval, Ngrids, GMODopt)
            CASE ('KNSED')
              Npts=load_r(Nval, Rval, Ngrids, KNSED)
            CASE ('KNWC')
              Npts=load_r(Nval, Rval, Ngrids, KNWC)
            CASE ('TOPT')
              Npts=load_r(Nval, Rval, Ngrids, TOPT)
            CASE ('THTA')
              Npts=load_r(Nval, Rval, Ngrids, THTA)
            CASE ('THTA2')
              Npts=load_r(Nval, Rval, Ngrids, THTA2)
            CASE ('SCL')
              Npts=load_r(Nval, Rval, Ngrids, SCL)
            CASE ('SCL2')
              Npts=load_r(Nval, Rval, Ngrids, SCL2)
            CASE ('KI')
              Npts=load_r(Nval, Rval, Ngrids, KI)
            CASE ('SR')
              Npts=load_r(Nval, Rval, Ngrids, SR)
            CASE ('LMBAMX')
              Npts=load_r(Nval, Rval, Ngrids, LMBAMX)
            CASE ('KMAG')
              Npts=load_r(Nval, Rval, Ngrids, KMAG)
            CASE ('ARSC')
              Npts=load_r(Nval, Rval, Ngrids, ARSC)
            CASE ('ARC')
              Npts=load_r(Nval, Rval, Ngrids, ARC)
            CASE ('BSRC')
              Npts=load_r(Nval, Rval, Ngrids, BSRC)
            CASE ('RC')
              Npts=load_r(Nval, Rval, Ngrids, RC)
            CASE ('RtStTL')
              Npts=load_r(Nval, Rval, Ngrids, RtStTL)
            CASE ('DOWNt')
              Npts=load_r(Nval, Rval, Ngrids, DOWNt)
            CASE ('TRNS')
              Npts=load_r(Nval, Rval, Ngrids, TRNS)
            CASE ('TCRIT')
              Npts=load_r(Nval, Rval, Ngrids, TCRIT)
            CASE ('KM')
              Npts=load_r(Nval, Rval, Ngrids, KM)
#endif 
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  nl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  nl_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu2(i,ng)=Rbio(itrc,ng)
                  tl_tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_tnu4(i,ng)=Rbio(itrc,ng)
                  ad_tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('LtracerSponge')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSponge(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Akt_bak(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('ad_AKT_fac')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  ad_Akt_fac(i,ng)=Rbio(itrc,ng)
                  tl_Akt_fac(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Tnudg(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idbio(iTrcStr), idbio(iTrcEnd),             &
     &                      Vname(1,idTvar(idbio(itracer))), LBC)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idbio(iTrcStr), idbio(iTrcEnd),             &
     &                      Vname(1,idTvar(idbio(itracer))), ad_LBC)
#endif
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('LnudgeTCLM')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LnudgeTCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTvar)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idTvar(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTsur)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,30)                          &
     &                                'idTsur(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            CASE ('Aout(idTvar)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idTTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idUTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idVTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idVTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHUTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=iHUTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHVTav)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=iHVTav(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
#ifdef DIAGNOSTICS_TS
            CASE ('Dout(iTrate)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTrate),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iThadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iThadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTxadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTyadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTvadv)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltrc(i,ng)
                END DO
              END DO
# if defined TS_DIF2 || defined TS_DIF4
            CASE ('Dout(iThdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iThdif),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTxdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltrc(i,ng)
                END DO
              END DO
            CASE ('Dout(iTydif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTydif),ng)=Ltrc(i,ng)
                END DO
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            CASE ('Dout(iTsdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltrc(i,ng)
                END DO
              END DO
#  endif
# endif
            CASE ('Dout(iTvdif)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO i=1,NBT
                  itrc=idbio(i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltrc(i,ng)
                END DO
              END DO
#endif
#ifdef DIAGNOSTICS_BIO
# ifdef CARBON
            CASE ('Dout(iCOfx)')
              IF (iDbio2(iCOfx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(iCOfx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(iCOfx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
# ifdef DENITRIFICATION
            CASE ('Dout(iDNIT)')
              IF (iDbio2(iDNIT).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(iDNIT)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(iDNIT)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
# ifdef CARBON
            CASE ('Dout(ipCO2)')
              IF (iDbio2(ipCO2).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(ipCO2)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(ipCO2)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
# ifdef OXYGEN
            CASE ('Dout(iO2fx)')
              IF (iDbio2(iO2fx).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio2(iO2fx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio2(iO2fx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
            CASE ('Dout(iPPro)')
              IF (iDbio3(iPPro).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iPPro)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iPPro)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iNO3u)')
              IF (iDbio3(iNO3u).eq.0) THEN
                IF (Master) WRITE (out,40) 'iDbio3(iNO3u)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iNO3u)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
#endif
#if defined SPECTRAL_LIGHT 
            CASE ('Hout(idPARo)') 
              IF (idPARo.eq.0) THEN
                IF (Master) WRITE (out,40) 'idPARo'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idPARo,1:Ngrids))
            CASE ('Hout(idPARs)') 
              IF (idPARs.eq.0) THEN
                IF (Master) WRITE (out,40) 'idPARs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idPARs,1:Ngrids))
            CASE ('Hout(idSpKd)') 
              IF (idSpKd.eq.0) THEN
                IF (Master) WRITE (out,40) 'idSpKd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSpKd,1:Ngrids))
#endif
#if defined SAV_BIOMASS
            CASE ('Hout(iddins)')
              IF (iddins.eq.0) THEN
                IF (Master) WRITE (out,40) 'iddins'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(iddins,1:Ngrids))
            CASE ('Hout(iddinw)')
              IF (iddinw.eq.0) THEN
                IF (Master) WRITE (out,40) 'iddinw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(iddinw,1:Ngrids))
            CASE ('Hout(iddowc)')
              IF (iddowc.eq.0) THEN
                IF (Master) WRITE (out,40) 'iddowc'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(iddowc,1:Ngrids))
            CASE ('Hout(idwsvl)')
              IF (idwsvl.eq.0) THEN
                IF (Master) WRITE (out,40) 'idwsvl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idwsvl,1:Ngrids))
            CASE ('Hout(idsagb)')
              IF (idsagb.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsagb'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsagb,1:Ngrids))
            CASE ('Hout(idsbgb)')
              IF (idsbgb.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsbgb'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsbgb,1:Ngrids))
            CASE ('Hout(idsvpp)')
              IF (idsvpp.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsvpp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsvpp,1:Ngrids))
            CASE ('Hout(idsvam)')
              IF (idsvam.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsvam'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsvam,1:Ngrids))
            CASE ('Hout(idsgar)')
              IF (idsgar.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsgar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsgar,1:Ngrids))
            CASE ('Hout(idsvbr)')
              IF (idsvbr.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsvbr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsvbr,1:Ngrids))
            CASE ('Hout(idsvrs)')
              IF (idsvrs.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsvrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsvrs,1:Ngrids))
            CASE ('Hout(idsvbg)')
              IF (idsvbg.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsvbg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsvbg,1:Ngrids))
            CASE ('Hout(idsvag)')
              IF (idsvag.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsvag'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsvag,1:Ngrids))
            CASE ('Hout(idsbgr)')
              IF (idsbgr.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsbgr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsbgr,1:Ngrids))
            CASE ('Hout(idsbgm)')
              IF (idsbgm.eq.0) THEN
                IF (Master) WRITE (out,40) 'idsbgm'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idsbgm,1:Ngrids))
#endif
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,60) ng
            WRITE (out,70) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,80) AttSW(ng), 'AttSW',                          &
     &            'Light attenuation of seawater (m-1).'
            WRITE (out,80) AttChl(ng), 'AttChl',                        &
     &            'Light attenuation by chlorophyll (1/(mg_Chl m-2)).'
            WRITE (out,90) PARfrac(ng), 'PARfrac',                      &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,90) Vp0(ng), 'Vp0',                              &
     &            'Eppley temperature-limited growth parameter',        &
     &            '(nondimensional).'
            WRITE (out,80) I_thNH4(ng), 'I_thNH4',                      &
     &            'Radiation threshold for nitrification (W/m2).'
            WRITE (out,80) D_p5NH4(ng), 'D_p5NH4',                      &
     &            'Half-saturation radiation for nitrification (W/m2).'
            WRITE (out,80) NitriR(ng), 'NitriR',                        &
     &            'Nitrification rate (day-1).'
            WRITE (out,90) K_NO3(ng), 'K_NO3',                          &
     &            'Inverse half-saturation for phytoplankton NO3',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,90) K_NH4(ng), 'K_NH4',                          &
     &            'Inverse half-saturation for phytoplankton NH4',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,90) K_Phy(ng), 'K_Phy',                          &
     &            'Zooplankton half-saturation constant for ingestion', &
     &            '(mmol_N m-3)^2.'
            WRITE (out,80) Chl2C_m(ng), 'Chl2C_m',                      &
     &            'Maximum chlorophyll to carbon ratio (mg_Chl/mg_C).'
            WRITE (out,80) ChlMin(ng), 'ChlMin',                        &
     &            'Chlorophyll minimum threshold (mg_Chl/m3).'
            WRITE (out,80) PhyCN(ng), 'PhyCN',                          &
     &            'Phytoplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,80) PhyIP(ng), 'PhyIP',                          &
     &            'Phytoplankton NH4 inhibition parameter (1/mmol_N).'
            WRITE (out,90) PhyIS(ng), 'PhyIS',                          &
     &            'Phytoplankton growth, initial slope of P-I curve',   &
     &            '(mg_C/(mg_Chl Watts m-2 day)).'
            WRITE (out,80) PhyMin(ng), 'PhyMin',                        &
     &            'Phytoplankton minimum threshold (mmol_N/m3).'
            WRITE (out,80) PhyMR(ng), 'PhyMR',                          &
     &            'Phytoplankton mortality rate (day-1).'
            WRITE (out,90) ZooAE_N(ng), 'ZooAE_N',                      &
     &            'Zooplankton nitrogen assimilation efficiency',       &
     &            '(nondimensional).'
            WRITE (out,80) ZooBM(ng), 'ZooBM',                          &
     &            'Rate for zooplankton basal metabolism (1/day).'
            WRITE (out,80) ZooCN(ng), 'ZooCN',                          &
     &            'Zooplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,80) ZooER(ng), 'ZooER',                          &
     &            'Zooplankton specific excretion rate (day-1).'
            WRITE (out,80) ZooGR(ng), 'ZooGR',                          &
     &            'Zooplankton maximum growth rate (day-1).'
            WRITE (out,80) ZooMin(ng), 'ZooMin',                        &
     &            'Zooplankton minimum threshold (mmol_N/m3).'
            WRITE (out,80) ZooMR(ng), 'ZooMR',                          &
     &            'Zooplankton mortality rate (day-1).'
            WRITE (out,80) LDeRRN(ng), 'LDeRRN',                        &
     &            'Large detritus N re-mineralization rate (day-1).'
            WRITE (out,80) LDeRRC(ng), 'LDeRRC',                        &
     &            'Large detritus C re-mineralization rate (day-1).'
            WRITE (out,80) CoagR(ng), 'CoagR',                          &
     &            'Coagulation rate (day-1).'
            WRITE (out,80) SDeRRN(ng), 'SDeRRN',                        &
     &            'Remineralization rate for small detritus N (day-1).'
            WRITE (out,80) SDeRRC(ng), 'SDeRRC',                        &
     &            'Remineralization rate for small detritus C (day-1).'
            WRITE (out,80) wPhy(ng), 'wPhy',                            &
     &            'Phytoplankton sinking velocity (m/day).'
            WRITE (out,80) wLDet(ng), 'wLDet',                          &
     &            'Large detritus sinking velocity (m/day).'
            WRITE (out,80) wSDet(ng), 'wSDet',                          &
     &            'Small detritus sinking velocity (m/day).'
            WRITE (out,80) pCO2air(ng), 'pCO2air',                      &
     &            'CO2 partial pressure in air (ppm by volume).'
#ifdef SPECTRAL_LIGHT
            WRITE (out,80) SIGATRB(ng), 'SIGATRB',                      &
     &            'NAP absorb. cross section at 440 nm (m2/g).'
            WRITE (out,80) STRB(ng), 'STRB',                            &
     &            'Spectral slope of NAP absorption.'
            WRITE (out,80) BLTRB(ng), 'BLTRB',                          &
     &            'Baseline NAP absorption (m2/g).'
            WRITE (out,80) SIGBTRB(ng), 'SIGBTRB',                      &
     &            'Scattering cross section of turbidity (m2/g).'
            WRITE (out,80) ETASPEC(ng), 'ETASPEC',                      &
     &            'Scattering spectral exponent.'
            WRITE (out,80) BB2B(ng), 'BB2B',                            &
     &            'Particulate backscattering ratio.'
#endif
#ifdef SAV_BIOMASS
            WRITE (out,70) GMODopt(ng), 'GMODopt',                      &
     &            'SAV growth rate options.'
            WRITE (out,80) KNSED(ng), 'KNSED',                          &
     &            'Half saturation coeff.-sediment N uptake (mu-M).'
            WRITE (out,80) KNWC(ng), 'KNWC',                            &
     &            'Half saturation coeff.-water-column N uptake (mu-M).'
            WRITE (out,80) TOPT(ng), 'TOPT',                            &
     &            'Optimum temperature for SAV growth (deg C).'
            WRITE (out,80) THTA(ng), 'THTA',                            &
     &            'Temp. growth theta for growth Model #1 (unitless).'
            WRITE (out,80) THTA2(ng), 'THTA2',                          &
     &            'Temp. growth theta for growth Model #2 (unitless).'
            WRITE (out,80) SCL(ng), 'SCL',                              &
     &            'SAV growth fraction for growth Model #1 (unitless).'
            WRITE (out,80) SCL2(ng), 'SCL2',                            &
     &            'SAV growth fraction for growth Model #2 (unitless).'
            WRITE (out,80) KI(ng), 'KI',                                &
     &            'Half saturation for light (muE m-2 s-1).'
            WRITE (out,80) SR(ng), 'SR',                                &
     &            'Surface reflectance of light.'
            WRITE (out,80) LMBAMX(ng), 'LMBAMX',                        &
     &            'Maximum AGB in self-shading formulation (g C m-2). '
            WRITE (out,80) KMAG(ng), 'KMAG',                            &
     &            'Above biomass ground mortality rate (1/d).'
            WRITE (out,80) ARSC(ng), 'ARSC',                            &
     &            'Maximum fraction of PP respired.'
            WRITE (out,80) ARC(ng), 'ARC',                              &
     &            'Active respiration coefficient.'
            WRITE (out,80) BSRC(ng), 'BSRC',                            &
     &            'Maximum fraction of biomass respired.'
            WRITE (out,80) RC(ng), 'RC',                                &
     &            'Basal respiration coefficient (1/d).'
            WRITE (out,80) RtStTL(ng), 'RtStTL',                        &
     &            'Seasonal root storage coefficient (1/d).'
            WRITE (out,80) DOWNt(ng), 'DOWNt',                          &
     &            'Downward translocation coefficient (1/d).'
            WRITE (out,80) TRNS(ng), 'TRNS',                            &
     &            'Upward translocation coefficient (1/d).'
            WRITE (out,80) TCRIT(ng), 'TCRIT',                          &
     &            'Crit. temp. development of abovegrd. biomass (deg C).'
            WRITE (out,80) KM(ng), 'KM',                                &
     &            'Below ground biomass mortality (1/d).'
#endif 
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) nl_tnu2(i,ng), 'nl_tnu2', i,              &
     &              'NLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,100) ad_tnu2(i,ng), 'ad_tnu2', i,              &
     &              'ADM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,100) tl_tnu2(i,ng), 'tl_tnu2', i,              &
     &              'TLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) nl_tnu4(i,ng), 'nl_tnu4', i,              &
     &              'NLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,100) ad_tnu4(i,ng), 'ad_tnu4', i,              &
     &              'ADM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,100) tl_tnu4(i,ng), 'tl_tnu4', i,              &
     &              'TLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,110) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,100) Akt_bak(i,ng), 'Akt_bak', i,               &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef FORWARD_MIXING
            DO itrc=1,NBT
              i=idbio(itrc)
# ifdef ADJOINT
              WRITE (out,100) ad_Akt_fac(i,ng), 'ad_Akt_fac', i,        &
     &              'ADM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,100) tl_Akt_fac(i,ng), 'tl_Akt_fac', i,        &
     &              'TLM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) Tnudg(i,ng), 'Tnudg', i,                  &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,110) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning ON  point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning OFF point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,110) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,110) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,110) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,120)                   &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,120)                   &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            WRITE (out,'(1x)')
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTvar(i),ng)) WRITE (out,120)                   &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                   &
     &            'Write out averaged tracer ', i,                      &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTTav(i),ng)) WRITE (out,120)                   &
     &            Aout(idTTav(i),ng), 'Aout(idTTav)',                   &
     &            'Write out averaged <t*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idUTav(i),ng)) WRITE (out,120)                   &
     &            Aout(idUTav(i),ng), 'Aout(idUTav)',                   &
     &            'Write out averaged <u*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idVTav(i),ng)) WRITE (out,120)                   &
     &            Aout(idVTav(i),ng), 'Aout(idVTav)',                   &
     &            'Write out averaged <v*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHUTav(i),ng)) WRITE (out,120)                   &
     &            Aout(iHUTav(i),ng), 'Aout(iHUTav)',                   &
     &            'Write out averaged <Huon*t> for tracer ', i,         &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHVTav(i),ng)) WRITE (out,120)                   &
     &            Aout(iHVTav(i),ng), 'Aout(iHVTav)',                   &
     &            'Write out averaged <Hvom*t> for tracer ', i,         &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
#endif
#ifdef DIAGNOSTICS_TS
            WRITE (out,'(1x)')
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTrate),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTrate)',                 &
     &              'Write out rate of change of tracer ', itrc,        &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThadv),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iThadv)',                 &
     &              'Write out horizontal advection, tracer ', itrc,    &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTxadv),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTxadv)',                 &
     &              'Write out horizontal X-advection, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTyadv),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTyadv)',                 &
     &              'Write out horizontal Y-advection, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvadv),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTvadv)',                 &
     &              'Write out vertical advection, tracer ', itrc,      &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
# if defined TS_DIF2 || defined TS_DIF4
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThdif),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iThdif)',                 &
     &              'Write out horizontal diffusion, tracer ', itrc,    &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(i,iTxdif),ng))                            &
     &          WRITE (out,120) .TRUE., 'Dout(iTxdif)',                 &
     &              'Write out horizontal X-diffusion, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTydif),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTydif)',                 &
     &              'Write out horizontal Y-diffusion, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTsdif),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTsdif)',                 &
     &              'Write out horizontal S-diffusion, tracer ', itrc,  &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
#  endif
# endif
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvdif),ng))                         &
     &          WRITE (out,120) .TRUE., 'Dout(iTvdif)',                 &
     &              'Write out vertical diffusion, tracer ', itrc,      &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
#endif
#ifdef DIAGNOSTICS_BIO
            IF (NDbio2d.gt.0) THEN
              DO itrc=1,NDbio2d
                i=iDbio2(itrc)
                IF (Dout(i,ng)) WRITE (out,130)                         &
     &              Dout(i,ng), 'Hout(iDbio2)',                         &
     &              'Write out diagnostics for', TRIM(Vname(1,i))
              END DO
            END IF
            DO itrc=1,NDbio3d
              i=iDbio3(itrc)
              IF (Dout(i,ng)) WRITE (out,130)                           &
     &            Dout(i,ng), 'Hout(iDbio3)',                           &
     &            'Write out diagnostics for', TRIM(Vname(1,i))
            END DO
#endif
#ifdef SPECTRAL_LIGHT
          IF (Hout(idPARo,ng)) WRITE (out,120) Hout(idPARo,ng),         &
     &       'Hout(idPARo)',                                            &
     &       'Write out Photosynthetical Available Radiation PAR.'
          IF (Hout(idPARs,ng)) WRITE (out,120) Hout(idPARs,ng),         &
     &       'Hout(idPARs)',                                            &
     &       'Write out Spectral PAR.'
          IF (Hout(idSpKd,ng)) WRITE (out,120) Hout(idSpKd,ng),         &
     &       'Hout(idSpKd)',                                            &
     &       'Write out Spectral Kd (light attenuation).'
#endif
#ifdef CDOM_VARIABLE
            WRITE (out,70) Ndom, 'Ndom',                      &
     &            'Number of DOM constituents.'
#endif
#ifdef SAV_BIOMASS 
          IF (Hout(iddins,ng)) WRITE (out,120) Hout(iddins,ng),         &
     &       'Hout(iddins)',                                            &
     &       'Write out sediment dissolved inorganic N in units (mu-M).'
          IF (Hout(iddinw,ng)) WRITE (out,120) Hout(iddinw,ng),         &
     &       'Hout(iddinw)',                                            &
     &       'Write out water column inorganic N in units (mu-M).'
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters.
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          nl_tnu4(i,ng)=SQRT(ABS(nl_tnu4(i,ng)))
#ifdef ADJOINT
          ad_tnu4(i,ng)=SQRT(ABS(ad_tnu4(i,ng)))
#endif
#if defined TANGENT || defined TL_IOMS
          tl_tnu4(i,ng)=SQRT(ABS(tl_tnu4(i,ng)))
#endif
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)
  40  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
  50  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Fennel Model Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  70  FORMAT (1x,i10,2x,a,t32,a)
  80  FORMAT (1p,e11.4,2x,a,t32,a)
  90  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
 100  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 120  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 130  FORMAT (10x,l1,2x,a,t32,a,1x,a)

      RETURN
      END SUBROUTINE read_BioPar
