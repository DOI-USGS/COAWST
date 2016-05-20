      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in EcoSim bio-optical model input parameters.    !
!  They are specified in input script "ecosim.in".                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_eclight
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
      integer :: i, ifield, igrid, is, itracer, itrc, ng, nline, status
      integer :: ibac, iband, ifec, iphy

      integer :: decode_line, load_i, load_l, load_lbc, load_r

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
!  Read in EcoSim bio-optical model parameters.
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
            CASE ('RtUVR_flag')
              Npts=load_l(Nval, Cval, Ngrids, RtUVR_flag)
            CASE ('NFIX_flag')
              Npts=load_l(Nval, Cval, Ngrids, NFIX_flag)
            CASE ('Regen_flag')
              Npts=load_l(Nval, Cval, Ngrids, Regen_flag)
            CASE ('HsNO3')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsNO3)
            CASE ('HsNH4')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsNH4)
            CASE ('HsSiO')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsSiO)
            CASE ('HsPO4')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsPO4)
            CASE ('HsFe')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsFe)
            CASE ('GtALG_max')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, GtALG_max)
            CASE ('PhyTbase')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, PhyTbase)
            CASE ('PhyTfac')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, PhyTfac)
            CASE ('BET_')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, BET_)
            CASE ('maxC2nALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2nALG)
            CASE ('minC2nALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2nALG)
            CASE ('C2nALGminABS')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, C2nALGminABS)
            CASE ('maxC2SiALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2SiALG)
            CASE ('minC2SiALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2SiALG)
            CASE ('C2SiALGminABS')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, C2SiALGminABS)
            CASE ('maxC2pALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2pALG)
            CASE ('minC2pALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2pALG)
            CASE ('C2pALGminABS')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, C2pALGminABS)
            CASE ('maxC2FeALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2FeALG)
            CASE ('minC2FeALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2FeALG)
            CASE ('C2FeALGminABS')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, C2FeALGminABS)
            CASE ('qu_yld')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, qu_yld)
            CASE ('E0_comp')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, E0_comp)
            CASE ('E0_inhib')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, E0_inhib)
            CASE ('inhib_fac')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, inhib_fac)
            CASE ('C2CHL_max')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, C2CHL_max)
            CASE ('mxC2Cl')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxC2Cl)
            CASE ('b_C2Cl')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_C2Cl)
            CASE ('mxC2Cn')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxC2Cn)
            CASE ('b_C2Cn')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_C2Cn)
            CASE ('mxPacEff')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPacEff)
            CASE ('b_PacEff')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PacEff)
            CASE ('mxChlB')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxChlB)
            CASE ('b_ChlB')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_ChlB)
            CASE ('mxChlC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxChlC)
            CASE ('b_ChlC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_ChlC)
            CASE ('mxPSC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPSC)
            CASE ('b_PSC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PSC)
            CASE ('mxPPC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPPC)
            CASE ('b_PPC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PPC)
            CASE ('mxLPUb')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxLPUb)
            CASE ('b_LPUb')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_LPUb)
            CASE ('mxHPUb')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, mxHPUb)
            CASE ('b_HPUb')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, b_HPUb)
            CASE ('FecDOC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, FecDOC)
            CASE ('FecPEL')
              Npts=load_r(Nval, Rval, Nphy*Nfec*Ngrids, FecPEL)
            CASE ('FecCYC')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, FecCYC)
            CASE ('ExALG')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, ExALG)
            CASE ('WS')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, WS)
            CASE ('HsGRZ')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsGRZ)
            CASE ('MinRefuge')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, MinRefuge)
            CASE ('RefugeDep')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, RefugeDep)
            CASE ('Norm_Vol')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, Norm_Vol)
            CASE ('Norm_Surf')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, Norm_Surf)
            CASE ('HsDOP')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsDOP)
            CASE ('C2pALKPHOS')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, C2pALKPHOS)
            CASE ('HsDON')
              Npts=load_r(Nval, Rval, Nphy*Ngrids, HsDON)
            CASE ('C2nNupDON')
              Npts=load_r(Nval, Rval, Ngrids, C2nNupDON)
            CASE ('C2nBAC')
              Npts=load_r(Nval, Rval, Ngrids, C2nBAC)
            CASE ('C2pBAC')
              Npts=load_r(Nval, Rval, Ngrids, C2pBAC)
            CASE ('C2FeBAC')
              Npts=load_r(Nval, Rval, Ngrids, C2FeBAC)
            CASE ('HsDOC_ba')
              Npts=load_r(Nval, Rval, Nbac*Ngrids, HsDOC_ba)
            CASE ('GtBAC_max')
              Npts=load_r(Nval, Rval, Nbac*Ngrids, GtBAC_max)
            CASE ('BacTbase')
              Npts=load_r(Nval, Rval, Nbac*Ngrids, BacTbase)
            CASE ('BacTfac')
              Npts=load_r(Nval, Rval, Nbac*Ngrids, BacTfac)
            CASE ('BacDOC')
              Npts=load_r(Nval, Rval, Ngrids, BacDOC)
            CASE ('BacPEL')
              Npts=load_r(Nval, Rval, Ngrids, BacPEL)
            CASE ('BacCYC')
              Npts=load_r(Nval, Rval, Ngrids, BacCYC)
            CASE ('ExBAC_c')
              Npts=load_r(Nval, Rval, Ngrids, ExBAC_c)
            CASE ('ExBacC2N')
              Npts=load_r(Nval, Rval, Ngrids, ExBacC2N)
            CASE ('Bac_Ceff')
              Npts=load_r(Nval, Rval, Ngrids, Bac_Ceff)
            CASE ('RtNIT')
              Npts=load_r(Nval, Rval, Ngrids, RtNIT)
            CASE ('HsNIT')
              Npts=load_r(Nval, Rval, Ngrids, HsNIT)
            CASE ('cDOCfrac_c')
              Npts=load_r(Nval, Rval, Ndom*Ngrids, cDOCfrac_c)
            CASE ('RtUVR_DIC')
              Npts=load_r(Nval, Rval, Ngrids, RtUVR_DIC)
            CASE ('RtUVR_DOC')
              Npts=load_r(Nval, Rval, Ngrids, RtUVR_DOC)
            CASE ('WF')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, WF)
            CASE ('RegTbase')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegTbase)
            CASE ('RegTfac')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegTfac)
            CASE ('RegCR')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegCR)
            CASE ('RegNR')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegNR)
            CASE ('RegSR')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegSR)
            CASE ('RegPR')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegPR)
            CASE ('RegFR')
              Npts=load_r(Nval, Rval, Nfec*Ngrids, RegFR)
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
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Initialize secondary parameters.
!-----------------------------------------------------------------------
!
!  Convert rates from day-1 to second-1.
!
      DO ng=1,Ngrids
        DO iphy=1,Nphy
          GtALG_max(iphy,ng)=GtALG_max(iphy,ng)*sec2day
          ExALG(iphy,ng)=ExALG(iphy,ng)*sec2day
          HsGRZ(iphy,ng)=HsGRZ(iphy,ng)*sec2day
          WS(iphy,ng)=WS(iphy,ng)*sec2day
        END DO
        DO ibac=1,Nbac
          GtBAC_max(ibac,ng)=GtBAC_max(ibac,ng)*sec2day
        END DO
        DO ifec=1,Nfec
          WF(ifec,ng)=WF(ifec,ng)*sec2day
        END DO
        RtNIT(ng)=RtNIT(ng)*sec2day
      END DO
!
!  Calculated reciprocal phytoplankton parameters.
!
      DO ng=1,Ngrids
        DO iphy=1,Nphy
          IF (maxC2nALG(iphy,ng).gt.SMALL) THEN
            ImaxC2nALG(iphy,ng)=1.0_r8/maxC2nALG(iphy,ng)
          ELSE
            ImaxC2nALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2SiALG(iphy,ng).gt.SMALL) THEN
            ImaxC2SiALG(iphy,ng)=1.0_r8/maxC2SiALG(iphy,ng)
          ELSE
            ImaxC2SiALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2pALG(iphy,ng).gt.SMALL) THEN
            ImaxC2pALG(iphy,ng)=1.0_r8/maxC2pALG(iphy,ng)
          ELSE
            ImaxC2pALG(iphy,ng)=0.0_r8
          END IF
          IF (maxC2FeALG(iphy,ng).gt.SMALL) THEN
            ImaxC2FeALG(iphy,ng)=1.0_r8/maxC2FeALG(iphy,ng)
          ELSE
            ImaxC2FeALG(iphy,ng)=0.0_r8
          END IF
        END DO
      END DO
!
!  Calculated bacterial parameters.
!
      DO ng=1,Ngrids
        DO ibac=1,Nbac
          HsNH4_ba(ibac,ng)=HsDOC_ba(ibac,ng)/C2nBAC(ng)
          HsPO4_ba(ibac,ng)=HsDOC_ba(ibac,ng)/C2pBAC(ng)
          HsFe_ba (ibac,ng)=HsDOC_ba(ibac,ng)/C2FeBAC(ng)
        END DO
      END DO
!
!  Inverse parameters for computational efficiency.
!
      DO ng=1,Ngrids
        N2cBAC(ng)=1.0_r8/C2nBAC(ng)
        P2cBAC(ng)=1.0_r8/C2pBAC(ng)
        Fe2cBAC(ng)=1.0_r8/C2FeBAC(ng)
        I_Bac_Ceff(ng)=1.0_r8/Bac_Ceff(ng)
      END DO
!
!  Reciprocal of non baterial recalcitran carbon excretion.
!
      DO ng=1,Ngrids
        R_ExBAC_c(ng)=1.0_r8/(1.0_r8-ExBAC_c(ng))
      END DO
!
!  Bacterial recalcitrant nitrogen excretion as a function of uptake.
!
      DO ng=1,Ngrids
        ExBAC_n(ng)=ExBAC_c(ng)*C2nBAC(ng)/ExBacC2N(ng)
        Frac_ExBAC_n(ng)=1.0_r8-ExBAC_n(ng)
      END DO
!
!  Scale UV degradation parameters.
!
      DO ng=1,Ngrids
        RtUVR_DIC(ng)=RtUVR_DIC(ng)/3600.0_r8
        RtUVR_DOC(ng)=RtUVR_DOC(ng)/3600.0_r8
      END DO
!
!  If applicable, zero-out fecal regeneration parameters.
!
      DO ng=1,Ngrids
        IF (Regen_flag(ng)) THEN
          DO ifec=1,Nfec
            RegCR(ifec,ng)=RegCR(ifec,ng)*sec2day
            RegNR(ifec,ng)=RegNR(ifec,ng)*sec2day
            RegPR(ifec,ng)=RegPR(ifec,ng)*sec2day
            RegFR(ifec,ng)=RegFR(ifec,ng)*sec2day
            RegSR(ifec,ng)=RegSR(ifec,ng)*sec2day
          END DO
        ELSE
          DO ifec=1,Nfec
            RegCR(ifec,ng)=0.0_r8
            RegNR(ifec,ng)=0.0_r8
            RegPR(ifec,ng)=0.0_r8
            RegFR(ifec,ng)=0.0_r8
            RegSR(ifec,ng)=0.0_r8
          END DO
        END IF
      END DO
!
!  Spectral dependency for scattering and backscattering.
!
      DO iband=1,NBands
        wavedp(iband)=(550.0_r8/(397.0_r8+REAL(iband,r8)*DLAM))
      END DO
!
!  Calculated IOP parameter values.
!
      aDOC410(ilab)=aDOC(ilab,1)*EXP(0.014_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC410(irct)=aDOC(irct,1)*EXP(0.025_r8*(ec_wave_ab(1)-410.0_r8))
      aDOC300(ilab)=EXP(0.0145_r8*(410.0_r8-300.0_r8))
      aDOC300(irct)=EXP(0.0145_r8*(410.0_r8-300.0_r8))
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,70) RtUVR_flag(ng), 'RtUVR_flag',                &
     &            'Switch to calculate CDOC UV photolysis.'
            WRITE (out,70) NFIX_flag(ng), 'NFIX_flag',                  &
     &            'Switch to calculate temperature based N fixation.'
            WRITE (out,70) Regen_flag(ng), 'Regen_flag',                &
     &            'Switch to calculate fecal matter regeneration.'
            WRITE (out,80) 'HsNO3',                                     &
     &            'Half-saturation for phytoplankton NO3 uptake',       &
     &            '(micromole_NO3/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsNO3(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsNH4',                                     &
     &            'Half-saturation for phytoplankton NH4 uptake',       &
     &            '(micromole_NH4/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsNH4(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsSiO',                                     &
     &            'Half-saturation for phytoplankton SiO uptake',       &
     &            '(micromole_SiO/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsSiO(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsPO4',                                     &
     &            'Half-saturation for phytoplankton PO4 uptake',       &
     &            '(micromole_PO4/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsPO4(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsFe',                                      &
     &            'Half-saturation for phytoplankton Fe uptake',        &
     &            '(micromole_Fe/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsFe(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'GtALG_max',                                &
     &            'Maximum 24 hour growth rate (1/day):'
              DO is=1,Nphy
                WRITE (out,90) GtALG_max(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'PhyTbase',                                  &
     &            'Temperature base for exponential response to',       &
     &            'temperature (Celsius):'
              DO is=1,Nphy
                WRITE (out,90) PhyTbase(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'PhyTfac',                                   &
     &            'Phytoplankton exponential temperature factor',       &
     &            '(1/Celsius):'
              DO is=1,Nphy
                WRITE (out,90) PhyTfac(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'BET_',                                     &
     &            'Nitrate uptake inhibition for NH4 (l/micromole):'
              DO is=1,Nphy
                WRITE (out,90) BET_(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'maxC2nALG',                                 &
     &            'Maximum phytoplankton C:N ratio',                    &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) maxC2nALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'minC2nALG',                                 &
     &            'Balanced phytoplankton C:N ratio',                   &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) minC2nALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2nALGminABS',                              &
     &            'Absolute minimum phytoplankton C:N ratio',           &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) C2nALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'maxC2SiALG',                                &
     &            'Maximum phytoplankton C:Si ratio',                   &
     &            '(micromole_C/micromole_Si)'
              DO is=1,Nphy
                WRITE (out,90) maxC2SiALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'minC2SiALG',                                &
     &            'Balanced phytoplankton C:Si ratio',                  &
     &            '(micromole_C/micromole_Si):'
              DO is=1,Nphy
                WRITE (out,90) minC2SiALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2SiALGminABS',                             &
     &            'Absolute minimum phytoplankton C:Si ratio',          &
     &            '(micromole_C/micromole_Si):'
              DO is=1,Nphy
                WRITE (out,90) C2SiALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'maxC2pALG',                                 &
     &            'Maximum phytoplankton C:P ratio',                    &
     &            '(micromole_C/micromole_P):'
              DO is=1,Nphy
                WRITE (out,90) maxC2pALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'minC2pALG',                                 &
     &            'Balanced phytoplankton C:P ratio',                   &
     &            '(micromole_C/micromole_P):'
              DO is=1,Nphy
                WRITE (out,90) minC2pALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2pALGminABS',                              &
     &            'Absolute minimum phytoplankton C:P ratio',           &
     &            '(micromole_C/micromole_P)'
              DO is=1,Nphy
                WRITE (out,90) C2pALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'maxC2FeALG',                                &
     &            'Maximum phytoplankton C:Fe ratio',                   &
     &            '(micromole_C/micromole_Fe):'
              DO is=1,Nphy
                WRITE (out,90) maxC2FeALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'minC2FeALG',                                &
     &            'Balanced phytoplankton C:Fe ratio',                  &
     &            '(micromole_C/micromole_Fe):'
              DO is=1,Nphy
                WRITE (out,90) minC2FeALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2FeALGminABS',                             &
     &            'Absolute minimum phytoplankton C:Fe ratio',          &
     &            '(micromole_C/micromole_Fe):'
              DO is=1,Nphy
                WRITE (out,90) C2FeALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'qu_yld',                                    &
     &            'Maximum quantum yield',                              &
     &            '(micromole_C/micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) qu_yld(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'E0_comp',                                  &
     &            'Compensation light level (micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) E0_comp(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'E0_inhib',                                  &
     &            'Light level for onset of photoinhibition',           &
     &            '(micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) E0_inhib(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'inhib_fac',                                 &
     &            'Exponential decay factor for light limited growth',  &
     &            '(1/micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) inhib_fac(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2CHL_max',                                 &
     &            'Maximum lighted limited C:Chl ratio',                &
     &            '(microgram_C/microgram_Chl):'
              DO is=1,Nphy
                WRITE (out,90) C2CHL_max(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxC2Cl',                                    &
     &            'Rate of change in light limited C:Chl ratio',        &
     &            '(microgram_C/microgram_Chl/micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) mxC2Cl(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_C2Cl',                                    &
     &            'Minimum lighted limited C:Chl ratio',                &
     &            '(microgram_C/microgram_Chl):'
              DO is=1,Nphy
                WRITE (out,90) b_C2Cl(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxC2Cn',                                    &
     &            'Rate of change in nutient limited C:Chl ratio',      &
     &            '[(ug_C/ug_Chl)/(umole_C/umole_N)]:'
              DO is=1,Nphy
                WRITE (out,90) mxC2Cn(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_C2Cn',                                    &
     &            'Minimum nutrient limited C:Chl ratio',               &
     &            '(microgram_C/microgram_Chl):'
              DO is=1,Nphy
                WRITE (out,90) b_C2Cn(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxPacEff',                                  &
     &            'Rate of change in package effect',                   &
     &            '[1/(microgram_C/microgram_Chl)]:'
              DO is=1,Nphy
                WRITE (out,90) mxPacEff(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_PacEff',                                  &
     &            'Maximum package effect',                             &
     &            '[1/(microgram_C/microgram_Chl)]:'
              DO is=1,Nphy
                WRITE (out,90) b_PacEff(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxChlB',                                    &
     &            'Rate of change in the Chl_b:Chl_a ratio',            &
     &            '[(ug_Chl_b/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxChlB(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_ChlB',                                    &
     &            'Maximum Chl_b:Chl_a ratio',                          &
     &            '(microgram_Chl_b/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_ChlB(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxChlC',                                    &
     &            'Rate of change in the Chl_c:Chl_a ratio',            &
     &            '[(ug_Chl_c/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxChlC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_ChlC',                                    &
     &            'Maximum Chl_c:Chl_a ratio',                          &
     &            '(microgram_Chl_c/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_ChlC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxPSC',                                     &
     &            'Rate of change in the PSC:Chl_a ratio',              &
     &            '[(ug_PSC/ug_Chl_a)/ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxPSC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_PSC',                                     &
     &            'Maximum PSC:Chl_a ratio',                            &
     &            '(microgram_PSC/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_PSC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxPPC',                                     &
     &            'Rate of change in the PPC:Chl_a ratio',              &
     &            '[(ug_PPC/ug_Chl_a)/(ug_C/ug_Chl_ a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxPPC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_PPC',                                     &
     &            'Maximum PPC:Chl_a ratio',                            &
     &            '(microgram_Chl_c/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_PPC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxLPUb',                                    &
     &            'Rate of change in the LPUb:Chl_a ratio',             &
     &            '[(ug_LPUb/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxLPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'b_LPUb',                                    &
     &            'Maximum LPUb:Chl_a ratio',                           &
     &            '(migrogram_HPUb/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_LPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'mxHPUb',                                    &
     &            'Rate of change in the HPUb:Chl_a ratio',             &
     &            '[(ug_HPUb/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxHPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80)'b_HPUb',                                     &
     &            'Maximum HPUb:Chl_a ratio',                           &
     &            '(microgram_HPUb/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_HPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'FecDOC',                                    &
     &            'Proportion of grazing stress apportioned to DOM',    &
     &            '(nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) FecDOC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'FecPEL',                                    &
     &            'Proportion of grazing stress apportioned to fecal',  &
     &            '(nondimensional):'
              DO i=1,Nfec
                DO is=1,Nphy
                  WRITE (out,110) FecPEL(is,i,ng), i, TRIM(PhyName(is))
                END DO
              END DO
            WRITE (out,80) 'FecCYC',                                    &
     &            'Proportion of grazing stress that is recycled',      &
     &            '(nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) FecCYC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'ExALG',                                     &
     &            'Proportion of daily production lost to excretion',   &
     &            '(nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) ExALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'WS',                                       &
     &            'Phytoplankton sinking speed (meters/day):'
              DO is=1,Nphy
                WRITE (out,90) WS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'HsGRZ',                                    &
     &            'Phytoplankton grazing parameter (nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) HsGRZ(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'MinRefuge',                                &
     &            'Refuge Phytoplankton population (micromole_C/liter):'
              DO is=1,Nphy
                WRITE (out,90) MinRefuge(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'RefugeDep',                                &
     &            'Maximum Refuge Phytoplankton depth (meters):'
              DO is=1,Nphy
                WRITE (out,90) RefugeDep(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'Norm_Vol',                                 &
     &            'Normalized Volume factor (nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) Norm_Vol(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,100) 'Norm_Surf',                                &
     &            'Normalized Surface Area factor (nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) Norm_Surf(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsDOP',                                     &
     &            'Half Saturation Constant for DOP uptake',            &
     &            '(micromole_DOP/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsDOP(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2pALKPHOS',                                &
     &            'C:P ratio where DOP uptake begins',                  &
     &            '(micromole_C/micromole_P):'
              DO is=1,Nphy
                WRITE (out,90) C2pALKPHOS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsDON',                                     &
     &            'Half Saturation Constant for DON uptake',            &
     &            '(micromole_DON/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsDON(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'C2nNupDON',                                 &
     &            'C:N ratio where DON uptake begins',                  &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) C2nNupDON(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsDOC_ba',                                  &
     &            'Half saturation constant for bacteria DOC uptake',   &
     &            '(micromole_DOC/liter):'
              DO is=1,Nbac
                WRITE (out,90) HsDOC_ba(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,100) 'GtBAC_max',                                &
     &            'Maximum 24 hour bacterial growth rate (1/day):'
              DO is=1,Nbac
                WRITE (out,90) GtBAC_max(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,80) 'BacTbase',                                  &
     &            'Temperature base for exponetial response to',        &
     &            'temperature (Celsius):'
              DO is=1,Nbac
                WRITE (out,90) BacTbase(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,80) 'BacTfac',                                   &
     &            'Bacteria exponential temperature factor',            &
     &            '(1/Celsius):'
              DO is=1,Nbac
                WRITE (out,90) BacTfac(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,120) C2nBAC(ng), 'C2nBAC',                       &
     &            'Carbon to Nitrogen ratio of Bacteria',               &
     &            '(micromole_C/micromole_N).'
            WRITE (out,120) C2pBAC(ng), 'C2pBAC',                       &
     &            'Carbon to Phosphorus ratio of Bacteria',             &
     &            '(micromole_C/micromole_P).'
            WRITE (out,120) C2FeBAC(ng), 'C2FeBAC',                     &
     &            'Carbon to Iron ratio of Bacteria',                   &
     &            '(micromole_C/micromole_Fe).'
            WRITE (out,120) BacDOC(ng), 'BacDOC',                       &
     &            'Proportion of bacteria grazing stress apportioned',  &
     &            'to DOM (nondimensional).'
            WRITE (out,120) BacPEL(ng), 'BacPEL',                       &
     &            'Proportion of bacteria grazing stress apportioned',  &
     &            'to fecal (nondimensional).'
            WRITE (out,120) BacCYC(ng), 'BacCYC',                       &
     &            'Proportion of bacteria grazing stress recycled',     &
     &            '(nondimensional).'
            WRITE (out,120) ExBAC_c(ng), 'ExBAC_c',                     &
     &            'Bacterial recalcitrant C excretion as proportion',   &
     &            'of uptake (nondimensional).'
            WRITE (out,120) ExBacC2N(ng), 'ExBacC2N',                   &
     &            'Bacterial recalcitrant excretion carbon:nitrogen',   &
     &            'ratio (micromole_C/micromole_N).'
            WRITE (out,120) Bac_Ceff(ng), 'Bac_Ceff',                   &
     &            'Bacterial gross growth carbon efficiency',           &
     &            '(nondimensional).'
            WRITE (out,130) RtNIT(ng), 'RtNIT',                         &
     &            'Maximum nitrification rate (1/day).'
            WRITE (out,120) HsNIT(ng), 'HsNIT',                         &
     &            'Half saturation constant for bacteria nitrification',&
     &            '(micromole_NH4/liter).'
            WRITE (out,80) 'cDOCfrac_c',                                &
     &            'Colored fraction of DOC from phytoplakton and',      &
     &            'bacterial losses (nondimensional):'
              DO is=1,Ndom
                WRITE (out,90) cDOCfrac_c(is,ng), TRIM(DomName(is))
              END DO
            WRITE (out,120) RtUVR_DIC(ng), 'RtUVR_DIC',                 &
     &            'UV degradation of DOC into DIC at 410 nm',           &
     &            '(micromole/meter/liter/hour).'
            WRITE (out,120) RtUVR_DOC(ng), 'RtUVR_DOC',                 &
     &            'UV degradation of DOC to colorless labile DOC',      &
     &            'at 410 nm (micromole/meter/liter/hour).'
            WRITE (out,100) 'WF',                                       &
     &            'Fecal sinking flux (meters/day):'
              DO is=1,Nfec
                WRITE (out,90) WF(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegTbase',                                  &
     &            'Fecal regeneration temperature base for exponential',&
     &            'response to temperature (Celsius):'
              DO is=1,Nfec
                WRITE (out,90) RegTbase(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegTfac',                                   &
     &            'Fecal regeneration exponential temperature factor',  &
     &            '(1/Celsius):'
              DO is=1,Nfec
                WRITE (out,90) RegTfac(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,100) 'RegCR',                                    &
     &            'Fecal carbon regeneration rate (1/day):'
              DO is=1,Nfec
                WRITE (out,90) RegCR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,100) 'RegNR',                                    &
     &            'Fecal nitrogen regeneration rate (1/day)'
              DO is=1,Nfec
                WRITE (out,90) RegNR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,100) 'RegSR',                                    &
     &            'Fecal silica regeneration rate (1/day):'
              DO is=1,Nfec
                WRITE (out,90) RegSR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,100) 'RegPR',                                    &
     &            'Fecal phosphorus regeneration rate (1/day):'
              DO is=1,Nfec
                WRITE (out,90) RegPR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,100) 'RegFR',                                    &
     &            'Fecal iron regeneration rate (1/day)'
              DO is=1,Nfec
                WRITE (out,90) RegFR(is,ng), TRIM(FecName(is))
              END DO
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,140) nl_tnu2(i,ng), 'nl_tnu2', i,              &
     &              'NLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,140) ad_tnu2(i,ng), 'ad_tnu2', i,              &
     &              'ADM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,140) tl_tnu2(i,ng), 'tl_tnu2', i,              &
     &              'TLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,140) nl_tnu4(i,ng), 'nl_tnu4', i,              &
     &              'NLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,140) ad_tnu4(i,ng), 'ad_tnu4', i,              &
     &              'ADM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,140) tl_tnu4(i,ng), 'tl_tnu4', i,              &
     &              'TLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,150) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,140) Akt_bak(i,ng), 'Akt_bak', i,              &
     &              'Background vertical mixing coefficient (m2/s)',    &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef FORWARD_MIXING
            DO itrc=1,NBT
              i=idbio(itrc)
# ifdef ADJOINT
              WRITE (out,140) ad_Akt_fac(i,ng), 'ad_Akt_fac', i,        &
     &              'ADM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,140) tl_Akt_fac(i,ng), 'tl_Akt_fac', i,        &
     &              'TLM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,140) Tnudg(i,ng), 'Tnudg', i,                  &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning ON  point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning OFF point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,150) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,160)                   &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,160)                   &
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
              IF (Aout(idTvar(i),ng)) WRITE (out,160)                   &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                   &
     &            'Write out averaged tracer ', i,                      &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTTav(i),ng)) WRITE (out,160)                   &
     &            Aout(idTTav(i),ng), 'Aout(idTTav)',                   &
     &            'Write out averaged <t*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idUTav(i),ng)) WRITE (out,160)                   &
     &            Aout(idUTav(i),ng), 'Aout(idUTav)',                   &
     &            'Write out averaged <u*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idVTav(i),ng)) WRITE (out,160)                   &
     &            Aout(idVTav(i),ng), 'Aout(idVTav)',                   &
     &            'Write out averaged <v*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHUTav(i),ng)) WRITE (out,160)                   &
     &            Aout(iHUTav(i),ng), 'Aout(iHUTav)',                   &
     &            'Write out averaged <Huon*t> for tracer ', i,         &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHVTav(i),ng)) WRITE (out,160)                   &
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
     &          WRITE (out,160) .TRUE., 'Dout(iTrate)',                 &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iThadv)',                 &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTxadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTxadv)',                 &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTyadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTyadv)',                 &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTvadv)',                 &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
# if defined TS_DIF2 || defined TS_DIF4
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThdif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iThdif)',                 &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(i,iTxdif),ng))                            &
     &          WRITE (out,160) .TRUE., 'Dout(iTxdif)',                 &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTydif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTydif)',                 &
     &            'Write out horizontal Y-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTsdif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTsdif)',                 &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  endif
# endif
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvdif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTvdif)',                 &
     &            'Write out vertical diffusion, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
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
  40  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' EcoSim Parameters, Grid: ',i2.2,                    &
     &        /,  ' ===========================',/)
  60  FORMAT (1x,i10,2x,a,t32,a)
  70  FORMAT (10x,l1,2x,a,t32,a)
  80  FORMAT ('...........',2x,a,t32,a,/,t34,a)
  90  FORMAT (1p,e11.4,t33,a)
 100  FORMAT ('...........',2x,a,t32,a)
 110  FORMAT (1p,e11.4,t33,'Fecal Group ',i1,', ',a)
 120  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
 130  FORMAT (1p,e11.4,2x,a,t32,a)
 140  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 160  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_BioPar
