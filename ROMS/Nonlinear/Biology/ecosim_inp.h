      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: ecosim_inp.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      integer :: Npts, Nval, i, is, itrc, ng, status
      integer :: ibac, iband, ifec, iphy

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in EcoSim bio-optical model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'RtUVR_flag') THEN
            Npts=load_l(Nval, Cval, Ngrids, RtUVR_flag)
          ELSE IF (TRIM(KeyWord).eq.'NFIX_flag') THEN
            Npts=load_l(Nval, Cval, Ngrids, NFIX_flag)
          ELSE IF (TRIM(KeyWord).eq.'Regen_flag') THEN
            Npts=load_l(Nval, Cval, Ngrids, Regen_flag)
          ELSE IF (TRIM(KeyWord).eq.'HsNO3') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsNO3)
          ELSE IF (TRIM(KeyWord).eq.'HsNH4') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsNH4)
          ELSE IF (TRIM(KeyWord).eq.'HsSiO') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsSiO)
          ELSE IF (TRIM(KeyWord).eq.'HsPO4') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsPO4)
          ELSE IF (TRIM(KeyWord).eq.'HsFe') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsFe)
          ELSE IF (TRIM(KeyWord).eq.'GtALG_max') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, GtALG_max)
          ELSE IF (TRIM(KeyWord).eq.'PhyTbase') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, PhyTbase)
          ELSE IF (TRIM(KeyWord).eq.'PhyTfac') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, PhyTfac)
          ELSE IF (TRIM(KeyWord).eq.'BET_') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, BET_)
          ELSE IF (TRIM(KeyWord).eq.'maxC2nALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2nALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2nALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2nALG)
          ELSE IF (TRIM(KeyWord).eq.'C2nALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2nALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'maxC2SiALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2SiALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2SiALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2SiALG)
          ELSE IF (TRIM(KeyWord).eq.'C2SiALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2SiALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'maxC2pALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2pALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2pALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2pALG)
          ELSE IF (TRIM(KeyWord).eq.'C2pALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2pALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'maxC2FeALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2FeALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2FeALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2FeALG)
          ELSE IF (TRIM(KeyWord).eq.'C2FeALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2FeALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'qu_yld') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, qu_yld)
          ELSE IF (TRIM(KeyWord).eq.'E0_comp') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, E0_comp)
          ELSE IF (TRIM(KeyWord).eq.'E0_inhib') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, E0_inhib)
          ELSE IF (TRIM(KeyWord).eq.'inhib_fac') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, inhib_fac)
          ELSE IF (TRIM(KeyWord).eq.'C2CHL_max') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2CHL_max)
          ELSE IF (TRIM(KeyWord).eq.'mxC2Cl') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxC2Cl)
          ELSE IF (TRIM(KeyWord).eq.'b_C2Cl') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_C2Cl)
          ELSE IF (TRIM(KeyWord).eq.'mxC2Cn') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxC2Cn)
          ELSE IF (TRIM(KeyWord).eq.'b_C2Cn') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_C2Cn)
          ELSE IF (TRIM(KeyWord).eq.'mxPacEff') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPacEff)
          ELSE IF (TRIM(KeyWord).eq.'b_PacEff') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PacEff)
          ELSE IF (TRIM(KeyWord).eq.'mxChlB') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxChlB)
          ELSE IF (TRIM(KeyWord).eq.'b_ChlB') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_ChlB)
          ELSE IF (TRIM(KeyWord).eq.'mxChlC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxChlC)
          ELSE IF (TRIM(KeyWord).eq.'b_ChlC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_ChlC)
          ELSE IF (TRIM(KeyWord).eq.'mxPSC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPSC)
          ELSE IF (TRIM(KeyWord).eq.'b_PSC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PSC)
          ELSE IF (TRIM(KeyWord).eq.'mxPPC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPPC)
          ELSE IF (TRIM(KeyWord).eq.'b_PPC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PPC)
          ELSE IF (TRIM(KeyWord).eq.'mxLPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxLPUb)
          ELSE IF (TRIM(KeyWord).eq.'b_LPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_LPUb)
          ELSE IF (TRIM(KeyWord).eq.'mxHPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxHPUb)
          ELSE IF (TRIM(KeyWord).eq.'b_HPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_HPUb)
          ELSE IF (TRIM(KeyWord).eq.'FecDOC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, FecDOC)
          ELSE IF (TRIM(KeyWord).eq.'FecPEL') THEN
            Npts=load_r(Nval, Rval, Nphy*Nfec*Ngrids, FecPEL)
          ELSE IF (TRIM(KeyWord).eq.'FecCYC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, FecCYC)
          ELSE IF (TRIM(KeyWord).eq.'ExALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, ExALG)
          ELSE IF (TRIM(KeyWord).eq.'WS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, WS)
          ELSE IF (TRIM(KeyWord).eq.'HsGRZ') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsGRZ)
          ELSE IF (TRIM(KeyWord).eq.'MinRefuge') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, MinRefuge)
          ELSE IF (TRIM(KeyWord).eq.'RefugeDep') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, RefugeDep)
          ELSE IF (TRIM(KeyWord).eq.'Norm_Vol') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, Norm_Vol)
          ELSE IF (TRIM(KeyWord).eq.'Norm_Surf') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, Norm_Surf)
          ELSE IF (TRIM(KeyWord).eq.'HsDOP') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsDOP)
          ELSE IF (TRIM(KeyWord).eq.'C2pALKPHOS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2pALKPHOS)
          ELSE IF (TRIM(KeyWord).eq.'HsDON') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsDON)
          ELSE IF (TRIM(KeyWord).eq.'C2nNupDON') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2nNupDON)
          ELSE IF (TRIM(KeyWord).eq.'C2nBAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2nBAC)
          ELSE IF (TRIM(KeyWord).eq.'C2pBAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2pBAC)
          ELSE IF (TRIM(KeyWord).eq.'C2FeBAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2FeBAC)
          ELSE IF (TRIM(KeyWord).eq.'HsDOC_ba') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, HsDOC_ba)
          ELSE IF (TRIM(KeyWord).eq.'GtBAC_max') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, GtBAC_max)
          ELSE IF (TRIM(KeyWord).eq.'BacTbase') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, BacTbase)
          ELSE IF (TRIM(KeyWord).eq.'BacTfac') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, BacTfac)
          ELSE IF (TRIM(KeyWord).eq.'BacDOC') THEN
            Npts=load_r(Nval, Rval, Ngrids, BacDOC)
          ELSE IF (TRIM(KeyWord).eq.'BacPEL') THEN
            Npts=load_r(Nval, Rval, Ngrids, BacPEL)
          ELSE IF (TRIM(KeyWord).eq.'BacCYC') THEN
            Npts=load_r(Nval, Rval, Ngrids, BacCYC)
          ELSE IF (TRIM(KeyWord).eq.'ExBAC_c') THEN
            Npts=load_r(Nval, Rval, Ngrids, ExBAC_c)
          ELSE IF (TRIM(KeyWord).eq.'ExBacC2N') THEN
            Npts=load_r(Nval, Rval, Ngrids, ExBacC2N)
          ELSE IF (TRIM(KeyWord).eq.'Bac_Ceff') THEN
            Npts=load_r(Nval, Rval, Ngrids, Bac_Ceff)
          ELSE IF (TRIM(KeyWord).eq.'RtNIT') THEN
            Npts=load_r(Nval, Rval, Ngrids, RtNIT)
          ELSE IF (TRIM(KeyWord).eq.'HsNIT') THEN
            Npts=load_r(Nval, Rval, Ngrids, HsNIT)
          ELSE IF (TRIM(KeyWord).eq.'cDOCfrac_c') THEN
            Npts=load_r(Nval, Rval, Ndom*Ngrids, cDOCfrac_c)
          ELSE IF (TRIM(KeyWord).eq.'RtUVR_DIC') THEN
            Npts=load_r(Nval, Rval, Ngrids, RtUVR_DIC)
          ELSE IF (TRIM(KeyWord).eq.'RtUVR_DOC') THEN
            Npts=load_r(Nval, Rval, Ngrids, RtUVR_DOC)
          ELSE IF (TRIM(KeyWord).eq.'WF') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, WF)
          ELSE IF (TRIM(KeyWord).eq.'RegTbase') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegTbase)
          ELSE IF (TRIM(KeyWord).eq.'RegTfac') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegTfac)
          ELSE IF (TRIM(KeyWord).eq.'RegCR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegCR)
          ELSE IF (TRIM(KeyWord).eq.'RegNR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegNR)
          ELSE IF (TRIM(KeyWord).eq.'RegSR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegSR)
          ELSE IF (TRIM(KeyWord).eq.'RegPR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegPR)
          ELSE IF (TRIM(KeyWord).eq.'RegFR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegFR)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                nl_tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                nl_tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                ad_tnu2(i,ng)=Rbio(itrc,ng)
                tl_tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                ad_tnu4(i,ng)=Rbio(itrc,ng)
                ad_tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_AKT_fac') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                ad_Akt_fac(i,ng)=Rbio(itrc,ng)
                tl_Akt_fac(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
#ifdef TS_PSOURCE
          ELSE IF (TRIM(KeyWord).eq.'LtracerSrc') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                LtracerSrc(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                            &
     &                              'idTvar(idbio(', itrc, '))'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,30)                            &
     &                              'idTsur(idbio(', itrc, '))'
                  exit_flag=5
                  RETURN
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          END IF
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
#ifdef TS_PSOURCE
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',           &
     &              i, 'Processing point sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
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
  60  FORMAT (1x,i10,2x,a,t30,a)
  70  FORMAT (10x,l1,2x,a,t30,a)
  80  FORMAT ('...........',2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,t33,a)
 100  FORMAT ('...........',2x,a,t30,a)
 110  FORMAT (1p,e11.4,t33,'Fecal Group ',i1,', ',a)
 120  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 130  FORMAT (1p,e11.4,2x,a,t30,a)
 140  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 160  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_BioPar
