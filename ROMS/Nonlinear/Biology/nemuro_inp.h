      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: nemuro_inp.h 790 2016-05-05 19:09:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in Nemuro ecosystem model input parameters.      !
!  They are specified in input script "nemuro.in".                     !
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
!  Read in Nemuro biological model parameters.
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
            CASE ('AttPS')
              Npts=load_r(Nval, Rval, Ngrids, AttPS)
            CASE ('AttPL')
              Npts=load_r(Nval, Rval, Ngrids, AttPL)
            CASE ('PARfrac')
              Npts=load_r(Nval, Rval, Ngrids, PARfrac)
            CASE ('AlphaPS')
              Npts=load_r(Nval, Rval, Ngrids, AlphaPS)
            CASE ('AlphaPL')
              Npts=load_r(Nval, Rval, Ngrids, AlphaPL)
            CASE ('BetaPS')
              Npts=load_r(Nval, Rval, Ngrids, BetaPS)
            CASE ('BetaPL')
              Npts=load_r(Nval, Rval, Ngrids, BetaPL)
            CASE ('VmaxS')
              Npts=load_r(Nval, Rval, Ngrids, VmaxS)
            CASE ('VmaxL')
              Npts=load_r(Nval, Rval, Ngrids, VmaxL)
            CASE ('KNO3S')
              Npts=load_r(Nval, Rval, Ngrids, KNO3S)
            CASE ('KNO3L')
              Npts=load_r(Nval, Rval, Ngrids, KNO3L)
            CASE ('KNH4S')
              Npts=load_r(Nval, Rval, Ngrids, KNH4S)
            CASE ('KNH4L')
              Npts=load_r(Nval, Rval, Ngrids, KNH4L)
            CASE ('KSiL')
              Npts=load_r(Nval, Rval, Ngrids, KSiL)
            CASE ('PusaiS')
              Npts=load_r(Nval, Rval, Ngrids, PusaiS)
            CASE ('PusaiL')
              Npts=load_r(Nval, Rval, Ngrids, PusaiL)
            CASE ('KGppS')
              Npts=load_r(Nval, Rval, Ngrids, KGppS)
            CASE ('KGppL')
              Npts=load_r(Nval, Rval, Ngrids, KGppL)
            CASE ('ResPS0')
              Npts=load_r(Nval, Rval, Ngrids, ResPS0)
            CASE ('ResPL0')
              Npts=load_r(Nval, Rval, Ngrids, ResPL0)
            CASE ('KResPS')
              Npts=load_r(Nval, Rval, Ngrids, KResPS)
            CASE ('KResPL')
              Npts=load_r(Nval, Rval, Ngrids, KResPL)
            CASE ('GammaS')
              Npts=load_r(Nval, Rval, Ngrids, GammaS)
            CASE ('GammaL')
              Npts=load_r(Nval, Rval, Ngrids, GammaL)
            CASE ('MorPS0')
              Npts=load_r(Nval, Rval, Ngrids, MorPS0)
            CASE ('MorPL0')
              Npts=load_r(Nval, Rval, Ngrids, MorPL0)
            CASE ('KMorPS')
              Npts=load_r(Nval, Rval, Ngrids, KMorPS)
            CASE ('KMorPL')
              Npts=load_r(Nval, Rval, Ngrids, KMorPL)
            CASE ('GRmaxSps')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxSps)
            CASE ('GRmaxSpl')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxSpl)
            CASE ('GRmaxLps')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxLps)
            CASE ('GRmaxLpl')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxLpl)
            CASE ('GRmaxLzs')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxLzs)
            CASE ('GRmaxPpl')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxPpl)
            CASE ('GRmaxPzs')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxPzs)
            CASE ('GRmaxPzl')
              Npts=load_r(Nval, Rval, Ngrids, GRmaxPzl)
            CASE ('KGraS')
              Npts=load_r(Nval, Rval, Ngrids, KGraS)
            CASE ('KGraL')
              Npts=load_r(Nval, Rval, Ngrids, KGraL)
            CASE ('KGraP')
              Npts=load_r(Nval, Rval, Ngrids, KGraP)
            CASE ('LamS')
              Npts=load_r(Nval, Rval, Ngrids, LamS)
            CASE ('LamL')
              Npts=load_r(Nval, Rval, Ngrids, LamL)
            CASE ('LamP')
              Npts=load_r(Nval, Rval, Ngrids, LamP)
            CASE ('KPS2ZS')
              Npts=load_r(Nval, Rval, Ngrids, KPS2ZS)
            CASE ('KPL2ZS')
              Npts=load_r(Nval, Rval, Ngrids, KPL2ZS)
            CASE ('KPS2ZL')
              Npts=load_r(Nval, Rval, Ngrids, KPS2ZL)
            CASE ('KPL2ZL')
              Npts=load_r(Nval, Rval, Ngrids, KPL2ZL)
            CASE ('KZS2ZL')
              Npts=load_r(Nval, Rval, Ngrids, KZS2ZL)
            CASE ('KPL2ZP')
              Npts=load_r(Nval, Rval, Ngrids, KPL2ZP)
            CASE ('KZS2ZP')
              Npts=load_r(Nval, Rval, Ngrids, KZS2ZP)
            CASE ('KZL2ZP')
              Npts=load_r(Nval, Rval, Ngrids, KZL2ZP)
            CASE ('PS2ZSstar')
              Npts=load_r(Nval, Rval, Ngrids, PS2ZSstar)
            CASE ('PL2ZSstar')
              Npts=load_r(Nval, Rval, Ngrids, PL2ZSstar)
            CASE ('PS2ZLstar')
              Npts=load_r(Nval, Rval, Ngrids, PS2ZLstar)
            CASE ('PL2ZLstar')
              Npts=load_r(Nval, Rval, Ngrids, PL2ZLstar)
            CASE ('ZS2ZLstar')
              Npts=load_r(Nval, Rval, Ngrids, ZS2ZLstar)
            CASE ('PL2ZPstar')
              Npts=load_r(Nval, Rval, Ngrids, PL2ZPstar)
            CASE ('ZS2ZPstar')
              Npts=load_r(Nval, Rval, Ngrids, ZS2ZPstar)
            CASE ('ZL2ZPstar')
              Npts=load_r(Nval, Rval, Ngrids, ZL2ZPstar)
            CASE ('PusaiPL')
              Npts=load_r(Nval, Rval, Ngrids, PusaiPL)
            CASE ('PusaiZS')
              Npts=load_r(Nval, Rval, Ngrids, PusaiZS)
            CASE ('MorZS0')
              Npts=load_r(Nval, Rval, Ngrids, MorZS0)
            CASE ('MorZL0')
              Npts=load_r(Nval, Rval, Ngrids, MorZL0)
            CASE ('MorZP0')
              Npts=load_r(Nval, Rval, Ngrids, MorZP0)
            CASE ('KMorZS')
              Npts=load_r(Nval, Rval, Ngrids, KMorZS)
            CASE ('KMorZL')
              Npts=load_r(Nval, Rval, Ngrids, KMorZL)
            CASE ('KMorZP')
              Npts=load_r(Nval, Rval, Ngrids, KMorZP)
            CASE ('AlphaZS')
              Npts=load_r(Nval, Rval, Ngrids, AlphaZS)
            CASE ('AlphaZL')
              Npts=load_r(Nval, Rval, Ngrids, AlphaZL)
            CASE ('AlphaZP')
              Npts=load_r(Nval, Rval, Ngrids, AlphaZP)
            CASE ('BetaZS')
              Npts=load_r(Nval, Rval, Ngrids, BetaZS)
            CASE ('BetaZL')
              Npts=load_r(Nval, Rval, Ngrids, BetaZL)
            CASE ('BetaZP')
              Npts=load_r(Nval, Rval, Ngrids, BetaZP)
            CASE ('Nit0')
              Npts=load_r(Nval, Rval, Ngrids, Nit0)
            CASE ('VP2N0')
              Npts=load_r(Nval, Rval, Ngrids, VP2N0)
            CASE ('VP2D0')
              Npts=load_r(Nval, Rval, Ngrids, VP2D0)
            CASE ('VD2N0')
              Npts=load_r(Nval, Rval, Ngrids, VD2N0)
            CASE ('VO2S0')
              Npts=load_r(Nval, Rval, Ngrids, VO2S0)
            CASE ('KNit')
              Npts=load_r(Nval, Rval, Ngrids, KNit)
            CASE ('KP2D')
              Npts=load_r(Nval, Rval, Ngrids, KP2D)
            CASE ('KP2N')
              Npts=load_r(Nval, Rval, Ngrids, KP2N)
            CASE ('KD2N')
              Npts=load_r(Nval, Rval, Ngrids, KD2N)
            CASE ('KO2S')
              Npts=load_r(Nval, Rval, Ngrids, KO2S)
            CASE ('RSiN')
              Npts=load_r(Nval, Rval, Ngrids, RSiN)
            CASE ('setVPON')
              Npts=load_r(Nval, Rval, Ngrids, setVPON)
            CASE ('setVOpal')
              Npts=load_r(Nval, Rval, Ngrids, setVOpal)
# ifdef IRON_LIMIT
            CASE ('T_Fe')
              Npts=load_r(Nval, Rval, Ngrids, T_Fe)
            CASE ('A_Fe')
              Npts=load_r(Nval, Rval, Ngrids, A_Fe)
            CASE ('B_Fe')
              Npts=load_r(Nval, Rval, Ngrids, B_Fe)
            CASE ('SK_FeC')
              Npts=load_r(Nval, Rval, Ngrids, SK_FeC)
            CASE ('LK_FeC')
              Npts=load_r(Nval, Rval, Ngrids, LK_FeC)
            CASE ('FeRR')
              Npts=load_r(Nval, Rval, Ngrids, FeRR)
# endif
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
#ifdef NEMURO_SED1
            CASE ('Hout(idPONsed)')
              IF (idPONsed.eq.0) THEN
                IF (Master) WRITE (out,280) 'idPONsed'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idPONsed,:))
            CASE ('Hout(idOPALsed)')
              IF (idOPALsed.eq.0) THEN
                IF (Master) WRITE (out,280) 'idOPALsed'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idOPALsed,:))
            CASE ('Hout(idDENITsed)')
              IF (idDENITsed.eq.0) THEN
                IF (Master) WRITE (out,280) 'idDENITsed'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idDENITsed,:))
            CASE ('Hout(idPONbur)')
              IF (idPONbur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idPONbur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idPONbur,:))
            CASE ('Hout(idOPALbur)')
              IF (idOPALbur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idOPALbur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idOPALbur,:))
#endif
#ifdef PRIMARY_PROD
            CASE ('Hout(idNPP)')
              IF (idNPP.eq.0) THEN
                IF (Master) WRITE (out,280) 'idNPP'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idNPP,:))
#endif
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
# ifdef NEMURO_SED1
            CASE ('Aout(idPONsed)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idPONsed,:))
            CASE ('Aout(idOPALsed)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idOPALsed,:))
            CASE ('Aout(idDENITsed)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idDENITsed,:))
            CASE ('Aout(idPONbur)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idPONbur,:))
            CASE ('Aout(idOPALbur)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idOPALbur,:))
# endif
# ifdef PRIMARY_PROD
            CASE ('Aout(idNPP)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idNPP,:))
# endif
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
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,70) AttSW(ng), 'AttSW',                          &
     &            'Light attenuation due to seawater (m-1)'
            WRITE (out,80) AttPS(ng), 'AttPS',                          &
     &            'Light attenuation due to small phytoplankton',       &
     &            '(m2/mmole_N).'
            WRITE (out,80) AttPL(ng), 'AttPL',                          &
     &            'Light attenuation due to large phytoplankton',       &
     &            '(m2/mmole_N).'
            WRITE (out,80) PARfrac(ng), 'PARfrac',                      &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,80) AlphaPS(ng), 'AlphaPS',                      &
     &            'Small phytoplankton initial slope of the P-I curve', &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,80) AlphaPL(ng), 'AlphaPL',                      &
     &            'Small phytoplankton initial slope of the P-I curve', &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,80) BetaPS(ng), 'BetaPS',                        &
     &            'Small phytoplankton photoinhibition coefficient',    &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,80) BetaPL(ng), 'BetaPL',                        &
     &            'Large phytoplankton photoinhibition coefficient',    &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,80) VmaxS(ng), 'VmaxS',                          &
     &            'Small phytoplankton maximum photosynthetic rate',    &
     &            '(1/day).'
            WRITE (out,80) VmaxL(ng), 'VmaxL',                          &
     &            'Large phytoplankton maximum photosynthetic rate',    &
     &            '(1/day).'
            WRITE (out,80) KNO3S(ng), 'KNO3S',                          &
     &            'Small phytoplankton NO3 half saturation constant',   &
     &            '(mmole_N/m3).'
            WRITE (out,80) KNO3L(ng), 'KNO3L',                          &
     &            'Large phytoplankton NO3 half saturation constant',   &
     &            '(mmole_N/m3).'
            WRITE (out,80) KNH4S(ng), 'KNH4S',                          &
     &            'Small phytoplankton NH4 half saturation constant',   &
     &            '(mmole_N/m3).'
            WRITE (out,80) KNH4L(ng), 'KNH4L',                          &
     &            'Large phytoplankton NH4 half saturation constant',   &
     &            '(mmole_N/m3).'
            WRITE (out,80) KSiL(ng), 'KSiL',                            &
     &            'Small phytoplankton SiOH4 half saturation constant', &
     &            '(mmole_Si/m3).'
            WRITE (out,80) PusaiS(ng), 'PusaiS',                        &
     &            'Small phytoplankton NH4 inhibition coefficient',     &
     &            '(m3/mmole_N).'
            WRITE (out,80) PusaiL(ng), 'PusaiL',                        &
     &            'Large phytoplankton NH4 inhibition coefficient',     &
     &            '(m3/mmole_N).'
            WRITE (out,80) KGppS(ng), 'KGppS',                          &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'photosynthetic rate (1/Celsius).'
            WRITE (out,80) KGppL(ng), 'KGppL',                          &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'photosynthetic rate (1/Celsius).'
            WRITE (out,70) ResPS0(ng), 'ResPS0',                        &
     &            'Small phytoplankton respiration rate (1/day).'
            WRITE (out,70) ResPL0(ng), 'ResPL0',                        &
     &            'Large phytoplankton respiration rate (1/day).'
            WRITE (out,80) KResPS(ng), 'KResPS',                        &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'respiration (1/Celsius).'
            WRITE (out,80) KResPL(ng), 'KResPL',                        &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'respiration (1/Celsius).'
            WRITE (out,80) GammaS(ng), 'GammaS',                        &
     &            'Small phytoplankton ratio of extracellular',         &
     &            'excretion to photosynthesis (nondimensional).'
            WRITE (out,80) GammaL(ng), 'GammaL',                        &
     &            'Large phytoplankton ratio of extracellular',         &
     &            'excretion to photosynthesis (nondimensional).'
            WRITE (out,80) MorPS0(ng), 'MorPS0',                        &
     &            'Small phytoplankton mortality rate',                 &
     &            '(m3/mmole_N/day).'
            WRITE (out,80) MorPL0(ng), 'MorPL0',                        &
     &            'Large phytoplankton mortality rate',                 &
     &            '(m3/mmole_N/day).'
            WRITE (out,80) KMorPS(ng), 'KMorPS',                        &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'mortality (1/Celsius).'
            WRITE (out,80) KMorPL(ng), 'KMorPL',                        &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'mortality (1/Celsius).'
            WRITE (out,80) GRmaxSps(ng), 'GRmaxSps',                    &
     &            'Small zooplankton grazing rate on small',            &
     &            'phytoplankton (1/day).'
            WRITE (out,80) GRmaxSpl(ng), 'GRmaxSpl',                    &
     &            'Small zooplankton grazing rate on large',            &
     &            'phytoplankton (1/day).'
            WRITE (out,80) GRmaxLps(ng), 'GRmaxLps',                    &
     &            'Large zooplankton grazing rate on small',            &
     &            'phytoplankton (1/day).'
            WRITE (out,80) GRmaxLpl(ng), 'GRmaxLpl',                    &
     &            'Large zooplankton grazing rate on large',            &
     &            'phytoplankton (1/day).'
            WRITE (out,80) GRmaxLzs(ng), 'GRmaxLzs',                    &
     &            'Large zooplankton grazing rate on small',            &
     &            'zooplankton (1/day).'
            WRITE (out,80) GRmaxPpl(ng), 'GRmaxPpl',                    &
     &            'Predator zooplankton grazing rate on large',         &
     &            'phytoplankton (1/day).'
            WRITE (out,80) GRmaxPzs(ng), 'GRmaxPzs',                    &
     &            'Predator zooplankton grazing rate on small',         &
     &            'zooplankton (1/day).'
            WRITE (out,80) GRmaxPzl(ng), 'GRmaxPzl',                    &
     &            'Predator zooplankton grazing rate on large',         &
     &            'zooplankton (1/day).'
            WRITE (out,80) KGraS(ng), 'KGraS',                          &
     &            'Small zooplankton temperature coefficient for',      &
     &            'grazing (1/Celsius).'
            WRITE (out,80) KGraL(ng), 'KGraL',                          &
     &            'Large zooplankton temperature coefficient for',      &
     &            'grazing (1/Celsius).'
            WRITE (out,80) KGraP(ng), 'KGraP',                          &
     &            'Predator zooplankton temperature coefficient for',   &
     &            'grazing (1/Celsius).'
            WRITE (out,80) LamS(ng), 'LamS',                            &
     &            'Small zooplankton grazing Ivlev constant',           &
     &            '(m3/mmole_N).'
            WRITE (out,80) LamL(ng), 'LamL',                            &
     &            'Large zooplankton grazing Ivlev constant',           &
     &            '(m3/mmole_N).'
            WRITE (out,80) LamP(ng), 'LamP',                            &
     &            'Preditor zooplankton grazing Ivlev constant',        &
     &            '(m3/mmole_N).'
#ifdef HOLLING_GRAZING
            WRITE (out,80) KPS2ZS(ng), 'KPS2ZS',                        &
     &            'Half-saturation constant for small zooplankton',     &
     &            'grazing on small phytoplankton (mmole_N/m3)^2.'
            WRITE (out,80) KPL2ZS(ng), 'KPL2ZS',                        &
     &            'Half-saturation constant for small zooplankton',     &
     &            'grazing on large phytoplankton (mmole_N/m3)^2.'
            WRITE (out,80) KPS2ZL(ng), 'KPS2ZL',                        &
     &            'Half-saturation constant for large zooplankton',     &
     &            'grazing on small phytoplankton (mmole_N/m3)^2.'
            WRITE (out,80) KPL2ZL(ng), 'KPL2ZL',                        &
     &            'Half-saturation constant for large zooplankton',     &
     &            'grazing on large phytoplankton (mmole_N/m3)^2.'
            WRITE (out,80) KPL2ZP(ng), 'KPL2ZP',                        &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on large phytoplankton (mmole_N/m3)^2.'
            WRITE (out,80) KZS2ZP(ng), 'KZS2ZP',                        &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on small zooplankton (mmole_N/m3)^2.'
            WRITE (out,80) KZL2ZP(ng), 'KZL2ZP',                        &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on large zooplankton (mmole_N/m3)^2.'
#else
            WRITE (out,80) PS2ZSstar(ng), 'PS2ZSstar',                  &
     &            'Small zooplankton threshold for grazing on small',   &
     &            'phytoplankton (mmole_N/m3).'
            WRITE (out,80) PL2ZSstar(ng), 'PL2ZSstar',                  &
     &            'Small zooplankton threshold for grazing on large',   &
     &            'phytoplankton (mmole_N/m3).'
            WRITE (out,80) PS2ZLstar(ng), 'PS2ZLstar',                  &
     &            'Large zooplankton threshold for grazing on small',   &
     &            'phytoplankton (mmole_N/m3).'
            WRITE (out,80) PL2ZLstar(ng), 'PL2ZLstar',                  &
     &            'Large zooplankton threshold for grazing on large',   &
     &            'phytoplankton (mmole_N/m3).'
            WRITE (out,80) PL2ZLstar(ng), 'PL2ZLstar',                  &
     &            'Large zooplankton threshold for grazing on small',   &
     &            'zooplankton (mmole_N/m3).'
            WRITE (out,80) PL2ZPstar(ng), 'PL2ZPstar',                  &
     &            'Predator zooplankton threshold for grazing on large',&
     &            'phytoplankton (mmole_N/m3).'
            WRITE (out,80) ZS2ZPstar(ng), 'ZS2ZPstar',                  &
     &            'Predator zooplankton threshold for grazing on small',&
     &            'zooplankton (mmole_N/m3).'
            WRITE (out,80) ZL2ZPstar(ng), 'ZL2ZPstar',                  &
     &            'Predator zooplankton threshold for grazing on large',&
     &            'zooplankton (mmole_N/m3).'
#endif
            WRITE (out,80) PusaiPL(ng), 'PusauPL',                      &
     &            'Predator zooplankton grazing inhibition on large',   &
     &            'phytoplankton (mmole_N/m3).'
            WRITE (out,80) PusaiZS(ng), 'PusauZS',                      &
     &            'Predator zooplankton grazing inhibition on small',   &
     &            'zootoplankton (mmole_N/m3).'
            WRITE (out,80) MorZS0(ng), 'MorZS0',                        &
     &            'Small zooplankton mortality rate at 0 Celsius',      &
     &            '(m3/mmole_N/day).'
            WRITE (out,80) MorZL0(ng), 'MorZL0',                        &
     &            'Large zooplankton mortality rate at 0 Celsius',      &
     &            '(m3/mmole_N/day).'
            WRITE (out,80) MorZP0(ng), 'MorZP0',                        &
     &            'Predator zooplankton mortality rate at 0 Celsius',   &
     &            '(m3/mmole_N/day).'
            WRITE (out,80) KMorZS(ng), 'KMorZS',                        &
     &            'Small zooplankton temperature coefficient for',      &
     &            'mortality (1/Celsius).'
            WRITE (out,80) KMorZL(ng), 'KMorZL',                        &
     &            'Large zooplankton temperature coefficient for',      &
     &            'mortality (1/Celsius).'
            WRITE (out,80) KMorZP(ng), 'KMorZP',                        &
     &            'Predator zooplankton temperature coefficient for',   &
     &            'mortality (1/Celsius).'
            WRITE (out,80) AlphaZS(ng), 'AlphaZS',                      &
     &            'Small zooplankton assimilation efficiency',          &
     &            '(nondimensional).'
            WRITE (out,80) AlphaZL(ng), 'AlphaZL',                      &
     &            'Large zooplankton assimilation efficiency',          &
     &            '(nondimensional).'
            WRITE (out,80) AlphaZP(ng), 'AlphaZP',                      &
     &            'Predator zooplankton assimilation efficiency',       &
     &            '(nondimensional).'
            WRITE (out,80) BetaZS(ng), 'BetaZS',                        &
     &            'Small zooplankton growth efficiency',                &
     &            '(nondimensional).'
            WRITE (out,80) BetaZL(ng), 'BetaZL',                        &
     &            'Large zooplankton growth efficiency',                &
     &            '(nondimensional).'
            WRITE (out,80) BetaZP(ng), 'BetaZP',                        &
     &            'Predator zooplankton growth efficiency',             &
     &            '(nondimensional).'
            WRITE (out,70) Nit0(ng), 'Nit0',                            &
     &            'NH4 to NO3 decomposition rate (1/day).'
            WRITE (out,70) VP2N0(ng), 'VP2N0',                          &
     &            'PON to NH4 decomposition rate (1/day).'
            WRITE (out,70) VP2D0(ng), 'VP2D0',                          &
     &            'PON to DON decomposition rate (1/day).'
            WRITE (out,70) VD2N0(ng), 'VD2N0',                          &
     &            'DON to NH4 decomposition rate (1/day).'
            WRITE (out,70) VO2S0(ng), 'VO2S0',                          &
     &            'Opal to SiOH4 decomposition rate (1/day).'
            WRITE (out,80) KNit(ng), 'KNit',                            &
     &            'Temperature coefficient for NH4 to NO3',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,80) KP2D(ng), 'KP2D',                            &
     &            'Temperature coefficient for PON to DON',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,80) KP2N(ng), 'KP2N',                            &
     &            'Temperature coefficient for PON to NH4',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,80) KD2N(ng), 'KD2N',                            &
     &            'Temperature coefficient for DON to NH4',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,80) KO2S(ng), 'KO2S',                            &
     &            'Temperature coefficient for Opal to SiOH4',          &
     &            'decomposition (1/Celsius).'
            WRITE (out,70) RSiN(ng), 'RSiN',                            &
     &            'Si:N ratio (mmole_Si/mmole_N)'
            WRITE (out,70) setVPON(ng), 'setVPON',                      &
     &            'PON sinking velocity (m/day).'
            WRITE (out,70) setVOpal(ng), 'setVOpal',                    &
     &            'Opal sinking velocity (m/day).'
#ifdef IRON_LIMIT
            WRITE (out,70) T_Fe(ng), 'T_Fe',                            &
     &            'Iron uptake time scale (day-1).'
            WRITE (out,70) A_Fe(ng), 'A_Fe',                            &
     &            'Empirical Fe:C power (-).'
            WRITE (out,70) B_Fe(ng), 'B_Fe',                            &
     &            'Empirical Fe:C coefficient (1/M-C).'
            WRITE (out,70) SK_FeC(ng), 'SK_FeC',                        &
     &            'Small P Fe:C at F=0.5 (muM-Fe/M-C).'
            WRITE (out,70) LK_FeC(ng), 'LK_FeC',                        &
     &            'Large P Fe:C at F=0.5 (muM-Fe/M-C).'
            WRITE (out,70) FeRR(ng), 'FeRR',                            &
     &            'Fe remineralization rate (day-1).'
#endif
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu2(i,ng), 'nl_tnu2', i,               &
     &              'NLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu2(i,ng), 'ad_tnu2', i,               &
     &              'ADM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu2(i,ng), 'tl_tnu2', i,               &
     &              'TLM Horizontal, harmonic mixing coefficient',      &
     &              '(m2/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) nl_tnu4(i,ng), 'nl_tnu4', i,               &
     &              'NLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# ifdef ADJOINT
              WRITE (out,90) ad_tnu4(i,ng), 'ad_tnu4', i,               &
     &              'ADM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_tnu4(i,ng), 'tl_tnu4', i,               &
     &              'TLM Horizontal, biharmonic mixing coefficient',    &
     &              '(m4/s) for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,100) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef FORWARD_MIXING
            DO itrc=1,NBT
              i=idbio(itrc)
# ifdef ADJOINT
              WRITE (out,90) ad_Akt_fac(i,ng), 'ad_Akt_fac', i,         &
     &              'ADM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
# if defined TANGENT || defined TL_IOMS
              WRITE (out,90) tl_Akt_fac(i,ng), 'tl_Akt_fac', i,         &
     &              'TLM basic state vertical mixing scale factor',     &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
# endif
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning ON  point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',         &
     &              i, 'Turning OFF point sources/Sink on tracer ', i,  &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,100) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,100) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,100) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,110)                   &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,110)                   &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef PRIMARY_PROD
            IF (Hout(idNPP,ng)) WRITE (out,110)                         &
     &          Hout(idNPP,ng), 'Hout(idNPP)',                          &
     &          'Write out primary productivity', 0,                    &
     &          TRIM(Vname(1,idNPP))
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            WRITE (out,'(1x)')
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTvar(i),ng)) WRITE (out,110)                   &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                   &
     &            'Write out averaged tracer ', i,                      &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
# ifdef PRIMARY_PROD
            IF (Aout(idNPP,ng)) WRITE (out,110)                         &
     &          Aout(idNPP,ng), 'Aout(idNPP)',                          &
     &          'Write out primary productivity', 0,                    &
     &          TRIM(Vname(1,idNPP))
# endif
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTTav(i),ng)) WRITE (out,110)                   &
     &            Aout(idTTav(i),ng), 'Aout(idTTav)',                   &
     &            'Write out averaged <t*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idUTav(i),ng)) WRITE (out,110)                   &
     &            Aout(idUTav(i),ng), 'Aout(idUTav)',                   &
     &            'Write out averaged <u*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idVTav(i),ng)) WRITE (out,110)                   &
     &            Aout(idVTav(i),ng), 'Aout(idVTav)',                   &
     &            'Write out averaged <v*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHUTav(i),ng)) WRITE (out,110)                   &
     &            Aout(iHUTav(i),ng), 'Aout(iHUTav)',                   &
     &            'Write out averaged <Huon*t> for tracer ', i,         &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(iHVTav(i),ng)) WRITE (out,110)                   &
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
     &          WRITE (out,110) .TRUE., 'Dout(iTrate)',                 &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iThadv)',                 &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTxadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTxadv)',                 &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTyadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTyadv)',                 &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvadv),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTvadv)',                 &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
# if defined TS_DIF2 || defined TS_DIF4
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iThdif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iThdif)',                 &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(i,iTxdif),ng))                            &
     &          WRITE (out,110) .TRUE., 'Dout(iTxdif)',                 &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTydif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTydif)',                 &
     &            'Write out horizontal Y-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTsdif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTsdif)',                 &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  endif
# endif
            DO i=1,NBT
              itrc=idbio(i)
              IF (Dout(idDtrc(itrc,iTvdif),ng))                         &
     &          WRITE (out,110) .TRUE., 'Dout(iTvdif)',                 &
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
  50  FORMAT (/,/,' Nemuro Model Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  60  FORMAT (1x,i10,2x,a,t32,a)
  70  FORMAT (1p,e11.4,2x,a,t32,a)
  80  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 100  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 280  FORMAT (/,' READ_BioPar - variable info not yet loaded, ', a)

      RETURN
      END SUBROUTINE read_BioPar
