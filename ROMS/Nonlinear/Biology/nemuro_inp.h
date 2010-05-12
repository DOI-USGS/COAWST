      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: nemuro_inp.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      integer :: Npts, Nval, i, itrc, ng, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in Nemuro biological model parameters.
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
          ELSE IF (TRIM(KeyWord).eq.'AttSW') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSW)
          ELSE IF (TRIM(KeyWord).eq.'AttPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttPS)
          ELSE IF (TRIM(KeyWord).eq.'AttPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttPL)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'AlphaPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaPS)
          ELSE IF (TRIM(KeyWord).eq.'AlphaPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaPL)
          ELSE IF (TRIM(KeyWord).eq.'BetaPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaPS)
          ELSE IF (TRIM(KeyWord).eq.'BetaPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaPL)
          ELSE IF (TRIM(KeyWord).eq.'VmaxS') THEN
            Npts=load_r(Nval, Rval, Ngrids, VmaxS)
          ELSE IF (TRIM(KeyWord).eq.'VmaxL') THEN
            Npts=load_r(Nval, Rval, Ngrids, VmaxL)
          ELSE IF (TRIM(KeyWord).eq.'KNO3S') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNO3S)
          ELSE IF (TRIM(KeyWord).eq.'KNO3L') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNO3L)
          ELSE IF (TRIM(KeyWord).eq.'KNH4S') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNH4S)
          ELSE IF (TRIM(KeyWord).eq.'KNH4L') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNH4L)
          ELSE IF (TRIM(KeyWord).eq.'KSiL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KSiL)
          ELSE IF (TRIM(KeyWord).eq.'PusaiS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiS)
          ELSE IF (TRIM(KeyWord).eq.'PusaiL') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiL)
          ELSE IF (TRIM(KeyWord).eq.'KGppS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGppS)
          ELSE IF (TRIM(KeyWord).eq.'KGppL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGppL)
          ELSE IF (TRIM(KeyWord).eq.'ResPS0') THEN
            Npts=load_r(Nval, Rval, Ngrids, ResPS0)
          ELSE IF (TRIM(KeyWord).eq.'ResPL0') THEN
            Npts=load_r(Nval, Rval, Ngrids, ResPL0)
          ELSE IF (TRIM(KeyWord).eq.'KResPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KResPS)
          ELSE IF (TRIM(KeyWord).eq.'KResPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KResPL)
          ELSE IF (TRIM(KeyWord).eq.'GammaS') THEN
            Npts=load_r(Nval, Rval, Ngrids, GammaS)
          ELSE IF (TRIM(KeyWord).eq.'GammaL') THEN
            Npts=load_r(Nval, Rval, Ngrids, GammaL)
          ELSE IF (TRIM(KeyWord).eq.'MorPS0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorPS0)
          ELSE IF (TRIM(KeyWord).eq.'MorPL0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorPL0)
          ELSE IF (TRIM(KeyWord).eq.'KMorPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorPS)
          ELSE IF (TRIM(KeyWord).eq.'KMorPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorPL)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxSps') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxSps)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxLps') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxLps)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxLpl') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxLpl)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxLzs') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxLzs)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxPpl') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxPpl)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxPzs') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxPzs)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxPzl') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxPzl)
          ELSE IF (TRIM(KeyWord).eq.'KGraS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGraS)
          ELSE IF (TRIM(KeyWord).eq.'KGraL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGraL)
          ELSE IF (TRIM(KeyWord).eq.'KGraP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGraP)
          ELSE IF (TRIM(KeyWord).eq.'LamS') THEN
            Npts=load_r(Nval, Rval, Ngrids, LamS)
          ELSE IF (TRIM(KeyWord).eq.'LamL') THEN
            Npts=load_r(Nval, Rval, Ngrids, LamL)
          ELSE IF (TRIM(KeyWord).eq.'LamP') THEN
            Npts=load_r(Nval, Rval, Ngrids, LamP)
          ELSE IF (TRIM(KeyWord).eq.'KPS2ZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPS2ZS)
          ELSE IF (TRIM(KeyWord).eq.'KPS2ZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPS2ZL)
          ELSE IF (TRIM(KeyWord).eq.'KPL2ZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPL2ZL)
          ELSE IF (TRIM(KeyWord).eq.'KZS2ZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KZS2ZL)
          ELSE IF (TRIM(KeyWord).eq.'KPL2ZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPL2ZP)
          ELSE IF (TRIM(KeyWord).eq.'KZS2ZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KZS2ZP)
          ELSE IF (TRIM(KeyWord).eq.'KZL2ZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KZL2ZP)
          ELSE IF (TRIM(KeyWord).eq.'PS2ZSstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PS2ZSstar)
          ELSE IF (TRIM(KeyWord).eq.'PS2ZLstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PS2ZLstar)
          ELSE IF (TRIM(KeyWord).eq.'PL2ZLstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PL2ZLstar)
          ELSE IF (TRIM(KeyWord).eq.'ZS2ZLstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZS2ZLstar)
          ELSE IF (TRIM(KeyWord).eq.'PL2ZPstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PL2ZPstar)
          ELSE IF (TRIM(KeyWord).eq.'ZS2ZPstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZS2ZPstar)
          ELSE IF (TRIM(KeyWord).eq.'ZL2ZPstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZL2ZPstar)
          ELSE IF (TRIM(KeyWord).eq.'PusaiPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiPL)
          ELSE IF (TRIM(KeyWord).eq.'PusaiZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiZS)
          ELSE IF (TRIM(KeyWord).eq.'MorZS0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorZS0)
          ELSE IF (TRIM(KeyWord).eq.'MorZL0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorZL0)
          ELSE IF (TRIM(KeyWord).eq.'MorZP0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorZP0)
          ELSE IF (TRIM(KeyWord).eq.'KMorZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorZS)
          ELSE IF (TRIM(KeyWord).eq.'KMorZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorZL)
          ELSE IF (TRIM(KeyWord).eq.'KMorZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorZP)
          ELSE IF (TRIM(KeyWord).eq.'AlphaZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaZS)
          ELSE IF (TRIM(KeyWord).eq.'AlphaZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaZL)
          ELSE IF (TRIM(KeyWord).eq.'AlphaZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaZP)
          ELSE IF (TRIM(KeyWord).eq.'BetaZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaZS)
          ELSE IF (TRIM(KeyWord).eq.'BetaZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaZL)
          ELSE IF (TRIM(KeyWord).eq.'BetaZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaZP)
          ELSE IF (TRIM(KeyWord).eq.'Nit0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Nit0)
          ELSE IF (TRIM(KeyWord).eq.'VP2N0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VP2N0)
          ELSE IF (TRIM(KeyWord).eq.'VP2D0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VP2D0)
          ELSE IF (TRIM(KeyWord).eq.'VD2N0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VD2N0)
          ELSE IF (TRIM(KeyWord).eq.'VO2S0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VO2S0)
          ELSE IF (TRIM(KeyWord).eq.'KNit') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNit)
          ELSE IF (TRIM(KeyWord).eq.'KP2D') THEN
            Npts=load_r(Nval, Rval, Ngrids, KP2D)
          ELSE IF (TRIM(KeyWord).eq.'KP2N') THEN
            Npts=load_r(Nval, Rval, Ngrids, KP2N)
          ELSE IF (TRIM(KeyWord).eq.'KD2N') THEN
            Npts=load_r(Nval, Rval, Ngrids, KD2N)
          ELSE IF (TRIM(KeyWord).eq.'KO2S') THEN
            Npts=load_r(Nval, Rval, Ngrids, KO2S)
          ELSE IF (TRIM(KeyWord).eq.'RSiN') THEN
            Npts=load_r(Nval, Rval, Ngrids, RSiN)
          ELSE IF (TRIM(KeyWord).eq.'setVPON') THEN
            Npts=load_r(Nval, Rval, Ngrids, setVPON)
          ELSE IF (TRIM(KeyWord).eq.'setVOpal') THEN
            Npts=load_r(Nval, Rval, Ngrids, setVOpal)
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
     &            '(m2/millimole_N).'
            WRITE (out,80) AttPL(ng), 'AttPL',                          &
     &            'Light attenuation due to large phytoplankton',       &
     &            '(m2/millimole_N).'
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
     &            '(millimole_N/m3).'
            WRITE (out,80) KNO3L(ng), 'KNO3L',                          &
     &            'Large phytoplankton NO3 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,80) KNH4S(ng), 'KNH4S',                          &
     &            'Small phytoplankton NH4 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,80) KNH4L(ng), 'KNH4L',                          &
     &            'Large phytoplankton NH4 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,80) KSiL(ng), 'KSiL',                            &
     &            'Small phytoplankton SiOH4 half saturation constant', &
     &            '(millimole_Si/m3).'
            WRITE (out,80) PusaiS(ng), 'PusaiS',                        &
     &            'Small phytoplankton NH4 inhibition coefficient',     &
     &            '(m3/millimole_N).'
            WRITE (out,80) PusaiL(ng), 'PusaiL',                        &
     &            'Large phytoplankton NH4 inhibition coefficient',     &
     &            '(m3/millimole_N).'
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
     &            '(m3/millimole_N/day).'
            WRITE (out,80) MorPS0(ng), 'MorPL0',                        &
     &            'Large phytoplankton mortality rate',                 &
     &            '(m3/millimole_N/day).'
            WRITE (out,80) KMorPS(ng), 'KMorPS',                        &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'mortality (1/Celsius).'
            WRITE (out,80) KMorPL(ng), 'KMorPL',                        &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'mortality (1/Celsius).'
            WRITE (out,80) GRmaxSps(ng), 'GRmaxSps',                    &
     &            'Small zooplankton grazing rate on small',            &
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
     &            '(m3/millimole_N).'
            WRITE (out,80) LamL(ng), 'LamL',                            &
     &            'Large zooplankton grazing Ivlev constant',           &
     &            '(m3/millimole_N).'
            WRITE (out,80) LamP(ng), 'LamP',                            &
     &            'Preditor zooplankton grazing Ivlev constant',        &
     &            '(m3/millimole_N).'
#ifdef HOLLING_GRAZING
            WRITE (out,80) KPS2ZS(ng), 'KPS2ZS',                        &
     &            'Half-saturation constant for small zooplankton',     &
     &            'grazing on small phytoplankton (millimole_N/m3)^2.'
            WRITE (out,80) KPS2ZL(ng), 'KPS2ZL',                        &
     &            'Half-saturation constant for large zooplankton',     &
     &            'grazing on small phytoplankton (millimole_N/m3)^2.'
            WRITE (out,80) KPL2ZL(ng), 'KPL2ZL',                        &
     &            'Half-saturation constant for large zooplankton',     &
     &            'grazing on large phytoplankton (millimole_N/m3)^2.'
            WRITE (out,80) KPL2ZP(ng), 'KPL2ZP',                        &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on large phytoplankton (millimole_N/m3)^2.'
            WRITE (out,80) KZS2ZP(ng), 'KZS2ZP',                        &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on small zooplankton (millimole_N/m3)^2.'
            WRITE (out,80) KZL2ZP(ng), 'KZL2ZP',                        &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on large zooplankton (millimole_N/m3)^2.'
#else
            WRITE (out,80) PS2ZSstar(ng), 'PS2ZSstar',                  &
     &            'Small zooplankton threshold for grazing on small',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,80) PS2ZLstar(ng), 'PS2ZLstar',                  &
     &            'Large zooplankton threshold for grazing on small',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,80) PL2ZLstar(ng), 'PL2ZLstar',                  &
     &            'Large zooplankton threshold for grazing on large',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,80) PL2ZLstar(ng), 'PL2ZLstar',                  &
     &            'Large zooplankton threshold for grazing on small',   &
     &            'zooplankton (millimole_N/m3).'
            WRITE (out,80) PL2ZPstar(ng), 'PL2ZPstar',                  &
     &            'Predator zooplankton threshold for grazing on large',&
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,80) ZS2ZPstar(ng), 'ZS2ZPstar',                  &
     &            'Predator zooplankton threshold for grazing on small',&
     &            'zooplankton (millimole_N/m3).'
            WRITE (out,80) ZL2ZPstar(ng), 'ZL2ZPstar',                  &
     &            'Predator zooplankton threshold for grazing on large',&
     &            'zooplankton (millimole_N/m3).'
#endif
            WRITE (out,80) PusaiPL(ng), 'PusauPL',                      &
     &            'Predator zooplankton grazing inhibition on large',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,80) PusaiZS(ng), 'PusauZS',                      &
     &            'Predator zooplankton grazing inhibition on small',   &
     &            'zootoplankton (millimole_N/m3).'
            WRITE (out,80) MorZS0(ng), 'MorZS0',                        &
     &            'Small zooplankton mortality rate at 0 Celsius',      &
     &            '(m3/millimole_N/day).'
            WRITE (out,80) MorZL0(ng), 'MorZL0',                        &
     &            'Large zooplankton mortality rate at 0 Celsius',      &
     &            '(m3/millimole_N/day).'
            WRITE (out,80) MorZP0(ng), 'MorZP0',                        &
     &            'Predator zooplankton mortality rate at 0 Celsius',   &
     &            '(m3/millimole_N/day).'
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
            WRITE (out,80) AlphaZS(ng), 'AlphaZL',                      &
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
            WRITE (out,70) VP2N0(ng), 'VO2S0',                          &
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
            WRITE (out,80) KP2N(ng), 'KD2N',                            &
     &            'Temperature coefficient for DON to NH4',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,80) KO2S(ng), 'KO2S',                            &
     &            'Temperature coefficient for Opal to SiOH4',          &
     &            'decomposition (1/Celsius).'
            WRITE (out,70) RSiN(ng), 'RSiN',                            &
     &            'Si:N ratio (millimole_Si/millimole_N)'
            WRITE (out,70) setVPON(ng), 'setVPON',                      &
     &            'PON sinking velocity (m/day).'
            WRITE (out,70) setVOpal(ng), 'setVOpal',                    &
     &            'Opal sinking velocity (m/day).'
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
#ifdef TS_PSOURCE
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,100) LtracerSrc(i,ng), 'LtracerSrc',           &
     &              i, 'Processing point sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
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
     &        /,  ' =============================',/)
  60  FORMAT (1x,i10,2x,a,t30,a)
  70  FORMAT (1p,e11.4,2x,a,t30,a)
  80  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_BioPar
