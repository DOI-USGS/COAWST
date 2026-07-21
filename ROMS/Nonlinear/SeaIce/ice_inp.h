      SUBROUTINE read_IcePar (model, inp, out, Lwrite)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads and reports ice model input parameters.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ice
      USE mod_ncparam
      USE mod_scalars
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
!
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, icevar, igrid, ng, nline, status
      integer :: Bindex
#ifdef AVERAGES
      integer :: nf, ns
#endif
      integer :: tile, LBi, UBi, LBj, UBj
!
      real(r8) :: Rvalue(1)
      real(r8), dimension(200) :: Rval
!
      character (len=40) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Allocate ice model parameters that depend on Ngrids/
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
      tile=MyRank
#else
      tile=0
#endif
      DO ng=1,Ngrids
        LBi=BOUNDS(ng)%LBi(tile)
        UBi=BOUNDS(ng)%UBi(tile)
        LBj=BOUNDS(ng)%LBj(tile)
        UBj=BOUNDS(ng)%UBj(tile)
      END DO
!
      CALL allocate_ice (ng, LBi, UBi, LBj, UBj, .FALSE.)
!
!-----------------------------------------------------------------------
!  Read in ice model parameters.
!-----------------------------------------------------------------------
!
      iceVar=0
      igrid=1
      nline=0
      Nval=0
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('nEVP')
              Npts=load_i(Nval, Rval, Ngrids, nEVP)
            CASE ('AirRho')
              Npts=load_r(Nval, Rval, Ngrids, AirRho)
            CASE ('IceRho')
              Npts=load_r(Nval, Rval, Ngrids, IceRho)
            CASE ('SnowDryRho')
              Npts=load_r(Nval, Rval, Ngrids, SnowDryRho)
            CASE ('SnowWetRho')
              Npts=load_r(Nval, Rval, Ngrids, SnowWetRho)
            CASE ('Cd_ai')
              Npts=load_r(Nval, Rval, Ngrids, Cd_ai)
            CASE ('Cd_io')
              Npts=load_r(Nval, Rval, Ngrids, Cd_io)
            CASE ('Astrength')
              Npts=load_r(Nval, Rval, Ngrids, Astrength)
            CASE ('Pstar')
              Npts=load_r(Nval, Rval, Ngrids, Pstar)
            CASE ('zetaMin')
              Npts=load_r(Nval, Rval, Ngrids, zetaMin)
            CASE ('zetaMax')
              Npts=load_r(Nval, Rval, Ngrids, zetaMax)
            CASE ('ellip_sq')
              Npts=load_r(Nval, Rval, Ngrids, ellip_sq)
            CASE ('min_ai')
              Npts=load_r(Nval, Rval, Ngrids, min_ai)
            CASE ('max_ai')
              Npts=load_r(Nval, Rval, Ngrids, max_ai)
            CASE ('min_hi')
              Npts=load_r(Nval, Rval, Ngrids, min_hi)
            CASE ('max_hmelt')
              Npts=load_r(Nval, Rval, Ngrids, max_hmelt)
            CASE ('stressAng')
              Npts=load_r(Nval, Rval, Ngrids, stressAng)
              DO ng=1,Ngrids
                 stressAng(ng) = stressAng(ng)*deg2rad
              END DO
            CASE ('ice_emiss')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              ice_emiss=Rvalue(1)
            CASE ('spec_heat_air')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              spec_heat_air=Rvalue(1)
            CASE ('trans_coeff')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              trans_coeff=Rvalue(1)
            CASE ('sublimation')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              sublimation=Rvalue(1)
            CASE ('LBC(isAice)')
              Bindex=isAice+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idAice), LBC)
            CASE ('LBC(isHice)')
              Bindex=isHice+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idHice), LBC)
            CASE ('LBC(isHmel)')
              Bindex=isHmel+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idHmel), LBC)
            CASE ('LBC(isHsno)')
              Bindex=isHsno+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idHsno), LBC)
            CASE ('LBC(isIage)')
              Bindex=isIage+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idIage), LBC)
            CASE ('LBC(isISxx)')
              Bindex=isISxx+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idISxx), LBC)
            CASE ('LBC(isISxy)')
              Bindex=isISxy+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idISxy), LBC)
            CASE ('LBC(isISyy)')
              Bindex=isISyy+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idISyy), LBC)
            CASE ('LBC(isTice)')
              Bindex=isTice+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idTice), LBC)
            CASE ('LBC(isUice)')
              Bindex=isUice+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idUice), LBC)
            CASE ('LBC(isVice)')
              Bindex=isVice+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idVice), LBC)
            CASE ('LBC(isEnth)')
              Bindex=isEnth+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idEnth), LBC)
            CASE ('LBC(isHage)')
              Bindex=isHage+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idHage), LBC)
            CASE ('LBC(isUevp)')
              Bindex=isUevp+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idUevp), LBC)
            CASE ('LBC(isVevp)')
              Bindex=isVevp+(nLBCvar-nIceS)
              Npts=load_lbc(Nval, Cval, line, nline, Bindex, igrid,     &
     &                      0, 0, Vname(1,idVevp), LBC)
            CASE ('Hout(idUice)')
              IF (idUice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUice,:))
            CASE ('Hout(idVice)')
              IF (idVice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idVice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVice,:))
            CASE ('Hout(idUiER)')
              IF (idUiER.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUiER'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUiER,:))
            CASE ('Hout(idViNR)')
              IF (idViNR.eq.0) THEN
                IF (Master) WRITE (out,80) 'idViNR'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idViNR,:))
            CASE ('Hout(idAice)')
              IF (idAice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idAice,:))
            CASE ('Hout(idIage)')
              IF (idIage.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIage'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIage,:))
            CASE ('Hout(idHice)')
              IF (idHice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHice,:))
            CASE ('Hout(idHmel)')
              IF (idHmel.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHmel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHmel,:))
            CASE ('Hout(idHsno)')
              IF (idHsno.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHsno'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHsno,:))
            CASE ('Hout(idTice)')
              IF (idTice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTice,:))
            CASE ('Hout(idISxx)')
              IF (idISxx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idISxx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idISxx,:))
            CASE ('Hout(idISxy)')
              IF (idISxy.eq.0) THEN
                IF (Master) WRITE (out,80) 'idISxy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idISxy,:))
            CASE ('Hout(idISyy)')
              IF (idISyy.eq.0) THEN
                IF (Master) WRITE (out,80) 'idISyy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idISyy,:))
            CASE ('Hout(idEnth)')
              IF (idEnth.eq.0) THEN
                IF (Master) WRITE (out,80) 'idEnth'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idEnth,:))
            CASE ('Hout(idHage)')
              IF (idHage.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHage'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHage,:))
            CASE ('Hout(idUevp)')
              IF (idUevp.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUevp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUevp,:))
            CASE ('Hout(idVevp)')
              IF (idVevp.eq.0) THEN
                IF (Master) WRITE (out,80) 'idVevp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVevp,:))
            CASE ('Hout(idIsst)')
              IF (idIsst.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIsst'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIsst,:))
            CASE ('Hout(idIOmf)')
              IF (idIOmf.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIOmf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIOmf,:))
            CASE ('Hout(idIOfv)')
              IF (idIOfv.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIOfv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIOfv,:))
            CASE ('Hout(idIOmt)')
              IF (idIOmt.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIOmt'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIOmt,:))
            CASE ('Hout(idS0mk)')
              IF (idS0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idS0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idS0mk,:))
            CASE ('Hout(idT0mk)')
              IF (idT0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idT0mk,:))
            CASE ('Hout(idWdiv)')
              IF (idWdiv.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWdiv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdiv,:))
            CASE ('Hout(idW_fr)')
              IF (idW_fr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_fr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW_fr,:))
            CASE ('Hout(idW_ai)')
              IF (idW_ai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_ai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW_ai,:))
            CASE ('Hout(idW_ao)')
              IF (idW_ao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_ao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW_ao,:))
            CASE ('Hout(idW_io)')
              IF (idW_io.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_io'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW_io,:))
            CASE ('Hout(idW_ro)')
              IF (idW_ro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idW_ro'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW_ro,:))
            CASE ('Hout(idQcon)')
              IF (idQcon.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQcon'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idQcon,:))
            CASE ('Hout(idQrhs)')
              IF (idQrhs.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQrhs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idQrhs,:))
            CASE ('Qout(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idUice,:))
            CASE ('Qout(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idVice,:))
            CASE ('Qout(idUiER)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idUiER,:))
            CASE ('Qout(idViNR)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idViNR,:))
            CASE ('Qout(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idAice,:))
            CASE ('Qout(idIage)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIage,:))
            CASE ('Qout(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idHice,:))
            CASE ('Qout(idHmel)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idHmel,:))
            CASE ('Qout(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idHsno,:))
            CASE ('Qout(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTice,:))
            CASE ('Qout(idISxx)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idISxx,:))
            CASE ('Qout(idISxy)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idISxy,:))
            CASE ('Qout(idISyy)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idISyy,:))
            CASE ('Qout(idEnth)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idEnth,:))
            CASE ('Qout(idHage)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idHage,:))
            CASE ('Qout(idUevp)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idUevp,:))
            CASE ('Qout(idVevp)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idVevp,:))
            CASE ('Qout(idIsst)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIsst,:))
            CASE ('Qout(idIOmf)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIOmf,:))
            CASE ('Qout(idIOfv)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIOfv,:))
            CASE ('Qout(idIOmt)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idIOmt,:))
            CASE ('Qout(idQcon)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idQcon,:))
            CASE ('Qout(idQrhs)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idQrhs,:))
            CASE ('Qout(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idS0mk,:))
            CASE ('Qout(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idT0mk,:))
            CASE ('Qout(idWdiv)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idWdiv,:))
            CASE ('Qout(idW_fr)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idW_fr,:))
            CASE ('Qout(idW_ai)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idW_ai,:))
            CASE ('Qout(idW_ao)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idW_ao,:))
            CASE ('Qout(idW_io)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idW_io,:))
            CASE ('Qout(idW_ro)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idW_ro,:))
#ifdef AVERAGES
            CASE ('Aout(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUice,:))
            CASE ('Aout(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVice,:))
            CASE ('Aout(idUiER)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUiER,:))
            CASE ('Aout(idViNR)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idViNR,:))
            CASE ('Aout(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idAice,:))
            CASE ('Aout(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHice,:))
            CASE ('Aout(idHmel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHmel,:))
            CASE ('Aout(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHsno,:))
            CASE ('Aout(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idTice,:))
            CASE ('Aout(idISxx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idISxx,:))
            CASE ('Aout(idISxy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idISxy,:))
            CASE ('Aout(idISyy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idISyy,:))
            CASE ('Aout(idEnth)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idEnth,:))
            CASE ('Aout(idHage)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHage,:))
            CASE ('Aout(idUevp)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUevp,:))
            CASE ('Aout(idVevp)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVevp,:))
            CASE ('Aout(idIsst)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIsst,:))
            CASE ('Aout(idIOmf)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIOmf,:))
            CASE ('Aout(idIOfv)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIOfv,:))
            CASE ('Aout(idIOmt)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIOmt,:))
            CASE ('Aout(idIage)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIage,:))
            CASE ('Aout(idQcon)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idQcon,:))
            CASE ('Aout(idQrhs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idQrhs,:))
            CASE ('Aout(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idS0mk,:))
            CASE ('Aout(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idT0mk,:))
            CASE ('Aout(idWdiv)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdiv,:))
            CASE ('Aout(idW_fr)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW_fr,:))
            CASE ('Aout(idW_ai)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW_ai,:))
            CASE ('Aout(idW_ao)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW_ao,:))
            CASE ('Aout(idW_io)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW_io,:))
            CASE ('Aout(idW_ro)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW_ro,:))
#endif
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)

! Set ice time step to ocean time step

      DO ng = 1,Ngrids
        dtice(ng) = dt(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          WRITE (out,40) ng
          WRITE (out,50) nEVP(ng), 'nEVP',                              &
     &          'Number of elastic steps per plastic step in EVP.'
          WRITE (out,60) AirRho(ng), 'AirRho',                          &
     &          'Air density (kg/m3).'
          WRITE (out,60) IceRho(ng), 'IceRho',                          &
     &          'Density of sea ice (kg/m3).'
          WRITE (out,60) SnowDryRho(ng), 'SnowDryRho',                  &
     &          'Dry snow density (kg/m3).'
          WRITE (out,60) SnowWetRho(ng), 'SnowWetRho',                  &
     &          'Wet snow density (kg/m3).'
          WRITE (out,60) Cd_ai(ng), 'Cd_ai',                            &
     &          'Air-Ice drag coefficient (nondimensional).'
          WRITE (out,60) Cd_io(ng), 'Cd_io',                            &
     &          'Ice-Ocean drag coefficient (nondimensional).'
          WRITE (out,60) Astrength(ng), 'Astrength',                    &
     &          'Ice strength exponential weighting (nondimensional).'
          WRITE (out,60) Pstar(ng), 'Pstar',                            &
     &          'Compressive ice strength parameter, Hibler (1979).'
          WRITE (out,60) ZetaMin(ng), 'zetaMin',                        &
     &          'Minimum ice shear strength (N/m2) limiter.'
          WRITE (out,60) ZetaMax(ng), 'zetaMax',                        &
     &          'Maximum ice shear strength (N/m2) limiter.'
          WRITE (out,60) ellip_sq(ng), 'ellip_sq',                      &
     &          'Ellipticity squared of yield curve (nondimensional).'
          WRITE (out,60) min_ai(ng), 'min_ai',                          &
     &          'Minimum ice concentration (nondimensional) limiter.'
          WRITE (out,60) max_ai(ng), 'max_ai',                          &
     &          'Maximum ice concentration (nondimensional) limiter.'
          WRITE (out,60) min_hi(ng), 'min_hi',                          &
     &          'Minimum average ice thickness (m) limiter.'
          WRITE (out,60) max_hmelt(ng), 'max_hmelt',                    &
     &          'Maximum surface melt water thickness (m) limiter.'
          WRITE (out,60) stressAng(ng)*rad2deg, 'stressAng',            &
     &          'Turning angle for ice-water drag (degrees).'
          WRITE (out,60) ice_emiss, 'ice_emiss',                        &
     &          'Ice emissivity (nondimensional).'
          WRITE (out,60) spec_heat_air, 'spec_heat_air',                &
     &          'Specific heat of air (J/kg/K).'
          WRITE (out,60) trans_coeff, 'trans_coeff',                    &
     &          'Heat transfer coefficient (nondimensional).'
          WRITE (out,60) sublimation, 'sublimation',                    &
     &          'Latent heat of sublimation (J/kg).'
!
          IF (Hout(idUice,ng)) WRITE (out,70) Hout(idUice,ng),          &
     &       'Hout(idUice)',                                            &
     &       'Write out ice U-velocity component.'
          IF (Hout(idVice,ng)) WRITE (out,70) Hout(idVice,ng),          &
     &       'Hout(idVice)',                                            &
     &       'Write out ice V-velocity component.'
          IF (Hout(idUiER,ng)) WRITE (out,70) Hout(idUiER,ng),          &
     &       'Hout(idUiER)',                                            &
     &       'Write out ice Eastward velocity component at RHO-points.'
          IF (Hout(idViNR,ng)) WRITE (out,70) Hout(idViNR,ng),          &
     &       'Hout(idViNR)',                                            &
     &       'Write out ice Northward velocity component at RHO-points.'
          IF (Hout(idAice,ng)) WRITE (out,70) Hout(idAice,ng),          &
     &       'Hout(idAice)',                                            &
     &       'Write out ice concentration (fractional area coverage).'
          IF (Hout(idIage,ng)) WRITE (out,70) Hout(idIage,ng),          &
     &       'Hout(idIage)',                                            &
     &       'Write out age of ice.'
          IF (Hout(idHice,ng)) WRITE (out,70) Hout(idHice,ng),          &
     &       'Hout(idHice)',                                            &
     &       'Write out average ice thickness (ice mass per area).'
          IF (Hout(idHmel,ng)) WRITE (out,70) Hout(idHmel,ng),          &
     &       'Hout(idHmel)',                                            &
     &       'Write out melt pond water thickness on ice.'
          IF (Hout(idHsno,ng)) WRITE (out,70) Hout(idHsno,ng),          &
     &       'Hout(idHsno)',                                            &
     &       'Write out average snow coverage thickness.'
          IF (Hout(idTice,ng)) WRITE (out,70) Hout(idTice,ng),          &
     &       'Hout(idTice)',                                            &
     &       'Write out interior ice temperature.'
          IF (Hout(idISxx,ng)) WRITE (out,70) Hout(idISxx,ng),          &
     &       'Hout(idISxx)',                                            &
     &       'Write out internal ice stress tensor, xx-component.'
          IF (Hout(idISxy,ng)) WRITE (out,70) Hout(idISxy,ng),          &
     &       'Hout(idISxy)',                                            &
     &       'Write out internal ice stress tensor, xy-component.'
          IF (Hout(idISyy,ng)) WRITE (out,70) Hout(idISyy,ng),          &
     &       'Hout(idISyy)',                                            &
     &       'Write out internal ice stress tensor, yy-component.'
          IF (Hout(idEnth,ng)) WRITE (out,70) Hout(idEnth,ng),          &
     &       'Hout(idEnth)',                                            &
     &       'Write out ice-brine enthalpy.'
          IF (Hout(idHage,ng)) WRITE (out,70) Hout(idHage,ng),          &
     &       'Hout(idHage)',                                            &
     &       'Write out thickness linked sea ice age.'
          IF (Hout(idUevp,ng)) WRITE (out,70) Hout(idUevp,ng),          &
     &       'Hout(idUevp)',                                            &
     &       'Write out elastic-viscous-plastic rheology U-velocity.'
          IF (Hout(idVevp,ng)) WRITE (out,70) Hout(idVevp,ng),          &
     &       'Hout(idVevp)',                                            &
     &       'Write out elastic-viscous-plastic rheology V-velocity.'
          IF (Hout(idIsst,ng)) WRITE (out,70) Hout(idIsst,ng),          &
     &       'Hout(idIsst)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Hout(idIOmf,ng)) WRITE (out,70) Hout(idIOmf,ng),          &
     &       'Hout(idIOmf)',                                            &
     &       'Write out ice-ocean mass flux.'
          IF (Hout(idIOfv,ng)) WRITE (out,70) Hout(idIOfv,ng),          &
     &       'Hout(idIOfv)',                                            &
     &       'Write out ice-ocean friction velocity.'
          IF (Hout(idIOmt,ng)) WRITE (out,70) Hout(idIOmt,ng),          &
     &       'Hout(idIOmt)',                                            &
     &       'Write out ice-ocean momentum transfer coefficient.'
          IF (Hout(idQcon,ng)) WRITE (out,70) Hout(idQcon,ng),          &
     &       'Hout(idQcon)',                                            &
     &       'Write out ice-snow heat conductivity.'
          IF (Hout(idQrhs,ng)) WRITE (out,70) Hout(idQrhs,ng),          &
     &       'Hout(idQrhs)',                                            &
     &       'Write out RHS heat flux over ice/snow.'
          IF (Hout(idS0mk,ng)) WRITE (out,70) Hout(idS0mk,ng),          &
     &       'Hout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
          IF (Hout(idT0mk,ng)) WRITE (out,70) Hout(idT0mk,ng),          &
     &       'Hout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Hout(idWdiv,ng)) WRITE (out,70) Hout(idWdiv,ng),          &
     &       'Hout(idWdiv)',                                            &
     &       'Write out rate of ice divergence.'
          IF (Hout(idW_fr,ng)) WRITE (out,70) Hout(idW_fr,ng),          &
     &       'Hout(idW_fr)',                                            &
     &       'Write out rate of ice accretion by frazil ice growth.'
          IF (Hout(idW_ai,ng)) WRITE (out,70) Hout(idW_ai,ng),          &
     &       'Hout(idW_ai)',                                            &
     &       'Write out rate of melt/freeze at air/ice interface.'
          IF (Hout(idW_ao,ng)) WRITE (out,70) Hout(idW_ao,ng),          &
     &       'Hout(idW_ao)',                                            &
     &       'Write out rate of melt/freeze at air/ocean interface.'
          IF (Hout(idW_io,ng)) WRITE (out,70) Hout(idW_io,ng),          &
     &       'Hout(idW_io)',                                            &
     &       'Write out rate of melt/freeze at ice/ocean interface.'
          IF (Hout(idW_ro,ng)) WRITE (out,70) Hout(idW_ro,ng),          &
     &       'Hout(idW_ro)',                                            &
     &       'Write out rate of melt/freeze runoff into ocean.'
!
          IF (Qout(idUice,ng)) WRITE (out,70) Qout(idUice,ng),          &
     &       'Qout(idUice)',                                            &
     &       'Write out ice U-velocity component.'
          IF (Qout(idVice,ng)) WRITE (out,70) Qout(idVice,ng),          &
     &       'Qout(idVice)',                                            &
     &       'Write out ice V-velocity component.'
          IF (Qout(idUiER,ng)) WRITE (out,70) Qout(idUiER,ng),          &
     &       'Qout(idUiER)',                                            &
     &       'Write out ice Eastward velocity component at RHO-points.'
          IF (Qout(idViNR,ng)) WRITE (out,70) Qout(idViNR,ng),          &
     &       'Qout(idViNR)',                                            &
     &       'Write out ice Northward velocity component at RHO-points.'
          IF (Qout(idAice,ng)) WRITE (out,70) Qout(idAice,ng),          &
     &       'Qout(idAice)',                                            &
     &       'Write out ice concentration (fractional area coverage).'
          IF (Qout(idIage,ng)) WRITE (out,70) Qout(idIage,ng),          &
     &       'Qout(idIage)',                                            &
     &       'Write out age of ice.'
          IF (Qout(idHice,ng)) WRITE (out,70) Qout(idHice,ng),          &
     &       'Qout(idHice)',                                            &
     &       'Write out average ice thickness (ice mass per area).'
          IF (Qout(idHmel,ng)) WRITE (out,70) Qout(idHmel,ng),          &
     &       'Qout(idHmel)',                                            &
     &       'Write out surface melt water thickness on ice.'
          IF (Qout(idHsno,ng)) WRITE (out,70) Qout(idHsno,ng),          &
     &       'Qout(idHsno)',                                            &
     &       'Write out average snow coverage thickness.'
          IF (Qout(idTice,ng)) WRITE (out,70) Qout(idTice,ng),          &
     &       'Qout(idTice)',                                            &
     &       'Write out interior ice temperature.'
          IF (Qout(idISxx,ng)) WRITE (out,70) Qout(idISxx,ng),          &
     &       'Qout(idISxx)',                                            &
     &       'Write out internal ice stress tensor, xx-component.'
          IF (Qout(idISxy,ng)) WRITE (out,70) Qout(idISxy,ng),          &
     &       'Qout(idISxy)',                                            &
     &       'Write out internal ice stress tensor, xy-component.'
          IF (Qout(idISyy,ng)) WRITE (out,70) Qout(idISyy,ng),          &
     &       'Qout(idISyy)',                                            &
     &       'Write out internal ice stress tensor, yy-component.'
          IF (Qout(idEnth,ng)) WRITE (out,70) Qout(idEnth,ng),          &
     &       'Qout(idEnth)',                                            &
     &       'Write out ice-brine enthalpy.'
          IF (Qout(idHage,ng)) WRITE (out,70) Qout(idHage,ng),          &
     &       'Qout(idHage)',                                            &
     &       'Write out thickness linked sea ice age.'
          IF (Qout(idUevp,ng)) WRITE (out,70) Qout(idUevp,ng),          &
     &       'Qout(idUevp)',                                            &
     &       'Write out elastic-viscous-plastic rheology U-velocity.'
          IF (Qout(idVevp,ng)) WRITE (out,70) Qout(idVevp,ng),          &
     &       'Qout(idVevp)',                                            &
     &       'Write out elastic-viscous-plastic rheology V-velocity.'
          IF (Qout(idIsst,ng)) WRITE (out,70) Qout(idIsst,ng),          &
     &       'Qout(idIsst)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Qout(idIOmf,ng)) WRITE (out,70) Qout(idIOmf,ng),          &
     &       'Qout(idIOmf)',                                            &
     &       'Write out ice-ocean mass flux.'
          IF (Qout(idIOfv,ng)) WRITE (out,70) Qout(idIOfv,ng),          &
     &       'Qout(idIOfv)',                                            &
     &       'Write out ice-ocean friction velocity.'
          IF (Qout(idIOmt,ng)) WRITE (out,70) Qout(idIOmt,ng),          &
     &       'Qout(idIOmt)',                                            &
     &       'Write out ice-ocean momentum transfer coefficient.'
          IF (Qout(idQcon,ng)) WRITE (out,70) Qout(idQcon,ng),          &
     &       'Qout(idQcon)',                                            &
     &       'Write out ice-snow heat conductivity.'
          IF (Qout(idQrhs,ng)) WRITE (out,70) Qout(idQrhs,ng),          &
     &       'Qout(idQrhs)',                                            &
     &       'Write out RHS heat flux over ice/snow.'
          IF (Qout(idS0mk,ng)) WRITE (out,70) Qout(idS0mk,ng),          &
     &       'Qout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
          IF (Qout(idT0mk,ng)) WRITE (out,70) Qout(idT0mk,ng),          &
     &       'Qout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Qout(idWdiv,ng)) WRITE (out,70) Qout(idWdiv,ng),          &
     &       'Qout(idWdiv)',                                            &
     &       'Write out rate of ice divergence.'
          IF (Qout(idW_fr,ng)) WRITE (out,70) Qout(idW_fr,ng),          &
     &       'Qout(idW_fr)',                                            &
     &       'Write out rate of ice accretion by frazil ice growth.'
          IF (Qout(idW_ai,ng)) WRITE (out,70) Qout(idW_ai,ng),          &
     &       'Qout(idW_ai)',                                            &
     &       'Write out rate of melt/freeze at air/ice interface.'
          IF (Qout(idW_ao,ng)) WRITE (out,70) Qout(idW_ao,ng),          &
     &       'Qout(idW_ao)',                                            &
     &       'Write out rate of melt/freeze at air/ocean interface.'
          IF (Qout(idW_io,ng)) WRITE (out,70) Qout(idW_io,ng),          &
     &       'Qout(idW_io)',                                            &
     &       'Write out rate of melt/freeze at ice/ocean interface.'
          IF (Qout(idW_ro,ng)) WRITE (out,70) Qout(idW_ro,ng),          &
     &       'Qout(idW_ro)',                                            &
     &       'Write out rate of melt/freeze runoff into ocean.'
#ifdef AVERAGES
!
          IF (Aout(idUice,ng)) WRITE (out,70) Aout(idUice,ng),          &
     &       'Aout(idUice)',                                            &
     &       'Write out ice U-velocity component.'
          IF (Aout(idVice,ng)) WRITE (out,70) Aout(idVice,ng),          &
     &       'Aout(idVice)',                                            &
     &       'Write out ice V-velocity component.'
          IF (Aout(idUiER,ng)) WRITE (out,70) Aout(idUiER,ng),          &
     &       'Aout(idUiER)',                                            &
     &       'Write out ice Eastward velocity component at RHO-points.'
          IF (Aout(idViNR,ng)) WRITE (out,70) Aout(idViNR,ng),          &
     &       'Aout(idViNR)',                                            &
     &       'Write out ice Northward velocity component at RHO-points.'
          IF (Aout(idAice,ng)) WRITE (out,70) Aout(idAice,ng),          &
     &       'Aout(idAice)',                                            &
     &       'Write out ice concentration (fractional area coverage).'
          IF (Aout(idIage,ng)) WRITE (out,70) Aout(idIage,ng),          &
     &       'Aout(idIage)',                                            &
     &       'Write out age of ice.'
          IF (Aout(idHice,ng)) WRITE (out,70) Aout(idHice,ng),          &
     &       'Aout(idHice)',                                            &
     &       'Write out average ice thickness (ice mass per area).'
          IF (Aout(idHmel,ng)) WRITE (out,70) Aout(idHmel,ng),          &
     &       'Aout(idHmel)',                                            &
     &       'Write out surface melt water thickness on ice.'
          IF (Aout(idHsno,ng)) WRITE (out,70) Aout(idHsno,ng),          &
     &       'Aout(idHsno)',                                            &
     &       'Write out average snow coverage thickness.'
          IF (Aout(idTice,ng)) WRITE (out,70) Aout(idTice,ng),          &
     &       'Aout(idTice)',                                            &
     &       'Write out interior ice temperature.'
          IF (Aout(idISxx,ng)) WRITE (out,70) Aout(idISxx,ng),          &
     &       'Aout(idISxx)',                                            &
     &       'Write out internal ice stress tensor, xx-component.'
          IF (Aout(idISxy,ng)) WRITE (out,70) Aout(idISxy,ng),          &
     &       'Aout(idISxy)',                                            &
     &       'Write out internal ice stress tensor, xy-component.'
          IF (Aout(idISyy,ng)) WRITE (out,70) Aout(idISyy,ng),          &
     &       'Aout(idISyy)',                                            &
     &       'Write out internal ice stress tensor, yy-component.'
          IF (Aout(idEnth,ng)) WRITE (out,70) Aout(idEnth,ng),          &
     &       'Aout(idEnth)',                                            &
     &       'Write out ice-brine enthalpy.'
          IF (Aout(idHage,ng)) WRITE (out,70) Aout(idHage,ng),          &
     &       'Aout(idHage)',                                            &
     &       'Write out thickness linked sea ice age.'
          IF (Aout(idUevp,ng)) WRITE (out,70) Aout(idUevp,ng),          &
     &       'Aout(idUevp)',                                            &
     &       'Write out elastic-viscous-plastic rheology U-velocity.'
          IF (Aout(idVevp,ng)) WRITE (out,70) Aout(idVevp,ng),          &
     &       'Aout(idVevp)',                                            &
     &       'Write out elastic-viscous-plastic rheology V-velocity.'
          IF (Aout(idIsst,ng)) WRITE (out,70) Aout(idisst,ng),          &
     &       'Aout(idIsst)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Aout(idIOmf,ng)) WRITE (out,70) Aout(idIOmf,ng),          &
     &       'Aout(idIOmf)',                                            &
     &       'Write out ice-ocean mass flux.'
          IF (Aout(idIOfv,ng)) WRITE (out,70) Aout(idIOfv,ng),          &
     &       'Aout(idIOfv)',                                            &
     &       'Write out ice-ocean friction velocity.'
          IF (Aout(idIOmt,ng)) WRITE (out,70) Aout(idIOmt,ng),          &
     &       'Aout(idIOmt)',                                            &
     &       'Write out ice-ocean momentum transfer coefficient.'
          IF (Aout(idQcon,ng)) WRITE (out,70) Aout(idQcon,ng),          &
     &       'Aout(idQcon)',                                            &
     &       'Write out ice-snow heat conductivity.'
          IF (Aout(idQrhs,ng)) WRITE (out,70) Aout(idQrhs,ng),          &
     &       'Aout(idQrhs)',                                            &
     &       'Write out RHS heat flux over ice/snow.'
          IF (Aout(idS0mk,ng)) WRITE (out,70) Aout(idS0mk,ng),          &
     &       'Aout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
          IF (Aout(idT0mk,ng)) WRITE (out,70) Aout(idT0mk,ng),          &
     &       'Aout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Aout(idWdiv,ng)) WRITE (out,70) Aout(idWdiv,ng),          &
     &       'Aout(idWdiv)',                                            &
     &       'Write out rate of ice divergence.'
          IF (Aout(idW_fr,ng)) WRITE (out,70) Aout(idW_fr,ng),          &
     &       'Aout(idW_fr)',                                            &
     &       'Write out rate of ice accretion by frazil ice growth.'
          IF (Aout(idW_ai,ng)) WRITE (out,70) Aout(idW_ai,ng),          &
     &       'Aout(idW_ai)',                                            &
     &       'Write out rate of melt/freeze at air/ice interface.'
          IF (Aout(idW_ao,ng)) WRITE (out,70) Aout(idW_ao,ng),          &
     &       'Aout(idW_ao)',                                            &
     &       'Write out rate of melt/freeze at air/ocean interface.'
          IF (Aout(idW_io,ng)) WRITE (out,70) Aout(idW_io,ng),          &
     &       'Aout(idW_io)',                                            &
     &       'Write out rate of melt/freeze at ice/ocean interface.'
          IF (Aout(idW_ro,ng)) WRITE (out,70) Aout(idW_ro,ng),          &
     &       'Aout(idW_ro)',                                            &
     &       'Write out rate of melt/freeze runoff into ocean.'
#endif
        END DO
      END IF

#ifdef AVERAGES
!
!  Turn out switches to process time-averaged ice model state and
!  internal variables. It is necessary to avoid recursive dependency
!  between "mod_ncparam" and "mod_ice".
!
      LiceFavg=.FALSE.
      LiceSavg=.FALSE.
!
      DO ng=1,Ngrids
        DO nf=1,nIceF
          IF (iFice(nf).gt.0) THEN
            IF (Aout(iFice(nf),ng)) THEN
              LiceFavg(nf,ng)=.TRUE.
            END IF
          END IF
        END DO
!
        DO ns=1,nIceS
          IF (iSice(ns).gt.0) THEN
            IF (Aout(iSice(ns),ng)) THEN
              LiceSavg(ns,ng)=.TRUE.
            END IF
          END IF
        END DO
      END DO
#endif
!
  30  FORMAT (/,' READ_IcePar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Ice Parameters, Grid: ',i2.2,                       &
     &        /,  ' ========================',/)
  50  FORMAT (1x,i10,2x,a,t32,a)
  60  FORMAT (1p,e11.4,2x,a,t32,a)
  70  FORMAT (10x,l1,2x,a,t32,a)
  80  FORMAT (/,' READ_IcePar - variable index not yet loaded, ', a)
!
      RETURN
      END SUBROUTINE read_IcePar
