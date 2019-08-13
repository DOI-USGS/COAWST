#include "cppdefs.h"
#ifdef ICE_MODEL
      SUBROUTINE read_IcePar (inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports ice model input parameters.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, ng, status

      real(r8), dimension(nRval) :: Rval

      character (len=40) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
      logical, dimension(Ngrids) :: Lswitch
      real(r8) :: Rvalue(1)

!
!-----------------------------------------------------------------------
!  Read in ice model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lice')
              Npts=load_l(Nval, Cval, Ngrids, Lice)
            CASE ('nstrs')
              Npts=load_i(Nval, Rval, Ngrids, nstrs)
            CASE ('nevp')
              Npts=load_i(Nval, Rval, Ngrids, nevp)
            CASE ('rhoice')
              Npts=load_r(Nval, Rval, Ngrids, rhoice)
            CASE ('cdiw')
              Npts=load_r(Nval, Rval, Ngrids, cdiw)
            CASE ('cdai')
              Npts=load_r(Nval, Rval, Ngrids, cdai)
            CASE ('rho_air')
              Npts=load_r(Nval, Rval, Ngrids, rho_air)
            CASE ('rhosnow_dry')
              Npts=load_r(Nval, Rval, Ngrids, rhosnow_dry)
            CASE ('rhosnow_wet')
              Npts=load_r(Nval, Rval, Ngrids, rhosnow_wet)
            CASE ('pstar')
              Npts=load_r(Nval, Rval, Ngrids, pstar)
            CASE ('pstar_e')
              Npts=load_r(Nval, Rval, Ngrids, pstar_e)
            CASE ('astren')
              Npts=load_r(Nval, Rval, Ngrids, astren)
            CASE ('zetamax')
              Npts=load_r(Nval, Rval, Ngrids, zetamax)
            CASE ('zetamin')
              Npts=load_r(Nval, Rval, Ngrids, zetamin)
            CASE ('ellip_sq')
              Npts=load_r(Nval, Rval, Ngrids, ellip_sq)
            CASE ('alphai')
              Npts=load_r(Nval, Rval, Ngrids, alphai)
              do ng=1,Ngrids
                 alphai(ng) = alphai(ng)*deg2rad
              enddo
            CASE ('tol')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              tol=Rvalue(1)
            CASE ('min_h')
              Npts=load_r(Nval, Rval, Ngrids, min_h)
            CASE ('min_a')
              Npts=load_r(Nval, Rval, Ngrids, min_a)
            CASE ('max_a')
              Npts=load_r(Nval, Rval, Ngrids, max_a)
            CASE ('stressang')
              Npts=load_r(Nval, Rval, Ngrids, stressang)
              do ng=1,Ngrids
                 stressang(ng) = stressang(ng)*deg2rad
              enddo
            CASE ('ice_emiss')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              ice_emiss=Rvalue(1)
            CASE ('spec_heat_air')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              spec_heat_air=Rvalue(1)
            CASE ('trans_coeff')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              trans_coeff=Rvalue(1)
            CASE ('sublim_latent_heat')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              sublim_latent_heat=Rvalue(1)
            CASE ('t0deg')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              t0deg=Rvalue(1)
#ifdef ICE_LANDFAST
            CASE ('lf_k1')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              lf_k1=Rvalue(1)
            CASE ('lf_k2')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              lf_k2=Rvalue(1)
            CASE ('lf_u0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              lf_u0=Rvalue(1)
#endif
            CASE ('Hout(idUice)')
              IF (idUice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idVice)')
              IF (idVice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idVice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idUiceE)')
              IF (idUiceE.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUiceE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idUiceE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idViceN)')
              IF (idViceN.eq.0) THEN
                IF (Master) WRITE (out,80) 'idViceN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idViceN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idAice)')
              IF (idAice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHice)')
              IF (idHice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idHsno)')
              IF (idHsno.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHsno'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTice)')
              IF (idTice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTimid)')
              IF (idTimid.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTimid'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idTimid,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idT2ice)')
              IF (idT2ice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT2ice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idT2ice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idTauiw)')
              IF (idTauiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTauiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idTauiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idChuiw)')
              IF (idChuiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idChuiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idChuiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idAgeice)')
              IF (idAgeice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAgeice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAgeice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idIomflx)')
              IF (idIomflx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIomflx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idIomflx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSig11)')
              IF (idSig11.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig11'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSig11,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSig12)')
              IF (idSig12.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig12'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSig12,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSig22)')
              IF (idSig22.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig22'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSig22,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idT0mk)')
              IF (idT0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idS0mk)')
              IF (idS0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idS0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWfr)')
              IF (idWfr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWfr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWfr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWai)')
              IF (idWai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWao)')
              IF (idWao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWio)')
              IF (idWio.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWio'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWio,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWro)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWro'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWro,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idWdiv)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWdiv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
#ifdef ICE_DIAGS
            CASE ('Hout(idAlbe)')
              IF (idAlbe.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAlbe'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAlbe,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idAlbI)')
              IF (idAlbI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAlbI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idAlbI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idLWdn)')
              IF (idLWdn.eq.0) THEN
                IF (Master) WRITE (out,80) 'idLWdn'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idLWdn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSWdn)')
              IF (idSWdn.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWdn'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSWdn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idLatI)')
              IF (idLatI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idLatI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idLatI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSenI)')
              IF (idSenI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSenI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSenI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idLWupI)')
              IF (idLWupI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idLWupI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idLWupI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSWupI)')
              IF (idSWupI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWupI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSWupI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSWinI)')
              IF (idSWinI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWinI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSWinI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSWthruI)')
              IF (idSWthruI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWthruI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSWthruI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSaltI)')
              IF (idSaltI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSaltI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSaltI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSaltL)')
              IF (idSaltL.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSaltL'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSaltL,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idQao)')
              IF (idQao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idQao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idQai)')
              IF (idQai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idQai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idQio)')
              IF (idQio.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQio'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idQio,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idQi2)')
              IF (idQi2.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQi2'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idQi2,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Hout(idSnoI)')
              IF (idSnoI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSnoI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Hout(idSnoI,1:Ngrids)=Lswitch(1:Ngrids)
#endif
#ifdef AVERAGES
            CASE ('Aout(idUice)')
              IF (idUice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idVice)')
              IF (idVice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idVice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idUiceE)')
              IF (idUiceE.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUiceE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idUiceE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idViceN)')
              IF (idViceN.eq.0) THEN
                IF (Master) WRITE (out,80) 'idViceN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idViceN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idAice)')
              IF (idAice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idHice)')
              IF (idHice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idHsno)')
              IF (idHsno.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHsno'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idTice)')
              IF (idTice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idTimid)')
              IF (idTimid.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTimid'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idTimid,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idT2ice)')
              IF (idT2ice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT2ice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idT2ice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idTauiw)')
              IF (idTauiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTauiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idTauiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idChuiw)')
              IF (idChuiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idChuiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idChuiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idAgeice)')
              IF (idAgeice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAgeice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAgeice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idIomflx)')
              IF (idIomflx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIomflx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idIomflx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSig11)')
              IF (idSig11.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig11'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSig11,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSig12)')
              IF (idSig12.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig12'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSig12,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSig22)')
              IF (idSig22.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig22'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSig22,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idT0mk)')
              IF (idT0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idS0mk)')
              IF (idS0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idS0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWfr)')
              IF (idWfr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWfr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWfr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWai)')
              IF (idWai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWao)')
              IF (idWao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWio)')
              IF (idWio.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWio'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWio,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWro)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWro'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWro,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idWdiv)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWdiv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
# ifdef ICE_DIAGS
            CASE ('Aout(idAlbe)')
              IF (idAlbe.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAlbe'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAlbe,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idAlbI)')
              IF (idAlbI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAlbI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idAlbI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idLWdn)')
              IF (idLWdn.eq.0) THEN
                IF (Master) WRITE (out,80) 'idLWdn'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idLWdn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSWdn)')
              IF (idSWdn.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWdn'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSWdn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idLatI)')
              IF (idLatI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idLatI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idLatI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSenI)')
              IF (idSenI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSenI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSenI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idLWupI)')
              IF (idLWupI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idLWupI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idLWupI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSWupI)')
              IF (idSWupI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWupI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSWupI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSWinI)')
              IF (idSWinI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWinI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSWinI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSWthruI)')
              IF (idSWthruI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSWthruI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSWthruI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSaltI)')
              IF (idSaltI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSaltI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSaltI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSaltL)')
              IF (idSaltL.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSaltL'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSaltL,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idQao)')
              IF (idQao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idQao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idQai)')
              IF (idQai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idQai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idQio)')
              IF (idQio.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQio'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idQio,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idQi2)')
              IF (idQi2.eq.0) THEN
                IF (Master) WRITE (out,80) 'idQi2'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idQi2,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout(idSnoI)')
              IF (idSnoI.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSnoI'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout(idSnoI,1:Ngrids)=Lswitch(1:Ngrids)
# endif
#endif
#ifdef AVERAGES2
            CASE ('Aout2(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idUiceE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idUiceE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idViceN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idViceN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idTimid)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idTimid,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idT2ice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idT2ice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idTauiw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idTauiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idChuiw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idChuiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idAgeice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idAgeice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idIomflx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idIomflx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSig11)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSig11,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSig12)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSig12,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSig22)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSig22,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idWfr)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idWfr,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idWai)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idWai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idWao)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idWao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idWio)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idWio,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idWro)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idWro,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idWdiv)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idWdiv,1:Ngrids)=Lswitch(1:Ngrids)
# ifdef ICE_DIAGS
            CASE ('Aout2(idAlbe)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idAlbe,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idAlbI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idAlbI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idLWdn)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idLWdn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSWdn)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSWdn,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idLatI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idLatI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSenI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSenI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idLWupI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idLWupI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSWupI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSWupI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSWinI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSWinI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSWthruI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSWthruI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSaltI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSaltI,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSaltL)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSaltL,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idQao)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idQao,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idQai)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idQai,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idQio)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idQio,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idQi2)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idQi2,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Aout2(idSnoI)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Aout2(idSnoI,1:Ngrids)=Lswitch(1:Ngrids)
# endif
#endif
            CASE ('Qout(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idVice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idUiceE)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idUiceE,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idViceN)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idViceN,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idAice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idHsno,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idTice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTimid)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idTimid,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idT2ice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idT2ice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idTauiw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idTauiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idChuiw)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idChuiw,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idAgeice)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idAgeice,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idIomflx)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idIomflx,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSig11)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSig11,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSig12)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSig12,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idSig22)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idSig22,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idT0mk,1:Ngrids)=Lswitch(1:Ngrids)
            CASE ('Qout(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              Qout(idS0mk,1:Ngrids)=Lswitch(1:Ngrids)
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
          IF (.not. Lice(ng)) RETURN
          WRITE (out,40) ng
          WRITE(out,*) 'Ice time step = ocean time step'
          WRITE (out,60) dtice(ng), 'DTICE',                            &
     &          'Ice model time step (s).'
          WRITE (out,50) nstrs(ng), 'nstrs',                            &
     &          'Number of iterations for nonlinear ice dynamics.'
          WRITE (out,50) nevp(ng), 'nevp',                              &
     &          'Number of elastic steps per plastic step in EVP.'
          WRITE (out,60) rhoice(ng), 'rhoice',                          &
     &          'Density of sea ice (kg/m3).'
          WRITE (out,60) cdiw(ng), 'cdiw',                              &
     &          'Ice-water drag coefficient (nondimensional).'
          WRITE (out,60) cdai(ng), 'cdai',                              &
     &          'Air-ice drag coefficient (nondimensional).'
          WRITE (out,60) rho_air(ng), 'rho_air',                        &
     &          'Air density (kg/m3).'
          WRITE (out,60) rhosnow_dry(ng), 'rhosnow_dry',                &
     &          'Dry snow density (kg/m3).'
          WRITE (out,60) rhosnow_wet(ng), 'rhosnow_wet',                &
     &          'Wet snow density (kg/m3).'
          WRITE (out,60) alphai(ng)*rad2deg, 'alphai',                  &
     &          'Mohr-Coulomb stress angle (degrees).'
          WRITE (out,60) min_h(ng), 'min_h',                            &
     &          'Minimum average ice thickness (m).'
          WRITE (out,60) min_a(ng), 'min_a',                            &
     &          'Minimum ice concentration (nondimensional).'
          WRITE (out,60) max_a(ng), 'max_a',                            &
     &          'Maximum ice concentration (nondimensional).'
          WRITE (out,60) stressang(ng)*rad2deg, 'stressang',            &
     &          'Turning angle for ice-water drag (degrees).'
          WRITE (out,60) tol, 'tol',                                    &
     &          'Numerical tolerance in rheology calculations .'
          WRITE (out,60) ice_emiss, 'ice_emiss',                        &
     &          'Ice emissivity.'
          WRITE (out,60) spec_heat_air, 'spec_heat_air',                &
     &          'Specific heat of air.'
          WRITE (out,60) trans_coeff, 'trans_coeff',                    &
     &          'Transfer coefficient.'
          WRITE (out,60) sublim_latent_heat, 'sublim_latent_heat',      &
     &          'Latent_heat of sublimation.'
          WRITE (out,60) t0deg, 't0deg',                                &
     &          'Zero degrees Celsius in degrees Kelvin.'
#ifdef ICE_LANDFAST
          WRITE (out,60) lf_k1, 'lf_k1',                                &
     &          'Landfast ice tunable k1.'
          WRITE (out,60) lf_k2, 'lf_k2',                                &
     &          'Landfast ice tunable k2.'
          WRITE (out,60) lf_u0, 'lf_u0',                                &
     &          'Landfast ice minimal velocity.'
#endif
          IF (Hout(idUice,ng)) WRITE (out,70) Hout(idUice,ng),          &
     &       'Hout(idUice)',                                            &
     &       'Write out U-component ice velocity.'
          IF (Hout(idVice,ng)) WRITE (out,70) Hout(idVice,ng),          &
     &       'Hout(idVice)',                                            &
     &       'Write out V-component ice velocity.'
          IF (Hout(idUiceE,ng)) WRITE (out,70) Hout(idUiceE,ng),        &
     &       'Hout(idUiceE)',                                           &
     &       'Write out East component ice velocity.'
          IF (Hout(idViceN,ng)) WRITE (out,70) Hout(idViceN,ng),        &
     &       'Hout(idViceN)',                                           &
     &       'Write out North component ice velocity.'
          IF (Hout(idAice,ng)) WRITE (out,70) Hout(idAice,ng),          &
     &       'Hout(idAice)',                                            &
     &       'Write out ice concentration.'
          IF (Hout(idHice,ng)) WRITE (out,70) Hout(idHice,ng),          &
     &       'Hout(idHice)',                                            &
     &       'Write out average ice thickness.'
          IF (Hout(idHsno,ng)) WRITE (out,70) Hout(idHsno,ng),          &
     &       'Hout(idHsno)',                                            &
     &       'Write out snow thickness.'
          IF (Hout(idTice,ng)) WRITE (out,70) Hout(idTice,ng),          &
     &       'Hout(idTice)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Hout(idTimid,ng)) WRITE (out,70) Hout(idTimid,ng),        &
     &       'Hout(idTimid)',                                           &
     &       'Write out interior ice temperature.'
          IF (Hout(idT2ice,ng)) WRITE (out,70) Hout(idT2ice,ng),        &
     &       'Hout(idT2ice)',                                           &
     &       'Write out top of ice temperature.'
          IF (Hout(idAgeice,ng)) WRITE (out,70) Hout(idAgeice,ng),      &
     &       'Hout(idAgeice)',                                          &
     &       'Write out ice age.'
          IF (Hout(idIomflx,ng)) WRITE (out,70) Hout(idIomflx,ng),      &
     &       'Hout(idIomflx)',                                          &
     &       'Write out ice-ocean mass flux'
          IF (Hout(idSig11,ng)) WRITE (out,70) Hout(idSig11,ng),        &
     &       'Hout(idSig11)',                                           &
     &       'Write out internal ice stress component 11.'
          IF (Hout(idSig12,ng)) WRITE (out,70) Hout(idSig12,ng),        &
     &       'Hout(idSig12)',                                           &
     &       'Write out internal ice stress component 12.'
          IF (Hout(idSig22,ng)) WRITE (out,70) Hout(idSig22,ng),        &
     &       'Hout(idSig22)',                                           &
     &       'Write out internal ice stress component 22.'
          IF (Hout(idTauiw,ng)) WRITE (out,70) Hout(idTauiw,ng),        &
     &       'Hout(idTauiw)',                                           &
     &       'Write out ice-water friction velocity.'
          IF (Hout(idChuiw,ng)) WRITE (out,70) Hout(idChuiw,ng),        &
     &       'Hout(idChuiw)',                                           &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Hout(idT0mk,ng)) WRITE (out,70) Hout(idT0mk,ng),          &
     &       'Hout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Hout(idS0mk,ng)) WRITE (out,70) Hout(idS0mk,ng),          &
     &       'Hout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
          IF (Hout(idWfr,ng)) WRITE (out,70) Hout(idWfr,ng),            &
     &       'Hout(idWfr)',                                             &
     &       'Write out frazil ice growth rate.'
          IF (Hout(idWai,ng)) WRITE (out,70) Hout(idWai,ng),            &
     &       'Hout(idWai)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Hout(idWao,ng)) WRITE (out,70) Hout(idWao,ng),            &
     &       'Hout(idWao)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Hout(idWio,ng)) WRITE (out,70) Hout(idWio,ng),            &
     &       'Hout(idWio)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Hout(idWro,ng)) WRITE (out,70) Hout(idWro,ng),            &
     &       'Hout(idWro)',                                             &
     &       'Write out ice melt runoff rate.'
          IF (Hout(idWdiv,ng)) WRITE (out,70) Hout(idWdiv,ng),          &
     &       'Hout(idWdiv)',                                            &
     &       'Write out ice divergence rate.'
#ifdef ICE_DIAGS
          IF (Hout(idAlbe,ng)) WRITE (out,70) Hout(idAlbe,ng),          &
     &       'Hout(idAlbe)',                                            &
     &       'Write out ocean surface albedo.'
          IF (Hout(idAlbI,ng)) WRITE (out,70) Hout(idAlbI,ng),          &
     &       'Hout(idAlbI)',                                            &
     &       'Write out ice albedo.'
          IF (Hout(idLWdn,ng)) WRITE (out,70) Hout(idLWdn,ng),          &
     &       'Hout(idLWdn)',                                            &
     &       'Write out downward longwave radiation.'
          IF (Hout(idSWdn,ng)) WRITE (out,70) Hout(idSWdn,ng),          &
     &       'Hout(idSWdn)',                                            &
     &       'Write out downward shortwave radiation.'
          IF (Hout(idLatI,ng)) WRITE (out,70) Hout(idLatI,ng),          &
     &       'Hout(idLatI)',                                            &
     &       'Write out latent heat flux over ice.'
          IF (Hout(idSenI,ng)) WRITE (out,70) Hout(idSenI,ng),          &
     &       'Hout(idSenI)',                                            &
     &       'Write out sensible heat flux over ice.'
          IF (Hout(idLWupI,ng)) WRITE (out,70) Hout(idLWupI,ng),        &
     &       'Hout(idLWupI)',                                           &
     &       'Write out upward longwave over ice.'
          IF (Hout(idSWupI,ng)) WRITE (out,70) Hout(idSWupI,ng),        &
     &       'Hout(idSWupI)',                                           &
     &       'Write out upward shortwave over ice.'
          IF (Hout(idSWinI,ng)) WRITE (out,70) Hout(idSWinI,ng),        &
     &       'Hout(idSWinI)',                                           &
     &       'Write out upward shortwave over ice.'
          IF (Hout(idSWthruI,ng)) WRITE (out,70) Hout(idSWthruI,ng),    &
     &       'Hout(idSWthruI)',                                         &
     &       'Write out shortwave radiation through ice.'
          IF (Hout(idSaltI,ng)) WRITE (out,70) Hout(idSaltI,ng),        &
     &       'Hout(idSaltI)',                                           &
     &       'Write out salt flux under ice.'
          IF (Hout(idSaltL,ng)) WRITE (out,70) Hout(idSaltL,ng),        &
     &       'Hout(idSaltL)',                                           &
     &       'Write out salt flux in leads.'
          IF (Hout(idQao,ng)) WRITE (out,70) Hout(idQao,ng),            &
     &       'Hout(idQao)',                                             &
     &       'Write out heat flux in leads.'
          IF (Hout(idQai,ng)) WRITE (out,70) Hout(idQai,ng),            &
     &       'Hout(idQai)',                                             &
     &       'Write out heat flux on top of ice.'
          IF (Hout(idQio,ng)) WRITE (out,70) Hout(idQio,ng),            &
     &       'Hout(idQio)',                                             &
     &       'Write out heat flux into ice.'
          IF (Hout(idQi2,ng)) WRITE (out,70) Hout(idQi2,ng),            &
     &       'Hout(idQi2)',                                             &
     &       'Write out heat flux in ice.'
          IF (Hout(idSnoI,ng)) WRITE (out,70) Hout(idSnoI,ng),          &
     &       'Hout(idSnoI)',                                            &
     &       'Write out heat flux into ice.'
#endif
#ifdef AVERAGES
          IF (Aout(idUice,ng)) WRITE (out,70) Aout(idUice,ng),          &
     &       'Aout(idUice)',                                            &
     &       'Write out time-averaged U-component ice velocity.'
          IF (Aout(idVice,ng)) WRITE (out,70) Aout(idVice,ng),          &
     &       'Aout(idVice)',                                            &
     &       'Write out time-averaged V-component ice velocity.'
          IF (Aout(idUiceE,ng)) WRITE (out,70) Aout(idUiceE,ng),        &
     &       'Aout(idUiceE)',                                           &
     &       'Write out time-averaged East component ice velocity.'
          IF (Aout(idViceN,ng)) WRITE (out,70) Aout(idViceN,ng),        &
     &       'Aout(idViceN)',                                           &
     &       'Write out time-averaged North component ice velocity.'
          IF (Aout(idAice,ng)) WRITE (out,70) Aout(idAice,ng),          &
     &       'Aout(idAice)',                                            &
     &       'Write out time-averaged ice concentration.'
          IF (Aout(idHice,ng)) WRITE (out,70) Aout(idHice,ng),          &
     &       'Aout(idHice)',                                            &
     &       'Write out time-averaged average ice thickness.'
          IF (Aout(idHsno,ng)) WRITE (out,70) Aout(idHsno,ng),          &
     &       'Aout(idHsno)',                                            &
     &       'Write out time-averaged snow thickness.'
          IF (Aout(idTice,ng)) WRITE (out,70) Aout(idTice,ng),          &
     &       'Aout(idTice)',                                            &
     &       'Write out time-averaged ice/snow surface temperature.'
          IF (Aout(idTimid,ng)) WRITE (out,70) Aout(idTimid,ng),        &
     &       'Aout(idTimid)',                                           &
     &       'Write out time-averaged interior ice temperature.'
          IF (Aout(idT2ice,ng)) WRITE (out,70) Aout(idT2ice,ng),        &
     &       'Aout(idT2ice)',                                           &
     &       'Write out time-averaged top of ice temperature.'
          IF (Aout(idAgeice,ng)) WRITE (out,70) Aout(idAgeice,ng),      &
     &       'Aout(idAgeice)',                                          &
     &       'Write out time-averaged ice age.'
          IF (Aout(idIomflx,ng)) WRITE (out,70) Aout(idIomflx,ng),      &
     &       'Aout(idIomflx)',                                          &
     &       'Write out time-averaged ice-ocean mass flux'
          IF (Aout(idSig11,ng)) WRITE (out,70) Aout(idSig11,ng),        &
     &       'Aout(idSig11)',                                           &
     &       'Write out time-averaged internal ice stress component 11.'
          IF (Aout(idSig12,ng)) WRITE (out,70) Aout(idSig12,ng),        &
     &       'Aout(idSig12)',                                           &
     &       'Write out time-averaged internal ice stress component 12.'
          IF (Aout(idSig22,ng)) WRITE (out,70) Aout(idSig22,ng),        &
     &       'Aout(idSig22)',                                           &
     &       'Write out time-averaged internal ice stress component 22.'
          IF (Aout(idTauiw,ng)) WRITE (out,70) Aout(idTauiw,ng),        &
     &       'Aout(idTauiw)',                                           &
     &       'Write out time-averaged ice-water friction velocity.'
          IF (Aout(idChuiw,ng)) WRITE (out,70) Aout(idChuiw,ng),        &
     &       'Aout(idChuiw)',                                           &
     &       'Write out time-averaged ice-water transfer coefficient.'
          IF (Aout(idT0mk,ng)) WRITE (out,70) Aout(idT0mk,ng),          &
     &       'Aout(idT0mk)',                                            &
     &       'Write out time-averaged under ice temperature.'
          IF (Aout(idS0mk,ng)) WRITE (out,70) Aout(idS0mk,ng),          &
     &       'Aout(idS0mk)',                                            &
     &       'Write out time-averaged under ice salinity.'
          IF (Aout(idWfr,ng)) WRITE (out,70) Aout(idWfr,ng),            &
     &       'Aout(idWfr)',                                             &
     &       'Write out time-averaged frazil ice growth rate.'
          IF (Aout(idWai,ng)) WRITE (out,70) Aout(idWai,ng),            &
     &       'Aout(idWai)',                                             &
     &       'Write out time-averaged ice growth/melt rate.'
          IF (Aout(idWao,ng)) WRITE (out,70) Aout(idWao,ng),            &
     &       'Aout(idWao)',                                             &
     &       'Write out time-averaged ice growth/melt rate.'
          IF (Aout(idWio,ng)) WRITE (out,70) Aout(idWio,ng),            &
     &       'Aout(idWio)',                                             &
     &       'Write out time-averaged ice growth/melt rate.'
          IF (Aout(idWro,ng)) WRITE (out,70) Aout(idWro,ng),            &
     &       'Aout(idWro)',                                             &
     &       'Write out time-averaged ice melt runoff rate.'
          IF (Aout(idWdiv,ng)) WRITE (out,70) Aout(idWdiv,ng),          &
     &       'Aout(idWdiv)',                                            &
     &       'Write out time-averaged ice divergence rate.'
# ifdef ICE_DIAGS
          IF (Aout(idAlbe,ng)) WRITE (out,70) Aout(idAlbe,ng),          &
     &       'Aout(idAlbe)',                                            &
     &       'Write out ocean surface albedo.'
          IF (Aout(idAlbI,ng)) WRITE (out,70) Aout(idAlbI,ng),          &
     &       'Aout(idAlbI)',                                            &
     &       'Write out ice albedo.'
          IF (Aout(idLWdn,ng)) WRITE (out,70) Aout(idLWdn,ng),          &
     &       'Aout(idLWdn)',                                            &
     &       'Write out downward longwave radiation.'
          IF (Aout(idSWdn,ng)) WRITE (out,70) Aout(idSWdn,ng),          &
     &       'Aout(idSWdn)',                                            &
     &       'Write out downward shortwave radiation.'
          IF (Aout(idLatI,ng)) WRITE (out,70) Aout(idLatI,ng),          &
     &       'Aout(idLatI)',                                            &
     &       'Write out latent heat flux over ice.'
          IF (Aout(idSenI,ng)) WRITE (out,70) Aout(idSenI,ng),          &
     &       'Aout(idSenI)',                                            &
     &       'Write out sensible heat flux over ice.'
          IF (Aout(idLWupI,ng)) WRITE (out,70) Aout(idLWupI,ng),        &
     &       'Aout(idLWupI)',                                           &
     &       'Write out upward longwave over ice.'
          IF (Aout(idSWupI,ng)) WRITE (out,70) Aout(idSWupI,ng),        &
     &       'Aout(idSWupI)',                                           &
     &       'Write out upward shortwave over ice.'
          IF (Aout(idSWinI,ng)) WRITE (out,70) Aout(idSWinI,ng),        &
     &       'Aout(idSWinI)',                                           &
     &       'Write out upward shortwave over ice.'
          IF (Aout(idSWthruI,ng)) WRITE (out,70) Aout(idSWthruI,ng),    &
     &       'Aout(idSWthruI)',                                         &
     &       'Write out shortwave radiation through ice.'
          IF (Aout(idSaltI,ng)) WRITE (out,70) Aout(idSaltI,ng),        &
     &       'Aout(idSaltI)',                                           &
     &       'Write out salt flux under ice.'
          IF (Aout(idSaltL,ng)) WRITE (out,70) Aout(idSaltL,ng),        &
     &       'Aout(idSaltL)',                                           &
     &       'Write out salt flux in leads.'
          IF (Aout(idQao,ng)) WRITE (out,70) Aout(idQao,ng),            &
     &       'Aout(idQao)',                                             &
     &       'Write out heat flux in leads.'
          IF (Aout(idQai,ng)) WRITE (out,70) Aout(idQai,ng),            &
     &       'Aout(idQai)',                                             &
     &       'Write out heat flux on top of ice.'
          IF (Aout(idQio,ng)) WRITE (out,70) Aout(idQio,ng),            &
     &       'Aout(idQio)',                                             &
     &       'Write out heat flux into ice.'
          IF (Aout(idQi2,ng)) WRITE (out,70) Aout(idQi2,ng),            &
     &       'Aout(idQi2)',                                             &
     &       'Write out heat flux in ice.'
          IF (Aout(idSnoI,ng)) WRITE (out,70) Aout(idSnoI,ng),          &
     &       'Aout(idSnoI)',                                            &
     &       'Write out heat flux into ice.'
# endif
#endif
#ifdef AVERAGES2
          IF (Aout2(idUice,ng)) WRITE (out,70) Aout2(idUice,ng),        &
     &       'Aout2(idUice)',                                           &
     &       'Write out U-component ice velocity.'
          IF (Aout2(idVice,ng)) WRITE (out,70) Aout2(idVice,ng),        &
     &       'Aout2(idVice)',                                           &
     &       'Write out V-component ice velocity.'
          IF (Aout2(idUiceE,ng)) WRITE (out,70) Aout2(idUiceE,ng),      &
     &       'Aout2(idUiceE)',                                          &
     &       'Write out East component ice velocity.'
          IF (Aout2(idViceN,ng)) WRITE (out,70) Aout2(idViceN,ng),      &
     &       'Aout2(idViceN)',                                          &
     &       'Write out North component ice velocity.'
          IF (Aout2(idAice,ng)) WRITE (out,70) Aout2(idAice,ng),        &
     &       'Aout2(idAice)',                                           &
     &       'Write out ice concentration.'
          IF (Aout2(idHice,ng)) WRITE (out,70) Aout2(idHice,ng),        &
     &       'Aout2(idHice)',                                           &
     &       'Write out average ice thickness.'
          IF (Aout2(idHsno,ng)) WRITE (out,70) Aout2(idHsno,ng),        &
     &       'Aout2(idHsno)',                                           &
     &       'Write out snow thickness.'
          IF (Aout2(idTice,ng)) WRITE (out,70) Aout2(idTice,ng),        &
     &       'Aout2(idTice)',                                           &
     &       'Write out ice/snow surface temperature.'
          IF (Aout2(idTimid,ng)) WRITE (out,70) Aout2(idTimid,ng),      &
     &       'Aout2(idTimid)',                                          &
     &       'Write out interior ice temperature.'
          IF (Aout2(idT2ice,ng)) WRITE (out,70) Aout2(idT2ice,ng),      &
     &       'Aout2(idT2ice)',                                          &
     &       'Write out top of ice temperature.'
          IF (Aout2(idAgeice,ng)) WRITE (out,70) Aout2(idAgeice,ng),    &
     &       'Aout2(idAgeice)',                                         &
     &       'Write out ice age.'
          IF (Aout2(idIomflx,ng)) WRITE (out,70) Aout2(idIomflx,ng),    &
     &       'Aout2(idIomflx)',                                         &
     &       'Write out ice-ocean mass flux'
          IF (Aout2(idSig11,ng)) WRITE (out,70) Aout2(idSig11,ng),      &
     &       'Aout2(idSig11)',                                          &
     &       'Write out internal ice stress component 11.'
          IF (Aout2(idSig12,ng)) WRITE (out,70) Aout2(idSig12,ng),      &
     &       'Aout2(idSig12)',                                          &
     &       'Write out internal ice stress component 12.'
          IF (Aout2(idSig22,ng)) WRITE (out,70) Aout2(idSig22,ng),      &
     &       'Aout2(idSig22)',                                          &
     &       'Write out internal ice stress component 22.'
          IF (Aout2(idTauiw,ng)) WRITE (out,70) Aout2(idTauiw,ng),      &
     &       'Aout2(idTauiw)',                                          &
     &       'Write out ice-water friction velocity.'
          IF (Aout2(idChuiw,ng)) WRITE (out,70) Aout2(idChuiw,ng),      &
     &       'Aout2(idChuiw)',                                          &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Aout2(idT0mk,ng)) WRITE (out,70) Aout2(idT0mk,ng),        &
     &       'Aout2(idT0mk)',                                           &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Aout2(idS0mk,ng)) WRITE (out,70) Aout2(idS0mk,ng),        &
     &       'Aout2(idS0mk)',                                           &
     &       'Write out salinity of molecular sublayer under ice.'
#endif
          IF (Qout(idUice,ng)) WRITE (out,70) Qout(idUice,ng),          &
     &       'Qout(idUice)',                                            &
     &       'Write out U-component ice velocity.'
          IF (Qout(idVice,ng)) WRITE (out,70) Qout(idVice,ng),          &
     &       'Qout(idVice)',                                            &
     &       'Write out V-component ice velocity.'
          IF (Qout(idUiceE,ng)) WRITE (out,70) Qout(idUiceE,ng),        &
     &       'Qout(idUiceE)',                                           &
     &       'Write out East component ice velocity.'
          IF (Qout(idViceN,ng)) WRITE (out,70) Qout(idViceN,ng),        &
     &       'Qout(idViceN)',                                           &
     &       'Write out North component ice velocity.'
          IF (Qout(idAice,ng)) WRITE (out,70) Qout(idAice,ng),          &
     &       'Qout(idAice)',                                            &
     &       'Write out ice concentration.'
          IF (Qout(idHice,ng)) WRITE (out,70) Qout(idHice,ng),          &
     &       'Qout(idHice)',                                            &
     &       'Write out average ice thickness.'
          IF (Qout(idHsno,ng)) WRITE (out,70) Qout(idHsno,ng),          &
     &       'Qout(idHsno)',                                            &
     &       'Write out snow thickness.'
          IF (Qout(idTice,ng)) WRITE (out,70) Qout(idTice,ng),          &
     &       'Qout(idTice)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Qout(idTimid,ng)) WRITE (out,70) Qout(idTimid,ng),        &
     &       'Qout(idTimid)',                                           &
     &       'Write out interior ice temperature.'
          IF (Qout(idT2ice,ng)) WRITE (out,70) Qout(idT2ice,ng),        &
     &       'Qout(idT2ice)',                                           &
     &       'Write out top of ice temperature.'
          IF (Qout(idAgeice,ng)) WRITE (out,70) Qout(idAgeice,ng),      &
     &       'Qout(idAgeice)',                                          &
     &       'Write out ice age.'
          IF (Qout(idIomflx,ng)) WRITE (out,70) Qout(idIomflx,ng),      &
     &       'Qout(idIomflx)',                                          &
     &       'Write out ice-ocean mass flux'
          IF (Qout(idSig11,ng)) WRITE (out,70) Qout(idSig11,ng),        &
     &       'Qout(idSig11)',                                           &
     &       'Write out internal ice stress component 11.'
          IF (Qout(idSig12,ng)) WRITE (out,70) Qout(idSig12,ng),        &
     &       'Qout(idSig12)',                                           &
     &       'Write out internal ice stress component 12.'
          IF (Qout(idSig22,ng)) WRITE (out,70) Qout(idSig22,ng),        &
     &       'Qout(idSig22)',                                           &
     &       'Write out internal ice stress component 22.'
          IF (Qout(idTauiw,ng)) WRITE (out,70) Qout(idTauiw,ng),        &
     &       'Qout(idTauiw)',                                           &
     &       'Write out ice-water friction velocity.'
          IF (Qout(idChuiw,ng)) WRITE (out,70) Qout(idChuiw,ng),        &
     &       'Qout(idChuiw)',                                           &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Qout(idT0mk,ng)) WRITE (out,70) Qout(idT0mk,ng),          &
     &       'Qout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Qout(idS0mk,ng)) WRITE (out,70) Qout(idS0mk,ng),          &
     &       'Qout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
        END DO
      END IF

  30  FORMAT (/,' READ_IcePar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Ice Parameters, Grid: ',i2.2,                   &
     &        /,' ============================',/)
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (1p,e11.4,2x,a,t28,a)
  70  FORMAT (10x,l1,2x,a,t30,a)
  80  FORMAT (/,' READ_IcePar - variable info not yet loaded, ', a)
      RETURN
      END SUBROUTINE read_IcePar
#endif
