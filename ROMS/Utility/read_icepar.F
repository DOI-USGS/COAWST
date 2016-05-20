#include "cppdefs.h"
#ifdef ICE_MODEL
      SUBROUTINE read_IcePar (inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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

      integer :: decode_line, load_i, load_l, load_r

      real(r8), dimension(200) :: Rval

      character (len=40) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval

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
              Npts=load_l(Nval, Cval, Ngrids, Lfloats)
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
              Npts=load_r(Nval, Rval, 1, tol)
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
              Npts=load_r(Nval, Rval, 1, ice_emiss)
            CASE ('spec_heat_air')
              Npts=load_r(Nval, Rval, 1, spec_heat_air)
            CASE ('trans_coeff')
              Npts=load_r(Nval, Rval, 1, trans_coeff)
            CASE ('sublim_latent_heat')
              Npts=load_r(Nval, Rval, 1, sublim_latent_heat)
            CASE ('t0deg')
              Npts=load_r(Nval, Rval, 1, t0deg)
#ifdef ICE_LANDFAST
            CASE ('lf_k1')
              Npts=load_r(Nval, Rval, 1, lf_k1)
            CASE ('lf_k2')
              Npts=load_r(Nval, Rval, 1, lf_k2)
            CASE ('lf_u0')
              Npts=load_r(Nval, Rval, 1, lf_u0)
#endif
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
            CASE ('Hout(idUiceE)')
              IF (idUiceE.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUiceE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUiceE,:))
            CASE ('Hout(idViceN)')
              IF (idViceN.eq.0) THEN
                IF (Master) WRITE (out,80) 'idViceN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idViceN,:))
            CASE ('Hout(idAice)')
              IF (idAice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idAice,:))
            CASE ('Hout(idHice)')
              IF (idHice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHice,:))
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
            CASE ('Hout(idTimid)')
              IF (idTimid.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTimid'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTimid,:))
            CASE ('Hout(idTauiw)')
              IF (idTauiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTauiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTauiw,:))
            CASE ('Hout(idChuiw)')
              IF (idChuiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idChuiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idChuiw,:))
#ifdef MELT_PONDS
            CASE ('Hout(idApond)')
              IF (idApond.eq.0) THEN
                IF (Master) WRITE (out,80) 'idApond'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idApond,:))
            CASE ('Hout(idHpond)')
              IF (idHpond.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHpond'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHpond,:))
#endif
            CASE ('Hout(idAgeice)')
              IF (idAgeice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAgeice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idAgeice,:))
            CASE ('Hout(idIomflx)')
              IF (idIomflx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIomflx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idIomflx,:))
            CASE ('Hout(idSig11)')
              IF (idSig11.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig11'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSig11,:))
            CASE ('Hout(idSig12)')
              IF (idSig12.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig12'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSig12,:))
            CASE ('Hout(idSig22)')
              IF (idSig22.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig22'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSig22,:))
            CASE ('Hout(idT0mk)')
              IF (idT0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idT0mk,:))
            CASE ('Hout(idS0mk)')
              IF (idS0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idS0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idS0mk,:))
            CASE ('Hout(idWfr)')
              IF (idWfr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWfr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWfr,:))
            CASE ('Hout(idWai)')
              IF (idWai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWai,:))
            CASE ('Hout(idWao)')
              IF (idWao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWao,:))
            CASE ('Hout(idWio)')
              IF (idWio.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWio'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWio,:))
            CASE ('Hout(idWro)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWro'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWro,:))
            CASE ('Hout(idWdiv)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWdiv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdiv,:))
#ifdef AVERAGES
            CASE ('Aout(idUice)')
              IF (idUice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUice,:))
            CASE ('Aout(idVice)')
              IF (idVice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idVice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVice,:))
            CASE ('Aout(idUiceE)')
              IF (idUiceE.eq.0) THEN
                IF (Master) WRITE (out,80) 'idUiceE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUiceE,:))
            CASE ('Aout(idViceN)')
              IF (idViceN.eq.0) THEN
                IF (Master) WRITE (out,80) 'idViceN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idViceN,:))
            CASE ('Aout(idAice)')
              IF (idAice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idAice,:))
            CASE ('Aout(idHice)')
              IF (idHice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHice,:))
            CASE ('Aout(idHsno)')
              IF (idHsno.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHsno'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHsno,:))
            CASE ('Aout(idTice)')
              IF (idTice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idTice,:))
            CASE ('Aout(idTimid)')
              IF (idTimid.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTimid'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idTimid,:))
            CASE ('Aout(idTauiw)')
              IF (idTauiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idTauiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idTauiw,:))
            CASE ('Aout(idChuiw)')
              IF (idChuiw.eq.0) THEN
                IF (Master) WRITE (out,80) 'idChuiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idChuiw,:))
#ifdef MELT_PONDS
            CASE ('Aout(idApond)')
              IF (idApond.eq.0) THEN
                IF (Master) WRITE (out,80) 'idApond'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idApond,:))
            CASE ('Aout(idHpond)')
              IF (idHpond.eq.0) THEN
                IF (Master) WRITE (out,80) 'idHpond'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHpond,:))
#endif
            CASE ('Aout(idAgeice)')
              IF (idAgeice.eq.0) THEN
                IF (Master) WRITE (out,80) 'idAgeice'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idAgeice,:))
            CASE ('Aout(idIomflx)')
              IF (idIomflx.eq.0) THEN
                IF (Master) WRITE (out,80) 'idIomflx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idIomflx,:))
            CASE ('Aout(idSig11)')
              IF (idSig11.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig11'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSig11,:))
            CASE ('Aout(idSig12)')
              IF (idSig12.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig12'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSig12,:))
            CASE ('Aout(idSig22)')
              IF (idSig22.eq.0) THEN
                IF (Master) WRITE (out,80) 'idSig22'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSig22,:))
            CASE ('Aout(idT0mk)')
              IF (idT0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idT0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idT0mk,:))
            CASE ('Aout(idS0mk)')
              IF (idS0mk.eq.0) THEN
                IF (Master) WRITE (out,80) 'idS0mk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idS0mk,:))
            CASE ('Aout(idWfr)')
              IF (idWfr.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWfr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWfr,:))
            CASE ('Aout(idWai)')
              IF (idWai.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWai'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWai,:))
            CASE ('Aout(idWao)')
              IF (idWao.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWao'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWao,:))
            CASE ('Aout(idWio)')
              IF (idWio.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWio'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWio,:))
            CASE ('Aout(idWro)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWro'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWro,:))
            CASE ('Aout(idWdiv)')
              IF (idWro.eq.0) THEN
                IF (Master) WRITE (out,80) 'idWdiv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdiv,:))
#endif
#ifdef AVERAGES2
            CASE ('Aout2(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUice,:))
            CASE ('Aout2(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVice,:))
            CASE ('Aout2(idUiceE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUiceE,:))
            CASE ('Aout2(idViceN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idViceN,:))
            CASE ('Aout2(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idAice,:))
            CASE ('Aout2(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHice,:))
            CASE ('Aout2(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHsno,:))
            CASE ('Aout2(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idTice,:))
            CASE ('Aout2(idTimid)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idTimid,:))
            CASE ('Aout2(idTauiw)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idTauiw,:))
            CASE ('Aout2(idChuiw)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idChuiw,:))
#ifdef MELT_PONDS
            CASE ('Aout2(idApond)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idApond,:))
            CASE ('Aout2(idHpond)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHpond,:))
#endif
            CASE ('Aout2(idAgeice)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idAgeice,:))
            CASE ('Aout2(idIomflx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idIomflx,:))
            CASE ('Aout2(idSig11)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idSig11,:))
            CASE ('Aout2(idSig12)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idSig12,:))
            CASE ('Aout2(idSig22)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idSig22,:))
            CASE ('Aout2(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idT0mk,:))
            CASE ('Aout2(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idS0mk,:))
#endif
#ifdef HISTORY2
            CASE ('Hout2(idUice)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUice,:))
            CASE ('Hout2(idVice)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVice,:))
            CASE ('Hout2(idUiceE)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUiceE,:))
            CASE ('Hout2(idViceN)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idViceN,:))
            CASE ('Hout2(idAice)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idAice,:))
            CASE ('Hout2(idHice)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHice,:))
            CASE ('Hout2(idHsno)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHsno,:))
            CASE ('Hout2(idTice)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idTice,:))
            CASE ('Hout2(idTimid)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idTimid,:))
            CASE ('Hout2(idTauiw)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idTauiw,:))
            CASE ('Hout2(idChuiw)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idChuiw,:))
#ifdef MELT_PONDS
            CASE ('Hout2(idApond)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idApond,:))
            CASE ('Hout2(idHpond)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHpond,:))
#endif
            CASE ('Hout2(idAgeice)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idAgeice,:))
            CASE ('Hout2(idIomflx)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idIomflx,:))
            CASE ('Hout2(idSig11)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idSig11,:))
            CASE ('Hout2(idSig12)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idSig12,:))
            CASE ('Hout2(idSig22)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idSig22,:))
            CASE ('Hout2(idT0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idT0mk,:))
            CASE ('Hout2(idS0mk)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idS0mk,:))
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
#ifdef MELT_PONDS
          IF (Hout(idApond,ng)) WRITE (out,70) Hout(idApond,ng),        &
     &       'Hout(idApond)',                                           &
     &       'Write out surface water (on ice) fraction.'
          IF (Hout(idHpond,ng)) WRITE (out,70) Hout(idHpond,ng),        &
     &       'Hout(idHpond)',                                           &
     &       'Write out surface water (on ice) thickness.'
#endif
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
#ifdef MELT_PONDS
          IF (Aout(idApond,ng)) WRITE (out,70) Aout(idApond,ng),        &
     &       'Aout(idApond)',                                           &
     &       'Write out time-averaged surface water (on ice) fraction.'
          IF (Aout(idHpond,ng)) WRITE (out,70) Aout(idHpond,ng),        &
     &       'Aout(idHpond)',                                           &
     &       'Write out time-averaged surface water (on ice) thickness.'
#endif
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
#ifdef MELT_PONDS
          IF (Aout2(idApond,ng)) WRITE (out,70) Aout2(idApond,ng),      &
     &       'Aout2(idApond)',                                          &
     &       'Write out surface water (on ice) fraction.'
          IF (Aout2(idHpond,ng)) WRITE (out,70) Aout2(idHpond,ng),      &
     &       'Aout2(idHpond)',                                          &
     &       'Write out surface water (on ice) thickness.'
#endif
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
#ifdef HISTORY2
          IF (Hout2(idUice,ng)) WRITE (out,70) Hout2(idUice,ng),        &
     &       'Hout2(idUice)',                                           &
     &       'Write out U-component ice velocity.'
          IF (Hout2(idVice,ng)) WRITE (out,70) Hout2(idVice,ng),        &
     &       'Hout2(idVice)',                                           &
     &       'Write out V-component ice velocity.'
          IF (Hout2(idUiceE,ng)) WRITE (out,70) Hout2(idUiceE,ng),      &
     &       'Hout2(idUiceE)',                                          &
     &       'Write out East component ice velocity.'
          IF (Hout2(idViceN,ng)) WRITE (out,70) Hout2(idViceN,ng),      &
     &       'Hout2(idViceN)',                                          &
     &       'Write out North component ice velocity.'
          IF (Hout2(idAice,ng)) WRITE (out,70) Hout2(idAice,ng),        &
     &       'Hout2(idAice)',                                           &
     &       'Write out ice concentration.'
          IF (Hout2(idHice,ng)) WRITE (out,70) Hout2(idHice,ng),        &
     &       'Hout2(idHice)',                                           &
     &       'Write out average ice thickness.'
          IF (Hout2(idHsno,ng)) WRITE (out,70) Hout2(idHsno,ng),        &
     &       'Hout2(idHsno)',                                           &
     &       'Write out snow thickness.'
          IF (Hout2(idTice,ng)) WRITE (out,70) Hout2(idTice,ng),        &
     &       'Hout2(idTice)',                                           &
     &       'Write out ice/snow surface temperature.'
          IF (Hout2(idTimid,ng)) WRITE (out,70) Hout2(idTimid,ng),      &
     &       'Hout2(idTimid)',                                          &
     &       'Write out interior ice temperature.'
#ifdef MELT_PONDS
          IF (Hout2(idApond,ng)) WRITE (out,70) Hout2(idApond,ng),      &
     &       'Hout2(idApond)',                                          &
     &       'Write out surface water (on ice) fraction.'
          IF (Hout2(idHpond,ng)) WRITE (out,70) Hout2(idHpond,ng),      &
     &       'Hout2(idHpond)',                                          &
     &       'Write out surface water (on ice) thickness.'
#endif
          IF (Hout2(idAgeice,ng)) WRITE (out,70) Hout2(idAgeice,ng),    &
     &       'Hout2(idAgeice)',                                         &
     &       'Write out ice age.'
          IF (Hout2(idIomflx,ng)) WRITE (out,70) Hout2(idIomflx,ng),    &
     &       'Hout2(idIomflx)',                                         &
     &       'Write out ice-ocean mass flux'
          IF (Hout2(idSig11,ng)) WRITE (out,70) Hout2(idSig11,ng),      &
     &       'Hout2(idSig11)',                                          &
     &       'Write out internal ice stress component 11.'
          IF (Hout2(idSig12,ng)) WRITE (out,70) Hout2(idSig12,ng),      &
     &       'Hout2(idSig12)',                                          &
     &       'Write out internal ice stress component 12.'
          IF (Hout2(idSig22,ng)) WRITE (out,70) Hout2(idSig22,ng),      &
     &       'Hout2(idSig22)',                                          &
     &       'Write out internal ice stress component 22.'
          IF (Hout2(idTauiw,ng)) WRITE (out,70) Hout2(idTauiw,ng),      &
     &       'Hout2(idTauiw)',                                          &
     &       'Write out ice-water friction velocity.'
          IF (Hout2(idChuiw,ng)) WRITE (out,70) Hout2(idChuiw,ng),      &
     &       'Hout2(idChuiw)',                                          &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Hout2(idT0mk,ng)) WRITE (out,70) Hout2(idT0mk,ng),        &
     &       'Hout2(idT0mk)',                                           &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Hout2(idS0mk,ng)) WRITE (out,70) Hout2(idS0mk,ng),        &
     &       'Hout2(idS0mk)',                                           &
     &       'Write out salinity of molecular sublayer under ice.'
#endif
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
