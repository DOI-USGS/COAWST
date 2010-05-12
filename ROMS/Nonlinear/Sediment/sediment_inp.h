      SUBROUTINE read_SedPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in cohesive and non-cohesive sediment model      !
!  parameters.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
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
      integer :: Npts, Nval, i, ng, itrc, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(Ngrids) :: Lbed
      logical, dimension(MBOTP,Ngrids) :: Lbottom
      logical, dimension(NCS,Ngrids) :: Lmud
      logical, dimension(NNS,Ngrids) :: Lsand

      real(r8), dimension(Ngrids) :: Rbed
      real(r8), dimension(NCS,Ngrids) :: Rmud
      real(r8), dimension(NNS,Ngrids) :: Rsand

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in cohesive and non-cohesive model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lsediment') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lsediment)
          ELSE IF (TRIM(KeyWord).eq.'NEWLAYER_THICK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              newlayer_thick(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MINLAYER_THICK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              minlayer_thick(ng)=Rbed(ng)
            END DO
#ifdef MIXED_BED
          ELSE IF (TRIM(KeyWord).eq.'TRANSC') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              transC(ng)=Rbed(ng)
          END DO
          ELSE IF (TRIM(KeyWord).eq.'TRANSN') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              transN(ng)=Rbed(ng)
          END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'BEDLOAD_COEFF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              bedload_coeff(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_SD50') THEN
            IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Sd50(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_CSED') THEN
            IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud )
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Csed(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_SRHO') THEN
            IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Srho(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_WSED') THEN
            IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Wsed(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_ERATE') THEN
            IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Erate(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAU_CE') THEN
            IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                tau_ce(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAU_CD') THEN
            IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                tau_cd(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_POROS') THEN
            IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                poros(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNU2') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                nl_tnu2(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNU4') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                nl_tnu4(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_MUD_TNU2') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                ad_tnu2(i,ng)=Rmud(itrc,ng)
                tl_tnu2(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_MUD_TNU4') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                ad_tnu4(i,ng)=Rmud(itrc,ng)
                nl_tnu4(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                Akt_bak(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_AKT_fac') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                ad_Akt_fac(i,ng)=Rmud(itrc,ng)
                tl_Akt_fac(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNUDG') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                Tnudg(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_MORPH_FAC') THEN
            IF (.not.allocated(morph_fac)) THEN
              allocate (morph_fac(NST,Ngrids))
            END IF
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                morph_fac(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
#if defined COHESIVE_BED || defined MIXED_BED
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_MIN') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_min(ng)=Rbed(ng)
          END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_MAX') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_max(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_SLOPE') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_slp(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_OFF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_off(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_TIME') THEN
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_tim(ng)=Rbed(ng)
          END DO
#endif
#ifdef TS_PSOURCE
          ELSE IF (TRIM(KeyWord).eq.'MUD_Ltracer') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                LtracerSrc(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idmud)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idTvar(idsed(itrc))
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMfrac)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idfrac(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMmass)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idBmas(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
#ifdef BEDLOAD
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMUbld)') THEN
            DO ng=1,Ngrids
              DO itrc=1,NCS
                IF (idUbld(itrc).eq.0) THEN
                  IF (Master) WRITE (out,30) 'idUbld'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idUbld(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMVbld)') THEN
            DO ng=1,Ngrids
              DO itrc=1,NCS
                IF (idVbld(itrc).eq.0) THEN
                  IF (Master) WRITE (out,30) 'idVbld'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idVbld(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'SAND_SD50') THEN
            IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Sd50(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_CSED') THEN
            IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand )
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Csed(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_SRHO') THEN
            IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Srho(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_WSED') THEN
            IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Wsed(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_ERATE') THEN
            IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Erate(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TAU_CE') THEN
            IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                tau_ce(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TAU_CD') THEN
            IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                tau_cd(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_POROS') THEN
            IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                poros(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNU2') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                nl_tnu2(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNU4') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                nl_tnu4(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_SAND_TNU2') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                ad_tnu2(i,ng)=Rsand(itrc,ng)
                tl_tnu2(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_SAND_TNU4') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                ad_tnu4(i,ng)=Rsand(itrc,ng)
                tl_tnu4(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                Akt_bak(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_AKT_fac') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                ad_Akt_fac(i,ng)=Rsand(itrc,ng)
                tl_Akt_fac(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNUDG') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                Tnudg(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_MORPH_FAC') THEN
            IF (.not.allocated(morph_fac)) THEN
              allocate (morph_fac(NST,Ngrids))
            END IF
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                morph_fac(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
#ifdef TS_PSOURCE
          ELSE IF (TRIM(KeyWord).eq.'SAND_Ltracer') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                LtracerSrc(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idsand)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idTvar(idsed(NCS+itrc))
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSfrac)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idfrac(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSmass)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idBmas(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
#ifdef BEDLOAD
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSUbld)') THEN
            DO ng=1,Ngrids
              DO itrc=NCS+1,NST
                IF (idUbld(itrc).eq.0) THEN
                  IF (Master) WRITE (out,30) 'idUbld'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idUbld(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSVbld)') THEN
            DO ng=1,Ngrids
              DO itrc=NCS+1,NST
                IF (idVbld(itrc).eq.0) THEN
                  IF (Master) WRITE (out,30) 'idVbld'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idVbld(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(ithck)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(ithck)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iaged)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(iaged)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iporo)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(iporo)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
          ELSE IF (TRIM(KeyWord).eq.'Hout(ibtcr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(ibtcr)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idiff)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(idiff)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
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
          IF (Lsediment(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60)
            DO itrc=1,NST
              WRITE (out,70) itrc, Sd50(itrc,ng), Csed(itrc,ng),        &
     &                       Srho(itrc,ng), Wsed(itrc,ng),              &
     &                       Erate(itrc,ng), poros(itrc,ng)
            END DO
            WRITE (out,80)
            DO itrc=1,NST
              i=idsed(itrc)
              WRITE (out,70) itrc, tau_ce(itrc,ng), tau_cd(itrc,ng),    &
     &                       nl_tnu2(i,ng), nl_tnu4(i,ng),              &
     &                       Akt_bak(i,ng), Tnudg(i,ng)
            END DO
            WRITE (out,90)
            DO itrc=1,NST
              WRITE (out,70) itrc,  morph_fac(itrc,ng)
            END DO
            WRITE (out,100) newlayer_thick(ng)
            WRITE (out,110) minlayer_thick(ng)
            WRITE (out,120) bedload_coeff(ng)
#ifdef MIXED_BED
            WRITE (out,130) transC(ng)
            WRITE (out,140) transN(ng)
#endif
#ifdef TS_PSOURCE
            DO itrc=1,NST
              i=idsed(itrc)
              WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc',           &
     &              i, 'Processing point sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
            DO itrc=1,NST
              i=idTvar(idsed(itrc))
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &            'Hout(idTvar)',                                       &
     &            'Write out sediment', itrc, TRIM(Vname(1,i))
            END DO
            DO itrc=1,NST
              i=idfrac(itrc)
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &            'Hout(idfrac)',                                       &
     &            'Write out bed fraction, sediment ', itrc,            &
     &            TRIM(Vname(1,i))
            END DO
            DO itrc=1,NST
              i=idBmas(itrc)
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &            'Hout(idfrac)',                                       &
     &            'Write out mass, sediment ', itrc,                    &
     &            TRIM(Vname(1,i))
            END DO
#ifdef BEDLOAD
            DO itrc=1,NST
              i=idUbld(itrc)
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &            'Hout(idUbld)',                                       &
     &            'Write out bed load at U-points, sediment ', itrc,    &
     &            TRIM(Vname(1,i))
              i=idVbld(itrc)
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &            'Hout(idVbld)',                                       &
     &            'Write out bed load at V-points, sediment ', itrc,    &
     &            TRIM(Vname(1,i))
            END DO
#endif
            DO itrc=1,MBEDP
              i=idSbed(itrc)
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &            'Hout(idSbed)',                                       &
     &            'Write out BED property ', itrc, TRIM(Vname(1,i))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant input parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NST
          Sd50(i,ng)=Sd50(i,ng)*0.001_r8
          Wsed(i,ng)=Wsed(i,ng)*0.001_r8
          tau_ce(i,ng)=tau_ce(i,ng)/rho0
          tau_cd(i,ng)=tau_cd(i,ng)/rho0
          nl_tnu4(idsed(i),ng)=SQRT(ABS(nl_tnu4(idsed(i),ng)))
#ifdef ADJOINT
          ad_tnu4(idsed(i),ng)=SQRT(ABS(ad_tnu4(idsed(i),ng)))
#endif
#if defined TANGENT || defined TL_IOMS
          tl_tnu4(idsed(i),ng)=SQRT(ABS(tl_tnu4(idsed(i),ng)))
#endif
          IF (Tnudg(idsed(i),ng).gt.0.0_r8) THEN
            Tnudg(idsed(i),ng)=1.0_r8/(Tnudg(idsed(i),ng)*86400.0_r8)
          ELSE
            Tnudg(idsed(i),ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' READ_SedPar - variable info not yet loaded, ', a)
  40  FORMAT (/,' READ_SedPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' Sediment Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  60  FORMAT (/,1x,'Size',5x,'Sd50',8x,'Csed',8x,'Srho',8x,'Wsed',      &
     &        8x,'Erate',7x,'poros',/,1x,'Class',4x,'(mm)',7x,          &
     &        '(kg/m3)',5x,'(kg/m3)',5x,'(mm/s)',5x,'(kg/m2/s)',4x,     &
     &        '(nondim)',/)
  70  FORMAT (2x,i2,2x,6(1x,1p,e11.4))
  80  FORMAT (/,9x,'tau_ce',6x,'tau_cd',6x,'nl_tnu2',5x,'nl_tnu4',5x,   &
     &        'Akt_bak',6x,'Tnudg',/,9x,'(N/m2)',6x,'(N/m2)',6x,        &
     &        '(m2/s)',6x,'(m4/s)',7x,'(m2/s)',6x,'(day)',/)
  90  FORMAT (/,9x,'morph_fac',/,9x,'(nondim)',/)
 100  FORMAT (/,' New bed layer formed when deposition exceeds ',e11.5, &
     &        ' (m).')
 110  FORMAT (' Two first layers are combined when 2nd layer smaller ', &
     &         'than ',e11.5,' (m).')
 120  FORMAT (' Rate coefficient for bed load transport = ',e11.5,/)
 130  FORMAT (' Transition for mixed sediment =',e11.5,/)
 140  FORMAT (' Transition for cohesive sediment =',e11.5,/)
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 160  FORMAT (10x,l1,2x,a,t29,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_SedPar

