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
      USE mod_sedflocs
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

      logical, dimension(Ngrids) :: Lbed
      logical, dimension(Ngrids) :: Lbottom
      logical, dimension(NCS,Ngrids) :: Lmud
      logical, dimension(NNS,Ngrids) :: Lsand

      real(r8), dimension(Ngrids) :: Rbed
      real(r8), dimension(NCS,Ngrids) :: Rmud
      real(r8), dimension(NNS,Ngrids) :: Rsand

      real(r8), dimension(100) :: Rval

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
      iTrcEnd=NST                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in cohesive and non-cohesive model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lsediment')
              Npts=load_l(Nval, Cval, Ngrids, Lsediment)
            CASE ('NEWLAYER_THICK')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                newlayer_thick(ng)=Rbed(ng)
              END DO
            CASE ('MINLAYER_THICK')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                minlayer_thick(ng)=Rbed(ng)
              END DO
#ifdef MIXED_BED
            CASE ('TRANSC')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                transC(ng)=Rbed(ng)
              END DO
            CASE ('TRANSN')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                transN(ng)=Rbed(ng)
              END DO
#endif
            CASE ('BEDLOAD_COEFF')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                bedload_coeff(ng)=Rbed(ng)
              END DO
#ifdef SED_BIODIFF
	        CASE ('DBMAX')
	          Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
	            Dbmx(ng)=Rbed(ng)
              END DO
            CASE ('DBMIN')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
	            Dbmm(ng)=Rbed(ng)
              END DO
            CASE ('DBZS')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                Dbzs(ng)=Rbed(ng)
              END DO
            CASE ('DBZM')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                Dbzm(ng)=Rbed(ng)
              END DO
            CASE ('DBZP')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                Dbzp(ng)=Rbed(ng)
              END DO
#endif
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NST) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idsed(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idsed(iTrcStr), idsed(iTrcEnd),             &
     &                      Vname(1,idTvar(idsed(itracer))), LBC)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_LBC(isTvar)')
              IF (itracer.lt.NST) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idsed(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idsed(iTrcStr), idsed(iTrcEnd),             &
     &                      Vname(1,idTvar(idsed(itracer))), ad_LBC)
#endif
            CASE ('MUD_SD50')
              IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Sd50(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_CSED')
              IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud )
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Csed(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_SRHO')
              IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Srho(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_WSED')
              IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Wsed(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_ERATE')
              IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Erate(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TAU_CE')
              IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  tau_ce(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TAU_CD')
              IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  tau_cd(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_POROS')
              IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  poros(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNU2')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  nl_tnu2(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNU4')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  nl_tnu4(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('ad_MUD_TNU2')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_tnu2(i,ng)=Rmud(itrc,ng)
                  tl_tnu2(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('ad_MUD_TNU4')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_tnu4(i,ng)=Rmud(itrc,ng)
                  nl_tnu4(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Sponge')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerSponge(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_AKT_BAK')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  Akt_bak(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_AKT_fac')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_Akt_fac(i,ng)=Rmud(itrc,ng)
                  tl_Akt_fac(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNUDG')
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  Tnudg(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_MORPH_FAC')
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
            CASE ('MUD_TAUCR_MIN')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                tcr_min(ng)=Rbed(ng)
              END DO
            CASE ('MUD_TAUCR_MAX')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                tcr_max(ng)=Rbed(ng)
              END DO
            CASE ('MUD_TAUCR_SLOPE')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                tcr_slp(ng)=Rbed(ng)
              END DO
            CASE ('MUD_TAUCR_OFF')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                tcr_off(ng)=Rbed(ng)
              END DO
            CASE ('MUD_TAUCR_TIME')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                tcr_tim(ng)=Rbed(ng)
              END DO
#endif
            CASE ('MUD_Ltsrc', 'MUD_Ltracer')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerSrc(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Ltclm')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerCLM(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Tnudge')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LnudgeTCLM(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idmud)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMfrac)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idfrac(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMmass)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idBmas(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
#ifdef BEDLOAD
            CASE ('Hout(iMUbld)')
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
            CASE ('Hout(iMVbld)')
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
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            CASE ('Aout(idmud)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMTTav)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMUTav)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idUTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMVTav)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idVTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(MHUTav)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=iHUTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(MHVTav)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=iHVTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
# ifdef BEDLOAD
            CASE ('Aout(iMUbld)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idUbld(itrc)
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMVbld)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idVbld(itrc)
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
# endif
#endif
#ifdef DIAGNOSTICS_TS
            CASE ('Dout(MTrate)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTrate),ng)=Lmud(i,ng)
                END DO
              END DO
            CASE ('Dout(MThadv)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iThadv),ng)=Lmud(i,ng)
                END DO
              END DO
            CASE ('Dout(MTxadv)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Lmud(i,ng)
                END DO
              END DO
            CASE ('Dout(MTyadv)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Lmud(i,ng)
                END DO
              END DO
            CASE ('Dout(MTvadv)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Lmud(i,ng)
                END DO
              END DO
# if defined TS_DIF2 || defined TS_DIF4
            CASE ('Dout(MThdif)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iThdif),ng)=Lmud(i,ng)
                END DO
              END DO
            CASE ('Dout(MTxdif)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Lmud(i,ng)
                END DO
              END DO
            CASE ('Dout(MTydif)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTydif),ng)=Lmud(i,ng)
                END DO
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            CASE ('Dout(MTsdif)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Lmud(i,ng)
                END DO
              END DO
#  endif
# endif
            CASE ('Dout(MTvdif)')
              Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
              DO ng=1,Ngrids
                DO i=1,NCS
                  itrc=idsed(i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Lmud(i,ng)
                END DO
              END DO
#endif
            CASE ('SAND_SD50')
              IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Sd50(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_CSED')
              IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand )
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Csed(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_SRHO')
              IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Srho(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_WSED')
              IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Wsed(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_ERATE')
              IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Erate(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TAU_CE')
              IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  tau_ce(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TAU_CD')
              IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  tau_cd(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_POROS')
              IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  poros(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNU2')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  nl_tnu2(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNU4')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  nl_tnu4(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('ad_SAND_TNU2')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_tnu2(i,ng)=Rsand(itrc,ng)
                  tl_tnu2(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('ad_SAND_TNU4')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_tnu4(i,ng)=Rsand(itrc,ng)
                  tl_tnu4(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Sponge')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerSponge(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_AKT_BAK')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  Akt_bak(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_AKT_fac')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_Akt_fac(i,ng)=Rsand(itrc,ng)
                  tl_Akt_fac(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNUDG')
              Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  Tnudg(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_MORPH_FAC')
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
            CASE ('SAND_Ltsrc', 'SAND_Ltracer')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerSrc(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Ltclm')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerCLM(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Tnudge')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LnudgeTCLM(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsand)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSfrac)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idfrac(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSmass)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idBmas(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
#ifdef BEDLOAD
            CASE ('Hout(iSUbld)')
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
            CASE ('Hout(iSVbld)')
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
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            CASE ('Aout(idsand)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSTTav)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSUTav)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idUTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSVTav)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idVTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(SHUTav)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=iHUTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(SHVTav)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=iHVTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
# ifdef BEDLOAD
            CASE ('Aout(iSUbld)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idUbld(NCS+itrc)
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSVbld)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idVbld(NCS+itrc)
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
# endif
#endif
#ifdef DIAGNOSTICS_TS
            CASE ('Dout(STrate)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTrate),ng)=Lsand(i,ng)
                END DO
              END DO
            CASE ('Dout(SThadv)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iThadv),ng)=Lsand(i,ng)
                END DO
              END DO
            CASE ('Dout(STxadv)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Lsand(i,ng)
                END DO
              END DO
            CASE ('Dout(STyadv)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Lsand(i,ng)
                END DO
              END DO
            CASE ('Dout(STvadv)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Lsand(i,ng)
                END DO
              END DO
# if defined TS_DIF2 || defined TS_DIF4
            CASE ('Dout(SThdif)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iThdif),ng)=Lsand(i,ng)
                END DO
              END DO
            CASE ('Dout(STxdif)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Lsand(i,ng)
                END DO
              END DO
            CASE ('Dout(STydif)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTydif),ng)=Lsand(i,ng)
                END DO
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            CASE ('Dout(STsdif)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Lsand(i,ng)
                END DO
              END DO
#  endif
# endif
            CASE ('Dout(STvdif)')
              Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
              DO ng=1,Ngrids
                DO i=1,NNS
                  itrc=idsed(NCS+i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Lsand(i,ng)
                END DO
              END DO
#endif
            CASE ('Hout(ithck)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(ithck)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(iaged)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iaged)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(iporo)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iporo)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
            CASE ('Hout(ibtcr)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(ibtcr)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
#endif
            CASE ('Hout(isd50)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(isd50)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idens)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idens)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
             END DO
            CASE ('Hout(iwsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iwsed)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(itauc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(itauc)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(irlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irlen)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(irhgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irhgt)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(ibwav)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ibwav)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izdef)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izdef)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izapp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izapp)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izNik)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izNik)
             DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbio)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbio)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbfm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbfm)
              DO ng=1,Ngrids
               Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbld)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbld)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izwbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izwbl)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(iactv)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iactv)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(ishgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ishgt)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(imaxD)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(imaxD)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idnet)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idnet)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
#if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
            CASE ('Hout(idoff)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idoff)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idslp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idslp)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idtim)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idtim)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idbmx)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idbmx)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idbmm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idbmm)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idbzs)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idbzs)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idbzm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idbzm)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idbzp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idbzp)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
#endif 
#if defined MIXED_BED
            CASE ('Hout(idprp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idprp)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
#endif
#if defined SED_FLOCS
            CASE ('l_ADS')
              Npts=load_l(Nval, Cval, Ngrids, l_ADS)
            CASE ('l_ASH')
              Npts=load_l(Nval, Cval, Ngrids, l_ASH)
            CASE ('l_COLLFRAG')
              Npts=load_l(Nval, Cval, Ngrids, l_COLLFRAG)
            CASE ('f_dp0')
              Npts=load_r(Nval, Rval, Ngrids, f_dp0)
            CASE ('f_nf')
              Npts=load_r(Nval, Rval, Ngrids, f_nf)
            CASE ('f_dmax')
              Npts=load_r(Nval, Rval, Ngrids, f_dmax)
            CASE ('f_nb_frag')
              Npts=load_r(Nval, Rval, Ngrids, f_nb_frag)
            CASE ('f_alpha')
              Npts=load_r(Nval, Rval, Ngrids, f_alpha)
            CASE ('f_beta')
              Npts=load_r(Nval, Rval, Ngrids, f_beta)
            CASE ('f_ater')
              Npts=load_r(Nval, Rval, Ngrids, f_ater)
            CASE ('f_ero_frac')
              Npts=load_r(Nval, Rval, Ngrids, f_ero_frac)
            CASE ('f_ero_nbfrag')
              Npts=load_r(Nval, Rval, Ngrids, f_ero_nbfrag)
            CASE ('f_ero_iv')
              Npts=load_r(Nval, Rval, Ngrids, f_ero_iv)
            CASE ('f_collfragparam')
              Npts=load_r(Nval, Rval, Ngrids, f_collfragparam)
            CASE ('f_clim')
              Npts=load_r(Nval, Rval, Ngrids, f_clim)
            CASE ('l_testcase')
              Npts=load_l(Nval, Cval, Ngrids, l_testcase)
#endif
#if defined SED_FLOCS && defined SED_DEFLOC
            CASE ('MUD_FRAC_EQ')
              IF (.not.allocated(mud_frac_eq)) THEN 
                allocate (mud_frac_eq(NCS,Ngrids))
              ENDIF
              Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  mud_frac_eq(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
	    CASE ('MUD_T_DFLOC')
!              Npts=load_r(Nval, Rval, Ngrids, t_dfloc)
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  t_dfloc(ng)=Rbed(ng)
                END DO
#endif
            CASE ('Hout(idiff)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(idiff)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
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
            DO itrc=1,NST
              i=idsed(itrc)
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
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc', i,      &
     &              'Turning ON  point sources/Sink on tracer ', i,     &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc', i,      &
     &              'Turning OFF point sources/Sink on tracer ', i,     &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
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
            DO itrc=1,NST
              i=idsed(itrc)
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
            DO itrc=1,MBOTP
              i=idBott(itrc)
              IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),               &
     &           'Hout(idBott)',                                        &
     &           'Write out BOTTOM property', itrc, TRIM(Vname(1,i))
            END DO
!
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            WRITE (out,'(1x)')
            DO itrc=1,NST
              i=idTvar(idsed(itrc))
              IF (Aout(i,ng)) WRITE (out,160) Aout(i,ng),               &
     &            'Aout(idTvar)',                                       &
     &            'Write out averaged sediment', itrc, TRIM(Vname(1,i))
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (Aout(idTTav(i),ng)) WRITE (out,160)                   &
     &            Aout(idTTav(i),ng), 'Aout(idTTav)',                   &
     &            'Write out averaged <t*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (Aout(idUTav(i),ng)) WRITE (out,160)                   &
     &            Aout(idUTav(i),ng), 'Aout(idUTav)',                   &
     &            'Write out averaged <u*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (Aout(idVTav(i),ng)) WRITE (out,160)                   &
     &            Aout(idVTav(i),ng), 'Aout(idVTav)',                   &
     &            'Write out averaged <v*t> for tracer ', i,            &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (Aout(iHUTav(i),ng)) WRITE (out,160)                   &
     &            Aout(iHUTav(i),ng), 'Aout(iHUTav)',                   &
     &            'Write out averaged <Huon*t> for tracer ', i,         &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (Aout(iHVTav(i),ng)) WRITE (out,160)                   &
     &            Aout(iHVTav(i),ng), 'Aout(iHVTav)',                   &
     &            'Write out averaged <Hvom*t> for tracer ', i,         &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
# ifdef BEDLOAD
            DO itrc=1,NST
              i=idUbld(itrc)
              IF (Aout(i,ng)) WRITE (out,160) Aout(i,ng),               &
     &            'Aout(idUbld)',                                       &
     &            'Write out U-bedload, sediment ', itrc,               &
     &            TRIM(Vname(1,i))
              i=idVbld(itrc)
              IF (Aout(i,ng)) WRITE (out,160) Aout(i,ng),               &
     &            'Aout(idVbld)',                                       &
     &            'Write out V-bedload, sediment ', itrc,               &
     &            TRIM(Vname(1,i))
            END DO
# endif
#endif
#ifdef DIAGNOSTICS_TS
            WRITE (out,'(1x)')
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iTrate),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTrate)',                 &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iThadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iThadv)',                 &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iTxadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTxadv)',                 &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iTyadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTyadv)',                 &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iTvadv),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTvadv)',                 &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
# if defined TS_DIF2 || defined TS_DIF4
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iThdif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iThdif)',                 &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(i,iTxdif),ng))                            &
     &          WRITE (out,160) .TRUE., 'Dout(iTxdif)',                 &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iTydif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTydif)',                 &
     &            'Write out horizontal Y-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            DO i=1,NST
              itrc=idsed(i)
              IF (Dout(idDtrc(itrc,iTsdif),ng))                         &
     &          WRITE (out,160) .TRUE., 'Dout(iTsdif)',                 &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#  endif
# endif
            DO i=1,NST
              itrc=idsed(i)
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

