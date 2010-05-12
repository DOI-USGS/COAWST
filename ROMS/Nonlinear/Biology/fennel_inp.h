      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: fennel_inp.h 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in Fennel et al. (2006) ecosystem model input    !
!  parameters. They are specified in input script "fennel.in".         !
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

      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in Fennel et al. (2006) biological model parameters.
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
          ELSE IF (TRIM(KeyWord).eq.'AttChl') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttChl)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'Vp0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vp0)
          ELSE IF (TRIM(KeyWord).eq.'I_thNH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, I_thNH4)
          ELSE IF (TRIM(KeyWord).eq.'D_p5NH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, D_p5NH4)
          ELSE IF (TRIM(KeyWord).eq.'NitriR') THEN
            Npts=load_r(Nval, Rval, Ngrids, NitriR)
          ELSE IF (TRIM(KeyWord).eq.'K_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NO3)
          ELSE IF (TRIM(KeyWord).eq.'K_NH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NH4)
          ELSE IF (TRIM(KeyWord).eq.'K_Phy') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_Phy)
          ELSE IF (TRIM(KeyWord).eq.'Chl2C_m') THEN
            Npts=load_r(Nval, Rval, Ngrids, Chl2C_m)
          ELSE IF (TRIM(KeyWord).eq.'ChlMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, ChlMin)
          ELSE IF (TRIM(KeyWord).eq.'PhyCN') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyCN)
          ELSE IF (TRIM(KeyWord).eq.'PhyIP') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIP)
          ELSE IF (TRIM(KeyWord).eq.'PhyIS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIS)
          ELSE IF (TRIM(KeyWord).eq.'PhyMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMin)
          ELSE IF (TRIM(KeyWord).eq.'PhyMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMR)
          ELSE IF (TRIM(KeyWord).eq.'ZooAE_N') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooAE_N)
          ELSE IF (TRIM(KeyWord).eq.'ZooBM') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooBM)
          ELSE IF (TRIM(KeyWord).eq.'ZooCN') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooCN)
          ELSE IF (TRIM(KeyWord).eq.'ZooER') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooER)
          ELSE IF (TRIM(KeyWord).eq.'ZooGR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGR)
          ELSE IF (TRIM(KeyWord).eq.'ZooMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMin)
          ELSE IF (TRIM(KeyWord).eq.'ZooMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMR)
          ELSE IF (TRIM(KeyWord).eq.'LDeRRN') THEN
            Npts=load_r(Nval, Rval, Ngrids, LDeRRN)
          ELSE IF (TRIM(KeyWord).eq.'LDeRRC') THEN
            Npts=load_r(Nval, Rval, Ngrids, LDeRRC)
          ELSE IF (TRIM(KeyWord).eq.'CoagR') THEN
            Npts=load_r(Nval, Rval, Ngrids, CoagR)
          ELSE IF (TRIM(KeyWord).eq.'SDeRRN') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeRRN)
          ELSE IF (TRIM(KeyWord).eq.'SDeRRC') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeRRC)
          ELSE IF (TRIM(KeyWord).eq.'wPhy') THEN
            Npts=load_r(Nval, Rval, Ngrids, wPhy)
          ELSE IF (TRIM(KeyWord).eq.'wLDet') THEN
            Npts=load_r(Nval, Rval, Ngrids, wLDet)
          ELSE IF (TRIM(KeyWord).eq.'wSDet') THEN
            Npts=load_r(Nval, Rval, Ngrids, wSDet)
          ELSE IF (TRIM(KeyWord).eq.'pCO2air') THEN
            Npts=load_r(Nval, Rval, Ngrids, pCO2air)
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
#ifdef DIAGNOSTICS_BIO
# ifdef CARBON
          ELSE IF (TRIM(KeyWord).eq.'Hout(iCOfx)') THEN
            IF (iDbio2(iCOfx).eq.0) THEN
              IF (Master) WRITE (out,40) 'iDbio2(iCOfx)'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(iCOfx)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
# endif
# ifdef DENITRIFICATION
          ELSE IF (TRIM(KeyWord).eq.'Hout(iDNIT)') THEN
            IF (iDbio2(iDNIT).eq.0) THEN
              IF (Master) WRITE (out,40) 'iDbio2(iDNIT)'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(iDNIT)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
# endif
# ifdef CARBON
          ELSE IF (TRIM(KeyWord).eq.'Hout(ipCO2)') THEN
            IF (iDbio2(ipCO2).eq.0) THEN
              IF (Master) WRITE (out,40) 'iDbio2(ipCO2)'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(ipCO2)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
# endif
# ifdef OXYGEN
          ELSE IF (TRIM(KeyWord).eq.'Hout(iO2fx)') THEN
            IF (iDbio2(iO2fx).eq.0) THEN
              IF (Master) WRITE (out,40) 'iDbio2(iO2fx)'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(iO2fx)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
# endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(iPPro)') THEN
            IF (iDbio3(iPPro).eq.0) THEN
              IF (Master) WRITE (out,40) 'iDbio3(iPPro)'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio3(iPPro)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iNO3u)') THEN
            IF (iDbio3(iNO3u).eq.0) THEN
              IF (Master) WRITE (out,40) 'iDbio3(iNO3u)'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio3(iNO3u)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
#endif
          END IF
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
#ifdef TS_PSOURCE
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,110) LtracerSrc(i,ng), 'LtracerSrc',           &
     &              i, 'Processing point sources/Sink on tracer ', i,   &
     &              TRIM(Vname(1,idTvar(i)))
            END DO
#endif
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
#ifdef DIAGNOSTICS_BIO
# if !(defined CARBON || defined OXYGEN || defined DENITRIFICATION)
            DO itrc=1,NDbio2d
              i=iDbio2(itrc)
              IF (Hout(i,ng)) WRITE (out,130)                           &
     &            Hout(i,ng), 'Hout(iDbio2)', 'Diagnostics for',        &
     &            TRIM(Vname(1,i))
            END DO
# endif
            DO itrc=1,NDbio3d
              i=iDbio3(itrc)
              IF (Hout(i,ng)) WRITE (out,130)                           &
     &            Hout(i,ng), 'Hout(iDbio3)', 'Diagnostics for',        &
     &            TRIM(Vname(1,i))
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
  40  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
  50  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Fennel Model Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  70  FORMAT (1x,i10,2x,a,t30,a)
  80  FORMAT (1p,e11.4,2x,a,t30,a)
  90  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 100  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 110  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 120  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 130  FORMAT (10x,l1,2x,a,t30,a,1x,a)

      RETURN
      END SUBROUTINE read_BioPar
