      SUBROUTINE read_VegPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in vegetation model parameters.                  !
!  Equivalent of read_phypar.F so it gets read in that                 !
!  This routine also outputs vegetation model parameters.              !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE mod_vegetation
!
      USE inp_decode_mod
#if defined MARSH_SED_EROSION 
      USE mod_sediment 
#endif 
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
      integer :: iveg, ised, it
!      integer :: decode_line, load_i, load_l, load_lbc, load_r
!
      real(r8), dimension(200) :: Rval
#if defined MARSH_SED_EROSION  || defined MARSH_VERT_GROWTH
      real(r8), dimension(Ngrids) :: Rmarsh
#endif 
#ifdef VEG_DRAG 
      real(r8), allocatable :: Rveg(:,:)
#endif 
      logical, dimension(NNS,Ngrids) :: Lsand1
!
      character (len=40 ) :: KeyWord
      character (len=256) :: line 
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Read input parameters from vegetation.in
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
#ifdef VEG_DRAG 
            CASE ('NVEG') 
              Npts=load_i(Nval, Rval, Ngrids, NVEG)
                IF (NVEG.lt.0) THEN
                  IF (Master) WRITE (out,30) 'NVEG', ng,                &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
            IF (.not.allocated(Rveg)) allocate(Rveg(NVEG,Ngrids)) 
            CASE ('CD_VEG')
              IF (.not.allocated(CD_VEG)) allocate(CD_VEG(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  CD_VEG(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
            CASE ('E_VEG') 
              IF (.not.allocated(E_VEG)) allocate(E_VEG(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  E_VEG(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
            CASE ('VEG_MASSDENS') 
              IF (.not.allocated(VEG_MASSDENS))                         &
     &                 allocate(VEG_MASSDENS(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  VEG_MASSDENS(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
            CASE ('VEGHMIXCOEF') 
              IF (.not.allocated(VEGHMIXCOEF))                          &
     &                 allocate(VEGHMIXCOEF(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  VEGHMIXCOEF(iveg,ng)=Rveg(iveg,ng)
                END DO
              END DO
#endif 
#if defined MARSH_SED_EROSION
!           IF (.not.allocated(Rmarsh)) allocate(Rmarsh(Ngrids))
             CASE ('KFAC_MARSH')
               IF (.not.allocated(KFAC_MARSH))                          &
     &                allocate(KFAC_MARSH(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 KFAC_MARSH(ng)=Rmarsh(ng)
               END DO
# if defined MARSH_RETREAT
             CASE ('SCARP_HGHT')
               IF (.not.allocated(SCARP_HGHT))                          &
     &                allocate(SCARP_HGHT(Ngrids))
                 Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
                 DO ng=1,Ngrids
                   SCARP_HGHT(ng)=Rmarsh(ng)
                 END DO
# endif 
#endif
#if defined MARSH_TIDAL_RANGE_CALC
            CASE ('NTIMES_MARSH') 
              Npts=load_i(Nval, Rval, Ngrids, NTIMES_MARSH)
                IF (NTIMES_MARSH.lt.0) THEN
                  IF (Master) WRITE (out,30) 'NTIMES_MARSH', ng,        &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
#endif
#if defined MARSH_VERT_GROWTH 
            IF (.not.allocated(Rveg)) allocate(Rveg(NVEG,Ngrids)) 
             CASE ('PAR_FAC1')
               IF (.not.allocated(PAR_FAC1))                            &
     &                allocate(PAR_FAC1(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 PAR_FAC1(ng)=Rmarsh(ng)
               END DO
             CASE ('PAR_FAC2')
               IF (.not.allocated(PAR_FAC2))                            &
     &                allocate(PAR_FAC2(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 PAR_FAC2(ng)=Rmarsh(ng)
               END DO
             CASE ('TDAYS_MARSH_GROWTH')
               IF (.not.allocated(TDAYS_MARSH_GROWTH))                  &
     &                   allocate(TDAYS_MARSH_GROWTH(Ngrids))
                 Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
                 DO ng=1,Ngrids
                   TDAYS_MARSH_GROWTH(ng)=Rmarsh(ng)
                 END DO
                 IF (TDAYS_MARSH_GROWTH(ng).lt.0) THEN
                   IF (Master) WRITE (out,30) 'TDAYS_MARSH_GROWTH', ng,  &
     &                'must be greater than zero.'
                      exit_flag=5
                   RETURN
                 END IF
!             CASE ('MARSH_BULK_DENS')
!               IF (.not.allocated(MARSH_BULK_DENS))                     &
!     &                allocate(MARSH_BULK_DENS(Ngrids))
!               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
!               DO ng=1,Ngrids
!                 MARSH_BULK_DENS(ng)=Rmarsh(ng)
!               END DO
             CASE ('NUGP')
               IF (.not.allocated(NUGP))                                &
     &                allocate(NUGP(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 NUGP(ng)=Rmarsh(ng)
               END DO
             CASE ('BMAX')
               IF (.not.allocated(BMAX))                                &
     &                allocate(BMAX(Ngrids))
                 Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
                 DO ng=1,Ngrids
                   BMAX(ng)=Rmarsh(ng)
                 END DO
             CASE ('CHIREF')
               IF (.not.allocated(CHIREF))                              &
     &                allocate(CHIREF(Ngrids))
                 Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
                 DO ng=1,Ngrids
                   CHIREF(ng)=Rmarsh(ng)
                 END DO
# if defined MARSH_BIOMASS_VEG
             CASE ('ALPHA_PDENS')
               IF (.not.allocated(ALPHA_PDENS))                         &
     &                allocate(ALPHA_PDENS(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 ALPHA_PDENS(ng)=Rmarsh(ng)
               END DO
             CASE ('BETA_PDENS')
               IF (.not.allocated(BETA_PDENS))                          &
     &                allocate(BETA_PDENS(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 BETA_PDENS(ng)=Rmarsh(ng)
               END DO
             CASE ('ALPHA_PHGHT')
               IF (.not.allocated(ALPHA_PHGHT))                         &
     &                allocate(ALPHA_PHGHT(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 ALPHA_PHGHT(ng)=Rmarsh(ng)
               END DO
             CASE ('BETA_PHGHT')
               IF (.not.allocated(BETA_PHGHT))                          &
     &                allocate(BETA_PHGHT(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 BETA_PHGHT(ng)=Rmarsh(ng)
               END DO
             CASE ('ALPHA_PDIAM')
               IF (.not.allocated(ALPHA_PDIAM))                         &
     &                allocate(ALPHA_PDIAM(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 ALPHA_PDIAM(ng)=Rmarsh(ng)
               END DO
             CASE ('BETA_PDIAM')
               IF (.not.allocated(BETA_PDIAM))                          &
     &                allocate(BETA_PDIAM(Ngrids))
               Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
               DO ng=1,Ngrids
                 BETA_PDIAM(ng)=Rmarsh(ng)
               END DO
# endif 
#endif 
!
!-----------------------------------------------------------------------
!  Read output ids from vegetation.in
!-----------------------------------------------------------------------
!
#if defined VEG_DRAG || defined VEG_BIOMASS
            CASE ('Hout(ipdens)')
              IF (idvprp(pdens).eq.0) THEN
                IF (Master) WRITE (out,30) 'ipdens'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(pdens),:))
            CASE ('Hout(iphght)')
              IF (idvprp(phght).eq.0) THEN
                IF (Master) WRITE (out,30) 'iphght'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(phght),:))
            CASE ('Hout(ipdiam)')
              IF (idvprp(pdiam).eq.0) THEN
                IF (Master) WRITE (out,30) 'ipdiam'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(pdiam),:))
            CASE ('Hout(ipthck)')
              IF (idvprp(pthck).eq.0) THEN
                IF (Master) WRITE (out,30) 'ipthck'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(pthck),:))
#endif                                         
#ifdef VEG_STREAMING
            CASE ('Hout(idWdvg)')
              IF ((idWdvg).eq.0) THEN
                IF (Master) WRITE (out,30) 'idWdvg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdvg,:))
            CASE ('Hout(idCdvg)')
              IF ((idCdvg).eq.0) THEN
                IF (Master) WRITE (out,30) 'idCdvg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idCdvg,:))
#endif
#ifdef MARSH_DYNAMICS
            CASE ('Hout(idTims)')
              IF (idTims.eq.0) THEN
                IF (Master) WRITE (out,30) 'idTims'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTims,:))
# ifdef MARSH_WAVE_THRUST
            CASE ('Hout(idTtot)')
              IF (idTtot.eq.0) THEN
                IF (Master) WRITE (out,30) 'idTtot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTtot,:))
#  ifdef MARSH_SED_EROSION
            CASE ('Hout(idTmfo)')
              DO ng=1,Ngrids
                DO ised=1,NST
                  IF (idTmfo(ised).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idTmfo'
                    exit_flag=5
                    RETURN
                  END IF
                END DO
              END DO
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand1)
              DO ng=1,Ngrids
                DO ised=1,NST
                  i=idTmfo(ised)
                  Hout(i,ng)=Lsand1(ised,ng)
                END DO
              END DO
#   ifdef MARSH_RETREAT
            CASE ('Hout(idTmmr)')
              IF (idTmmr.eq.0) THEN 
                IF (Master) WRITE (out,30) 'idTmmr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmmr,:))
#   endif 
#  endif 
# endif 
# ifdef MARSH_TIDAL_RANGE_CALC 
            CASE ('Hout(idTmtr)')
              IF (idTmtr.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmtr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmtr,1:Ngrids))
# endif 
!
# ifdef MARSH_VERT_GROWTH 
            CASE ('Hout(idTmhw)')
              IF (idTmhw.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmhw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmhw,1:Ngrids))
            CASE ('Hout(idTmlw)')
              IF (idTmlw.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmlw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmlw,1:Ngrids))
            CASE ('Hout(idTmvg)')
              IF (idTmvg.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmvg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmvg,:))
            CASE ('Hout(idTmvt)')
              IF (idTmvt.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmvt'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmvt,:))
            CASE ('Hout(idTmbp)')
              IF (idTmbp.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmbp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmbp,1:Ngrids))
# endif 
#endif
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
! Print/Report input parameters (values specified in vegetation.in).
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
#if defined VEG_DRAG || defined VEG_BIOMASS 
            WRITE (out,50) ng
            WRITE (out,60)
            DO iveg=1,NVEG
              WRITE (out,70) NVEG, CD_VEG(iveg,ng), E_VEG(iveg,ng),     &
     &                       VEG_MASSDENS(iveg,ng), VEGHMIXCOEF(iveg,ng)
            END DO 
#endif 
#if defined MARSH_DYNAMICS
           WRITE (out,80)  ng
# if defined MARSH_SED_EROSION  
           WRITE (out,90)  KFAC_MARSH(ng)
# endif 
# if defined MARSH_RETREAT 
           WRITE (out,100)  SCARP_HGHT(ng)
# endif 
# ifdef MARSH_TIDAL_RANGE_CALC 
!           WRITE (out,110)
           WRITE (out,120) NTIMES_MARSH
# endif 
# ifdef MARSH_VERT_GROWTH
!	   WRITE(out,130) 
	   WRITE(out,130) PAR_FAC1(ng)
	   WRITE(out,140) PAR_FAC2(ng)
           WRITE(out,150) TDAYS_MARSH_GROWTH(ng)
           WRITE (out,160) NUGP(ng)
	   WRITE (out,170) BMAX(ng)
	   WRITE (out,180) CHIREF(ng)
#  ifdef MARSH_BIOMASS_VEG 
           WRITE (out,190)
           WRITE (out,200)
	   WRITE (out,210) ALPHA_PDENS(ng), BETA_PDENS(ng),             &
     &                     ALPHA_PHGHT(ng), BETA_PHGHT(ng),             & 
     &                     ALPHA_PDIAM(ng), BETA_PDIAM(ng) 
#  endif
# endif 	
#endif 
        END DO 
      ENDIF 
!
!-----------------------------------------------------------------------
!  Report output parameters (switched on in vegetation.in).
!-----------------------------------------------------------------------
! 
   30  FORMAT (/,' read_VegPar - variable info not yet loaded, ',a)
   40  FORMAT (/,' read_VegPar - Error while processing line: ',/,a)
#if defined VEG_DRAG || defined VEG_BIOMASS 
   50  FORMAT (/,/,' Vegetation Parameters, Grid: ',i2.2,               &
      &        /,  ' =====================================',/)
   60  FORMAT (/,1x,'Nveg(unitless)',2x,'Cd_veg(unitless)',2x,          &
      &        'E_veg(N/m2)',2x,'Veg_massdens(kg/m3)',2x,'VegHMixCoef'/)
   70  FORMAT (2x,i2,2(10x,1p,e11.4),2(5x,1p,e11.4))
#endif
#ifdef MARSH_DYNAMICS
   80  FORMAT (/,/,' Marsh dynamics model Parameters, Grid: ',i2.2,     &
      &        /,  ' =====================================',/)
# if defined MARSH_SED_EROSION
   90  FORMAT ('Marsh erosion coefficient (s/m)        = ',e11.3,/,a)
#  if defined MARSH_RETREAT
  100  FORMAT ('Marsh scarp height        (m)          = ',e11.3,/,a)
#  endif 
!  110  FORMAT (1x,l1,2x,a,t29,a,i2.2,':',1x,a)
# endif 
# ifdef MARSH_TIDAL_RANGE_CALC
  120  FORMAT ('Days after which MHW calc. starts     = ', i4,/,a)
# endif
# ifdef MARSH_VERT_GROWTH 
  130  FORMAT ('Parabolic growth factor 1              = ',e11.3,/,a)
  140  FORMAT ('Parabolic growth factor 2              = ',e11.3,/,a)
  150  FORMAT ('Number of growing days for marsh biomass  = ',e11.3,/,a)
!  160  FORMAT ('Marsh organic sed. bulk density (kg/m3)= ',e11.3,/,a)
!  130  FORMAT (/,1x,'par_fac1',5x,'par_fac2',7x,                        &
!      &        'tdays_marsh_growth(tdays)',3x,'marsh_bulk_dens(kg/m3)'/) 
!  140  FORMAT ((1x,1p,e11.4),(2x,1p,e11.4),(5x,1p,e11.3),(5x,1p,e11.4))
  160  FORMAT ('Fraction of below ground biomass       = ',e11.3,/,a)
  170  FORMAT ('Peak biomass (kg/m2)                   = ',e11.3,/,a)
  180  FORMAT ('Fraction of recalcitrant Carbon        = ',e11.3,/,a)
#  ifdef MARSH_BIOMASS_VEG
  190  FORMAT (/,'Marsh vegetation growth parameters: ',/,a)
  200  FORMAT (/,2x,'alpha_pdens', 4x,'beta_pdens',                     &
      &         4x, 'alpha_phght',4x,'beta_phght',                      &
      &         4x,'alpha_pdiam', 4x,'beta_pdiam'/)
  210  FORMAT (6(3x,1p,e11.3))
#  endif 
# endif 
#endif
      RETURN
      END SUBROUTINE read_VegPar

