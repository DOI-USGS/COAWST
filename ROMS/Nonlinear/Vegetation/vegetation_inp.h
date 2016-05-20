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
      integer :: iveg, ng, status
      integer :: decode_line, load_i, load_l, load_lbc, load_r
!
      real(r8), dimension(200) :: Rval
      real(r8), allocatable :: Rveg(:,:)
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
              Npts=load_r(Nval, Rval, NVEG*Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  CD_VEG(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
            CASE ('E_VEG') 
              IF (.not.allocated(E_VEG)) allocate(E_VEG(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG*Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  E_VEG(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
            CASE ('VEG_MASSDENS') 
              IF (.not.allocated(VEG_MASSDENS))                         &
     &                 allocate(VEG_MASSDENS(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG*Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  VEG_MASSDENS(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
            CASE ('VEGHMIXCOEF') 
              IF (.not.allocated(VEGHMIXCOEF))                         &
     &                 allocate(VEGHMIXCOEF(NVEG,Ngrids)) 
              Npts=load_r(Nval, Rval, NVEG*Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  VEGHMIXCOEF(iveg,ng)=Rveg(iveg,ng)
                END DO 
              END DO
!
!-----------------------------------------------------------------------
!  Read output ids from vegetation.in
!-----------------------------------------------------------------------
!          
#ifdef VEG_DRAG
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
#ifdef VEG_BIOMASS
            CASE ('Hout(ipupbm)')
              IF (idvprp(pupbm).eq.0) THEN 
                IF (Master) WRITE (out,30) 'ipupbm'
                exit_flag=5
                RETURN
              END IF 
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(pupbm),:))
            CASE ('Hout(ipdwbm)')
              IF (idvprp(pdwbm).eq.0) THEN 
                IF (Master) WRITE (out,30) 'ipdwbm'
                exit_flag=5
                RETURN
              END IF 
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(pdwbm),:))
#endif 
#ifdef VEG_STREAMING 
            CASE ('Hout(idWdvg)')
              IF ((idWdvg).eq.0) THEN 
                IF (Master) WRITE (out,30) 'idWdvg'
                exit_flag=5
                RETURN
              END IF 
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdvg,:))
#endif 
#ifdef MARSH_WAVE_THRUST
            CASE ('Hout(idTims)')
              IF (idTims.eq.0) THEN 
                IF (Master) WRITE (out,30) 'idTims'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTims,:))
            CASE ('Hout(idTmsk)')
              IF (idTmsk.eq.0) THEN 
                IF (Master) WRITE (out,30) 'idTmsk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmsk,:))
            CASE ('Hout(idTmax)')
              IF (idTmax.eq.0) THEN 
                IF (Master) WRITE (out,30) 'idTmax'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmax,:))
            CASE ('Hout(idTton)')
              IF (idTton.eq.0) THEN 
                IF (Master) WRITE (out,30) 'idTton'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTton,:))
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
            WRITE (out,50) ng
            WRITE (out,60)
            DO iveg=1,NVEG
              WRITE (out,70) NVEG, CD_VEG(iveg,ng), E_VEG(iveg,ng),     &
     &                       VEG_MASSDENS(iveg,ng), VEGHMIXCOEF(iveg,ng)
            END DO 
        END DO 
!     END IF 
!
!-----------------------------------------------------------------------
!  Report output parameters (switched on in vegetation.in).
!-----------------------------------------------------------------------
! 
#ifdef VEG_DRAG
!      IF (Hout(idvprp(pdens),ng)WRITE (out,90) Hout(idvprp(pdens),ng),  & 
!     &  'Hout(ipdens)',                                                 &
!     &  'Write out plant density for each veg type.'
!      IF (Hout(idvprp(pdiam),:)) WRITE (out,90) Hout(idvprp(pdiam),ng),  & 
!     &  'Hout(ipdiam)',                                                 &
!     &  'Write out plant diameter for each veg type.'
!      IF (Hout(idvprp(pthck),ng)) WRITE (out,90) Hout(idvprp(pthck),:), & 
!     &  'Hout(ipthck)',                                                 &
!     &  'Write out plant thickness for each veg type.'
#endif 
#ifdef VEG_BIOMASS
!      IF (Hout(idvprp(pupbm),ng) WRITE (out,90) Hout(idvprp(pupbm),ng), & 
!     &  'Hout(ipupbm)',                                                 &
!     &  'Write out above ground plant biomass.'
!      IF (Hout(idvprp(pdwbm),ng) WRITE (out,90) Hout(idvprp(pdwbm),ng), & 
!     &  'Hout(idwbm)',                                                  &
!     &  'Write out below ground plant biomass.'
#endif 
#ifdef VEG_STREAMING 
!      IF (Hout(idWdvg,ng) WRITE (out,90) Hout(idWdvg,ng), & 
!     &  'Hout(idWdvg).    
!     &  'Write out wave dissipation due to vegetation.'
#endif
#ifdef MARSH_WAVE_THRUST
          IF (Hout(idTims,ng)) WRITE (out,90) Hout(idTims,ng),          &
     &       'Hout(idTims)',                                            &
     &       'Write out initial masking for marshes.'
          IF (Hout(idTmsk,ng)) WRITE (out,90) Hout(idTmsk,ng),          &
     &       'Hout(idTmsk)',                                            &
     &       'Write out masking for getting thrust due to waves.'
          IF (Hout(idTmax,ng)) WRITE (out,90) Hout(idTmax,ng),          &
     &       'Hout(idTmax)',                                            &
     &       'Write out maximum thrust due to waves.'
         IF (Hout(idTton,ng)) WRITE (out,90) Hout(idTton,ng),           &
     &       'Hout(idTton)',                                            &
     &       'Write out Tonelli masking based thrust due to waves.'
#endif 
       END IF 
   30  FORMAT (/,' read_VegPar - variable info not yet loaded, ',a)
   40  FORMAT (/,' read_VegPar - Error while processing line: ',/,a)
   50  FORMAT (/,/,' Vegetation Parameters, Grid: ',i2.2,               &
      &        /,  ' =================================',/)
   60  FORMAT (/,1x,'Nveg(unitless)',2x,'Cd_veg(unitless)',2x,          &
      &        'E_veg(N/m2)',2x,'Veg_massdens(kg/m3)',2x,'VegHMixCoef'/)
   70  FORMAT (2x,i2,2(10x,1p,e11.4),2(5x,1p,e11.4))
   80  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
!   90  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
   90  FORMAT (10x,l1,2x,a,t32,a)
        
  100  FORMAT (10x,l1,2x,a,t32,a,1x,a)

      RETURN
      END SUBROUTINE read_VegPar

