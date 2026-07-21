      SUBROUTINE read_VegPar (model, inp, out, Lwrite)
!
!git $Id$
!=======================================================================
!                                                                      !
!  This routine reads in submerged aquatic vegetation model parameters !
!  from its standard input file, usually "vegetation.in".              !
!                                                                      !
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
      logical, dimension(NST,Ngrids) :: Lsed
      logical, dimension(Ngrids) :: Lswitch
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status
      integer :: iveg, ised, it
      integer :: Ivalue(1)
!
      real(r8), dimension(200) :: Rval
#if defined MARSH_SED_EROSION  || defined MARSH_VERT_GROWTH
      real(r8), dimension(Ngrids) :: Rmarsh
#endif
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
#ifdef VEG_DRAG
            CASE ('NVEG')
              Npts=load_i(Nval, Rval, 1, NVEG)
              IF (NVEG.lt.0) THEN
                IF (Master) WRITE (out,30) 'NVEG', ng,                  &
     &                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              IF (.not.allocated(Rveg)) THEN
                allocate ( Rveg(NVEG,Ngrids) )    ! local variable
              END IF
              DO ng=1,Ngrids                      ! allocate variables
                CALL allocate_vegetation (ng, 0)  ! solely depending on
              END DO                              ! NVEG and Ngrids
            CASE ('CD_VEG')
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  CD_VEG(iveg,ng)=Rveg(iveg,ng)
                END DO
              END DO
            CASE ('E_VEG')
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  E_VEG(iveg,ng)=Rveg(iveg,ng)
                END DO
              END DO
            CASE ('VEG_MASSDENS')
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  VEG_MASSDENS(iveg,ng)=Rveg(iveg,ng)
                END DO
              END DO
            CASE ('VEGHMIXCOEF')
              Npts=load_r(Nval, Rval, NVEG, Ngrids, Rveg)
              DO ng=1,Ngrids
                DO iveg=1,NVEG
                  VEGHMIXCOEF(iveg,ng)=Rveg(iveg,ng)
                END DO
              END DO
#endif
#if defined MARSH_SED_EROSION
            CASE ('KFAC_MARSH')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                KFAC_MARSH(ng)=Rmarsh(ng)
              END DO
# if defined MARSH_RETREAT
            CASE ('SCARP_HGHT')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                SCARP_HGHT(ng)=Rmarsh(ng)
              END DO
# endif
#endif
#if defined MARSH_TIDAL_RANGE
            CASE ('NTIMES_MARSH')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NTIMES_MARSH=Ivalue(1)
              IF (NTIMES_MARSH.lt.0) THEN
                IF (Master) WRITE (out,30) 'NTIMES_MARSH', ng,          &
     &                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
#endif
#if defined MARSH_VERT_GROWTH
            CASE ('PAR_FAC1')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                PAR_FAC1(ng)=Rmarsh(ng)
              END DO
            CASE ('PAR_FAC2')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                 PAR_FAC2(ng)=Rmarsh(ng)
              END DO
            CASE ('TDAYS_MARSH_GROWTH')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                TDAYS_MARSH_GROWTH(ng)=Rmarsh(ng)
                IF (TDAYS_MARSH_GROWTH(ng).lt.0) THEN
                  IF (Master) WRITE (out,30) 'TDAYS_MARSH_GROWTH', ng,  &
     &                        'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
!!          CASE ('MARSH_BULK_DENS')
!!            Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
!!            DO ng=1,Ngrids
!!              MARSH_BULK_DENS(ng)=Rmarsh(ng)
!!            END DO
            CASE ('NUGP')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                NUGP(ng)=Rmarsh(ng)
              END DO
            CASE ('BMAX')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                BMAX(ng)=Rmarsh(ng)
              END DO
            CASE ('CHIREF')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                CHIREF(ng)=Rmarsh(ng)
              END DO
# if defined MARSH_BIOMASS_VEG
            CASE ('ALPHA_PDENS')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                 ALPHA_PDENS(ng)=Rmarsh(ng)
              END DO
            CASE ('BETA_PDENS')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                BETA_PDENS(ng)=Rmarsh(ng)
              END DO
            CASE ('ALPHA_PHGHT')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                ALPHA_PHGHT(ng)=Rmarsh(ng)
              END DO
            CASE ('BETA_PHGHT')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                BETA_PHGHT(ng)=Rmarsh(ng)
              END DO
            CASE ('ALPHA_PDIAM')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                ALPHA_PDIAM(ng)=Rmarsh(ng)
              END DO
            CASE ('BETA_PDIAM')
              Npts=load_r(Nval, Rval, Ngrids, Rmarsh)
              DO ng=1,Ngrids
                BETA_PDIAM(ng)=Rmarsh(ng)
              END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Read in output logical switches.
!-----------------------------------------------------------------------
!
#if defined VEG_DRAG || defined VEG_BIOMASS
            CASE ('Hout(isDens)', 'Hout(ipdens)')
              IF (idvprp(isDens).eq.0) THEN
                IF (Master) WRITE (out,30) 'idvprp(isDens)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(isDens),:))
            CASE ('Hout(isHght)', 'Hout(iphght)')
              IF (idvprp(isHght).eq.0) THEN
                IF (Master) WRITE (out,30) 'idvprp(isHght)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(isHght),:))
            CASE ('Hout(isDiam)', 'Hout(ipdiam)')
              IF (idvprp(isDiam).eq.0) THEN
                IF (Master) WRITE (out,30) 'idvprp(isDiam)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(isDiam),:))
            CASE ('Hout(isThck)', 'Hout(ipthck)')
              IF (idvprp(isThck).eq.0) THEN
                IF (Master) WRITE (out,30) 'idvprp(isThck)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idvprp(isThck),:))
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
              Npts=load_l(Nval, Cval, NST, Ngrids, Lsed)
              DO ng=1,Ngrids
                DO ised=1,NNS
                  i=idTmfo(ised)
                  Hout(i,ng)=Lsed(ised,ng)
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
# ifdef MARSH_TIDAL_RANGE
            CASE ('Hout(idTmtr)')
              IF (idTmtr.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmtr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmtr,:))
# endif
# ifdef MARSH_VERT_GROWTH
            CASE ('Hout(idTmhw)')
              IF (idTmhw.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmhw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmhw,:))
            CASE ('Hout(idTmlw)')
              IF (idTmlw.eq.0) THEN
                IF (Master) WRITE (out,40) 'idTmlw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmlw,:))
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
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTmbp,:))
# endif
#endif
#if defined VEG_DRAG || defined VEG_BIOMASS
            CASE ('Qout(isDens)', 'Qout(ipdens)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idvprp(isDens),:))
            CASE ('Qout(isHght)', 'Qout(iphght)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idvprp(isHght),:))
            CASE ('Qout(isDiam)', 'Qout(ipdiam)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idvprp(isDiam),:))
            CASE ('Qout(isThck)', 'Qout(ipthck)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idvprp(isThck),:))
#endif
#ifdef VEG_STREAMING
            CASE ('Qout(idWdvg)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idWdvg,:))
            CASE ('Qout(idCdvg)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idCdvg,:))
#endif
#ifdef MARSH_DYNAMICS
            CASE ('Qout(idTims)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTims,:))
# ifdef MARSH_WAVE_THRUST
            CASE ('Qout(idTtot)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTtot,:))
#  ifdef MARSH_SED_EROSION
            CASE ('Qout(idTmfo)')
              Npts=load_l(Nval, Cval, NST, Ngrids, Lsed)
              DO ng=1,Ngrids
                DO ised=1,NST
                  i=idTmfo(ised)
                  Qout(i,ng)=Lsed(ised,ng)
                END DO
              END DO
#   ifdef MARSH_RETREAT
            CASE ('Qout(idTmmr)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmmr,:))
#   endif
#  endif
# endif
# ifdef MARSH_TIDAL_RANGE
            CASE ('Qout(idTmtr)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmtr,:))
# endif
# ifdef MARSH_VERT_GROWTH
            CASE ('Qout(idTmhw)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmhw,:))
            CASE ('Qout(idTmlw)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmlw,:))
            CASE ('Qout(idTmvg)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmvg,:))
            CASE ('Qout(idTmvt)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmvt,:))
            CASE ('Qout(idTmbp)')
              Npts=load_l(Nval, Cval, Ngrids, Qout(idTmbp,:))
# endif
#endif
#if defined DIAGNOSTICS_UV && defined VEG_DRAG
            CASE ('Dout(M3fveg)')
              IF (M3fveg.le.0) THEN
                IF (Master) WRITE (out,40) 'M3fveg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              DO ng=1,Ngrids
                Dout(idDu3d(M3fveg),ng)=Lswitch(ng)
                Dout(idDv3d(M3fveg),ng)=Lswitch(ng)
              END DO
            CASE ('Dout(M2fveg)')
              IF (M2fveg.le.0) THEN
                IF (Master) WRITE (out,40) 'M2fveg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lswitch)
              DO ng=1,Ngrids
                Dout(idDu2d(M2fveg),ng)=Lswitch(ng)
                Dout(idDv2d(M2fveg),ng)=Lswitch(ng)
              END DO
#endif
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
      IF (allocated(Rveg)) deallocate (Rveg)
!
!-----------------------------------------------------------------------
! Print/Report input parameters (values specified in vegetation.in).
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          WRITE (out,50) ng
#if defined VEG_DRAG || defined VEG_BIOMASS
          WRITE (out,60) 'CD_VEG',                                      &
     &          'Plant drag coefficient (nondimensional)'
          DO iveg=1,NVEG
            WRITE (out,70) CD_VEG(iveg,ng), iveg
          END DO
          WRITE (out,60) 'E_VEG',                                       &
     &          'Plant Young modulus: stress/tension (nondimensional)'
          DO iveg=1,NVEG
            WRITE (out,70) E_VEG(iveg,ng), iveg
          END DO
          WRITE (out,60) 'VEG_MassDens',                                &
     &          'Plant stems mass density (kg/m3)'
          DO iveg=1,NVEG
            WRITE (out,70) VEG_MASSDENS(iveg,ng), iveg
          END DO
          WRITE (out,60) 'VegHmixCoef',                                 &
     &          'Plant patch horizontal viscosity (m2/s)'
          DO iveg=1,NVEG
            WRITE (out,70) VEGHMIXCOEF(iveg,ng), iveg
          END DO
#endif
#if defined MARSH_DYNAMICS
# if defined MARSH_SED_EROSION
          WRITE (out,80) KFAC_MARSH(ng), 'KFAC_MARCH',                  &
     &          'Marsh sediment erodibility coefficient (m2/s)'
# endif
# if defined MARSH_RETREAT
          WRITE (out,80) SCARP_HGHT(ng), 'SCARP_HGHT',                  &
     &          'Scarp elevation (m) criteria for marsh retreat'
# endif
# ifdef MARSH_TIDAL_RANGE
          WRITE (out,90) NTIMES_MARSH, 'NTIMES_MARSH',                  &
     &          'Number of days for mean higher high water (MHHW)'
# endif
# ifdef MARSH_VERT_GROWTH
          WRITE (out,80) PAR_FAC1(ng), 'PAR_FAC1',                      &
     &          'Marsh parabolic curve growth biomass parameter 1'
          WRITE (out,80) PAR_FAC2(ng), 'PAR_FAC2',                      &
     &          'Marsh parabolic curve growth biomass parameter 2'
          WRITE (out,80) TDAYS_MARSH_GROWTH(ng), 'TDAYS_MARSH_GROWTH',  &
     &          'Number of growing days for marsh biomass'
          WRITE (out,80) NUGP(ng), 'NuGP',                              &
     &          'Fraction of below ground biomass (1/days)'
          WRITE (out,80) BMAX(ng), 'Bmax',                              &
     &          'Optimal marsh biomass growth rate (kg/m2/year)'
          WRITE (out,80) CHIREF(ng), 'ChiRef',                          &
     &          'Fraction of recalcitrant Carbon in effective biomass'
#  ifdef MARSH_BIOMASS_VEG
          WRITE (out,80) ALPHA_PDENS(ng), 'alpha_pdens',                &
     &          'Marsh alpha growth parameter for plant density'
          WRITE (out,80) BETA_PDENS(ng), 'beta_pdens',                  &
     &          'Marsh beta  growth parameter for plant density'
          WRITE (out,80) ALPHA_PHGHT(ng), 'alpha_phght',                &
     &          'Marsh alpha growth parameter for plant height'
          WRITE (out,80) BETA_PHGHT(ng), 'beta_pdens',                  &
     &          'Marsh beta  growth parameter for plant height'
          WRITE (out,80) ALPHA_PDIAM(ng), 'alpha_pdiam',                &
     &          'Marsh alpha growth parameter for plant diameter'
          WRITE (out,80) BETA_PDIAM(ng), 'beta_pdens',                  &
     &          'Marsh beta  growth parameter for plant diameter'
#  endif
# endif
#endif
#if defined VEG_DRAG || defined VEG_BIOMASS
          IF (Hout(idvprp(isDens),ng)) WRITE (out,100)                  &
     &        Hout(idvprp(isDens),ng), 'Hout(isDens)',                  &
     &        'Plant density, individuals per area'
          IF (Hout(idvprp(isHght),ng)) WRITE (out,100)                  &
     &        Hout(idvprp(isHght),ng), 'Hout(isHght)',                  &
     &        'Dominant plant mean height'
          IF (Hout(idvprp(isDiam),ng)) WRITE (out,100)                  &
     &        Hout(idvprp(isDiam),ng), 'Hout(isDiam)',                  &
     &        'Dominant plant mean diameter'
          IF (Hout(idvprp(isThck),ng)) WRITE (out,100)                  &
     &        Hout(idvprp(isThck),ng), 'Hout(isThck)',                  &
     &        'Dominant plant mean thickness'
#endif
#ifdef VEG_STREAMING
          IF (Hout(idWdvg,ng)) WRITE (out,100)                          &
     &        Hout(idWdvg,ng), 'Hout(idWdvg)',                          &
     &        'Wave dissipation due to vegetation'
          IF (Hout(idCdvg,ng)) WRITE (out,100)                          &
     &        Hout(idCdvg,ng), 'Hout(idCdvg)',                          &
     &        'Spectral drag coefficient due to waves and vegetation'
#endif
#ifdef MARSH_DYNAMICS
          IF (Hout(idTims,ng)) WRITE (out,100)                          &
     &        Hout(idTims,ng), 'Hout(idTims)',                          &
     &        'Marsh cell mask'
# ifdef MARSH_WAVE_THRUST
          IF (Hout(idTtot,ng)) WRITE (out,100)                          &
     &        Hout(idTtot,ng), 'Hout(idTtot)',                          &
     &        'Total wave thrust on marsh cells'
#  ifdef MARSH_SED_EROSION
          DO i=1,NST
            itrc=idsed(i)
            IF (Hout(idTmfo(i),ng)) WRITE (out,110)                     &
     &        Hout(idTmfo(i),ng), 'Hout(idTmfo)',                       &
     &        'Sediment flux out of marsh cells, tracer ', itrc,        &
     &        TRIM(Vname(1,idTvar(itrc)))
          END DO
#   ifdef MARSH_RETREAT
          IF (Hout(idTmmr,ng)) WRITE (out,100)                          &
     &        Hout(idTmmr,ng), 'Hout(idTmmr)',                          &
     &        'Amount of marsh retreat'
#   endif
#  endif
# endif
# ifdef MARSH_TIDAL_RANGE
          IF (Hout(idTmtr,ng)) WRITE (out,100)                          &
     &        Hout(idTmtr,ng), 'Hout(idTmtr)',                          &
     &        'Tidal range for marsh growth'
# endif
# ifdef MARSH_VERT_GROWTH
          IF (Hout(idTmhw,ng)) WRITE (out,100)                          &
     &        Hout(idTmhw,ng), 'Hout(idTmhw)',                          &
     &        'Mean high water for marsh cells'
          IF (Hout(idTmlw,ng)) WRITE (out,100)                          &
     &        Hout(idTmlw,ng), 'Hout(idTmlw)',                          &
     &        'Mean low water for marsh cells'
          IF (Hout(idTmvg,ng)) WRITE (out,100)                          &
     &        Hout(idTmvg,ng), 'Hout(idTmvg)',                          &
     &        'Rate of marsh vertical growth'
          IF (Hout(idTmbp,ng)) WRITE (out,100)                          &
     &        Hout(idTmbp,ng), 'Hout(idTmbp)',                          &
     &        'Marsh biomass peak production'
# endif
#endif
#if defined VEG_DRAG || defined VEG_BIOMASS
          IF (Qout(idvprp(isDens),ng)) WRITE (out,100)                  &
     &        Qout(idvprp(isDens),ng), 'Qout(isDens)',                  &
     &        'Plant density, individuals per area'
          IF (Qout(idvprp(isHght),ng)) WRITE (out,100)                  &
     &        Qout(idvprp(isHght),ng), 'Qout(isHght)',                  &
     &        'Dominant plant mean height'
          IF (Qout(idvprp(isDiam),ng)) WRITE (out,100)                  &
     &        Qout(idvprp(isDiam),ng), 'Qout(isDiam)',                  &
     &        'Dominant plant mean diameter'
          IF (Qout(idvprp(isThck),ng)) WRITE (out,100)                  &
     &        Qout(idvprp(isThck),ng), 'Qout(isThck)',                  &
     &        'Dominant plant mean thickness'
#endif
#ifdef VEG_STREAMING
          IF (Qout(idWdvg,ng)) WRITE (out,100)                          &
     &        Qout(idWdvg,ng), 'Qout(idWdvg)',                          &
     &        'Wave dissipation due to vegetation'
          IF (Qout(idCdvg,ng)) WRITE (out,100)                          &
     &        Qout(idCdvg,ng), 'Qout(idCdvg)',                          &
     &        'Spectral drag coefficient due to waves and vegetation'
#endif
#ifdef MARSH_DYNAMICS
          IF (Qout(idTims,ng)) WRITE (out,100)                          &
     &        Qout(idTims,ng), 'Qout(idTims)',                          &
     &        'Marsh cell mask'
# ifdef MARSH_WAVE_THRUST
          IF (Qout(idTtot,ng)) WRITE (out,100)                          &
     &        Qout(idTtot,ng), 'Qout(idTtot)',                          &
     &        'Total wave thrust on marsh cells'
#  ifdef MARSH_SED_EROSION
          DO i=1,NST
            itrc=idsed(i)
            IF (Qout(idTmfo(i),ng)) WRITE (out,110)                     &
     &        Qout(idTmfo(i),ng), 'Qout(idTmfo)',                       &
     &        'Sediment flux out of marsh cells, tracer ', itrc,        &
     &        TRIM(Vname(1,idTvar(itrc)))
          END DO
#   ifdef MARSH_RETREAT
          IF (Qout(idTmmr,ng)) WRITE (out,100)                          &
     &        Qout(idTmmr,ng), 'Qout(idTmmr)',                          &
     &        'Amount of marsh retreat'
#   endif
#  endif
# endif
# ifdef MARSH_TIDAL_RANGE
          IF (Qout(idTmtr,ng)) WRITE (out,100)                          &
     &        Qout(idTmtr,ng), 'Qout(idTmtr)',                          &
     &        'Tidal range for marsh growth'
# endif
# ifdef MARSH_VERT_GROWTH
          IF (Qout(idTmhw,ng)) WRITE (out,100)                          &
     &        Qout(idTmhw,ng), 'Qout(idTmhw)',                          &
     &        'Mean high water for marsh cells'
          IF (Qout(idTmlw,ng)) WRITE (out,100)                          &
     &        Qout(idTmlw,ng), 'Qout(idTmlw)',                          &
     &        'Mean low water for marsh cells'
          IF (Qout(idTmvg,ng)) WRITE (out,100)                          &
     &        Qout(idTmvg,ng), 'Qout(idTmvg)',                          &
     &        'Rate of marsh vertical growth'
          IF (Qout(idTmbp,ng)) WRITE (out,100)                          &
     &        Qout(idTmbp,ng), 'Qout(idTmbp)',                          &
     &        'Marsh biomass peak production'
# endif
#endif
#if defined DIAGNOSTICS_UV && defined VEG_DRAG
          IF (Dout(idDu2d(M2fveg),ng)) WRITE (out,100)                  &
     &        Dout(idDu2d(M2fveg),ng), 'Dout(M2fveg)',                  &
     &        '2D u-momentum drag force term due to vegetation'
          IF (Dout(idDv2d(M2fveg),ng)) WRITE (out,100)                  &
     &        Dout(idDv2d(M2fveg),ng), 'Dout(M2fveg)',                  &
     &        '2D v-momentum drag force term due to vegetation'
          IF (Dout(idDu3d(M3fveg),ng)) WRITE (out,100)                  &
     &        Dout(idDu3d(M3fveg),ng), 'Dout(M3fveg)',                  &
     &        '3D u-momentum drag force term due to vegetation'
          IF (Dout(idDv3d(M3fveg),ng)) WRITE (out,100)                  &
     &        Dout(idDv3d(M3fveg),ng), 'Dout(M3fveg)',                  &
     &        '3D v-momentum drag force term due to vegetation'
#endif
        END DO
      ENDIF
!
  30  FORMAT (/,' read_VegPar - variable info not yet loaded, ',a)
  40  FORMAT (/,' read_VegPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' Vegetation Model Parameters, Grid: ',i2.2,          &
     &        /,  ' =====================================',/)
  60  FORMAT ('...........',2x,a,t34,a)
  70  FORMAT (1p,e11.4,t35,'vegetation type ',i0)
  80  FORMAT (1p,e11.4,2x,a,t34,a)
  90  FORMAT (1x,i10,2x,a,t34,a)
 100  FORMAT (10x,l1,2x,a,t34,a,1x,a)
 110  FORMAT (10x,l1,2x,a,t34,a,i2.2,':',1x,a)
!
      RETURN
      END SUBROUTINE read_VegPar
