#undef GROWTH_ONLY

      SUBROUTINE biology_floats (ng, Lstr, Lend, Predictor, my_thread)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Diego A. Narvaez   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine sets behavior for Lagrangian particles that simulates  !
!  planktonic larvae.  The  planktonic behavior is based on the model  !
!  of Dekshenieks et al. (1993; 1996; 1997).  It  calculates the size  !
!  (length) and development of oyster larvae. Results about this model !
!  can be found in Narvaez et al. 2012 a,b (in review)                 !
!                                                                      !
!  The governing equation are:                                         !
!                                                                      !
!    d(Lsize)/d(t) = growth(food,Lsize) * Gfactor(T,S) * turb_ef       !
!                                                                      !
!    w_bio = TS * SW - (1 - TS) * SR                                   !
!                                                                      !
!  where                                                               !
!                                                                      !
!      turbef = m * turb + c,                for turbidity < 0.1 g/l   !
!  or                                                                  !
!      turbef = b * EXP(-beta*(turb-turb0)), for turbidity > 0.1 g/l   !
!                                                                      !
!      TS =  c * DS + d,    for increasing salinity gradient DS        !
!  or                                                                  !
!      TS = -e * DS + f,    for decreasing salinity gradient DS        !
!                                                                      !
!      SR = 2.665 * EXP(0.0058*(Lsize-220))                            !
!                                                                      !
!  TS: fraction of active larvae: fractional swimming time             !
!  SW: larval swimming rate (mm/s)                                     !
!  SR: larval sinking rate (mm/s)                                      !
!  DS: salinity change rate (1/s)                                      !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Dekshenieks, M.M., E.E. Hofmann, and E.N. Powell, 1993:           !
!      Environmental effect on the growth and development of           !
!      Eastern oyster, Crassostrea virginica (Gmelin, 1791),           !
!      larvae: A modeling study, J. Shelfish Res, 12, 241-254.         !
!                                                                      !
!    Dekshenieks, M.M., E.E. Hofmann, J.M. Klinck, and E.N.            !
!      Powell, 1996: Modeling the vertical distribution of             !
!      oyster larvae in response to environmental conditions,          !
!      Mar. Ecol. Prog. Ser., 136, 97-110.                             !
!                                                                      !
!    Dekshenicks, M.M., E.E. Hofmann, J.M. Klinck, and E.N.            !
!      Powell, 1997: A modeling study of the effect of size-           !
!      and depth_dependent predation on larval survival, J.            !
!      Plankton Res., 19, (11), 1583-1598.                             !
!                                                                      !
!    Narvaez, D.A, J.M. Klinck, E.N. Powell, E.E. Hofmann, J. Wilkin   !
!      and D. B. Haidvogel. Modeling the dispersal of Eastern          !
!      oyster (Crassostrea virginica) larvae in Delaware Bay, J. Mar.  !
!      Res., in review.                                                !
!                                                                      !
!    Narvaez, D.A, J.M. Klinck, E.N. Powell, E.E. Hofmann, J. Wilkin   !
!      and D. B. Haidvogel. Circulation and behavior controls on       !
!      dispersal of Eastern oyster (Crassostrea virginica) larvae in   !
!      Delaware Bay, J. Mar. Res., in review.                          !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend

      logical, intent(in) :: Predictor
#ifdef ASSUMED_SHAPE
      logical, intent(in) :: my_thread(Lstr:)
#else
      logical, intent(in) :: my_thread(Lstr:Lend)
#endif
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
#endif
      CALL biology_floats_tile (ng, Lstr, Lend,                         &
     &                          nfm3(ng), nfm2(ng), nfm1(ng), nf(ng),   &
     &                          nfp1(ng),                               &
     &                          Predictor, my_thread,                   &
     &                          DRIFTER(ng) % bounded,                  &
     &                          DRIFTER(ng) % Tinfo,                    &
     &                          DRIFTER(ng) % track)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
#endif

      RETURN
      END SUBROUTINE biology_floats
!
!***********************************************************************
      SUBROUTINE biology_floats_tile (ng, Lstr, Lend,                   &
     &                                nfm3, nfm2, nfm1, nf, nfp1,       &
     &                                Predictor, my_thread, bounded,    &
     &                                Tinfo, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_behavior
#ifdef BIOLOGY
      USE mod_biology
#endif
      USE mod_floats
      USE mod_grid
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: nfm3, nfm2, nfm1, nf, nfp1
      logical, intent(in) :: Predictor
!
#ifdef ASSUMED_SHAPE
      logical, intent(in) :: bounded(:)
      logical, intent(in) :: my_thread(Lstr:)

      real(r8), intent(in) :: Tinfo(0:,:)
      real(r8), intent(inout) :: track(:,0:,:)
#else
      logical, intent(in) :: bounded(Nfloats(ng))
      logical, intent(in) :: my_thread(Lstr:Lend)

      real(r8), intent(in) :: Tinfo(0:izrhs,Nfloats(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
#endif
!
!  Local variable declarations.
!
      integer :: i, i1, i2, j1, j2, l

      real(r8) :: dsalt, temp, salt
      real(r8) :: Lfood, Lturb, Lsize, LsizeNew
      real(r8) :: Grate, Gfactor, turb_ef
      real(r8) :: SwimRate, SwimTime, SwimTimeNew
      real(r8) :: bottom, brhs, sink, w_bio
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8) :: p1, p2, q1, q2
      real(r8) :: my_food, my_salt, my_size, my_temp
      real(r8) :: oGfactor_DS, oGfactor_DT
      real(r8) :: oGrate_DF, oGrate_DL
      real(r8) :: oswim_DL, oswim_DT
      real(r8) :: HalfDT
!
!-----------------------------------------------------4-----------------
!  Estimate larval growth, as length (um), based on food, salinity,
!  temperature and turbidity.  Then, estimate swimming time (s),
!  larvae sinking velocity (m/s) and larvae vertical velocity (m/s).
!-----------------------------------------------------------------------
!
!  Compute look tables inverse axis increments to avoid repetitive
!  divisions.
!
      oGrate_DF=1.0_r8/Grate_DF
      oGrate_DL=1.0_r8/Grate_DL

      oGfactor_DS=1.0_r8/Gfactor_DS
      oGfactor_DT=1.0_r8/Gfactor_DT

      oswim_DT=1.0_r8/swim_DT
      oswim_DL=1.0_r8/swim_DL
!
!  Assign predictor/corrector weights.
!
      IF (Predictor) THEN
        cff1=8.0_r8/3.0_r8
        cff2=4.0_r8/3.0_r8
      ELSE
        cff1=9.0_r8/8.0_r8
        cff2=1.0_r8/8.0_r8
        cff3=3.0_r8/8.0_r8
        cff4=6.0_r8/8.0_r8
      END IF
!
!  Compute biological behavior fields.
!
      HalfDT=0.5_r8*dt(ng)

      DO l=Lstr,Lend
        IF (my_thread(l).and.bounded(l)) THEN
!
!  If newly relased float, initialize biological behavior fields. Note
!  that since we need temperature and salinity, we need to initialize
!  their values to all time levels. Otherwise, we will have a parallel
!  bug.
!
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                    &
     &        time(ng)+HalfDT.gt.Tinfo(itstr,l)) THEN
            temp=track(ifTvar(itemp),nfp1,l)
            salt=track(ifTvar(isalt),nfp1,l)
            DO i=0,NFT
              track(isizf,i,l)=Larvae_size0(ng)
              track(iswim,i,l)=0.5_r8*(swim_Tmin(ng)+swim_Tmax(ng))
              track(ifTvar(itemp),i,l)=temp
              track(ifTvar(isalt),i,l)=salt
            END DO
          END IF
!
!  Get temperature, salinity, larvae size (length), swimming time, and
!  food supply. For now, assume constant food and turbidity.  This
!  can be changed in the future with spatial and temporal variability
!  for food and/or turbidity.  If this is the case, we need to have
!  something like:
!
!         IF (Predictor) THEN
!           Lfood=track(ifood,nf,l)
!           Lturb=track(iturb,nf,l)
!           ...
!         ELSE
!           Lfood=track(ifood,nfp1,l)
!           Lturb=track(iturb,nfp1,l)
!           ...
!         END IF
!
!  The food variability may be from data, ecosystem model, or analytical
!  functions. Similarly, the turbidity may be from data, sediment model,
!  or analytical functions.
!
          temp=track(ifTvar(itemp),nfp1,l)
          salt=track(ifTvar(isalt),nfp1,l)
          dsalt=track(ifTvar(isalt),nfp1,l)-                            &
     &          track(ifTvar(isalt),nf  ,l)
          IF (Predictor) THEN
            Lsize=track(isizf,nf,l)
            SwimTime=track(iswim,nf,l)
            Lfood=food_supply(ng)
            Lturb=turb_ambi(ng)
          ELSE
            Lsize=track(isizf,nfp1,l)
            SwimTime=track(iswim,nfp1,l)
            Lfood=food_supply(ng)
            Lturb=turb_ambi(ng)
          END IF
!
!  Determine larval growth rate (um/day) contribution as function of
!  food supply (mg_Carbon/l) and larval size (um). Linearly interpolate
!  growth rate from the look table of food concentration versus larval
!  size. Notice that extrapolation is suppresed by bounding "Lfood" and
!  "Lsize" to the range values in the look table.
!
          my_food=MIN(MAX(Grate_F0,Lfood),                              &
     &                Grate_F0+Grate_DF*REAL(Grate_Im-1,r8))
          my_size=MIN(MAX(Grate_L0,Lsize),                              &
     &                Grate_L0+Grate_DL*REAL(Grate_Jm-1,r8))

          i1=INT(1.0_r8+(my_food-Grate_F0)*oGrate_DF)
          i2=MIN(i1+1,Grate_Im)
          j1=INT(1.0_r8+(my_size-Grate_L0)*oGrate_DL)
          j2=MIN(j1+1,Grate_Jm)

          p2=(my_food-(Grate_F0+REAL(i1-1,r8)*Grate_DF))*oGrate_DF
          q2=(my_size-(Grate_L0+REAL(j1-1,r8)*Grate_DL))*oGrate_DL
          p1=1.0_r8-p2
          q1=1.0_r8-q2

          Grate=p1*q1*Grate_table(i1,j1)+                               &
     &          p2*q1*Grate_table(i2,j1)+                               &
     &          p1*q2*Grate_table(i1,j2)+                               &
     &          p2*q2*Grate_table(i2,j2)
!
!  Determine larval growth rate factor (nondimensional) as function of
!  salinity and temperature (Celsius). Linearly interpolate growth rate
!  factor from the look table of salinity versus temperature. Notice
!  that extrapolation is suppresed by bounding "salt" and "temp" to the
!  range values in the look table.
!
          IF (temp.lt.Gfactor_T0) THEN
            Gfactor=0.0_r8
          ELSE
            my_salt=MIN(MAX(Gfactor_S0,salt),                           &
     &                  Gfactor_S0+Gfactor_DS*REAL(Gfactor_Im-1,r8))
            my_temp=MIN(MAX(Gfactor_T0,temp),                           &
     &                  Gfactor_T0+Gfactor_DT*REAL(Gfactor_Jm-1,r8))

            i1=INT(1.0_r8+(my_salt-Gfactor_S0)*oGfactor_DS)
            i2=MIN(i1+1,Gfactor_Im)
            j1=INT(1.0_r8+(my_temp-Gfactor_T0)*oGfactor_DT)
            j2=MIN(j1+1,Gfactor_Jm)

            p2=(my_salt-(Gfactor_S0+REAL(i1-1,r8)*Gfactor_DS))*         &
     &         oGfactor_DS
            q2=(my_temp-(Gfactor_T0+REAL(j1-1,r8)*Gfactor_DT))*         &
     &         oGfactor_DT
            p1=1.0_r8-p2
            q1=1.0_r8-q2

            Gfactor=p1*q1*Gfactor_table(i1,j1)+                         &
     &              p2*q1*Gfactor_table(i2,j1)+                         &
     &              p1*q2*Gfactor_table(i1,j2)+                         &
     &              p2*q2*Gfactor_table(i2,j2)
          END IF
!
!  Determine turbidity effect (linear or exponential) on larval growth.
!  Then, compute new larvae size (um) as function of growth rate (um/s)
!  which is loaded in track(ibrhs,:,:).
!
          IF (Lsize.gt.turb_size(ng)) THEN
            IF (Lturb.gt.turb_crit(ng)) THEN
              turb_ef=turb_base(ng)*                                    &
     &                EXP(-turb_rate(ng)*(Lturb-turb_mean(ng)))
            ELSE
              turb_ef=turb_slop(ng)*Lturb+turb_axis(ng)
            END IF
            track(ibrhs,nfp1,l)=Grate*Gfactor*turb_ef*sec2day
          ELSE
            track(ibrhs,nfp1,l)=Larvae_GR0(ng)*Gfactor*sec2day
          END IF
          IF (Predictor) THEN
            track(isizf,nfp1,l)=track(isizf,nfm3,l)+                    &
     &                          dt(ng)*(cff1*track(ibrhs,nf  ,l)-       &
     &                                  cff2*track(ibrhs,nfm1,l)+       &
     &                                  cff1*track(ibrhs,nfm2,l))
          ELSE
            track(isizf,nfp1,l)=cff1*track(isizf,nf  ,l)-               &
     &                          cff2*track(isizf,nfm2,l)+               &
     &                          dt(ng)*(cff3*track(ibrhs,nfp1,l)+       &
     &                                  cff4*track(ibrhs,nf  ,l)-       &
     &                                  cff3*track(ibrhs,nfm1,l))
          END IF
          Lsize=track(isizf,nfp1,l)
!
!  Estimate the fraction of time that the larvae spend swimming.
!
          IF (ABS(dsalt).lt.0.00001_r8) THEN
            dsalt=0.0_r8
          END IF
          IF (dsalt.gt.0.0_r8) THEN
            SwimTimeNew=MIN(SwimTime+dsalt*slope_Sinc(ng),swim_Tmax(ng))
          ELSE
            SwimTimeNew=MAX(SwimTime+dsalt*slope_Sdec(ng),swim_Tmin(ng))
          END IF
!
!  Compute swim behavior as function of larval size and temperature.
!  Linearly interpolate swimming rate (mm/s) from the look table of
!  larval size (um) versus temperature (Celsius).  Notice that
!  extrapolation is suppresed by bounding "Lsize" and "temp" to
!  the range values in the look table.
!
          IF ((temp.lt.swim_T0).or.(Lsize.lt.swim_L0)) THEN
            SwimRate=0.0_r8
          ELSE
            my_size=MIN(MAX(swim_L0,Lsize),                             &
     &                  swim_L0+swim_DL*REAL(swim_Im-1,r8))
            my_temp=MIN(MAX(swim_T0,temp),                              &
     &                  swim_T0+swim_DT*REAL(swim_Jm-1,r8))

            i1=INT(1.0_r8+(my_size-swim_L0)*oswim_DL)
            i2=MIN(i1+1,swim_Im)
            j1=INT(1.0_r8+(my_temp-swim_T0)*oswim_DT)
            j2=MIN(j1+1,swim_Jm)

            p2=(my_size-(swim_L0+REAL(i1-1,r8)*swim_DL))*oswim_DL
            q2=(my_temp-(swim_T0+REAL(j1-1,r8)*swim_DT))*oswim_DT
            p1=1.0_r8-p2
            q1=1.0_r8-q2

            SwimRate=p1*q1*swim_table(i1,j1)+                           &
     &               p2*q1*swim_table(i2,j1)+                           &
     &               p1*q2*swim_table(i1,j2)+                           &
     &               p2*q2*swim_table(i2,j2)

            SwimRate=SwimRate*0.001_r8        ! convert from mm/s to m/s
          END IF
!
!  Compute larvae sinking velocity (m/s).
!
          sink=sink_base(ng)*(EXP(sink_rate(ng)*(Lsize-sink_size(ng))))

          sink=sink*0.001_r8                  ! convert from mm/s to m/s
!
!  Compute larvae vertical velocity (m/s).
!
#ifdef GROWTH_ONLY
          w_bio=0.0_r8
#else
          w_bio=SwimTime*SwimRate-(1.0_r8-SwimTime)*sink
#endif
!
!  Load behavior into track array. Apply settlement condition: larvae
!  greater or equal than SETTLE_SIZE, settle on the bottom.
!
          IF (track(isizf,nfp1,l).lt.settle_size(ng)) THEN
            track(iwbio,nfp1,l)=w_bio
            track(iwsin,nfp1,l)=sink
            track(iswim,nfp1,l)=SwimTimeNew
          ELSE
            i1=MIN(MAX(0,INT(track(ixgrd,nfp1,l))),Lm(ng)+1)
            i2=MIN(i1+1,Lm(ng)+1)
            j1=MIN(MAX(0,INT(track(iygrd,nfp1,l))),Mm(ng)+1)
            j2=MIN(j1+1,Mm(ng)+1)

            p2=REAL(i2-i1,r8)*(track(ixgrd,nfp1,l)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(track(iygrd,nfp1,l)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2

            bottom=p1*q1*GRID(ng)%h(i1,j1)+                             &
     &             p2*q1*GRID(ng)%h(i2,j1)+                             &
     &             p1*q2*GRID(ng)%h(i1,j2)+                             &
     &             p2*q2*GRID(ng)%h(i2,j2)

            track(idpth,nfp1,l)=-bottom
            track(isizf,nfp1,l)=track(isizf,nf,l)
            track(iwbio,nfp1,l)=0.0_r8
            track(iwsin,nfp1,l)=0.0_r8
            track(iswim,nfp1,l)=0.0_r8
          END IF
!
!  If newly relased float, set vertical migration fields for all time
!  levels.
!
          IF (time(ng)-HalfDT.le.Tinfo(itstr,l).and.                    &
     &        time(ng)+HalfDT.gt.Tinfo(itstr,l)) THEN
            brhs=track(ibrhs,nfp1,l)
            sink=track(iwsin,nfp1,l)
            w_bio=track(iwbio,nfp1,l)
            DO i=0,NFT
              track(ibrhs,i,l)=brhs
              track(iwsin,i,l)=sink
              track(iwbio,i,l)=w_bio
            END DO
          END IF

        END IF
      END DO

      RETURN
      END SUBROUTINE biology_floats_tile
