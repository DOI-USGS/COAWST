      SUBROUTINE read_FltBioPar (model, inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in input biological floats parameters.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_behavior
      USE mod_iounits
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
      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status
      integer :: Ivalue(1)

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: Rvalue(1)

      real(r8), dimension(200) :: Rval

      character (len=35) :: frmt
      character (len=40) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval

      character (len=1 ), parameter :: blank = ' '
!
!-----------------------------------------------------------------------
!  Read in initial float locations.
!-----------------------------------------------------------------------
!
!  Allocate oyster model parameter.
!
      CALL allocate_behavior
!
!  Notice I added one when allocating local scratch arrays to avoid
!  out of bounds in some compilers when reading the last blank line
!  which signal termination of input data.
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Larvae_size0')
              Npts=load_r(Nval, Rval, Ngrids, Larvae_size0)
            CASE ('Larvae_GR0')
              Npts=load_r(Nval, Rval, Ngrids, Larvae_GR0)
            CASE ('settle_size')
              Npts=load_r(Nval, Rval, Ngrids, settle_size)
            CASE ('food_supply')
              Npts=load_r(Nval, Rval, Ngrids, food_supply)
            CASE ('turb_ambi')
              Npts=load_r(Nval, Rval, Ngrids, turb_ambi)
            CASE ('turb_crit')
              Npts=load_r(Nval, Rval, Ngrids, turb_crit)
            CASE ('turb_slop')
              Npts=load_r(Nval, Rval, Ngrids, turb_slop)
            CASE ('turb_axis')
              Npts=load_r(Nval, Rval, Ngrids, turb_axis)
            CASE ('turb_base')
              Npts=load_r(Nval, Rval, Ngrids, turb_base)
            CASE ('turb_rate')
              Npts=load_r(Nval, Rval, Ngrids, turb_rate)
            CASE ('turb_mean')
              Npts=load_r(Nval, Rval, Ngrids, turb_mean)
            CASE ('turb_size')
              Npts=load_r(Nval, Rval, Ngrids, turb_size)
            CASE ('swim_Tmin')
              Npts=load_r(Nval, Rval, Ngrids, swim_Tmin)
            CASE ('swim_Tmax')
              Npts=load_r(Nval, Rval, Ngrids, swim_Tmax)
            CASE ('swim_Sinc')
              Npts=load_r(Nval, Rval, Ngrids, swim_Sinc)
            CASE ('swim_Sdec')
              Npts=load_r(Nval, Rval, Ngrids, swim_Sdec)
            CASE ('slope_Sinc')
              Npts=load_r(Nval, Rval, Ngrids, slope_Sinc)
            CASE ('slope_Sdec')
              Npts=load_r(Nval, Rval, Ngrids, slope_Sdec)
            CASE ('sink_base')
              Npts=load_r(Nval, Rval, Ngrids, sink_base)
            CASE ('sink_rate')
              Npts=load_r(Nval, Rval, Ngrids, sink_rate)
            CASE ('sink_size')
              Npts=load_r(Nval, Rval, Ngrids, sink_size)
            CASE ('swim_Im')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              swim_Im=Ivalue(1)
            CASE ('swim_Jm')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              swim_Jm=Ivalue(1)
            CASE ('swim_L0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              swim_L0=Rvalue(1)
            CASE ('swim_T0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              swim_T0=Rvalue(1)
            CASE ('swim_DL')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              swim_DL=Rvalue(1)
            CASE ('swim_DT')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              swim_DT=Rvalue(1)
            CASE ('swim_table')
              IF (.not.allocated(swim_table)) THEN
                allocate ( swim_table(swim_Im,swim_Jm) )
                swim_table=0.0_r8
              END IF
              READ (inp,*,ERR=20,END=30)                                &
                   ((swim_table(i,j),i=1,swim_Im),j=1,swim_Jm)
            CASE ('Gfactor_Im')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Gfactor_Im=Ivalue(1)
            CASE ('Gfactor_Jm')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Gfactor_Jm=Ivalue(1)
            CASE ('Gfactor_S0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Gfactor_S0=Rvalue(1)
            CASE ('Gfactor_T0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Gfactor_T0=Rvalue(1)
            CASE ('Gfactor_DS')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Gfactor_DS=Rvalue(1)
            CASE ('Gfactor_DT')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Gfactor_DT=Rvalue(1)
            CASE ('Gfactor_table')
              IF (.not.allocated(Gfactor_table)) THEN
                allocate ( Gfactor_table(Gfactor_Im,Gfactor_Jm) )
                Gfactor_table=0.0_r8
              END IF
              READ (inp,*,ERR=20,END=30)                                &
                   ((Gfactor_table(i,j),i=1,Gfactor_Im),j=1,Gfactor_Jm)
            CASE ('Grate_Im')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Grate_Im=Ivalue(1)
            CASE ('Grate_Jm')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Grate_Jm=Ivalue(1)
            CASE ('Grate_F0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Grate_F0=Rvalue(1)
            CASE ('Grate_L0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Grate_L0=Rvalue(1)
            CASE ('Grate_DF')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Grate_DF=Rvalue(1)
            CASE ('Grate_DL')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Grate_DL=Rvalue(1)
            CASE ('Grate_table')
              IF (.not.allocated(Grate_table)) THEN
                allocate ( Grate_table(Grate_Im,Grate_Jm) )
                Grate_table=0.0_r8
              END IF
              READ (inp,*,ERR=20,END=30)                                &
                   ((Grate_table(i,j),i=1,Grate_Im),j=1,Grate_Jm)
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  20  IF (Master) WRITE (out,50) TRIM(KeyWord)
      exit_flag=4
      RETURN
  30  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          WRITE (out,60) ng
          WRITE (out,70) Larvae_size0(ng), 'Larvae_size0',              &
     &          'Initial larvae size (um).'
          WRITE (out,70) Larvae_GR0(ng), 'Larvae_GR0',                  &
     &          'Initial larvae growth rate (um/day).'
          WRITE (out,70) settle_size(ng), 'settle_size',                &
     &          'Larvae settlement size (um).'
          WRITE (out,70) food_supply(ng), 'food_supply',                &
     &          'Food supply (mg Carbon/l).'
          WRITE (out,70) turb_ambi(ng), 'turb_ambi',                    &
     &          'Ambient turbidity level (g/l).'
          WRITE (out,70) turb_crit(ng), 'turb_crit',                    &
     &          'Critical turbidity value (g/l).'
          WRITE (out,70) turb_slop(ng), 'turb_slop',                    &
     &          'Turbidity linear slope (l/g).'
          WRITE (out,70) turb_axis(ng), 'turb_axis',                    &
     &          'Turbidity linear axis crossing (g/l).'
          WRITE (out,70) turb_base(ng), 'turb_base',                    &
     &          'Turbidity exponential base factor (g/l).'
          WRITE (out,70) turb_rate(ng), 'turb_rate',                    &
     &          'Turbidity exponential rate (l/g).'
          WRITE (out,70) turb_mean(ng), 'turb_mean',                    &
     &          'Turbidity exponential mean (g/l).'
          WRITE (out,70) turb_size(ng), 'turb_size',                    &
     &          'Minimum larvae size (um) affected by turbidity.'
          WRITE (out,70) swim_Tmin(ng), 'swim_Tmin',                    &
     &          'Minimum swimming time fraction.'
          WRITE (out,70) swim_Tmax(ng), 'swim_Tmax',                    &
     &          'Maximum swimming time fraction.'
          WRITE (out,70) swim_Sinc(ng), 'swim_Sinc',                    &
     &          'Swimming, active fraction due to increasing salinity.'
          WRITE (out,70) swim_Sdec(ng), 'swim_Sdec',                    &
     &          'Swimming, active fraction due to decreasing salinity.'
          WRITE (out,70) slope_Sinc(ng), 'slope_Sinc',                  &
     &          'Swimming, coefficient due to increasing salinity.'
          WRITE (out,70) slope_Sdec(ng), 'slope_Sdec',                  &
     &          'Swimming, coefficient due to increasing salinity.'
          WRITE (out,70) sink_base(ng), 'sink_base',                    &
     &          'Sinking, exponential base factor (mm/s).'
          WRITE (out,70) sink_rate(ng), 'sink_rate',                    &
     &          'Sinking, exponential rate factor (1/um).'
          WRITE (out,70) sink_size(ng), 'sink_mean',                    &
     &          'Sinking, exponential mean size (um).'
          WRITE (out,80) swim_Im, 'swim_Im',                            &
     &          'Swim table, number of values in larval size I-axis.'
          WRITE (out,80) swim_Jm, 'swim_Jm',                            &
     &          'Swim table, number of values in temperature J-axis.'
          WRITE (out,70) swim_L0, 'swim_L0',                            &
     &          'Swim table, starting value for larval size I-axis.'
          WRITE (out,70) swim_T0, 'swim_T0',                            &
     &          'Swim table, starting value for temperature J-axis.'
          WRITE (out,70) swim_DL, 'swim_DL',                            &
     &          'Swim table, larval size I-axis increment.'
          WRITE (out,70) swim_DT, 'swim_DT',                            &
     &          'Swim table, temperature J-axis increment.'
          WRITE (out,80) Gfactor_Im, 'Gfactor_Im',                      &
     &          'Gfactor table, number of values in salinity I-axis.'
          WRITE (out,80) Gfactor_Jm, 'Gfactor_Jm',                      &
     &          'Gfactor table, number of values in temperature J-axis.'
          WRITE (out,70) Gfactor_S0, 'Gfactor_S0',                      &
     &          'Gfactor table, starting value for salinity I-axis.'
          WRITE (out,70) Gfactor_T0, 'Gfactor_T0',                      &
     &          'Gfactor table, starting value for temperature J-axis.'
          WRITE (out,70) Gfactor_DS, 'Gfactor_DS',                      &
     &          'Gfactor table, starting value for salinity I-axis.'
          WRITE (out,70) Gfactor_DT, 'Gfactor_DT',                      &
     &          'Gfactor table, starting value for temperature J-axis.'
          WRITE (out,80) Grate_Im, 'Grate_Im',                          &
     &          'Grate table, number of values in food supply I-axis.'
          WRITE (out,80) Grate_Jm, 'Grate_Jm',                          &
     &          'Grate table, number of values in larval size J-axis.'
          WRITE (out,70) Grate_F0, 'Grate_F0',                          &
     &          'Grate table, starting value for food supply I-axis.'
          WRITE (out,70) Grate_L0, 'Grate_L0',                          &
     &          'Grate table, starting value for larval size J-axis.'
          WRITE (out,70) Grate_DF, 'Grate_DF',                          &
     &          'Grate table, food supply I-axis increment.'
          WRITE (out,70) Grate_DL, 'Grate_DL',                          &
     &          'Grate table, larval size J-axis increment.'
        END DO
      END IF

  40  FORMAT (/,' READ_FloatsBioPar - Error while processing line: ',/, &
     &        a)
  50  FORMAT (/,' READ_FloatsBioPar - Error reading look table: ',a)
  60  FORMAT (/,/,' Biological Floats Behavior Parameters, Grid: ',i2.2, &
     &        /,  ' ===============================================',/)
  70  FORMAT (1p,e11.4,2x,a,t32,a)
  80  FORMAT (1x,i10,2x,a,t32,a)

      RETURN
      END SUBROUTINE read_FltBioPar
