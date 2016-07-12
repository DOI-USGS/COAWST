#include "cppdefs.h"
#define LONG_NUMS
#ifdef NONLINEAR
      SUBROUTINE output (ng)
!
!svn $Id: output.F 1597 2008-07-25 02:31:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine manages nonlinear model output. It creates output   !
!  NetCDF files and writes out data into NetCDF files. If requested,   !
!  it can create several history and/or time-averaged files to avoid   !
!  generating too large files during a single model run.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FILTERED
      USE mod_filter, ONLY: FILN
# endif
# ifdef FLOATS
      USE mod_floats
# endif
# if defined FOUR_DVAR || defined VERIFICATION
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasts
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, NewFile

      integer :: Fcount, ifile, status, tile
# ifdef FILTERED_RST
      integer :: i
# endif
!
      SourceFile='output.F'

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 8)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process nonlinear history NetCDF file.
!-----------------------------------------------------------------------
!
!  Set tile for local array manipulations in output routines.
!
#ifdef DISTRIBUTE
      tile=MyRank
#else
      tile=-1
#endif
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  Create output history NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
# ifndef NO_HIS
      IF (LdefHIS(ng)) THEN
        IF (ndefHIS(ng).gt.0) THEN
          IF (idefHIS(ng).lt.0) THEN
            idefHIS(ng)=((ntstart(ng)-1)/ndefHIS(ng))*ndefHIS(ng)
            IF (idefHIS(ng).lt.iic(ng)-1) THEN
              idefHIS(ng)=idefHIS(ng)+ndefHIS(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefHIS(ng)) THEN
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              Ldefine=.TRUE.
              NewFile=.FALSE.                 ! unfinished file, inquire
            END IF                            ! content for appending
            idefHIS(ng)=idefHIS(ng)+nHIS(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefHIS(ng)) THEN
            idefHIS(ng)=idefHIS(ng)+ndefHIS(ng)
            IF (nHIS(ng).ne.ndefHIS(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefHIS(ng)=idefHIS(ng)+nHIS(ng)  ! multiple record offset
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN                     ! create new file or
            Fcount=HIS(ng)%Fcount               ! inquire existing file
            HIS(ng)%Nrec(Fcount)=0
            ifile=(iic(ng)-1)/ndefHIS(ng)+1
            IF (Master) THEN
              WRITE (HIS(ng)%name,10) TRIM(HIS(ng)%base), ifile
#ifdef LONG_NUMS
  10          FORMAT (a,'_',i5.5,'.nc')
#else
  10          FORMAT (a,'_',i4.4,'.nc')
#endif
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, HIS(ng)%name)
# endif
            IF (HIS(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, HIS(ng)%ncid)
            END IF
            CALL def_his (ng, NewFile)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtHIS(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtHIS(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_his (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtHIS(ng)=.TRUE.
            LdefHIS(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into history NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtHIS(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nHIS(ng)).eq.0)) THEN
            IF (nrrec(ng).eq.0.or.iic(ng).ne.ntstart(ng)) THEN
              CALL wrt_his (ng, tile)
            END IF
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nHIS(ng)).eq.0) THEN
            CALL wrt_his (ng, tile)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
      END IF
# endif
# if defined HISTORY2
!
!-----------------------------------------------------------------------
!  If appropriate, process secondary surface NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged NetCDF file or prepare existing file
!  to append new data to it. Also, notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefHIS2(ng)) THEN
        IF (ndefHIS2(ng).gt.0) THEN
          IF (idefHIS2(ng).lt.0) THEN
            idefHIS2(ng)=((ntstart(ng)-1)/ndefHIS2(ng))*ndefHIS2(ng)
            IF ((ndefHIS2(ng).eq.nHIS2(ng)).and.(idefHIS2(ng).le.0)) THEN
              idefHIS2(ng)=ndefHIS2(ng)         ! one file per record
            ELSE IF (idefHIS2(ng).lt.iic(ng)-1) THEN
              idefHIS2(ng)=idefHIS2(ng)+ndefHIS2(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefHIS2(ng)) THEN
              Ldefine=.FALSE.                    ! finished file, delay
            ELSE                                 ! creation of next file
              NewFile=.FALSE.
              Ldefine=.TRUE.                     ! unfinished file, inquire
            END IF                               ! content for appending
            idefHIS2(ng)=idefHIS2(ng)+nHIS2(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefHIS2(ng)) THEN
            idefHIS2(ng)=idefHIS2(ng)+ndefHIS2(ng)
            IF (nHIS2(ng).ne.ndefHIS2(ng).and.iic(ng).eq.ntstart(ng))   &
     &              THEN
              idefHIS2(ng)=idefHIS2(ng)+nHIS2(ng)
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN
            Fcount=HIS2(ng)%Fcount
            HIS2(ng)%Nrec(Fcount)=0
            IF (nHIS2(ng).eq.ndefHIS2(ng)) THEN
              ifile=(iic(ng)-1)/ndefHIS2(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefHIS2(ng)+1
            END IF
            IF (Master) THEN
              WRITE (HIS2(ng)%name,12) TRIM(HIS2(ng)%base), ifile
            END IF
#ifdef LONG_NUMS
  12          FORMAT (a,'_',i5.5,'.nc')
#else
  12          FORMAT (a,'_',i4.4,'.nc')
#endif
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, HIS2(ng)%name)
# endif
            IF (HIS2(ng)%ncid.ne.-1) THEN
              CALL netcdf_close(ng, iNLM, HIS2(ng)%ncid)
            END IF
            CALL def_his2 (ng, Newfile)
            IF (exit_flag.ne.NoError) RETURN
            LwrtHIS2(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_his2 (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtHIS2(ng)=.TRUE.
            LdefHIS2(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF (LwrtHIS2(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nHIS2(ng)).eq.0)) THEN
          CALL wrt_his2 (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifdef AVERAGES
!
!-----------------------------------------------------------------------
!  If appropriate, process time-averaged NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged NetCDF file or prepare existing file
!  to append new data to it. Also, notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefAVG(ng)) THEN
        IF (ndefAVG(ng).gt.0) THEN
          IF (idefAVG(ng).lt.0) THEN
            idefAVG(ng)=((ntstart(ng)-1)/ndefAVG(ng))*ndefAVG(ng)
            IF ((ndefAVG(ng).eq.nAVG(ng)).and.(idefAVG(ng).le.0)) THEN
              idefAVG(ng)=ndefAVG(ng)         ! one file per record
            ELSE IF (idefAVG(ng).lt.iic(ng)-1) THEN
              idefAVG(ng)=idefAVG(ng)+ndefAVG(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefAVG(ng)) THEN
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              NewFile=.FALSE.
              Ldefine=.TRUE.                  ! unfinished file, inquire
            END IF                            ! content for appending
            idefAVG(ng)=idefAVG(ng)+nAVG(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefAVG(ng)) THEN
            idefAVG(ng)=idefAVG(ng)+ndefAVG(ng)
            IF (nAVG(ng).ne.ndefAVG(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefAVG(ng)=idefAVG(ng)+nAVG(ng)
            END IF
            Ldefine=.TRUE.
            Newfile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN
            Fcount=AVG(ng)%Fcount
            AVG(ng)%Nrec(Fcount)=0
            IF (ndefAVG(ng).eq.nAVG(ng)) THEN
              ifile=(iic(ng)-1)/ndefAVG(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefAVG(ng)+1
            END IF
            IF (Master) THEN
              WRITE (AVG(ng)%name,20) TRIM(AVG(ng)%base), ifile
#ifdef LONG_NUMS
  20          FORMAT (a,'_',i5.5,'.nc')
#else
  20          FORMAT (a,'_',i4.4,'.nc')
#endif
            END IF
#  ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, AVG(ng)%name)
#  endif
            IF (AVG(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, AVG(ng)%ncid)
            END IF
            IF (ifile .ne. 0) CALL def_avg (ng, Newfile)
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_avg (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG(ng)=.TRUE.
            LdefAVG(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF (LwrtAVG(ng)) THEN
        IF (((iic(ng).gt.ntstart(ng)).and.                              &
     &       (MOD(iic(ng)-1,nAVG(ng)).eq.0)).or.                        &
     &      ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1))) THEN
          CALL wrt_avg (ng)
          IF (exit_flag.ne.NoError) RETURN
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          CALL wrt_tides (ng)
          IF (exit_flag.ne.NoError) RETURN
#  endif
        END IF
      END IF
# endif
# if defined AVERAGES2
!
!-----------------------------------------------------------------------
!  If appropriate, process secondary time-averaged NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged NetCDF file or prepare existing file
!  to append new data to it. Also, notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefAVG2(ng)) THEN
        IF (ndefAVG2(ng).gt.0) THEN
          IF (idefAVG2(ng).lt.0) THEN
            idefAVG2(ng)=((ntstart(ng)-1)/ndefAVG2(ng))*ndefAVG2(ng)
            IF ((ndefAVG2(ng).eq.nAVG2(ng)).and.(idefAVG2(ng).le.0)) THEN
              idefAVG2(ng)=ndefAVG2(ng)         ! one file per record
            ELSE IF (idefAVG2(ng).lt.iic(ng)-1) THEN
              idefAVG2(ng)=idefAVG2(ng)+ndefAVG2(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefAVG2(ng)) THEN
              Ldefine=.FALSE.                    ! finished file, delay
            ELSE                                 ! creation of next file
              NewFile=.FALSE.
              Ldefine=.TRUE.                     ! unfinished file, inquire
            END IF                               ! content for appending
            idefAVG2(ng)=idefAVG2(ng)+nAVG2(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefAVG2(ng)) THEN
            idefAVG2(ng)=idefAVG2(ng)+ndefAVG2(ng)
            IF (nAVG2(ng).ne.ndefAVG2(ng).and.iic(ng).eq.ntstart(ng))   &
     &              THEN
              idefAVG2(ng)=idefAVG2(ng)+nAVG2(ng)
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN
            Fcount=AVG2(ng)%Fcount
            AVG2(ng)%Nrec(Fcount)=0
            IF (nAVG2(ng).eq.ndefAVG2(ng)) THEN
              ifile=(iic(ng)-1)/ndefAVG2(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefAVG2(ng)+1
            END IF
            IF (Master) THEN
              WRITE (AVG2(ng)%name,22) TRIM(AVG2(ng)%base), ifile
            END IF
#ifdef LONG_NUMS
  22          FORMAT (a,'_',i5.5,'.nc')
#else
  22          FORMAT (a,'_',i4.4,'.nc')
#endif
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, AVG2(ng)%name)
# endif
            IF (AVG2(ng)%ncid.ne.-1) THEN
              CALL netcdf_close(ng, iNLM, AVG2(ng)%ncid)
            END IF
            CALL def_avg2 (ng, Newfile)
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG2(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_avg2 (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG2(ng)=.TRUE.
            LdefAVG2(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF (LwrtAVG2(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nAVG2(ng)).eq.0)) THEN
          CALL wrt_avg2 (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  If appropriate, process time-averaged diagnostics NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged diagnostics NetCDF file or prepare
!  existing file to append new data to it. Also, notice that it is
!  possible to create several files during a single model run.
!
      IF (LdefDIA(ng)) THEN
        IF (ndefDIA(ng).gt.0) THEN
          IF (idefDIA(ng).lt.0) THEN
            idefDIA(ng)=((ntstart(ng)-1)/ndefDIA(ng))*ndefDIA(ng)
            IF ((ndefDIA(ng).eq.nDIA(ng)).and.(idefDIA(ng).le.0)) THEN
              idefDIA(ng)=ndefDIA(ng)         ! one file per record
            ELSE IF (idefDIA(ng).lt.iic(ng)-1) THEN
              idefDIA(ng)=idefDIA(ng)+ndefDIA(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefDIA(ng)) THEN
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              NewFile=.FALSE.
              Ldefine=.TRUE.                  ! unfinished file, inquire
            END IF                            ! content for appending
            idefDIA(ng)=idefDIA(ng)+nDIA(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefDIA(ng)) THEN
            idefDIA(ng)=idefDIA(ng)+ndefDIA(ng)
            IF (nDIA(ng).ne.ndefDIA(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefDIA(ng)=idefDIA(ng)+nDIA(ng)
            END IF
            Ldefine=.TRUE.
            Newfile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN
            Fcount=DIA(ng)%Fcount
            DIA(ng)%Nrec(Fcount)=0
            IF (ndefDIA(ng).eq.nDIA(ng)) THEN
              ifile=(iic(ng)-1)/ndefDIA(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefDIA(ng)+1
            END IF
            IF (Master) THEN
              WRITE (DIA(ng)%name,30) TRIM(DIA(ng)%base), ifile
#ifdef LONG_NUMS
  30          FORMAT (a,'_',i5.5,'.nc')
#else
  30          FORMAT (a,'_',i4.4,'.nc')
#endif
            END IF
#  ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, DIA(ng)%name)
#  endif
            IF (DIA(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, DIA(ng)%ncid)
            END IF
            CALL def_diags (ng, Newfile)
            IF (exit_flag.ne.NoError) RETURN
            LwrtDIA(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_diags (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtDIA(ng)=.TRUE.
            LdefDIA(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged diagnostics NetCDF file.
!
      IF (LwrtDIA(ng)) THEN
        IF (((iic(ng).gt.ntstart(ng)).and.                              &
     &       (MOD(iic(ng)-1,nDIA(ng)).eq.0)).or.                        &
     &      ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1))) THEN
          CALL wrt_diags (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifdef STATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process stations NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Lstations(ng).and.                                            &
     &    (Nstation(ng).gt.0).and.(nSTA(ng).gt.0)) THEN
!
!  Create output station NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefSTA(ng).and.(iic(ng).eq.ntstart(ng))) THEN
#if defined NEP5 || defined NEP6 || defined CASPIAN || defined CHUKCHI \
           || defined CORAL || defined NWA || defined ARC_NATL \
           || defined ARCTIC || defined NWGOA
          CALL def_station (ng, .true.)
#else
          CALL def_station (ng, ldefout(ng))
#endif
          IF (exit_flag.ne.NoError) RETURN
          LdefSTA(ng)=.FALSE.
        END IF
!
!  Write out data into stations NetCDF file.
!
        IF (MOD(iic(ng)-1,nSTA(ng)).eq.0) THEN
          CALL wrt_station (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  If appropriate, process floats NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Lfloats(ng).and.                                              &
     &    (Nfloats(ng).gt.0).and.(nFLT(ng).gt.0)) THEN
!
!  Create output floats NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefFLT(ng)) THEN
          IF (frrec(ng).eq.0) THEN
            NewFile=.TRUE.
          ELSE
            NewFile=.FALSE.
          END IF
          CALL def_floats (ng, NewFile)
          IF (exit_flag.ne.NoError) RETURN
          LdefFLT(ng)=.FALSE.
        END IF
!
!  Write out data into floats NetCDF file.
!
        IF ((MOD(iic(ng)-1,nFLT(ng)).eq.0).and.                         &
     &      ((frrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
          CALL wrt_floats (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process restart NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output restart NetCDF file or prepare existing file to
!  append new data to it.
!
      IF (LdefRST(ng)) THEN
        CALL def_rst (ng)
        IF (exit_flag.ne.NoError) RETURN
        LwrtRST(ng)=.TRUE.
        LdefRST(ng)=.FALSE.
#ifdef FILTERED_RST
        CALL def_filt (ng)
#endif
      END IF
!
!  Write out data into restart NetCDF file.
!
      IF (LwrtRST(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     & ((MOD(iic(ng)-1,nRST(ng)).eq.0) .or. ((iic(ng)-1)==ntimes(ng)))) THEN
          CALL wrt_rst (ng)
          IF (exit_flag.ne.NoError) RETURN
# ifdef FILTERED_RST
          DO i=1,FILN
            CALL wrt_filt (ng, i)
            IF (exit_flag.ne.NoError) RETURN
          END DO
# endif
        END IF
      END IF

# if (defined FOUR_DVAR || defined VERIFICATION) && \
     !defined IS4DVAR_SENSITIVITY
#  ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process and write model state at observation
!  locations. Compute misfit (model-observations) cost function.
!-----------------------------------------------------------------------
!
      IF (((time(ng)-0.5_r8*dt(ng)).le.ObsTime(ng)).and.                &
     &    (ObsTime(ng).lt.(time(ng)+0.5_r8*dt(ng)))) THEN
        ProcessObs=.TRUE.
        CALL obs_read (ng, iNLM, .FALSE.)
        CALL obs_write (ng, tile, iNLM)
#   if !(defined VERIFICATION || defined W4DVAR)
        CALL obs_cost (ng, iNLM)
#   endif
      ELSE
        ProcessObs=.FALSE.
      END IF
#  endif
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 8)
# endif
      RETURN
      END SUBROUTINE output
#else
      SUBROUTINE output
      RETURN
      END SUBROUTINE output
#endif
