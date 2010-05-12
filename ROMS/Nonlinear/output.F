#include "cppdefs.h"
#ifdef NONLINEAR
      SUBROUTINE output (ng)
!
!svn $Id: output.F 1597 2008-07-25 02:31:35Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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

      integer :: ifile, lstr, status, tile
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
            NrecHIS(ng)=0                       ! inquire existing file
            ifile=(iic(ng)-1)/ndefHIS(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(HISbase(ng))
              WRITE (HISname(ng),10) HISbase(ng)(1:lstr-3),ifile
  10          FORMAT (a,'_',i4.4,'.nc')
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, HISname(ng))
# endif
            IF (ncHISid(ng).ne.-1) THEN
              CALL netcdf_close (ng, iNLM, ncHISid(ng))
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
              CALL wrt_his (ng)
            END IF
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nHIS(ng)).eq.0) THEN
            CALL wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
      END IF

# if defined AVERAGES && !defined ADJOINT
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
            NrecAVG(ng)=0
            IF (ndefAVG(ng).eq.nAVG(ng)) THEN
              ifile=(iic(ng)-1)/ndefAVG(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefAVG(ng)+1
            END IF
            IF (Master) THEN
              lstr=LEN_TRIM(AVGbase(ng))
              WRITE (AVGname(ng),20) AVGbase(ng)(1:lstr-3),ifile
  20          FORMAT (a,'_',i4.4,'.nc')
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, AVGname(ng))
# endif
            IF (ncAVGid(ng).ne.-1) THEN
              CALL netcdf_close (ng, iNLM, ncAVGid(ng))
            END IF
            CALL def_avg (ng, Newfile)
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
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nAVG(ng)).eq.0)) THEN
          CALL wrt_avg (ng)
          IF (exit_flag.ne.NoError) RETURN
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          CALL wrt_tides (ng)
          IF (exit_flag.ne.NoError) RETURN
# endif
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
            NrecDIA(ng)=0
            IF (ndefDIA(ng).eq.nDIA(ng)) THEN
              ifile=(iic(ng)-1)/ndefDIA(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefDIA(ng)+1
            END IF
            IF (Master) THEN
              lstr=LEN_TRIM(DIAbase(ng))
              WRITE (DIAname(ng),30) DIAbase(ng)(1:lstr-3),ifile
  30          FORMAT (a,'_',i4.4,'.nc')
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, DIAname(ng))
# endif
            IF (ncDIAid(ng).ne.-1) THEN
              CALL netcdf_close (ng, iNLM, ncDIAid(ng))
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
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nDIA(ng)).eq.0)) THEN
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
          CALL def_station (ng, ldefout(ng))
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
      END IF
!
!  Write out data into restart NetCDF file.
!
      IF (LwrtRST(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nRST(ng)).eq.0)) THEN
          CALL wrt_rst (ng)
          IF (exit_flag.ne.NoError) RETURN
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
#   ifdef DISTRIBUTE
        tile=MyRank
#   else
        tile=-1
#   endif
        CALL obs_read (ng, iNLM, .FALSE.)
        CALL obs_write (ng, tile, iNLM)
#   ifndef W4DVAR
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
