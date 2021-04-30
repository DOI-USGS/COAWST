#include "cppdefs.h"
#define LONG_NUMS
#ifdef NONLINEAR
      SUBROUTINE output (ng)
!
!svn $Id: output.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasts
# endif
      USE strings_mod,    ONLY : FoundError
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
!
      integer :: Fcount, ifile, status, tile
# ifdef FILTERED_RST
      integer :: i
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 8, __LINE__, MyFile)
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
              HIS(ng)%load=0                  ! restart, reset counter
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
            IF (iic(ng).eq.ntstart(ng)) THEN    ! inquire existing file
              HIS(ng)%load=0                    ! reset filename counter
            END IF
            ifile=(iic(ng)-1)/ndefHIS(ng)+1     ! next filename suffix
            HIS(ng)%load=HIS(ng)%load+1
            IF (HIS(ng)%load.gt.HIS(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'HIS(ng)%load = ', HIS(ng)%load,      &
     &                             HIS(ng)%Nfiles, TRIM(HIS(ng)%base),  &
     &                             ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            Fcount=HIS(ng)%load
            HIS(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (HIS(ng)%name,20) TRIM(HIS(ng)%base), ifile
            END IF
#  ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, HIS(ng)%name)
#  endif
            HIS(ng)%files(Fcount)=TRIM(HIS(ng)%name)
            IF (HIS(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, HIS(ng)%ncid)
            END IF
            CALL def_his (ng, NewFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtHIS(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtHIS(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_his (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nHIS(ng)).eq.0) THEN
            CALL wrt_his (ng, tile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process nonlinear quicksave NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output quicksave NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefQCK(ng)) THEN
        IF (ndefQCK(ng).gt.0) THEN
          IF (idefQCK(ng).lt.0) THEN
            idefQCK(ng)=((ntstart(ng)-1)/ndefQCK(ng))*ndefQCK(ng)
            IF (idefQCK(ng).lt.iic(ng)-1) THEN
              idefQCK(ng)=idefQCK(ng)+ndefQCK(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefQCK(ng)) THEN
              QCK(ng)%load=0                  ! restart, reset counter
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              Ldefine=.TRUE.
              NewFile=.FALSE.                 ! unfinished file, inquire
            END IF                            ! content for appending
            idefQCK(ng)=idefQCK(ng)+nQCK(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefQCK(ng)) THEN
            idefQCK(ng)=idefQCK(ng)+ndefQCK(ng)
            IF (nQCK(ng).ne.ndefQCK(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefQCK(ng)=idefQCK(ng)+nQCK(ng)  ! multiple record offset
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN                     ! create new file or
            IF (iic(ng).eq.ntstart(ng)) THEN    ! inquire existing file
              QCK(ng)%load=0                    ! reset filename counter
            END IF
            ifile=(iic(ng)-1)/ndefQCK(ng)+1     ! next filename suffix
            QCK(ng)%load=QCK(ng)%load+1
            IF (QCK(ng)%load.gt.QCK(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'QCK(ng)%load = ', QCK(ng)%load,      &
     &                             QCK(ng)%Nfiles, TRIM(QCK(ng)%base),  &
     &                             ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            Fcount=QCK(ng)%load
            QCK(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (QCK(ng)%name,20) TRIM(QCK(ng)%base), ifile
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, QCK(ng)%name)
# endif
            QCK(ng)%files(Fcount)=TRIM(QCK(ng)%name)
            IF (QCK(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, QCK(ng)%ncid)
            END IF
            CALL def_quick (ng, NewFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtQCK(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtQCK(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_quick (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtQCK(ng)=.TRUE.
            LdefQCK(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into quicksave NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtQCK(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nQCK(ng)).eq.0)) THEN
            IF (nrrec(ng).eq.0.or.iic(ng).ne.ntstart(ng)) THEN
              CALL wrt_quick (ng, tile)
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nQCK(ng)).eq.0) THEN
            CALL wrt_quick (ng, tile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END IF

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
              AVG(ng)%load=0                  ! restart, reset counter
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
            IF (iic(ng).eq.ntstart(ng)) THEN
              AVG(ng)%load=0                  ! reset filename counter
            END IF
            IF (ndefAVG(ng).eq.nAVG(ng)) THEN ! next filename suffix
              ifile=(iic(ng)-1)/ndefAVG(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefAVG(ng)+1
            END IF
            AVG(ng)%load=AVG(ng)%load+1
            IF (AVG(ng)%load.gt.AVG(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'AVG(ng)%load = ', AVG(ng)%load,      &
     &                             AVG(ng)%Nfiles, TRIM(AVG(ng)%base),  &
     &                             ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            Fcount=AVG(ng)%load
            AVG(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (AVG(ng)%name,20) TRIM(AVG(ng)%base), ifile
            END IF
#  ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, AVG(ng)%name)
#  endif
            AVG(ng)%files(Fcount)=TRIM(AVG(ng)%name)
            IF (AVG(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, AVG(ng)%ncid)
            END IF
            CALL def_avg (ng, Newfile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtAVG(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_avg (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          CALL wrt_tides (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
              DIA(ng)%load=0                  ! restart, reset counter
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
            IF (iic(ng).eq.ntstart(ng)) THEN
              DIA(ng)%load=0                  ! reset filename counter
            END IF
            IF (ndefDIA(ng).eq.nDIA(ng)) THEN ! next filename suffix
              ifile=(iic(ng)-1)/ndefDIA(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefDIA(ng)+1
            END IF
            DIA(ng)%load=DIA(ng)%load+1
            IF (DIA(ng)%load.gt.DIA(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'DIA(ng)%load = ', DIA(ng)%load,      &
     &                             DIA(ng)%Nfiles, TRIM(DIA(ng)%base),  &
     &                             ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            Fcount=DIA(ng)%load
            DIA(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (DIA(ng)%name,20) TRIM(DIA(ng)%base), ifile
            END IF
#  ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, DIA(ng)%name)
#  endif
            DIA(ng)%files(Fcount)=TRIM(DIA(ng)%name)
            IF (DIA(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iNLM, DIA(ng)%ncid)
            END IF
            CALL def_diags (ng, Newfile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtDIA(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_diags (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          LdefSTA(ng)=.FALSE.
        END IF
!
!  Write out data into stations NetCDF file.
!
        IF (MOD(iic(ng)-1,nSTA(ng)).eq.0) THEN
          CALL wrt_station (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          LdefFLT(ng)=.FALSE.
        END IF
!
!  Write out data into floats NetCDF file.
!
        IF ((MOD(iic(ng)-1,nFLT(ng)).eq.0).and.                         &
     &      ((frrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
          CALL wrt_floats (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
     &      (MOD(iic(ng)-1,nRST(ng)).eq.0)) THEN
          CALL wrt_rst (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# ifdef FILTERED_RST
          DO i=1,FILN
            CALL wrt_filt (ng, i)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END DO
# endif
        END IF
      END IF

# if (defined FOUR_DVAR              || \
      defined VERIFICATION)          && \
     !defined I4DVAR_ANA_SENSITIVITY
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
#   ifdef SP4DVAR
        IF (iic(ng).ne.ntend(ng)+1.or.Lsadd(ng)) THEN
          CALL obs_write (ng, tile, iNLM)
          CALL obs_cost (ng, iNLM)
        END IF
#   else
        CALL obs_write (ng, tile, iNLM)
#    if !(defined R4DVAR || defined VERIFICATION)
        CALL obs_cost (ng, iNLM)
#    endif
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
      CALL wclock_off (ng, iNLM, 8, __LINE__, MyFile)
# endif
!
 10   FORMAT (/,' OUTPUT - multi-file counter ',a,i0,                   &
     &          ', is greater than Nfiles = ',i0,1x,'dimension',        &
     &        /,10x,'in structure when creating next file: ',           &
     &           a,'_',i4.4,'.nc',                                      &
     &        /,10x,'Incorrect OutFiles logic in ''read_phypar''.')
#  ifdef LONG_NUMS
 20   FORMAT (a,'_',i5.5,'.nc')
#  else
 20   FORMAT (a,'_',i4.4,'.nc')
#  endif
!
      RETURN
      END SUBROUTINE output
#else
      SUBROUTINE output
      RETURN
      END SUBROUTINE output
#endif
