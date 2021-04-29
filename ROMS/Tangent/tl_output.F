#include "cppdefs.h"
#ifdef TANGENT
      SUBROUTINE tl_output (ng)
!
!svn $Id: tl_output.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine manages tangent linear model output. It creates output !
!  NetCDF files and writes out data into NetCDF files.  If  requested, !
!  it can create several tangent history files to avoid generating too !
!  large files during a single model run.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FOUR_DVAR
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
      CALL wclock_on (ng, iTLM, 8, __LINE__, MyFile)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process tangent linear history NetCDF file.
!-----------------------------------------------------------------------
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  Create output tangent NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefTLM(ng)) THEN
        IF (ndefTLM(ng).gt.0) THEN
          IF (idefTLM(ng).lt.0) THEN
            idefTLM(ng)=((ntstart(ng)-1)/ndefTLM(ng))*ndefTLM(ng)
            IF (idefTLM(ng).lt.iic(ng)-1) THEN
              idefTLM(ng)=idefTLM(ng)+ndefTLM(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefTLM(ng)) THEN
              TLM(ng)%load=0                  ! restart, reset counter
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              Ldefine=.TRUE.
              NewFile=.FALSE.                 ! unfinished file, inquire
            END IF                            ! content for appending
            idefTLM(ng)=idefTLM(ng)+nTLM(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefTLM(ng)) THEN
            idefTLM(ng)=idefTLM(ng)+ndefTLM(ng)
            IF (nTLM(ng).ne.ndefTLM(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefTLM(ng)=idefTLM(ng)+nTLM(ng)  ! multiple record offset
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN                     ! create new file or
            IF (iic(ng).eq.ntstart(ng)) THEN    ! inquire existing file
              TLM(ng)%load=0                    ! reset filename counter
            END IF
            ifile=(iic(ng)-1)/ndefTLM(ng)+1     ! next filename suffix
            TLM(ng)%load=TLM(ng)%load+1
            IF (TLM(ng)%load.gt.TLM(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'TLM(ng)%load = ', TLM(ng)%load,      &
     &                            TLM(ng)%Nfiles, TRIM(TLM(ng)%base),   &
     &                            ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            Fcount=TLM(ng)%load
            TLM(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (TLM(ng)%name,20) TRIM(TLM(ng)%base), ifile
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iTLM, TLM(ng)%name)
# endif
            TLM(ng)%files(Fcount)=TRIM(TLM(ng)%name)
            IF (TLM(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
            END IF
            CALL tl_def_his (ng, NewFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtTLM(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtTLM(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL tl_def_his (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtTLM(ng)=.TRUE.
            LdefTLM(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into tangent NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtTLM(ng)) THEN
# if defined STOCHASTIC_OPT || defined FORCING_SV || \
     defined HESSIAN_SO     || defined HESSIAN_FSV
!
!  Write only on first time step.
!
        IF (iic(ng).eq.1) THEN
          CALL tl_wrt_his (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# else
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nTLM(ng)).eq.0)) THEN
            CALL tl_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        ELSE
          IF ((MOD(iic(ng)-1,nTLM(ng)).eq.0).and.                       &
     &        ((nrrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
            CALL tl_wrt_his (ng)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
# endif
      END IF

# ifdef TL_AVERAGES
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
            AVG(ng)%load=AVG(ng)%load+1
            IF (ndefAVG(ng).eq.nAVG(ng)) THEN ! next filename suffix
              ifile=(iic(ng)-1)/ndefAVG(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefAVG(ng)+1
            END IF
            IF (AVG(ng)%load.gt.AVG(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'AVG(ng)%load = ', AVG(ng)%load,      &
     &                            AVG(ng)%Nfiles, TRIM(AVG(ng)%base),   &
     &                            ifile
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
            CALL mp_bcasts (ng, iTLM, AVG(ng)%name)
#  endif
            AVG(ng)%files(Fcount)=TRIM(AVG(ng)%name)
            IF (AVG(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iTLM, AVG(ng)%ncid)
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
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nAVG(ng)).eq.0)) THEN
          CALL wrt_avg (ng)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
# endif

# if defined FOUR_DVAR
#  ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process and write model state at observation
!  locations. Compute misfit (model-observations) cost function.
!-----------------------------------------------------------------------
!
        IF (((time(ng)-0.5_r8*dt(ng)).le.ObsTime(ng)).and.              &
     &      (ObsTime(ng).lt.(time(ng)+0.5_r8*dt(ng)))) THEN
          ProcessObs(ng)=.TRUE.
#   ifdef DISTRIBUTE
          tile=MyRank
#   else
          tile=-1
#   endif
          CALL obs_read (ng, iTLM, .FALSE.)
#   ifdef SP4DVAR
!
!  Do not assimilate obs on the last time step unless inter=Nsaddle.
!
          IF (iic(ng).ne.ntend(ng)+1.or.Lsadd(ng)) THEN
            CALL obs_write (ng, tile, iTLM)
          END IF
#   else
          CALL obs_write (ng, tile, iTLM)
#   endif
#   ifndef WEAK_CONSTRAINT
          CALL obs_cost (ng, iTLM)
#   endif
        ELSE
          ProcessObs(ng)=.FALSE.
        END IF
#  endif
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iTLM, 8, __LINE__, MyFile)
# endif
!
 10   FORMAT (/,' TL_OUTPUT - multi-file counter ',a,i0,                &
     &          ', is greater than Nfiles = ',i0,1x,'dimension',        &
     &        /,13x,'in structure when creating next file: ',           &
     &           a,'_',i4.4,'.nc',                                      &
     &        /,13x,'Incorrect OutFiles logic in ''read_phypar''.')
 20   FORMAT (a,'_',i4.4,'.nc')
!
      RETURN
      END SUBROUTINE tl_output
#else
      SUBROUTINE tl_output
      RETURN
      END SUBROUTINE tl_output
#endif
