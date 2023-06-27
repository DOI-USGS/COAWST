#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_output (ng)
!
!git $Id$
!svn $Id: ad_output.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine manages adjoint model output. It creates output NetCDF !
!  files and writes out data into NetCDF files.  If requested,  it can !
!  create several adjoint history files  to avoid generating too large !
!  files during a single model run.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE ad_def_his_mod, ONLY : ad_def_his
      USE ad_wrt_his_mod, ONLY : ad_wrt_his
      USE close_io_mod,   ONLY : close_file
# ifdef AD_AVERAGES
      USE def_avg_mod,    ONLY : def_avg
# endif
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasts
# endif
      USE strings_mod,    ONLY : FoundError
# ifdef AD_AVERAGES
      USE wrt_avg_mod,    ONLY : wrt_avg
# endif
!
      implicit none
!
!  Imported variables declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, Lupdate, NewFile, wrtHIS

      integer :: Fcount, ifile, status
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
      CALL wclock_on (ng, iADM, 8, __LINE__, MyFile)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process adjoint history NetCDF file.
!-----------------------------------------------------------------------
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  If appropriate, set switch for updating biology header file global
!  attribute in output NetCDF files.
!
#ifdef BIOLOGY
      Lupdate=.TRUE.
#else
      Lupdate=.FALSE.
#endif
!
!  Create output adjoint NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefADJ(ng)) THEN
        IF (ndefADJ(ng).gt.0) THEN
          IF (idefADJ(ng).lt.0) THEN
            idefADJ(ng)=((ntstart(ng)-1)/ndefADJ(ng))*ndefADJ(ng)
            IF (idefADJ(ng).lt.iic(ng)-1) THEN
              idefADJ(ng)=idefADJ(ng)+ndefADJ(ng)
            END IF
          END IF
          IF ((nrrec(ng).ne.0).and.(iic(ng).eq.ntstart(ng))) THEN
            IF ((iic(ng)-1).eq.idefADJ(ng)) THEN
              ADM(ng)%load=0                  ! restart, reset counter
              Ldefine=.FALSE.                 ! finished file, delay
            ELSE                              ! creation of next file
              Ldefine=.TRUE.
              NewFile=.FALSE.                 ! unfinished file, inquire
            END IF                            ! content for appending
            idefADJ(ng)=idefADJ(ng)+nADJ(ng)  ! restart offset
          ELSE IF ((iic(ng)-1).eq.idefADJ(ng)) THEN
            idefADJ(ng)=idefADJ(ng)+ndefADJ(ng)
            IF (nADJ(ng).ne.ndefADJ(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefADJ(ng)=idefADJ(ng)+nADJ(ng)  ! multiple record offset
            END IF
            Ldefine=.TRUE.
            NewFile=.TRUE.
          ELSE
            Ldefine=.FALSE.
          END IF
          IF (Ldefine) THEN                     ! create new file or
            IF (iic(ng).eq.ntstart(ng)) THEN    ! inquire existing file
              ADM(ng)%load=0                    ! reset filename counter
            END IF
            ifile=(iic(ng)-1)/ndefADJ(ng)+1     ! next filename suffix
            ADM(ng)%load=ADM(ng)%load+1
            IF (ADM(ng)%load.gt.ADM(ng)%Nfiles) THEN
              IF (Master) THEN
                WRITE (stdout,10) 'TLM(ng)%load = ', ADM(ng)%load,      &
     &                            ADM(ng)%Nfiles, TRIM(ADM(ng)%base),   &
     &                            ifile
              END IF
              exit_flag=4
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            Fcount=ADM(ng)%load
            ADM(ng)%Nrec(Fcount)=0
            IF (Master) THEN
              WRITE (ADM(ng)%name,20) TRIM(ADM(ng)%base), ifile
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iADM, ADM(ng)%name)
# endif
            ADM(ng)%files(Fcount)=TRIM(ADM(ng)%name)
            CALL close_file (ng, iADM, ADM(ng), ADM(ng)%name, Lupdate)
            CALL ad_def_his (ng, NewFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtADJ(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtADJ(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL ad_def_his (ng, ldefout(ng))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            LwrtADJ(ng)=.TRUE.
            LdefADJ(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into adjoint NetCDF file.
!
      IF (LwrtADJ(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).ne.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nADJ(ng)).eq.0)) THEN
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        ELSE
# ifdef HESSIAN_SO
          wrtHIS=(MOD(iic(ng)-1,nADJ(ng)).eq.0)   ! otherwise
# else
          IF (nADJ(ng).eq.ntimes(ng)) THEN
            wrtHIS=(iic(ng).ne.ntstart(ng)).and.                        &
     &             (MOD(iic(ng)-1,nADJ(ng)).eq.0) ! avoid ntstart rec
          ELSE
#  if defined WEAK_CONSTRAINT && !defined WEAK_NOINTERP
#   ifdef SP4DVAR
            wrtHIS=(MOD(iic(ng)-1,nADJ(ng)).eq.0)
#   endif
            wrtHIS=(iic(ng).ne.ntstart(ng)).and.                        &
     &             (MOD(iic(ng)-1,nADJ(ng)).eq.0) ! avoid ntstart-1 rec
#  else
            wrtHIS=(MOD(iic(ng)-1,nADJ(ng)).eq.0) ! otherwise
#  endif
          END IF
# endif
          IF (wrtHIS) THEN
# ifdef DISTRIBUTE
            CALL ad_wrt_his (ng, MyRank)
# else
            CALL ad_wrt_his (ng, -1)
# endif
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
      END IF

# ifdef AD_AVERAGES
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
            CALL mp_bcasts (ng, iADM, AVG(ng)%name)
#  endif
            AVG(ng)%files(Fcount)=TRIM(AVG(ng)%name)
            CALL close_file (ng, iADM, ADM(ng), ADM(ng)%name, Lupdate)
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
        IF ((iic(ng).ne.ntstart(ng)).and.                               &
     &      (MOD(iic(ng),nAVG(ng)).eq.1)) THEN
#  ifdef DISTRIBUTE
          CALL wrt_avg (ng, MyRank)
#  else
          CALL wrt_avg (ng, -1)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iADM, 8, __LINE__, MyFile)
# endif
!
 10   FORMAT (/,' AD_OUTPUT - multi-file counter ',a,i0,                &
     &          ', is greater than Nfiles = ',i0,1x,'dimension',        &
     &        /,13x,'in structure when creating next file: ',           &
     &           a,'_',i4.4,'.nc',                                      &
     &        /,13x,'Incorrect OutFiles logic in ''read_phypar''.')
 20   FORMAT (a,'_',i4.4,'.nc')
!
      RETURN
      END SUBROUTINE ad_output
#else
      SUBROUTINE ad_output
      RETURN
      END SUBROUTINE ad_output
#endif
