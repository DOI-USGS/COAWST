#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_output (ng)
!
!svn $Id: ad_output.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasts
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
      logical :: Ldefine, NewFile, wrtHIS

      integer :: Fcount, ifile, status
!
      SourceFile='ad_output.F'

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iADM, 8)
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
            Fcount=ADM(ng)%Fcount               ! inquire existing file
            ADM(ng)%Nrec(Fcount)=0
            ifile=(iic(ng)-1)/ndefADJ(ng)+1
            IF (Master) THEN
              WRITE (ADM(ng)%name,10) TRIM(ADM(ng)%base), ifile
  10          FORMAT (a,'_',i4.4,'.nc')
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iADM, ADM(ng)%name)
# endif
            IF (ADM(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iADM, ADM(ng)%ncid)
            END IF
            CALL ad_def_his (ng, NewFile)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          IF ((iic(ng).eq.ntstart(ng)).and.(nrrec(ng).ne.0)) THEN
            LwrtADJ(ng)=.FALSE.                 ! avoid writing initial
          ELSE                                  ! fields during restart
            LwrtADJ(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL ad_def_his (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
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
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
# ifdef HESSIAN_SO
          wrtHIS=(MOD(iic(ng)-1,nADJ(ng)).eq.0)   ! otherwise
# else
          IF (nADJ(ng).eq.ntimes(ng)) THEN
            wrtHIS=(iic(ng).ne.ntstart(ng)).and.                        &
     &             (MOD(iic(ng)-1,nADJ(ng)).eq.0) ! avoid ntstart rec
          ELSE
#  ifdef WEAK_CONSTRAINT
            wrtHIS=(iic(ng).ne.ntstart(ng)).and.                        &
     &             (MOD(iic(ng)-1,nADJ(ng)).eq.0) ! avoid ntstart-1 rec
#  else
            wrtHIS=(MOD(iic(ng)-1,nADJ(ng)).eq.0) ! otherwise
#  endif
          END IF
# endif
          IF (wrtHIS) THEN
            CALL ad_wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
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
            AVG(ng)%Nrec(Fcount)=0
            IF (ndefAVG(ng).eq.nAVG(ng)) THEN
              ifile=(iic(ng)-1)/ndefAVG(ng)
            ELSE
              ifile=(iic(ng)-1)/ndefAVG(ng)+1
            END IF
            IF (Master) THEN
              WRITE (AVG(ng)%name,20) TRIM(AVG(ng)%base), ifile
  20          FORMAT (a,'_',i4.4,'.nc')
            END IF
#  ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iADM, AVG(ng)%name)
#  endif
            IF (AVG(ng)%ncid.ne.-1) THEN
              CALL netcdf_close (ng, iADM, AVG(ng)%ncid)
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
        IF ((iic(ng).ne.ntstart(ng)).and.                               &
     &      (MOD(iic(ng),nAVG(ng)).eq.1)) THEN
          CALL wrt_avg (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iADM, 8)
# endif
      RETURN
      END SUBROUTINE ad_output
#else
      SUBROUTINE ad_output
      RETURN
      END SUBROUTINE ad_output
#endif
