#include "cppdefs.h"
      MODULE close_io_mod
!
!git $Id$
!svn $Id: close_io.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module closes input and output files using either the standard !
!  NetCDF library or the Parallel-IO (PIO) library.                    !
!                                                                      !
!  During initialization, the input input files need to be in closed   !
!  state to facilitate multi-file processing. This is important in     !
!  iterative algorithms that run the model kernels repetitevely.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
!
      USE dateclock_mod, ONLY : get_date
      USE strings_mod,   ONLY : FoundError
!
      implicit none
!
      PUBLIC :: close_file
      PUBLIC :: close_inp
      PUBLIC :: close_out
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE close_file (ng, model, S, ncname, Lupdate)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      TYPE(T_IO), intent(inout) :: S
!
      logical, intent(in), optional :: Lupdate
!
      character (len=*), intent(in), optional :: ncname
!
!  Local variable declarations.
!
      integer :: ClosedState = -1
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", close_file_nf90"
!
!-----------------------------------------------------------------------
!  Close specified NetCDF file.
!-----------------------------------------------------------------------
!
      SELECT CASE (S%IOtype)
        CASE (io_nf90)
          IF (S%ncid.ne.ClosedState) THEN
            CALL netcdf_close (ng, model, S%ncid,                       &
     &                         TRIM(ncname), Lupdate)
            S%ncid=ClosedState
          END IF

#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          IF (ASSOCIATED(S%pioFile%iosystem)) THEN
            IF (S%pioFile%fh.ne.ClosedState) THEN
              CALL pio_netcdf_close (ng, model, S%pioFile,              &
     &                               TRIM(ncname), Lupdate)
              S%pioFile%fh=ClosedState
            END IF
          END IF
#endif
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE close_file
!
!***********************************************************************
      SUBROUTINE close_inp (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: Fcount, i
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", close_inp"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  If multi-file input fields, close several input files.
!-----------------------------------------------------------------------
!
!  Skip if configuration error.
!
      IF ((exit_flag.eq.5).or.(exit_flag.eq.6)) RETURN

#ifdef FRC_FILE
!
!  If appropriate, close input forcing files and set several parameter
!  to closed state.
!
      DO i=1,nFfiles(ng)
        IF ((FRC(i,ng)%Nfiles.gt.0).and.(FRC(i,ng)%ncid.ne.-1)) THEN
          Fcount=FRC(i,ng)%Fcount
          CALL close_file (ng, model, FRC(i,ng),                        &
     &                     FRC(i,ng)%files(Fcount), .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          FRCids=-1
          FRCncid=-1
          Fcount=1
          FRC(i,ng)%Fcount=Fcount
          FRC(i,ng)%name=TRIM(FRC(i,ng)%files(Fcount))
        END IF
      END DO
#endif
!
!  If appropriate, close boundary files.
!
      IF (ObcData(ng)) THEN
        DO i=1,nBCfiles(ng)
          IF ((BRY(i,ng)%Nfiles.gt.0).and.(BRY(i,ng)%ncid.ne.-1)) THEN
            Fcount=BRY(i,ng)%Fcount
            CALL close_file (ng, model, BRY(i,ng),                      &
     &                       BRY(i,ng)%files(Fcount), .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            BRYids=-1
            BRYncid=-1
            Fcount=1
            BRY(i,ng)%Fcount=Fcount
            BRY(i,ng)%name=TRIM(BRY(i,ng)%files(Fcount))
          END IF
        END DO
      END IF
!
!  If appropriate, close climatology files.
!
      IF (CLM_FILE(ng)) THEN
        DO i=1,nCLMfiles(ng)
          IF ((CLM(i,ng)%Nfiles.gt.0).and.(CLM(i,ng)%ncid.ne.-1)) THEN
            Fcount=CLM(i,ng)%Fcount
            CALL close_file (ng, model, CLM(i,ng),                      &
     &                       CLM(i,ng)%files(Fcount), .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            CLMids=-1
            CLMncid=-1
            Fcount=1
            CLM(i,ng)%Fcount=Fcount
            CLM(i,ng)%name=TRIM(CLM(i,ng)%files(Fcount))
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE close_inp
!
      SUBROUTINE close_out
!
!=======================================================================
!                                                                      !
! This subroutine flushes and closes all output files.                 !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      logical :: First, Lupdate
!
      integer :: Fcount, MyError, i, ivalue, ng
!
      character (len=256) :: ana_string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", close_out"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
!  Skip if configuration error.
!
      IF ((exit_flag.eq.5).or.(exit_flag.eq.6)) RETURN
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
!  Close output NetCDF files.
!
      DO ng=1,Ngrids
        CALL close_file (ng, iNLM, RST(ng), RST(ng)%name, Lupdate)
#if defined FOUR_DVAR || defined ENKF_RESTART || defined VERIFICATION
        CALL close_file (ng, iNLM, DAI(ng), DAI(ng)%name, Lupdate)
        CALL close_file (ng, iNLM, DAV(ng), DAV(ng)%name, Lupdate)
#endif
#if defined FORWARD_READ || defined FORWARD_WRITE
        IF (FWD(ng)%IOtype.eq.io_nf90) THEN
          IF ((FWD(ng)%ncid.ne.-1).and.                                 &
     &        (FWD(ng)%ncid.eq.HIS(ng)%ncid)) THEN
            FWD(ng)%ncid=-1
          END IF
# if defined PIO_LIB && defined DISTRIBUTE
        ELSE IF (FWD(ng)%IOtype.eq.io_pio) THEN
          IF ((FWD(ng)%pioFile%fh.ne.-1).and.                           &
     &        (FWD(ng)%pioFile%fh.eq.HIS(ng)%pioFile%fh)) THEN
            FWD(ng)%pioFile%fh=-1
          END IF
# endif
        END IF
        CALL close_file (ng, iNLM, FWD(ng), FWD(ng)%name, Lupdate)
#endif
        CALL close_file (ng, iNLM, HIS(ng), HIS(ng)%name, Lupdate)
        CALL close_file (ng, iNLM, QCK(ng), QCK(ng)%name, Lupdate)
#ifdef SP4DVAR
        CALL close_file (ng, iTLM, SPT(ng), SPT(ng)%name, Lupdate)
        CALL close_file (ng, iTLM, SCT(ng), SCT(ng)%name, Lupdate)
        CALL close_file (ng, iADM, SPA(ng), SPA(ng)%name, Lupdate)
#endif
#ifdef ADJOINT
        CALL close_file (ng, iADM, ADM(ng), ADM(ng)%name, Lupdate)
#endif
#ifdef TANGENT
        CALL close_file (ng, iTLM, ITL(ng), ITL(ng)%name, Lupdate)
        CALL close_file (ng, iTLM, TLM(ng), TLM(ng)%name, Lupdate)
#endif
#if defined TL_IOMS && defined FOUR_DVAR
        CALL close_file (ng, iRPM, IRP(ng), IRP(ng)%name, Lupdate)
#endif
#ifdef WEAK_CONSTRAINT
        CALL close_file (ng, iTLM, TLF(ng), TLF(ng)%name, Lupdate)
#endif
#ifdef FOUR_DVAR
        CALL close_file (ng, iADM, HSS(ng), HSS(ng)%name, Lupdate)
        CALL close_file (ng, iADM, LCZ(ng), LCZ(ng)%name, Lupdate)
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        CALL close_file (ng, iNLM, AVG(ng), AVG(ng)%name, Lupdate)
#endif
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
        CALL close_file (ng, iNLM, HAR(ng), HAR(ng)%name, Lupdate)
#endif
#ifdef DIAGNOSTICS
        CALL close_file (ng, iNLM, DIA(ng), DIA(ng)%name, Lupdate)
#endif
#ifdef FLOATS
        CALL close_file (ng, iNLM, FLT(ng), FLT(ng)%name, Lupdate)
#endif
#ifdef STATIONS
        CALL close_file (ng, iNLM, STA(ng), STA(ng)%name, Lupdate)
#endif
#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
        CALL close_file (ng, iTLM, ERR(ng), ERR(ng)%name, Lupdate)
#endif
!
!  Report number of time records written.
!
        IF (Master) THEN
          WRITE (stdout,10) ng

          IF (associated(HIS(ng)%Nrec)) THEN
            IF (ANY(HIS(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'HISTORY', SUM(HIS(ng)%Nrec)
            END IF
          END IF

          IF (associated(RST(ng)%Nrec)) THEN
            Fcount=RST(ng)%load
            IF (RST(ng)%Nrec(Fcount).gt.0) THEN
              IF (LcycleRST(ng)) THEN
                IF (RST(ng)%Nrec(Fcount).gt.1) THEN
                  RST(ng)%Nrec(Fcount)=2
                ELSE
                  RST(ng)%Nrec(Fcount)=1
                END IF
              END IF
              WRITE (stdout,20) 'RESTART', RST(ng)%Nrec(Fcount)
            END IF
          END IF

#if defined FOUR_DVAR || defined ENKF_RESTART
          IF (associated(DAI(ng)%Nrec)) THEN
            IF (ANY(DAI(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'DA IC  ', SUM(DAI(ng)%Nrec)
            END IF
          END IF
#endif

#ifdef ADJOINT
          IF (associated(ADM(ng)%Nrec)) THEN
            IF (ANY(ADM(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'ADJOINT', SUM(ADM(ng)%Nrec)
            END IF
          END IF
#endif

#ifdef TANGENT
          IF (associated(TLM(ng)%Nrec)) THEN
            IF (ANY(TLM(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'TANGENT', SUM(TLM(ng)%Nrec)
            END IF
          END IF
#endif

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          IF (associated(AVG(ng)%Nrec)) THEN
            IF (ANY(AVG(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'AVERAGE', SUM(AVG(ng)%Nrec)
            END IF
          END IF
#endif

#ifdef STATIONS
          IF (associated(STA(ng)%Nrec)) THEN
            IF (ANY(STA(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'STATION', SUM(STA(ng)%Nrec)
            END IF
          END IF
#endif

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
          IF (associated(ERR(ng)%Nrec)) THEN
            IF (ANY(ERR(ng)%Nrec.gt.0)) THEN
              WRITE (stdout,20) 'ERROR  ', SUM(ERR(ng)%Nrec)
            END IF
          END IF
#endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report analytical header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,39
          ana_string=TRIM(ANANAME(i))
          ivalue=ICHAR(ana_string(1:1))
          IF (ivalue.ge.47) THEN          ! decimal value for characters
            IF (exit_flag.ne.5) THEN
              IF (First) THEN
                First=.FALSE.
                WRITE (stdout,30) ' Analytical header files used:'
              END IF
              WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(ANANAME(i)))
            END IF
          END IF
        END DO
      END IF

#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Report biology model header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,4
          ana_string=TRIM(BIONAME(i))
          ivalue=ICHAR(ana_string(1:1))
          IF (ivalue.ge.47) THEN          ! decimal value for characters
            IF (exit_flag.ne.5) THEN
              IF (First) THEN
                First=.FALSE.
                WRITE (stdout,30) ' Biology model header files used:'
              END IF
              WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(BIONAME(i)))
            END IF
          END IF
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  If applicable, report internal exit errors.
!-----------------------------------------------------------------------
!
      IF (Master.and.                                                   &
     &    (FoundError(exit_flag, NoError, __LINE__, MyFile))) THEN
        WRITE (stdout,40) Rerror(exit_flag), exit_flag
      END IF
      IF (blowup.ne.0) THEN
        IF (Master) WRITE (stdout,50) TRIM(blowup_string)
      ELSE IF (exit_flag.eq.NoError) THEN
        CALL get_date (date_str)
        IF (Master) WRITE (stdout,60) TRIM(date_str)
      ELSE IF (exit_flag.eq.2) THEN
        IF (Master) WRITE (stdout,70) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.3) THEN
        IF (Master) WRITE (stdout,80) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.4) THEN
        IF (Master) WRITE (stdout,90)
      ELSE IF (exit_flag.eq.5) THEN
        IF (Master) WRITE (stdout,100)
      ELSE IF (exit_flag.eq.6) THEN
        IF (Master) WRITE (stdout,110)
      ELSE IF (exit_flag.eq.7) THEN
        IF (Master) WRITE (stdout,120)
      ELSE IF (exit_flag.eq.8) THEN
        IF (Master) WRITE (stdout,130)
      END IF

#ifdef ROMS_STDOUT
!
!-----------------------------------------------------------------------
!  Close ROMS standard outpu file.
!-----------------------------------------------------------------------
!
      CALL my_flush (stdout)
      CLOSE (stdout)
#endif
!
 10   FORMAT (/,' ROMS/TOMS - Output NetCDF summary for Grid ',         &
     &        i2.2,':')
 20   FORMAT (13x,'number of time records written in ',                 &
     &        a,' file = ',i0)
 30   FORMAT (/,a,/)
 40   FORMAT (/,a,i0,/)
 50   FORMAT (/,' MAIN: Abnormal termination: BLOWUP.',/,               &
     &          ' REASON: ',a)
 60   FORMAT (/,' ROMS/TOMS: DONE... ',a)
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.',/,        &
     &          ' REASON: ',a)
 80   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.',/,       &
     &          ' REASON: ',a)
 90   FORMAT (/,' ERROR: I/O related problem.')
100   FORMAT (/,' ERROR: Illegal model configuration.')
110   FORMAT (/,' ERROR: Illegal domain partition.')
120   FORMAT (/,' ERROR: Illegal input parameter.')
130   FORMAT (/,' ERROR: Fatal algorithm result.')
!
      RETURN
      END SUBROUTINE close_out

      END MODULE close_io_mod
