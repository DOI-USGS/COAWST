#include "cppdefs.h"
      SUBROUTINE close_inp (ng, model)
!
!svn $Id: close_io.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This subroutine checks if input files are in close state.  It is     !
! used during initialization to force all multi-file input fields to   !
! in close state. This is important in iterative algorithms that run   !
! the full model repetitevely.                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
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
     &  __FILE__//", close_inp.F"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  If multi-file input fields, close several input files.
!-----------------------------------------------------------------------

#ifdef FRC_FILE
!
!  If appropriate, close input forcing files and set several parameter
!  to closed state.
!
      DO i=1,nFfiles(ng)
        IF ((FRC(i,ng)%Nfiles.gt.0).and.(FRC(i,ng)%ncid.ne.-1)) THEN
          Fcount=FRC(i,ng)%Fcount
          CALL netcdf_close (ng, model, FRC(i,ng)%ncid,                 &
     &                       FRC(i,ng)%files(Fcount), .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          FRCids=-1
          ncFRCid=-1
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
            CALL netcdf_close (ng, model, BRY(i,ng)%ncid,               &
     &                         BRY(i,ng)%files(Fcount),  .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            BRYids=-1
            ncBRYid=-1
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
            CALL netcdf_close (ng, model, CLM(i,ng)%ncid,               &
     &                         CLM(i,ng)%files(Fcount),  .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            CLMids=-1
            ncCLMid=-1
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
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#if defined BIOFLUX
      USE mod_biology
#endif
      USE strings_mod, ONLY : FoundError
!
      USE dateclock_mod, ONLY : get_date
!
      implicit none
!
!  Local variable declarations.
!
      logical :: First
!
      integer :: Fcount, MyError, i, ng
#if defined FILTERED_RST
      integer :: ifile
#endif

#if defined BIOFLUX
      integer :: j
      real(r8) :: dtdays
!
!---------------------------------------------------------------------
!  Write biological fluxes to standard out.
!---------------------------------------------------------------------
!
      
      write(stdout,*) ' '
      write(stdout,*) ' '
      write(stdout,134) ' Cumulative biological fluxes'
      DO ng=1,Ngrids
        IF (Master) THEN
# if defined BIO_GOANPZ || (defined BIOFLUX && defined BEST_NPZ)
!         Scale fluxes by maximum
          xi = 0.0_r8
          DO i = itemp,iDet
            DO j = itemp,iDet
              xi = MAX(xi, bflx(i,j))
            END DO
          END DO
          xi = 10.0_r8/xi
          write(stdout,134) ' To  :  From:'
          write(stdout,131) 'NO3','NH4','PhS','PhL','MZS','MZL','Cop',  &
     &                      'NCa','Eup','Det'
          write(stdout,132) 'NO3 : ', (xi*bflx(i,iNO3), i=iNO3,IDet)
          write(stdout,132) 'NH4 : ', (xi*bflx(i,iNH4), i=iNO3,IDet)
          write(stdout,132) 'PhS : ', (xi*bflx(i,iPhS), i=iNO3,IDet)
          write(stdout,132) 'PhL : ', (xi*bflx(i,iPhL), i=iNO3,IDet)
          write(stdout,132) 'MZS : ', (xi*bflx(i,iMZS), i=iNO3,IDet)
          write(stdout,132) 'MZL : ', (xi*bflx(i,iMZL), i=iNO3,IDet)
          write(stdout,132) 'Cop : ', (xi*bflx(i,iCop), i=iNO3,IDet)
          write(stdout,132) 'NCa : ', (xi*bflx(i,iNCa), i=iNO3,IDet)
          write(stdout,132) 'Eup : ', (xi*bflx(i,iEup), i=iNO3,IDet)
          write(stdout,132) 'Det : ', (xi*bflx(i,iDet), i=iNO3,IDet)
          write(stdout,133) 'pred: ', (xi*bflx(i,itemp), i=iCop,iEup)
          write(stdout,*) ' '
 131  format('BIOFLUX  ',7x,10a6)      
 132  format('BIOFLUX  ',a7,10f6.1)      
 133  format('BIOFLUX  ',a7,36x,3f6.1)      
 134  format('BIOFLUX  ',a)      
# endif
        END IF
      END DO
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", close_out"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (RST(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, RST(ng)%ncid, RST(ng)%name)
        END IF
#if defined FOUR_DVAR || defined ENKF_RESTART || defined VERIFICATION
        IF (DAI(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, DAI(ng)%ncid, DAI(ng)%name)
        END IF
        IF (DAV(ng)%ncid.ne.-1) THEN
          CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
          CALL netcdf_close (ng, iNLM, DAV(ng)%ncid, DAV(ng)%name)
        END IF
#endif
#if defined FORWARD_READ || defined FORWARD_WRITE
        IF ((FWD(ng)%ncid.ne.-1).and.                                   &
     &      (FWD(ng)%ncid.eq.HIS(ng)%ncid)) THEN
          FWD(ng)%ncid=-1
        END IF
        IF (FWD(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, FWD(ng)%ncid, FWD(ng)%name)
        END IF
#endif
        IF (HIS(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, HIS(ng)%ncid, HIS(ng)%name)
        END IF
        IF (QCK(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, QCK(ng)%ncid, QCK(ng)%name)
        END IF
#ifdef SP4DVAR
        IF (SPT(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, SPT(ng)%ncid, SPT(ng)%name)
        END IF
        IF (SCT(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, SCT(ng)%ncid, SCT(ng)%name)
        END IF
        IF (SPA(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, SPA(ng)%ncid, SPA(ng)%name)
        END IF
#endif
#ifdef ADJOINT
        IF (ADM(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, ADM(ng)%ncid, ADM(ng)%name)
        END IF
#endif
#ifdef TANGENT
        IF (ITL(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, ITL(ng)%ncid, ITL(ng)%name)
        END IF
        IF (TLM(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, TLM(ng)%ncid, TLM(ng)%name)
        END IF
#endif
#if defined TL_IOMS && defined FOUR_DVAR
        IF (IRP(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iRPM, IRP(ng)%ncid, IRP(ng)%name)
        END IF
#endif
#ifdef WEAK_CONSTRAINT
        IF (TLF(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, TLF(ng)%ncid, TLF(ng)%name)
        END IF
#endif
#ifdef FOUR_DVAR
        IF (HSS(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, HSS(ng)%ncid, HSS(ng)%name)
        END IF
        IF (LCZ(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, LCZ(ng)%ncid, LCZ(ng)%name)
        END IF
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        IF (AVG(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, AVG(ng)%ncid, AVG(ng)%name)
        END IF
#endif
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
        IF (HAR(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, HAR(ng)%ncid, HAR(ng)%name)
        END IF
#endif
#ifdef FILTERED_RST
        DO ifile=1,nfile
          IF (FIL(ifile,ng)%ncid.ne.-1) THEN
            CALL netcdf_close (ng, iNLM, FIL(ifile,ng)%ncid)
          END IF
        END DO
#endif
#ifdef DIAGNOSTICS
        IF (DIA(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, DIA(ng)%ncid, DIA(ng)%name)
        END IF
#endif
#ifdef FLOATS
        IF (FLT(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, FLT(ng)%ncid, FLT(ng)%name)
        END IF
#endif
#ifdef STATIONS
        IF (STA(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, STA(ng)%ncid)
        END IF
#endif
#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
        IF (ERR(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, ERR(ng)%ncid, ERR(ng)%name)
        END IF
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
#ifdef FILTERED_RST
          DO ifile=1,nfile
            IF (associated(FIL(ifile,ng)%Nrec)) THEN
              Fcount=FIL(ifile,ng)%Fcount
              IF (FIL(ifile,ng)%Nrec(Fcount).gt.0) THEN
                WRITE (stdout,20) 'FILTER RESTART',                     &
     &                     FIL(ifile,ng)%Nrec(Fcount)
              END IF
            END IF
          END DO
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
        DO i=1,51
          IF ((LEN_TRIM(ANANAME(i)).gt.0).and.(exit_flag.ne.5)) THEN
            IF (First) THEN
              First=.FALSE.
              WRITE (stdout,30) ' Analytical header files used:'
            END IF
            WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(ANANAME(i)))
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
          IF ((LEN_TRIM(BIONAME(i)).gt.0).and.(exit_flag.ne.5)) THEN
            IF (First) THEN
              First=.FALSE.
              WRITE (stdout,30) ' Biology model header files used:'
            END IF
            WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(BIONAME(i)))
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
#ifndef NO_DATE_STAMP
        CALL get_date (date_str)
#endif
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
