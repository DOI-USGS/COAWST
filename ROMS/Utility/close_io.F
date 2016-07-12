#include "cppdefs.h"
      SUBROUTINE close_inp (ng, model)
!
!svn $Id: close_io.F 798 2016-05-21 17:04:58Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This subroutine checks some input files are in close state.  It is   !
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
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: Fcount, i, j, lstr
!
      SourceFile='close_inp, close_io.F'
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
        IF ((FRC(i,ng)%Nfiles.gt.1).and.(FRC(i,ng)%ncid.ge.0)) THEN
          FRCids(i,ng)=-1
          DO j=1,NV
            IF (ncFRCid(j,ng).eq.FRC(i,ng)%ncid) THEN
              ncFRCid(j,ng)=-1
            END IF
          END DO
          IF (model.eq.iADM) THEN
            DO j=1,FRC(i,ng)%Nfiles
              IF ((FRC(i,ng)%time_min(j).le.tdays(ng)).and.              &
     &            (tdays(ng).le.FRC(i,ng)%time_max(j))) THEN
                Fcount=j
                EXIT
              END IF
            END DO
          ELSE
            Fcount=1
          END IF
          FRC(i,ng)%Fcount=Fcount
          FRC(i,ng)%name=TRIM(FRC(i,ng)%files(Fcount))
          lstr=LEN_TRIM(FRC(i,ng)%name)
          FRC(i,ng)%base=FRC(i,ng)%name(1:lstr-3)
          CALL netcdf_close (ng, model, FRC(i,ng)%ncid,                 &
     &                       FRC(i,ng)%files(i), .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
#endif
!
!  If appropriate, close boundary files.
!
      IF (ObcData(ng)) THEN
        DO i=1,nBCfiles(ng)
          IF ((BRY(i,ng)%Nfiles.gt.1).and.(BRY(i,ng)%ncid.ge.0)) THEN
            IF (model.eq.iADM) THEN
              DO j=1,BRY(i,ng)%Nfiles
                IF ((BRY(i,ng)%time_min(j).le.tdays(ng)).and.           &
     &              (tdays(ng).le.BRY(i,ng)%time_max(j))) THEN
                  Fcount=j
                  EXIT
                END IF
              END DO
            ELSE
              Fcount=1
            END IF
            BRY(i,ng)%Fcount=Fcount
            BRY(i,ng)%name=TRIM(BRY(i,ng)%files(Fcount))
            lstr=LEN_TRIM(BRY(i,ng)%name)
            BRY(i,ng)%base=BRY(i,ng)%name(1:lstr-3)
            CALL netcdf_close (ng, model, BRY(i,ng)%ncid,               &
     &                         BRY(i,ng)%files(i),  .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
      END IF
!
!  If appropriate, close climatology files.
!
      IF (CLM_FILE(ng)) THEN
        DO i=1,nCLMfiles(ng)
          IF ((CLM(i,ng)%Nfiles.gt.1).and.(CLM(i,ng)%ncid.ge.0)) THEN
            IF (model.eq.iADM) THEN
              DO j=1,CLM(i,ng)%Nfiles
                IF ((CLM(i,ng)%time_min(i).le.tdays(ng)).and.             &
     &            (tdays(ng).le.CLM(i,ng)%time_max(i))) THEN
                  Fcount=j
                  EXIT
                END IF
              END DO
            ELSE
              Fcount=1
            END IF
            CLM(i,ng)%Fcount=Fcount
            CLM(i,ng)%name=TRIM(CLM(i,ng)%files(Fcount))
            lstr=LEN_TRIM(CLM(i,ng)%name)
            CLM(i,ng)%base=CLM(i,ng)%name(1:lstr-3)
            CALL netcdf_close (ng, model, CLM(i,ng)%ncid,                 &
     &                     CLM(i,ng)%files(i),  .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
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
!
      implicit none
!
!  Local variable declarations.
!
      logical :: First

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
      SourceFile='close_out, close_io.F'
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (RST(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, RST(ng)%ncid)
        END IF
#if defined FOUR_DVAR || defined ENKF_RESTART
        IF (DAI(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, DAI(ng)%ncid)
        END IF
#endif
#if defined FORWARD_READ || defined FORWARD_WRITE
        IF ((FWD(ng)%ncid.ne.-1).and.(FWD(ng)%ncid.eq.HIS(ng)%ncid)) THEN
          FWD(ng)%ncid=-1
        END IF
        IF (FWD(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, FWD(ng)%ncid)
        END IF
#endif
        IF (HIS(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, HIS(ng)%ncid)
        END IF
#ifdef ADJOINT
        IF (ADM(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iADM, ADM(ng)%ncid)
        END IF
#endif
#ifdef TANGENT
        IF (TLM(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iTLM, TLM(ng)%ncid)
        END IF
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        IF (AVG(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, AVG(ng)%ncid)
        END IF
#endif
#ifdef AVERAGES2
        IF (AVG2(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, AVG2(ng)%ncid)
        END IF
#endif
#ifdef HISTORY2
        IF (HIS2(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, HIS2(ng)%ncid)
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
          CALL netcdf_close (ng, iNLM, DIA(ng)%ncid)
        END IF
#endif
#ifdef FLOATS
        IF (FLT(ng)%ncid.ne.-1) THEN
          CALL netcdf_close (ng, iNLM, FLT(ng)%ncid)
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
          CALL netcdf_close (ng, iTLM, ERR(ng)%ncid)
        END IF
#endif
!
!  Report number of time records written.
!
        IF (Master) THEN
          WRITE (stdout,10) ng

          IF (associated(HIS(ng)%Nrec)) THEN
            Fcount=HIS(ng)%Fcount
            IF (HIS(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'HISTORY', HIS(ng)%Nrec(Fcount)
            END IF
          END IF

          IF (associated(RST(ng)%Nrec)) THEN
            Fcount=RST(ng)%Fcount
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
            Fcount=DAI(ng)%Fcount
            IF (DAI(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'DA IC  ', DAI(ng)%Nrec(Fcount)
            END IF
          END IF
#endif

#ifdef ADJOINT
          IF (associated(ADM(ng)%Nrec)) THEN
            Fcount=ADM(ng)%Fcount
            IF (ADM(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'ADJOINT', ADM(ng)%Nrec(Fcount)
            END IF
          END IF
#endif

#ifdef TANGENT
          IF (associated(TLM(ng)%Nrec)) THEN
            Fcount=TLM(ng)%Fcount
            IF (TLM(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'TANGENT', TLM(ng)%Nrec(Fcount)
            END IF
          END IF
#endif

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          IF (associated(AVG(ng)%Nrec)) THEN
            Fcount=AVG(ng)%Fcount
            IF (AVG(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'AVERAGE', AVG(ng)%Nrec(Fcount)
            END IF
          END IF
#endif
#ifdef AVERAGES2
          IF (associated(AVG2(ng)%Nrec)) THEN
            Fcount=AVG2(ng)%Fcount
            IF (AVG2(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'AVERAGE2', AVG2(ng)%Nrec(Fcount)
            END IF
          END IF
#endif
#ifdef HISTORY2
          IF (associated(HIS2(ng)%Nrec)) THEN
            Fcount=HIS2(ng)%Fcount
            IF (HIS2(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'HISTORY2', HIS2(ng)%Nrec(Fcount)
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
            Fcount=STA(ng)%Fcount
            IF (STA(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'STATION', STA(ng)%Nrec(Fcount)
            END IF
          END IF
#endif

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
          IF (associated(ERR(ng)%Nrec)) THEN
            Fcount=ERR(ng)%Fcount
            IF (ERR(ng)%Nrec(Fcount).gt.0) THEN
              WRITE (stdout,20) 'ERROR  ', ERR(ng)%Nrec(Fcount)
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
      IF (Master.and.(exit_flag.ne.NoError)) THEN
        WRITE (stdout,40) Rerror(exit_flag), exit_flag
      END IF
      IF (exit_flag.eq.NoError) THEN
#ifndef NO_DATE_STAMP
        CALL get_date (date_str)
#endif
        IF (Master) WRITE (stdout,50) TRIM(date_str)
      ELSE IF ((exit_flag.eq.1).or.(blowup.ne.0)) THEN
        IF (Master) WRITE (stdout,60)
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
      ELSE IF (exit_flag.eq.9) THEN
        IF (Master) WRITE (stdout,140)
      END IF
!
 10   FORMAT (/,' ROMS/TOMS - Output NetCDF summary for Grid ',         &
     &        i2.2,':')
 20   FORMAT (13x,'number of time records written in ',                 &
     &        a,' file = ',i8.8)
 30   FORMAT (/,a,/)
 40   FORMAT (/,a,i3,/)
 50   FORMAT (/,' ROMS/TOMS: DONE... ',a)
 60   FORMAT (/,' MAIN: Abnormal termination: BLOWUP.')
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.',/,        &
     &          ' REASON: ',a)
 80   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.',/,       &
     &          ' REASON: ',a)
 90   FORMAT (/,' ERROR: I/O related problem.')
100   FORMAT (/,' ERROR: Illegal model configuration.')
110   FORMAT (/,' ERROR: Illegal domain partition.')
120   FORMAT (/,' ERROR: Illegal input parameter.')
130   FORMAT (/,' ERROR: Fatal algorithm result.')
140   FORMAT (/,' ERROR: Fatal frazil ice check.')

      RETURN
      END SUBROUTINE close_out
