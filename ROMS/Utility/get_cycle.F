#include "cppdefs.h"
      SUBROUTINE get_cycle (ng, ncid, tfldid, ntime, smday, cycle,      &
     &                      clength, Tindex, Tstr, Tend, Tmin, Tmax,    &
     &                      Tscale)
!
!svn $Id: get_cycle.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine determines relevant parameters for time cycling        !
!  of data from a input NetCDF file.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number.                                     !
!     ncid     NetCDF file ID.                                         !
!     tfldid   NetCDF field time variable ID.                          !
!     ntime    Size of time dimension.                                 !
!     smday    Starting model day.                                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     cycle    Logical switch indicating cycling of input fields.      !
!     clength  Length of field time cycle.                             !
!     Tindex   Starting field time index to read.                      !
!     Tstr     Data time lower bound enclosing "smday".                !
!     Tend     Data time upper bound enclosing "smday".                !
!     Tmin     Data starting (first record) day.                       !
!     Tmax     Data ending (last record) day.                          !
!     Tscale   Scale to convert time coordinate to day units.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: cycle

      integer, intent(in) :: ng, ncid, tfldid, ntime

      integer, intent(out) :: Tindex

      real(r8), intent(in) :: smday

      real(r8), intent(out) :: Tmax, Tmin, Tend, Tscale, Tstr, clength
!
!  Local variable declarations.
!
      integer :: i, nvatts, nvdims, status, vartype

      integer, dimension(5) :: start, total, vdims

      real(r8) :: mday, tstart
      real(r8) :: Tval(1)

      character (len=15) :: attnam, varnam
      character (len=40) :: tunits
!
!-----------------------------------------------------------------------
!  Find time cycling parameters, if any.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      cycle=.FALSE.
      Tindex=0
      clength=0.0_r8
      Tstr=0.0_r8
      Tscale=1.0_r8
!
!  Check if time cycling attribute is present and then read in time
!  cycle length.  Check time coordinate units and determine time
!  scale.  The internal processing of all fields requires time in
!  day units.  Check if more than one time record is available.
!
      IF (ntime.ge.1) THEN
        status=nf90_inquire_variable(ncid, tfldid, varnam, vartype,       &
     &                               nvdims, vdims, nvatts)
        IF (status.eq.nf90_noerr) THEN
          DO i=1,nvatts
            status=nf90_inq_attname(ncid, tfldid, i, attnam)
            IF (status.eq.nf90_noerr) THEN
              IF (TRIM(attnam).eq.'cycle_length') THEN
                cycle=.TRUE.
                status=nf90_get_att(ncid, tfldid, TRIM(attnam),           &
     &                              clength)
                IF (status.ne.nf90_noerr) THEN
                  WRITE (stdout,10) TRIM(attnam)
                  exit_flag=2
                  ioerror=status
                  RETURN
                END IF
              ELSE IF (TRIM(attnam).eq.'units') THEN
                status=nf90_get_att(ncid, tfldid, TRIM(attnam),           &
     &                              tunits)
                IF (status.eq.nf90_noerr) THEN
                  IF (tunits(1:6).eq.'second') THEN
                    Tscale=sec2day
                  END IF
                ELSE
                  WRITE (stdout,10) TRIM(attnam)
                  exit_flag=2
                  RETURN
                END IF
              END IF
            ELSE
              WRITE (stdout,20) TRIM(varnam)
              exit_flag=2
              RETURN
            END IF
          END DO
        ELSE
          WRITE (stdout,30)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Read in minimum (Tmin) and maximum (Tmax) times available for input
!  field.
!
      start(1)=1
      total(1)=1
      status=nf90_get_var(ncid, tfldid, Tval, start, total)
      Tmin=Tval(1)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,40) TRIM(varnam)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      start(1)=ntime
      total(1)=1
      status=nf90_get_var(ncid, tfldid, Tval, start, total)
      Tmax=Tval(1)
      IF (status.ne.nf90_noerr) THEN
        WRITE (stdout,40) TRIM(varnam)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
      Tmin=Tmin*Tscale
      Tmax=Tmax*Tscale
!
!  Find lower time-snapshot (Tstr) and its associated record (Tindex).
!
      IF (cycle) THEN
        mday=MOD(smday,clength)
      ELSE
        mday=smday
      END IF
      IF ((mday.lt.Tmin).or.(mday.ge.Tmax)) THEN
        Tindex=ntime
        Tstr=Tmax
      ELSE
        tstart=Tmin
        DO i=2,ntime
          start(1)=i
          total(1)=1
          status=nf90_get_var(ncid, tfldid, Tval, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(varnam)
            exit_flag=2
            ioerror=status
            RETURN
          END IF
          Tval(1)=Tval(1)*Tscale
          IF ((tstart.le.mday).and.(mday.lt.Tval(1))) THEN
            Tindex=i-1
            Tstr=tstart
            GO TO 50
          END IF
          tstart=Tval(1)
        END DO
  50    CONTINUE
      END IF
!
!  Find upper time-snapshot (Tend).
!
      IF (cycle.and.(Tindex.eq.ntime)) THEN
        Tend=Tmin
      ELSE
        start(1)=MIN(ntime,Tindex+1)
        total(1)=1
        status=nf90_get_var(ncid, tfldid, Tval, start, total)
        Tend=Tval(1)
        IF (status.ne.nf90_noerr) THEN
          write(stdout,40) TRIM(varnam)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        Tend=Tend*Tscale
      END IF
!
!  If not cycling, stop execution if there is not field data
!  available for current model time.
!
      IF (.not.cycle.and.(ntime.gt.1)) THEN
        IF (smday.lt.Tmin) THEN
          WRITE (stdout,60) TRIM(varnam), Tmin, smday
          exit_flag=2
          RETURN
        END IF
      END IF
!
  10  FORMAT (/,' GET_CYCLE - error while reading attribute: ',a)
  20  FORMAT (/,' GET_CYCLE - error while inquiring attributes for',    &
     &          ' variable: ',a)
  30  FORMAT (/,' GET_CYCLE - cannot inquire about time variable from', &
     &          ' input NetCDF file.')
  40  FORMAT (/,' GET_CYCLE - error while reading variable: ',a)
  60  FORMAT (/,' GET_CYCLE - starting time for variable: ',a,          &
     &        /,13x,'is greater than current model time. ',             &
     &        /,13x,'TMIN = ',f15.4,2x,'TDAYS = ',f15.4)
      RETURN
      END SUBROUTINE get_cycle
