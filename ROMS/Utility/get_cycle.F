#include "cppdefs.h"
      SUBROUTINE get_cycle (ng, model, ncname, TvarName, ntime,         &
     &                      smday, Tid, cycle, clength, Tindex,         &
     &                      Tstr, Tend, Tmin, Tmax, Tscale)
!
!svn $Id: get_cycle.F 1336 2008-01-24 02:45:56Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine determines relevant parameters for time cycling        !
!  of data from a input NetCDF file.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     TvarName   Requested time variable name (string)                 !
!     ntime      Size of time dimension (integer)                      !
!     smday      Starting model day (real)                             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Tid        NetCDF field time variable ID (integer)               !
!     cycle      Switch indicating cycling of input fields (logical)   !
!     clength    Length of field time cycle (real)                     !
!     Tindex     Starting field time index to read (integer)           !
!     Tstr       Data time lower bound enclosing "smday" (real)        !
!     Tend       Data time upper bound enclosing "smday" (real)        !
!     Tmin       Data starting (first record) day (real)               !
!     Tmax       Data ending (last record) day (real)                  !
!     Tscale     Scale to convert time coordinate to day units (real)  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: cycle

      integer, intent(in) :: ng, model, ntime

      integer, intent(out) :: Tindex, Tid

      real(r8), intent(in) :: smday

      real(r8), intent(out) :: Tmax, Tmin, Tend, Tscale, Tstr, clength

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: TvarName
!
!  Local variable declarations.
!
      logical :: TimeLess

      integer :: i, nvatt, nvdim

      real(r8) :: mday, tstart
      real(r8) :: Tval(ntime)

      character (len=40) :: tunits
!
      SourceFile='get_cycle.F'
!
!-----------------------------------------------------------------------
!  Find time cycling parameters, if any.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      cycle=.FALSE.
      TimeLess=.FALSE.
      Tindex=0
      clength=0.0_r8
      Tstr=0.0_r8
      Tscale=1.0_r8
!
!  Determine if the associated time variable is actually a time
!  variable.
!
      IF ((INDEX(TRIM(TvarName),'period').gt.0).or.                     &
     &    (TRIM(TvarName).eq.'river')) THEN
        TimeLess=.TRUE.
      END IF
!
!  Inquire input NetCDF file about the time variable.
!
      IF (ntime.ge.1) THEN
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       MyVarName = TRIM(TvarName),                &
     &                       VarID = Tid,                               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (exit_flag.ne.NoError) RETURN
!
!  Check if time cycling attribute is present and then read in time
!  cycle length.  Check time coordinate units and determine time
!  scale.  The internal processing of all fields requires time in
!  day units.  Check if more than one time record is available.
!
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'cycle_length') THEN
            cycle=.TRUE.
            IF (var_Afloat(i).gt.0.0_r8) THEN
              clength=var_Afloat(i)
            ELSE IF (var_Aint(i).gt.0) THEN    ! no CF compliance
              clength=REAL(var_Aint(i),r8)     ! attribute is an integer
            ELSE
              IF (Master) WRITE (stdout,10) TRIM(TvarName)
              exit_flag=2
              RETURN
            END IF
          ELSE IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:6).eq.'second') THEN
              Tscale=sec2day
            END IF
          END IF
        END DO
      END IF
!
!  Read in time variable values.
!
      CALL netcdf_get_fvar (ng, model, ncname, TvarName,                &
     &                      Tval,                                       &
     &                      start = (/1/),                              &
     &                      total = (/ntime/),                          &
     &                      min_val = Tmin,                             &
     &                      max_val = Tmax)
!
!  Scale time variable to days.  Determine the minimum and maximun time
!  values available
!
      DO i=1,ntime
        Tval(i)=Tval(i)*Tscale
      END DO
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
          IF ((tstart.le.mday).and.(mday.lt.Tval(i))) THEN
            Tindex=i-1
            Tstr=tstart
            EXIT
          END IF
          tstart=Tval(i)
        END DO
      END IF
!
!  Find upper time-snapshot (Tend).
!
      IF (cycle.and.(Tindex.eq.ntime)) THEN
        Tend=Tmin
      ELSE
        i=MIN(ntime,Tindex+1)
        Tend=Tval(i)
      END IF
!
!  If not cycling, stop execution if there is not field data
!  available for current model time. Avoid check on tidal data
!  since time is in terms of frequencies.
!
      IF (.not.TimeLess) THEN
        IF (.not.cycle.and.(ntime.gt.1)) THEN
          IF (smday.lt.Tmin) THEN
            IF (Master) WRITE (stdout,20) TRIM(TvarName), Tmin, smday
            exit_flag=2
            RETURN
          END IF
        END IF
      END IF
!
  10  FORMAT (/,' GET_CYCLE - unable to get value for attribute: ',a,   &
     &        /,13x,'in variable: ',a,                                  &
     &        /,13x,'This attribute value is expected to be of',        &
     &        /,13x,'the same external type as the variable.')
  20  FORMAT (/,' GET_CYCLE - starting time for variable: ',a,          &
     &        /,13x,'is greater than current model time. ',             &
     &        /,13x,'TMIN = ',f15.4,2x,'TDAYS = ',f15.4)

      RETURN
      END SUBROUTINE get_cycle
