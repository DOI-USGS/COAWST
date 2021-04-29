#include "cppdefs.h"
      SUBROUTINE get_cycle (ng, model, ifield, job, Lmulti,             &
     &                      ncname, ncid, TvarName, ntime, smday,       &
     &                      Tid, Lcycle, clength, Tindex,               &
     &                      Tstr, Tend, Tmin, Tmax, Tscale)
!
!svn $Id: get_cycle.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!     ifield     Field ID.                                             !
!     job        Processing job (integer):                             !
!                  job < 0        backward time logic, adjoint model   !
!                  job > 0        forward time logic                   !
!     Lmulti     Switch to process a multi-file field. That is, the    !
!                  time records are split in several NetCDF files.     !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
!     TvarName   Requested time variable name (string)                 !
!     ntime      Size of time dimension (integer)                      !
!     smday      Starting model day (real)                             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Tid        NetCDF field time variable ID (integer)               !
!     Lcycle     Switch indicating cycling of input fields (logical)   !
!     clength    Length of field time cycle (real)                     !
!     Tindex     Starting field time index plus/minus one to read      !
!                  according job flag (integer)                        !
!                    job > 0       read Tindex+1   (NLM, TLM, RPM)     !
!                    job < 0       read Tindex-1   (ADM)               !
!     Tstr       Data time lower bound enclosing "smday" (real)        !
!     Tend       Data time upper bound enclosing "smday" (real)        !
!     Tmin       Data starting (first record) day (real)               !
!     Tmax       Data ending (last record) day (real)                  !
!     Tscale     Scale to convert time coordinate to day units (real)  !
!                                                                      !
!                                                                      !
! Notice that during initialization (initial, tl_initial, rp_initial,  !
! or ad_initial) this routine process the information for the first    !
! snapshot of the time interpolation. In the forward time stepping     !
! (NLM, TLM, and RPM), the first data snapshot corresponds to the      !
! LOWER interpolant.  Contrarily, in backward time stepping (ADM),     !
! the first data snapshot corresponds to the UPPER interpolant.  At    !
! initialization, there is no information yet if a particular field is !
! from a multifile or not.  So the strategy is to treat such field     !
! differently and resolve it multifile status latter elsewhere when    !
! processing the second snapshot (UPPER interpolant for the forward    !
! and LOWER interpolant for the backward time stepping).               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
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
      logical, intent(in) :: Lmulti
      logical, intent(out) :: Lcycle

      integer, intent(in) :: ng, model, ifield, job, ncid, ntime

      integer, intent(out) :: Tindex, Tid

      real(dp), intent(in) :: smday

      real(dp), intent(out) :: Tmax, Tmin, Tend, Tscale, Tstr, clength

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: TvarName
!
!  Local variable declarations.
!
      logical :: TimeLess
      logical :: Linside, LowerBound, Upperbound
!
      integer :: i, nvatt, nvdim
!
      real(dp) :: mday, tstart
      real(dp) :: Tval(ntime)
!
      character (len=40) :: tunits

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Find time cycling parameters, if any.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      Lcycle=.FALSE.
      Linside=.FALSE.
      LowerBound=.FALSE.
      TimeLess=.FALSE.
      UpperBound=.FALSE.
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
     &                       ncid = ncid,                               &
     &                       MyVarName = TRIM(TvarName),                &
     &                       VarID = Tid,                               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Check if time cycling attribute is present and then read in time
!  cycle length.  Check time coordinate units and determine time
!  scale.  The internal processing of all fields requires time in
!  day units.  Check if more than one time record is available.
!
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'cycle_length') THEN
            Lcycle=.TRUE.
            IF (var_Afloat(i).gt.0.0_r8) THEN
              clength=var_Afloat(i)
            ELSE IF (var_Aint(i).gt.0) THEN    ! no CF compliance
              clength=REAL(var_Aint(i),r8)     ! attribute is an integer
            ELSE
              IF (Master) WRITE (stdout,10) TRIM(var_Aname(i)),         &
     &                                      TRIM(TvarName)
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
      CALL netcdf_get_time (ng, model, ncname, TvarName,                &
     &                      Rclock%DateNumber, Tval,                    &
     &                      ncid = ncid,                                &
     &                      start = (/1/),                              &
     &                      total = (/ntime/),                          &
     &                      min_val = Tmin,                             &
     &                      max_val = Tmax)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Scale time variable to days.  Determine the minimum and maximum time
!  values available
!
      DO i=1,ntime
        Tval(i)=Tval(i)*Tscale
      END DO
      Tmin=Tmin*Tscale
      Tmax=Tmax*Tscale
      IF (Lcycle) THEN
        mday=MOD(smday,clength)
      ELSE
        mday=smday
      END IF
!
!  Is the model time inside the data time range? If not, check if the
!  data just has the LOWER- or the UPPER-snapshot interpolant.
!
      IF ((Tmin.le.mday).and.(mday.le.Tmax)) THEN
        Linside=.TRUE.
      ELSE IF (mday.ge.Tmax) THEN
        LowerBound=.TRUE.
      ELSE IF (mday.le.Tmin) THEN
        UpperBound=.TRUE.
      END IF
!
!  If processing timeless variable, set processing parameter.
!
      IF (TimeLess) THEN
        Tstr=Tmin
        Tend=Tmax
        Tindex=ntime
        RETURN
      END IF
!
!  If processing field split in several files, find UPPER time-snapshot
!  and its associated time record (Tindex).
!
      IF (Lmulti) THEN
        IF (job.gt.0) THEN          ! forward  time-stepping logic
          DO i=1,ntime
            IF (Tval(i).gt.mday) THEN
              Tindex=i-1
              Tend=Tval(i)
              EXIT
            END IF
          END DO
        ELSE                        ! backward time-stepping logic
          DO i=ntime,1,-1
            IF (Tval(i).le.mday) THEN
              Tindex=i+1
              Tstr=Tval(i)
              EXIT
            END IF
          END DO
        END IF
!
!  If not processing a multi-file field, find LOWER time-snapshot
!  and its associated time record (Tindex).
!
      ELSE
        IF (job.gt.0) THEN     ! forward:   Tval(i) =< mday =< Tval(i+1)
          IF (Linside) THEN
            tstart=Tmin
            IF (ntime.eq.1) THEN
              Tstr=Tmin
              Tindex=1
            ELSE
              DO i=2,ntime
                IF ((tstart.le.mday).and.(mday.le.Tval(i))) THEN
                  Tstr=tstart
                  Tindex=i-1   ! one is added when processing
                  EXIT
                END IF
                tstart=Tval(i)
              END DO
            END IF
          ELSE
            Tstr=Tmax          ! LowerBound for next multifile or
            Tindex=ntime       ! time cycling
          END IF
        ELSE                   ! backward:  Tval(i) =< mday =< Tval(i+1)
          IF (Linside) THEN
            tstart=Tmin
            IF (ntime.eq.1) THEN
              Tstr=Tmin
              Tindex=1
            ELSE
              DO i=2,ntime
                IF ((tstart.le.mday).and.(mday.le.Tval(i))) THEN
                  Tstr=Tval(i)
                  Tindex=i
                  EXIT
                END IF
                tstart=Tval(i)
              END DO
            END IF
          ELSE IF (LowerBound) THEN
            Tstr=Tmax
            Tindex=ntime
          ELSE IF (UpperBound) THEN
            Tstr=Tmin
            Tindex=1
          END IF
        END IF
      END IF
!
!  If processing a multi-file field, set LOWER time-snapshot. It
!  is the last value from previous file. Otherwise, set UPPER
!  time-snapshot.
!
      IF (Lmulti) THEN
        IF (job.gt.0) THEN
          Tstr=Finfo(2,ifield,ng)      ! Tmax from previous file
        ELSE
          Tend=Finfo(1,ifield,ng)      ! Tmin from previous file
        END IF
      ELSE
        IF (Lcycle.and.(Tindex.eq.ntime)) THEN
          Tend=Tmin
        ELSE
          IF (job.gt.0) THEN
            i=MIN(ntime,Tindex+1)
            Tend=Tval(i)
          ELSE
            i=MAX(1,Tindex-1)
            Tend=Tval(i)
          END IF
        END IF
      END IF
!
!  If not cycling, stop execution if there is not field data
!  available for current model time. Avoid check on tidal data
!  since time is in terms of frequencies.
!
      IF (.not.TimeLess) THEN
        IF (.not.Lcycle.and.(ntime.gt.1)) THEN
          IF (Lmulti) THEN
            IF (job.gt.0) THEN
              IF (smday.gt.Tmax) THEN
                IF (Master) WRITE (stdout,20) TRIM(TvarName),           &
     &                                        Tmax, smday
                exit_flag=2
                RETURN
              ELSE IF (Linfo(6,ifield,ng)) THEN
                IF (Tmin.lt.smday) THEN
                  IF (Master) THEN
                    WRITE (stdout,30)                                   &
     &                  'Upper snapshot time for multi-file variable:', &
     &                                TRIM(TvarName),                   &
     &                                TRIM(Vname(1,ifield)),            &
     &                  'is less than current model time.',             &
     &                                'Tmin = ', Tmin, smday
                  END IF
                  exit_flag=2
                  RETURN
                END IF
              END IF
            ELSE
              IF (Linfo(5,ifield,ng)) THEN
                IF (Tmin.gt.smday) THEN
                  IF (Master) THEN
                    WRITE (stdout,30)                                   &
     &                  'Lower snapshot time for multi-file variable:', &
     &                                TRIM(TvarName),                   &
     &                                TRIM(Vname(1,ifield)),            &
     &                  'is greater than current model time.',          &
     &                                'Tmin = ', Tmin, smday
                  END IF
                  exit_flag=2
                  RETURN
                END IF
              ELSE
                IF (smday.lt.Tmax) THEN
                  IF (Master) THEN
                    WRITE (stdout,30)                                   &
     &                    'starting time for multi-file variable:',     &
     &                                TRIM(TvarName),                   &
     &                                TRIM(Vname(1,ifield)),            &
     &                    'is greater than current model time.',        &
     &                                'Tmax = ', Tmax, smday
                  END IF
                  exit_flag=2
                  RETURN
                END IF
              END IF
            END IF
          ELSE
            IF (.not.UpperBound.and.(smday.lt.Tmin)) THEN
              IF (Master) WRITE (stdout,40) TRIM(TvarName),             &
     &                                      TRIM(Vname(1,ifield)),      &
     &                                      Tmin, smday
              exit_flag=2
              RETURN
            END IF
          END IF
        END IF
      END IF
!
  10  FORMAT (/,' GET_CYCLE - unable to get value for attribute: ',a,   &
     &        /,13x,'in variable: ',a,                                  &
     &        /,13x,'This attribute value is expected to be of',        &
     &        /,13x,'the same external type as the variable.')
  20  FORMAT (/,' GET_CYCLE - ending time for multi-file variable: ',a, &
     &        /,13x,'is less than current model time. ',                &
     &        /,13x,'TMAX = ',f15.4,2x,'TDAYS = ',f15.4)
  30  FORMAT (/,' GET_CYCLE - ',a,1x,a,2x,'(',a,')',/,13x,a,            &
     &        /,13x,a,f15.4,2x,'TDAYS = ',f15.4)
  40  FORMAT (/,' GET_CYCLE - starting time for variable: ',a,2x,       &
     &        '(',a,')',/,13x,'is greater than current model time.',    &
     &        /,13x,'TMIN = ',f15.4,2x,'TDAYS = ',f15.4)

      RETURN
      END SUBROUTINE get_cycle
