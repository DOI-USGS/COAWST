#include "cppdefs.h"
      SUBROUTINE get_ngfld (ng, model, ifield, ncid, nfiles, fname,     &
     &                      update, LBi, UBi, UBj, UBk,                 &
     &                      Istr, Iend, Jrec,                           &
     &                      Fout)
!
!svn $Id: get_ngfld.F 1598 2008-07-25 03:29:09Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested non-grided field from specified     !
!  NetCDF file.  A non-grided field has different dimensions  than     !
!  model spatial dimensions. Forward time processing.                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ifield     Field ID.                                             !
!     ncid       NetCDF file ID.                                       !
!     nfiles     Number of input NetCDF files.                         !
!     fname      NetCDF file name(s).                                  !
!     LBi        "Fout" 1st dimension lower-bound value.               !
!     UBi        "Fout" 1st dimension upper-bound value.               !
!     UBj        "Fout" 2nd dimension upper-bound value, if any.       !
!                  Otherwise, a value of one is expected.              !
!     UBk        "Fout" time dimension upper-bound value, if any.      !
!                  Otherwise, a value of one is expected.              !
!     Istr       Starting location of read data in the 1st dimension.  !
!     Iend       Ending location of read data in the 1st dimension;    !
!                  Number of records read is: Iend-Istr+1.             !
!     Jrec       Number of records read in the 2st dimension.          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Read field.                                           !
!     update     Switch indicating reading of the requested field      !
!                  the current time step.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
      USE distribute_mod, ONLY : mp_bcastl, mp_bcasts
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, model, ifield, nfiles
      integer, intent(in) :: LBi, UBi, UBj, UBk, Istr, Jrec
      integer, intent(inout) :: ncid, Iend

      character (len=*), intent(in) :: fname(nfiles)

      real(r8), intent(inout) :: Fout(LBi:UBi,UBj,UBk)
!
!  Local variable declarations.
!
      logical :: Lgrided, Liocycle, Lonerec, got_var, got_time

      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ic, j, lend, lfvar, lstr, ltvar, npts, nvatts
      integer :: nvd, status, vartype

      integer :: start(3), total(3), Vsize(4)

      real(r8) :: Aval, Clength, Fmax, Fmin, Tdelta, Tend
      real(r8) :: Tmax, Tmin, Tmono, Tscale, Tstr, scale
      real(r8) :: Tval(1)

      real(r8), dimension((UBi-LBi+1)*UBj) :: A

      character (len=14) :: t_code
      character (len=45) :: attnam
!
!-----------------------------------------------------------------------
!  On first call, inquire about the contents of input NetCDF file.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
      IF (iic(ng).eq.0) THEN
!
!  Intialize local variables.
!
        Vid=-1
        Tid=-1        
        Liocycle=.FALSE.
        Lgrided=.FALSE.
        Lonerec=.FALSE.
        got_var=.FALSE.
        got_time=.FALSE.
        Vtype=Iinfo(1,ifield,ng)
        lfvar=LEN_TRIM(Vname(1,ifield))
        ltvar=LEN_TRIM(Tname(ifield))
!
!  Inquire about the dimensions and variables. Check for consistency.
!
        IF (InpThread) THEN
          CALL opencdf (ng, nfiles, fname, ncfile, N(ng), ifield, Nrec, &
     &                  nvd, Vsize)
          IF (exit_flag.ne.NoError) RETURN
          IF ((UBk.eq.1).and.(Nrec.gt.UBi)) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), UBi, Nrec
            exit_flag=4
            RETURN
          END IF
        END IF
        IF (Iend.eq.0) THEN
          Iend=Nrec
#ifdef DISTRIBUTE
          CALL mp_bcasti (ng, model, Iend, 1)
#endif
        END IF
!
!  Scan variable list from input NetCDF and check for requested
!  variables.
!
        IF (InpThread) THEN
          DO i=1,nvars
            IF (lfvar.gt.0) THEN
              IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
                Vid=i
                got_var=.TRUE.
              END IF
            END IF
            IF (ltvar.gt.0) THEN
              IF (TRIM(varnam(i)).eq.TRIM(Tname(ifield))) THEN
                Tid=i
                got_time=.TRUE.
              END IF
            END IF
          END DO
          Linfo(1,ifield,ng)=Lgrided
          Iinfo(2,ifield,ng)=Vid
          Iinfo(3,ifield,ng)=Tid
          Iinfo(4,ifield,ng)=Nrec
          Iinfo(5,ifield,ng)=Vsize(1)
          Iinfo(6,ifield,ng)=Vsize(2)
          Iinfo(7,ifield,ng)=Vsize(3)
!
!  Terminate execution requested variables are not found.
!
          IF (.not.got_var) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)), TRIM(ncfile)
            exit_flag=2
            RETURN
          END IF
          IF (.not.got_time) THEN
            WRITE (stdout,20) TRIM(Tname(ifield)), TRIM(ncfile)
            exit_flag=2
            RETURN
          END IF
!
!  If appropriate, open input NetCDF file for reading.  If processing
!  model forcing (multiple files allowed), check if file for requested
!  field has been already opened and get/save its ID from/to the
!  association table. 
!
          CALL openids (nFfiles(ng), FRCname(1,ng), FRCids(1,ng),       &
     &                  ncfile, ncid)
          IF (ncid.eq.-1) THEN
            status=nf90_open(TRIM(ncfile), nf90_nowrite, ncid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(ncfile)
              exit_flag=2
              ioerror=status
              RETURN
            END IF
            CALL openids (nFfiles(ng), FRCname(1,ng), FRCids(1,ng),     &
     &                    ncfile, ncid)
          END IF
          Cinfo(ifield,ng)=TRIM(ncfile)
!
!  If "scale_factor" attribute is present for a variable, the data are
!  to be multiplied by this factor.  Check if only water points are
!  available.
!
          IF (got_var) THEN
            status=nf90_inquire_variable(ncid, Vid, varnam(ifield),     &
     &                                   vartype, nvdims(ifield),       &
     &                                   vdims(:,ifield), nvatts)
            IF (status.eq.nf90_noerr) THEN
              DO i=1,nvatts
                status=nf90_inq_attname(ncid, Vid, i, attnam)
                IF (status.eq.nf90_noerr) THEN
                  IF (TRIM(attnam).eq.'scale_factor') THEN
                    status=nf90_get_att(ncid, Vid, TRIM(attnam),        &
     &                                  scale)
                    IF (status.eq.nf90_noerr) THEN
                      Fscale(ifield,ng)=Fscale(ifield,ng)*scale
                    ELSE
                      WRITE (stdout,30) TRIM(attnam)
                      exit_flag=2
                      RETURN
                    END IF
                  ELSE IF (TRIM(attnam).eq.'water_points') THEN
                    Iinfo(1,ifield,ng)=-ABS(Iinfo(1,ifield,ng))
                    Vtype=Iinfo(1,ifield,ng)
                  END IF
                ELSE
                  WRITE (stdout,40) TRIM(varnam(ifield))
                  exit_flag=2
                  RETURN
                END IF
              END DO
            ELSE
              WRITE (stdout,50)
              exit_flag=2
              RETURN
            END IF
          END IF
!
!  Determine initial time record to read and cycling switch.
!
          IF (UBk.eq.1) THEN
            Tstr=1
            Tend=1
            Trec=1
            Clength=0.0_r8
            Tscale=1.0_r8
          ELSE
            CALL get_cycle (ng, ncid, Tid, Nrec, tdays(ng), Liocycle,   &
     &                      Clength, Trec, Tstr, Tend, Tmin, Tmax,      &
     &                      Tscale)
            IF (exit_flag.ne.NoError) RETURN
          END IF
          Linfo(2,ifield,ng)=Liocycle
          Finfo(1,ifield,ng)=Tmin
          Finfo(2,ifield,ng)=Tmax
          Finfo(3,ifield,ng)=Tstr
          Finfo(4,ifield,ng)=Tend
          Finfo(5,ifield,ng)=Clength
          Finfo(6,ifield,ng)=Tscale
!
!  The strategy here is to create a local, monotonically increasing
!  time variable so the interpolation between snapshots is trivial
!  when cycling forcing fields. Subtract one to time record counter
!  "Trec" to avoid doing special case at initialization.
!
          IF (Iend.eq.1) THEN
            Tindex=UBk
          ELSE
            Tindex=1
          END IF
          IF (Liocycle) THEN
            IF (Trec.eq.Nrec) THEN
              IF (tdays(ng).lt.Tmax) THEN
                Tmono=Tstr-Clength
              ELSE
                Tmono=tdays(ng)+(Tstr-Clength)
                IF (Tstr.eq.Tmax) THEN
                  Tmono=Tmono+(Tmin-MOD(tdays(ng)+Tmin,Clength))
                ELSE
                  Tmono=Tmono+(Tstr-MOD(tdays(ng)+Tstr,Clength))
                END IF
              END IF
            ELSE
              IF (tdays(ng).gt.Clength) THEN
                Tmono=tdays(ng)-MOD(tdays(ng)-Tstr,Clength)
              ELSE
                Tmono=Tstr
              END IF
            END IF
          ELSE
            Tmono=Tstr
          END IF
          Tmono=Tmono*day2sec
          Trec=Trec-1
          Iinfo(8,ifield,ng)=Tindex
          Iinfo(9,ifield,ng)=Trec
          Finfo(7,ifield,ng)=Tmono
!
!  Set switch for one time record dataset. In this case, the input
!  data is always the same and time interpolation is not performed.
!
          IF (Nrec.eq.1) Lonerec=.TRUE.
          Linfo(3,ifield,ng)=Lonerec
        END IF
#ifdef DISTRIBUTE
!
!  Broadcast variable information to all nodes in the group.
!
        CALL mp_bcastl (ng, model, Linfo(:,ifield,ng), 3)
        CALL mp_bcasti (ng, model, Iinfo(:,ifield,ng), 9)
        CALL mp_bcastf (ng, model, Finfo(:,ifield,ng), 7)
        CALL mp_bcasts (ng, model, Cinfo(ifield,ng), 80)
!
#endif
        Tindex=Iinfo(8,ifield,ng)
        Vtime(Tindex,ifield,ng)=Finfo(3,ifield,ng)
      END IF         
!
!-----------------------------------------------------------------------
!  Get requested field information from global storage.
!-----------------------------------------------------------------------
!
      Liocycle=Linfo(2,ifield,ng)
      Lonerec =Linfo(3,ifield,ng)
      Vtype   =Iinfo(1,ifield,ng)
      Vid     =Iinfo(2,ifield,ng)
      Tid     =Iinfo(3,ifield,ng)
      Nrec    =Iinfo(4,ifield,ng)
      Tindex  =Iinfo(8,ifield,ng)
      Trec    =Iinfo(9,ifield,ng)
      Tmin    =Finfo(1,ifield,ng)
      Tmax    =Finfo(2,ifield,ng)
      Clength =Finfo(5,ifield,ng)
      Tscale  =Finfo(6,ifield,ng)
      Tmono   =Finfo(7,ifield,ng)
      ncfile  =Cinfo(ifield,ng)
!
!-----------------------------------------------------------------------
!  If appropriate, read in new data.
!-----------------------------------------------------------------------
!
      update=.FALSE.
      IF ((Tmono.lt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart(ng))) THEN
        IF (Liocycle) THEN
          Trec=MOD(Trec,Nrec)+1
        ELSE
          Trec=Trec+1
        END IF
        Iinfo(9,ifield,ng)=Trec
        IF (Trec.le.Nrec) THEN
!
!  Set rolling index for two-time record storage of input data.  If
!  "UBk" is unity, input data is stored in recordless array by the
!  calling program.
!
         IF (UBk.eq.1) THEN
           Tindex=1
         ELSE
           Tindex=3-Tindex
         END IF
         Iinfo(8,ifield,ng)=Tindex
!
!  Read in time coordinate.
!
          IF (InpThread.and.((Tid.ge.0).and.(Tid.ne.Vid))) THEN
            start(1)=Trec
            total(1)=1
            status=nf90_get_var(ncid, Tid, Tval, start, total)
            IF (InpThread.and.(status.ne.nf90_noerr)) THEN
              WRITE (stdout,40) TRIM(Tname(ifield)), Trec
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
          IF ((Tid.ge.0).and.(Tid.ne.Vid)) THEN
#ifdef DISTRIBUTE
            CALL mp_bcasti (ng, model, exit_flag, 1)
            CALL mp_bcastf (ng, model, Tval, 1)
#endif
            Tval(1)=Tval(1)*Tscale
            Vtime(Tindex,ifield,ng)=Tval(1)
          END IF
!
!  Read in non-grided data.
!
          IF (Vid.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF (Jrec.gt.1) THEN
              start(1)=1
              total(1)=(Iend-Istr)+1
              start(2)=1
              total(2)=Jrec
              start(3)=Trec
              total(3)=1
              npts=total(1)*total(2)*total(3)
            ELSE
              start(1)=1
              total(1)=(Iend-Istr)+1
              start(2)=Trec
              total(2)=1
              npts=total(1)*total(2)
            END IF
            IF (InpThread) THEN
              status=nf90_get_var(ncid, Vid, A, start, total)
              IF (status.ne.nf90_noerr) THEN
                WRITE (stdout,40) TRIM(Vname(1,ifield)), Trec
                exit_flag=2
                ioerror=status
                RETURN
              END IF
            END IF
#ifdef DISTRIBUTE
            CALL mp_bcasti (ng, model, exit_flag, 1)
            CALL mp_bcastf (ng, model, A, npts)
#endif
            Fmin=A(1)*Fscale(ifield,ng)
            Fmax=A(1)*Fscale(ifield,ng)
            ic=0
            DO j=1,Jrec
              DO i=Istr,Iend
                ic=ic+1
                Aval=A(ic)*Fscale(ifield,ng)
                Fmin=MIN(Fmin,Aval)
                Fmax=MAX(Fmax,Aval)
                Fout(i,j,Tindex)=Aval
              END DO
            END DO
            Finfo(8,ifield,ng)=Fmin
            Finfo(9,ifield,ng)=Fmax
            IF (Master) THEN
              IF (UBk.eq.1) THEN
                WRITE (stdout,50) TRIM(Vname(2,ifield)), Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                Tval(1)=Tval(1)*day2sec
                CALL time_string (Tval(1), t_code)
                WRITE (stdout,60) TRIM(Vname(2,ifield)), t_code,        &
     &                            ncfile(lstr:lend), Trec, Tindex,      &
     &                            Tmin, Tmax, Fmin, Fmax
              END IF
            END IF
            update=.TRUE.
          END IF
        END IF
!
!  Increment the local time variable "Tmono" by the interval between
!  snapshots. If the interval is negative, indicating cycling, add in
!  a cycle length.  Load time value (sec) into "Tintrp" which used
!  during interpolation between snapshots.
!
        IF (.not.Lonerec) THEN
          Tdelta=Vtime(Tindex,ifield,ng)-Vtime(3-Tindex,ifield,ng)
          IF (Liocycle.and.(Tdelta.lt.0.0_r8)) THEN
            Tdelta=Tdelta+Clength
          END IF
          Tmono=Tmono+Tdelta*day2sec
          Finfo(7,ifield,ng)=Tmono
          Tintrp(Tindex,ifield,ng)=Tmono
        END IF
      END IF
!
  10  FORMAT (/,' GET_NGFLD   - too small dimension for variable ',a,   &
     &        ': ',2i5)
  20  FORMAT (/,' GET_NGFLD   - unable to find requested variable: ',a, &
     &       /,15x,'in input NetCDF file: ',a)
  30  FORMAT (/,' GET_NGFLD   - unable to open input NetCDF file: ',a)
  40  FORMAT (/,' GET_NGFLD   - error while reading variable: ',a,2x,   &
     &          ' at TIME index = ',i4)
  50  FORMAT (3x,' GET_NGFLD   - ',a,/,19x,'(Min = ',1pe15.8,' Max = ', &
     &        1pe15.8,')')
  60  FORMAT (3x,' GET_NGFLD   - ',a,',',t62,'t = ',a,/,19x,            &
     &        '(File: ',a,', Rec=',i4.4,', Index=',i1,')',/,19x,        &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',/, 19x,            &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
      RETURN
      END SUBROUTINE get_ngfld
