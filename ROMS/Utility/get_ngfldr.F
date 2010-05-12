#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE get_ngfldr (ng, model, ifield, ncid, nfiles, fname,    &
     &                       update, LBi, UBi, UBj, UBk,                &
     &                       Istr, Iend, Jrec,                          &
     &                       Fout)
!
!svn $Id: get_ngfldr.F 1598 2008-07-25 03:29:09Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested non-grided field from specified     !
!  NetCDF file.  A non-grided field has different dimensions  than     !
!  model spatial dimensions. Backward time processing.                 !
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
      logical :: Lgridded, Liocycle, Lonerec
      logical :: foundit, got_var, got_time

      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ic, ifile, j, lend, lstr, npts
      integer :: nvatt, nvdim, status
      integer :: Vsize(4)

      real(r8) :: Aval, Clength, Fmax, Fmin, Tdelta, Tend
      real(r8) :: Tmax, Tmin, Tmono, Tscale, Tstr, scale
      real(r8) :: Tval

      real(r8), dimension((UBi-LBi+1)*UBj) :: A

      character (len=14) :: t_code
      character (len=45) :: attnam
!
      SourceFile='get_ngfldr.F'
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
        Nrec=0
        Liocycle=.FALSE.
        Lgridded=.FALSE.
        Lonerec=.FALSE.
        got_var=.FALSE.
        got_time=.FALSE.
        Vtype=Iinfo(1,ifield,ng)
!
!  If several input NetCDF files (nfiles>1), scan files until the
!  requested variable is found.
!
        foundit=.FALSE.
        QUERY: DO ifile=1,nfiles
!
!  Inquire about the dimensions and check for consistency.
!
          CALL netcdf_check_dim (ng, model, fname(ifile))
          IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about requested variable.
!
          CALL netcdf_inq_var (ng, model, fname(ifile),                 &
     &                         MyVarName = TRIM(Vname(1,ifield)),       &
     &                         SearchVar = foundit,                     &
     &                         VarID = Vid,                             &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
!
!  Determine if gridded or point data.  Set variable dimensions.
!
          IF (foundit) THEN
            got_var=.TRUE.
            ncfile=fname(ifile)
            Vsize=0
            DO i=1,nvdim
              Vsize(i)=var_Dsize(i)
            END DO
          END IF
!
!  If "scale_factor" attribute is present for a variable, the data are
!  to be multiplied by this factor.  Check if only water points are
!  available.
!
          IF (foundit) THEN
            DO i=1,nvatt
              IF (TRIM(var_Aname(i)).eq.'scale_factor') THEN
                scale=var_Afloat(i)
                Finfo(10,ifield,ng)=scale
              ELSE IF (TRIM(var_Aname(i)).eq.'water_points') THEN
                Iinfo(1,ifield,ng)=-ABS(Iinfo(1,ifield,ng))
                Vtype=Iinfo(1,ifield,ng)
              END IF
            END DO
          END IF
!
!  Inquire about associated time dimension, if any, and get number of
!  available time records.
!
          IF (foundit) THEN
            IF (LEN_TRIM(Vname(5,ifield)).gt.0) THEN
              Tname(ifield)=TRIM(Vname(5,ifield))
              DO i=1,nvdim
                IF (var_Dname(i).eq.TRIM(Tname(ifield))) THEN
                  Nrec=var_Dsize(i)
                  got_time=.TRUE.
                END IF
              END DO
            END IF
!
!  If the associated time variable is different to that specified in
!  "varinfo.dat" (Nrec still zero), reset associated time-variable name
!  to that specified in the "time" attribute.
!
            IF (Nrec.eq.0) THEN
              DO i=1,nvatt
                IF (TRIM(var_Aname(i)).eq.'time') THEN
                  Tname(ifield)=TRIM(var_Achar(i))
                  got_time=.TRUE.
                END IF
              END DO
              IF (got_time) THEN
                DO i=1,n_dim
                  IF (TRIM(dim_name(i)).eq.TRIM(Tname(ifield))) THEN
                    Nrec=dim_size(i)
                  END IF
                END DO
              END IF
            END IF
!
!  If Nrec=0, input file is not CF compliant, check variable dimension
!  to see if the dimension contains the "time" string.
!
            IF (got_time.and.(Nrec.eq.0)) THEN
              DO i=1,n_vdim
                IF (INDEX(TRIM(var_Dname(i)),'time').ne.0) THEN
                  Nrec=var_Dsize(i)
                END IF
              END DO
            END IF
            IF (got_time.and.(Nrec.eq.0)) THEN
              IF (Master) WRITE (stdout,10) TRIM(Tname(ifield)),        &
     &                                      TRIM(Vname(1,ifield)),      &
     &                                      TRIM(fname(ifile))
              exit_flag=4
              RETURN
            END IF
          END IF
          IF ((UBk.eq.1).and.(Nrec.gt.UBi)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,ifield)),        &
     &                                    UBi, Nrec
            exit_flag=4
            RETURN
          END IF
!
!  Determine initial time record to read and cycling switch.
!
          IF (foundit) THEN
            CALL get_cycle (ng, model, ncfile, Tname(ifield), Nrec,     &
     &                      tdays(ng), Tid, Liocycle, Clength,          &
     &                      Trec, Tstr, Tend, Tmin, Tmax,  Tscale)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Store variable information into global information arrays.
!
          IF (foundit) THEN
            Linfo(1,ifield,ng)=Lgridded
            Linfo(2,ifield,ng)=Liocycle
            Iinfo(2,ifield,ng)=Vid
            Iinfo(3,ifield,ng)=Tid
            Iinfo(4,ifield,ng)=Nrec
            Iinfo(5,ifield,ng)=Vsize(1)
            Iinfo(6,ifield,ng)=Vsize(2)
            Iinfo(7,ifield,ng)=Vsize(3)
            Finfo(1,ifield,ng)=Tmin
            Finfo(2,ifield,ng)=Tmax
            Finfo(3,ifield,ng)=Tstr
            Finfo(4,ifield,ng)=Tend
            Finfo(5,ifield,ng)=Clength
            Finfo(6,ifield,ng)=Tscale
            EXIT QUERY
          END IF
        END DO QUERY
!
!  Terminate execution requested variables are not found.
!
        IF (.not.got_var) THEN
          WRITE (stdout,30) TRIM(Vname(1,ifield)), TRIM(ncfile)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_time) THEN
          WRITE (stdout,30) TRIM(Tname(ifield)), TRIM(ncfile)
          exit_flag=2
          RETURN
        END IF
!
!  If appropriate, open input NetCDF file for reading.  If processing
!  model forcing (multiple files allowed), check if file for requested
!  field has been already opened and get/save its ID from/to the
!  association table.
!
        CALL netcdf_openid (nFfiles(ng), FRCname(1,ng), FRCids(1,ng),   &
     &                      ncfile, ncid)
        IF (ncid.eq.-1) THEN
          CALL netcdf_open (ng, model, ncfile, 0, ncid)
          IF (exit_flag.ne.NoError) THEN
            WRITE (stdout,40) TRIM(ncfile)
            RETURN
          END IF
          CALL netcdf_openid (nFfiles(ng), FRCname(1,ng), FRCids(1,ng), &
     &                        ncfile, ncid)
        END IF
        Cinfo(ifield,ng)=TRIM(ncfile)
!
!  The strategy here is to create a local, monotonically decreasing
!  time variable so the interpolation between snapshots is trivial
!  when cycling forcing fields.
!
        IF (Iend.eq.1) THEN
          Tindex=UBk
        ELSE
          Tindex=1
        END IF
        IF (Liocycle) THEN
          IF (Trec.eq.Nrec) THEN
            IF (tdays(ng).lt.Tmax) THEN
              Tmono=Tend
            ELSE
              Tmono=Tend+(tdays(ng)-MOD(tdays(ng),Clength))
            END IF
          ELSE
            IF (tdays(ng).gt.Clength) THEN
              Tmono=Tend+(tdays(ng)-MOD(tdays(ng),Clength))
            ELSE
              Tmono=Tend
            END IF
          END IF
          Trec=Trec+2
        ELSE
          Tmono=Tend
          Trec=Trec+1
          CALL netcdf_get_fvar (ng, model, ncfile, Tname(ifield),       &
     &                          Tval,                                   &
     &                          ncid = ncid,                            &
     &                          start = (/Trec-1/),                     &
     &                          total = (/1/))
          IF (exit_flag.ne.NoError) THEN
            IF (Master) WRITE (stdout,50) TRIM(Tname(ifield)), Trec
            RETURN
          END IF
          Tval=Tval*Tscale
          IF (Tval.lt.Tend) THEN
            Trec=Trec+1
          END IF
        END IF
        Tmono=Tmono*day2sec
        Iinfo(8,ifield,ng)=Tindex
        Iinfo(9,ifield,ng)=Trec
        Finfo(7,ifield,ng)=Tmono
!
!  Set switch for one time record dataset. In this case, the input
!  data is always the same and time interpolation is not performed.
!
        IF (Nrec.eq.1) Lonerec=.TRUE.
        Linfo(3,ifield,ng)=Lonerec
        Tindex=Iinfo(8,ifield,ng)
        Vtime(Tindex,ifield,ng)=Finfo(4,ifield,ng)
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
      IF ((Tmono.gt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart(ng))) THEN
        IF (Liocycle) THEN
          Trec=MOD(Trec,Nrec)-1
          IF (Trec.le.0) Trec=Nrec+Trec
        ELSE
          Trec=Trec-1
        END IF
        Iinfo(9,ifield,ng)=Trec

        IF ((1.le.Trec).and.(Trec.le.Nrec)) THEN
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
          IF ((Tid.ge.0).and.(Tid.ne.Vid)) THEN
            CALL netcdf_get_fvar (ng, model, ncfile, Tname(ifield),     &
     &                            Tval,                                 &
     &                            ncid = ncid,                          &
     &                            start = (/Trec/),                     &
     &                            total = (/1/))
            IF (exit_flag.ne.NoError) THEN
              IF (Master) WRITE (stdout,50) TRIM(Tname(ifield)), Trec
              RETURN
            END IF
          END IF
          IF ((Tid.ge.0).and.(Tid.ne.Vid)) THEN
            Tval=Tval*Tscale
            Vtime(Tindex,ifield,ng)=Tval
          END IF
!
!  Read in non-grided data. The conditional statement on Jrec is to
!  differentiate between reading a 3D and 2D array.
!
          IF (Vid.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF ((Jrec.gt.1).or.                                         &
     &          ((Jrec.eq.1).and.(Iinfo(7,ifield,ng).gt.0))) THEN
              npts=(Iend-Istr+1)*Jrec
              CALL netcdf_get_fvar (ng, model, ncfile, Vname(1,ifield), &
     &                              A,                                  &
     &                              ncid = ncid,                        &
     &                              start = (/1,1,Trec/),               &
     &                              total = (/Iend-Istr+1,Jrec,1/))
            ELSE
              npts=Iend-Istr+1
              CALL netcdf_get_fvar (ng, model, ncfile, Vname(1,ifield), &
     &                              A,                                  &
     &                              ncid = ncid,                        &
     &                              start = (/1,Trec/),                 &
     &                              total = (/Iend-Istr+1,1/))
            END IF
            IF (exit_flag.ne.NoError) THEN
              IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)), Trec
              RETURN
            END IF
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
                WRITE (stdout,60) TRIM(Vname(2,ifield)), Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                Tval=Tval*day2sec
                CALL time_string (Tval, t_code)
                WRITE (stdout,70) TRIM(Vname(2,ifield)), t_code,        &
     &                            Trec, Tindex, ncfile(lstr:lend),      &
     &                            Tmin, Tmax, Fmin, Fmax
              END IF
            END IF
            update=.TRUE.
          END IF
        END IF
!
!  Decrease the local time variable "Tmono" by the interval between
!  snapshots. If the interval is negative, indicating cycling, add in
!  a cycle length.  Load time value (sec) into "Tintrp" which used
!  during interpolation between snapshots.
!
        IF (.not.Lonerec) THEN
          Tdelta=Vtime(3-Tindex,ifield,ng)-Vtime(Tindex,ifield,ng)
          IF (Liocycle.and.(Tdelta.lt.0.0_r8)) THEN
            Tdelta=Tdelta+Clength
          END IF
          Tmono=Tmono-Tdelta*day2sec
          Finfo(7,ifield,ng)=Tmono
          Tintrp(Tindex,ifield,ng)=Tmono
        END IF
      END IF
!
  10  FORMAT (/,' GET_NGFLDR  - unable to find dimension ',a,           &
     &        /,15x,'for variable: ',a,/,15x,'in file: ',a,             &
     &        /,15x,'file is not CF compliant...')
  20  FORMAT (/,' GET_NGFLDR  - too small dimension for variable ',a,   &
     &        ': ',2i5)
  30  FORMAT (/,' GET_NGFLDR  - unable to find requested variable: ',a, &
     &       /,15x,'in file: ',a)
  40  FORMAT (/,' GET_NGFLDR  - unable to open input NetCDF file: ',a)
  50  FORMAT (/,' GET_NGFLDR  - error while reading variable: ',a,2x,   &
     &          ' at TIME index = ',i4)
  60  FORMAT (3x,' GET_NGFLDR  - ',a,/,19x,'(Min = ',1pe15.8,' Max = ', &
     &        1pe15.8,')')
  70  FORMAT (3x,' GET_NGFLDR  - ',a,',',t62,'t = ',a,/,19x,            &
     &        '(Rec=',i4.4,', Index=',i1,', File: ',a,')',/,19x,        &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',/, 19x,            &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE get_ngfldr
#else
      SUBROUTINE get_ngfldr
      RETURN
      END SUBROUTINE get_ngfldr
#endif
