#include "cppdefs.h"
      SUBROUTINE inquire (ng, model, job, Iout, Irec, Mrec, ifield,     &
     &                    ncid, Lmulti, nfiles, S)
!
!svn $Id: inquire.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine inquires about the contents of input data from NetCDF  !
!  files.                                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     job        Processing job (adjoint model tasks are negative):    !
!                  job = -1, 1    process non-grided field             !
!                  job = -2, 2    process 2D field                     !
!                  job = -3, 3    process 3D field                     !
!     Iout       Size of the outer dimension. Usually, when processing !
!                  data between time snapshots Iout=2. Otherwise,      !
!                  Iout=1 in the calling program.                      !
!     Irec       Number of field records to process. Usually, Irec=1   !
!                  except when processing special fields.              !
!     Mrec       Size of the record dimension when processing          !
!                  non-gridded fields, ABS(job)=1.                     !
!     ifield     Field ID.                                             !
!     ncid       NetCDF file ID.                                       !
!     Lmulti     Switch to process a multi-file field. That is, the    !
!                  time records are split in several NetCDF files.     !
!     nfiles     Number of input NetCDF files.                         !
!     S          I/O derived type structure, TYPE(T_IO).               !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid       NetCDF file ID.                                       !
!     S          Updated I/O derived type structure, TYPE(T_IO).       !
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
      logical, intent(in) :: Lmulti

      integer, intent(in) :: ng, model, job, ifield, nfiles
      integer, intent(in) :: Iout, Irec, Mrec
      integer, intent(inout) :: ncid

      TYPE(T_IO), intent(inout) :: S(nfiles)
!
!  Local variable declarations.
!
      logical :: Liocycle, Lgridded, Lonerec
      logical :: foundit, got_var, got_time, special

      integer :: Fcount, Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ifile, lstr, nvatt, nvdim
      integer :: Vsize(4)

      real(r8) :: Clength, Tend, Tmax, Tmin, Tmono, Tscale, Tstr, Tval
      real(r8) :: scale

      character (len=1  ), parameter :: blank = ' '
      character (len=3  ) :: label
      character (len=256) :: Fname
!
      SourceFile='inquire.F'
!
!-----------------------------------------------------------------------
!  On first call, inquire about the contents of input NetCDF file.
!-----------------------------------------------------------------------
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
      label=S(1)%label(1:3)
      Vtype=Iinfo(1,ifield,ng)
!
!  Clear ncfile module variable.
!
      DO i=1,LEN(ncfile)
        ncfile(i:i)=blank
      END DO
!
!  If multi-files, increase (decrease if backward logic) file counter
!  and set new file names.
!
      IF (Lmulti) THEN
        DO ifile=1,nfiles
          IF (TRIM(Cinfo(ifield,ng)).eq.TRIM(S(ifile)%name)) THEN
            IF (job.gt.0) THEN
              Fcount=S(ifile)%Fcount+1
            ELSE
              Fcount=S(ifile)%Fcount-1
            END IF
            IF ((1.gt.Fcount).and.(Fcount.gt.S(ifile)%Nfiles)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,ifield)),                &
     &                            Fcount, S(ifile)%Nfiles
              END IF
              exit_flag=4
              RETURN
            END IF
            S(ifile)%Fcount=Fcount
            S(ifile)%name=TRIM(S(ifile)%files(Fcount))
            lstr=LEN_TRIM(S(ifile)%name)
            S(ifile)%base=S(ifile)%name(1:lstr-3)
            CALL netcdf_close (ng, model, ncid)
            IF (exit_flag.ne.NoError) RETURN
            IF (label(1:3).eq.'FRC') THEN
              FRCids(ifile,ng)=-1
            END IF
            IF (label(1:3).eq.'BRY') THEN
              BRYids(ifile,ng)=-1
            END IF
            IF (label(1:3).eq.'CLM') THEN
              CLMids(ifile,ng)=-1
            END IF
            EXIT
          ELSE                       ! S(ifile)%name and Fcount already
            Fcount=S(ifile)%Fcount   ! updated in first field processed
          END IF                     ! for current new file
        END DO
      ELSE
        Fcount=S(1)%Fcount
      END IF
      IF (Fcount.eq.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) Fcount, label, TRIM(Vname(1,ifield))
        END IF
        exit_flag=4
      END IF
!
!  If several input NetCDF files (nfiles>1), scan files until the
!  requested variable is found.
!
      foundit=.FALSE.
      QUERY: DO ifile=1,nfiles
!! Not sure if this is needed
        Fcount=S(ifile)%Fcount
        Fname=S(ifile)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, Fname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about requested variable.
!
        CALL netcdf_inq_var (ng, model, Fname,                          &
     &                       MyVarName = TRIM(Vname(1,ifield)),         &
     &                       SearchVar = foundit,                       &
     &                       VarID = Vid,                               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (exit_flag.ne.NoError) RETURN
!
!  Determine if gridded or point data.  Set variable dimensions.
!
        IF (foundit) THEN
          got_var=.TRUE.
          ncfile=Fname
          IF ((nvdim.gt.1).and.(ABS(job).gt.1)) THEN
            Lgridded=.TRUE.
          END IF
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
!  Check if processing special 2D fields with additional dimensions
!  (for example, tides).
!
        IF (foundit.and.(ABS(job).eq.2)) THEN
          special=.FALSE.
          DO i=1,nvdim
            IF (INDEX(TRIM(var_Dname(i)),'period').ne.0) THEN
              special=.TRUE.
            END IF
          END DO
          Linfo(4,ifield,ng)=special
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
!  to that specified in the "time" attribute.  If founded, also check
!  for a C-language null termination character in the time variable
!  string, which is possible in NetCDF files generated from C-programs.
!
          IF (Nrec.eq.0) THEN
            DO i=1,nvatt
              IF (TRIM(var_Aname(i)).eq.'time') THEN
                Tname(ifield)=TRIM(var_Achar(i))
                got_time=.TRUE.
              END IF
            END DO
            IF (got_time) THEN
              lstr=LEN_TRIM(Tname(ifield))
              DO i=1,n_dim
                IF (TRIM(dim_name(i)).eq.Tname(ifield)(1:lstr)) THEN
                  Nrec=dim_size(i)
                ELSE IF (TRIM(dim_name(i)).eq.                          &
     &                   Tname(ifield)(1:lstr-1)) THEN  ! There is a C
                  Nrec=dim_size(i)                      ! language null
                  Tname(ifield)=Tname(ifield)(1:lstr-1) ! termination
                END IF                                  ! character
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
            IF (Master) WRITE (stdout,30) TRIM(Tname(ifield)),          &
     &                                    TRIM(Vname(1,ifield)),        &
     &                                    TRIM(Fname)
            exit_flag=4
            RETURN
          END IF
        END IF
        IF (ABS(job).eq.1) THEN
          IF ((Iout.eq.1).and.(Nrec.gt.Mrec)) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,ifield)),        &
     &                                    Mrec, Nrec
            exit_flag=4
            RETURN
          END IF
        END IF
!
!  Determine initial time record to read and cycling switch.
!
        IF (foundit) THEN
          CALL get_cycle (ng, model, ifield, job, Lmulti,               &
     &                    ncfile, Tname(ifield), Nrec, tdays(ng),       &
     &                    Tid, Liocycle, Clength, Trec,                 &
     &                    Tstr, Tend, Tmin, Tmax,  Tscale)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Store variable information into global information arrays.  Also
!  fill some structure values.
!
        IF (foundit) THEN
          Linfo( 1,ifield,ng)=Lgridded
          Linfo( 2,ifield,ng)=Liocycle
          Iinfo( 2,ifield,ng)=Vid
          Iinfo( 3,ifield,ng)=Tid
          Iinfo( 4,ifield,ng)=Nrec
          Iinfo( 5,ifield,ng)=Vsize(1)
          Iinfo( 6,ifield,ng)=Vsize(2)
          Iinfo(10,ifield,ng)=S(ifile)%Nfiles
          IF (ABS(job).eq.3) THEN
            Iinfo(7,ifield,ng)=Vsize(3)
          END IF
          Finfo(1,ifield,ng)=Tmin
          Finfo(2,ifield,ng)=Tmax
          Finfo(3,ifield,ng)=Tstr
          Finfo(4,ifield,ng)=Tend
          Finfo(5,ifield,ng)=Clength
          Finfo(6,ifield,ng)=Tscale
          S(ifile)%Nrec(Fcount)=Nrec
          S(ifile)%time_min(Fcount)=Tmin
          S(ifile)%time_max(Fcount)=Tmax
          EXIT QUERY
        END IF
      END DO QUERY
!
!  Terminate execution requested variables are not found.
!
      IF (.not.got_var) THEN
        IF ((nfiles.gt.1).and.(label(1:3).eq.'FRC')) THEN
          WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
          DO i=1,nfiles
            WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
          END DO
        ELSE
          lstr=LEN_TRIM(ncfile)
          WRITE (stdout,50) TRIM(Vname(1,ifield)), 'file:'
          IF (lstr.gt.0) THEN
            WRITE (stdout,'(15x,a)') TRIM(ncfile)
          ELSE
            WRITE (stdout,'(15x,a,a)') 'file name is blank, ',          &
     &                                 'cannot be determined.'
          END IF
        END IF
        exit_flag=2
        RETURN
      END IF
      IF (.not.got_time) THEN
        IF ((nfiles.gt.1).and.(label(1:3).eq.'FRC')) THEN
          WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
          DO i=1,nfiles
            WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
          END DO
        ELSE
          WRITE (stdout,50) TRIM(Tname(ifield)), 'file:'
          WRITE (stdout,'(15x,a)') TRIM(ncfile)
        END IF
        exit_flag=2
        RETURN
      END IF
!
!  If processing model boundaries (multiple files allowed), check if file
!  for requested field (based on ncid) has been already opened and get
!  its ID from the association table (BRYids).
!
      IF (label(1:3).eq.'BRY') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (ncid.ne.BRYids(ifile,ng))) THEN
            ncid=BRYids(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
!
!  If processing model climatology (multiple files allowed), check if file
!  for requested field (based on ncid) has been already opened and get
!  its ID from the association table (CLMids).
!
      IF (label(1:3).eq.'CLM') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (ncid.ne.CLMids(ifile,ng))) THEN
            ncid=CLMids(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
!
!  If processing model forcing (multiple files allowed), check if file
!  for requested field (based on ncid) has been already opened and get
!  its ID from the association table (FRCids).
!
      IF (label(1:3).eq.'FRC') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (ncid.ne.FRCids(ifile,ng))) THEN
            ncid=FRCids(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
!
!  If appropriate, open input NetCDF file for reading.
!
      IF (ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncfile, 0, ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,60) TRIM(ncfile)
          RETURN
        END IF
!
!  If processing model boundaries, load ID into association table (BRYids).
!
        IF (label(1:3).eq.'BRY') THEN
          DO ifile=1,nfiles
            IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
              BRYids(ifile,ng)=ncid
              S(ifile)%ncid=ncid
              EXIT
            END IF
          END DO
!
!  If processing model climatology, load ID into association table (CLMids).
!
        ELSE IF (label(1:3).eq.'CLM') THEN
          DO ifile=1,nfiles
            IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
              CLMids(ifile,ng)=ncid
              S(ifile)%ncid=ncid
              EXIT
            END IF
          END DO
!
!  If processing model forcing, load ID into association table (FRCids).
!
        ELSE IF (label(1:3).eq.'FRC') THEN
          DO ifile=1,nfiles
            IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
              FRCids(ifile,ng)=ncid
              S(ifile)%ncid=ncid
              EXIT
            END IF
          END DO
        ELSE
          S(1)%ncid=ncid
        END IF
      END IF
      Cinfo(ifield,ng)=TRIM(ncfile)
!
!  The strategy here is to create a local, monotonically increasing
!  time variable so the interpolation between snapshots is trivial
!  when cycling forcing fields. Subtract one to time record counter
!  "Trec" to avoid doing special case at initialization. If the
!  job task is negative, the logic is backaward in time for the
!  adjoint model.
!
      IF (.not.Lmulti) THEN
        IF (job.lt.0) THEN            ! backward time logic, adjoint
          IF (Irec.eq.1) THEN
            Tindex=Iout
          ELSE
            Tindex=1
          END IF
          IF (Liocycle) THEN
            IF (Trec.eq.nrec) THEN
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
            CALL netcdf_get_fvar (ng, model, ncfile, Tname(ifield),     &
     &                            Tval,                                 &
     &                            ncid = ncid,                          &
     &                            start = (/Trec-1/),                   &
     &                            total = (/1/))
            IF (exit_flag.ne.NoError) THEN
              IF (Master) WRITE (stdout,70) TRIM(Tname(ifield)), Trec
              RETURN
            END IF
            Tval=Tval*Tscale
            IF (Tval.lt.Tend) THEN
              Trec=Trec+1
            END IF
          END IF
        ELSE                          ! forward time logic
          IF (Irec.eq.1) THEN
            Tindex=Iout
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
          Trec=Trec-1
        END IF
        Tmono=Tmono*day2sec
        Iinfo(8,ifield,ng)=Tindex
        Iinfo(9,ifield,ng)=Trec
        Finfo(7,ifield,ng)=Tmono
      ELSE
        Iinfo(9,ifield,ng)=Trec
      END IF
!
!  Set switch for one time record dataset. In this case, the input
!  data is always the same and time interpolation is not performed.
!
      IF (Nrec.eq.1) Lonerec=.TRUE.
      Linfo(3,ifield,ng)=Lonerec
      Tindex=Iinfo(8,ifield,ng)
      IF (job.lt.0) THEN
        Vtime(Tindex,ifield,ng)=Finfo(4,ifield,ng)
      ELSE
        Vtime(Tindex,ifield,ng)=Finfo(3,ifield,ng)
      END IF
!
  10  FORMAT (/,' INQUIRE     - out of range multi-files counter for ', &
     &        'variable; ',a,/,15x,'Fcount = ',i2.2,                    &
     &        ',   Expected range: 1 - ',i2.2)
  20  FORMAT (/,' INQUIRE     - unable to assign file counter, ',       &
     &        'Fcount = ',i4,/,15x,'while processing structure: ',a,    &
     &        /,15x,'and variable; ',a)
  30  FORMAT (/,' INQUIRE     - unable to find dimension ',a,           &
     &        /,15x,'for variable: ',a,/,15x,'in file: ',a,             &
     &        /,15x,'file is not CF compliant...')
  40  FORMAT (/,' INQUIRE     - too small dimension for variable ',a,   &
     &        ': ',2i5)
  50  FORMAT (/,' INQUIRE     - unable to find requested variable: ',a, &
     &        /,15x,'in ',a)
  60  FORMAT (/,' INQUIRE     - unable to open input NetCDF file: ',a)
  70  FORMAT (/,' INQUIRE     - error while reading variable: ',a,2x,   &
     &        ' at TIME index = ',i4)

      RETURN
      END SUBROUTINE inquire
