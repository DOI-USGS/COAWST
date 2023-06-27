#include "cppdefs.h"
      MODULE inquiry_mod
!
!git $Id$
!svn $Id: inquiry.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
#if defined PIO_LIB && defined DISTRIBUTE
!or   pioFile    PIO file descriptor structure, TYPE(file_desc_t)      !
!                  pioFile%fh         file handler                     !
!                  pioFile%iosystem   IO system descriptor (struct)    !
#endif
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
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
!
      USE get_cycle_mod, ONLY : get_cycle
      USE strings_mod,   ONLY : FoundError, lowercase
!
      implicit none
!
!  Interface for same name routine overloading.
!
      INTERFACE inquiry
        MODULE PROCEDURE inquiry_nf90
#if defined PIO_LIB && defined DISTRIBUTE
        MODULE PROCEDURE inquiry_pio
#endif
      END INTERFACE inquiry
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE inquiry_nf90 (ng, model, job, Iout, Irec, Mrec,        &
     &                         ifield, ncid, Lmulti, nfiles, S)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lmulti
!
      integer, intent(in) :: ng, model, job, ifield, nfiles
      integer, intent(in) :: Iout, Irec, Mrec
      integer, intent(inout) :: ncid
!
      TYPE (T_IO), intent(inout) :: S(nfiles)
!
!  Local variable declarations.
!
      logical :: CloseFile, Liocycle, Lgridded, Lonerec
      logical :: foundit, got_var, got_time, special
!
      integer :: Fcount, Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ifile, lstr, my_ncid, nvatt, nvdim
      integer :: Vsize(4)
!
      real(dp) :: Clength, Tend, Tmax, Tmin, Tmono, Tscale, Tstr, Tval
      real(dp) :: scale
!
      character (len=1  ), parameter :: blank = ' '
      character (len=3  ) :: label
      character (len=256) :: Fname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", inquiry_nf90"
!
      SourceFile=MyFile
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
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            S(ifile)%Fcount=Fcount
            S(ifile)%name=TRIM(S(ifile)%files(Fcount))
            CALL netcdf_close (ng, model, ncid)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (label(1:3).eq.'FRC') THEN
              FRCids(ifile,ng)=-1
              S(ifile)%ncid=-1
            END IF
            IF (label(1:3).eq.'BRY') THEN
              BRYids(ifile,ng)=-1
              S(ifile)%ncid=-1
            END IF
            IF (label(1:3).eq.'CLM') THEN
              CLMids(ifile,ng)=-1
              S(ifile)%ncid=-1
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  If several input NetCDF files (nfiles>1), scan files until the
!  requested variable is found.
!
      foundit=.FALSE.
      QUERY: DO ifile=1,nfiles
        Fname=S(ifile)%name
!
!  Open NetCDF file for reading.
!
        IF (S(ifile)%ncid.eq.-1) THEN
          CALL netcdf_open (ng, model, Fname, 0, my_ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Fname)
            RETURN
          END IF
          CloseFile=.TRUE.
        ELSE
          my_ncid=S(ifile)%ncid
          CloseFile=.FALSE.
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, Fname,                        &
     &                         ncid = my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about requested variable.
!
        CALL netcdf_inq_var (ng, model, Fname,                          &
     &                       ncid = my_ncid,                            &
     &                       MyVarName = TRIM(Vname(1,ifield)),         &
     &                       SearchVar = foundit,                       &
     &                       VarID = Vid,                               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)

        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        ELSE
          IF (.not.Lmulti) THEN
            ncfile=Fname                    ! need for error report
          END IF
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
!  "varinfo.yaml" (Nrec still zero), reset associated time-variable name
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
              IF (INDEX(TRIM(lowercase(var_Dname(i))),'time').ne.0) THEN
                Nrec=var_Dsize(i)
              END IF
            END DO
          END IF
          IF (got_time.and.(Nrec.eq.0)) THEN
            IF (Master) WRITE (stdout,30) TRIM(Tname(ifield)),          &
     &                                    TRIM(Vname(1,ifield)),        &
     &                                    TRIM(Fname)
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
        IF (ABS(job).eq.1) THEN
          IF ((Iout.eq.1).and.(Nrec.gt.Mrec)) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,ifield)),        &
     &                                    Mrec, Nrec
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
!
!  Determine initial time record to read and cycling switch.
!
        IF (foundit) THEN
          CALL get_cycle (ng, model, ifield, job, Lmulti,               &
     &                    ncfile, my_ncid, Tname(ifield), Nrec,         &
     &                    tdays(ng), Tid, Liocycle, Clength, Trec,      &
     &                    Tstr, Tend, Tmin, Tmax,  Tscale)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          SourceFile=MyFile
        END IF
!
!  Store variable information into global information arrays.  Also
!  fill some structure values.
!
        IF (foundit) THEN
          ncid=my_ncid
          Linfo( 1,ifield,ng)=Lgridded
          Linfo( 2,ifield,ng)=Liocycle
          Iinfo( 2,ifield,ng)=Vid
          Iinfo( 3,ifield,ng)=Tid
          Iinfo( 4,ifield,ng)=Nrec
          Iinfo( 5,ifield,ng)=Vsize(1)
          Iinfo( 6,ifield,ng)=Vsize(2)
          Iinfo( 7,ifield,ng)=Vsize(3)
          Iinfo(10,ifield,ng)=S(ifile)%Nfiles
          Iinfo(11,ifield,ng)=nvdim
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
!
!  Close input NetCDF file if opened during the query.  Files opened
!  outside the query loop remain open.  This is done to avoid opening
!  too many files.
!
        IF (CloseFile) THEN
          CALL netcdf_close (ng, model, my_ncid, Fname, .FALSE.)
        END IF
      END DO QUERY
!
!  Terminate execution requested variables are not found.
!
      IF (.not.got_var) THEN
        IF ((nfiles.gt.1).and.(label(1:3).eq.'FRC')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'BRY')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'CLM')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE
          lstr=LEN_TRIM(ncfile)
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'file:'
            IF (lstr.gt.0) THEN
              WRITE (stdout,'(15x,a)') TRIM(ncfile)
            ELSE
              WRITE (stdout,'(15x,a,a)') 'file name is blank, ',        &
     &                                   'cannot be determined.'
            END IF
          END IF
        END IF
        exit_flag=2
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
      IF (.not.got_time) THEN
        IF ((nfiles.gt.1).and.(label(1:3).eq.'FRC')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'BRY')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'CLM')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE
          lstr=LEN_TRIM(ncfile)
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'file:'
            IF (lstr.gt.0) THEN
               WRITE (stdout,'(15x,a)') TRIM(ncfile)
            ELSE
              WRITE (stdout,'(15x,a,a)') 'file name is blank, ',        &
     &                                   'cannot be determined.'
            END IF
          END IF
        END IF
        exit_flag=2
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  If processing model forcing (multiple files allowed), check if file
!  for requested field (based on ncid) has been already opened and get
!  its ID from the association table (FRCids).
!
      IF (label(1:3).eq.'FRC') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (FRCids(ifile,ng).ne.-1).and.                             &
     &        (ncid.ne.FRCids(ifile,ng))) THEN
            ncid=FRCids(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
      IF (label(1:3).eq.'BRY') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (BRYids(ifile,ng).ne.-1).and.                             &
     &        (ncid.ne.BRYids(ifile,ng))) THEN
            ncid=BRYids(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
      IF (label(1:3).eq.'CLM') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (CLMids(ifile,ng).ne.-1).and.                             &
     &        (ncid.ne.CLMids(ifile,ng))) THEN
            ncid=CLMids(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
!
!  If appropriate, open input NetCDF file for reading.
!  (HGA: it is not needed since file was already open and
!        ncid=my_ncid was assigned above.
!
      IF (ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncfile, 0, ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
           IF (Master) WRITE (stdout,60) TRIM(ncfile)
          RETURN
        END IF
      END IF
!
!  If processing model forcing, load ID into association table (FRCids).
!
      IF (label(1:3).eq.'FRC') THEN
        DO ifile=1,nfiles
          IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
            FRCids(ifile,ng)=ncid
            S(ifile)%ncid=ncid
            EXIT
          END IF
        END DO
      ELSE IF (label(1:3).eq.'BRY') THEN
        DO ifile=1,nfiles
          IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
            BRYids(ifile,ng)=ncid
            S(ifile)%ncid=ncid
            EXIT
          END IF
        END DO
      ELSE IF (label(1:3).eq.'CLM') THEN
        DO ifile=1,nfiles
          IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
            CLMids(ifile,ng)=ncid
            S(ifile)%ncid=ncid
            EXIT
          END IF
        END DO
      ELSE
        S(1)%ncid=ncid
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
            Tindex=Iinfo(8,ifield,ng)
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
            CALL netcdf_get_time (ng, model, ncfile, Tname(ifield),     &
     &                            Rclock%DateNumber, Tval,              &
     &                            ncid = ncid,                          &
     &                            start = (/Trec-1/),                   &
     &                            total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
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
  10  FORMAT (/,' INQUIRY_NF90 - out of range multi-files counter for', &
     &        ' variable: ',a,/,16x,'Fcount = ',i0,                     &
     &        ',   Expected range: 1 - ',i0)
  20  FORMAT (/,' INQUIRY_NF90 - unable to assign file counter, ',      &
     &        'Fcount = ',i0,/,16x,'while processing structure: ',a,    &
     &        /,16x,'and variable; ',a)
  30  FORMAT (/,' INQUIRY_NF90 - unable to find dimension ',a,          &
     &        /,16x,'for variable: ',a,/,16x,'in file: ',a,             &
     &        /,16x,'file is not CF compliant...')
  40  FORMAT (/,' INQUIRY_NF90 - too small dimension for variable ',a,  &
     &        ': ',i0,2x,i0)
  50  FORMAT (/,' INQUIRY_NF90 - unable to find requested variable: ',  &
     &        a,/,16x,'in ',a)
  60  FORMAT (/,' INQUIRY_NF90 - unable to open input NetCDF file: ',a)
  70  FORMAT (/,' INQUIRY_NF90 - error while reading variable: ',a,2x,  &
     &        ' at TIME index = ',i0)
!
      RETURN
      END SUBROUTINE inquiry_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE inquiry_pio (ng, model, job, Iout, Irec, Mrec,         &
     &                        ifield, pioFile, Lmulti, nfiles, S)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lmulti
!
      integer, intent(in) :: ng, model, job, ifield, nfiles
      integer, intent(in) :: Iout, Irec, Mrec
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE (T_IO),        intent(inout) :: S(nfiles)
!
!  Local variable declarations.
!
      logical :: CloseFile, Liocycle, Lgridded, Lonerec
      logical :: foundit, got_var, got_time, special
!
      integer :: Fcount, Nrec, Tindex, Trec, Vtype
      integer :: i, ifile, lstr, nvatt, nvdim
      integer :: Vsize(4)
!
      real(dp) :: Clength, Tend, Tmax, Tmin, Tmono, Tscale, Tstr, Tval
      real(dp) :: scale
!
      character (len=1  ), parameter :: blank = ' '
      character (len=3  ) :: label
      character (len=256) :: Fname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", inquiry_pio"
!
      TYPE (File_desc_t) :: my_pioFile
      TYPE (Var_desc_t)  :: pioVar, TpioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  On first call, inquire about the contents of input NetCDF file.
!-----------------------------------------------------------------------
!
!  Intialize local variables.
!
      pioVar%varID=-1
      TpioVar%varID=-1
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
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
            S(ifile)%Fcount=Fcount
            S(ifile)%name=TRIM(S(ifile)%files(Fcount))
            CALL pio_netcdf_close (ng, model, pioFile)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
            IF (label(1:3).eq.'FRC') THEN
              FRCdesc(ifile,ng)%fh=-1
              S(ifile)%pioFile%fh=-1
            END IF
            IF (label(1:3).eq.'BRY') THEN
              BRYdesc(ifile,ng)%fh=-1
              S(ifile)%pioFile%fh=-1
            END IF
            IF (label(1:3).eq.'CLM') THEN
              CLMdesc(ifile,ng)%fh=-1
              S(ifile)%pioFile%fh=-1
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  If several input NetCDF files (nfiles>1), scan files until the
!  requested variable is found.
!
      foundit=.FALSE.
      QUERY: DO ifile=1,nfiles
        Fname=S(ifile)%name
!
!  Open NetCDF file for reading.
!
        IF (S(ifile)%pioFile%fh.eq.-1) THEN
          CALL pio_netcdf_open (ng, model, Fname, 0, my_pioFile)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Fname)
            RETURN
          END IF
          CloseFile=.TRUE.
        ELSE
          my_pioFile=S(ifile)%pioFile
          CloseFile=.FALSE.
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, Fname,                    &
     &                             pioFile = my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about requested variable.
!
        CALL pio_netcdf_inq_var (ng, model, Fname,                      &
     &                           pioFile = my_pioFile,                  &
     &                           MyVarName = TRIM(Vname(1,ifield)),     &
     &                           SearchVar = foundit,                   &
     &                           pioVar = pioVar,                       &
     &                           nVarDim = nvdim,                       &
     &                           nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        ELSE
          IF (.not.Lmulti) THEN
            ncfile=Fname                    ! need for error report
          END IF
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
!  "varinfo.yaml" (Nrec still zero), reset associated time-variable name
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
              IF (INDEX(TRIM(lowercase(var_Dname(i))),'time').ne.0) THEN
                Nrec=var_Dsize(i)
              END IF
            END DO
          END IF
          IF (got_time.and.(Nrec.eq.0)) THEN
            IF (Master) WRITE (stdout,30) TRIM(Tname(ifield)),          &
     &                                    TRIM(Vname(1,ifield)),        &
     &                                    TRIM(Fname)
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
        IF (ABS(job).eq.1) THEN
          IF ((Iout.eq.1).and.(Nrec.gt.Mrec)) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,ifield)),        &
     &                                    Mrec, Nrec
            exit_flag=4
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END IF
!
!  Determine initial time record to read and cycling switch.
!
        IF (foundit) THEN
          CALL get_cycle (ng, model, ifield, job, Lmulti,               &
     &                    ncfile, my_pioFile, Tname(ifield), Nrec,      &
     &                    tdays(ng), TpioVar, Liocycle, Clength, Trec,  &
     &                    Tstr, Tend, Tmin, Tmax,  Tscale)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          SourceFile=MyFile
        END IF
!
!  Store variable information into global information arrays.  Also
!  fill some structure values.
!
        IF (foundit) THEN
          pioFile=my_pioFile
          Dinfo( 1,ifield,ng)%vd=pioVar
          Dinfo( 1,ifield,ng)%dkind=PIO_TYPE
          Dinfo( 1,ifield,ng)%gtype=Vtype
          Dinfo( 2,ifield,ng)%vd=TpioVar
          Dinfo( 1,ifield,ng)%dkind=PIO_TOUT
          Dinfo( 2,ifield,ng)%gtype=0
          Linfo( 1,ifield,ng)=Lgridded
          Linfo( 2,ifield,ng)=Liocycle
          Iinfo( 4,ifield,ng)=Nrec
          Iinfo( 5,ifield,ng)=Vsize(1)
          Iinfo( 6,ifield,ng)=Vsize(2)
          Iinfo( 7,ifield,ng)=Vsize(3)
          Iinfo(10,ifield,ng)=S(ifile)%Nfiles
          Iinfo(11,ifield,ng)=nvdim
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
!
!  Close input NetCDF file if opened during the query.  Files opened
!  outside the query loop remain open.  This is done to avoid opening
!  too many files.
!
        IF (CloseFile) THEN
          CALL pio_netcdf_close (ng, model, my_pioFile, Fname, .FALSE.)
        END IF
      END DO QUERY
!
!  Terminate execution requested variables are not found.
!
      IF (.not.got_var) THEN
        IF ((nfiles.gt.1).and.(label(1:3).eq.'FRC')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'BRY')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'CLM')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE
          lstr=LEN_TRIM(ncfile)
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Vname(1,ifield)), 'file:'
            IF (lstr.gt.0) THEN
              WRITE (stdout,'(15x,a)') TRIM(ncfile)
            ELSE
              WRITE (stdout,'(15x,a,a)') 'file name is blank, ',        &
     &                                   'cannot be determined.'
            END IF
          END IF
        END IF
        exit_flag=2
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
      IF (.not.got_time) THEN
        IF ((nfiles.gt.1).and.(label(1:3).eq.'FRC')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'BRY')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE IF ((nfiles.gt.1).and.(label(1:3).eq.'CLM')) THEN
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'files:'
            DO i=1,nfiles
              WRITE (stdout,'(15x,a)') TRIM(S(i)%name)
            END DO
          END IF
        ELSE
          lstr=LEN_TRIM(ncfile)
          IF (Master) THEN
            WRITE (stdout,50) TRIM(Tname(ifield)), 'file:'
            IF (lstr.gt.0) THEN
               WRITE (stdout,'(15x,a)') TRIM(ncfile)
            ELSE
              WRITE (stdout,'(15x,a,a)') 'file name is blank, ',        &
     &                                   'cannot be determined.'
            END IF
          END IF
        END IF
        exit_flag=2
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  If processing model forcing (multiple files allowed), check if file
!  for requested field (based on pioFile) has been already opened and
!  get its file descriptor from the association table (FRCdesc).
!
      IF (label(1:3).eq.'FRC') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (FRCdesc(ifile,ng)%fh.ne.-1).and.                         &
     &        (pioFile%fh.ne.FRCdesc(ifile,ng)%fh)) THEN
            pioFile=FRCdesc(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
      IF (label(1:3).eq.'BRY') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (BRYdesc(ifile,ng)%fh.ne.-1).and.                         &
     &        (pioFile%fh.ne.BRYdesc(ifile,ng)%fh)) THEN
            pioFile=BRYdesc(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
      IF (label(1:3).eq.'CLM') THEN
        DO ifile=1,nfiles
          IF ((TRIM(ncfile).eq.TRIM(S(ifile)%name)).and.                &
     &        (CLMdesc(ifile,ng)%fh.ne.-1).and.                         &
     &        (pioFile%fh.ne.CLMdesc(ifile,ng)%fh)) THEN
            pioFile=CLMdesc(ifile,ng)
            EXIT
          END IF
        END DO
      END IF
!
!  If appropriate, open input NetCDF file for reading.
!  (HGA: it is not needed since file was already open and
!        pioFile=my_pioFile was assigned above.
!
      IF (pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncfile, 0, pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
           IF (Master) WRITE (stdout,60) TRIM(ncfile)
          RETURN
        END IF
      END IF
!
!  If processing model forcing, load ID into association table (FRCids).
!
      IF (label(1:3).eq.'FRC') THEN
        DO ifile=1,nfiles
          IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
            FRCdesc(ifile,ng)=pioFile
            S(ifile)%pioFile=pioFile
            EXIT
          END IF
        END DO
      ELSE IF (label(1:3).eq.'BRY') THEN
        DO ifile=1,nfiles
          IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
            BRYdesc(ifile,ng)=pioFile
            S(ifile)%pioFile=pioFile
            EXIT
          END IF
        END DO
      ELSE IF (label(1:3).eq.'CLM') THEN
        DO ifile=1,nfiles
          IF (TRIM(ncfile).eq.TRIM(S(ifile)%name)) THEN
            CLMdesc(ifile,ng)=pioFile
            S(ifile)%pioFile=pioFile
            EXIT
          END IF
        END DO
      ELSE
        S(1)%pioFile=piofile
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
            Tindex=Iinfo(8,ifield,ng)
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
            CALL pio_netcdf_get_time (ng, model, ncfile, Tname(ifield), &
     &                                Rclock%DateNumber, Tval,          &
     &                                pioFile = pioFile,                &
     &                                start = (/Trec-1/),               &
     &                                total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
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
  10  FORMAT (/,' INQUIRY_PIO  - out of range multi-files counter for', &
     &        ' variable: ',a,/,16x,'Fcount = ',i0,                     &
     &        ',   Expected range: 1 - ',i0)
  20  FORMAT (/,' INQUIRY_PIO  - unable to assign file counter, ',      &
     &        'Fcount = ',i0,/,16x,'while processing structure: ',a,    &
     &        /,16x,'and variable; ',a)
  30  FORMAT (/,' INQUIRY_PIO  - unable to find dimension ',a,          &
     &        /,16x,'for variable: ',a,/,16x,'in file: ',a,             &
     &        /,16x,'file is not CF compliant...')
  40  FORMAT (/,' INQUIRY_PIO  - too small dimension for variable ',a,  &
     &        ': ',i0,2x,i0)
  50  FORMAT (/,' INQUIRY_PIO  - unable to find requested variable: ',  &
     &        a,/,16x,'in ',a)
  60  FORMAT (/,' INQUIRY_PIO  - unable to open input NetCDF file: ',a)
  70  FORMAT (/,' INQUIRY_PIO  - error while reading variable: ',a,2x,  &
     &        ' at TIME index = ',i0)
!
      RETURN
      END SUBROUTINE inquiry_pio
# endif

      END MODULE inquiry_mod
