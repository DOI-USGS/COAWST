#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE get_ngfldr (ng, model, ifield, ncid,                   &
# if defined PIO_LIB && defined DISTRIBUTE
     &                       pioFile,                                   &
# endif
     &                       nfiles, S, recordless, update,             &
     &                       LBi, UBi, UBj, UBk, Istr, Iend, Jrec,      &
     &                       Fout)
!
!git $Id$
!svn $Id: get_ngfldr.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
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
!     S          I/O derived type structure, TYPE(T_IO).               !
!     recordless Switch for time invariant field (logical).            !
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
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in)  :: recordless
      logical, intent(out) :: update
!
      integer, intent(in) :: ng, model, ifield, nfiles
      integer, intent(in) :: LBi, UBi, UBj, UBk, Istr, Iend, Jrec
      integer, intent(inout) :: ncid
!
# if defined PIO_LIB && defined DISTRIBUTE
      TYPE (File_desc_t), intent(inout) :: pioFile
# endif
      TYPE(T_IO), intent(inout) :: S(nfiles)
!
      real(r8), intent(inout) :: Fout(LBi:UBi,UBj,UBk)
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Reads in requested non-grided field according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (S(1)%IOtype)
        CASE (io_nf90)
          CALL get_ngfldr_nf90 (ng, model, ifield, ncid,                &
     &                          nfiles, S, recordless, update,          &
     &                          LBi, UBi, UBj, UBk, Istr, Iend, Jrec,   &
     &                          Fout)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_ngfldr_pio (ng, model, ifield, pioFile,              &
     &                         nfiles, S, recordless, update,           &
     &                         LBi, UBi, UBj, UBk, Istr, Iend, Jrec,    &
     &                         Fout)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) S(1)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_NGFLDR - Illegal input file type, io_type = ',i0,   &
     &        /,12x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_ngfldr
!
!***********************************************************************
      SUBROUTINE get_ngfldr_nf90 (ng, model, ifield, ncid,              &
     &                            nfiles, S, recordless, update,        &
     &                            LBi, UBi, UBj, UBk, Istr, Iend, Jrec, &
     &                            Fout)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE dateclock_mod, ONLY : time_string
# ifdef CHECKSUM
      USE get_hash_mod,  ONLY : get_hash
# endif
      USE inquiry_mod,   ONLY : inquiry
      USE strings_mod,   ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in)  :: recordless
      logical, intent(out) :: update
!
      integer, intent(in) :: ng, model, ifield, nfiles
      integer, intent(in) :: LBi, UBi, UBj, UBk, Istr, Iend, Jrec
      integer, intent(inout) :: ncid
!
      TYPE(T_IO), intent(inout) :: S(nfiles)
!
      real(r8), intent(inout) :: Fout(LBi:UBi,UBj,UBk)
!
!  Local variable declarations.
!
      logical :: Linquire, Liocycle, Lmulti, Lonerec
!
      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ic, j, job, lend, lstr, npts, nvdim, status
# ifdef CHECKSUM
      integer(i8b) :: Fhash
# endif
!
      real(r8) :: Aval, Fmax, Fmin

      real(dp) :: Clength, Tdelta, Tend
      real(dp) :: Tmax, Tmin, Tmono, Tscale, Tstr
      real(dp) :: Tsec, Tval

      real(r8), dimension((UBi-LBi+1)*UBj) :: Awrk
!
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_ngfldr_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Determine if inquiring about the field to process in input NetCDF
!  file(s).  This usually happens on first call or when the field
!  time records are split (saved) in several multi-files.
!
      Linquire=.FALSE.
      Lmulti=.FALSE.
      IF (iic(ng).eq.0) Linquire=.TRUE.
      IF (.not.Linquire.and.                                            &
     &    ((Iinfo(10,ifield,ng).gt.1).and.                              &
     &     (Linfo( 5,ifield,ng).or.                                     &
     &     (Finfo( 1,ifield,ng)*day2sec.gt.time(ng))))) THEN
        Linquire=.TRUE.
        Lmulti=.TRUE.
      END IF
!
!  If appropriate, inquire about the contents of input NetCDF file and
!  fill information arrays.
!
!  Also, if appropriate, deactivate the Linfo(5,ifield,ng) switch after
!  the query for the LOWER snapshot interpolant from previous multifile
!  in the list. The switch was activated previously to indicate the
!  processing of the FIRST record of the file for the UPPER snapshot
!  interpolant.
!
      IF (Linquire) THEN
        job=-1
        CALL inquiry (ng, model, job, UBk, Iend, UBi, ifield, ncid,     &
     &                Lmulti, nfiles, S)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        IF (Linfo(5,ifield,ng)) THEN
          Linfo(5,ifield,ng)=.FALSE.
        END IF
        IF (Linfo(6,ifield,ng)) THEN
          Linfo(6,ifield,ng)=.FALSE.
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, read in new data.
!-----------------------------------------------------------------------
!
      update=.FALSE.
      Tmono=Finfo(7,ifield,ng)
!
      IF ((Tmono.gt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart(ng))) THEN
!
!  Load properties for requested field from information arrays.
!
        Liocycle=Linfo( 2,ifield,ng)
        Lonerec =Linfo( 3,ifield,ng)
        Vtype   =Iinfo( 1,ifield,ng)
        Vid     =Iinfo( 2,ifield,ng)
        Tid     =Iinfo( 3,ifield,ng)
        Nrec    =Iinfo( 4,ifield,ng)
        Tindex  =Iinfo( 8,ifield,ng)
        Trec    =Iinfo( 9,ifield,ng)
        nvdim   =Iinfo(11,ifield,ng)
        Tmin    =Finfo( 1,ifield,ng)
        Tmax    =Finfo( 2,ifield,ng)
        Clength =Finfo( 5,ifield,ng)
        Tscale  =Finfo( 6,ifield,ng)
        ncfile  =Cinfo(ifield,ng)
!
        IF (Liocycle) THEN
          Trec=MOD(Trec,Nrec)-1
          IF (Trec.le.0) Trec=Nrec+Trec
        ELSE
          Trec=Trec-1
        END IF
        Iinfo(9,ifield,ng)=Trec
!
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
          IF (.not.recordless.and.(Tid.ge.0).and.(Tid.ne.Vid)) THEN
            CALL netcdf_get_time (ng, model, ncfile, Tname(ifield),     &
     &                            Rclock%DateNumber, Tval,              &
     &                            ncid = ncid,                          &
     &                            start = (/Trec/),                     &
     &                            total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,50) TRIM(Tname(ifield)), Trec
              RETURN
            END IF
            Tval=Tval*Tscale
            Vtime(Tindex,ifield,ng)=Tval
!
!  Activate switch Linfo(5,ifield,ng) if processing the FIRST record of
!  the file for the UPPER time snapshot. We need to get the LOWER time
!  snapshot from previous multifile.
!
            IF ((Trec.eq.1).and.(Tval*day2sec.ge.time(ng))) THEN
              Linfo(5,ifield,ng)=.TRUE.
            END IF
          END IF
!
!  Read in non-grided data. The conditional statement on Jrec is to
!  differentiate between reading a 3D and 2D array.
!
          IF (Vid.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF (nvdim.eq.1) THEN
              npts=Iend-Istr+1
              CALL netcdf_get_fvar (ng, model, ncfile,                  &
     &                              Vname(1,ifield), Awrk,              &
     &                              ncid = ncid,                        &
     &                              start = (/1/),                      &
     &                              total = (/Iend-Istr+1/))
            ELSE IF (nvdim.eq.2) THEN
              IF (recordless) THEN
                npts=(Iend-Istr+1)*Jrec
                CALL netcdf_get_fvar (ng, model, ncfile,                &
     &                                Vname(1,ifield), Awrk,            &
     &                                ncid = ncid,                      &
     &                                start = (/1,1/),                  &
     &                                total = (/Iend-Istr+1,Jrec/))
              ELSE
                npts=Iend-Istr+1
                CALL netcdf_get_fvar (ng, model, ncfile,                &
     &                                Vname(1,ifield), Awrk,            &
     &                                ncid = ncid,                      &
     &                                start = (/1,Trec/),               &
     &                                total = (/Iend-Istr+1,1/))
              END IF
            ELSE IF (nvdim.eq.3) THEN
              npts=(Iend-Istr+1)*Jrec
              CALL netcdf_get_fvar (ng, model, ncfile,                  &
     &                              Vname(1,ifield), Awrk,              &
     &                              ncid = ncid,                        &
     &                              start = (/1,1,Trec/),               &
     &                              total = (/Iend-Istr+1,Jrec,1/))
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)), Trec
              RETURN
            END IF
            Fmin=Awrk(1)*Fscale(ifield,ng)
            Fmax=Awrk(1)*Fscale(ifield,ng)
            ic=0
            DO j=1,Jrec
              DO i=Istr,Iend
                ic=ic+1
                Aval=Awrk(ic)*Fscale(ifield,ng)
                Fmin=MIN(Fmin,Aval)
                Fmax=MAX(Fmax,Aval)
                Fout(i,j,Tindex)=Aval
              END DO
            END DO
            Finfo(8,ifield,ng)=Fmin
            Finfo(9,ifield,ng)=Fmax
#ifdef CHECKSUM
            CALL get_hash (Awrk, npts, Fhash)
#endif
            IF (Master) THEN
              IF (UBk.eq.1) THEN
                WRITE (stdout,60) TRIM(Vname(2,ifield)), ng, Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                Tsec=Tval*day2sec
                CALL time_string (Tsec, t_code)
                WRITE (stdout,70) TRIM(Vname(2,ifield)), t_code,        &
     &                            ng, Trec, Tindex, ncfile(lstr:lend),  &
     &                            Tmin, Tmax, Tval, Fmin, Fmax
              END IF
# ifdef CHECKSUM
              WRITE (stdout,80) Fhash
# endif
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
        IF (.not.Lonerec.and.(.not.recordless)) THEN
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
  10  FORMAT (/,' GET_NGFLDR_NF90 - unable to find dimension ',a,       &
     &        /,19x,'for variable: ',a,/,19x,'in file: ',a,             &
     &        /,19x,'file is not CF compliant...')
  20  FORMAT (/,' GET_NGFLDR_NF90 - too small dimension for variable',  &
     &        1x,a,': ',i0,2x,i0)
  30  FORMAT (/,' GET_NGFLDR_NF90 - unable to find requested',          &
     &        ' variable: ',a,/,19x,'in file: ',a)
  40  FORMAT (/,' GET_NGFLDR_NF90 - unable to open input NetCDF',       &
     &        ' file: ',a)
  50  FORMAT (/,' GET_NGFLDR_NF90 - error while reading variable: ',a,  &
     &        2x,' at TIME index = ',i0)
  60  FORMAT (2x,'GET_NGFLDR_NF90  - ',a,/,22x,'(Grid = ',i2.2,         &
     &        ', Min = ',1pe15.8,' Max = ',1pe15.8,')')
  70  FORMAT (2x,'GET_NGFLDR_NF90  - ',a,',',t75,a,/,22x,               &
     &        '(Grid=',i2.2,', Rec=',i0,', Index=',i1,                  &
     &        ', File: ',a,')',/,22x,                                   &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',                   &
     &        t71, 't = ', f15.4 ,/,22x,                                &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
# ifdef CHECKSUM
  80  FORMAT (22x,'(CheckSum = ',i0,')')
# endif
!
      RETURN
      END SUBROUTINE get_ngfldr_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_ngfldr_pio (ng, model, ifield, pioFile,            &
     &                           nfiles, S, recordless, update,         &
     &                           LBi, UBi, UBj, UBk, Istr, Iend, Jrec,  &
     &                           Fout)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_pio_netcdf
      USE mod_scalars
!
      USE dateclock_mod, ONLY : time_string
#  ifdef CHECKSUM
      USE get_hash_mod,  ONLY : get_hash
#  endif
      USE inquiry_mod,   ONLY : inquiry
      USE strings_mod,   ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in)  :: recordless
      logical, intent(out) :: update
!
      integer, intent(in) :: ng, model, ifield, nfiles
      integer, intent(in) :: LBi, UBi, UBj, UBk, Istr, Iend, Jrec
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE(T_IO), intent(inout) :: S(nfiles)
!
      real(r8), intent(inout) :: Fout(LBi:UBi,UBj,UBk)
!
!  Local variable declarations.
!
      logical :: Linquire, Liocycle, Lmulti, Lonerec
!
      integer :: Nrec, Tindex, Trec, Vtype
      integer :: i, ic, j, job, lend, lstr, npts, nvdim, status
#  ifdef CHECKSUM
      integer(i8b) :: Fhash
#  endif
!
      real(r8) :: Aval, Fmax, Fmin

      real(dp) :: Clength, Tdelta, Tend
      real(dp) :: Tmax, Tmin, Tmono, Tscale, Tstr
      real(dp) :: Tsec, Tval

      real(r8), dimension((UBi-LBi+1)*UBj) :: Awrk
!
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_ngfldr_pio"
!
      TYPE (My_VarDesc) :: TpioVar, VpioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Determine if inquiring about the field to process in input NetCDF
!  file(s).  This usually happens on first call or when the field
!  time records are split (saved) in several multi-files.
!
      Linquire=.FALSE.
      Lmulti=.FALSE.
      IF (iic(ng).eq.0) Linquire=.TRUE.
      IF (.not.Linquire.and.                                            &
     &    ((Iinfo(10,ifield,ng).gt.1).and.                              &
     &     (Linfo( 5,ifield,ng).or.                                     &
     &     (Finfo( 1,ifield,ng)*day2sec.gt.time(ng))))) THEN
        Linquire=.TRUE.
        Lmulti=.TRUE.
      END IF
!
!  If appropriate, inquire about the contents of input NetCDF file and
!  fill information arrays.
!
!  Also, if appropriate, deactivate the Linfo(5,ifield,ng) switch after
!  the query for the LOWER snapshot interpolant from previous multifile
!  in the list. The switch was activated previously to indicate the
!  processing of the FIRST record of the file for the UPPER snapshot
!  interpolant.
!
      IF (Linquire) THEN
        job=-1
        CALL inquiry (ng, model, job, UBk, Iend, UBi, ifield, pioFile,  &
     &                Lmulti, nfiles, S)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        IF (Linfo(5,ifield,ng)) THEN
          Linfo(5,ifield,ng)=.FALSE.
        END IF
        IF (Linfo(6,ifield,ng)) THEN
          Linfo(6,ifield,ng)=.FALSE.
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, read in new data.
!-----------------------------------------------------------------------
!
      update=.FALSE.
      Tmono=Finfo(7,ifield,ng)
!
      IF ((Tmono.gt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart(ng))) THEN
!
!  Load properties for requested field from information arrays.
!
        Liocycle=Linfo( 2,ifield,ng)
        Lonerec =Linfo( 3,ifield,ng)
        Vtype   =Iinfo( 1,ifield,ng)
        VpioVar =Dinfo( 1,ifield,ng)
        TpioVar =Dinfo( 2,ifield,ng)
        Nrec    =Iinfo( 4,ifield,ng)
        Tindex  =Iinfo( 8,ifield,ng)
        Trec    =Iinfo( 9,ifield,ng)
        nvdim   =Iinfo(11,ifield,ng)
        Tmin    =Finfo( 1,ifield,ng)
        Tmax    =Finfo( 2,ifield,ng)
        Clength =Finfo( 5,ifield,ng)
        Tscale  =Finfo( 6,ifield,ng)
        ncfile  =Cinfo(ifield,ng)
!
        IF (Liocycle) THEN
          Trec=MOD(Trec,Nrec)-1
          IF (Trec.le.0) Trec=Nrec+Trec
        ELSE
          Trec=Trec-1
        END IF
        Iinfo(9,ifield,ng)=Trec
!
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
          IF (.not.recordless.and.(TpioVar%vd%varID.ge.0).and.          &
     &        (TpioVar%vd%varID.ne.VpioVar%vd%varID)) THEN
            CALL pio_netcdf_get_time (ng, model, ncfile, Tname(ifield), &
     &                                Rclock%DateNumber, Tval,          &
     &                                pioFile = pioFile,                &
     &                                start = (/Trec/),                 &
     &                                total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,50) TRIM(Tname(ifield)), Trec
              RETURN
            END IF
            Tval=Tval*Tscale
            Vtime(Tindex,ifield,ng)=Tval
!
!  Activate switch Linfo(5,ifield,ng) if processing the FIRST record of
!  the file for the UPPER time snapshot. We need to get the LOWER time
!  snapshot from previous multifile.
!
            IF ((Trec.eq.1).and.(Tval*day2sec.ge.time(ng))) THEN
              Linfo(5,ifield,ng)=.TRUE.
            END IF
          END IF
!
!  Read in non-grided data. The conditional statement on Jrec is to
!  differentiate between reading a 3D and 2D array.
!
          IF (VpioVar%vd%varID.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF (nvdim.eq.1) THEN
              npts=Iend-Istr+1
              CALL pio_netcdf_get_fvar (ng, model, ncfile,              &
     &                                  Vname(1,ifield), Awrk,          &
     &                                  pioFile = pioFile,              &
     &                                  start = (/1/),                  &
     &                                  total = (/Iend-Istr+1/))
            ELSE IF (nvdim.eq.2) THEN
              IF (recordless) THEN
                npts=(Iend-Istr+1)*Jrec
                CALL pio_netcdf_get_fvar (ng, model, ncfile,            &
     &                                    Vname(1,ifield), Awrk,        &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1,1/),              &
     &                                    total = (/Iend-Istr+1,Jrec/))
              ELSE
                npts=Iend-Istr+1
                CALL pio_netcdf_get_fvar (ng, model, ncfile,            &
     &                                    Vname(1,ifield), Awrk,        &
     &                                    pioFile = pioFile,            &
     &                                    start = (/1,Trec/),           &
     &                                    total = (/Iend-Istr+1,1/))
              END IF
            ELSE IF (nvdim.eq.3) THEN
              npts=(Iend-Istr+1)*Jrec
              CALL pio_netcdf_get_fvar (ng, model, ncfile,              &
     &                                  Vname(1,ifield), Awrk,          &
     &                                  pioFile = pioFile,              &
     &                                  start = (/1,1,Trec/),           &
     &                                  total = (/Iend-Istr+1,Jrec,1/))
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,50) TRIM(Vname(1,ifield)), Trec
              RETURN
            END IF
            Fmin=Awrk(1)*Fscale(ifield,ng)
            Fmax=Awrk(1)*Fscale(ifield,ng)
            ic=0
            DO j=1,Jrec
              DO i=Istr,Iend
                ic=ic+1
                Aval=Awrk(ic)*Fscale(ifield,ng)
                Fmin=MIN(Fmin,Aval)
                Fmax=MAX(Fmax,Aval)
                Fout(i,j,Tindex)=Aval
              END DO
            END DO
            Finfo(8,ifield,ng)=Fmin
            Finfo(9,ifield,ng)=Fmax
#  ifdef CHECKSUM
            CALL get_hash (Awrk, npts, Fhash)
#  endif
            IF (Master) THEN
              IF (UBk.eq.1) THEN
                WRITE (stdout,60) TRIM(Vname(2,ifield)), ng, Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                Tsec=Tval*day2sec
                CALL time_string (Tsec, t_code)
                WRITE (stdout,70) TRIM(Vname(2,ifield)), t_code,        &
     &                            ng, Trec, Tindex, ncfile(lstr:lend),  &
     &                            Tmin, Tmax, Tval, Fmin, Fmax
              END IF
#  ifdef CHECKSUM
              WRITE (stdout,80) Fhash
#  endif
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
        IF (.not.Lonerec.and.(.not.recordless)) THEN
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
  10  FORMAT (/,' GET_NGFLDR_PIO - unable to find dimension ',a,        &
     &        /,18x,'for variable: ',a,/,18x,'in file: ',a,             &
     &        /,18x,'file is not CF compliant...')
  20  FORMAT (/,' GET_NGFLDR_PIO - too small dimension for variable',   &
     &        1x,a,': ',i0,2x,i0)
  30  FORMAT (/,' GET_NGFLDR_PIO - unable to find requested',           &
     &        ' variable: ',a,/,18x,'in file: ',a)
  40  FORMAT (/,' GET_NGFLDR_PIO - unable to open input NetCDF',        &
     &        ' file: ',a)
  50  FORMAT (/,' GET_NGFLDR_PIO - error while reading variable: ',a,   &
     &        2x,' at TIME index = ',i0)
  60  FORMAT (2x,'GET_NGFLDR_PIO   - ',a,/,22x,'(Grid = ',i2.2,         &
     &        ', Min = ',1pe15.8,' Max = ',1pe15.8,')')
  70  FORMAT (2x,'GET_NGFLDR_PIO   - ',a,',',t75,a,/,22x,               &
     &        '(Grid=',i2.2,', Rec=',i0,', Index=',i1,                  &
     &        ', File: ',a,')',/,22x,                                   &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',                   &
     &        t71, 't = ', f15.4 ,/,22x,                                &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
#  ifdef CHECKSUM
  80  FORMAT (22x,'(CheckSum = ',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE get_ngfldr_pio
# endif
#else
      SUBROUTINE get_ngfldr
      RETURN
      END SUBROUTINE get_ngfldr
#endif
