#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE get_2dfldr (ng, model, ifield, ncid,                   &
# if defined PIO_LIB && defined DISTRIBUTE
     &                       pioFile,                                   &
# endif
     &                       nfiles, S, update,                         &
     &                       LBi, UBi, LBj, UBj, Iout, Irec,            &
# ifdef MASKING
     &                       Fmask,                                     &
# endif
     &                       Fout)
!
!git $Id$
!svn $Id: get_2dfldr.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested 2D field (point or grided) from     !
!  specified NetCDF file. Backward time processing.                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ifield     Field ID.                                             !
!     ncid       NetCDF file ID.                                       !
# if defined PIO_LIB && defined DISTRIBUTE
!     pioFile    PIO file descriptor structure, TYPE(file_desc_t)      !
!                  pioFile%fh         file handler                     !
!                  pioFile%iosystem   IO system descriptor (struct)    !
# endif
!     nfiles     Number of input NetCDF files.                         !
!     S          I/O derived type structure, TYPE(T_IO).               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     Iout       Size of the outer dimension,  if any.  Otherwise,     !
!                  Iout must be set to one by the calling program.     !
!     Irec       Number of 2D field records to read.                   !
!     Fmask      Land/Sea mask, if any.                                !
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
      logical, intent(out) :: update
!
      integer, intent(in) :: ng, model, ifield, nfiles, Iout, Irec
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(inout) :: ncid
!
# if defined PIO_LIB && defined DISTRIBUTE
      TYPE (File_desc_t), intent(inout) :: pioFile
# endif
      TYPE(T_IO), intent(inout) :: S(nfiles)
!
# ifdef MASKING
      real(r8), intent(in) :: Fmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: Fout(LBi:UBi,LBj:UBj,Iout)
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Read in requested 2D field for backward time processing according
!  to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (S(1)%IOtype)
        CASE (io_nf90)
          CALL get_2dfldr_nf90 (ng, model, ifield, ncid,                &
     &                          nfiles, S, update,                      &
     &                          LBi, UBi, LBj, UBj, Iout, Irec,         &
# ifdef MASKING
     &                          Fmask,                                  &
# endif
     &                          Fout)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_2dfldr_pio (ng, model, ifield, pioFile,              &
     &                         nfiles, S, update,                       &
     &                         LBi, UBi, LBj, UBj, Iout, Irec,          &
#  ifdef MASKING
     &                         Fmask,                                   &
#  endif
     &                         Fout)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) S(1)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_2DFLDR - Illegal input file type, io_type = ',i0,   &
     &        /,14x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_2dfldr
!
!***********************************************************************
      SUBROUTINE get_2dfldr_nf90 (ng, model, ifield, ncid,              &
     &                            nfiles, S, update,                    &
     &                            LBi, UBi, LBj, UBj, Iout, Irec,       &
# ifdef MASKING
     &                            Fmask,                                &
# endif
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
      USE dateclock_mod,  ONLY : time_string
      USE inquiry_mod,    ONLY : inquiry
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE nf_fread3d_mod, ONLY : nf_fread3d
      USE strings_mod,    ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update
!
      integer, intent(in) :: ng, model, ifield, nfiles, Iout, Irec
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(inout) :: ncid
!
      TYPE(T_IO), intent(inout) :: S(nfiles)
!
# ifdef MASKING
      real(r8), intent(in) :: Fmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: Fout(LBi:UBi,LBj:UBj,Iout)
!
!  Local variable declarations.
!
      logical :: Lgridded, Linquire, Liocycle, Lmulti, Lonerec, Lregrid
      logical :: special
!
      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: gtype, job, lend, lstr, lvar, status
      integer :: Vsize(4)
# ifdef CHECKSUM
      integer(i8b) :: Fhash
# endif
!
      real(r8) :: Fmax, Fmin, Fval

      real(dp) :: Clength, Tdelta, Tend
      real(dp) :: Tmax, Tmin, Tmono, Tscale, Tstr
      real(dp) :: Tsec, Tval
!
      character (len= 1) :: Rswitch
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
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
        job=-2
        CALL inquiry (ng, model, job, Iout, Irec, 1, ifield, ncid,      &
     &                Lmulti, nfiles, S)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        IF (Linfo(5,ifield,ng)) THEN
          Linfo(5,ifield,ng)=.FALSE.
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, read in new data.
!-----------------------------------------------------------------------
!
      update=.FALSE.
      Lregrid=.FALSE.
      Tmono=Finfo(7,ifield,ng)
!
      IF ((Tmono.gt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart(ng))) THEN
!
!  Load properties for requested field from information arrays.
!
        Lgridded=Linfo(1,ifield,ng)
        Liocycle=Linfo(2,ifield,ng)
        Lonerec =Linfo(3,ifield,ng)
        special =Linfo(4,ifield,ng)
        Vtype   =Iinfo(1,ifield,ng)
        Vid     =Iinfo(2,ifield,ng)
        Tid     =Iinfo(3,ifield,ng)
        Nrec    =Iinfo(4,ifield,ng)
        Vsize(1)=Iinfo(5,ifield,ng)
        Vsize(2)=Iinfo(6,ifield,ng)
        Tindex  =Iinfo(8,ifield,ng)
        Trec    =Iinfo(9,ifield,ng)
        Tmin    =Finfo(1,ifield,ng)
        Tmax    =Finfo(2,ifield,ng)
        Clength =Finfo(5,ifield,ng)
        Tscale  =Finfo(6,ifield,ng)
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
!  "Iout" is unity, input data is stored in timeless array by the
!  calling program.  If Irec > 1, this routine is used to read a 2D
!  field varying in another non-time dimension like tidal constituent
!  component.
!
          IF (.not.special.and.(Irec.eq.1)) THEN
            IF (Iout.eq.1) THEN
              Tindex=1
            ELSE
              Tindex=3-Tindex
            END IF
            Iinfo(8,ifield,ng)=Tindex
          END IF
!
!  Read in time coordinate and scale it to day units.
!
          IF (.not.special.and.(Tid.ge.0)) THEN
            CALL netcdf_get_time (ng, model, ncfile, Tname(ifield),     &
     &                            Rclock%DateNumber, Tval,              &
     &                            ncid = ncid,                          &
     &                            start = (/Trec/),                     &
     &                            total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,40) TRIM(Tname(ifield)), Trec
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
!  Read in 2D-grided or point data. Notice for special 2D fields, Vtype
!  is augmented by four indicating reading a 3D field. This rational is
!  used to read fields like tide data.
!
          IF (Vid.ge.0) THEN
            IF (Lgridded) THEN
              IF (special) THEN
                Vsize(3)=Irec
                gtype=Vtype+4
                status=nf_fread3d(ng, model, ncfile, ncid,              &
     &                            Vname(1,ifield), Vid,                 &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj, 1, Irec,          &
     &                            Fscale(ifield,ng), Fmin, Fmax,        &
# ifdef MASKING
     &                            Fmask,                                &
# endif
# ifdef CHECKSUM
     &                            Fout,                                 &
     &                            checksum = Fhash)
# else
     &                            Fout)
# endif
              ELSE
                status=nf_fread2d(ng, model, ncfile, ncid,              &
     &                            Vname(1,ifield), Vid,                 &
     &                            Trec, Vtype, Vsize,                   &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscale(ifield,ng), Fmin, Fmax,        &
# ifdef MASKING
     &                            Fmask,                                &
# endif
     &                            Fout(:,:,Tindex),                     &
# ifdef CHECKSUM
     &                            checksum = Fhash,                     &
# endif
     &                            Lregrid = Lregrid)
              END IF
            ELSE
              CALL netcdf_get_fvar (ng, model, ncfile,                  &
     &                              Vname(1,ifield), Fval,              &
     &                              ncid = ncid,                        &
     &                              start = (/Trec/),                   &
     &                              total = (/1/))
              Fval=Fval*Fscale(ifield,ng)
              Fpoint(Tindex,ifield,ng)=Fval
              Fmin=Fval
              Fmax=Fval
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,40) TRIM(Vname(1,ifield)), Trec
              RETURN
            END IF
            Finfo(8,ifield,ng)=Fmin
            Finfo(9,ifield,ng)=Fmax
            IF (Master) THEN
              IF (special) THEN
                WRITE (stdout,50) TRIM(Vname(2,ifield)), ng, Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                lvar=MIN(46,LEN_TRIM(Vname(2,ifield)))
                Tsec=Tval*day2sec
                CALL time_string (Tsec, t_code)
                IF (Lregrid) THEN
                  Rswitch='T'
                ELSE
                  Rswitch='F'
                END IF
                WRITE (stdout,60) Vname(2,ifield)(1:lvar), t_code,      &
     &                            ng, Trec, Tindex, ncfile(lstr:lend),  &
     &                            Tmin, Tmax, Tval, Fmin, Fmax, Rswitch
              END IF
# ifdef CHECKSUM
              WRITE (stdout,70) Fhash
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
        IF (.not.Lonerec.and.(.not.special)) THEN
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
  10  FORMAT (/,' GET_2DFLDR_NF90 - unable to find dimension ',a,       &
     &        /,19x,'for variable: ',a,/,19x,'in file: ',a,             &
     &        /,19x,'file is not CF compliant...')
  20  FORMAT (/,' GET_2DFLDR_NF90 - unable to find requested',          &
     &        ' variable:',a,/,19x,'in ',a)
  30  FORMAT (/,' GET_2DFLDR_NF90 - unable to open input NetCDF',       &
     &        ' file: ',a)
  40  FORMAT (/,' GET_2DFLDR_NF90 - error while reading variable: ',a,  &
     &        2x,' at TIME index = ',i0)
  50  FORMAT (2x,'GET_2DFLDR_NF90  - ',a,/,22x,'(Grid = ',              &
     &        ', Min = ',1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  60  FORMAT (2x,'GET_2DFLDR_NF90  - ',a,',',t75,a,/,22x,               &
     &        '(Grid=',i2.2,', Rec=',i0,', Index=',i1,                  &
     &        ', File: ',a,')',/,22x,                                   &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',                   &
     &        t71, 't = ', f15.4 ,/,22x,                                &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')',         &
     &        t71, 'regrid = ',a)
# ifdef CHECKSUM
  70  FORMAT (22x,'(CheckSum = ',i0,')')
# endif
!
      RETURN
      END SUBROUTINE get_2dfldr_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_2dfldr_pio (ng, model, ifield, pioFile,            &
     &                           nfiles, S, update,                     &
     &                           LBi, UBi, LBj, UBj, Iout, Irec,        &
#  ifdef MASKING
     &                           Fmask,                                 &
#  endif
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
      USE dateclock_mod,  ONLY : time_string
      USE inquiry_mod,    ONLY : inquiry
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE nf_fread3d_mod, ONLY : nf_fread3d
      USE strings_mod,    ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update
!
      integer, intent(in) :: ng, model, ifield, nfiles, Iout, Irec
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE (File_desc_t), intent(inout) :: pioFile
      TYPE(T_IO), intent(inout) :: S(nfiles)
!
#  ifdef MASKING
      real(r8), intent(in) :: Fmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(inout) :: Fout(LBi:UBi,LBj:UBj,Iout)
!
!  Local variable declarations.
!
      logical :: Lgridded, Linquire, Liocycle, Lmulti, Lonerec, Lregrid
      logical :: special
!
      integer :: Nrec, Tindex, Trec, Vtype
      integer :: i, job, lend, lstr, lvar, status
      integer :: Vsize(4)
#  ifdef CHECKSUM
      integer(i8b) :: Fhash
#  endif
!
      real(r8) :: Fmax, Fmin, Fval

      real(dp) :: Clength, Tdelta, Tend
      real(dp) :: Tmax, Tmin, Tmono, Tscale, Tstr
      real(dp) :: Tsec, Tval
!
      character (len= 1) :: Rswitch
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_3dfldr_pio"
!
      TYPE (IO_Desc_t), pointer :: ioDesc
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
        job=-2
        CALL inquiry (ng, model, job, Iout, Irec, 1, ifield, pioFile,   &
     &                Lmulti, nfiles, S)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        IF (Linfo(5,ifield,ng)) THEN
          Linfo(5,ifield,ng)=.FALSE.
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  If appropriate, read in new data.
!-----------------------------------------------------------------------
!
      update=.FALSE.
      Lregrid=.FALSE.
      Tmono=Finfo(7,ifield,ng)
!
      IF ((Tmono.gt.time(ng)).or.(iic(ng).eq.0).or.                     &
     &    (iic(ng).eq.ntstart(ng))) THEN
!
!  Load properties for requested field from information arrays.
!
        Lgridded=Linfo(1,ifield,ng)
        Liocycle=Linfo(2,ifield,ng)
        Lonerec =Linfo(3,ifield,ng)
        special =Linfo(4,ifield,ng)
        Vtype   =Iinfo(1,ifield,ng)
        VpioVar =Dinfo(1,ifield,ng)
        TpioVar =Dinfo(2,ifield,ng)
        Nrec    =Iinfo(4,ifield,ng)
        Vsize(1)=Iinfo(5,ifield,ng)
        Vsize(2)=Iinfo(6,ifield,ng)
        Tindex  =Iinfo(8,ifield,ng)
        Trec    =Iinfo(9,ifield,ng)
        Tmin    =Finfo(1,ifield,ng)
        Tmax    =Finfo(2,ifield,ng)
        Clength =Finfo(5,ifield,ng)
        Tscale  =Finfo(6,ifield,ng)
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
!  "Iout" is unity, input data is stored in timeless array by the
!  calling program.  If Irec > 1, this routine is used to read a 2D
!  field varying in another non-time dimension like tidal constituent
!  component.
!
          IF (.not.special.and.(Irec.eq.1)) THEN
            IF (Iout.eq.1) THEN
              Tindex=1
            ELSE
              Tindex=3-Tindex
            END IF
            Iinfo(8,ifield,ng)=Tindex
          END IF
!
!  Read in time coordinate and scale it to day units.
!
          IF (.not.special.and.(TpioVar%vd%varID.ge.0)) THEN
            CALL pio_netcdf_get_time (ng, model, ncfile, Tname(ifield), &
     &                                Rclock%DateNumber, Tval,          &
     &                                pioFile = pioFile,                &
     &                                start = (/Trec/),                 &
     &                                total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,40) TRIM(Tname(ifield)), Trec
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
!  Read in 2D-grided or point data. Notice for special 2D fields, Vtype
!  is augmented by four indicating reading a 3D field (tidal data with
!  tidal constituents in the third dimension. Currently, all input tidal
!  constituents data is at RHO-points.
!
          IF (special) Vtype=Vtype+4
!
          IF (VpioVar%vd%varID.ge.0) THEN
            SELECT CASE (Vtype)
              CASE (p2dvar)
                IF (KIND(Fout).eq.8) THEN
                  ioDesc => ioDesc_dp_p2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_p2dvar(ng)
                END IF
              CASE (r2dvar)
                IF (KIND(Fout).eq.8) THEN
                  ioDesc => ioDesc_dp_r2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_r2dvar(ng)
                END IF
              CASE (u2dvar)
                IF (KIND(Fout).eq.8) THEN
                  ioDesc => ioDesc_dp_u2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_u2dvar(ng)
                END IF
              CASE (v2dvar)
                IF (KIND(Fout).eq.8) THEN
                  ioDesc => ioDesc_dp_v2dvar(ng)
                ELSE
                  ioDesc => ioDesc_sp_v2dvar(ng)
                END IF
#  if defined SSH_TIDES || defined UV_TIDES

              CASE (r2dvar+4)
                IF (KIND(Fout).eq.8) THEN
                  ioDesc => ioDesc_dp_rtides(ng)
                ELSE
                  ioDesc => ioDesc_sp_rtides(ng)
                END IF
                VpioVar%gtype=r2dvar+4
#  endif
            END SELECT
!
            IF (Lgridded) THEN
              Fmin=spval
              Fmax=-spval
              IF (special) THEN
                Vsize(3)=Irec
                status=nf_fread3d(ng, model, ncfile, pioFile,           &
     &                            Vname(1,ifield), VpioVar,             &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj, 1, Irec,          &
     &                            Fscale(ifield,ng), Fmin, Fmax,        &
#  ifdef MASKING
     &                            Fmask,                                &
#  endif
#  ifdef CHECKSUM
     &                            Fout,                                 &
     &                            checksum = Fhash)
#  else
     &                            Fout)
#  endif
              ELSE
                status=nf_fread2d(ng, model, ncfile, pioFile,           &
     &                            Vname(1,ifield), VpioVar,             &
     &                            Trec, ioDesc, Vsize,                  &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            Fscale(ifield,ng), Fmin, Fmax,        &
#  ifdef MASKING
     &                            Fmask,                                &
#  endif
     &                            Fout(:,:,Tindex),                     &
#  ifdef CHECKSUM
     &                            checksum = Fhash,                     &
#  endif
     &                            Lregrid = Lregrid)
              END IF
            ELSE
              CALL pio_netcdf_get_fvar (ng, model, ncfile,              &
     &                                  Vname(1,ifield), Fval,          &
     &                                  pioFile = pioFile,              &
     &                                  start = (/Trec/),               &
     &                                  total = (/1/))
              Fval=Fval*Fscale(ifield,ng)
              Fpoint(Tindex,ifield,ng)=Fval
              Fmin=Fval
              Fmax=Fval
            END IF
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,40) TRIM(Vname(1,ifield)), Trec
              RETURN
            END IF
            Finfo(8,ifield,ng)=Fmin
            Finfo(9,ifield,ng)=Fmax
            IF (Master) THEN
              IF (special) THEN
                WRITE (stdout,50) TRIM(Vname(2,ifield)), ng, Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                lvar=MIN(46,LEN_TRIM(Vname(2,ifield)))
                Tsec=Tval*day2sec
                CALL time_string (Tsec, t_code)
                IF (Lregrid) THEN
                  Rswitch='T'
                ELSE
                  Rswitch='F'
                END IF
                WRITE (stdout,60) Vname(2,ifield)(1:lvar), t_code,      &
     &                            ng, Trec, Tindex, ncfile(lstr:lend),  &
     &                            Tmin, Tmax, Tval, Fmin, Fmax, Rswitch
              END IF
#  ifdef CHECKSUM
              WRITE (stdout,70) Fhash
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
        IF (.not.Lonerec.and.(.not.special)) THEN
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
  10  FORMAT (/,' GET_2DFLDR_PIO - unable to find dimension ',a,        &
     &        /,18x,'for variable: ',a,/,18x,'in file: ',a,             &
     &        /,18x,'file is not CF compliant...')
  20  FORMAT (/,' GET_2DFLDR_PIO - unable to find requested',           &
     &        ' variable:',a,/,18x,'in ',a)
  30  FORMAT (/,' GET_2DFLDR_PIO - unable to open input NetCDF',        &
     &        ' file: ',a)
  40  FORMAT (/,' GET_2DFLDR_PIO - error while reading variable: ',a,   &
     &        2x,' at TIME index = ',i0)
  50  FORMAT (2x,'GET_2DFLDR_PIO   - ',a,/,22x,'(Grid = ',              &
     &        ', Min = ',1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
  60  FORMAT (2x,'GET_2DFLDR_PIO   - ',a,',',t75,a,/,22x,               &
     &        '(Grid=',i2.2,', Rec=',i0,', Index=',i1,                  &
     &        ', File: ',a,')',/,22x,                                   &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',                   &
     &        t71, 't = ', f15.4 ,/,22x,                                &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')',         &
     &        t71, 'regrid = ',a)
#  ifdef CHECKSUM
  70  FORMAT (22x,'(CheckSum = ',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE get_2dfldr_pio
# endif
#else
      SUBROUTINE get_2dfldr
      RETURN
      END SUBROUTINE get_2dfldr
#endif
