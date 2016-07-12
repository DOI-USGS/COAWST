#include "cppdefs.h"
#if defined ADJOINT && defined SOLVE3D
      SUBROUTINE get_3dfldr (ng, model, ifield, ncid, nfiles, S,        &
     &                       update, LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                       Iout, Irec,                                &
# ifdef MASKING
     &                       Fmask,                                     &
# endif
     &                       Fout)
!
!svn $Id: get_3dfldr.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested 3D field (point or grided) from     !
!  specified NetCDF file. Backward time processing.                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ifield     Field ID.                                             !
!     ncid       NetCDF file ID.                                       !
!     nfiles     Number of input NetCDF files.                         !
!     S          I/O derived type structure, TYPE(T_IO).               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Iout       Size of the outer dimension,  if any.  Otherwise,     !
!                  Iout must be set to one by the calling program.     !
!     Irec       Number of 3D field records to read.                   !
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
      USE mod_netcdf
      USE mod_scalars
!
      USE nf_fread3d_mod, ONLY : nf_fread3d
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, model, ifield, nfiles
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Iout, Irec
      integer, intent(inout) :: ncid

      TYPE(T_IO), intent(inout) :: S(nfiles)

# ifdef MASKING
      real(r8), intent(in) :: Fmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: Fout(LBi:UBi,LBj:UBj,LBk:UBk,Iout)
!
!  Local variable declarations.
!
      logical :: Lgridded, Linquire, Liocycle, Lmulti, Lonerec

      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, job, lend, lstr, lvar, status
      integer :: Vsize(4)

      real(r8) :: Clength, Fmax, Fmin, Tdelta, Tend
      real(r8) :: Tmax, Tmin, Tmono, Tscale, Tstr
      real(r8) :: Fval, Tval

      character (len=14) :: t_code
!
      SourceFile='get_3dfldr.F'
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
     &     (Finfo( 1,ifield,ng)*day2sec.gt.time(ng)))) THEN
        Linquire=.TRUE.
        Lmulti=.TRUE.
      END IF
!
!  On first call, inquire about the contents of input NetCDF file and
!  fill information arrays.
!
      IF (Linquire) THEN
        job=-3
        CALL inquire (ng, model, job, Iout, Irec, 1, ifield, ncid,      &
     &                Lmulti, nfiles, S)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Load properties for requested field from information arrays.
!
      Lgridded=Linfo(1,ifield,ng)
      Liocycle=Linfo(2,ifield,ng)
      Lonerec =Linfo(3,ifield,ng)
      Vtype   =Iinfo(1,ifield,ng)
      Vid     =Iinfo(2,ifield,ng)
      Tid     =Iinfo(3,ifield,ng)
      Nrec    =Iinfo(4,ifield,ng)
      Vsize(1)=Iinfo(5,ifield,ng)
      Vsize(2)=Iinfo(6,ifield,ng)
      Vsize(3)=Iinfo(7,ifield,ng)
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
!  "Iout" is unity, input data is stored in timeless array by the
!  calling program.  If Irec > 1, this routine is used to read a 3D
!  field varying in another non-time dimension.
!
          IF (Irec.eq.1) THEN
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
          IF (Tid.ge.0) THEN
            CALL netcdf_get_fvar (ng, model, ncfile, Tname(ifield),     &
     &                            Tval,                                 &
     &                            ncid = ncid,                          &
     &                            start = (/Trec/),                     &
     &                            total = (/1/))
            IF (exit_flag.ne.NoError) THEN
              IF (Master) WRITE (stdout,40) TRIM(Tname(ifield)), Trec
              RETURN
            END IF
          END IF
          Tval=Tval*Tscale
          Vtime(Tindex,ifield,ng)=Tval
!
!  Read in 3D-grided or point data.
!
          IF (Vid.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF (Lgridded) THEN
              IF (Irec.gt.1) THEN
                Finfo(8,ifield,ng)=1.0E+35_r8
                Finfo(9,ifield,ng)=-1.0E+35_r8
                DO i=1,Irec
                  status=nf_fread3d(ng, model, ncfile, ncid,            &
     &                              Vname(1,ifield), Vid,               &
     &                              i, Vtype, Vsize,                    &
     &                              LBi, UBi, LBj, UBj, LBk, UBk,       &
     &                              Fscale(ifield,ng), Fmin, Fmax,      &
# ifdef MASKING
     &                              Fmask,                              &
# endif
     &                              Fout(:,:,:,i))
                  Finfo(8,ifield,ng)=MIN(Fmin,Finfo(8,ifield,ng))
                  Finfo(9,ifield,ng)=MAX(Fmax,Finfo(9,ifield,ng))
                END DO
              ELSE
                status=nf_fread3d(ng, model, ncfile, ncid,              &
     &                            Vname(1,ifield), Vid,                 &
     &                            Trec, Vtype, Vsize,                   &
     &                            LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                            Fscale(ifield,ng), Fmin, Fmax,        &
# ifdef MASKING
     &                            Fmask,                                &
# endif
     &                            Fout(:,:,:,Tindex))
                Finfo(8,ifield,ng)=Fmin
                Finfo(9,ifield,ng)=Fmax
              END IF
            ELSE
              CALL netcdf_get_fvar (ng, model, ncfile, Vname(1,ifield), &
     &                              Fval,                               &
     &                              ncid = ncid,                        &
     &                              start = (/Trec/),                   &
     &                              total = (/1/))
              Fval=Fval*Fscale(ifield,ng)
              Fpoint(Tindex,ifield,ng)=Fval
              Fmin=Fval
              Fmax=Fval
            END IF
            IF (exit_flag.ne.NoError) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,ifield)), Trec
              END IF
              RETURN
            END IF
            IF (Master) THEN
              IF (Irec.gt.1) THEN
                WRITE (stdout,50) TRIM(Vname(2,ifield)), Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                lvar=MIN(43,LEN_TRIM(Vname(2,ifield)))
                Tval=Tval*day2sec
                CALL time_string (Tval, t_code)
                WRITE (stdout,60) Vname(2,ifield)(1:lvar), t_code,      &
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
  10  FORMAT (/,' GET_3DFLDR  - unable to find dimension ',a,           &
     &        /,15x,'for variable: ',a,/,15x,'in file: ',a,             &
     &        /,15x,'file is not CF compliant...')
  20  FORMAT (/,' GET_3DFLDR  - unable to find requested variable: ',a, &
     &        /,15x,'in input NetCDF file: ',a)
  30  FORMAT (/,' GET_3DFLDR  - unable to open input NetCDF file: ',a)
  40  FORMAT (/,' GET_3DFLDR  - error while reading variable: ',a,2x,   &
     &       ' at TIME index = ',i7)
  50  FORMAT (3x,' GET_3DFLDR  - ',a,/,19x,'(Min = ',1p,e15.8,0p,       &
     &        ' Max = ',1p,e15.8,0p,')')
  60  FORMAT (3x,' GET_3DFLDR  - ',a,',',t62,'t = ',a,/,19x,            &
     &        '(Rec=',i7.7,', Index=',i1,', File: ',a,')',/,19x,        &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',/, 19x,            &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE get_3dfldr
#else
      SUBROUTINE get_3dfldr
      RETURN
      END SUBROUTINE get_3dfldr
#endif
