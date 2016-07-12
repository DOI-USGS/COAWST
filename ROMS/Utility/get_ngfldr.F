#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE get_ngfldr (ng, model, ifield, ncid, nfiles, S,        &
     &                       update, LBi, UBi, UBj, UBk,                &
     &                       Istr, Iend, Jrec,                          &
     &                       Fout)
!
!svn $Id: get_ngfldr.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      integer, intent(in) :: LBi, UBi, UBj, UBk, Istr, Iend, Jrec
      integer, intent(inout) :: ncid

      TYPE(T_IO), intent(inout) :: S(nfiles)

      real(r8), intent(inout) :: Fout(LBi:UBi,UBj,UBk)
!
!  Local variable declarations.
!
      logical :: Linquire, Liocycle, Lmulti, Lonerec

      integer :: Nrec, Tid, Tindex, Trec, Vid, Vtype
      integer :: i, ic, j, job, lend, lstr, npts, status

      real(r8) :: Aval, Clength, Fmax, Fmin, Tdelta, Tend
      real(r8) :: Tmax, Tmin, Tmono, Tscale, Tstr
      real(r8) :: Tval

      real(r8), dimension((UBi-LBi+1)*UBj) :: A

      character (len=14) :: t_code
!
      SourceFile='get_ngfldr.F'
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
!  If appropriate, inquire about the contents of input NetCDF file and
!  fill information arrays.
!
      IF (Linquire) THEN
        job=-1
        CALL inquire (ng, model, job, UBk, Iend, UBi, ifield, ncid,     &
     &                Lmulti, nfiles, S)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Load properties for requested field from information arrays.
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
     &          ' at TIME index = ',i7)
  60  FORMAT (3x,' GET_NGFLDR  - ',a,/,19x,'(Min = ',1pe15.8,' Max = ', &
     &        1pe15.8,')')
  70  FORMAT (3x,' GET_NGFLDR  - ',a,',',t62,'t = ',a,/,19x,            &
     &        '(Rec=',i7.7,', Index=',i1,', File: ',a,')',/,19x,        &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',/, 19x,            &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')

      RETURN
      END SUBROUTINE get_ngfldr
#else
      SUBROUTINE get_ngfldr
      RETURN
      END SUBROUTINE get_ngfldr
#endif
