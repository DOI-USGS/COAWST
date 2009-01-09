#include "cppdefs.h"
#if defined ADJOINT && defined SOLVE3D
      SUBROUTINE get_3dfldr (ng, model, ifield, ncid, nfiles, fname,    &
     &                       update, LBi, UBi, LBj, UBj, LBk, UBk,      &
     &                       Iout, Irec,                                &
# ifdef MASKING
     &                       Fmask,                                     &
# endif
     &                       Fout)
!
!svn $Id: get_3dfldr.F 1598 2008-07-25 03:29:09Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!     fname      NetCDF file name(s).                                  !
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

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
      USE distribute_mod, ONLY : mp_bcastl, mp_bcasts
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(out) :: update

      integer, intent(in) :: ng, model, ifield, nfiles, Iout, Irec
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(inout) :: ncid

      character (len=*), intent(in) :: fname(nfiles)

# ifdef MASKING
      real(r8), intent(in) :: Fmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: Fout(LBi:UBi,LBj:UBj,LBk:UBk,Iout)
!
!  Local variable declarations.
!
      logical :: Liocycle, Lgrided, Lonerec, got_var, got_time

      integer :: Tid, Tindex, Trec, Vid, Vtype
      integer :: i, lend, lfvar, lstr, ltvar, lvar, nrec, nvatts
      integer :: nvd, status, vartype
      integer :: Vsize(4), start(4), total(4)

      integer :: nf_fread3d

      real(r8) :: Clength, Fmax, Fmin, Tdelta, Tend
      real(r8) :: Tmax, Tmin, Tmono, Tscale, Tstr, scale
      real(r8) :: Fval(1), Tval(1)

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
          CALL opencdf (ng, nfiles, fname, ncfile, N(ng), ifield, nrec, &
     &                  nvd, Vsize)
          IF (exit_flag.ne.NoError) RETURN
!
!  Scan variable list from input NetCDF and check for requested
!  variables.
!
          DO i=1,nvars
            IF (lfvar.gt.0) THEN
              IF (TRIM(varnam(i)).eq.TRIM(Vname(1,ifield))) THEN
                Vid=i
                got_var=.TRUE.
                IF (nvdims(i).gt.1) Lgrided=.TRUE.
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
          Iinfo(4,ifield,ng)=nrec
          Iinfo(5,ifield,ng)=Vsize(1)
          Iinfo(6,ifield,ng)=Vsize(2)
          Iinfo(7,ifield,ng)=Vsize(3)
!
!  Terminate execution requested variables are not found.
!
          IF (.not.got_var.and.(lfvar.gt.0)) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), TRIM(ncfile)
            exit_flag=2
            RETURN
          END IF
          IF (.not.got_time.and.(ltvar.gt.0)) THEN
            WRITE (stdout,10) TRIM(Tname(ifield)), TRIM(ncfile)
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
          CALL get_cycle (ng, ncid, Tid, nrec, tdays(ng), Liocycle,     &
     &                    Clength, Trec, Tstr, Tend, Tmin, Tmax, Tscale)
          IF (exit_flag.ne.NoError) RETURN
          Linfo(2,ifield,ng)=Liocycle
          Finfo(1,ifield,ng)=Tmin
          Finfo(2,ifield,ng)=Tmax
          Finfo(3,ifield,ng)=Tstr
          Finfo(4,ifield,ng)=Tend
          Finfo(5,ifield,ng)=Clength
          Finfo(6,ifield,ng)=Tscale
!
!  The strategy here is to create a local, monotonically decreasing
!  time variable so the interpolation between snapshots is trivial
!  when cycling forcing fields.
!
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
            start(1)=Trec-1
            total(1)=1
            status=nf90_get_var(ncid, Tid, Tval, start, total)
            Tval(1)=Tval(1)*Tscale
            IF (Tval(1).lt.Tend) THEN
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
          IF (nrec.eq.1) Lonerec=.TRUE.
          Linfo(3,ifield,ng)=Lonerec
        END IF
# ifdef DISTRIBUTE
!
!  Broadcast variable information to all nodes in the group.
!
        CALL mp_bcastl (ng, model, Linfo(:,ifield,ng), 3)
        CALL mp_bcasti (ng, model, Iinfo(:,ifield,ng), 9)
        CALL mp_bcastf (ng, model, Finfo(:,ifield,ng), 7)
        CALL mp_bcasts (ng, model, Cinfo(ifield,ng), 80)
# endif
        Tindex=Iinfo(8,ifield,ng)
        Vtime(Tindex,ifield,ng)=Finfo(4,ifield,ng)
      END IF         
!
!-----------------------------------------------------------------------
!  Get requested field information from global storage.
!-----------------------------------------------------------------------
!
      Lgrided =Linfo(1,ifield,ng)
      Liocycle=Linfo(2,ifield,ng)
      Lonerec =Linfo(3,ifield,ng)
      Vtype   =Iinfo(1,ifield,ng)
      Vid     =Iinfo(2,ifield,ng)
      Tid     =Iinfo(3,ifield,ng)
      nrec    =Iinfo(4,ifield,ng)
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
          Trec=MOD(Trec,nrec)-1
          IF (Trec.le.0) Trec=nrec-Trec
        ELSE
          Trec=Trec-1
        END IF
        Iinfo(9,ifield,ng)=Trec
        IF ((1.le.Trec).and.(Trec.le.nrec)) THEN
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
          IF (InpThread.and.(Tid.ge.0)) THEN
            start(1)=Trec
            total(1)=1
            status=nf90_get_var(ncid, Tid, Tval, start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,60) TRIM(Tname(ifield)), Trec
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END IF
# ifdef DISTRIBUTE
          CALL mp_bcasti (ng, model, exit_flag, 1)
          CALL mp_bcastf (ng, model, Tval, 1)
# endif
          Tval(1)=Tval(1)*Tscale
          Vtime(Tindex,ifield,ng)=Tval(1)
!
!  Read in 3D-grided or point data.
!
          IF (Vid.ge.0) THEN
            Fmin=0.0_r8
            Fmax=0.0_r8
            IF (Lgrided) THEN
              IF (Irec.gt.1) THEN
                Finfo(8,ifield,ng)=1.0E+35_r8
                Finfo(9,ifield,ng)=-1.0E+35_r8
                DO i=1,Irec
                  status=nf_fread3d(ng, model, ncid, Vid, i, Vtype,     &
     &                              Vsize, LBi, UBi, LBj, UBj, LBk, UBk,&
     &                              Fscale(ifield,ng), Fmin, Fmax,      &
# ifdef MASKING
     &                              Fmask(LBi,LBj),                     &
# endif
     &                              Fout(LBi,LBj,LBk,i))
                  Finfo(8,ifield,ng)=MIN(Fmin,Finfo(8,ifield,ng))
                  Finfo(9,ifield,ng)=MAX(Fmax,Finfo(9,ifield,ng))
                END DO
              ELSE
                status=nf_fread3d(ng, model, ncid, Vid, Trec, Vtype,    &
     &                            Vsize, LBi, UBi, LBj, UBj, LBk, UBk,  &
     &                            Fscale(ifield,ng), Fmin, Fmax,        &
# ifdef MASKING
     &                            Fmask(LBi,LBj),                       &
# endif
     &                            Fout(LBi,LBj,LBk,Tindex))
                Finfo(8,ifield,ng)=Fmin
                Finfo(9,ifield,ng)=Fmax
              END IF
            ELSE
              IF (InpThread) THEN
                start(1)=Trec
                total(1)=1
                status=nf90_get_var(ncid, Vid, Fval, start, total)
              END IF
# ifdef DISTRIBUTE
              CALL mp_bcasti (ng, model, status, 1)
              CALL mp_bcastf (ng, model, Fval, 1)
# endif
              Fpoint(Tindex,ifield,ng)=Fval(1)*Fscale(ifield,ng)
              Fmin=Fval(1)
              Fmax=Fval(1)
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) TRIM(Vname(1,ifield)), Trec
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
            IF (Master) THEN
              IF (Irec.gt.1) THEN
                WRITE (stdout,70) TRIM(Vname(2,ifield)), Fmin, Fmax
              ELSE
                lstr=SCAN(ncfile,'/',BACK=.TRUE.)+1
                lend=LEN_TRIM(ncfile)
                lvar=MIN(43,LEN_TRIM(Vname(2,ifield)))
                Tval(1)=Tval(1)*day2sec
                CALL time_string (Tval(1), t_code)
                WRITE (stdout,80) Vname(2,ifield)(1:lvar),              &
     &                            t_code, ncfile(lstr:lend),            &
     &                            Trec, Tindex, Tmin, Tmax, Fmin, Fmax
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
  10  FORMAT (/,' GET_3DFLDR  - unable to find requested variable: ',a, &
     &        /,15x,'in input NetCDF file: ',a)
  20  FORMAT (/,' GET_3DFLDR  - unable to open input NetCDF file: ',a)
  30  FORMAT (/,' GET_3DFLDR  - error while reading attribute: ',a)
  40  FORMAT (/,' GET_3DFLDR  - error while inquiring attributes for',  &
     &        ' variable: ',a)
  50  FORMAT (/,' GET_3DFLDR  - cannot inquire about time variable',    &
     &          ' from input NetCDF file.')
  60  FORMAT (/,' GET_3DFLDR  - error while reading variable: ',a,2x,   &
     &       ' at TIME index = ',i4)
  70  FORMAT (3x,' GET_3DFLDR  - ',a,/,19x,'(Min = ',1p,e15.8,0p,       &
     &        ' Max = ',1p,e15.8,0p,')')
  80  FORMAT (3x,' GET_3DFLDR  - ',a,',',t62,'t = ',a,/,19x,            &
     &        '(File: ',a,', Rec=',i4.4,', Index=',i1,')',/,19x,        &
     &        '(Tmin= ', f15.4, ' Tmax= ', f15.4,')',/, 19x,            &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
      RETURN
      END SUBROUTINE get_3dfldr
#else
      SUBROUTINE get_3dfldr
      RETURN
      END SUBROUTINE get_3dfldr
#endif
