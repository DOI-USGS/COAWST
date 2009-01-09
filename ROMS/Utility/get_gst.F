#include "cppdefs.h"
#if defined PROPAGATOR && defined CHECKPOINTING
      SUBROUTINE get_gst (ng, model)
!
!svn $Id: get_gst.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in GST checkpointing restart NetCDF file.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_storage

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
      USE distribute_mod, ONLY : mp_ncread
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: i, ivar, status, tile, varid

      integer, dimension(2) :: start, total

# ifdef DISTRIBUTE
      integer :: vrecord = -1

      real(r8) :: scale =1.0_r8
# endif
      real(r8) :: rval

      character (len=1 ) :: char1, lchar(5)
      character (len=2 ) :: char2
      character (len=6 ) :: var
      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Read GST checkpointing restart variables.  Chack for consistency.
!-----------------------------------------------------------------------
!
!  Open checkpointing NetCDF file for reading and writing.
!
      ncname=GSTname(ng)
      IF (InpThread.and.(ncGSTid(ng).eq.-1)) THEN
        status=nf90_open(TRIM(ncname), nf90_write, ncGSTid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Read in number of eigenvalues to compute.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'NEV', varid)
        status=nf90_get_var(ncGSTid(ng), varid, ivar)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'NEV', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (ivar.ne.NEV) THEN
          WRITE (stdout,30) ', NEV = ', ivar, NEV
          exit_flag=6
          RETURN
        END IF
      END IF
!
!  Read in number of Lanczos vectors to compute.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'NCV', varid)
        status=nf90_get_var(ncGSTid(ng), varid, ivar)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'NCV', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (ivar.ne.NCV) THEN
          WRITE (stdout,30) ', NCV = ', ivar, NCV
          exit_flag=6
          RETURN
        END IF
      END IF
!
!  Read in size of the eigenvalue problem.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'Mstate', varid)
        status=nf90_get_var(ncGSTid(ng), varid, ivar)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'Mstate', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (ivar.ne.Mstate(ng)) THEN
          WRITE (stdout,30) ', Mstate = ', ivar, Mstate(ng)
          exit_flag=6
          RETURN
        END IF
      END IF

# ifdef DISTRIBUTE
!
!  Read in number of Lanczos vectors to compute.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'Nnodes', varid)
        status=nf90_get_var(ncGSTid(ng), varid, ivar)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'Nnodes', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (ivar.ne.numnodes) THEN
          WRITE (stdout,30) ', Nnodes = ', ivar, numnodes
          exit_flag=6
          RETURN
        END IF
      END IF
# endif
!
!  Read in iteration number.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'iter', varid)
        status=nf90_get_var(ncGSTid(ng), varid, Nrun)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'iter', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, Nrun, 1)
# endif
!
!  Read in reverse communications flag.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'ido', varid)
        status=nf90_get_var(ncGSTid(ng), varid, ido)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'ido', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, ido, 1)
# endif
!
!  Read in information and error flag.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'info', varid)
        status=nf90_get_var(ncGSTid(ng), varid, info)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'info', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, info, 1)
# endif
!
!  Read in eigenvalue problem type.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'bmat', varid)
        status=nf90_get_var(ncGSTid(ng), varid, char1)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'bmat', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (char1.ne.bmat) THEN
          WRITE (stdout,40) ', bmat = ', char1, bmat
          exit_flag=6
          RETURN
        END IF
      END IF
!
!  Read in Ritz eigenvalues to compute.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'which', varid)
        start(1)=1
        total(1)=2
        status=nf90_get_var(ncGSTid(ng), varid, char2, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'which', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (char2(1:2).ne.which(1:2)) THEN
          WRITE (stdout,40) ', which = ', char2, which
          exit_flag=6
          RETURN
        END IF
      END IF
!
!  Read in form of basis function.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'howmany', varid)
        status=nf90_get_var(ncGSTid(ng), varid, char1)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'howmany', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (char1.ne.howmany) THEN
          WRITE (stdout,40) ', howmany = ', char1, howmany
          exit_flag=6
          RETURN
        END IF
      END IF
!
!  Read in relative accuracy of computed Ritz values.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'Ritz_tol', varid)
        status=nf90_get_var(ncGSTid(ng), varid, rval)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'Ritz_tol', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        IF (rval.ne.Ritz_tol) THEN
          WRITE (stdout,50) ', Ritz_tol = ', rval, Ritz_tol
        END IF
        Ritz_tol=rval
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, Ritz_tol, 1)
# endif
!
!  Read in eigenproblem parameters.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'iparam', varid)
        start(1)=1
        total(1)=SIZE(iparam)
        status=nf90_get_var(ncGSTid(ng), varid, iparam, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'iparam', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, iparam, SIZE(iparam))
# endif
!
!  Read in pointers to mark starting location in work arrays.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'ipntr', varid)
        start(1)=1
        total(1)=SIZE(ipntr)
        status=nf90_get_var(ncGSTid(ng), varid, ipntr, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'ipntr', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, ipntr, SIZE(ipntr))
# endif
!
!  Read in ARPACK internal integer parameters to _aupd routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'iaupd', varid)
        start(1)=1
        total(1)=SIZE(iaupd)
        status=nf90_get_var(ncGSTid(ng), varid, iaupd, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'iaupd', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, iaupd, SIZE(iaupd))
# endif
!
!  Read in ARPACK internal integer parameters to _aitr routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'iaitr', varid)
        start(1)=1
        total(1)=SIZE(iaitr)
        status=nf90_get_var(ncGSTid(ng), varid, iaitr, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'iaitr', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, iaitr, SIZE(iaitr))
# endif
!
!  Read in ARPACK internal integer parameters to _aup2 routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'iaup2', varid)
        start(1)=1
        total(1)=SIZE(iaup2)
        status=nf90_get_var(ncGSTid(ng), varid, iaup2, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'iaup2', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, iaup2, SIZE(iaup2))
# endif
!
!  Read in ARPACK internal logical parameters to _aup2 routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'laitr', varid)
        start(1)=1
        total(1)=SIZE(laitr)
        status=nf90_get_var(ncGSTid(ng), varid, lchar, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'laitr', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        DO i=1,SIZE(laitr)
          IF (lchar(i).eq.'T') THEN
            laitr(i)=.TRUE.
          ELSE
            laitr(i)=.FALSE.
          END IF
        END DO
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastl (ng, model, laitr, SIZE(laitr))
# endif
!
!  Read in ARPACK internal logical parameters to _aup2 routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'laup2', varid)
        start(1)=1
        total(1)=SIZE(laup2)
        status=nf90_get_var(ncGSTid(ng), varid, lchar, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'laup2', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        DO i=1,SIZE(laup2)
          IF (lchar(i).eq.'T') THEN
            laup2(i)=.TRUE.
          ELSE
            laup2(i)=.FALSE.
          END IF
        END DO

      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastl (ng, model, laup2, SIZE(laup2))
# endif
!
!  Read in ARPACK internal real parameters to _aup2 routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'raitr', varid)
        start(1)=1
        total(1)=SIZE(raitr)
        status=nf90_get_var(ncGSTid(ng), varid, raitr, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'raitr', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, raitr, SIZE(raitr))
# endif
!
!  Read in ARPACK internal real parameters to _aup2 routines.
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncGSTid(ng), 'raup2', varid)
        start(1)=1
        total(1)=SIZE(raup2)
        status=nf90_get_var(ncGSTid(ng), varid, raup2, start, total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'raup2', TRIM(ncname)
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, raup2, SIZE(raup2))
# endif
!
!-----------------------------------------------------------------------
!  Read in checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Read in Lanczos/Arnoldi basis vectors.
!
# ifdef DISTRIBUTE
      var='Bvec'
      status=mp_ncread(ng, model, ncGSTid(ng), var, ncname,             &
     &                 vrecord, Nstr(ng), Nend(ng), 1, NCV, scale,      &
     &                 Bvec(Nstr(ng):,1))
# else
      start(1)=Nstr(ng)
      total(1)=Nend(ng)-Nstr(ng)+1
      start(2)=1
      total(2)=NCV
      status=nf90_inq_varid(ncGSTid(ng), 'Bvec', varid)
      status=nf90_get_var(ncGSTid(ng), varid, Bvec(Nstr(ng):,1),        &
     &                    start, total)
# endif
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,60) 'Bvec', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in eigenproblem residual vector.
!
# ifdef DISTRIBUTE
      var='resid'
      status=mp_ncread(ng, model, ncGSTid(ng), var, ncname,             &
     &                 vrecord, Nstr(ng), Nend(ng), 1, 1, scale,        &
     &                 resid(Nstr(ng):))
# else
      start(1)=Nstr(ng)
      total(1)=Nend(ng)-Nstr(ng)+1
      status=nf90_inq_varid(ncGSTid(ng), 'resid', varid)
      status=nf90_get_var(ncGSTid(ng), varid, resid(Nstr(ng):),         &
     &                    start, total)
# endif
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,60) 'resid', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in state reverse communication work array.
# ifdef DISTRIBUTE
!  Notice zero arguments indicating node dimension in NetCDF file.
# endif
!
# ifdef DISTRIBUTE
      var='SworkD'
      status=mp_ncread(ng, model, ncGSTid(ng), var, ncname,             &
     &                 vrecord, 1, 3*Nstate(ng), 0, 0, scale,           &
     &                 SworkD)
# else
      start(1)=1
      total(1)=3*Nstate(ng)
      status=nf90_inq_varid(ncGSTid(ng), 'SworkD', varid)
      status=nf90_get_var(ncGSTid(ng), varid, SworkD,                   &
     &                    start, total)
# endif
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,60) 'SworkD', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
!  Read in eigenproblem work array.
# ifdef DISTRIBUTE
!  Notice zero arguments indicating node dimension in NetCDF file.
# endif
!
# ifdef DISTRIBUTE
      var='SworkL'
      status=mp_ncread(ng, model, ncGSTid(ng), var, ncname,             &
     &                 vrecord, 1, LworkL, 0, 0, scale,                 &
     &                 SworkL)
# else
      start(1)=1
      total(1)=LworkL
      status=nf90_inq_varid(ncGSTid(ng), 'SworkL', varid)
      status=nf90_put_var(ncGSTid(ng), varid, SworkL,                   &
     &                    start, total)
# endif
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,20) 'SworkL', TRIM(ncname)
        exit_flag=2
        ioerror=status
        RETURN
      END IF
!
  10  FORMAT (/,' GET_GST - unable to open checkpointing NetCDF',       &
     &          ' file:', a)
  20  FORMAT (/,' GET_GST - error while reading variable: ', a,         &
     &        /, 11x, 'in GST checkpointing NetCDF file: ', a)
  30  FORMAT (/,' GET_GST - inconsistent input parameter', a, 2i4)
  40  FORMAT (/,' GET_GST - inconsistent input parameter', a, a, a)
  50  FORMAT (/,' GET_GST - input parameter', a, 1pe10.4,0p,            &
     &        /, 11x,'has been reset to: ', 1pe10.4)
  60  FORMAT (/,' GET_GST - error while reading variable: ', a,         &
     &        /, 11x, 'in GST checkpointing NetCDF file: ', a)

      RETURN
      END SUBROUTINE get_gst
#else
      SUBROUTINE get_gst
      RETURN
      END SUBROUTINE get_gst
#endif

