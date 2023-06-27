#include "cppdefs.h"
      MODULE get_gst_mod
#if defined PROPAGATOR && defined CHECKPOINTING
!
!git $Id$
!svn $Id: get_gst.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in GST checkpointing restart file using either   !
!  the standard NetCDF library or the Parallel-IO (PIO) library.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_storage
!
      USE strings_mod,    ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: get_gst
      PRIVATE :: get_gst_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: get_gst_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_gst (ng, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Read in GST checkpointing restart file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (GST(ng)%IOtype)
        CASE (io_nf90)
          CALL get_gst_nf90 (ng, model)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_gst_pio (ng, model)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) GST(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_GST - Illegal input file type, io_type = ',i0,      &
     &        /,11x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_gst
!
!***********************************************************************
      SUBROUTINE get_gst_nf90 (ng, model)
!***********************************************************************
!
      USE mod_netcdf

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
      USE distribute_mod, ONLY : mp_ncread1d, mp_ncread2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: i, ivar, status

# ifdef DISTRIBUTE
      integer :: vrecord = -1

      real(r8) :: scale = 1.0_r8
# endif
      real(r8) :: rval
!
      character (len=1  ) :: char1
      character (len=2  ) :: char2
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_gst_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read GST checkpointing restart variables.  Check for consistency.
!-----------------------------------------------------------------------
!
!  Open checkpointing NetCDF file for reading and writing.
!
      ncname=GST(ng)%name
      IF (GST(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 1, GST(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Read in number of eigenvalues to compute.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'NEV', ivar,             &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.NEV) THEN
        IF (Master) WRITE (stdout,20) ', NEV = ', ivar, NEV
        exit_flag=6
        RETURN
      END IF
!
!  Read in number of Lanczos vectors to compute.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'NCV', ivar,             &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.NCV) THEN
        IF (Master)  WRITE (stdout,20) ', NCV = ', ivar, NCV
        exit_flag=6
        RETURN
      END IF
!
!  Read in size of the eigenvalue problem.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'Mstate', ivar,          &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.Mstate(ng)) THEN
        IF (Master) WRITE (stdout,20) ', Mstate = ', ivar, Mstate(ng)
        exit_flag=6
        RETURN
      END IF

# ifdef DISTRIBUTE
!
!  Read in number of Lanczos vectors to compute.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'Nnodes', ivar,          &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.numthreads) THEN
        IF (Master) WRITE (stdout,20) ', Nnodes = ', ivar, numthreads
        exit_flag=6
        RETURN
      END IF
# endif
!
!  Read in iteration number.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'iter', Nrun,            &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in reverse communications flag.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'ido', ido,              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in information and error flag.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'info', ido,             &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in eigenvalue problem type.
!
      CALL netcdf_get_svar (ng, model, ncname, 'bmat', char1,           &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (char1.ne.bmat) THEN
        IF (Master) WRITE (stdout,30) ', bmat = ', char1, bmat
        exit_flag=6
        RETURN
      END IF
!
!  Read in Ritz eigenvalues to compute.
!
      CALL netcdf_get_svar (ng, model, ncname, 'which', char2,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/2/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (char2(1:2).ne.which(1:2)) THEN
        IF (Master) WRITE (stdout,30) ', which = ', char2, which
        exit_flag=6
        RETURN
      END IF
!
!  Read in form of basis function.
!
      CALL netcdf_get_svar (ng, model, ncname, 'howmany', char1,        &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (char1.ne.howmany) THEN
        IF (Master) WRITE (stdout,30) ', howmany = ', char1, howmany
        exit_flag=6
        RETURN
      END IF
!
!  Read in relative accuracy of computed Ritz values.
!
      CALL netcdf_get_fvar (ng, model, ncname, 'Ritz_tol', rval,        &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (rval.ne.Ritz_tol) THEN
        IF (Master) WRITE (stdout,40) ', Ritz_tol = ', rval, Ritz_tol
      END IF
      Ritz_tol=rval
!
!  Read in eigenproblem parameters.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'iparam', iparam,        &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(iparam)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in pointers to mark starting location in work arrays.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'ipntr', ipntr,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(ipntr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal integer parameters to _aupd routines.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'iaupd', iaupd,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(iaupd)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal integer parameters to _aitr routines.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'iaitr', iaitr,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(iaitr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal integer parameters to _aup2 routines.
!
      CALL netcdf_get_ivar (ng, model, ncname, 'iaup2', iaup2,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(iaup2)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal logical parameters to _aup2 routines.
!
      CALL netcdf_get_lvar (ng, model, ncname, 'laitr', laitr,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(laitr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal logical parameters to _aup2 routines.
!
      CALL netcdf_get_lvar (ng, model, ncname, 'laup2', laup2,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(laup2)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal real parameters to _aup2 routines.
!
      CALL netcdf_get_fvar (ng, model, ncname, 'raitr', raitr,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(raitr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal real parameters to _aup2 routines.
!
      CALL netcdf_get_fvar (ng, model, ncname, 'raup2', raup2,          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/SIZE(raup2)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Read in checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Read in Lanczos/Arnoldi basis vectors.
!
# ifdef DISTRIBUTE
      status=mp_ncread2d(ng, model, GST(ng)%ncid, 'Bvec',               &
     &                   TRIM(ncname), vrecord,                         &
     &                   Nstr(ng), Nend(ng), 1, NCV, scale,             &
     &                   STORAGE(ng)%Bvec(Nstr(ng):,:))
# else
      CALL netcdf_get_fvar (ng, model, ncname, 'Bvec',                  &
     &                      STORAGE(ng)%Bvec,                           &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1,1/),                            &
     &                      total = (/Nend(ng)-Nstr(ng)+1,NCV/))
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in eigenproblem residual vector.
!
# ifdef DISTRIBUTE
      status=mp_ncread1d(ng, model, GST(ng)%ncid, 'resid',              &
     &                   TRIM(ncname), vrecord,                         &
     &                   Nstr(ng), Nend(ng), scale,                     &
     &                   STORAGE(ng)%resid(Nstr(ng):))
# else
      CALL netcdf_get_fvar (ng, model, ncname, 'resid',                 &
     &                      STORAGE(ng)%resid,                          &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/Nend(ng)-Nstr(ng)+1/))
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in state reverse communication work array.
!
# ifdef DISTRIBUTE
      status=mp_ncread1d(ng, model, GST(ng)%ncid, 'SworkD',             &
     &                   TRIM(ncname), vrecord,                         &
     &                   1, 3*Nstate(ng), scale,                        &
     &                   STORAGE(ng)%SworkD)
# else
      CALL netcdf_get_fvar (ng, model, ncname, 'SworkD',                &
     &                      STORAGE(ng)%SworkD,                         &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/3*Nstate(ng)/))
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in eigenproblem work array.
!
      CALL netcdf_get_fvar (ng, model, ncname, 'SworkL',                &
     &                      SworkL(:,ng),                               &
     &                      ncid = GST(ng)%ncid,                        &
     &                      start = (/1/),                              &
     &                      total = (/LworkL/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' GET_GST_NF90 - unable to open checkpointing NetCDF',  &
     &          ' file:', a)
  20  FORMAT (/,' GET_GST_NF90 - inconsistent input parameter', a, 2i4)
  30  FORMAT (/,' GET_GST_NF90 - inconsistent input parameter', a, a, a)
  40  FORMAT (/,' GET_GST_NF90 - input parameter', a, 1pe10.2,0p,       &
     &        /, 16x,'has been reset to: ', 1pe10.2)
!
      RETURN
      END SUBROUTINE get_gst_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_gst_pio (ng, model)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: Is, Ie
      integer :: i, ivar, status

# ifdef DISTRIBUTE
      integer :: vrecord = -1

      real(r8) :: scale = 1.0_r8
# endif
      real(r8) :: rval
!
      character (len=1  ) :: char1
      character (len=2  ) :: char2
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_gst_pio"
!
      TYPE (Var_desc_t) :: pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read GST checkpointing restart variables.  Check for consistency.
!-----------------------------------------------------------------------
!
!  Open checkpointing NetCDF file for reading and writing.
!
      ncname=GST(ng)%name
      IF (GST(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 1, GST(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Read in number of eigenvalues to compute.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'NEV', ivar,         &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.NEV) THEN
        IF (Master) WRITE (stdout,20) ', NEV = ', ivar, NEV
        exit_flag=6
        RETURN
      END IF
!
!  Read in number of Lanczos vectors to compute.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'NCV', ivar,         &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.NCV) THEN
        IF (Master)  WRITE (stdout,20) ', NCV = ', ivar, NCV
        exit_flag=6
        RETURN
      END IF
!
!  Read in size of the eigenvalue problem.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'Mstate', ivar,      &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.Mstate(ng)) THEN
        IF (Master) WRITE (stdout,20) ', Mstate = ', ivar, Mstate(ng)
        exit_flag=6
        RETURN
      END IF
!
!  Read in number of Lanczos vectors to compute.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'Nnodes', ivar,      &
     &                          piofile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (ivar.ne.numthreads) THEN
        IF (Master) WRITE (stdout,20) ', Nnodes = ', ivar, numthreads
        exit_flag=6
        RETURN
      END IF
!
!  Read in iteration number.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'iter', Nrun,        &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in reverse communications flag.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'ido', ido,          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in information and error flag.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'info', ido,         &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in eigenvalue problem type.
!
      CALL pio_netcdf_get_svar (ng, model, ncname, 'bmat', char1,       &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (char1.ne.bmat) THEN
        IF (Master) WRITE (stdout,30) ', bmat = ', char1, bmat
        exit_flag=6
        RETURN
      END IF
!
!  Read in Ritz eigenvalues to compute.
!
      CALL pio_netcdf_get_svar (ng, model, ncname, 'which', char2,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/2/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (char2(1:2).ne.which(1:2)) THEN
        IF (Master) WRITE (stdout,30) ', which = ', char2, which
        exit_flag=6
        RETURN
      END IF
!
!  Read in form of basis function.
!
      CALL pio_netcdf_get_svar (ng, model, ncname, 'howmany', char1,    &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (char1.ne.howmany) THEN
        IF (Master) WRITE (stdout,30) ', howmany = ', char1, howmany
        exit_flag=6
        RETURN
      END IF
!
!  Read in relative accuracy of computed Ritz values.
!
      CALL pio_netcdf_get_fvar (ng, model, ncname, 'Ritz_tol', rval,    &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (rval.ne.Ritz_tol) THEN
        IF (Master) WRITE (stdout,40) ', Ritz_tol = ', rval, Ritz_tol
      END IF
      Ritz_tol=rval
!
!  Read in eigenproblem parameters.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'iparam', iparam,    &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(iparam)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in pointers to mark starting location in work arrays.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'ipntr', ipntr,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(ipntr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal integer parameters to _aupd routines.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'iaupd', iaupd,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(iaupd)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal integer parameters to _aitr routines.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'iaitr', iaitr,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(iaitr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal integer parameters to _aup2 routines.
!
      CALL pio_netcdf_get_ivar (ng, model, ncname, 'iaup2', iaup2,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(iaup2)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal logical parameters to _aup2 routines.
!
      CALL pio_netcdf_get_lvar (ng, model, ncname, 'laitr', laitr,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(laitr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal logical parameters to _aup2 routines.
!
      CALL pio_netcdf_get_lvar (ng, model, ncname, 'laup2', laup2,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(laup2)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal real parameters to _aup2 routines.
!
      CALL pio_netcdf_get_fvar (ng, model, ncname, 'raitr', raitr,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(raitr)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in ARPACK internal real parameters to _aup2 routines.
!
      CALL pio_netcdf_get_fvar (ng, model, ncname, 'raup2', raup2,      &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/SIZE(raup2)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Read in checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Read in Lanczos/Arnoldi basis vectors.
!
      status=PIO_inq_varid(GST(ng)%pioFile, 'Bvec', pioVar)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
!
      IF (KIND(STORAGE(ng)%Bvec).eq.8) THEN
        CALL PIO_read_darray (GST(ng)%pioFile, pioVar,                  &
     &                        ioDesc_dp_Bvec(ng),                       &
     &                        STORAGE(ng)%Bvec(Nstr(ng):,:), status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
      ELSE
        CALL PIO_read_darray (GST(ng)%pioFile, pioVar,                  &
     &                        ioDesc_sp_Bvec(ng),                       &
     &                        STORAGE(ng)%Bvec(Nstr(ng):,:), status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
      END IF
!
!  Read in eigenproblem residual vector.
!
      status=PIO_inq_varid(GST(ng)%pioFile, 'resid', pioVar)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
!
      IF (KIND(STORAGE(ng)%Bvec).eq.8) THEN
        CALL PIO_read_darray (GST(ng)%pioFile, pioVar,                  &
     &                        ioDesc_dp_resid(ng),                      &
     &                        STORAGE(ng)%resid(Nstr(ng):), status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
      ELSE
        CALL PIO_read_darray (GST(ng)%pioFile, pioVar,                  &
     &                        ioDesc_dp_resid(ng),                      &
     &                        STORAGE(ng)%resid(Nstr(ng):), status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
      END IF
!
!  Read in state reverse communication work array.
!
      Is=MyRank*3*Nstate(ng)+1
      Ie=MIN(Is+3*Nstate(ng)-1, 3*Mstate(ng))
!
      status=PIO_inq_varid(GST(ng)%pioFile, 'SworkD', pioVar)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
!
      IF (KIND(STORAGE(ng)%SworkD).eq.8) THEN
        CALL PIO_read_darray (GST(ng)%pioFile, pioVar,                  &
     &                        ioDesc_dp_SworkD(ng),                     &
     &                        STORAGE(ng)%SworkD(Nstr(ng):), status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
      ELSE
        CALL PIO_read_darray (GST(ng)%pioFile, pioVar,                  &
     &                        ioDesc_sp_SworkD(ng),                     &
     &                        STORAGE(ng)%SworkD(Nstr(ng):), status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
      END IF
!
!  Read in eigenproblem work array.
!
      CALL pio_netcdf_get_fvar (ng, model, ncname, 'SworkL',            &
     &                          SworkL(:,ng),                           &
     &                          pioFile = GST(ng)%pioFile,              &
     &                          start = (/1/),                          &
     &                          total = (/LworkL/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' GET_GST_PIO - unable to open checkpointing NetCDF',   &
     &          ' file:', a)
  20  FORMAT (/,' GET_GST_PIO - inconsistent input parameter', a, 2i4)
  30  FORMAT (/,' GET_GST_PIO - inconsistent input parameter', a, a, a)
  40  FORMAT (/,' GET_GST_PIO - input parameter', a, 1pe10.2,0p,        &
     &        /, 16x,'has been reset to: ', 1pe10.2)
!
      RETURN
      END SUBROUTINE get_gst_pio
# endif
#endif
      END MODULE get_gst_mod


