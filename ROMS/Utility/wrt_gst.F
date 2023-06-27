#include "cppdefs.h"
      MODULE wrt_gst_mod
#if defined PROPAGATOR && defined CHECKPOINTING
!
!git $Id$
!svn $Id: wrt_gst.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes checkpointing fields into GST restart file       !
!  using either the standard NetCDF library or the Parallel-IO (PIO)   !
!  library.                                                            !
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
      PUBLIC  :: wrt_gst
      PRIVATE :: wrt_gst_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_gst_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_gst (ng, model)
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
!  Write out GST checkpointing fields according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (GST(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_gst_nf90 (ng, model)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_gst_pio (ng, model)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) GST(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_GST - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_gst
!
!***********************************************************************
      SUBROUTINE wrt_gst_nf90 (ng, model)
!***********************************************************************
!
      USE mod_netcdf

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
      USE distribute_mod, ONLY : mp_ncwrite1d, mp_ncwrite2d
# endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: status

# ifdef DISTRIBUTE
      integer :: Is, Ie
      integer :: vrecord = -1

      real(r8) :: scale = 1.0_r8
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_gst_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out checkpointing information variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out number of eigenvalues to compute.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'NEV',             &
     &                      NEV, (/0/), (/0/),                          &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out number of Lanczos vectors to compute.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'NCV',             &
     &                      NCV, (/0/), (/0/),                          &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out size of the eigenvalue problem.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'Mstate',          &
     &                      Mstate(ng), (/0/), (/0/),                   &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out iteration number.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iter',            &
     &                      Nrun, (/0/), (/0/),                         &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out reverse communications flag.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'ido',             &
     &                      ido(ng), (/0/), (/0/),                      &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out information and error flag.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'info',            &
     &                      info(ng), (/0/), (/0/),                     &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out eigenvalue problem type.
!
      CALL netcdf_put_svar (ng, model, GST(ng)%name, 'bmat',            &
     &                      bmat, (/1/), (/1/),                         &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out Ritz eigenvalues to compute.
!
      CALL netcdf_put_svar (ng, model, GST(ng)%name, 'which',           &
     &                      which, (/1/), (/2/),                        &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out form of basis function.
!
      CALL netcdf_put_svar (ng, model, GST(ng)%name, 'howmany',         &
     &                      howmany, (/1/), (/1/),                      &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out relative accuracy of computed Ritz values.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'Ritz_tol',        &
     &                      Ritz_tol, (/0/), (/0/),                     &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out eigenproblem parameters.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iparam',          &
     &                      iparam(:,ng), (/1/), (/SIZE(iparam)/),      &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out pointers to mark starting location in work arrays.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'ipntr',           &
     &                      ipntr(:,ng), (/1/), (/SIZE(ipntr)/),        &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal integer parameters to _aupd routines.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iaupd',           &
     &                      iaupd, (/1/), (/SIZE(iaupd)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal integer parameters to _aitr routines.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iaitr',           &
     &                      iaitr, (/1/), (/SIZE(iaitr)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal integer parameters to _aup2 routines.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iaup2',           &
     &                      iaup2, (/1/), (/SIZE(iaup2)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal logical parameters to _aitr routines.
!
      CALL netcdf_put_lvar (ng, model, GST(ng)%name, 'laitr',           &
     &                      laitr, (/1/), (/SIZE(laitr)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal logical parameters to _aupd routines.
!
      CALL netcdf_put_lvar (ng, model, GST(ng)%name, 'laup2',           &
     &                      laup2, (/1/), (/SIZE(laup2)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'raitr',           &
     &                      raitr, (/1/), (/SIZE(raitr)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'raup2',           &
     &                      raup2, (/1/), (/SIZE(raup2)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Write out Lanczos/Arnoldi basis vectors.
!
# ifdef DISTRIBUTE
      status=mp_ncwrite2d(ng, model, GST(ng)%ncid, 'Bvec',              &
     &                    GST(ng)%name, vrecord,                        &
     &                    Nstr(ng), Nend(ng), 1, NCV, scale,            &
     &                    STORAGE(ng)%Bvec(Nstr(ng):,:))
# else
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'Bvec',            &
     &                      STORAGE(ng)%Bvec(Nstr(ng):,1),              &
     &                      (/Nstr(ng),1/),                             &
     &                      (/Nend(ng)-Nstr(ng)+1,NCV/),                &
     &                      ncid = GST(ng)%ncid)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out eigenproblem residual vector.
!
# ifdef DISTRIBUTE
      status=mp_ncwrite1d(ng, model, GST(ng)%ncid, 'resid',             &
     &                    GST(ng)%name, vrecord,                        &
     &                    Nstr(ng), Nend(ng), scale,                    &
     &                    STORAGE(ng)%resid(Nstr(ng):))
# else
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'resid',           &
     &                      STORAGE(ng)%resid(Nstr(ng):),               &
     &                      (/Nstr(ng)/), (/Nend(ng)-Nstr(ng)+1/),      &
     &                      ncid = GST(ng)%ncid)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out state reverse communication work array.
!
# ifdef DISTRIBUTE
      Is=MyRank*3*Nstate(ng)+1
      Ie=MIN(Is+3*Nstate(ng)-1, 3*Mstate(ng))
      status=mp_ncwrite1d(ng, model, GST(ng)%ncid, 'SworkD',            &
     &                    GST(ng)%name, vrecord,                        &
     &                    Is, Ie, scale,                                &
     &                    STORAGE(ng)%SworkD)
# else
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'SworkD',          &
     &                      STORAGE(ng)%SworkD,                         &
     &                      (/1/), (/3*Nstate(ng)/),                    &
     &                      ncid = GST(ng)%ncid)
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out eigenproblem work array.  In distributed-memory
!  applications, this array is identical in all the nodes.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'SworkL',          &
     &                      SworkL(:,ng), (/1/), (/LworkL/),            &
     &                      ncid = GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize GST checkpointing NetCDF file to disk so the file
!  is available to other processes.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, GST(ng)%name, GST(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,10) Nrun+1
!
  10  FORMAT (2x,'WRT_GST_NF90     - wrote GST checkpointing fields',   &
     &        ' at iteration: ', i0)

      RETURN
      END SUBROUTINE wrt_gst_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_gst_pio (ng, model)
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
      integer :: status
      integer :: Is, Ie, i, j
!
      real(r8) :: scale = 1.0_r8

      real(r4), pointer :: A1d_4(:), A2d_4(:,:)  ! single precision
      real(r8), pointer :: A1d_8(:), A2d_8(:,:)  ! double precision
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_gst_pio"
!
      TYPE (Var_desc_t) :: pioVar
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out checkpointing information variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out number of eigenvalues to compute.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'NEV',         &
     &                          NEV, (/0/), (/0/),                      &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out number of Lanczos vectors to compute.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'NCV',         &
     &                          NCV, (/0/), (/0/),                      &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out size of the eigenvalue problem.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'Mstate',      &
     &                          Mstate(ng), (/0/), (/0/),               &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out iteration number.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'iter',        &
     &                          Nrun, (/0/), (/0/),                     &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out reverse communications flag.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'ido',         &
     &                          ido(ng), (/0/), (/0/),                  &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out information and error flag.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'info',        &
     &                          info(ng), (/0/), (/0/),                 &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out eigenvalue problem type.
!
      CALL pio_netcdf_put_svar (ng, model, GST(ng)%name, 'bmat',        &
     &                          bmat, (/1/), (/1/),                     &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out Ritz eigenvalues to compute.
!
      CALL pio_netcdf_put_svar (ng, model, GST(ng)%name, 'which',       &
     &                          which, (/1/), (/2/),                    &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out form of basis function.
!
      CALL pio_netcdf_put_svar (ng, model, GST(ng)%name, 'howmany',     &
     &                          howmany, (/1/), (/1/),                  &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out relative accuracy of computed Ritz values.
!
      CALL pio_netcdf_put_fvar (ng, model, GST(ng)%name, 'Ritz_tol',    &
     &                          Ritz_tol, (/0/), (/0/),                 &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out eigenproblem parameters.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'iparam',      &
     &                          iparam(:,ng), (/1/), (/SIZE(iparam)/),  &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out pointers to mark starting location in work arrays.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'ipntr',       &
     &                          ipntr(:,ng), (/1/), (/SIZE(ipntr)/),    &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal integer parameters to _aupd routines.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'iaupd',       &
     &                          iaupd, (/1/), (/SIZE(iaupd)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal integer parameters to _aitr routines.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'iaitr',       &
     &                          iaitr, (/1/), (/SIZE(iaitr)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal integer parameters to _aup2 routines.
!
      CALL pio_netcdf_put_ivar (ng, model, GST(ng)%name, 'iaup2',       &
     &                          iaup2, (/1/), (/SIZE(iaup2)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal logical parameters to _aitr routines.
!
      CALL pio_netcdf_put_lvar (ng, model, GST(ng)%name, 'laitr',       &
     &                          laitr, (/1/), (/SIZE(laitr)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write ARPACK internal logical parameters to _aupd routines.
!
      CALL pio_netcdf_put_lvar (ng, model, GST(ng)%name, 'laup2',       &
     &                          laup2, (/1/), (/SIZE(laup2)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
      CALL pio_netcdf_put_fvar (ng, model, GST(ng)%name, 'raitr',       &
     &                          raitr, (/1/), (/SIZE(raitr)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
      CALL pio_netcdf_put_fvar (ng, model, GST(ng)%name, 'raup2',       &
     &                          raup2, (/1/), (/SIZE(raup2)/),          &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Write out Lanczos/Arnoldi basis vectors.
!
      status=PIO_inq_varid(GST(ng)%pioFile, 'Bvec', pioVar)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
!
      IF (PIO_FRST.eq.PIO_double) THEN
        allocate ( A2d_8(Nstr(ng):Nend(ng),1:NCV) )
        DO j=1,NCV
          DO i=Nstr(ng),Nend(ng)
            A2d_8(i,j)=scale*STORAGE(ng)%Bvec(i,j)
          END DO
        END DO
        CALL PIO_write_darray (GST(ng)%pioFile, pioVar,                 &
     &                         ioDesc_dp_Bvec(ng),                      &
     &                         A2d_8, status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        deallocate (A2d_8)
      ELSE
        allocate ( A2d_4(Nstr(ng):Nend(ng),1:NCV) )
        DO j=1,NCV
          DO i=Nstr(ng),Nend(ng)
            A2d_4(i,j)=REAL(scale*STORAGE(ng)%Bvec(i,j), r4)
          END DO
        END DO
        CALL PIO_write_darray (GST(ng)%pioFile, pioVar,                 &
     &                         ioDesc_sp_Bvec(ng),                      &
     &                         A2d_4, status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        deallocate (A2d_4)
      END IF
!
!  Write out eigenproblem residual vector.
!
      status=PIO_inq_varid(GST(ng)%pioFile, 'resid', pioVar)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
!
      IF (PIO_FRST.eq.PIO_double) THEN
        allocate ( A1d_8(Nstr(ng):Nend(ng)) )
        DO i=Nstr(ng),Nend(ng)
          A1d_8(i)=scale*STORAGE(ng)%resid(i)
        END DO
        CALL PIO_write_darray (GST(ng)%pioFile, pioVar,                 &
     &                         ioDesc_dp_resid(ng),                     &
     &                         A1d_8, status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        deallocate (A1d_8)
      ELSE
        allocate ( A1d_4(Nstr(ng):Nend(ng)) )
        DO i=Nstr(ng),Nend(ng)
          A1d_4(i)=REAL(scale*STORAGE(ng)%resid(i), r4)
        END DO
        CALL PIO_write_darray (GST(ng)%pioFile, pioVar,                 &
     &                         ioDesc_sp_resid(ng),                     &
     &                         A1d_4, status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        deallocate (A1d_4)
      END IF
!
!  Write out state reverse communication work array.
!
      Is=MyRank*3*Nstate(ng)+1
      Ie=MIN(Is+3*Nstate(ng)-1, 3*Mstate(ng))
!
      status=PIO_inq_varid(GST(ng)%pioFile, 'SworkD', pioVar)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
!
      IF (PIO_FRST.eq.PIO_double) THEN
        allocate ( A1d_8(Is:Ie) )
        DO i=Is,Ie
          A1d_8(i)=scale*STORAGE(ng)%SworkD(i)
        END DO
        CALL PIO_write_darray (GST(ng)%pioFile, pioVar,                 &
     &                         ioDesc_dp_SworkD(ng),                    &
     &                         A1d_8, status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        deallocate (A1d_8)
      ELSE
        allocate ( A1d_4(Is:Ie) )
        DO i=Is,Ie
          A1d_4(i)=REAL(scale*STORAGE(ng)%SworkD(i), r4)
        END DO
        CALL PIO_write_darray (GST(ng)%pioFile, pioVar,                 &
     &                         ioDesc_sp_SworkD(ng),                    &
     &                         A1d_4, status)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) RETURN
        deallocate (A1d_4)
      END IF
!
!  Write out eigenproblem work array.  In distributed-memory
!  applications, this array is identical in all the nodes.
!
      CALL pio_netcdf_put_fvar (ng, model, GST(ng)%name, 'SworkL',      &
     &                          SworkL(:,ng), (/1/), (/LworkL/),        &
     &                          pioFile = GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize GST checkpointing NetCDF file to disk so the file
!  is available to other processes.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, GST(ng)%name, GST(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) WRITE (stdout,10) Nrun+1
!
  10  FORMAT (2x,'WRT_GST_PIO      - wrote GST checkpointing fields',   &
     &        ' at iteration: ', i0)

      RETURN
      END SUBROUTINE wrt_gst_pio
# endif
#endif
      END MODULE wrt_gst_mod

