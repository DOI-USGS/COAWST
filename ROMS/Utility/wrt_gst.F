#include "cppdefs.h"
#if defined PROPAGATOR && defined CHECKPOINTING
      SUBROUTINE wrt_gst (ng, model)
!
!svn $Id: wrt_gst.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes checkpointing fields into GST restart NetCDF    !
!  file.                                                               !
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

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
      USE distribute_mod, ONLY : mp_ncwrite1d
      USE distribute_mod, ONLY : mp_ncwrite2d
#endif
!
      implicit none
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
      SourceFile='wrt_gst.F'
!
!-----------------------------------------------------------------------
!  Write out checkpointing information variables.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of eigenvalues to compute.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'NEV',             &
     &                      NEV, (/0/), (/0/),                          &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of Lanczos vectors to compute.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'NCV',             &
     &                      NCV, (/0/), (/0/),                          &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out size of the eigenvalue problem.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'Mstate',          &
     &                      Mstate(ng), (/0/), (/0/),                   &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out iteration number.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iter',            &
     &                      Nrun, (/0/), (/0/),                         &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out reverse communications flag.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'ido',             &
     &                      ido(ng), (/0/), (/0/),                      &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out information and error flag.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'info',            &
     &                      info(ng), (/0/), (/0/),                     &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenvalue problem type.
!
      CALL netcdf_put_svar (ng, model, GST(ng)%name, 'bmat',            &
     &                      bmat, (/1/), (/1/),                         &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out Ritz eigenvalues to compute.
!
      CALL netcdf_put_svar (ng, model, GST(ng)%name, 'which',           &
     &                      which, (/1/), (/2/),                        &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out form of basis function.
!
      CALL netcdf_put_svar (ng, model, GST(ng)%name, 'howmany',         &
     &                      howmany, (/1/), (/1/),                      &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out relative accuracy of computed Ritz values.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'Ritz_tol',        &
     &                      Ritz_tol, (/0/), (/0/),                     &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenproblem parameters.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iparam',          &
     &                      iparam(:,ng), (/1/), (/SIZE(iparam)/),      &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pointers to mark starting location in work arrays.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'ipntr',           &
     &                      ipntr(:,ng), (/1/), (/SIZE(ipntr)/),        &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal integer parameters to _aupd routines.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iaupd',           &
     &                      iaupd, (/1/), (/SIZE(iaupd)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal integer parameters to _aitr routines.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iaitr',           &
     &                      iaitr, (/1/), (/SIZE(iaitr)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal integer parameters to _aup2 routines.
!
      CALL netcdf_put_ivar (ng, model, GST(ng)%name, 'iaup2',           &
     &                      iaup2, (/1/), (/SIZE(iaup2)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal logical parameters to _aitr routines.
!
      CALL netcdf_put_lvar (ng, model, GST(ng)%name, 'laitr',           &
     &                      laitr, (/1/), (/SIZE(laitr)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal logical parameters to _aupd routines.
!
      CALL netcdf_put_lvar (ng, model, GST(ng)%name, 'laup2',           &
     &                      laup2, (/1/), (/SIZE(laup2)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'raitr',           &
     &                      raitr, (/1/), (/SIZE(raitr)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'raup2',           &
     &                      raup2, (/1/), (/SIZE(raup2)/),              &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
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
     &                    STORAGE(ng)%Bvec(Nstr(ng):,1))
# else
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'Bvec',            &
     &                      STORAGE(ng)%Bvec(Nstr(ng):,1),              &
     &                      (/Nstr(ng),1/),                             &
     &                      (/Nend(ng)-Nstr(ng)+1,NCV/),                &
     &                      ncid = GST(ng)%ncid)
# endif
      IF (exit_flag.ne.NoError) RETURN
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
      IF (exit_flag.ne.NoError) RETURN
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
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenproblem work array.  In distributed-memory
!  applications, this array is identical in all the nodes.
!
      CALL netcdf_put_fvar (ng, model, GST(ng)%name, 'SworkL',          &
     &                      SworkL, (/1/), (/LworkL/),                  &
     &                      ncid = GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize GST checkpointing NetCDF file to disk so the file
!  is available to other processes.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, GST(ng)%name, GST(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,10) Nrun+1
!
  10  FORMAT (6x,'WRT_GST   - wrote GST checkpointing fields at ',      &
     &        'iteration: ', i5.5)

      RETURN
      END SUBROUTINE wrt_gst
#else
      SUBROUTINE wrt_gst
      RETURN
      END SUBROUTINE wrt_gst
#endif
