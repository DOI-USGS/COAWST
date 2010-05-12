#include "cppdefs.h"
#if defined PROPAGATOR && defined CHECKPOINTING
      SUBROUTINE wrt_gst (ng, model)
!
!svn $Id: wrt_gst.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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
      USE distribute_mod, ONLY : mp_bcasti, mp_ncwrite
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
      integer :: i, status, varid
      integer :: start(4), total(4)

# ifdef DISTRIBUTE
      integer :: vrecord = -1

      real(r8) :: scale = 1.0_r8

      character (len=6) :: var
# endif
      character (len=1) :: lchar(5)
!
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
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'NEV',              &
     &                      NEV, (/0/), (/0/),                          &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of Lanczos vectors to compute.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'NCV',              &
     &                      NCV, (/0/), (/0/),                          &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out size of the eigenvalue problem.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'Mstate',           &
     &                      Mstate(ng), (/0/), (/0/),                   &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef DISTRIBUTE
!
!  Write out number of distributed-memory nodes.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'Nnodes',           &
     &                      numnodes, (/0/), (/0/),                     &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out iteration number.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'iter',             &
     &                      Nrun, (/0/), (/0/),                         &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out reverse communications flag.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'ido',              &
     &                      ido, (/0/), (/0/),                          &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out information and error flag.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'info',             &
     &                      info, (/0/), (/0/),                         &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenvalue problem type.
!
      CALL netcdf_put_svar (ng, model, GSTname(ng), 'bmat',             &
     &                      bmat, (/1/), (/1/),                         &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out Ritz eigenvalues to compute.
!
      CALL netcdf_put_svar (ng, model, GSTname(ng), 'which',            &
     &                      which, (/1/), (/2/),                        &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out form of basis function.
!
      CALL netcdf_put_svar (ng, model, GSTname(ng), 'howmany',          &
     &                      howmany, (/1/), (/1/),                      &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out relative accuracy of computed Ritz values.
!
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'Ritz_tol',         &
     &                      Ritz_tol, (/0/), (/0/),                     &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenproblem parameters.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'iparam',           &
     &                      iparam, (/1/), (/SIZE(iparam)/),            &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pointers to mark starting location in work arrays.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'ipntr',            &
     &                      ipntr, (/1/), (/SIZE(ipntr)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal integer parameters to _aupd routines.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'iaupd',            &
     &                      iaupd, (/1/), (/SIZE(iaupd)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal integer parameters to _aitr routines.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'iaitr',            &
     &                      iaitr, (/1/), (/SIZE(iaitr)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal integer parameters to _aup2 routines.
!
      CALL netcdf_put_ivar (ng, model, GSTname(ng), 'iaup2',            &
     &                      iaup2, (/1/), (/SIZE(iaup2)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal logical parameters to _aitr routines.
!
      DO i=1,SIZE(laitr)
        IF (laitr(i)) THEN
          lchar(i)='T'
        ELSE
          lchar(i)='F'
        END IF
      END DO
      CALL netcdf_put_svar (ng, model, GSTname(ng), 'laitr',            &
     &                      lchar, (/1/), (/SIZE(laitr)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write ARPACK internal logical parameters to _aupd routines.
!
      DO i=1,SIZE(laup2)
        IF (laup2(i)) THEN
          lchar(i)='T'
        ELSE
          lchar(i)='F'
        END IF
      END DO
      CALL netcdf_put_svar (ng, model, GSTname(ng), 'laup2',            &
     &                      lchar, (/1/), (/SIZE(laup2)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal real parameters to _aitr routines.
!
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'raitr',            &
     &                      raitr, (/1/), (/SIZE(raitr)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'raup2',            &
     &                      raup2, (/1/), (/SIZE(raup2)/),              &
     &                      ncid = ncGSTid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Write out Lanczos/Arnoldi basis vectors.
!
# ifdef DISTRIBUTE
      var='Bvec'
      status=mp_ncwrite (ng, model, ncGSTid(ng), var, GSTname(ng),      &
     &                   vrecord, Nstr(ng), Nend(ng), 1, NCV, scale,    &
     &                   Bvec(Nstr(ng):,1))
# else
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'Bvec',             &
     &                      Bvec(Nstr(ng):,1),                          &
     &                      (/Nstr(ng),1/),                             &
     &                      (/Nend(ng)-Nstr(ng)+1,NCV/),                &
     &                      ncid = ncGSTid(ng))
# endif
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenproblem residual vector.
!
# ifdef DISTRIBUTE
      var='resid'
      status=mp_ncwrite(ng, model, ncGSTid(ng), var, GSTname(ng),       &
     &                  vrecord, Nstr(ng), Nend(ng), 1, 1, scale,       &
     &                  resid(Nstr(ng):))
# else
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'resid',            &
     &                      resid(Nstr(ng):),                           &
     &                      (/Nstr(ng)/), (/Nend(ng)-Nstr(ng)+1/),      &
     &                      ncid = ncGSTid(ng))
# endif
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out state reverse communication work array.
# ifdef DISTRIBUTE
!  Notice zero arguments indicating node dimension in NetCDF file.
# endif
!
# ifdef DISTRIBUTE
      var='SworkD'
      status=mp_ncwrite(ng, model, ncGSTid(ng), var, GSTname(ng),       &
     &                  vrecord, 1, 3*Nstate(ng), 0, 0, scale,          &
     &                  SworkD)
# else
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'SworkD',           &
     &                      SworkD, (/1/), (/3*Nstate(ng)/),            &
     &                      ncid = ncGSTid(ng))
# endif
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out eigenproblem work array.
# ifdef DISTRIBUTE
!  Notice zero arguments indicating node dimension in NetCDF file.
# endif
!
# ifdef DISTRIBUTE
      var='SworkL'
      status=mp_ncwrite(ng, model, ncGSTid(ng), var, GSTname(ng),       &
     &                  vrecord, 1, LworkL, 0, 0, scale,                &
     &                  SworkL)
# else
      CALL netcdf_put_fvar (ng, model, GSTname(ng), 'SworkL',           &
     &                      SworkL, (/1/), (/LworkL/),                  &
     &                      ncid = ncGSTid(ng))
# endif
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize GST checkpointing NetCDF file to disk so the file
!  is available to other processes.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, GSTname(ng), ncGSTid(ng))
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
