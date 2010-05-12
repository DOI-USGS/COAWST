#include "cppdefs.h"
#if defined FOUR_DVAR || defined VERIFICATION
      SUBROUTINE def_mod (ng)
!
!svn $Id: def_mod.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine create model/observation output NetCDF which contains  !
!  model fields  processed at observations points.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings
!
      USE def_var_mod, ONLY : def_var
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(NV) :: got_var(NV)

      integer, parameter :: Natt = 25

      integer :: iterDim, datumDim, recordDim, stateDim
# if defined IS4DVAR || defined WEAK_CONSTRAINT
      integer :: MinnerDim, MouterDim, NinnerDim, NouterDim, threeDim
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      integer :: RetaDim, RxiDim
# endif
      integer :: i, j, status, varid
      integer :: CostDim(2), Vsize(4), vardim(3)
# ifdef DISTRIBUTE
      integer :: ibuffer(2)
# endif
      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname, type
!
      SourceFile='def_mod.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=MODname(ng)
!
      IF (Master) THEN
        IF (LdefMOD(ng)) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!=======================================================================
!  Create a new model/observation file.
!=======================================================================
!
      DEFINE : IF (LdefMOD(ng)) THEN

        CALL netcdf_create (ng, iNLM, TRIM(ncname), ncMODid(ng))
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define dimensions.
!-----------------------------------------------------------------------
!
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'xi_rho',         &
     &                 IOBOUNDS(ng)%xi_rho, RxiDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'eta_rho',        &
     &                 IOBOUNDS(ng)%eta_rho, RetaDim)
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'record',         &
     &                 2, recordDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'state_var',      &
     &                 NstateVar(ng)+1, stateDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'datum',          &
     &                 Ndatum(ng), datumDim)
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
#  if defined IS4DVAR || defined WEAK_CONSTRAINT
        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'Ninner',         &
     &                 Ninner, NinnerDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'Minner',         &
     &                 Ninner+1, MinnerDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'Nouter',         &
     &                 Nouter, NouterDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'Mouter',         &
     &                 Nouter+1, MouterDim)
        IF (exit_flag.ne.NoError) RETURN
#   ifdef IS4DVAR
        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'three',          &
     &                 3, threeDim)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

        status=def_dim(ng, iNLM, ncMODid(ng), ncname, 'iteration',      &
     &                 nf90_unlimited, iterDim)
        IF (exit_flag.ne.NoError) RETURN
!
        CostDim(1)=stateDim
        CostDim(2)=iterDim
# endif
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
        IF (OutThread) THEN
!
!  File type.
!
          IF (exit_flag.eq.NoError) THEN
            type='ROMS/TOMS 4DVAR output observation processing file'
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'type', TRIM(type))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'type', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  Input observations file.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'obs_file', TRIM(OBSname(ng)))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'obs_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  SVN repository information.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'svn_url', TRIM(svn_url))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'svn_url', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

# ifndef DEBUGGING
#  ifdef SVN_REV
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'svn_rev', TRIM(svn_rev))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'svn_rev', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif

#  ifdef ROOT_DIR
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'code_dir', TRIM(Rdir))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'code_dir', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif

#  ifdef HEADER_DIR
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'header_dir', TRIM(Hdir))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'header_dir', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif

#  ifdef ROMS_HEADER
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'header_file', TRIM(Hfile))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'header_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif
!
!  Attributes describing platform and compiler
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'os', TRIM(my_os))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'os', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'cpu', TRIM(my_cpu))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'cpu', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'compiler_system', TRIM(my_fort))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'compiler_system', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng),nf90_global,                &
     &                          'compiler_command', TRIM(my_fc))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'compiler_command', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'compiler_flags', TRIM(my_fflags))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'compiler_flags', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  History attribute.
!
          IF (exit_flag.eq.NoError) THEN
            IF (LEN_TRIM(date_str).gt.0) THEN
              WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',   &
     &                                          TRIM(version),          &
     &                                          TRIM(date_str)
            ELSE
              WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version',          &
     &                                   TRIM(version)
            END IF
            status=nf90_put_att(ncMODid(ng), nf90_global,               &
     &                          'history', TRIM(history))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) 'history', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
# endif
        END IF
# ifdef DISTRIBUTE
        ibuffer(1)=exit_flag
        ibuffer(2)=ioerror
        CALL mp_bcasti (ng, iNLM, ibuffer)
        exit_flag=ibuffer(1)
        ioerror=ibuffer(2)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------

# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Outer and inner loop contours.
!
        Vinfo( 1)='outer'
        Vinfo( 2)='outer loop counter'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_int,          &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='inner'
        Vinfo( 2)='inner loop counter'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_int,          &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifndef IS4DVAR_SENSITIVITY
!
!  Define model-observation comparison statistics.
!
        Vinfo( 1)='Nobs'
        Vinfo( 2)='Number of usable observations'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_int,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='obs_mean'
        Vinfo( 2)='observations mean'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='obs_std'
        Vinfo( 2)='observations standard deviation'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='model_mean'
        Vinfo( 2)='model mean'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='model_std'
        Vinfo( 2)='model standard deviation'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='model_bias'
        Vinfo( 2)='model bias'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo,ncname,             &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='SDE'
        Vinfo( 2)='model-observations standard deviation error'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='CC'
        Vinfo( 2)='model-observations cross-correlation'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='MSE'
        Vinfo( 2)='model-observations mean squared error'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='Number of converged Ritz eigenvalues'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_int,          &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='Number of converged Ritz eigenvalues'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_int,          &
     &                 1, (/Nouterdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Converged Ritz eigenvalues.
!
        Vinfo( 1)='Ritz'
        Vinfo( 2)='converged Ritz eigenvalues to approximate Hessian'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/Ninnerdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Converged Ritz eigenvalues.
!
        Vinfo( 1)='Ritz'
        Vinfo( 2)='converged Ritz eigenvalues to approximate Hessian'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Define conjugate gradient norm.
!
        Vinfo( 1)='cg_beta'
        Vinfo( 2)='conjugate gradient beta coefficient'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Define Lanczos algorithm coefficients.
!
        Vinfo( 1)='cg_delta'
        Vinfo( 2)='Lanczos algorithm delta coefficient'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#  ifdef WEAK_CONSTRAINT
        Vinfo( 1)='cg_dla'
        Vinfo( 2)='normalization coefficients for Lanczos vectors'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef IS4DVAR
        Vinfo( 1)='cg_gamma'
        Vinfo( 2)='Lanczos algorithm gamma coefficient'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# if defined IS4DVAR
!
!  Initial gradient vector normalization factor.
!
        Vinfo( 1)='cg_Gnorm'
        Vinfo( 2)='initial gradient normalization factor'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Initial gradient vector normalization factor.
!
        Vinfo( 1)='cg_Gnorm_v'
        Vinfo( 2)='initial gradient normalization factor, v-space'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='cg_Gnorm_y'
        Vinfo( 2)='initial gradient normalization factor, y-space'
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Lanczos vector normalization factor.
!
        Vinfo( 1)='cg_QG'
        Vinfo( 2)='Lanczos vector normalization factor'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Reduction in the gradient norm.
!
        Vinfo( 1)='cg_Greduc'
        Vinfo( 2)='reduction in the gradient norm'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Reduction in the gradient norm.
!
        Vinfo( 1)='cg_Greduc_v'
        Vinfo( 2)='reduction in the gradient norm, v-space'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='cg_Greduc_y'
        Vinfo( 2)='reduction in the gradient norm, y-space'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Lanczos recurrence tridiagonal matrix.
!
        Vinfo( 1)='cg_Tmatrix'
        Vinfo( 2)='Lanczos recurrence tridiagonal matrix'
        vardim(1)=NinnerDim
        vardim(2)=threeDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Lanczos tridiagonal matrix, upper diagonal elements.
!
        Vinfo( 1)='cg_zu'
        Vinfo( 2)='tridiagonal matrix, upper diagonal elements'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Eigenvalues of Lanczos recurrence relationship.
!
        Vinfo( 1)='cg_Ritz'
        Vinfo( 2)='Lanczos recurrence eigenvalues'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Eigenvalues relative error.
!
        Vinfo( 1)='cg_RitzErr'
        Vinfo( 2)='Ritz eigenvalues relative error'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Eigenvectors of Lanczos recurrence relationship.
!
        Vinfo( 1)='cg_zv'
        Vinfo( 2)='Lanczos recurrence eigenvectors'
        vardim(1)=NinnerDim
        vardim(2)=NinnerDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Eigenvectors of Lanczos recurrence relationship.
!
        Vinfo( 1)='cg_zv'
        Vinfo( 2)='Lanczos recurrence eigenvectors'
        vardim(1)=NinnerDim
        vardim(2)=NinnerDim
        vardim(3)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, NF_FRST,           &
     &                 3, vardim, Aval, Vinfo, ncname,                  &
     &                 SetFillVal = .FALSE.,                            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined TL_W4DVAR || defined W4DVAR || defined W4DVAR_SENSITIVITY
!
!  Define RPM data penalty function.
!
        Vinfo( 1)='RPcost_function'
        Vinfo( 2)='representer model data penalty function'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef WEAK_CONSTRAINT
!
!  Define first guess initial data misfit.
!
        Vinfo( 1)='Jf'
        Vinfo( 2)='first guess initial data misfit'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define state estimate data misfit.
!
        Vinfo( 1)='Jdata'
        Vinfo( 2)='state estimate data misfit'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model penalty function.
!
        Vinfo( 1)='Jmod'
        Vinfo( 2)='model penalty function'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define optimal penalty function.
!
        Vinfo( 1)='Jopt'
        Vinfo( 2)='optimal penalty function'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define actual model penalty function.
!
        Vinfo( 1)='Jb'
        Vinfo( 2)='actual model penalty function'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define actual data penalty function.
!
        Vinfo( 1)='Jobs'
        Vinfo( 2)='actual data penalty function'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define actual data penalty function.
!
        Vinfo( 1)='Jact'
        Vinfo( 2)='actual total penalty function'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Observations screening/normalization scale.
!
        Vinfo( 1)=Vname(1,idObsS)
        Vinfo( 2)=Vname(2,idObsS)
        status=def_var(ng, iNLM, ncMODid(ng), modVid(idObsS,ng),        &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# if defined FOUR_DVAR && !defined IS4DVAR_SENSITIVITY
!
!  Initial nonlinear model at observation locations.
!
        Vinfo( 1)=Vname(1,idNLmi)
        Vinfo( 2)=Vname(2,idNLmi)
        status=def_var(ng, iNLM, ncMODid(ng), modVid(idNLmi,ng),        &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifndef IS4DVAR_SENSITIVITY
!
!  Nonlinear model at observation points.
!
        haveNLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idNLmo)
        Vinfo( 2)=Vname(2,idNLmo)
        status=def_var(ng, iNLM, ncMODid(ng), modVid(idNLmo,ng),        &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR         || defined IS4DVAR_SENSITIVITY || \
     defined WEAK_CONSTRAINT
!
!  Tangent linear or representer model at observation points.
!
        haveTLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idTLmo)
#  ifdef IS4DVAR_SENSITIVITY
        Vinfo( 2)='4DVAR sensitivity analysis at observations location'
#  else
        Vinfo( 2)=Vname(2,idTLmo)
#  endif
        status=def_var(ng, iNLM, ncMODid(ng), modVid(idTLmo,ng),        &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Initial model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMi)
        Vinfo( 2)=Vname(2,idMOMi)
        status=def_var(ng, iNLM, ncMODid(ng), modVid(idMOMi,ng),        &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Final model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMf)
        Vinfo( 2)=Vname(2,idMOMf)
        status=def_var(ng, iNLM, ncMODid(ng), modVid(idMOMf,ng),        &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Define model minus observations misfit NLM cost function.
!
        Vinfo( 1)='NLcost_function'
        Vinfo( 2)='nonlinear model misfit cost function'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, (/stateDim,MouterDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model minus observations misfit TLM cost function.
!
        Vinfo( 1)='TLcost_function'
        Vinfo( 2)='tangent linear model misfit cost function'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, (/iterDim/), Aval, Vinfo, ncname,             &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef BACKGROUND
!
!  Define model minus background misfit cost function.
!
        Vinfo( 1)='back_function'
        Vinfo( 2)='model minus background misfit cost function'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, (/iterDim/), Aval, Vinfo, ncname,             &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Define optimality property that measures the consistency between
!  background and observation errors hypotheses (Chi-square).
!
        Vinfo( 1)='Jmin'
        Vinfo( 2)='normalized, optimal cost function minimum'
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, (/iterDim/), Aval, Vinfo, ncname,             &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef WEAK_CONSTRAINT
!
!  Define initial gradient for minimization.
!
        Vinfo( 1)='zgrad0'
        Vinfo( 2)='initial gradient for minimization, observation space'
        vardim(1)=datumDim
        vardim(2)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Lanczos vectors in observation space.
!
        Vinfo( 1)='zcglwk'
        Vinfo( 2)='Lanczos vectors, observation space'
        vardim(1)=datumDim
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define previous values of TLmodVal.
!
        Vinfo( 1)='TLmodVal_S'
        Vinfo( 2)='tangent linear model at observation locations'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        vardim(3)=NouterDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined TL_W4DVAR          || defined W4DVAR             || \
     defined W4DVAR_SENSITIVITY
!
!  Define initial values of RPmodVal.
!
        Vinfo( 1)='RPmodel_initial'
        Vinfo( 2)='initial representer model at observation locations'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined W4DVAR_SENSITIVITY  || defined W4DPSAS_SENSITIVITY || \
     defined IS4DVAR_SENSITIVITY
#  ifdef OBS_IMPACT
!
!  Define total observations impact.
!
        Vinfo( 1)='ObsImpact_total'
        Vinfo( 2)='total observation impact'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  else
!
!  Define total observation sensitivity.
!
        Vinfo( 1)='ObsSens_total'
        Vinfo( 2)='total observation sensitivity'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# if defined OBS_IMPACT_SPLIT     && \
     (defined W4DVAR_SENSITIVITY  || defined W4DPSAS_SENSITIVITY || \
      defined IS4DVAR_SENSITIVITY)
!
!  Define observation impact due to initial condition increments.
!
        Vinfo( 1)='ObsImpact_IC'
        Vinfo( 2)='observation impact due to initial conditions'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Define observation impact due to surface forcing increments.
!
        Vinfo( 1)='ObsImpact_FC'
        Vinfo( 2)='observation impact due to surface forcing'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  if defined ADJUST_BOUNDARY
!
!  Define observation impact due to boundary condition increments.
!
        Vinfo( 1)='ObsImpact_BC'
        Vinfo( 2)='observation impact due to open boundary conditions'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, ncMODid(ng), varid, nf90_double,       &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
!
!  Define reference free-surface used in the balance operator.
!
        IF (balance(isFsur)) THEN
          Vinfo( 1)='zeta_ref'
          Vinfo( 2)='reference free-surface, balance operator'
          Vinfo( 3)=Vname(3,idFsur)
          status=def_var(ng, iNLM, ncMODid(ng), modVid(idFsur,ng),      &
     &                   NF_FOUT, 2, (/RxiDim, RetaDim/), Aval, Vinfo,  &
     &                   ncname)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing model/observation file and check its contents.
!=======================================================================
!
      QUERY : IF (.not.LdefMOD(ng)) THEN
        ncname=MODname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open model/observation for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, ncMODid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from model/observation NetCDF and activate
!  switches for required variables.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsS))) THEN
            got_var(idObsS)=.TRUE.
            modVid(idObsS,ng)=var_id(i)
# ifdef FOUR_DVAR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmi))) THEN
            got_var(idNLmi)=.TRUE.
            modVid(idNLmi,ng)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmo))) THEN
            got_var(idNLmo)=.TRUE.
            haveNLmod(ng)=.TRUE.
            modVid(idNLmo,ng)=var_id(i)
# if defined IS4DVAR || defined WEAK_CONSTRAINT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTLmo))) THEN
            got_var(idTLmo)=.TRUE.
            haveTLmod(ng)=.TRUE.
            modVid(idTLmo,ng)=var_id(i)
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMOMi))) THEN
            got_var(idMOMi)=.TRUE.
            modVid(idMOMi,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMOMf))) THEN
            got_var(idMOMf)=.TRUE.
            modVid(idMOMf,ng)=var_id(i)
# endif
          END IF
        END DO
!
!  Check if needed variables are available.
!
        IF (.not.got_var(idObsS)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idObsS)),          &
     &                                  TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# ifdef FOUR_DVAR
        IF (.not.got_var(idNLmi)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idNLmi)),          &
     &                                  TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idNLmo)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idNLmo)),          &
     &                                  TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# if defined IS4DVAR || defined WEAK_CONSTRAINT
        IF (.not.got_var(idTLmo)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTLmo)),          &
     &                                  TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
        IF (.not.got_var(idMOMi)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMOMi)),          &
     &                                  TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idMOMf)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idMOMf)),          &
     &                                  TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
      END IF QUERY

  10  FORMAT (/,6x,'DEF_MOD   - creating model/observation data file:', &
     &        1x,a)
  20  FORMAT (/,6x,'DEF_MOD   - inquiring model/observation data file:',&
     &        1x,a)
  30  FORMAT (/,' DEF_MOD - unable to create model/observation file:',  &
     &        1x,a)
  40  FORMAT (/,' DEF_MOD - unable to create globat attribute: ',       &
     &        a,/,11x,a)
  50  FORMAT (/,' DEF_MOD - unable to open observation/model file: ',a)

  60  FORMAT (/,' DEF_MOD - unable to find model/observation variable:',&
     &        1x,a,/,11x,'in file: ',a)

      RETURN
      END SUBROUTINE def_mod
#else
      SUBROUTINE def_mod
      RETURN
      END SUBROUTINE def_mod
#endif
