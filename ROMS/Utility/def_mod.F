#include "cppdefs.h"
#if defined FOUR_DVAR || defined VERIFICATION
      SUBROUTINE def_mod (ng)
!
!svn $Id: def_mod.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine create model/observation output NetCDF which contains  !
!  model fields processed at observations locations.                   !
!                                                                      !
!  For completeness and to allow the Ensemble Kalman Filter (EnKF)     !
!  with the First Guess at Appropriate Time (FGAT), several variables  !
!  from the input observation file is also written.                    !
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
      USE strings_mod, ONLY : find_string
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
      logical :: foundit

      integer, parameter :: Natt = 25

      integer :: iterDim, datumDim, recordDim, stateDim, surveyDim
# ifdef RPCG
      integer :: datumDimp, NouterDimp
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
      integer :: MinnerDim, MouterDim, NinnerDim, NouterDim, threeDim
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      integer :: RetaDim, RxiDim
# endif
      integer :: i, j, nvatt, nvdim, status, varid, vindex
      integer :: OBSncid
      integer :: CostDim(2), Vsize(4), vardim(3)
# ifdef DISTRIBUTE
      integer :: ibuffer(2)
# endif
      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=80 )  :: type
      character (len=120)  :: Vinfo(Natt)
      character (len=256)  :: ncname
!
      SourceFile='def_mod.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=DAV(ng)%name
!
      IF (Master) THEN
        IF (LdefMOD(ng)) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
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
!  Open input observations NetCDF.
!
      CALL netcdf_open (ng, iNLM, OBS(ng)%name, 0, OBSncid)
      IF (exit_flag.ne.NoError) THEN
        WRITE (stdout,30) TRIM(OBS(ng)%name)
        RETURN
      END IF
!
!  Inquire about input observations variables.
!
      CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid)
      IF (exit_flag.ne.NoError) RETURN
!
!=======================================================================
!  Create a new model/observation file.
!=======================================================================
!
      DEFINE : IF (LdefMOD(ng)) THEN

        CALL netcdf_create (ng, iNLM, TRIM(ncname), DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,40) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define dimensions.
!-----------------------------------------------------------------------
!
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, RxiDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, RetaDim)
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'record',        &
     &                 2, recordDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'survey',        &
     &                 Nsurvey(ng), surveyDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'state_var',     &
     &                 NstateVar(ng)+1, stateDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'datum',         &
     &                 Ndatum(ng), datumDim)
        IF (exit_flag.ne.NoError) RETURN

# ifdef RPCG
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'datum+1',       &
     &                 Ndatum(ng)+1, datumDimp)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Nouter+1',      &
     &                 Nouter+1, NouterDimp)
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef FOUR_DVAR
#  if defined IS4DVAR || defined WEAK_CONSTRAINT
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Ninner',        &
     &                 Ninner, NinnerDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Minner',        &
     &                 Ninner+1, MinnerDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Nouter',        &
     &                 Nouter, NouterDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Mouter',        &
     &                 Nouter+1, MouterDim)
        IF (exit_flag.ne.NoError) RETURN
#   ifdef IS4DVAR
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'three',         &
     &                 3, threeDim)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
#  if defined IS4DVAR || defined BACKGROUND
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'iteration',     &
     &                 nf90_unlimited, iterDim)
        IF (exit_flag.ne.NoError) RETURN
!
        CostDim(1)=stateDim
        CostDim(2)=iterDim
#  endif
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
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'type', TRIM(type))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'type', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  Input observations file.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'obs_file', TRIM(OBS(ng)%name))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'obs_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  State variables IDs for observations. Copy global attribute from
!  observation file.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_copy_att(OBSncid, nf90_global,                  &
     &                           'state_variables',                     &
     &                           DAV(ng)%ncid, nf90_global)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,60) 'state_variables', 'nf90_global',       &
     &                          TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Observations provenance IDs. Copy global attribute from
!  observation file.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_copy_att(OBSncid, nf90_global,                  &
     &                           'obs_provenance',                      &
     &                           DAV(ng)%ncid, nf90_global)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,60) 'obs_provenance', 'nf90_global',        &
     &                          TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  SVN repository information.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'svn_url', TRIM(svn_url))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'svn_url', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

# ifndef DEBUGGING
#  ifdef SVN_REV
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'svn_rev', TRIM(svn_rev))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'svn_rev', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif

#  ifdef ROOT_DIR
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'code_dir', TRIM(Rdir))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'code_dir', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif

#  ifdef HEADER_DIR
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'header_dir', TRIM(Hdir))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'header_dir', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif

#  ifdef ROMS_HEADER
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'header_file', TRIM(Hfile))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'header_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
#  endif
!
!  Attributes describing platform and compiler
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'os', TRIM(my_os))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'os', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'cpu', TRIM(my_cpu))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'cpu', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'compiler_system', TRIM(my_fort))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'compiler_system', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid,nf90_global,               &
     &                          'compiler_command', TRIM(my_fc))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'compiler_command', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'compiler_flags', TRIM(my_fflags))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'compiler_flags', TRIM(ncname)
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
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'history', TRIM(history))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,50) 'history', TRIM(ncname)
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
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='inner'
        Vinfo( 2)='inner loop counter'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
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
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='obs_mean'
        Vinfo( 2)='observations mean'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='obs_std'
        Vinfo( 2)='observations standard deviation'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='model_mean'
        Vinfo( 2)='model mean'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='model_std'
        Vinfo( 2)='model standard deviation'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='model_bias'
        Vinfo( 2)='model bias'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo,ncname,             &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='SDE'
        Vinfo( 2)='model-observations standard deviation error'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='CC'
        Vinfo( 2)='model-observations cross-correlation'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='MSE'
        Vinfo( 2)='model-observations mean squared error'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='number of converged Ritz eigenvalues'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='number of converged Ritz eigenvalues'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/Nouterdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined WEAK_CONSTRAINT && \
    (defined ARRAY_MODES     || defined CLIPPING)
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='Nvct'
#  if defined ARRAY_MODES
        Vinfo( 2)='representer matrix array mode eigenvector '//        &
     &            'processed'
#  elif defined CLIPPING
        Vinfo( 2)='representer matric cut-off eigenvectors'
#  endif
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Converged Ritz eigenvalues.
!
        Vinfo( 1)='Ritz'
        Vinfo( 2)='converged Ritz eigenvalues to approximate Hessian'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
!
        Vinfo( 1)='cg_dla'
        Vinfo( 2)='normalization coefficients for Lanczos vectors'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef IS4DVAR
!
        Vinfo( 1)='cg_gamma'
        Vinfo( 2)='Lanczos algorithm gamma coefficient'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Initial gradient vector normalization factor.
!
        Vinfo( 1)='cg_Gnorm_v'
        Vinfo( 2)='initial gradient normalization factor, v-space'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='cg_Gnorm_y'
        Vinfo( 2)='initial gradient normalization factor, y-space'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='cg_Greduc_y'
        Vinfo( 2)='reduction in the gradient norm, y-space'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 3, vardim, Aval, Vinfo, ncname,                  &
     &                 SetFillVal = .FALSE.,                            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined TL_W4DPSAS          || defined W4DPSAS || \
     defined W4DPSAS_SENSITIVITY
!
!  Define NLM initial and final data penalty function.
!
        Vinfo( 1)='NL_iDataPenalty'
        Vinfo( 2)='nonlinear model initial data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='NL_fDataPenalty'
        Vinfo( 2)='nonlinear model final data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDim,NouterDim/), Aval, Vinfo, ncname,  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if (defined RECOMPUTE_4DVAR   && \
     (defined ARRAY_MODES       || defined CLIPPING)) || \
     defined TL_W4DVAR          || defined W4DVAR     || \
     defined W4DVAR_SENSITIVITY
!
!  Define RPM initial and final data penalty function.
!
        Vinfo( 1)='RP_iDataPenalty'
        Vinfo( 2)='representer model initial data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDim,NouterDim/), Aval, Vinfo, ncname,  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
        Vinfo( 1)='RP_fDataPenalty'
        Vinfo( 2)='representer model final data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDim,NouterDim/), Aval, Vinfo, ncname,  &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname,                  &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Observations survey time.
!
        IF (find_string(var_name,n_var,Vname(1,idOday),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idOday)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOday))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'units') THEN
              Vinfo(3)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'calendar') THEN
              Vinfo(4)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOday),   &
     &                   NF_FRST, 1, (/surveyDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observation type.
!
        IF (find_string(var_name,n_var,Vname(1,idOtyp),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idOTyp)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOTyp))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOTyp),   &
     &                   nf90_int, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations provenance.
!
        IF (find_string(var_name,n_var,Vname(1,idOpro),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idOpro)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOpro))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOpro),   &
     &                   nf90_int, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!                                       copy attribute from observations
          IF (OutThread) THEN
            status=nf90_copy_att(OBSncid, vindex,                       &
     &                           'flag_values',                         &
     &                           DAV(ng)%ncid, DAV(ng)%Vid(idOpro))
            IF (status.ne.nf90_noerr) THEN
               WRITE (stdout,50) 'flag_values', TRIM(Vname(1,idOpro)),  &
     &                           TRIM(ncname)
               exit_flag=3
               ioerror=status
            END IF
          END IF
# ifdef DISTRIBUTE
          ibuffer(1)=exit_flag
          ibuffer(2)=ioerror
          CALL mp_bcasti (ng, iNLM, ibuffer)
          exit_flag=ibuffer(1)
          ioerror=ibuffer(2)
# endif
          IF (exit_flag.ne.NoError) RETURN
!                                       copy attribute from observations
          IF (OutThread) THEN
            status=nf90_copy_att(OBSncid, vindex,                       &
     &                           'flag_meanings',                       &
     &                           DAV(ng)%ncid, DAV(ng)%Vid(idOpro))
            IF (status.ne.nf90_noerr) THEN
               WRITE (stdout,50) 'flag_meanings', TRIM(Vname(1,idOpro)),&
     &                           TRIM(ncname)
               exit_flag=3
               ioerror=status
            END IF
          END IF
# ifdef DISTRIBUTE
          ibuffer(1)=exit_flag
          ibuffer(2)=ioerror
          CALL mp_bcasti (ng, iNLM, ibuffer)
          exit_flag=ibuffer(1)
          ioerror=ibuffer(2)
# endif
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations time.
!
        IF (find_string(var_name,n_var,Vname(1,idObsT),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idObsT)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsT))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'units') THEN
              Vinfo(3)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'calendar') THEN
              Vinfo(4)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsT),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations longitude.
!
        IF (find_string(var_name,n_var,Vname(1,idOlon),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         MyVarName = TRIM(Vname(1,idOlon)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOlon))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'units') THEN
              Vinfo(3)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOlon),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations latitude.
!
        IF (find_string(var_name,n_var,Vname(1,idOlat),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         MyVarName = TRIM(Vname(1,idOlat)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOlat))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'units') THEN
              Vinfo(3)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOlat),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations depth.
!
        IF (find_string(var_name,n_var,Vname(1,idObsD),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         MyVarName = TRIM(Vname(1,idObsD)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsD))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'units') THEN
              Vinfo(3)=TRIM(var_Achar(i))
            ELSE IF (TRIM(var_Aname(i)).eq.'negative') THEN
              Vinfo(11)='downwards'
            END IF
          END DO
          Vinfo(17)='missing_value'
          Aval(4)=spval
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsD),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations X-fractional coordinate.
!
        IF (find_string(var_name,n_var,Vname(1,idObsX),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idObsX)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsX))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsX),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations Y-fractional coordinate.
!
        IF (find_string(var_name,n_var,Vname(1,idObsY),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idObsY)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsY))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsY),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations Z-fractional coordinate.
!
        IF (find_string(var_name,n_var,Vname(1,idObsZ),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idObsZ)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsZ))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsZ),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations total error (instrument + sampling + representation).
!
        IF (find_string(var_name,n_var,Vname(1,idOerr),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idOerr)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOerr))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOerr),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations value.
!
        IF (find_string(var_name,n_var,Vname(1,idOval),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name, OBSncid,         &
     &                         MyVarName = TRIM(Vname(1,idOval)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (exit_flag.ne.NoError) RETURN
          Vinfo(1)=TRIM(Vname(1,idOval))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOval),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Observations screening/normalization scale.
!
        Vinfo( 1)=Vname(1,idObsS)
        Vinfo( 2)=Vname(2,idObsS)
        Vinfo(24)='_FillValue'
        Aval(6)=0.0_r8
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsS),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# if defined FOUR_DVAR && !defined IS4DVAR_SENSITIVITY
!
!  Initial nonlinear model at observation locations.
!
        Vinfo( 1)=Vname(1,idNLmi)
        Vinfo( 2)=Vname(2,idNLmi)
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmi),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  ifdef BGQC
!
!  Background error at observation locations.
!
        Vinfo( 1)=Vname(1,idBgEr)
        Vinfo( 2)=Vname(2,idBgEr)
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idBgEr),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Threshold for background quality control check of observations.
!
        Vinfo( 1)=Vname(1,idBgTh)
        Vinfo( 2)=Vname(2,idBgTh)
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idBgTh),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
# if !defined IS4DVAR_SENSITIVITY && \
     (defined IS4DVAR             || defined VERIFICATION || \
      defined WEAK_CONSTRAINT)
!
!  Nonlinear model at observation points.
!
        haveNLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idNLmo)
        Vinfo( 2)=Vname(2,idNLmo)
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  ifdef VERIFICATION
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmo),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  else
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmo),     &
     &                 NF_FRST, 2, (/datumDim,NouterDim/), Aval, Vinfo, &
     &                 ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idTLmo),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  Initial model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMi)
        Vinfo( 2)=Vname(2,idMOMi)
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idMOMi),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Final model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMf)
        Vinfo( 2)=Vname(2,idMOMf)
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idMOMf),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR
!
!  Define model minus observations misfit NLM cost function.
!
        Vinfo( 1)='NLcost_function'
        Vinfo( 2)='nonlinear model misfit cost function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDim,MouterDim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define model minus observations misfit TLM cost function.
!
        Vinfo( 1)='TLcost_function'
        Vinfo( 2)='tangent linear model misfit cost function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
# ifdef RPCG
        vardim(1)=datumDimp
# else
        vardim(1)=datumDim
# endif
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# ifdef RPCG
!
!  Define initial gradient for minimization.
!
        Vinfo( 1)='vgrad0'
        Vinfo( 2)='initial gradient for minimization, v space'
        vardim(1)=datumDimp
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define sum of evolved outer-loop increments in observation space.
!
        Vinfo( 1)='Hbk'
        Vinfo( 2)='evolved sum of increments in observation space'
        vardim(1)=datumDim
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define outer-loop background cost function.
!
        Vinfo( 1)='Jb0'
        Vinfo( 2)='Outer-loop background cost function'
        vardim(1)=NouterDimp
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define Lanczos vectors in observation space.
!
# ifdef RPCG
        Vinfo( 1)='vcglwk'
        Vinfo( 2)='Preconditioned Lanczos vectors, observation space'
        vardim(1)=datumDimp
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='zcglwk'
        Vinfo( 2)='Lanczos vectors, observation space'
        vardim(1)=datumDimp
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# else
        Vinfo( 1)='zcglwk'
        Vinfo( 2)='Lanczos vectors, observation space'
        vardim(1)=datumDim
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define previous values of TLmodVal.
!
        Vinfo( 1)='TLmodVal_S'
        Vinfo( 2)='tangent linear model at observation locations'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  else
!
!  Define total observation sensitivity.
!
        Vinfo( 1)='ObsSens_total'
        Vinfo( 2)='total observation sensitivity'
        vardim(1)=datumDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Define observation impact due to surface forcing increments.
!
        Vinfo( 1)='ObsImpact_FC'
        Vinfo( 2)='observation impact due to surface forcing'
        vardim(1)=datumDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
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
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idFsur),   &
     &                   NF_FOUT, 2, (/RxiDim, RetaDim/), Aval, Vinfo,  &
     &                   ncname)
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!  Close input observations NetCDF file.
!
        CALL netcdf_close (ng, iNLM, OBSncid, OBS(ng)%name, .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing model/observation file and check its contents.
!=======================================================================
!
      QUERY : IF (.not.LdefMOD(ng)) THEN
        ncname=DAV(ng)%name
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
        CALL netcdf_open (ng, iNLM, ncname, 1, DAV(ng)%ncid)
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
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOday))) THEN
            got_var(idOday)=.TRUE.
            DAV(ng)%Vid(idOday)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOTyp))) THEN
            got_var(idOTyp)=.TRUE.
            DAV(ng)%Vid(idOTyp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOpro))) THEN
            got_var(idOpro)=.TRUE.
            DAV(ng)%Vid(idOpro)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsT))) THEN
            got_var(idObsT)=.TRUE.
            DAV(ng)%Vid(idObsT)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOlon))) THEN
            got_var(idOlon)=.TRUE.
            DAV(ng)%Vid(idOlon)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOlat))) THEN
            got_var(idOlat)=.TRUE.
            DAV(ng)%Vid(idOlat)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsD))) THEN
            got_var(idObsD)=.TRUE.
            DAV(ng)%Vid(idObsD)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsX))) THEN
            got_var(idObsX)=.TRUE.
            DAV(ng)%Vid(idObsX)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsY))) THEN
            got_var(idObsY)=.TRUE.
            DAV(ng)%Vid(idObsY)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsZ))) THEN
            got_var(idObsZ)=.TRUE.
            DAV(ng)%Vid(idObsZ)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOerr))) THEN
            got_var(idOerr)=.TRUE.
            DAV(ng)%Vid(idOerr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOval))) THEN
            got_var(idOval)=.TRUE.
            DAV(ng)%Vid(idOval)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idObsS))) THEN
            got_var(idObsS)=.TRUE.
            DAV(ng)%Vid(idObsS)=var_id(i)
# ifdef FOUR_DVAR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmi))) THEN
            got_var(idNLmi)=.TRUE.
            DAV(ng)%Vid(idNLmi)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmo))) THEN
            got_var(idNLmo)=.TRUE.
            haveNLmod(ng)=.TRUE.
            DAV(ng)%Vid(idNLmo)=var_id(i)
# if defined IS4DVAR || defined WEAK_CONSTRAINT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTLmo))) THEN
            got_var(idTLmo)=.TRUE.
            haveTLmod(ng)=.TRUE.
            DAV(ng)%Vid(idTLmo)=var_id(i)
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMOMi))) THEN
            got_var(idMOMi)=.TRUE.
            DAV(ng)%Vid(idMOMi)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMOMf))) THEN
            got_var(idMOMf)=.TRUE.
            DAV(ng)%Vid(idMOMf)=var_id(i)
# endif
          END IF
        END DO
!
!  Check if needed variables are available.
!
        IF (.not.got_var(idOday)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOday)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOpro)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOpro)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsT)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idObsT)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOlon)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOlon)),          &
     &                                  TRIM(DAV(ng)%name)
!!        exit_flag=2
!!        RETURN
        END IF
        IF (.not.got_var(idOlat)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOlat)),          &
     &                                  TRIM(DAV(ng)%name)
!!        exit_flag=2
!!        RETURN
        END IF
        IF (.not.got_var(idObsD)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idObsD)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsX)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idObsX)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsY)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idObsY)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOerr)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOerr)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOval)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOval)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF

        IF (.not.got_var(idObsS)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idObsS)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# ifdef FOUR_DVAR
        IF (.not.got_var(idNLmi)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idNLmi)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idNLmo)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idNLmo)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# if defined IS4DVAR || defined WEAK_CONSTRAINT
        IF (.not.got_var(idTLmo)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTLmo)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
# if defined IS4DVAR || defined WEAK_CONSTRAINT
        IF (.not.got_var(idMOMi)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMOMi)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idMOMf)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMOMf)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
      END IF QUERY

  10  FORMAT (/,6x,'DEF_MOD   - creating model/observation data file,', &
     &        ' Grid ',i2.2,': ',a)
  20  FORMAT (/,6x,'DEF_MOD   - inquiring model/observation data file,',&
     &        ' Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_MOD - unable to open observation/model file: ',a)
  40  FORMAT (/,' DEF_MOD - unable to create model/observation file:',  &
     &        1x,a)
  50  FORMAT (/,' DEF_MOD - unable to create global attribute: ',       &
     &        a,/,11x,a)
  60  FORMAT (/,' DEF_MOD - unable to copy attribute: ',1x,a,2x,        &
     &        'for variable: ',1x,a,/,11x,a)
  70  FORMAT (/,' DEF_MOD - unable to find model/observation variable:',&
     &        1x,a,/,11x,'in file: ',a)

      RETURN
      END SUBROUTINE def_mod
#else
      SUBROUTINE def_mod
      RETURN
      END SUBROUTINE def_mod
#endif
