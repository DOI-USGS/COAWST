#include "cppdefs.h"
#if defined FOUR_DVAR || defined VERIFICATION
      SUBROUTINE def_mod (ng)
!
!svn $Id: def_mod.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
      USE dateclock_mod,  ONLY : time_string
      USE def_var_mod,    ONLY : def_var
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE strings_mod,    ONLY : find_string
      USE strings_mod,    ONLY : FoundError
      USE strings_mod,    ONLY : uppercase
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
      logical :: foundAtt(2), foundit
!
      integer, parameter :: Natt = 25

      integer :: iterDim, datumDim, recordDim
      integer :: stateDim, stateDimp, surveyDim
      integer :: Fcount
# ifdef RPCG
      integer :: datumDimp, NouterDimp
# endif
# if defined I4DVAR || defined WEAK_CONSTRAINT
      integer :: MinnerDim, MouterDim, NinnerDim, NouterDim, threeDim
# endif
# ifdef SP4DVAR
      integer :: NinnerDim1
      integer :: NinnerDim2
# endif
# if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
      integer :: RetaDim, RxiDim
# endif
      integer :: i, j, lstr, nvatt, nvdim, status, varid, vindex
      integer :: OBSncid
      integer :: vardim(3)
# ifdef DISTRIBUTE
      integer :: ibuffer(2)
# endif
      integer :: def_dim
!
      real(dp) :: Tstart, Tfinal
      real(r8) :: Aval(6)
!
      character (len=22 )  :: str_date, end_date
      character (len=40 )  :: Aname, AttName(2)
      character (len=80 )  :: string
      character (len=120)  :: Vinfo(Natt)
      character (len=256)  :: ncname
      character (len=2048) :: AttValue(2)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
!=======================================================================
!  Create a new model/observation file.
!=======================================================================
!
      DEFINE : IF (LdefMOD(ng)) THEN
!
!  Open input observations NetCDF.
!
        IF (OBS(ng)%ncid.eq.-1) THEN
          CALL netcdf_open (ng, iNLM, OBS(ng)%name, 1, OBSncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
            WRITE (stdout,30) TRIM(OBS(ng)%name)
            RETURN
          END IF
        ELSE
          OBSncid=OBS(ng)%ncid
        END IF
!
!  Inquire about input observations variables.
!
        CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                    &
     &                       ncid = OBSncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire if the 'state_variables' and 'obs_provenance' attributes
!  are available in the observations file.
!
        AttName(1)='state_variables'
        AttName(2)='obs_provenance'
!
        CALL netcdf_get_satt (ng, iNLM, OBS(ng)%name, nf90_global,      &
     &                        AttName, AttValue, foundAtt,              &
     &                        ncid = OBSncid)
!
!  Create model/observation (DAV) file
!
        CALL netcdf_create (ng, iNLM, TRIM(ncname), DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, RetaDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'record',        &
     &                 2, recordDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'survey',        &
     &                 Nsurvey(ng), surveyDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'state_var',     &
     &                 NobsVar(ng), stateDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'cost_var',      &
     &                 NobsVar(ng)+1, stateDimp)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'datum',         &
     &                 Ndatum(ng), datumDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef RPCG
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'datum+1',       &
     &                 Ndatum(ng)+1, datumDimp)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Nouter+1',      &
     &                 Nouter+1, NouterDimp)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef FOUR_DVAR
#  if defined I4DVAR || defined WEAK_CONSTRAINT
#   ifdef SP4DVAR
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Ninner1',       &
     &                 Ninner+1, NinnerDim1)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Ninner2',       &
     &                 Ninner+2, NinnerDim2)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   endif
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Ninner',        &
     &                 Ninner, NinnerDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Minner',        &
     &                 Ninner+1, MinnerDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Nouter',        &
     &                 Nouter, NouterDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'Mouter',        &
     &                 Nouter+1, MouterDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   ifdef I4DVAR
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'three',         &
     &                 3, threeDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
#  if defined I4DVAR || defined BACKGROUND
        status=def_dim(ng, iNLM, DAV(ng)%ncid, ncname, 'iteration',     &
     &                 nf90_unlimited, iterDim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
            string='ROMS/TOMS 4D-Var output observation processing file'
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'type', TRIM(string))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'type', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  Algorithm.
!
          IF (exit_flag.eq.NoError) THEN
# if defined ARRAY_MODES
            string=uppercase('array_modes')
# elif defined SPLIT_I4DVAR
            string=uppercase('split_i4dvar')
# elif defined I4DVAR
            string=uppercase('i4dvar')
# elif defined I4DVAR_ANA_SENSITIVITY
            string=uppercase('i4dvar_ana_sensitivity')
# elif defined SPLIT_RBL4DVAR
            string=uppercase('split_rbl4dvar')
# elif defined RBL4DVAR
            string=uppercase('rbl4dvar')
# elif defined RBL4DVAR_ANA_SENSITIVITY
            string=uppercase('rbl4dvar_ana_sensitivity')
# elif defined RBL4DVAR_FCT_SENSITIVITY
            string=uppercase('rbl4dvar_fct_sensitivity')
# elif defined SPLIT_R4DVAR
            string=uppercase('split_r4dvar')
# elif defined R4DVAR
            string=uppercase('r4dvar')
# elif defined R4DVAR_ANA_SENSITIVITY
            string=uppercase('r4dvar_ana_sensitivity')
# elif defined TL_RBL4DVAR
            string=uppercase('tl_rbl4dvar')
# elif defined TL_R4DVAR
            string=uppercase('tl_r4dvar')
# elif defined VERIFICATION
            string=uppercase('verification')
# else
            string=uppercase('four_dvar')
# endif
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'Algorithm', TRIM(string))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'type', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
!  Set cycle start and final time and dates attributes.
!
          IF (INItime(ng).lt.0.0_dp) THEN
            Tstart=time(ng)+dt(ng)             ! called after 'initial'
          ELSE
            Tstart=INItime(ng)
          END IF
          Tfinal=Tstart+ntimes(ng)*dt(ng)
          CALL time_string (Tstart, str_date)
          CALL time_string (Tfinal, end_date)
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'str_day', Tstart*sec2day)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'str_day', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'end_day', Tfinal*sec2day)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'end_day', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'str_date', str_date)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'str_date', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'end_date', end_date)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'end_date', TRIM(ncname)
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (foundAtt(1)) THEN
              status=nf90_copy_att(OBSncid, nf90_global,                &
     &                             'state_variables',                   &
     &                             DAV(ng)%ncid, nf90_global)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                WRITE (stdout,50) 'state_variables', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
!
!  Observations provenance IDs. Copy global attribute from
!  observation file.
!
          IF (exit_flag.eq.NoError) THEN
            IF (foundAtt(2)) THEN
              status=nf90_copy_att(OBSncid, nf90_global,                &
     &                             'obs_provenance',                    &
     &                             DAV(ng)%ncid, nf90_global)
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                WRITE (stdout,50) 'obs_provenance', TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
              IF (FoundError(exit_flag, NoError,                        &
     &                       __LINE__, MyFile)) RETURN
            END IF
          END IF
!
!  SVN repository information.
!
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'svn_url', TRIM(svn_url))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'os', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'cpu', TRIM(my_cpu))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'cpu', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'compiler_system', TRIM(my_fort))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'compiler_system', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(DAV(ng)%ncid,nf90_global,               &
     &                          'compiler_command', TRIM(my_fc))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              WRITE (stdout,50) 'compiler_command', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            lstr=INDEX(my_fflags, 'free')-2
            IF (lstr.le.0) lstr=LEN_TRIM(my_fflags)
            status=nf90_put_att(DAV(ng)%ncid, nf90_global,              &
     &                          'compiler_flags', my_fflags(1:lstr))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------

# if defined I4DVAR || defined WEAK_CONSTRAINT
!
!  Outer and inner loop counters.
!
        Vinfo( 1)='outer'
        Vinfo( 2)='outer loop counter'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='inner'
        Vinfo( 2)='inner loop counter'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined ARRAY_MODES              || \
     defined I4DVAR_ANA_SENSITIVITY   || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined R4DVAR_ANA_SENSITIVITY
!
!  Observations impact/sensitivity outer loop beeing processed.
!
      Vinfo( 1)='Nimpact'
      Vinfo( 2)='Observations impact/sensitivity outer loop to use'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifndef I4DVAR_ANA_SENSITIVITY
!
!  Define model-observation comparison statistics.
!
        Vinfo( 1)='Nobs'
        Vinfo( 2)='number of observations with the same survey time'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/surveyDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='Nused_obs'
        Vinfo( 2)='Number of usable observations'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='obs_mean'
        Vinfo( 2)='observations mean'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='obs_std'
        Vinfo( 2)='observations standard deviation'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='model_mean'
        Vinfo( 2)='model mean'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='model_std'
        Vinfo( 2)='model standard deviation'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='model_bias'
        Vinfo( 2)='model bias'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo,ncname,             &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='SDE'
        Vinfo( 2)='model-observations standard deviation error'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='CC'
        Vinfo( 2)='model-observations cross-correlation'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='MSE'
        Vinfo( 2)='model-observations mean squared error'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/stateDim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='number of converged Ritz eigenvalues'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Number of converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='number of converged Ritz eigenvalues'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_int,         &
     &                 1, (/Nouterdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR
!
!  Converged Ritz eigenvalues.
!
        Vinfo( 1)='Ritz'
        Vinfo( 2)='converged Ritz eigenvalues to approximate Hessian'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/Ninnerdim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR || defined WEAK_CONSTRAINT
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR || defined WEAK_CONSTRAINT
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef I4DVAR
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
# if defined I4DVAR
!
!  Initial gradient vector normalization factor.
!
        Vinfo( 1)='cg_Gnorm'
        Vinfo( 2)='initial gradient normalization factor'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# elif defined WEAK_CONSTRAINT
!
!  Initial gradient vector normalization factor.
!
        Vinfo( 1)='cg_Gnorm_v'
        Vinfo( 2)='initial gradient normalization factor, v-space'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='cg_Gnorm_y'
        Vinfo( 2)='initial gradient normalization factor, y-space'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, NF_FRST,          &
     &                 1, (/NouterDim/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR || defined WEAK_CONSTRAINT
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR || defined WEAK_CONSTRAINT
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined I4DVAR
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if defined RBL4DVAR                 || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined SP4DVAR                  || \
     defined TL_RBL4DVAR
!
!  Define NLM initial and final data penalty function.
!
        Vinfo( 1)='NL_iDataPenalty'
        Vinfo( 2)='nonlinear model initial data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, (/stateDimp/), Aval, Vinfo, ncname,           &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='NL_fDataPenalty'
        Vinfo( 2)='nonlinear model final data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDimp,NouterDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# if (defined RECOMPUTE_4DVAR       && \
     (defined ARRAY_MODES           || defined CLIPPING)) || \
     defined R4DVAR_ANA_SENSITIVITY || defined R4DVAR     || \
     defined TL_R4DVAR
!
!  Define RPM initial and final data penalty function.
!
        Vinfo( 1)='RP_iDataPenalty'
        Vinfo( 2)='representer model initial data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDimp,NouterDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='RP_fDataPenalty'
        Vinfo( 2)='representer model final data penalty function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDimp,NouterDim/), Aval, Vinfo, ncname, &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  ifndef SP4DVAR
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
# endif
!
!  Observations survey time.
!
        IF (find_string(var_name,n_var,Vname(1,idOday),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOday)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
     &                   NF_TOUT, 1, (/surveyDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observation type.
!
        IF (find_string(var_name,n_var,Vname(1,idOtyp),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOTyp)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idOTyp))
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOTyp),   &
     &                   nf90_int, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!                                      copy attributes from observations
          IF (OutThread) THEN
            DO i=1,nvatt
              status=nf90_inq_attname(OBSncid, vindex, i, Aname)
              IF (status.eq.nf90_noerr) THEN
                status=nf90_copy_att(OBSncid, vindex, TRIM(Aname),      &
     &                               DAV(ng)%ncid, DAV(ng)%Vid(idOTyp))
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  WRITE (stdout,60) TRIM(Aname), TRIM(Vname(1,idOTyp)), &
     &                              TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                WRITE (stdout,70) i, TRIM(Vname(1,idOpro)),             &
     &                           TRIM(OBS(ng)%name)
              END IF
            END DO
          END IF
# ifdef DISTRIBUTE
          ibuffer(1)=exit_flag
          ibuffer(2)=ioerror
          CALL mp_bcasti (ng, iNLM, ibuffer)
          exit_flag=ibuffer(1)
          ioerror=ibuffer(2)
# endif
        END IF
!
!  Observations provenance.
!
        IF (find_string(var_name,n_var,Vname(1,idOpro),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOpro)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idOpro))
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOpro),   &
     &                   nf90_int, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!                                      copy attributes from observations
          IF (OutThread) THEN
            DO i=1,nvatt
              status=nf90_inq_attname(OBSncid, vindex, i, Aname)
              IF (status.eq.nf90_noerr) THEN
                status=nf90_copy_att(OBSncid, vindex, TRIM(Aname),      &
     &                               DAV(ng)%ncid, DAV(ng)%Vid(idOpro))
                IF (FoundError(status, nf90_noerr,                      &
     &                         __LINE__, MyFile)) THEN
                  WRITE (stdout,60) TRIM(Aname), TRIM(Vname(1,idOpro)), &
     &                              TRIM(ncname)
                  exit_flag=3
                  ioerror=status
                END IF
              ELSE
                WRITE (stdout,70) i, TRIM(Vname(1,idOpro)),             &
     &                           TRIM(OBS(ng)%name)
              END IF
            END DO
          END IF
# ifdef DISTRIBUTE
          ibuffer(1)=exit_flag
          ibuffer(2)=ioerror
          CALL mp_bcasti (ng, iNLM, ibuffer)
          exit_flag=ibuffer(1)
          ioerror=ibuffer(2)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations time.
!
        IF (find_string(var_name,n_var,Vname(1,idObsT),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idObsT)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
     &                   NF_TOUT, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations longitude.
!
        IF (find_string(var_name,n_var,Vname(1,idOlon),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOlon)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations latitude.
!
        IF (find_string(var_name,n_var,Vname(1,idOlat),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOlat)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations depth.
!
        IF (find_string(var_name,n_var,Vname(1,idObsD),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idObsD)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations X-fractional coordinate.
!
        IF (find_string(var_name,n_var,Vname(1,idObsX),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idObsX)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsX))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsX),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations Y-fractional coordinate.
!
        IF (find_string(var_name,n_var,Vname(1,idObsY),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idObsY)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsY))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsY),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations Z-fractional coordinate.
!
        IF (find_string(var_name,n_var,Vname(1,idObsZ),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idObsZ)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idObsZ))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsZ),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations total error (instrument + sampling + representation).
!
        IF (find_string(var_name,n_var,Vname(1,idOerr),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOerr)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idOerr))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          Vinfo(3)=TRIM(Vname(3,idOerr))
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOerr),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations value.
!
        IF (find_string(var_name,n_var,Vname(1,idOval),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOval)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idOval))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          Vinfo(3)=TRIM(Vname(3,idOval))
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOval),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations meta value.
!
        IF (find_string(var_name,n_var,Vname(1,idOmet),vindex)) THEN
          CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                  &
     &                         ncid = OBSncid,                          &
     &                         MyVarName = TRIM(Vname(1,idOmet)),       &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          Vinfo(1)=TRIM(Vname(1,idOmet))
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'long_name') THEN
              Vinfo(2)=TRIM(var_Achar(i))
            END IF
          END DO
          status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idOmet),   &
     &                   NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Observations screening/normalization scale.
!
        Vinfo( 1)=Vname(1,idObsS)
        Vinfo( 2)=Vname(2,idObsS)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idObsS),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined FOUR_DVAR && !defined I4DVAR_ANA_SENSITIVITY
!
!  Unvetted prior nonlinear model at observation locations.
!
        Vinfo( 1)=Vname(1,idNLmp)
        Vinfo( 2)=Vname(2,idNLmp)
        Vinfo( 3)=Vname(3,idNLmp)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmp),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initial nonlinear model at observation locations.
!
        Vinfo( 1)=Vname(1,idNLmi)
        Vinfo( 2)=Vname(2,idNLmi)
        Vinfo( 3)=Vname(3,idNLmi)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmi),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Final nonlinear model at observation locations.
!
        Vinfo( 1)=Vname(1,idNLmf)
        Vinfo( 2)=Vname(2,idNLmf)
        Vinfo( 3)=Vname(3,idNLmf)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmf),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  4D-Var background error standard deviation at observation locations.
!
        Vinfo( 1)=Vname(1,idBgEr)
        Vinfo( 2)=Vname(2,idBgEr)
        Vinfo( 3)=Vname(3,idBgEr)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idBgEr),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  4D-Var innovation vector: observation minus background.
!
        Vinfo( 1)=Vname(1,idInno)
        Vinfo( 2)=Vname(2,idInno)
        Vinfo( 3)=Vname(3,idInno)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idInno),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  4D-Var increment vector: analysis minus background.
!
        Vinfo( 1)=Vname(1,idIncr)
        Vinfo( 2)=Vname(2,idIncr)
        Vinfo( 3)=Vname(3,idIncr)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idIncr),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  4D-Var residual vector: observation minus analysis.
!
        Vinfo( 1)=Vname(1,idResi)
        Vinfo( 2)=Vname(2,idResi)
        Vinfo( 3)=Vname(3,idResi)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idResi),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef BGQC
!
!  Threshold for background quality control check of observations.
!
        Vinfo( 1)=Vname(1,idBgTh)
        Vinfo( 2)=Vname(2,idBgTh)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idBgTh),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if !defined I4DVAR_ANA_SENSITIVITY && \
     (defined I4DVAR                 || defined WEAK_CONSTRAINT)
!
!  Unvetted nonlinear model at observation points per outer-loop.
!
        haveNLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idNLmu)
        Vinfo( 2)=Vname(2,idNLmu)
        Vinfo( 3)=Vname(3,idNLmu)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmu),     &
     &                 NF_FRST, 2, (/datumDim,NouterDim/), Aval, Vinfo, &
     &                 ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !defined I4DVAR_ANA_SENSITIVITY && \
     (defined I4DVAR                 || defined VERIFICATION || \
      defined WEAK_CONSTRAINT)
!
!  Nonlinear model at observation points.
!
        haveNLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idNLmo)
        Vinfo( 2)=Vname(2,idNLmo)
        Vinfo( 3)=Vname(3,idNLmo)
#  ifdef VERIFICATION
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmo),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  else
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idNLmo),     &
     &                 NF_FRST, 2, (/datumDim,NouterDim/), Aval, Vinfo, &
     &                 ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined I4DVAR          || defined I4DVAR_ANA_SENSITIVITY || \
     defined WEAK_CONSTRAINT
!
!  Tangent linear or representer model at observation points.
!
        haveTLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idTLmo)
#  ifdef I4DVAR_ANA_SENSITIVITY
        Vinfo( 2)='4D-Var sensitivity analysis at observations location'
#  else
        Vinfo( 2)=Vname(2,idTLmo)
#  endif
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idTLmo),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined I4DVAR || defined WEAK_CONSTRAINT
!
!  Initial model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMi)
        Vinfo( 2)=Vname(2,idMOMi)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idMOMi),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Final model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMf)
        Vinfo( 2)=Vname(2,idMOMf)
        status=def_var(ng, iNLM, DAV(ng)%ncid, DAV(ng)%Vid(idMOMf),     &
     &                 NF_FRST, 1, (/datumDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined I4DVAR
!
!  Define model minus observations misfit NLM cost function.
!
        Vinfo( 1)='NLcost_function'
        Vinfo( 2)='nonlinear model misfit cost function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, (/stateDimp,MouterDim/), Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define model minus observations misfit TLM cost function.
!
        Vinfo( 1)='TLcost_function'
        Vinfo( 2)='tangent linear model misfit cost function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, (/iterDim/), Aval, Vinfo, ncname,             &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef BACKGROUND
!
!  Define model minus background misfit cost function.
!
        Vinfo( 1)='back_function'
        Vinfo( 2)='model minus background misfit cost function'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, (/iterDim/), Aval, Vinfo, ncname,             &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined I4DVAR
!
!  Define optimality property that measures the consistency between
!  background and observation errors hypotheses (Chi-square).
!
        Vinfo( 1)='Jmin'
        Vinfo( 2)='normalized, optimal cost function minimum'
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, (/iterDim/), Aval, Vinfo, ncname,             &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef WEAK_CONSTRAINT
#  ifndef RPCG
!
!  Representer coefficients estimate. Needed for restarting "congrad".
!
        Vinfo( 1)='cg_pxsave'
        Vinfo( 2)='representer coefficients estimatate'
        vardim(1)=datumDim
        vardim(2)=NouterDim
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Define initial gradient for minimization.
!
        Vinfo( 1)='zgrad0'
        Vinfo( 2)='initial gradient for minimization, observation space'
#  ifdef RPCG
        vardim(1)=datumDimp
#  else
        vardim(1)=datumDim
#  endif
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef RPCG
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Define Lanczos vectors in observation space.
!
#  ifdef RPCG
        Vinfo( 1)='vcglwk'
        Vinfo( 2)='Preconditioned Lanczos vectors, observation space'
        vardim(1)=datumDimp
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='vcglev'
        Vinfo( 2)='converged Lanczos vectors, observation space'
        vardim(1)=datumDimp
        vardim(2)=NinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='zcglwk'
        Vinfo( 2)='Lanczos vectors, observation space'
        vardim(1)=datumDimp
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  else
        Vinfo( 1)='zcglwk'
        Vinfo( 2)='Lanczos vectors, observation space'
        vardim(1)=datumDim
        vardim(2)=MinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Define previous values of TLmodVal.
!
        Vinfo( 1)='TLmodVal_S'
        Vinfo( 2)='tangent linear model at observation locations'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        vardim(3)=NouterDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SP4DVAR
!
!  Define saved values of ADmodVal.
!
        Vinfo( 1)='ADmodVal_S'
        Vinfo( 2)='adjoint model at observation locations'
        vardim(1)=datumDim
        vardim(2)=NinnerDim2
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define saved values of harnoldi.
!
        Vinfo( 1)='harnoldi'
        Vinfo( 2)='GMRES upper Hessenberg matrix'
        vardim(1)=NinnerDim1
        vardim(2)=NinnerDim
        vardim(3)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 3, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define saved values of gmze, the Arnoldi vector weights.
!
        Vinfo( 1)='gmze'
        Vinfo( 2)='GMRES Arnoldi vector weights'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define saved values of cg_beta0.
!
        Vinfo( 1)='cg_beta0'
        Vinfo( 2)='Residual norm'
        vardim(1)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define saved values of Jobs.
!
        Vinfo( 1)='Jobs'
        Vinfo( 2)='Obs cost function'
        vardim(1)=NinnerDim1
        vardim(2)=NouterDim
        Vinfo(24)='_FillValue'
        Aval(6)=spval
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
     defined TL_R4DVAR
!
!  Define initial values of RPmodVal.
!
        Vinfo( 1)='RPmodel_initial'
        Vinfo( 2)='initial representer model at observation locations'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined I4DVAR_ANA_SENSITIVITY   || \
     defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined R4DVAR_ANA_SENSITIVITY
#  ifdef OBS_IMPACT
!
!  Define total observations impact.
!
#   ifdef IMPACT_INNER
        Vinfo( 1)='ObsImpact_total'
        Vinfo( 2)='total observation impact'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   else
        Vinfo( 1)='ObsImpact_total'
        Vinfo( 2)='total observation impact'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  else
!
!  Define total observation sensitivity.
!
        Vinfo( 1)='ObsSens_total'
        Vinfo( 2)='total observation sensitivity'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined OBS_IMPACT_SPLIT          && \
     (defined I4DVAR_ANA_SENSITIVITY   || \
      defined RBL4DVAR_ANA_SENSITIVITY || \
      defined RBL4DVAR_FCT_SENSITIVITY || \
      defined R4DVAR_ANA_SENSITIVITY)
!
!  Define observation impact due to initial condition increments.
!
#  ifdef IMPACT_INNER
        Vinfo( 1)='ObsImpact_IC'
        Vinfo( 2)='observation impact due to initial conditions'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  else
        Vinfo( 1)='ObsImpact_IC'
        Vinfo( 2)='observation impact due to initial conditions'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
!
!  Define observation impact due to surface forcing increments.
!
#   ifdef IMPACT_INNER

        Vinfo( 1)='ObsImpact_FC'
        Vinfo( 2)='observation impact due to surface forcing'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   else
        Vinfo( 1)='ObsImpact_FC'
        Vinfo( 2)='observation impact due to surface forcing'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  if defined ADJUST_BOUNDARY
!
!  Define observation impact due to boundary condition increments.
!
#   ifdef IMPACT_INNER
        Vinfo( 1)='ObsImpact_BC'
        Vinfo( 2)='observation impact due to open boundary conditions'
        vardim(1)=datumDim
        vardim(2)=NinnerDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   else
        Vinfo( 1)='ObsImpact_BC'
        Vinfo( 2)='observation impact due to open boundary conditions'
        vardim(1)=datumDim
        status=def_var(ng, iNLM, DAV(ng)%ncid, varid, nf90_double,      &
     &                 1, vardim, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
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
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Close observations NetCDF.
!
        IF (OBS(ng)%ncid.eq.-1) THEN
          CALL netcdf_close (ng, iNLM, OBSncid, OBS(ng)%name, .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF DEFINE
!
!=======================================================================
!  Open an existing model/observation file and check its contents.
!=======================================================================
!
      QUERY : IF (.not.LdefMOD(ng)) THEN
        ncname=DAV(ng)%name
!
!  Open model/observation for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname,                        &
     &                         ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname,                          &
     &                       ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmp))) THEN
            got_var(idNLmp)=.TRUE.
            DAV(ng)%Vid(idNLmp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmu))) THEN
            got_var(idNLmu)=.TRUE.
            DAV(ng)%Vid(idNLmu)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmi))) THEN
            got_var(idNLmi)=.TRUE.
            DAV(ng)%Vid(idNLmi)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmf))) THEN
            got_var(idNLmf)=.TRUE.
            DAV(ng)%Vid(idNLmf)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNLmo))) THEN
            got_var(idNLmo)=.TRUE.
            haveNLmod(ng)=.TRUE.
            DAV(ng)%Vid(idNLmo)=var_id(i)
# if defined I4DVAR || defined WEAK_CONSTRAINT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTLmo))) THEN
            got_var(idTLmo)=.TRUE.
            haveTLmod(ng)=.TRUE.
            DAV(ng)%Vid(idTLmo)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBgEr))) THEN
            got_var(idBgEr)=.TRUE.
            DAV(ng)%Vid(idBgEr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idInno))) THEN
            got_var(idInno)=.TRUE.
            DAV(ng)%Vid(idInno)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIncr))) THEN
            got_var(idIncr)=.TRUE.
            DAV(ng)%Vid(idIncr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idResi))) THEN
            got_var(idResi)=.TRUE.
            DAV(ng)%Vid(idResi)=var_id(i)
# if defined I4DVAR || defined WEAK_CONSTRAINT
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
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idOday)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOpro)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idOpro)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsT)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idObsT)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOlon)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idOlon)),          &
     &                                  TRIM(DAV(ng)%name)
!!        exit_flag=2
!!        RETURN
        END IF
        IF (.not.got_var(idOlat)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idOlat)),          &
     &                                  TRIM(DAV(ng)%name)
!!        exit_flag=2
!!        RETURN
        END IF
        IF (.not.got_var(idObsD)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idObsD)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsX)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idObsX)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsY)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idObsY)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOerr)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idOerr)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOval)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idOval)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF

        IF (.not.got_var(idObsS)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idObsS)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# ifdef FOUR_DVAR
        IF (.not.got_var(idNLmp)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idNLmp)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF

        IF (.not.got_var(idNLmu)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idNLmu)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF

        IF (.not.got_var(idNLmi)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idNLmi)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF

        IF (.not.got_var(idNLmf)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idNLmf)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idNLmo)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idNLmo)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# if defined I4DVAR || defined WEAK_CONSTRAINT
        IF (.not.got_var(idTLmo)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idTLmo)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idBgEr)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idBgEr)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idInno)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idInno)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idIncr)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idIncr)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idResi)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idResi)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# if defined I4DVAR || defined WEAK_CONSTRAINT
        IF (.not.got_var(idMOMi)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idMOMi)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idMOMf)) THEN
          IF (Master) WRITE (stdout,90) TRIM(Vname(1,idMOMf)),          &
     &                                  TRIM(DAV(ng)%name)
          exit_flag=2
          RETURN
        END IF
# endif

# if defined I4DVAR || defined BACKGROUND
!
!  Set unlimited time record dimension to the appropriate value.
!
        DAV(ng)%Rindex=rec_size
        Fcount=DAV(ng)%Fcount
        DAV(ng)%Nrec(Fcount)=rec_size
# endif

      END IF QUERY

  10  FORMAT (/,6x,'DEF_MOD   - creating model/observation data',       &
     &        ' file, Grid ',i2.2,': ',a)
  20  FORMAT (/,6x,'DEF_MOD   - inquiring model/observation data',      &
     &        ' file Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_MOD - unable to open observation/model file: ',a)
  40  FORMAT (/,' DEF_MOD - unable to create model/observation file:',  &
     &        1x,a)
  50  FORMAT (/,' DEF_MOD - unable to create global attribute: ',       &
     &        a,/,11x,a)
  60  FORMAT (/,' DEF_MOD - unable to copy attribute; ',a,1x,           &
     &        'for variable: ',a,/,11x,'in file: 'a)
  70  FORMAT (/,' DEF_MOD - unable to inquire attribute ',i2.2,1x,      &
     &        'name for variable: ',a,/,11x,'in file: 'a)
  80  FORMAT (/,' DEF_MOD - unable to copy attribute: ',1x,a,2x,        &
     &        'for variable: ',1x,a,/,11x,a)
  90  FORMAT (/,' DEF_MOD - unable to find model/observation variable:',&
     &        1x,a,/,11x,'in file: ',a)
!
      RETURN
      END SUBROUTINE def_mod
#else
      SUBROUTINE def_mod
      RETURN
      END SUBROUTINE def_mod
#endif
