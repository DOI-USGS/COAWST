#include "cppdefs.h"
#if defined FOUR_DVAR || defined VERIFICATION
      SUBROUTINE def_mod (ng)
!
!svn $Id: def_mod.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine create 4DVAR output NetCDF which contains model fields !
!  processed at observations points.                                   !
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
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, dimension(NV) :: got_var(NV)

      integer, parameter :: Natt = 24

      integer :: iterDim, datumDim, recordDim, stateDim
# ifdef IS4DVAR
      integer :: MinnerDim, NinnerDim, NouterDim, threeDim
# endif
      integer :: i, j, nvd, recdim, status, varid
      integer :: CostDim(2), Vsize(4), vardim(2)

      integer :: def_var

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
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
          DO j=1,80
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!=======================================================================
!  Create a new 4DVAR model data NetCDF file.
!=======================================================================
!
      ncname=MODname(ng)
      IF (LdefMOD(ng).and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncMODid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Define dimensions.
!-----------------------------------------------------------------------
!
      IF (LdefMOD(ng).and.OutThread) THEN
        status=nf90_def_dim(ncMODid(ng),'record',2,                     &
     &                      recordDim)
        status=nf90_def_dim(ncMODid(ng),'state_var',NstateVar(ng)+1,    &
     &                      stateDim)
        status=nf90_def_dim(ncMODid(ng),'datum'    ,Ndatum(ng),         &
     &                      datumDim)
# ifdef FOUR_DVAR
#  ifdef IS4DVAR
        status=nf90_def_dim(ncMODid(ng),'Ninner'   ,Ninner,             &
     &                      NinnerDim)
        status=nf90_def_dim(ncMODid(ng),'Minner'   ,Ninner+1,           &
     &                      MinnerDim)
        status=nf90_def_dim(ncMODid(ng),'Nouter'   ,Nouter,             &
     &                      NouterDim)
#   ifdef LANCZOS
        status=nf90_def_dim(ncMODid(ng),'three '   ,3,                  &
     &                      threeDim)
#   endif
#  endif
        status=nf90_def_dim(ncMODid(ng),'iteration',nf90_unlimited,     &
     &                      iterDim)
!
        CostDim(1)=stateDim
        CostDim(2)=iterDim
# endif
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
!  File type.
!
        type='ROMS/TOMS 4DVAR output observation processing file'
        status=nf90_put_att(ncMODid(ng),nf90_global,'type',             &
     &                      TRIM(type))
!
!  Input observations file.
!
        status=nf90_put_att(ncMODid(ng),nf90_global,'obs_file',         &
     &                      TRIM(OBSname(ng)))
!
!  SVN repository information.
!
      status=nf90_put_att(ncMODid(ng),nf90_global,'svn_url',            &
     &                    TRIM(svn_url))
# ifdef SVN_REV
      status=nf90_put_att(ncMODid(ng),nf90_global,'svn_rev',            &
     &                    TRIM(svn_rev))
# endif
# ifdef ROOT_DIR
      status=nf90_put_att(ncMODid(ng),nf90_global,'code_dir',           &
     &                    TRIM(Rdir))
# endif
# ifdef HEADER_DIR
      status=nf90_put_att(ncMODid(ng),nf90_global,'header_dir',         &
     &                    TRIM(Hdir))
# endif
# ifdef ROMS_HEADER
      status=nf90_put_att(ncMODid(ng),nf90_global,'header_file',        &
     &                    TRIM(Hfile))
# endif
!
!  Attributes describing platform and compiler
!
        status=nf90_put_att(ncMODid(ng),nf90_global,'os',               &
     &                      TRIM(my_os))
        status=nf90_put_att(ncMODid(ng),nf90_global,'cpu',              &
     &                      TRIM(my_cpu))
        status=nf90_put_att(ncMODid(ng),nf90_global,'compiler_system',  &
     &                      TRIM(my_fort))
        status=nf90_put_att(ncMODid(ng),nf90_global,'compiler_command', &
     &                      TRIM(my_fc))
        status=nf90_put_att(ncMODid(ng),nf90_global,'compiler_flags',   &
     &                      TRIM(my_fflags))
!
!  History attribute.
!
        IF (LEN_TRIM(date_str).gt.0) THEN
          WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',       &
     &                                    TRIM( version), TRIM(date_str)
        ELSE
          WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version', TRIM(version)
        END IF
        status=nf90_put_att(ncMODid(ng),nf90_global,'history',          &
     &                      TRIM(history))
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------

# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
!
!  Outer and inner loop contours.
!
        Vinfo( 1)='outer'
        Vinfo( 2)='outer loop counter'
        status=def_var(ncMODid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='inner'
        Vinfo( 2)='inner loop counter'
        status=def_var(ncMODid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)
# endif
# ifndef OBS_SENSITIVITY
!
!  Define model-observations comparison statistics.
!
        Vinfo( 1)='Nobs'
        Vinfo( 2)='Number of usable observations'
        status=def_var(ncMODid(ng),varid,nf90_int,1,stateDim,           &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='obs_mean'
        Vinfo( 2)='observations mean'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='obs_std'
        Vinfo( 2)='observations standard deviation'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='model_mean'
        Vinfo( 2)='model mean'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='model_std'
        Vinfo( 2)='model standard deviation'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='model_bias'
        Vinfo( 2)='model bias'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='SDE'
        Vinfo( 2)='model-observations standard deviation error'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='CC'
        Vinfo( 2)='model-observations cross-correlation'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,            &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='MSE'
        Vinfo( 2)='model-observations mean squared error'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,stateDim,Aval,       &
     &                 Vinfo,ncname)
# endif
# if defined IS4DVAR
!
!  Converged Ritz eigenvalues.
!
        Vinfo( 1)='nConvRitz'
        Vinfo( 2)='Number of converged Ritz eigenvalues'
        status=def_var(ncMODid(ng),varid,nf90_int,1,0,                  &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='Ritz'
        Vinfo( 2)='converged Ritz eigenvalues to approximate Hessian'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,Ninnerdim,           &
     &                 Aval,Vinfo,ncname)
!
!  Define conjugate gradient norms.
!
#  ifndef LANCZOS
        Vinfo( 1)='cg_alpha'
        Vinfo( 2)='conjugate gradient alpha coefficient'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
#  endif
        Vinfo( 1)='cg_beta'
        Vinfo( 2)='conjugate gradient beta coefficient'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='cg_tau'
        Vinfo( 2)='conjugate gradient tau coefficient'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)

#  ifdef LANCZOS
!
!  Define Lanczos algorithm coefficients.
!
        Vinfo( 1)='cg_delta'
        Vinfo( 2)='Lanczos algorithm delta coefficient'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='cg_gamma'
        Vinfo( 2)='Lanczos algorithm gamma coefficient'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Initial gradient vector normalization factor.
!
        Vinfo( 1)='cg_Gnorm'
        Vinfo( 2)='initial gradient normalization factor'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,NouterDim,           &
     &                 Aval,Vinfo,ncname)
!
!  Lanczos vector normalization factor.
!
        Vinfo( 1)='cg_QG'
        Vinfo( 2)='Lanczos vector normalization factor'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Reduction in the gradient norm.
!
        Vinfo( 1)='cg_Greduc'
        Vinfo( 2)='reduction in the gradient norm'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Lanczos recurrence tridiagonal matrix.
!
        Vinfo( 1)='cg_Tmatrix'
        Vinfo( 2)='Lanczos recurrence tridiagonal matrix'
        vardim(1)=NinnerDim
        vardim(2)=threeDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Lanczos tridiagonal matrix, upper diagonal elements.
!
        Vinfo( 1)='cg_zu'
        Vinfo( 2)='tridiagonal matrix, upper diagonal elements'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Eigenvalues of Lanczos recurrence relationship.
!
        Vinfo( 1)='cg_Ritz'
        Vinfo( 2)='Lanczos recurrence eigenvalues'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Eigenvalues relative error.
!
        Vinfo( 1)='cg_RitzErr'
        Vinfo( 2)='Ritz eigenvalues relative error'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
!
!  Eigenvectors of Lanczos recurrence relationship.
!
        Vinfo( 1)='cg_zv'
        Vinfo( 2)='Lanczos recurrence eigenvectors'
        vardim(1)=NinnerDim
        vardim(2)=NinnerDim
        status=def_var(ncMODid(ng),varid,NF_FRST,2,vardim,              &
     &                 Aval,Vinfo,ncname)
#  endif
# endif
# if defined WEAK_CONSTRAINT || defined IOM
!
!  Define conjugate gradient norms
!
        Vinfo( 1)='cg_gamma'
        Vinfo( 2)='conjugate gradient gamma norm'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,0,                   &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='cg_sigma'
        Vinfo( 2)='conjugate gradient sigma norm'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,0,                   &
     &                 Aval,Vinfo,ncname)

        Vinfo( 1)='cg_rnorm'
        Vinfo( 2)='conjugate gradient right-hand-side norm'
        status=def_var(ncMODid(ng),varid,NF_FRST,1,0,                   &
     &                 Aval,Vinfo,ncname)
# endif
!
!  Observations screening/normalization scale.
!
        Vinfo( 1)=Vname(1,idObsS)
        Vinfo( 2)=Vname(2,idObsS)
        status=def_var(ncMODid(ng),modVid(idObsS,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
# if defined FOUR_DVAR && !defined OBS_SENSITIVITY
!
!  Initial nonlinear model at observation locations.
!
        Vinfo( 1)=Vname(1,idNLmi)
        Vinfo( 2)=Vname(2,idNLmi)
        Vinfo(24)='_FillValue'
        Aval(6)=1.0E-6_r8
        status=def_var(ncMODid(ng),modVid(idNLmi,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
# endif
# ifndef OBS_SENSITIVITY
!
!  Nonlinear model at observation points.
!
        haveNLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idNLmo)
        Vinfo( 2)=Vname(2,idNLmo)
        Vinfo(24)='_FillValue'
        Aval(6)=1.0E-6_r8
#  ifdef S4DVAR
        vardim(1)=datumDim
        vardim(2)=recordDim
        status=def_var(ncMODid(ng),modVid(idNLmo,ng),NF_FRST,           &
     &                 2,vardim,Aval,Vinfo,ncname)
#  else
        status=def_var(ncMODid(ng),modVid(idNLmo,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
#  endif
# endif
# if defined IS4DVAR         || defined IS4DVAR_OLD     || \
     defined IOM             || defined OBS_SENSITIVITY || \
     defined WEAK_CONSTRAINT 
!
!  Tangent linear or representer model at observation points.
!
        haveTLmod(ng)=.FALSE.
        Vinfo( 1)=Vname(1,idTLmo)
#  ifdef OBS_SENSITIVITY
        Vinfo( 2)='4DVAR sensitivity analysis at observations location'
#  else
        Vinfo( 2)=Vname(2,idTLmo)
#  endif
        Vinfo(24)='_FillValue'
        Aval(6)=1.0E-6_r8
#  ifdef IS4DVAR_OLD
        vardim(1)=datumDim
        vardim(2)=recordDim
        status=def_var(ncMODid(ng),modVid(idTLmo,ng),NF_FRST,           &
     &                 2,vardim,Aval,Vinfo,ncname)
#  else
        status=def_var(ncMODid(ng),modVid(idTLmo,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
#  endif
# endif
# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
!
!  Initial model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMi)
        Vinfo( 2)=Vname(2,idMOMi)
        Vinfo(24)='_FillValue'
        Aval(6)=0.0_r8
        status=def_var(ncMODid(ng),modVid(idMOMi,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
!
!  Final model-observation misfit (innovation) vector.
!
        Vinfo( 1)=Vname(1,idMOMf)
        Vinfo( 2)=Vname(2,idMOMf)
        Vinfo(24)='_FillValue'
        Aval(6)=0.0_r8
        status=def_var(ncMODid(ng),modVid(idMOMf,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
# endif
# if defined WEAK_CONSTRAINT || defined IOM
!
!  Representer coeficients (or their approximation PSI) at observation
!  points.
!
        Vinfo( 1)=Vname(1,idRepC)
        Vinfo( 2)=Vname(2,idRepC)
        status=def_var(ncMODid(ng),modVid(idRepC,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
!
!  Conjugate gradient vectors.
!
        Vinfo( 1)=Vname(1,idCG_P)
        Vinfo( 2)=Vname(2,idCG_P)
        status=def_var(ncMODid(ng),modVid(idCG_P,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)

        Vinfo( 1)=Vname(1,idCG_R)
        Vinfo( 2)=Vname(2,idCG_R)
        status=def_var(ncMODid(ng),modVid(idCG_R,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)

        Vinfo( 1)=Vname(1,idCG_S)
        Vinfo( 2)=Vname(2,idCG_S)
        status=def_var(ncMODid(ng),modVid(idCG_S,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)

        Vinfo( 1)=Vname(1,idCG_V)
        Vinfo( 2)=Vname(2,idCG_V)
        status=def_var(ncMODid(ng),modVid(idCG_V,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
       
        Vinfo( 1)=Vname(1,idCG_X)
        Vinfo( 2)=Vname(2,idCG_X)
        status=def_var(ncMODid(ng),modVid(idCG_X,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)

        Vinfo( 1)=Vname(1,idCG_Z)
        Vinfo( 2)=Vname(2,idCG_Z)
        status=def_var(ncMODid(ng),modVid(idCG_Z,ng),NF_FRST,           &
     &                 1,datumDim,Aval,Vinfo,ncname)
# endif
# if defined IS4DVAR || defined IS4DVAR_OLD
!
!  Define cost function gradient norm.
!
        Vinfo( 1)='cost_gradient'
        Vinfo( 2)='cost function gradient norm'
        Vinfo(17)='missing_value'
        Aval(1)=0.0_r8
        status=def_var(ncMODid(ng),varid,nf90_double,                   &
     &                 2,CostDim,Aval,Vinfo,ncname)
# endif
# if defined IS4DVAR || defined IS4DVAR_OLD || defined S4DVAR
!
!  Define model minus observations misfit cost function.
!
        Vinfo( 1)='cost_function'
        Vinfo( 2)='model minus observation misfit cost function'
        Vinfo(17)='missing_value'
        Aval(1)=0.0_r8
        status=def_var(ncMODid(ng),varid,nf90_double,                   &
     &                 2,CostDim,Aval,Vinfo,ncname)
# endif
# ifdef BACKGROUND
!
!  Define model minus background misfit cost function.
!
        Vinfo( 1)='back_function'
        Vinfo( 2)='model minus background misfit cost function'
        Vinfo(17)='missing_value'
        Aval(1)=0.0_r8
        status=def_var(ncMODid(ng),varid,nf90_double,                   &
     &                 2,CostDim,Aval,Vinfo,ncname)
# endif
# if defined IS4DVAR || defined IS4DVAR_OLD
!
!  Define optimality property that measures the consistency between
!  background and observation errors hypotheses (Chi-square). 
!
        Vinfo( 1)='Jmin'
        Vinfo( 2)='normalized, optimal cost function minimum'
        status=def_var(ncMODid(ng),varid,nf90_double,                   &
     &                 1,iterDim,Aval,Vinfo,ncname)
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncMODid(ng))
      END IF
!
!=======================================================================
!  Open an existing 4DVAR file and check its contents.
!=======================================================================
!
      IF (.not.LdefMOD(ng).and.Master) THEN
        CALL opencdf (ng, 1, MODname(ng), ncname, N(ng), 0, recdim, nvd, &
     &                Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from observation NetCDF and activate switches for
!  required variables.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idObsS))) THEN
            got_var(idObsS)=.TRUE.
            modVid(idObsS,ng)=i
          END IF
# ifdef FOUR_DVAR
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idNLmi))) THEN
            got_var(idNLmi)=.TRUE.
            modVid(idNLmi,ng)=i
          END IF
# endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idNLmo))) THEN
            got_var(idNLmo)=.TRUE.
            haveNLmod(ng)=.TRUE.
            modVid(idNLmo,ng)=i
          END IF
# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTLmo))) THEN
            got_var(idTLmo)=.TRUE.
            haveTLmod(ng)=.TRUE.
            modVid(idTLmo,ng)=i
          END IF
# endif
# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMOMi))) THEN
            got_var(idMOMi)=.TRUE.
            modVid(idMOMi,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMOMf))) THEN
            got_var(idMOMf)=.TRUE.
            modVid(idMOMf,ng)=i
          END IF
# endif
# if defined WEAK_CONSTRAINT || defined IOM
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRepC))) THEN
            got_var(idRepC)=.TRUE.
            haveADmod(ng)=.TRUE.
            modVid(idRepC,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCG_P))) THEN
            got_var(idCG_P)=.TRUE.
            modVid(idCG_P,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCG_R))) THEN
            got_var(idCG_R)=.TRUE.
            modVid(idCG_R,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCG_S))) THEN
            got_var(idCG_S)=.TRUE.
            modVid(idCG_S,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCG_V))) THEN
            got_var(idCG_V)=.TRUE.
            modVid(idCG_V,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCG_X))) THEN
            got_var(idCG_X)=.TRUE.
            modVid(idCG_X,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idCG_Z))) THEN
            got_var(idCG_Z)=.TRUE.
            modVid(idCG_Z,ng)=i
          END IF
# endif
        END DO
!
!  Check if needed variables are available.
!
        IF (.not.got_var(idObsS)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idObsS)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# ifdef FOUR_DVAR
        IF (.not.got_var(idNLmi)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idNLmi)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
        IF (.not.got_var(idNLmo)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idNLmo)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
        IF (.not.got_var(idTLmo)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idTLmo)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
# if defined IS4DVAR         || defined IS4DVAR_OLD || \
     defined WEAK_CONSTRAINT || defined IOM
        IF (.not.got_var(idMOMi)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idMOMi)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idMOMf)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idMOMf)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
# if defined WEAK_CONSTRAINT || defined IOM
        IF (.not.got_var(idRepC)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idRepC)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idCG_P)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCG_P)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idCG_R)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCG_R)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idCG_S)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCG_S)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idCG_V)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCG_V)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idCG_X)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCG_X)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idCG_Z)) THEN
          WRITE (stdout,30) TRIM(Vname(1,idCG_Z)), TRIM(MODname(ng))
          exit_flag=2
          RETURN
        END IF
# endif
      END IF

  10  FORMAT (/,6x,'DEF_MOD   - creating 4DVAR model data file: ', a)
  20  FORMAT (/,6x,'DEF_MOD   - inquiring 4DVAR model data file: ', a)
  30  FORMAT (/,' DEF_MOD - unable to find 4DVAR variable: ',a,         &
     &        /,11x,'in input NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_mod
#else
      SUBROUTINE def_mod
      RETURN
      END SUBROUTINE def_mod
#endif
