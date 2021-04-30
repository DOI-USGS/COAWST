#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE stats_modobs (ng)
!
!svn $Id: stats_modobs.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine computes several statistical quantities between     !
!  model and observations:                                             !
!                                                                      !
!     CC         Cross-Correlation                                     !
!     MB         Model Bias                                            !
!     MSE        Mean Squared Error                                    !
!     SDE        Standard Deviation Error                              !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Oke, P.R., J.S. Allen, R.N. Miller, G.D. Egbert, J.A. Austin,       !
!    J.A. Barth, T.J. Boyd, P.M. Kosro, and M.D. Levine, 2002: A       !
!    Modeling Study of the Three-Dimensional Continental Shelf         !
!    Circulation off Oregon. Part I: Model-Data Comparison, J.         !
!    Phys. Oceanogr., 32, 1360-1382.                                   !
!                                                                      !
!  Notice that this routine is never called inside of a parallel       !
!  region. Therefore, parallel reductions are not required.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY :  mp_collect
# endif
      USE obs_k2z_mod, ONLY : obs_k2z
      USE strings_mod, ONLY : FoundError, find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, ic, iobs, status, varid, vindex
      integer :: LBi, UBi, LBj, UBj
      integer :: saved_exit_flag

      integer, dimension(2) :: start, total
      integer, dimension(NobsVar(ng)) :: Ncount, is

      integer, dimension(Ndatum(ng)) :: obs_type
      integer, dimension(Ndatum(ng)) :: provenance
      integer, dimension(Nsurvey(ng)) :: survey_Nobs
!
      real(r8) :: cff1, cff2
      real(r8), parameter :: IniVal = 0.0_r8

      real(r8), dimension(NobsVar(ng)) :: CC, MB, MSE, SDE
      real(r8), dimension(NobsVar(ng)) :: mod_min, mod_max
      real(r8), dimension(NobsVar(ng)) :: mod_mean, mod_std
      real(r8), dimension(NobsVar(ng)) :: obs_min, obs_max
      real(r8), dimension(NobsVar(ng)) :: obs_mean, obs_std

      real(r8), dimension(Ndatum(ng)) :: mod_value
      real(r8), dimension(Ndatum(ng)) :: obs_depths
      real(r8), dimension(Ndatum(ng)) :: obs_scale
      real(dp), dimension(Ndatum(ng)) :: obs_time
      real(r8), dimension(Ndatum(ng)) :: obs_value
      real(r8), dimension(Ndatum(ng)) :: obs_Xgrid
      real(r8), dimension(Ndatum(ng)) :: obs_Ygrid
      real(r8), dimension(Ndatum(ng)) :: obs_work
      real(dp), dimension(Nsurvey(ng)) :: survey
!
      character (len=11), dimension(NobsVar(ng)) :: text, svar_name

      character (len=30) :: F1, F2, F3, F4

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Read in model and observations data.
!-----------------------------------------------------------------------
!
!  If exit_flag > 0, save and reset to allow completion of the DAV
!  when ROMS blows-up.
!
      IF ((exit_flag.eq.1).or.(blowup.gt.0)) THEN
        saved_exit_flag=exit_flag;
        exit_flag=NoError
      ELSE
        saved_exit_flag=NoError;
      END IF
!
!  Inquire about input observations variables.
!
      CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name,                      &
     &                     ncid = OBS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out number of observations per survey.
!
      CALL netcdf_get_ivar (ng, iNLM, OBS(ng)%name, Vname(1,idNobs),    &
     &                      survey_Nobs,                                &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Nsurvey(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, Vname(1,idNobs),    &
     &                      survey_Nobs, (/1/), (/Nsurvey(ng)/),        &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out survey time.
!
      CALL netcdf_get_time (ng, iNLM, OBS(ng)%name, Vname(1,idOday),    &
     &                      Rclock%DateNumber, survey,                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Nsurvey(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOday),    &
     &                      survey, (/1/), (/Nsurvey(ng)/),             &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observation type identifier.
!
      CALL netcdf_get_ivar (ng, iNLM, OBS(ng)%name, Vname(1,idOtyp),    &
     &                      obs_type,                                   &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, Vname(1,idOtyp),    &
     &                      obs_type, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observation provenance.
!
      CALL netcdf_get_ivar (ng, iNLM, OBS(ng)%name, Vname(1,idOpro),    &
     &                      provenance,                                 &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, Vname(1,idOpro),    &
     &                      provenance, (/1/), (/Ndatum(ng)/),          &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observation time.
!
      CALL netcdf_get_time (ng, iNLM, OBS(ng)%name, Vname(1,idObsT),    &
     &                      Rclock%DateNumber, obs_time,                &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsT),    &
     &                      obs_time, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observation longitude.
!
      IF (find_string(var_name,n_var,Vname(1,idOlon),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOlon),  &
     &                        obs_work,                                 &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOlon),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Read in and write out observation latitude.
!
      IF (find_string(var_name,n_var,Vname(1,idOlon),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOlat),  &
     &                        obs_work,                                 &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOlat),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Read in observation depth.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsD),    &
     &                      obs_depths,                                 &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out X-fractional coordinate.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsX),    &
     &                      obs_Xgrid,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsX),    &
     &                      obs_Xgrid, (/1/), (/Ndatum(ng)/),           &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out Y-fractional coordinate.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsY),    &
     &                      obs_Ygrid,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsY),    &
     &                      obs_Ygrid, (/1/), (/Ndatum(ng)/),           &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observations total error covariance
!  (instrument + sampling + representation).
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOerr),    &
     &                      obs_work,                                   &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOerr),    &
     &                      obs_work, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observation values.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOval),    &
     &                      obs_value,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOval),    &
     &                      obs_value, (/1/), (/Ndatum(ng)/),           &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in and write out observation meta values.
!
      IF (haveObsMeta(ng)) THEN
        CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOmet),  &
     &                        obs_work,                                 &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOmet),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Read in observation screening flag.
!
      CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsS),    &
     &                      obs_scale,                                  &
     &                      ncid = DAV(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in model values at observation locations.
!
# if defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
     defined TL_R4DVAR
      CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idTLmo),    &
     &                      mod_value,                                  &
     &                      ncid = DAV(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
# else
#  ifdef VERIFICATION
      CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmo),    &
     &                      mod_value,                                  &
     &                      ncid = DAV(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
#  else
      CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmo),    &
     &                      mod_value,                                  &
     &                      ncid = DAV(ng)%ncid,                        &
     &                      start = (/1,MAX(1,outer)/),                 &
     &                      total = (/Ndatum(ng),1/))
#  endif
# endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  If appropriate, convert depth in vertical fractional coordinates to
!  actual depths in meters (negative: downwards).
!-----------------------------------------------------------------------
!
      obs_work(1:Ndatum(ng))=IniVal
!
      CALL obs_k2z (ng, 0, Lm(ng)+1, 0, Mm(ng)+1,                       &
     &              LBi, UBi, LBj, UBj, 0, N(ng),                       &
     &              rXmin(ng), rXmax(ng),                               &
     &              rYmin(ng), rYmax(ng),                               &
     &              Ndatum(ng),                                         &
     &              obs_Xgrid, obs_Ygrid, obs_depths, obs_scale,        &
     &              GRID(ng) % z0_w,                                    &
     &              obs_work)

# ifdef DISTRIBUTE
!
!  Collect all the computed depths.
!
      CALL mp_collect (ng, iNLM, Ndatum(ng), IniVal, obs_work)
# endif
!
!  Write out to depths to NetCDF file.
!
      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsD),    &
     &                      obs_work, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined RBL4DVAR_ANA_SENSITIVITY || \
     defined RBL4DVAR_FCT_SENSITIVITY || \
     defined R4DVAR_ANA_SENSITIVITY
!
!-----------------------------------------------------------------------
!  Write out several variables from the 4D-Var estimation for
!  completeness.
!-----------------------------------------------------------------------
!
!  Inquire about input observations variables.
!
      CALL netcdf_inq_var (ng, iNLM, LCZ(ng)%name)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  4D-Var innovation (observation minus background).
!
      IF (find_string(var_name,n_var,Vname(1,idInno),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idInno),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idInno),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  4D-Var increment (analysis minus background).
!
      IF (find_string(var_name,n_var,Vname(1,idIncr),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idIncr),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idIncr),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  4D-Var residual (observation minus analysis).
!
      IF (find_string(var_name,n_var,Vname(1,idResi),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idResi),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idResi),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  4D-Var initial model-observation misfit.
!
      IF (find_string(var_name,n_var,Vname(1,idMOMi),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idMOMi),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idMOMi),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  4D-Var final model-observation misfit.
!
      IF (find_string(var_name,n_var,Vname(1,idMOMf),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idMOMf),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idMOMf),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Unvetted initial NLM at observation locations.
!
      IF (find_string(var_name,n_var,Vname(1,idNLmp),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idNLmp),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmp),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Final NLM at observation locations.
!
      IF (find_string(var_name,n_var,Vname(1,idNLmf),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idNLmf),  &
     &                        obs_work,                                 &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmf),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Unvetted NLM at observation locations per outer-loop.
!
      IF (find_string(var_name,n_var,Vname(1,idNLmu),vindex)) THEN
        DO i=1,Nouter
          CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idNLmu),&
     &                          obs_work,                               &
     &                          start = (/1,i/),                        &
                                total = (/Ndatum(ng),1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmu),&
     &                          obs_work, (/1,i/), (/Ndatum(ng),1/),    &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
!
!  NLM at observation locations per outer-loop
!
      IF (find_string(var_name,n_var,Vname(1,idNLmo),vindex)) THEN
        DO i=1,Nouter
          CALL netcdf_get_fvar (ng, iNLM, LCZ(ng)%name, Vname(1,idNLmo),&
     &                          obs_work,                               &
     &                          start = (/1,i/),                        &
                                total = (/Ndatum(ng),1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmo),&
     &                          obs_work, (/1,i/), (/Ndatum(ng),1/),    &
     &                          ncid = DAV(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
      END IF
# endif

# if defined I4DVAR       || defined RBL4DVAR       || \
     defined R4DVAR       || defined SP4DVAR        || \
     defined SPLIT_I4DVAR || defined SPLIT_RBL4DVAR || \
     defined SPLIT_R4DVAR || defined SPLIT_SP4DVAR  || \
     defined VERIFICATION
!
!-----------------------------------------------------------------------
!  Compute model and observations comparison statistics.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      IF (saved_exit_flag.eq.NoError) THEN
        DO i=1,NobsVar(ng)
          CC(i)=0.0_r8
          MB(i)=0.0_r8
          MSE(i)=0.0_r8
          SDE(i)=0.0_r8
          mod_min(i)=Large
          mod_max(i)=-Large
          mod_mean(i)=0.0_r8
          obs_min(i)=Large
          obs_max(i)=-Large
          obs_mean(i)=0.0_r8
          mod_std(i)=0.0_r8
          obs_std(i)=0.0_r8
          Ncount(i)=0
        END DO
!
!  Compute model and observations mean per each state variable.
!
        DO iobs=1,Ndatum(ng)
          IF (obs_scale(iobs).eq.1.0_r8) THEN
            i=ObsType2State(obs_type(iobs))
            Ncount(i)=Ncount(i)+1
            mod_min(i)=MIN(mod_min(i),mod_value(iobs))
            obs_min(i)=MIN(obs_min(i),obs_value(iobs))
            mod_max(i)=MAX(mod_max(i),mod_value(iobs))
            obs_max(i)=MAX(obs_max(i),obs_value(iobs))
            mod_mean(i)=mod_mean(i)+mod_value(iobs)
            obs_mean(i)=obs_mean(i)+obs_value(iobs)
          END IF
        END DO
        DO i=1,NobsVar(ng)
          IF (Ncount(i).gt.0) THEN
            mod_mean(i)=mod_mean(i)/REAL(Ncount(i),r8)
            obs_mean(i)=obs_mean(i)/REAL(Ncount(i),r8)
          END IF
        END DO
!
!  Compute standard deviation and cross-correlation between model and
!  observations (CC).
!
        DO iobs=1,Ndatum(ng)
          IF (obs_scale(iobs).eq.1.0_r8) THEN
            i=ObsType2State(obs_type(iobs))
            cff1=mod_value(iobs)-mod_mean(i)
            cff2=obs_value(iobs)-obs_mean(i)
            mod_std(i)=mod_std(i)+cff1*cff1
            obs_std(i)=obs_std(i)+cff2*cff2
            CC(i)=CC(i)+cff1*cff2
          END IF
        END DO
        DO i=1,NobsVar(ng)
          IF (Ncount(i).gt.1) THEN
            mod_std(i)=SQRT(mod_std(i)/REAL(Ncount(i)-1,r8))
            obs_std(i)=SQRT(obs_std(i)/REAL(Ncount(i)-1,r8))
            CC(i)=(CC(i)/REAL(Ncount(i),r8))/(mod_std(i)*obs_std(i))
          END IF
        END DO
!
!  Compute model bias (MB), standard deviation error (SDE), and mean
!  squared error (MSE).
!
        DO i=1,NobsVar(ng)
          IF (Ncount(i).gt.0) THEN
            MB(i)=mod_mean(i)-obs_mean(i)
            SDE(i)=mod_std(i)-obs_std(i)
            MSE(i)=MB(i)*MB(i)+                                         &
     &             SDE(i)*SDE(i)+                                       &
     &             2.0_r8*mod_std(i)*obs_std(i)*(1.0_r8-CC(i))
          END IF
        END DO
!
!  Report model and observations comparison statistics.
!

        IF (Master) THEN
          ic=0
          DO i=1,NobsVar(ng)
            svar_name(i)='           '
            text(i)='           '
            IF (Ncount(i).gt.0) THEN
              ic=ic+1
              is(ic)=i
              svar_name(ic)=TRIM(ObsName(i))
              text(ic)='-----------'
            END IF
          END DO
          WRITE(F1,10) MAX(1,ic)
          WRITE(F2,20) MAX(1,ic)
          WRITE(F3,30) MAX(1,ic)
          WRITE(F4,40) MAX(1,ic)
          WRITE(stdout,50)
          WRITE(stdout,F1) (svar_name(i),i=1,ic)
          WRITE(stdout,F2) (text(i),i=1,ic)
          WRITE(stdout,F3) 'Observation Min   ',(obs_min (is(i)),i=1,ic)
          WRITE(stdout,F3) 'Observation Max   ',(obs_max (is(i)),i=1,ic)
          WRITE(stdout,F3) 'Observation Mean  ',(obs_mean(is(i)),i=1,ic)
          WRITE(stdout,F3) 'Observation STD   ',(obs_std (is(i)),i=1,ic)
          WRITE(stdout,F3) 'Model Min         ',(mod_min (is(i)),i=1,ic)
          WRITE(stdout,F3) 'Model Max         ',(mod_max (is(i)),i=1,ic)
          WRITE(stdout,F3) 'Model Mean        ',(mod_mean(is(i)),i=1,ic)
          WRITE(stdout,F3) 'Model STD         ',(mod_std (is(i)),i=1,ic)
          WRITE(stdout,F3) 'Model Bias        ',(MB(is(i)),i=1,ic)
          WRITE(stdout,F3) 'STD Error         ',(SDE(is(i)),i=1,ic)
          WRITE(stdout,F3) 'Cross-Correlation ',(CC(is(i)),i=1,ic)
          WRITE(stdout,F3) 'Mean Squared Error',(MSE(is(i)),i=1,ic)
          WRITE(stdout,F4) 'Observation Count ',(Ncount(is(i)),i=1,ic)
        END IF
!
!  Write comparison statistics to NetCDF file.
!
        CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, 'Nused_obs',      &
     &                        Ncount, (/1/), (/NobsVar(ng)/),           &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'obs_mean',       &
     &                        obs_mean, (/1/), (/NobsVar(ng)/),         &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'obs_std',        &
     &                        obs_std, (/1/), (/NobsVar(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'model_mean',     &
     &                        mod_mean, (/1/), (/NobsVar(ng)/),         &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'model_std',      &
     &                        mod_std, (/1/), (/NobsVar(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'model_bias',     &
     &                        MB, (/1/), (/NobsVar(ng)/),               &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'SDE',            &
     &                        SDE, (/1/), (/NobsVar(ng)/),              &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'CC',             &
     &                        CC, (/1/), (/NobsVar(ng)/),               &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'MSE',            &
     &                        MSE, (/1/), (/NobsVar(ng)/),              &
     &                        ncid = DAV(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
!
!  Synchronize NetCDF file to disk.
!
      CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Reset "exit_flag" to its saved value when entering this routine.
!
      exit_flag=saved_exit_flag
!
 10   FORMAT ('(t22,',i2,'(a11,1x))')
 20   FORMAT ('(t22,',i2,'(a11,1x),/)')
 30   FORMAT ('(a,3x,',i2,'(1p,e11.4,0p,1x))')
 40   FORMAT ('(a,3x,',i2,'(i11,1x))')
 50   FORMAT (/,' Model-Observations Comparison Statistics:',/)

      RETURN
      END SUBROUTINE stats_modobs
#else
      SUBROUTINE stats_modobs
      RETURN
      END SUBROUTINE stats_modobs
#endif

