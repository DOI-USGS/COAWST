#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE stats_modobs (ng)
!
!svn $Id: stats_modobs.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      integer :: i, ic, iobs, status, varid, vindex
      integer :: LBi, UBi, LBj, UBj
      integer, dimension(2) :: start, total
      integer, dimension(0:NstateVar(ng)) :: Ncount, is

      integer, dimension(Ndatum(ng)) :: obs_type
      integer, dimension(Ndatum(ng)) :: provenance

      real(r8) :: cff1, cff2
      real(r8), parameter :: IniVal = 0.0_r8
      real(r8), parameter :: LARGE = 1.0e+10_r8

      real(r8), dimension(0:NstateVar(ng)) :: CC, MB, MSE, SDE
      real(r8), dimension(0:NstateVar(ng)) :: mod_min, mod_max
      real(r8), dimension(0:NstateVar(ng)) :: mod_mean, mod_std
      real(r8), dimension(0:NstateVar(ng)) :: obs_min, obs_max
      real(r8), dimension(0:NstateVar(ng)) :: obs_mean, obs_std

      real(r8), dimension(Ndatum(ng)) :: mod_value
      real(r8), dimension(Ndatum(ng)) :: obs_depths
      real(r8), dimension(Ndatum(ng)) :: obs_scale
      real(r8), dimension(Ndatum(ng)) :: obs_value
      real(r8), dimension(Ndatum(ng)) :: obs_Xgrid
      real(r8), dimension(Ndatum(ng)) :: obs_Ygrid
      real(r8), dimension(Ndatum(ng)) :: obs_work
      real(r8), dimension(Nsurvey(ng)) :: survey

      character (len=11), dimension(NstateVar(ng)) :: text, svar_name
!
      SourceFile='stats_modobs.F'
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
!  Inquire about input observations variables.
!
      CALL netcdf_inq_var (ng, iNLM, OBS(ng)%name)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out survey time.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOday),    &
     &                      survey,                                     &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Nsurvey(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOday),    &
     &                      survey, (/1/), (/Nsurvey(ng)/),             &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out observation type identifier.
!
      CALL netcdf_get_ivar (ng, iNLM, OBS(ng)%name, Vname(1,idOtyp),    &
     &                      obs_type,                                   &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, Vname(1,idOtyp),    &
     &                      obs_type, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out observation provenance.
!
      CALL netcdf_get_ivar (ng, iNLM, OBS(ng)%name, Vname(1,idOpro),    &
     &                      provenance,                                 &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, Vname(1,idOpro),    &
     &                      provenance, (/1/), (/Ndatum(ng)/),          &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out observation time.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsT),    &
     &                      obs_work,                                   &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsT),    &
     &                      obs_work, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out observation longitude.
!
      IF (find_string(var_name,n_var,Vname(1,idOlon),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOlon),  &
     &                        obs_work,                                 &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (exit_flag.ne. NoError) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOlon),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Read in and write out observation latitude.
!
      IF (find_string(var_name,n_var,Vname(1,idOlon),vindex)) THEN
        CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOlat),  &
     &                        obs_work,                                 &
     &                        ncid = OBS(ng)%ncid,                      &
     &                        start = (/1/), total = (/Ndatum(ng)/))
        IF (exit_flag.ne. NoError) RETURN

        CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOlat),  &
     &                        obs_work, (/1/), (/Ndatum(ng)/),          &
     &                        ncid = DAV(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Read in and write out observation depth.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsD),    &
     &                      obs_depths,                                 &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN
!
!  Read in and write out X-fractional coordinate.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsX),    &
     &                      obs_Xgrid,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsX),    &
     &                      obs_Xgrid, (/1/), (/Ndatum(ng)/),           &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out Y-fractional coordinate.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idObsY),    &
     &                      obs_Ygrid,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsY),    &
     &                      obs_Ygrid, (/1/), (/Ndatum(ng)/),           &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out observations total error
!  (instrument + sampling + representation).
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOerr),    &
     &                      obs_work,                                   &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOerr),    &
     &                      obs_work, (/1/), (/Ndatum(ng)/),            &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in and write out observation values.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOval),    &
     &                      obs_value,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idOval),    &
     &                      obs_value, (/1/), (/Ndatum(ng)/),           &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Read in observation screening flag.
!
      CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idObsS),    &
     &                      obs_scale,                                  &
     &                      ncid = DAV(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN
!
!  Read in model values at observation locations.
!
# if defined TL_W4DVAR          || defined W4DVAR || \
     defined W4DVAR_SENSITIVITY
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
     &                      start = (/1,outer/),                        &
     &                      total = (/Ndatum(ng),1/))
#  endif
# endif
      IF (exit_flag.ne. NoError) RETURN
!
!-----------------------------------------------------------------------
!  If appropriate, convert depth in vertical fractional coordinates to
!  actual depths in meted (negative: downwards).
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
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Compute model and observations comparison statistics.
!-----------------------------------------------------------------------
!
!  Initialize.
!
      DO i=0,NstateVar(ng)
        CC(i)=0.0_r8
        MB(i)=0.0_r8
        MSE(i)=0.0_r8
        SDE(i)=0.0_r8
        mod_min(i)=LARGE
        mod_max(i)=-LARGE
        mod_mean(i)=0.0_r8
        obs_min(i)=LARGE
        obs_max(i)=-LARGE
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
          i=obs_type(iobs)
          Ncount(i)=Ncount(i)+1
          mod_min(i)=MIN(mod_min(i),mod_value(iobs))
          obs_min(i)=MIN(obs_min(i),obs_value(iobs))
          mod_max(i)=MAX(mod_max(i),mod_value(iobs))
          obs_max(i)=MAX(obs_max(i),obs_value(iobs))
          mod_mean(i)=mod_mean(i)+mod_value(iobs)
          obs_mean(i)=obs_mean(i)+obs_value(iobs)
        END IF
      END DO
      DO i=1,NstateVar(ng)
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
          i=obs_type(iobs)
          cff1=mod_value(iobs)-mod_mean(i)
          cff2=obs_value(iobs)-obs_mean(i)
          mod_std(i)=mod_std(i)+cff1*cff1
          obs_std(i)=obs_std(i)+cff2*cff2
          CC(i)=CC(i)+cff1*cff2
        END IF
      END DO
      DO i=1,NstateVar(ng)
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
      DO i=1,NstateVar(ng)
        IF (Ncount(i).gt.0) THEN
          MB(i)=mod_mean(i)-obs_mean(i)
          SDE(i)=mod_std(i)-obs_std(i)
          MSE(i)=MB(i)*MB(i)+                                           &
     &           SDE(i)*SDE(i)+                                         &
     &           2.0_r8*mod_std(i)*obs_std(i)*(1.0_r8-CC(i))
        END IF
      END DO
!
!  Report model and observations comparison statistics.
!
      IF (Master) THEN
        ic=0
        DO i=1,NstateVar(ng)
          svar_name(i)='           '
          text(i)='           '
          IF (Ncount(i).gt.0) THEN
            ic=ic+1
            is(ic)=i
            svar_name(ic)=TRIM(Vname(1,idSvar(i)))
            text(ic)='-----------'
          END IF
        END DO
        WRITE (stdout,10)
        WRITE (stdout,20) (svar_name(i),i=1,ic)
        WRITE (stdout,30) (text(i),i=1,ic)
        WRITE (stdout,40) 'Observation Min   ', (obs_min (is(i)),i=1,ic)
        WRITE (stdout,40) 'Observation Max   ', (obs_max (is(i)),i=1,ic)
        WRITE (stdout,40) 'Observation Mean  ', (obs_mean(is(i)),i=1,ic)
        WRITE (stdout,40) 'Observation STD   ', (obs_std (is(i)),i=1,ic)
        WRITE (stdout,40) 'Model Min         ', (mod_min (is(i)),i=1,ic)
        WRITE (stdout,40) 'Model Max         ', (mod_max (is(i)),i=1,ic)
        WRITE (stdout,40) 'Model Mean        ', (mod_mean(is(i)),i=1,ic)
        WRITE (stdout,40) 'Model STD         ', (mod_std (is(i)),i=1,ic)
        WRITE (stdout,40) 'Model Bias        ', (MB(is(i)),i=1,ic)
        WRITE (stdout,40) 'STD Error         ', (SDE(is(i)),i=1,ic)
        WRITE (stdout,40) 'Cross-Correlation ', (CC(is(i)),i=1,ic)
        WRITE (stdout,40) 'Mean Squared Error', (MSE(is(i)),i=1,ic)
        WRITE (stdout,50) 'Observation Count ', (Ncount(is(i)),i=1,ic)
      END IF
!
!  Write comparison statistics to NetCDF file.
!
      CALL netcdf_put_ivar (ng, iNLM, DAV(ng)%name, 'Nobs',             &
     &                      Ncount(0:), (/1/), (/NstateVar(ng)+1/),     &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'obs_mean',         &
     &                      obs_mean(0:), (/1/), (/NstateVar(ng)+1/),   &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'obs_std',          &
     &                      obs_std(0:), (/1/), (/NstateVar(ng)+1/),    &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'model_mean',       &
     &                      mod_mean(0:), (/1/), (/NstateVar(ng)+1/),   &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'model_std',        &
     &                      mod_std(0:), (/1/), (/NstateVar(ng)+1/),    &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'model_bias',       &
     &                      MB(0:), (/1/), (/NstateVar(ng)+1/),         &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'SDE',              &
     &                      SDE(0:), (/1/), (/NstateVar(ng)+1/),        &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'CC',               &
     &                      CC(0:), (/1/), (/NstateVar(ng)+1/),         &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_fvar (ng, iNLM, DAV(ng)%name, 'MSE',              &
     &                      MSE(0:), (/1/), (/NstateVar(ng)+1/),        &
     &                      ncid = DAV(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Synchronize NetCDF file to disk.
!
      CALL netcdf_sync (ng, iNLM, DAV(ng)%name, DAV(ng)%ncid)
!
 10   FORMAT (/,' Model-Observations Comparison Statistics:',/)
 20   FORMAT (t22,5(a11,1x))
 30   FORMAT (t22,5(a11,1x),/)
 40   FORMAT (a,3x,5(1p,e11.4,0p,1x))
 50   FORMAT (a,3x,5(i11,1x))
 60   FORMAT (/,' STATS_4DVAR - unable to synchronize 4DVAR',           &
     &        1x,'NetCDF file to disk.')

      RETURN
      END SUBROUTINE stats_modobs
#else
      SUBROUTINE stats_modobs
      RETURN
      END SUBROUTINE stats_modobs
#endif

