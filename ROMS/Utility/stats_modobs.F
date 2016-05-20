#include "cppdefs.h"
#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
      SUBROUTINE stats_modobs (ng)
!
!svn $Id$
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
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: i, ic, iobs, status, varid
      integer, dimension(2) :: start, total
      integer, dimension(0:NstateVar(ng)) :: Ncount, is

      integer, allocatable :: Obs_type(:)

      real(r8) :: cff1, cff2
      real(r8), parameter :: LARGE = 1.0e+10_r8

      real(r8), dimension(0:NstateVar(ng)) :: CC, MB, MSE, SDE
      real(r8), dimension(0:NstateVar(ng)) :: mod_min, mod_max
      real(r8), dimension(0:NstateVar(ng)) :: mod_mean, mod_std
      real(r8), dimension(0:NstateVar(ng)) :: obs_min, obs_max
      real(r8), dimension(0:NstateVar(ng)) :: obs_mean, obs_std

      real(r8), allocatable :: mod_value(:)
      real(r8), allocatable :: Obs_scale(:)
      real(r8), allocatable :: Obs_value(:)

      character (len=11), dimension(NstateVar(ng)) :: text, svar_name
!
      SourceFile='stats_modobs.F'
!
!-----------------------------------------------------------------------
!  Read in model and observations data.
!-----------------------------------------------------------------------
!
!  Allocate working arrays.
!
      IF (.not.allocated(mod_value)) allocate (mod_value(Ndatum(ng)))
      IF (.not.allocated(obs_scale)) allocate (obs_scale(Ndatum(ng)))
      IF (.not.allocated(obs_type))  allocate (obs_type (Ndatum(ng)))
      IF (.not.allocated(obs_value)) allocate (obs_value(Ndatum(ng)))
!
!  Read in observation type identifier.
!
      CALL netcdf_get_ivar (ng, iNLM, OBS(ng)%name, Vname(1,idOtyp),    &
     &                      obs_type,                                   &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN
!
!  Read in observation values.
!
      CALL netcdf_get_fvar (ng, iNLM, OBS(ng)%name, Vname(1,idOval),    &
     &                      obs_value,                                  &
     &                      ncid = OBS(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
      IF (exit_flag.ne. NoError) RETURN
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
      CALL netcdf_get_fvar (ng, iNLM, DAV(ng)%name, Vname(1,idNLmo),    &
     &                      mod_value,                                  &
     &                      ncid = DAV(ng)%ncid,                        &
     &                      start = (/1/), total = (/Ndatum(ng)/))
# endif
      IF (exit_flag.ne. NoError) RETURN
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
!  Deallocate working arrays.
!
      IF (allocated(mod_value)) deallocate (mod_value)
      IF (allocated(obs_scale)) deallocate (obs_scale)
      IF (allocated(obs_type))  deallocate (obs_type)
      IF (allocated(obs_value)) deallocate (obs_value)
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
 10   FORMAT (/,' 4DVAR Model-Observations Comparison Statistics:',/)
 20   FORMAT (t22,5(a11,1x))
 30   FORMAT (t22,5(a11,1x),/)
 40   FORMAT (a,3x,5(1p,e11.4,0p,1x))
 50   FORMAT (a,3x,5(i11,1x))
 60   FORMAT (/,' STATS_4DVAR - unable to synchronize 4DVAR',           &
     &        1x,'NetCDF file to disk.')
#else
      SUBROUTINE stats_modobs
#endif
      RETURN
      END SUBROUTINE stats_modobs
