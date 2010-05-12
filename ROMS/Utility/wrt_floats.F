#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE wrt_floats (ng)
!
!svn $Id: wrt_floats.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes simulated drifter trajectories into floats   !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: itrc, l, status

      real(r8), dimension(Nfloats(ng)) :: Tout
!
      SourceFile='wrt_floats.F'
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tFLTindx(ng)=tFLTindx(ng)+1
      NrecFLT(ng)=NrecFLT(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                      &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tFLTindx(ng)/), (/1/),                    &
     &                      ncid = ncFLTid(ng),                         &
     &                      varid = fltVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out floats X-grid locations.
!
      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(ixgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                      &
     &                      'Xgrid', Tout,                              &
     &                      (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),      &
     &                      ncid = ncFLTid(ng),                         &
     &                      varid = fltVid(idXgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out floats Y-grid locations.
!
      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(iygrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                      &
     &                      'Ygrid', Tout,                              &
     &                      (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),      &
     &                      ncid = ncFLTid(ng),                         &
     &                      varid = fltVid(idYgrd,ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out floats Z-grid locations.
!
      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(izgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                      &
     &                      'Zgrid', Tout,                              &
     &                      (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),      &
     &                      ncid = ncFLTid(ng),                         &
     &                      varid = fltVid(idZgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out floats (lon,lat) or (x,y) locations.
!
      DO l=1,Nfloats(ng)
        Tout(l)=FLT(ng)%track(iflon,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'lon', Tout,                              &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idglon,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'x', Tout,                                &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idglon,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Nfloats(ng)
        Tout(l)=FLT(ng)%track(iflat,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'lat', Tout,                              &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idglat,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'y', Tout,                                &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idglat,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out floats depths.
!
      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(idpth,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                      &
     &                      'depth', Tout,                              &
     &                      (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),      &
     &                      ncid = ncFLTid(ng),                         &
     &                      varid = fltVid(iddpth,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out density anomaly.
!
      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(ifden,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                      &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),      &
     &                      ncid = ncFLTid(ng),                         &
     &                      varid = fltVid(idDano,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Tout(l)=FLT(ng)%track(ifTvar(itrc),nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltTid(itrc,ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize floats NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, FLTname(ng), ncFLTid(ng))

#else
      SUBROUTINE wrt_floats
#endif
      RETURN
      END SUBROUTINE wrt_floats
