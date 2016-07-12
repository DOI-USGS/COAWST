#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE wrt_floats (ng)
!
!svn $Id: wrt_floats.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
      integer :: Fcount, itrc, l, status

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
      FLT(ng)%Rindex=FLT(ng)%Rindex+1
      Fcount=FLT(ng)%Fcount
      FLT(ng)%Nrec(Fcount)=FLT(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/FLT(ng)%Rindex/), (/1/),                  &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out floats X-grid locations.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(ixgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'Xgrid', Tout,                              &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idXgrd))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out floats Y-grid locations.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iygrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'Ygrid', Tout,                              &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idYgrd))
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out floats Z-grid locations.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(izgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'Zgrid', Tout,                              &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idZgrd))
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef OFFLINE_FLOATS_LATLON
!
!  Write out floats (lon,lat) and (x,y) locations.
!
    IF (.not.spherical) THEN
        WRITE(stdout,*) 'OFFLINE FLOATS MODULE NEEDS SPHERICAL GRID!'
        exit_flag=2
        RETURN
    END IF
    IF (spherical) THEN
      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(iflon,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'lon', Tout,                              &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idglon,ng))
      IF (exit_flag.ne.NoError) RETURN

      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(ixspc,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'x', Tout,                                &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idxspc,ng))
      IF (exit_flag.ne.NoError) RETURN

      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(iflat,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'lat', Tout,                              &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idglat,ng))
      IF (exit_flag.ne.NoError) RETURN

      DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(iyspc,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'y', Tout,                                &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idyspc,ng))
      IF (exit_flag.ne.NoError) RETURN

     DO l=1,Nfloats(ng)
        IF (FLT(ng)%bounded(l)) THEN
          Tout(l)=FLT(ng)%track(iwdph,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
        CALL netcdf_put_fvar (ng, iNLM, FLTname(ng),                    &
     &                        'wd', Tout,                               &
     &                        (/1,tFLTindx(ng)/), (/Nfloats(ng),1/),    &
     &                        ncid = ncFLTid(ng),                       &
     &                        varid = fltVid(idwdph,ng))
      IF (exit_flag.ne.NoError) RETURN
     END IF
# else
!
!  Write out floats (lon,lat) or (x,y) locations.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iflon,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                   &
     &                        'lon', Tout,                              &
     &                        (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),  &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        varid = FLT(ng)%Vid(idglon))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                   &
     &                        'x', Tout,                                &
     &                        (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),  &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        varid = FLT(ng)%Vid(idglon))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iflat,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                   &
     &                        'lat', Tout,                              &
     &                        (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),  &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        varid = FLT(ng)%Vid(idglat))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                   &
     &                        'y', Tout,                                &
     &                        (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),  &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        varid = FLT(ng)%Vid(idglat))
      END IF
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
!
!  Write out floats depths.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(idpth,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'depth', Tout,                              &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(iddpth))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out density anomaly.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(ifden,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idDano))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Nfloats(ng)
          IF (DRIFTER(ng)%bounded(l)) THEN
#  ifdef OFFLINE_FLOATS
            Tout(l)=DRIFTER(ng)%track(itrc+(NFV(ng)-NT(ng)),nf(ng),l)
#  else
            Tout(l)=DRIFTER(ng)%track(ifTvar(itrc),nf(ng),l)
#  endif
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                   &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),  &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        varid = FLT(ng)%Tid(itrc))
        IF (exit_flag.ne.NoError) RETURN
      END DO

#  ifdef OFFLINE_FLOATS
      IF (Master) print *, '  WRT_FLOATS, model time (days), t=', tdays
#  endif
# endif

# ifdef FLOAT_OYSTER
!
!  Write out biological float swimming time.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iswim,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'swim_time', Tout,                          &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = Flt(ng)%Vid(idswim))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out biological float vertical velocity.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iwbio,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'w_bio', Tout,                              &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idwbio))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out biological float size (length).
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
                Tout(l)=DRIFTER(ng)%track(isizf,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'bio_size', Tout,                           &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idsize))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out biological float sinking velocity.
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iwsin,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                     &
     &                      'bio_sink', Tout,                           &
     &                      (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),    &
     &                      ncid = FLT(ng)%ncid,                        &
     &                      varid = FLT(ng)%Vid(idwsin))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Synchronize floats NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, FLT(ng)%name, FLT(ng)%ncid)

#else
      SUBROUTINE wrt_floats
#endif
      RETURN
      END SUBROUTINE wrt_floats
