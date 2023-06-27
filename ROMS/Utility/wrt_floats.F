#include "cppdefs.h"
      MODULE wrt_floats_mod
#ifdef FLOATS
!
!git $Id$
!svn $Id: wrt_floats.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes simulated drifter trajectories into output       !
!  file using the standard NetCDF library or the Parallel-IO (PIO)     !
!  library.                                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_stepping
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_floats
      PRIVATE :: wrt_floats_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_floats_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_floats (ng)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (FLT(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_floats_nf90 (ng)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_floats_pio (ng)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) FLT(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_FLOATS - Illegal output file type, io_type = ',i0,  &
     &        /,14x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_floats
!
!***********************************************************************
      SUBROUTINE wrt_floats_nf90 (ng)
!***********************************************************************
!
      USE mod_netcdf
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
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_floats_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Nfloats(ng)
          IF (DRIFTER(ng)%bounded(l)) THEN
            Tout(l)=DRIFTER(ng)%track(ifTvar(itrc),nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                   &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,FLT(ng)%Rindex/), (/Nfloats(ng),1/),  &
     &                        ncid = FLT(ng)%ncid,                      &
     &                        varid = FLT(ng)%Tid(itrc))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Synchronize floats NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, FLT(ng)%name, FLT(ng)%ncid)
!
      RETURN
      END SUBROUTINE wrt_floats_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_floats_pio (ng)
!***********************************************************************
!
      USE mod_pio_netcdf
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
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_floats_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      FLT(ng)%Rindex=FLT(ng)%Rindex+1
      Fcount=FLT(ng)%Fcount
      FLT(ng)%Nrec(Fcount)=FLT(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/FLT(ng)%Rindex/), (/1/),              &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'Xgrid', Tout,                          &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                           pioFile = FLT(ng)%pioFile,             &
     &                           pioVar = FLT(ng)%pioVar(idXgrd)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'Ygrid', Tout,                          &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idYgrd)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'Zgrid', Tout,                          &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idZgrd)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
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
        CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,               &
     &                            'lon', Tout,                          &
     &                            (/1,FLT(ng)%Rindex/),                 &
     &                            (/Nfloats(ng),1/),                    &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            pioVar = FLT(ng)%pioVar(idglon)%vd)
      ELSE
        CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,               &
     &                            'x', Tout,                            &
     &                            (/1,FLT(ng)%Rindex/),                 &
     &                            (/Nfloats(ng),1/),                    &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            pioVar = FLT(ng)%pioVar(idglon)%vd)
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      DO l=1,Nfloats(ng)
        IF (DRIFTER(ng)%bounded(l)) THEN
          Tout(l)=DRIFTER(ng)%track(iflat,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      IF (spherical) THEN
        CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,               &
     &                            'lat', Tout,                          &
     &                            (/1,FLT(ng)%Rindex/),                 &
     &                            (/Nfloats(ng),1/),                    &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            pioVar = FLT(ng)%pioVar(idglat)%vd)
      ELSE
        CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,               &
     &                            'y', Tout,                            &
     &                            (/1,FLT(ng)%Rindex/),                 &
     &                            (/Nfloats(ng),1/),                    &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            pioVar = FLT(ng)%pioVar(idglat)%vd)
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'depth', Tout,                          &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(iddpth)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          TRIM(Vname(1,idDano)), Tout,            &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idDano)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Nfloats(ng)
          IF (DRIFTER(ng)%bounded(l)) THEN
            Tout(l)=DRIFTER(ng)%track(ifTvar(itrc),nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,               &
     &                            TRIM(Vname(1,idTvar(itrc))), Tout,    &
     &                            (/1,FLT(ng)%Rindex/),                 &
     &                            (/Nfloats(ng),1/),                    &
     &                            pioFile = FLT(ng)%pioFile,            &
     &                            pioVar = FLT(ng)%pioTrc(itrc)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#  endif

#  ifdef FLOAT_OYSTER
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'swim_time', Tout,                      &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = Flt(ng)%pioVar(idswim)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'w_bio', Tout,                          &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idwbio))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'bio_size', Tout,                       &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idsize)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      CALL pio_netcdf_put_fvar (ng, iNLM, FLT(ng)%name,                 &
     &                          'bio_sink', Tout,                       &
     &                          (/1,FLT(ng)%Rindex/),                   &
     &                          (/Nfloats(ng),1/),                      &
     &                          pioFile = FLT(ng)%pioFile,              &
     &                          pioVar = FLT(ng)%pioVar(idwsin)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize floats NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, FLT(ng)%name, FLT(ng)%pioFile)
!
      RETURN
      END SUBROUTINE wrt_floats_pio
# endif
#endif
      END MODULE wrt_floats_mod
