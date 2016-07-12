#include "cppdefs.h"
#if defined FOUR_DVAR || defined ENKF_RESTART
      SUBROUTINE wrt_dai (ng, tile)
!
!svn $Id: wrt_dai.F 798 2016-05-21 17:04:58Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out Data Assimilation initial conditions        !
!  (4D-Var analysis) or Ensemble Kalman Filter (EnKF) restart file.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: i, j, k, itrc
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, status, varid
      real(r8) :: scale

# include "set_bounds.h"
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile='wrt_dai.F'
!
!-----------------------------------------------------------------------
!  Write out Data Assimilation initial/restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields.
!
      gfactor=1
!
!  Set time record index.
!
      DAI(ng)%Rindex=DAI(ng)%Rindex+1
      Fcount=DAI(ng)%Fcount
      DAI(ng)%Nrec(Fcount)=DAI(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in restart
!  file.
!
      DAI(ng)%Rindex=MOD(DAI(ng)%Rindex-1,2)+1

# ifdef SOLVE3D
!
!  Write out time independent (unperturb) depths of RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idpthR),   &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   GRID(ng) % z0_r,                               &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idpthR))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out time independent (unperturb) depths of U-points.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
            GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z0_r(i-1,j,k)+         &
     &                                  GRID(ng)%z0_r(i  ,j,k))
          END DO
        END DO
      END DO
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idpthU),   &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   GRID(ng) % z_v,                                &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idpthU))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out time independent (unperturb) depths of V-points.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
            GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z0_r(i,j-1,k)+         &
     &                                  GRID(ng)%z0_r(i,j  ,k))
          END DO
        END DO
      END DO
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idpthV),   &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % z_v,                                &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idpthV))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out time independent (unperturb) depths of W-points.
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idpthW),   &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   GRID(ng) % z0_w,                               &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idpthW))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, DAI(ng)%name,                     &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/DAI(ng)%Rindex/), (/1/),                  &
     &                      ncid = DAI(ng)%ncid,                        &
     &                      varid = DAI(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idFsur),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef WET_DRY
     &                   OCEAN(ng) % zeta(:,:,KOUT),                    &
     &                   SetFillVal = .FALSE.)
# else
     &                   OCEAN(ng) % zeta(:,:,KOUT))
# endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idUbar),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % ubar(:,:,KOUT))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idVbar),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % vbar(:,:,KOUT))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idUvel),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % u(:,:,:,NOUT))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idVvel),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % v(:,:,:,NOUT))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Tid(itrc),   &
     &                     DAI(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % t(:,:,:,NOUT,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))), DAI(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, DAI(ng)%name, DAI(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,30) KOUT, NOUT, DAI(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,30) KOUT, NOUT, DAI(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,30) KOUT, DAI(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,30) KOUT, DAI(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' WRT_DAI - error while writing variable: ',a,/,11x,    &
     &        'into DA initial/restart NetCDF file.')
  20  FORMAT (/,' WRT_DAI - error while writing variable: ',a,/,11x,    &
     &        'into DA initial/rstart NetCDF file for time record: ',i4)
# ifdef SOLVE3D
  30  FORMAT (6x,'WRT_DAI   - wrote DA IC    fields (Index=', i1,       &
#  ifdef NESTING
     &        ',',i1,') into time record = ',i7.7,2x,i2.2)
#  else
     &        ',',i1,') into time record = ',i7.7)
#  endif
# else
  30  FORMAT (6x,'WRT_DAI   - wrote DA IC    fields (Index=', i1,       &
#  ifdef NESTING
     &        ') into time record = ',i7.7,2x,i2.2)
#  else
     &        ') into time record = ',i7.7)
#  endif
# endif

      RETURN
      END SUBROUTINE wrt_dai
#else
      SUBROUTINE wrt_dai
      RETURN
      END SUBROUTINE wrt_dai
#endif
