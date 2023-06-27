#include "cppdefs.h"
      MODULE wrt_dai_mod
#if defined FOUR_DVAR || defined ENKF_RESTART
!
!git $Id$
!svn $Id: wrt_dai.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out Data Assimilation initial conditions        !
!  (4D-Var analysis) or Ensemble Kalman Filter (EnKF) restart file     !
!  using the standard NetCDF library or the Parallel-IO (PIO) library. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_dai
      PRIVATE :: wrt_dai_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_dai_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_dai (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (DIA(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_dai_nf90 (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_dai_pio (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) DAI(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_DAI - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_dai
!
!***********************************************************************
      SUBROUTINE wrt_dai_nf90 (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, j, k, itrc
      integer :: Fcount, gfactor, gtype, status, varid
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_dai_nf90"

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Data Assimilation initial/restart fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
      scale=1.0_dp
      gtype=gfactor*r3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idpthR),   &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   GRID(ng) % z0_r,                               &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
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
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
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
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idpthW),   &
     &                   0, gtype,                                      &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   GRID(ng) % z0_w,                               &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/DAI(ng)%Rindex/), (/1/),                  &
     &                      ncid = DAI(ng)%ncid,                        &
     &                      varid = DAI(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idFsur),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# if defined R4DVAR || defined SPLIT_R4DVAR
#  ifdef WET_DRY
     &                   OCEAN(ng) % tl_zeta(:,:,KOUT),                 &
     &                   SetFillVal = .FALSE.)
#  else
     &                   OCEAN(ng) % tl_zeta(:,:,KOUT))
#  endif
# else
#  ifdef WET_DRY
     &                   OCEAN(ng) % zeta(:,:,KOUT),                    &
     &                   SetFillVal = .FALSE.)
#  else
     &                   OCEAN(ng) % zeta(:,:,KOUT))
#  endif
# endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idUbar),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
# if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_ubar(:,:,KOUT))
# else
     &                   OCEAN(ng) % ubar(:,:,KOUT))
# endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idVbar),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
# if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_vbar(:,:,KOUT))
# else
     &                   OCEAN(ng) % vbar(:,:,KOUT))
# endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idUvel),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
#  if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_u(:,:,:,NOUT))
#  else
     &                   OCEAN(ng) % u(:,:,:,NOUT))
#  endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idVvel),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
#  if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_v(:,:,:,NOUT))
#  else
     &                   OCEAN(ng) % v(:,:,:,NOUT))
#  endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Tid(itrc),   &
     &                     DAI(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
#  if defined R4DVAR || defined SPLIT_R4DVAR
     &                     OCEAN(ng) % tl_t(:,:,:,NOUT,itrc))
#  else
     &                     OCEAN(ng) % t(:,:,:,NOUT,itrc))
#  endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))), DAI(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idVvis),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvis)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idTdif),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idTdif)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%ncid, DAI(ng)%Vid(idSdif),   &
     &                   DAI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idSdif)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, DAI(ng)%name, DAI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

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
  10  FORMAT (/,' WRT_DAI_NF90 - error while writing variable: ',a,     &
     &        /,11x,'into DA initial/restart NetCDF file.')
  20  FORMAT (/,' WRT_DAI_NF90 - error while writing variable: ',a,     &
     &        /,11x,'into DA initial/rstart NetCDF file for time ',     &
     &         'record: ',i0)
# ifdef SOLVE3D
  30  FORMAT (2x,'WRT_DAI_NF90     - wrote DA INI/RST', t40,            &
#  ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#  endif
# else
  30  FORMAT (2x,'WRT_DAI_NF90     - wrote DA INI/RST', t40,            &
#  ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,')   in record = ',i0)
#  endif
# endif
!
      RETURN
      END SUBROUTINE wrt_dai_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_dai_pio (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, j, k, itrc
      integer :: Fcount, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_dai_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc

#  include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Data Assimilation initial/restart fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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

#  ifdef SOLVE3D
!
!  Write out time independent (unperturb) depths of RHO-points.
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idpthR)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idpthR), 0,                     &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   GRID(ng) % z0_r,                               &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=IstrU-1,Iend+1
            GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z0_r(i-1,j,k)+         &
     &                                  GRID(ng)%z0_r(i  ,j,k))
          END DO
        END DO
      END DO
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idpthU)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idpthU), 0,                     &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   GRID(ng) % z_v,                                &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
      DO k=1,N(ng)
        DO j=JstrV-1,Jend+1
          DO i=Istr-1,Iend+1
            GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z0_r(i,j-1,k)+         &
     &                                  GRID(ng)%z0_r(i,j  ,k))
          END DO
        END DO
      END DO
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idpthV)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idpthV), 0,                     &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   GRID(ng) % z_v,                                &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idpthW)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idpthW), 0,                     &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   GRID(ng) % z0_w,                               &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idpthW))
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, iNLM, DAI(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/DAI(ng)%Rindex/), (/1/),              &
     &                          pioFile = DAI(ng)%pioFile,              &
     &                          pioVar = DAI(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idFsur),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
#  if defined R4DVAR || defined SPLIT_R4DVAR
#   ifdef WET_DRY
     &                   OCEAN(ng) % tl_zeta(:,:,KOUT),                 &
     &                   SetFillVal = .FALSE.)
#   else
     &                   OCEAN(ng) % tl_zeta(:,:,KOUT))
#   endif
#  else
#   ifdef WET_DRY
     &                   OCEAN(ng) % zeta(:,:,KOUT),                    &
     &                   SetFillVal = .FALSE.)
#   else
     &                   OCEAN(ng) % zeta(:,:,KOUT))
#   endif
#  endif
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idUbar),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
#  if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_ubar(:,:,KOUT))
#  else
     &                   OCEAN(ng) % ubar(:,:,KOUT))
#  endif
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idVbar),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
#  if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_vbar(:,:,KOUT))
#  else
     &                   OCEAN(ng) % vbar(:,:,KOUT))
#  endif
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idUvel),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
#   if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_u(:,:,:,NOUT))
#   else
     &                   OCEAN(ng) % u(:,:,:,NOUT))
#   endif
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idVvel),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
#   if defined R4DVAR || defined SPLIT_R4DVAR
     &                   OCEAN(ng) % tl_v(:,:,:,NOUT))
#   else
     &                   OCEAN(ng) % v(:,:,:,NOUT))
#   endif
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        IF (DAI(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                   &
     &                     DAI(ng)%pioTrc(itrc),                        &
     &                     DAI(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
#   if defined R4DVAR || defined SPLIT_R4DVAR
     &                     OCEAN(ng) % tl_t(:,:,:,NOUT,itrc))
#   else
     &                     OCEAN(ng) % t(:,:,:,NOUT,itrc))
#   endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))), DAI(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idVvis),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvis)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idTdif),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idTdif)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_dp
      IF (DAI(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_w3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_w3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, DAI(ng)%pioFile,                     &
     &                   DAI(ng)%pioVar(idSdif),                        &
     &                   DAI(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#    ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#    endif
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idSdif)), DAI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize restart NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, DAI(ng)%name, DAI(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,30) KOUT, NOUT, DAI(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,30) KOUT, NOUT, DAI(ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,30) KOUT, DAI(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,30) KOUT, DAI(ng)%Rindex
#   endif
#  endif
!
  10  FORMAT (/,' WRT_DAI_PIO - error while writing variable: ',a,      &
     &        /,11x,'into DA initial/restart NetCDF file.')
  20  FORMAT (/,' WRT_DAI_PIO - error while writing variable: ',a,      &
     &        /,11x,'into DA initial/rstart NetCDF file for time ',     &
     &         'record: ',i0)
#  ifdef SOLVE3D
  30  FORMAT (2x,'WRT_DAI_PIO      - wrote DA INI/RST', t40,            &
#   ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#   endif
#  else
  30  FORMAT (2x,'WRT_DAI_PIO      - wrote DA INI/RST', t40,            &
#   ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,')   in record = ',i0)
#   endif
#  endif
!
      RETURN
      END SUBROUTINE wrt_dai_pio
# endif
#endif
      END MODULE wrt_dai_mod
