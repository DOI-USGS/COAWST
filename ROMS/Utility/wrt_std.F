#include "cppdefs.h"
      MODULE wrt_std_mod

#if defined STD_MODEL
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes standard deviation output file, which is to      !
!  model the 4D-Var Bacground Error Covariance, B. It defines its      !
!  dimensions, attributes, and variables using either the standard     !
!  NetCDF library or the Parallel-IO (PIO) library.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
# endif
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_std
      PRIVATE :: wrt_std_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_std_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_std (ng, tile, Lstd)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstd
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out Background standard deviation fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (STD(5,ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_std_nf90 (ng, tile, Lstd,                            &
     &                       LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_std_pio (ng, tile, Lstd,                             &
     &                      LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) STD(5,ng)%IOtype
          exit_flag=3
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
      LwrtSTD(ng)=.FALSE.
!
  10  FORMAT (' WRT_STD - Illegal output file type, io_type = ',i0,     &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_std
!
!***********************************************************************
      SUBROUTINE wrt_std_nf90 (ng, tile, Lstd,                          &
     &                         LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstd
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, i, j, gfactor, gtype, status
# ifdef SOLVE3D
      integer :: itrc, k
# endif
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_std_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out full posterior error covariance (diagonal) matrix.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!  Set time record index.
!
      STD(5,ng)%Rindex=STD(5,ng)%Rindex+1
      Fcount=STD(5,ng)%Fcount
      STD(5,ng)%Nrec(Fcount)=STD(5,ng)%Nrec(Fcount)+1
!
!  Report.
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) Lstd, Lstd, STD(5,ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) Lstd, Lstd, STD(5,ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) Lstd, STD(5,ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) Lstd, STD(5,ng)%Rindex
#  endif
# endif
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, STD(5,ng)%name,                   &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/STD(5,ng)%Rindex/), (/1/),                &
     &                      ncid = STD(5,ng)%ncid,                      &
     &                      varid = STD(5,ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface standard deviation.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, STD(5,ng)%ncid, idFsur,              &
     &                   STD(5,ng)%Vid(idFsur),                         &
     &                   STD(5,ng)%Rindex, gtype,                       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng)% e_zeta(:,:,Lstd))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D U-momentum component standard deviation.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, STD(5,ng)%ncid, idUbar,              &
     &                   STD(5,ng)%Vid(idUbar),                         &
     &                   STD(5,ng)%Rindex, gtype,                       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % e_ubar(:,:,Lstd))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D V-momentum component standard deviation.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, STD(5,ng)%ncid, idVbar,              &
     &                   STD(5,ng)%Vid(idVbar),                         &
     &                   STD(5,ng)%Rindex, gtype,                       &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % e_vbar(:,:,Lstd))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component standard deviation.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, STD(5,ng)%ncid, idUvel,              &
     &                   STD(5,ng)%Vid(idUvel),                         &
     &                   STD(5,ng)%Rindex, gtype,                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % e_u(:,:,:,Lstd))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D V-momentum component standard deviation.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, STD(5,ng)%ncid, idVvel,              &
     &                   STD(5,ng)%Vid(idVvel),                         &
     &                   STD(5,ng)%Rindex, gtype,                       &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % e_v(:,:,:,Lstd))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables standard deviation.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, STD(5,ng)%ncid, idTvar(itrc),      &
     &                     STD(5,ng)%Tid(itrc),                         &
     &                     STD(5,ng)%Rindex, gtype,                     &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % e_t(:,:,:,Lstd,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),              &
     &                        STD(5,ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
# if !(defined OUT_NETCDF4 && defined DELAYED_SYNC_NF90)
!
!-----------------------------------------------------------------------
!  Synchronize Background standard deviation NetCDF file to disk to
!  allow other processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, STD(5,ng)%name, STD(5,ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
  10  FORMAT (2x,'WRT_STD_NF90     - writing standard deviation', t42,  &
#  ifdef SOLVE3D
#   ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#   endif
#  else
#   ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,')   in record = ',i0)
#   endif
#  endif
  20  FORMAT (/,' WRT_STD_NF90 - error while writing variable: ',a,     &
     &        /,16x,'into 4DVar error NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_std_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_std_pio (ng, tile, Lstd,                           &
     &                        LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstd
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, i, ifield, j, status
#  ifdef SOLVE3D
      integer :: itrc, k
#  endif
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_std_nf90"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out full posterior error covariance (diagonal) matrix.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      STD(5,ng)%Rindex=STD(5,ng)%Rindex+1
      Fcount=STD(5,ng)%Fcount
      STD(5,ng)%Nrec(Fcount)=STD(5,ng)%Nrec(Fcount)+1
!
!  Report.
!
#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) Lstd, Lstd, STD(5,ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) Lstd, Lstd, STD(5,ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) Lstd, STD(5,ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) Lstd, STD(5,ng)%Rindex
#   endif
#  endif
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, iNLM, STD(5,ng)%name,               &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/STD(5,ng)%Rindex/), (/1/),            &
     &                          pioFile = STD(5,ng)%pioFile,            &
     &                          pioVar = STD(5,ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface standard deviation.
!
      scale=1.0_dp
      IF (STD(5,ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, STD(5,ng)%pioFile, idFsur,           &
     &                   STD(5,ng)%pioVar(idFsur),                      &
     &                   STD(5,ng)%Rindex,                              &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   OCEAN(ng)% e_zeta(:,:,Lstd))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifndef SOLVE3D
!
!  Write out 2D U-momentum component standard deviation.
!
      scale=1.0_dp
      IF (STD(5,ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, STD(5,ng)%pioFile, idUbar,           &
     &                   STD(5,ng)%pioVar(idUbar),                      &
     &                   STD(5,ng)%Rindex,                              &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % e_ubar(:,:,Lstd))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D V-momentum component standard deviation.
!
      scale=1.0_dp
      IF (STD(5,ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iNLM, STD(5,ng)%pioFile, idVbar,           &
     &                   STD(5,ng)%pioVar(idVbar),                      &
     &                   STD(5,ng)%Rindex,                              &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % e_vbar(:,:,Lstd))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component standard deviation.
!
      scale=1.0_dp
      IF (STD(5,ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, STD(5,ng)%pioFile, idUvel,           &
     &                   STD(5,ng)%pioVar(idUvel),                      &
     &                   STD(5,ng)%Rindex,                              &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % e_u(:,:,:,Lstd))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D V-momentum component error variance.
!
      scale=1.0_dp
      IF (STD(5,ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iNLM, STD(5,ng)%pioFile, idVvel,           &
     &                   STD(5,ng)%pioVar(idVvel),                      &
     &                   STD(5,ng)%Rindex,                              &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % e_v(:,:,:,Lstd))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), STD(5,ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables standard deviation.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (STD(5,ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iNLM, STD(5,ng)%pioFile, idTvar(itrc),   &
     &                     STD(5,ng)%pioTrc(itrc),                      &
     &                     STD(5,ng)%Rindex,                            &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % e_t(:,:,:,Lstd,itrc))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),              &
     &                        STD(5,ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
#  ifndef DELAYED_SYNC_PIO
!
!-----------------------------------------------------------------------
!  Synchronize Background standard deviation NetCDF file to disk to
!  allow other processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iNLM, STD(5,ng)%name,                   &
     &                      STD(5,ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
  10  FORMAT (2x,'WRT_STD_PIO      - writing standard deviation', t42,  &
#  ifdef SOLVE3D
#   ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#   endif
#  else
#   ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,')   in record = ',i0)
#   endif
#  endif
  20  FORMAT (/,' WRT_STD_PIO - error while writing variable: ',a,      &
     &        /,15x,'into 4DVar error NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_std_pio
# endif
#endif
      END MODULE wrt_std_mod

