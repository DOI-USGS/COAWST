#include "cppdefs.h"
      MODULE wrt_aug_imp_mod
#if defined WEAK_CONSTRAINT && defined RPCG
!
!git $Id$
!svn $Id: wrt_aug_imp.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes augmented impulse state variables using the     !
!  standard NetCDF library or the Parallel-IO (PIO) library.           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      PUBLIC  :: wrt_aug_imp
      PRIVATE :: wrt_aug_imp_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_aug_imp_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_aug_imp (ng, tile, model, Iinp, Iout)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Iinp, Iout
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
      SELECT CASE (TLF(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_aug_imp_nf90 (ng, tile, model, Iinp, Iout,           &
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_aug_imp_pio (ng, tile, model, Iinp, Iout,            &
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TLF(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_AUG_IMP - Illegal output tile type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.' )
!
      RETURN
      END SUBROUTINE wrt_aug_imp
!
!***********************************************************************
      SUBROUTINE wrt_aug_imp_nf90 (ng, tile, model, Iinp, Iout,         &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Iinp, Iout
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: itrc, gfactor, gtype, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_aug_imp_nf90"

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
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
!  Process free-surface weak-constraint impulse forcing.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                       &
     &                   TLF(ng)%Vid(idZtlf),                           &
     &                   Iout, gtype,                                   &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng) % tl_zeta(:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idZtlf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifndef SOLVE3D
!
!  Process 2D U-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                       &
     &                   TLF(ng)%Vid(idUbtf),                           &
     &                   Iout, gype,                                    &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_ubar(:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbtf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Process 2D V-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                       &
     &                   TLF(ng)%Vid(idVbtf),                           &
     &                   Iout, gtype,                                   &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_vbar(:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbtf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
!
!  Process 3D U-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, model, TLF(ng)%ncid,                       &
     &                   TLF(ng)%Vid(idUtlf),                           &
     &                   Iout, gtype,                                   &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_u(:,:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUtlf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Process 3D V-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, model, TLF(ng)%ncid,                       &
     &                   TLF(ng)%Vid(idVtlf),                           &
     &                   Iout, gtype,                                   &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_v(:,:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVtlf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Process tracer type variables impulses.
!
      scale=1.0_dp
      gtype=gfactor*r3dvar
      DO itrc=1,NT(ng)
        status=nf_fwrite3d(ng, model, TLF(ng)%ncid,                     &
     &                     TLF(ng)%Tid(itrc),                           &
     &                     Iout, gtype,                                 &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % tl_t(:,:,:,Iinp,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTtlf(itrc))), Iout,        &
     &                        TRIM(TLF(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, TLF(ng)%name, TLF(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master)  WRITE (stdout,20) Iout, TRIM(TLF(ng)%name)
!
  10  FORMAT (/,' WRT_AUG_IMP_NF90 - error while writing variable: ',a, &
     &        2x,'at time record = ',i0,                                &
     &        /,20x,'into NetCDF file: ',a)
  20  FORMAT (2x,'WRT_AUG_IMP_NF90 - wrote augmented adjoint',          &
     &        ' impulses, record: ',i0,/,22x,'file: ',a)
!
      RETURN
      END SUBROUTINE wrt_aug_imp_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_aug_imp_pio (ng, tile, model, Iinp, Iout,          &
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, Iinp, Iout
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: itrc, status
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__", wrt_aug_imp_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc

#  include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Process free-surface weak-constraint impulse forcing.
!
      scale=1.0_dp
      IF (TLF(ng)%pioVar(idZtlf)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                    &
     &                   TLF(ng)%pioVar(idZtlf), Iout,                  &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   OCEAN(ng) % tl_zeta(:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idZtlf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifndef SOLVE3D
!
!  Process 2D U-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      IF (TLF(ng)%pioVar(idUbtf)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                    &
     &                   TLF(ng)%pioVar(idUbtf), Iout,                  &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_ubar(:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbtf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Process 2D V-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      IF (TLF(ng)%pioVar(idVbtf)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                    &
     &                   TLF(ng)%pioVar(idVbtf), Iout,                  &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_vbar(:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbtf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef SOLVE3D
!
!  Process 3D U-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      IF (TLF(ng)%pioVar(idUtlf)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                    &
     &                   TLF(ng)%pioVar(idUtlf), Iout,                  &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_u(:,:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUtlf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Process 3D V-momentum weak-constraint impulse forcing.
!
      scale=1.0_dp
      IF (TLF(ng)%pioVar(idVtlf)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                    &
     &                   TLF(ng)%pioVar(idVtlf), Iout,                  &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_v(:,:,:,Iinp))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVtlf)), Iout,                &
     &                      TRIM(TLF(ng)%name)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Process tracer type variables impulses.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (TLF(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                  &
     &                     TLF(ng)%pioTrc(itrc), Iout,                  &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_t(:,:,:,Iinp,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTtlf(itrc))), Iout,        &
     &                        TRIM(TLF(ng)%name)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, TLF(ng)%name, TLF(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master)  WRITE (stdout,20) Iout, TRIM(TLF(ng)%name)
!
  10  FORMAT (/,' WRT_AUG_IMP_PIO - error while writing variable: ',a,  &
     &        2x,'at time record = ',i0,                                &
     &        /,20x,'into NetCDF file: ',a)
  20  FORMAT (2x,'WRT_AUG_IMP_PIO  - wrote augmented adjoint',          &
     &        ' impulses, record: ',i0,/,22x,'file: ',a)
!
      RETURN
      END SUBROUTINE wrt_aug_imp_pio
# endif
#endif
      END MODULE wrt_aug_imp_mod
