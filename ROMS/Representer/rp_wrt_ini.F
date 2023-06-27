#include "cppdefs.h"
      MODULE rp_wrt_ini_mod
#if (defined TANGENT || defined TL_IOMS) && defined FOUR_DVAR
!
!git $Id$
!svn $Id: rp_wrt_ini.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes initial conditions into representer model file   !
!  using either the standard NetCDF library or the Parallel-IO (PIO)   !
!  library.                                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Tindex     State variables time index to write.                  !
!     OutRec     NetCDF file unlimited dimension record to write.      !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
# endif
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT_NOT_YET || defined BBL_MODEL_NOT_YET
      USE mod_sediment
# endif
      USE mod_stepping
!
      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
# endif
# ifdef SOLVE3D
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
#  ifdef ADJUST_BOUNDARY
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
#  endif
# endif
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: rp_wrt_ini
      PRIVATE :: rp_wrt_ini_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: rp_wrt_ini_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rp_wrt_ini (ng, tile, Tindex, OutRec)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
!
!  Local variable declarations.
!
# ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
# endif
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
# ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
# endif
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (IRP(ng)%IOtype)
        CASE (io_nf90)
          CALL rp_wrt_ini_nf90 (ng, tile, Tindex, OutRec,               &
# ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
# endif
     &                          LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL rp_wrt_ini_pio (ng, tile, Tindex, OutRec,                &
#  ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
#  endif
     &                         LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IRP(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' RP_WRT_INI - Illegal output file type, io_type = ',i0,  &
     &        /,14x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE rp_wrt_ini
!
!***********************************************************************
      SUBROUTINE rp_wrt_ini_nf90 (ng, tile, Tindex, OutRec,             &
# ifdef ADJUST_BOUNDARY
     &                            LBij, UBij,                           &
# endif
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: gfactor, gtype, i, itrc, status, varid
!
      real(dp) :: my_time, scale
!
      character (len=15) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", rp_wrt_ini_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out tangent linear initial conditions.
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
!  Write out model time (s).  Use the "tdays" variable here because of
!  the management of the "time" variable due to nesting.
!
      my_time=tdays(ng)*day2sec

      CALL netcdf_put_fvar (ng, iRPM, IRP(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), my_time,             &
     &                      (/OutRec/), (/1/),                          &
     &                      ncid = IRP(ng)%ncid,                        &
     &                      varid = IRP(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m)
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idFsur),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef WET_DRY
     &                   OCEAN(ng) % tl_zeta(:,:,Tindex),               &
     &                   SetFillVal = .FALSE.)
# else
     &                   OCEAN(ng) % tl_zeta(:,:,Tindex))
# endif
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iRPM, IRP(ng)%name, IRP(ng)%ncid,    &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         IRP(ng)%Vid(idSbry(isFsur)),             &
     &                         OutRec, r2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D U-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idUbar),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % tl_ubar(:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iRPM, IRP(ng)%name, IRP(ng)%ncid,    &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         IRP(ng)%Vid(idSbry(isUbar)),             &
     &                         OutRec, u2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D V-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idVbar),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % tl_vbar(:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iRPM, IRP(ng)%name, IRP(ng)%ncid,    &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         IRP(ng)%Vid(idSbry(isVbar)),             &
     &                         OutRec, v2dvar,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!!    scale=rho0
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idUsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % tl_ustr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
!!    scale=rho0
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idVsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % tl_vstr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idUvel),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_u(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iRPM, IRP(ng)%name, IRP(ng)%ncid,    &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         IRP(ng)%Vid(idSbry(isUvel)),             &
     &                         OutRec, u3dvar,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D V-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Vid(idVvel),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_v(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iRPM, IRP(ng)%name, IRP(ng)%ncid,    &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         IRP(ng)%Vid(idSbry(isVvel)),             &
     &                         OutRec, v3dvar,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iRPM, IRP(ng)%ncid, IRP(ng)%Tid(itrc),   &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_t(:,:,:,Tindex,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          status=nf_fwrite3d_bry(ng, iRPM, IRP(ng)%name, IRP(ng)%ncid,  &
     &                           Vname(1,idSbry(isTvar(itrc))),         &
     &                           IRP(ng)%Vid(idSbry(isTvar(itrc))),     &
     &                           OutRec, r3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,   &
     &                                                   Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iRPM, IRP(ng)%ncid,                    &
     &                       IRP(ng)%Vid(idTsur(itrc)),                 &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tl_tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent linear initial NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iRPM, IRP(ng)%name, IRP(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) THEN
        IF (OutRec.eq.2) THEN
          string='initial  fields'
        END IF
# ifdef SOLVE3D
        WRITE (stdout,20) string, outer, inner, Tindex, Tindex, OutRec
# else
        WRITE (stdout,20) string, outer, inner, Tindex, OutRec
# endif
      END IF
!
  10  FORMAT (/,' RP_WRT_INI_NF90 - error while writing variable: ',a,  &
     &        /,19x,'into representers initial file for time record:',  &
     &        1x,i4)
  20  FORMAT (2x,'RP_WRT_INI_NF90  - wrote ',a,                         &
     &        ' (Outer=',i2.2,', Inner=',i3.3,', Index=',i0,            &
# ifdef SOLVE3D
     &        ',',i0,', Rec=',i0,')')
# else
     &        ', Rec=',i0,')')
# endif
!
      RETURN
      END SUBROUTINE rp_wrt_ini_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE rp_wrt_ini_pio (ng, tile, Tindex, OutRec,              &
#  ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
#  endif
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Tindex, OutRec
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: i, ifield, itrc, status, varid
!
      real(dp) :: my_time, scale
!
      character (len=15) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", rp_wrt_ini_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out tangent linear initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out model time (s).  Use the "tdays" variable here because of
!  the management of the "time" variable due to nesting.
!
      my_time=tdays(ng)*day2sec

      CALL pio_netcdf_put_fvar (ng, iRPM, IRP(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), my_time,         &
     &                          (/OutRec/), (/1/),                      &
     &                          pioFile = IRP(ng)%pioFile,              &
     &                          pioVar = IRP(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m)
!
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idFsur),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
#  ifdef WET_DRY
     &                   OCEAN(ng) % tl_zeta(:,:,Tindex),               &
     &                   SetFillVal = .FALSE.)
#  else
     &                   OCEAN(ng) % tl_zeta(:,:,Tindex))
#  endif
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        IF (IRP(ng)%pioVar(idSBry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iRPM, IRP(ng)%name,                  &
     &                         IRP(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         IRP(ng)%pioVar(idSbry(isFsur)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D U-momentum component (m/s).
!
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idUbar),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_ubar(:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        IF (IRP(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iRPM, IRP(ng)%name,                  &
     &                         IRP(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         IRP(ng)%pioVar(idSbry(isUbar)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D V-momentum component (m/s).
!
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idVbar),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_vbar(:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        IF (IRP(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iRPM, IRP(ng)%name,                  &
     &                         IRP(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         IRP(ng)%pioVar(idSbry(isVbar)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,    &
     &                                                    Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!!    scale=rho0
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idUsms),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % tl_ustr(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
!!    scale=rho0
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idVsms),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % tl_vstr(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idUvel),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_u(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        IF (IRP(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iRPM, IRP(ng)%name,                  &
     &                         IRP(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         IRP(ng)%pioVar(idSbry(isUvel)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out 3D V-momentum component (m/s).
!
      scale=1.0_dp
      IF (IRP(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iRPM, IRP(ng)%pioFile,                     &
     &                   IRP(ng)%pioVar(idVvel),                        &
     &                   OutRec, ioDesc,                                &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_v(:,:,:,Tindex))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        IF (IRP(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iRPM, IRP(ng)%name,                  &
     &                         IRP(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         IRP(ng)%pioVar(idSbry(isVvel)),          &
     &                         OutRec, ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,     &
     &                                                 Tindex))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (IRP(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iRPM, IRP(ng)%pioFile,                   &
     &                     IRP(ng)%pioTrc(itrc),                        &
     &                     OutRec, ioDesc,                              &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % tl_t(:,:,:,Tindex,itrc))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#   ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          IF (IRP(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry(ng, iRPM, IRP(ng)%name,                &
     &                           IRP(ng)%pioFile,                       &
     &                           Vname(1,ifield),                       &
     &                           IRP(ng)%pioVar(ifield),                &
     &                           OutRec, ioDesc,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,   &
     &                                                   Tindex,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif

#   ifdef ADJUST_STFLUX
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp
          IF (IRP(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iRPM, IRP(ng)%pioFile,                 &
     &                       IRP(ng)%pioVar(idTsur(itrc)),              &
     &                       OutRec, ioDesc,                            &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % tl_tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent linear initial NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iRPM, IRP(ng)%name, IRP(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) THEN
        IF (OutRec.eq.2) THEN
          string='initial  fields'
        END IF
#  ifdef SOLVE3D
        WRITE (stdout,20) string, outer, inner, Tindex, Tindex, OutRec
#  else
        WRITE (stdout,20) string, outer, inner, Tindex, OutRec
#  endif
      END IF
!
  10  FORMAT (/,' RP_WRT_INI_PIO - error while writing variable: ',a,   &
     &        /,18x,'into representers initial file for time record:',  &
     &        1x,i4)
  20  FORMAT (2x,'RP_WRT_INI_PIO   - wrote ',a,                         &
     &        ' (Outer=',i2.2,', Inner=',i3.3,', Index=',i0,            &
#  ifdef SOLVE3D
     &        ',',i0,', Rec=',i0,')')
#  else
     &        ', Rec=',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE rp_wrt_ini_pio
# endif
#endif
      END MODULE rp_wrt_ini_mod
