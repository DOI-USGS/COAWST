#include "cppdefs.h"
      MODULE wrt_evolved_mod
#if (defined I4DVAR && defined EVOLVED_LCZ) || defined LCZ_FINAL
!
!git $Id$
!svn $Id: wrt_evolved.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes the evolved Lanczos vectors into the LZE file   !
!  using the standard NetCDF library or the Parallel-IO (PIO) library. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
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
      PUBLIC  :: wrt_evolved
      PRIVATE :: wrt_evolved_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_evolved_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_evolved (ng, tile, kout, nout)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout, nout
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
!  Write out time-averaged fields according to IO type.
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
      SELECT CASE (LZE(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_evolved_nf90 (ng, tile, kout, nout,                  &
# ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
# endif
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_evolved_pio (ng, tile, kout, nout,                   &
#  ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
#  endif
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) LZE(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_EVOLVED - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_evolved
!
!***********************************************************************
      SUBROUTINE wrt_evolved_nf90 (ng, tile, kout, nout,                &
# ifdef ADJUST_BOUNDARY
     &                             LBij, UBij,                          &
# endif
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout, nout
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
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
     &  __FILE__//", wrt_evolved_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out evolved Lanczos vectors.
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
      LZE(ng)%Rindex=LZE(ng)%Rindex+1
      Fcount=LZE(ng)%Fcount
      LZE(ng)%Nrec(Fcount)=LZE(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iADM, LZE(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/LZE(ng)%Rindex/), (/1/),                  &
     &                      ncid = LZE(ng)%ncid,                        &
     &                      varid = LZE(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idFsur),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng)% tl_zeta(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), LZE(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, LZE(ng)%name, LZE(ng)%ncid,   &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          LZE(ng)%Vid(idSbry(isFsur)),            &
     &                          LZE(ng)%Rindex, r2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D U-momentum component.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idUbar),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % tl_ubar(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), LZE(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, LZE(ng)%name, LZE(ng)%ncid,   &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          LZE(ng)%Vid(idSbry(isUbar)),            &
     &                          LZE(ng)%Rindex, u2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D V-momentum component.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idVbar),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % tl_vbar(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), LZE(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, LZE(ng)%name, LZE(ng)%ncid,   &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          LZE(ng)%Vid(idSbry(isVbar)),            &
     &                          LZE(ng)%Rindex, v2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idUvel),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_u(:,:,:,nout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), LZE(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, iADM, LZE(ng)%name, LZE(ng)%ncid,   &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          LZE(ng)%Vid(idSbry(isUvel)),            &
     &                          LZE(ng)%Rindex, u3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D V-momentum component.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idVvel),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_v(:,:,:,nout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), LZE(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, iADM, LZE(ng)%name, LZE(ng)%ncid,   &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          LZE(ng)%Vid(idSbry(isVvel)),            &
     &                          LZE(ng)%Rindex, v3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        LZE(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Tid(itrc),   &
     &                     LZE(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_t(:,:,:,nout,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        LZE(ng)%Rindex
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
          status=nf_fwrite3d_bry (ng, iADM, LZE(ng)%name, LZE(ng)%ncid, &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            LZE(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                            LZE(ng)%Rindex, r3dvar,               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,  &
     &                                                    nout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          LZE(ng)%Rindex
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
          status=nf_fwrite3d(ng, iADM, LZE(ng)%ncid,                    &
     &                       LZE(ng)%Vid(idTsur(itrc)),                 &
     &                       LZE(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tl_tflux(:,:,:,kout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          LZE(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idUsms),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   FORCES(ng) % tl_ustr(:,:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), LZE(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iADM, LZE(ng)%ncid, LZE(ng)%Vid(idVsms),   &
     &                   LZE(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   FORCES(ng) % tl_vstr(:,:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), LZE(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize evolved vectors NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iADM, LZE(ng)%name, LZE(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, LZE(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, nout, LZE(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, LZE(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, LZE(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' WRT_EVOLVED_NF90 - error while writing variable: ',a, &
     &        /,20x,'into evolved LCZ NetCDF file for time record: ',i0)
# ifdef SOLVE3D
  20  FORMAT (2x,'WRT_EVOLVED_NF90 - wrote evolved LCZ',t40,            &
#  ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#  endif
# else
  20  FORMAT (2x,'WRT_EVOLVED_NF90 - wrote evolved LCZ',t40,            &
#  ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,')   in record = ',i0)
#  endif
# endif
!
      RETURN
      END SUBROUTINE wrt_evolved_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_evolved_pio (ng, tile, kout, nout,                 &
#  ifdef ADJUST_BOUNDARY
     &                            LBij, UBij,                           &
#  endif
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, kout, nout
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
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
     &  __FILE__//", wrt_evolved_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out evolved Lanczos vectors.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      LZE(ng)%Rindex=LZE(ng)%Rindex+1
      Fcount=LZE(ng)%Fcount
      LZE(ng)%Nrec(Fcount)=LZE(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, iADM, LZE(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/LZE(ng)%Rindex/), (/1/),              &
     &                          pioFile = LZE(ng)%pioFile,              &
     &                          pioVar = LZE(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idFsur),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   OCEAN(ng)% tl_zeta(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), LZE(ng)%Rindex
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
        IF (LZE(ng)%pioVar(idSbry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iADM, LZE(ng)%name,                 &
     &                          LZE(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          LZE(ng)%pioVar(idSbry(isFsur)),         &
     &                          LZE(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D U-momentum component.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idUbar),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_ubar(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), LZE(ng)%Rindex
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
        IF (LZE(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iADM, LZE(ng)%name,                 &
     &                          LZE(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          LZE(ng)%pioVar(idSbry(isUbar)),         &
     &                          LZE(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D V-momentum component.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idVbar),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_vbar(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), LZE(ng)%Rindex
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
        IF (LZE(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iADM, LZE(ng)%name,                 &
     &                          LZE(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          LZE(ng)%pioVar(idSbry(isVbar)),         &
     &                          LZE(ng)%Rindex, v2dvar,                 &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idUvel),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_u(:,:,:,nout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), LZE(ng)%Rindex
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
        IF (LZE(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry (ng, iADM, LZE(ng)%name,                 &
     &                          LZE(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          LZE(ng)%pioVar(idSbry(isUvel)),         &
     &                          LZE(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        LZE(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out 3D V-momentum component.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idVvel),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_v(:,:,:,nout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), LZE(ng)%Rindex
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
        IF (LZE(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry (ng, iADM, LZE(ng)%name,                 &
     &                          LZE(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          LZE(ng)%pioVar(idSbry(isVvel)),         &
     &                          LZE(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        LZE(ng)%Rindex
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
        IF (LZE(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iADM, LZE(ng)%pioFile,                   &
     &                     LZE(ng)%pioTrc(itrc),                        &
     &                     LZE(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % tl_t(:,:,:,nout,itrc))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        LZE(ng)%Rindex
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
          IF (LZE(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry (ng, iADM, LZE(ng)%name,               &
     &                            LZE(ng)%pioFile,                      &
     &                            Vname(1,ifield),                      &
     &                            LZE(ng)%pioVar(ifield),               &
     &                            LZE(ng)%Rindex,                       &
     &                            ioDesc,                               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,  &
     &                                                    nout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield), LZE(ng)%Rindex
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
          IF (LZE(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iADM, LZE(ng)%pioFile,                 &
     &                       LZE(ng)%pioVar(idTsur(itrc)),              &
     &                       LZE(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % tl_tflux(:,:,:,kout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          LZE(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idUsms),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   FORCES(ng) % tl_ustr(:,:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), LZE(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      IF (LZE(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, LZE(ng)%pioFile,                     &
     &                   LZE(ng)%pioVar(idVsms),                        &
     &                   LZE(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   FORCES(ng) % tl_vstr(:,:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), LZE(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize evolved vectors NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iADM, LZE(ng)%name, LZE(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, LZE(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, nout, LZE(ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, LZE(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, LZE(ng)%Rindex
#   endif
#  endif
!
  10  FORMAT (/,' WRT_EVOLVED_PIO - error while writing variable: ',a,  &
     &        /,20x,'into evolved LCZ NetCDF file for time record: ',i0)
#  ifdef SOLVE3D
  20  FORMAT (2x,'WRT_EVOLVED_PIO  - wrote evolved LCZ',t40,            &
#   ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#   endif
#  else
  20  FORMAT (2x,'WRT_EVOLVED_PIO  - wrote evolved LCZ',t40,            &
#   ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,')   in record = ',i0)
#   endif
#  endif
!
      RETURN
      END SUBROUTINE wrt_evolved_pio
# endif
#endif
      END MODULE wrt_evolved_mod

