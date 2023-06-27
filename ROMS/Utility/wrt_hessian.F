#include "cppdefs.h"
      MODULE wrt_hessian_mod
#if defined I4DVAR                || \
    (defined WEAK_CONSTRAINT      && \
     (defined POSTERIOR_EOFS      || defined POSTERIOR_ERROR_I ||\
      defined POSTERIOR_ERROR_F)) || defined LCZ_FINAL
!
!git $Id$
!svn $Id: wrt_hessian.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes Hessian eigenvectors into output file using the  !
!  standard NetCDF library or the Parallel-IO (PIO) library.           !
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
      PUBLIC  :: wrt_hessian
      PRIVATE :: wrt_hessian_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_hessian_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_hessian (ng, tile, kout, nout)
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
      SELECT CASE (HSS(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_hessian_nf90 (ng, tile, kout, nout,                  &
# ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
# endif
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_hessian_pio (ng, tile, kout, nout,                   &
#  ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
#  endif
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) HSS(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_HESSIAN - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_hessian
!
!***********************************************************************
      SUBROUTINE wrt_hessian_nf90 (ng, tile, kout, nout,                &
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
     &  __FILE__//", wrt_hessian_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Hessian eigenvectors.
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
      HSS(ng)%Rindex=HSS(ng)%Rindex+1
      Fcount=HSS(ng)%Fcount
      HSS(ng)%Nrec(Fcount)=HSS(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iADM, HSS(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/HSS(ng)%Rindex/), (/1/),                  &
     &                      ncid = HSS(ng)%ncid,                        &
     &                      varid = HSS(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idFsur),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng)% ad_zeta(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), HSS(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, HSS(ng)%name, HSS(ng)%ncid,   &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          HSS(ng)%Vid(idSbry(isFsur)),            &
     &                          HSS(ng)%Rindex, r2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        HSS(ng)%Rindex
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
      status=nf_fwrite2d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idUbar),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % ad_ubar(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), HSS(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, HSS(ng)%name, HSS(ng)%ncid,   &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          HSS(ng)%Vid(idSbry(isUbar)),            &
     &                          HSS(ng)%Rindex, u2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        HSS(ng)%Rindex
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
      status=nf_fwrite2d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idVbar),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % ad_vbar(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), HSS(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, iADM, HSS(ng)%name, HSS(ng)%ncid,   &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          HSS(ng)%Vid(idSbry(isVbar)),            &
     &                          HSS(ng)%Rindex, v2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        HSS(ng)%Rindex
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
      status=nf_fwrite3d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idUvel),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % ad_u(:,:,:,nout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), HSS(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, iADM, HSS(ng)%name, HSS(ng)%ncid,   &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          HSS(ng)%Vid(idSbry(isUvel)),            &
     &                          HSS(ng)%Rindex, u3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        HSS(ng)%Rindex
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
      status=nf_fwrite3d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idVvel),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % ad_v(:,:,:,nout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), HSS(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, iADM, HSS(ng)%name, HSS(ng)%ncid,   &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          HSS(ng)%Vid(idSbry(isVvel)),            &
     &                          HSS(ng)%Rindex, v3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        HSS(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Tid(itrc),   &
     &                     HSS(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % ad_t(:,:,:,nout,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        HSS(ng)%Rindex
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
          status=nf_fwrite3d_bry (ng, iADM, HSS(ng)%name, HSS(ng)%ncid, &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            HSS(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                            HSS(ng)%Rindex, r3dvar,               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,  &
     &                                                    nout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          HSS(ng)%Rindex
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
          status=nf_fwrite3d(ng, iADM, HSS(ng)%ncid,                    &
     &                       HSS(ng)%Vid(idTsur(itrc)),                 &
     &                       HSS(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_tflux(:,:,:,kout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HSS(ng)%Rindex
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
      status=nf_fwrite3d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idUsms),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % ad_ustr(:,:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), HSS(ng)%Rindex
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
      status=nf_fwrite3d(ng, iADM, HSS(ng)%ncid, HSS(ng)%Vid(idVsms),   &
     &                   HSS(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % ad_vstr(:,:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), HSS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize Hessian eigenvectors NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iADM, HSS(ng)%name, HSS(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, HSS(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, nout, HSS(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, HSS(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, HSS(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' WRT_HESSIAN_NF90 - error while writing variable: ',a, &
     &        /,20x,'into Hessian NetCDF file for time record: ',i0)
# ifdef SOLVE3D
  20  FORMAT (2x,'WRT_HESSIAN_NF90 - wrote Hessian', t40,               &
#  ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#  endif
# else
  20  FORMAT (2x,'WRT_HESSIAN_NF90 - wrote Hessian', t40,               &
#  ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,')   in record = ',i0)
#  endif
# endif
!
      RETURN
      END SUBROUTINE wrt_hessian_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_hessian_pio (ng, tile, kout, nout,                 &
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
     &  __FILE__//", wrt_hessian_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out Hessian eigenvectors.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      HSS(ng)%Rindex=HSS(ng)%Rindex+1
      Fcount=HSS(ng)%Fcount
      HSS(ng)%Nrec(Fcount)=HSS(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, iADM, HSS(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/HSS(ng)%Rindex/), (/1/),              &
     &                          pioFile = HSS(ng)%pioFile,              &
     &                          pioVar = HSS(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface.
!
      scale=1.0_dp
      IF (HSS(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idFsur),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   OCEAN(ng)% ad_zeta(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), HSS(ng)%Rindex
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
        IF (HSS(ng)%pioVar(idSbry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iADM, HSS(ng)%name,                 &
     &                          HSS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          HSS(ng)%pioVar(idSbry(isFsur)),         &
     &                          HSS(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        HSS(ng)%Rindex
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
      IF (HSS(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idUbar),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % ad_ubar(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), HSS(ng)%Rindex
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
        IF (HSS(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iADM, HSS(ng)%name,                 &
     &                          HSS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          HSS(ng)%pioVar(idSbry(isUbar)),         &
     &                          HSS(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        HSS(ng)%Rindex
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
      IF (HSS(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idVbar),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % ad_vbar(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), HSS(ng)%Rindex
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
        IF (HSS(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iADM, HSS(ng)%name,                 &
     &                          HSS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          HSS(ng)%pioVar(idSbry(isVbar)),         &
     &                          HSS(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        HSS(ng)%Rindex
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
      IF (HSS(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idUvel),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % ad_u(:,:,:,nout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), HSS(ng)%Rindex
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
        IF (HSS(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry (ng, iADM, HSS(ng)%name,                 &
     &                          HSS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          HSS(ng)%pioVar(idSbry(isUvel)),         &
     &                          HSS(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        HSS(ng)%Rindex
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
      IF (HSS(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idVvel),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % ad_v(:,:,:,nout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), HSS(ng)%Rindex
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
        IF (HSS(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry (ng, iADM, HSS(ng)%name,                 &
     &                          HSS(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          HSS(ng)%pioVar(idSbry(isVvel)),         &
     &                          HSS(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        HSS(ng)%Rindex
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
        IF (HSS(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iADM, HSS(ng)%pioFile,                   &
     &                     HSS(ng)%pioTrc(itrc),                        &
     &                     HSS(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % ad_t(:,:,:,nout,itrc))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        HSS(ng)%Rindex
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
          IF (HSS(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry (ng, iADM, HSS(ng)%name,               &
     &                            HSS(ng)%pioFile,                      &
     &                            Vname(1,ifield),                      &
     &                            HSS(ng)%pioVar(ifield),               &
     &                            HSS(ng)%Rindex,                       &
     &                            ioDesc,                               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,  &
     &                                                    nout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), HSS(ng)%Rindex
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
          IF (HSS(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iADM, HSS(ng)%pioFile,                 &
     &                       HSS(ng)%pioVar(idTsur(itrc)),              &
     &                       HSS(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % ad_tflux(:,:,:,kout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          HSS(ng)%Rindex
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
      IF (HSS(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idUsms),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ad_ustr(:,:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), HSS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      IF (HSS(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, HSS(ng)%pioFile,                     &
     &                   HSS(ng)%pioVar(idVsms),                        &
     &                   HSS(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % ad_vstr(:,:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), HSS(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize Hessian eigenvectors NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iADM, HSS(ng)%name, HSS(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, HSS(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, nout, HSS(ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, HSS(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, HSS(ng)%Rindex
#   endif
#  endif
!
  10  FORMAT (/,' WRT_HESSIAN_PIO - error while writing variable: ',a,  &
     &        /,19x,'into Hessian NetCDF file for time record: ',i0)
#  ifdef SOLVE3D
  20  FORMAT (2x,'WRT_HESSIAN_PIO  - wrote Hessian', t40,               &
#   ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#   endif
#  else
  20  FORMAT (2x,'WRT_HESSIAN_PIO - wrote Hessian', t40,                &
#   ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,')   in record = ',i0)
#   endif
#  endif
!
      RETURN
      END SUBROUTINE wrt_hessian_pio
# endif
#endif
      END MODULE wrt_hessian_mod
