#include "cppdefs.h"
      MODULE wrt_error_mod
#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
!
!git $Id$
!svn $Id: wrt_error.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes the full posterior error covariance (diagonal)   !
!  matrix for weak constraint 4DVar data assimilation using the        !
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
      PUBLIC  :: wrt_error
      PRIVATE :: wrt_error_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_error_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_error (ng, tile, kout, nout)
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
      SELECT CASE (ERR(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_error_nf90 (ng, tile, kout, nout,                    &
# ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
# endif
     &                         LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_error_pio (ng, tile, kout, nout,                     &
#  ifdef ADJUST_BOUNDARY
     &                        LBij, UBij,                               &
#  endif
     &                        LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) ERR(ng)%IOtype
          exit_flag=3
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
!
  10  FORMAT (' WRT_ERROR - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_error
!
!***********************************************************************
      SUBROUTINE wrt_error_nf90 (ng, tile, kout, nout,                  &
# ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
# endif
     &                           LBi, UBi, LBj, UBj)
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
     &  __FILE__//", wrt_error_nf90"
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
      ERR(ng)%Rindex=ERR(ng)%Rindex+1
      Fcount=ERR(ng)%Fcount
      ERR(ng)%Nrec(Fcount)=ERR(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/ERR(ng)%Rindex/), (/1/),                  &
     &                      ncid = ERR(ng)%ncid,                        &
     &                      varid = ERR(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out inner-loop Lanczos vectors tridiagonal system.
!
      CALL netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                     &
     &                      'zLanczos_coef', zLanczos_coef,             &
     &                      (/1,1/), (/Ninner,Ninner/),                 &
     &                      ncid = ERR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                     &
     &                      'zLanczos_inv', zLanczos_inv,               &
     &                      (/1,1/), (/Ninner,Ninner/),                 &
     &                      ncid = ERR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                     &
     &                      'zLanczos_err', zLanczos_err,               &
     &                      (/1,1/), (/Ninner,Ninner/),                 &
     &                      ncid = ERR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface error variance.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idFsur),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng)% tl_zeta(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries error variance.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid,   &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          ERR(ng)%Vid(idSbry(isFsur)),            &
     &                          ERR(ng)%Rindex, r2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D U-momentum component error variance.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idUbar),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % tl_ubar(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid,   &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          ERR(ng)%Vid(idSbry(isUbar)),            &
     &                          ERR(ng)%Rindex, u2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D V-momentum component error variance.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idVbar),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % tl_vbar(:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid,   &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          ERR(ng)%Vid(idSbry(isVbar)),            &
     &                          ERR(ng)%Rindex, v2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component error variance.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idUvel),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_u(:,:,:,nout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid,   &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          ERR(ng)%Vid(idSbry(isUvel)),            &
     &                          ERR(ng)%Rindex, u3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D V-momentum component error variance.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idVvel),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_v(:,:,:,nout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid,   &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          ERR(ng)%Vid(idSbry(isVvel)),            &
     &                          ERR(ng)%Rindex, v3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out tracer type variables error variance.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Tid(itrc),   &
     &                     ERR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_t(:,:,:,nout,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries error variance.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          status=nf_fwrite3d_bry (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid, &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            ERR(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                            ERR(ng)%Rindex, r3dvar,               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,  &
     &                                                    nout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          ERR(ng)%Rindex
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
!  Write out surface net tracers fluxes error variance. Notice that
!  fluxes have their own fixed time-dimension (of size Nfrec) to allow
!  4DVar adjustments at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iTLM, ERR(ng)%ncid,                    &
     &                       ERR(ng)%Vid(idTsur(itrc)),                 &
     &                       ERR(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tl_tflux(:,:,:,kout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          ERR(ng)%Rindex
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
!  Write out surface U-momentum stress error variance.  Notice that the
!  stress has its own fixed time-dimension (of size Nfrec) to allow
!  4DVar adjustments at other times in addition to initialization time.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idUsms),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % tl_ustr(:,:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress error variance.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%ncid, ERR(ng)%Vid(idVsms),   &
     &                   ERR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % tl_vstr(:,:,:,kout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize posterior error covariance NetCDF file to disk to allow
!  other processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iTLM, ERR(ng)%name, ERR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, ERR(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, nout, ERR(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, ERR(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, ERR(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' WRT_ERROR_NF90 - error while writing variable: ',a,   &
     &        /,18x,'into 4DVar error NetCDF file for time record: ',i0)
# ifdef SOLVE3D
  20  FORMAT (2x,'WRT_ERROR_NF90   - wrote error    fields (Index=', i1,&
#  ifdef NESTING
     &        ',',i1,') into time record = ',i0,2x,i2.2)
#  else
     &        ',',i1,') into time record = ',i0)
#  endif
# else
  20  FORMAT (2x,'WRT_ERROR_NF90   - wrote error    fields (Index=', i1,&
#  ifdef NESTING
     &        ') into time record = ',i0,2x,i2.2)
#  else
     &        ') into time record = ',i0)
#  endif
# endif
!
      RETURN
      END SUBROUTINE wrt_error_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_error_pio (ng, tile, kout, nout,                   &
#  ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
#  endif
     &                          LBi, UBi, LBj, UBj)
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
     &  __FILE__//", wrt_error_nf90"
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
      ERR(ng)%Rindex=ERR(ng)%Rindex+1
      Fcount=ERR(ng)%Fcount
      ERR(ng)%Nrec(Fcount)=ERR(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/ERR(ng)%Rindex/), (/1/),              &
     &                          pioFile = ERR(ng)%pioFile,              &
     &                          pioVar = ERR(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out inner-loop Lanczos vectors tridiagonal system.
!
      CALL pio_netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                 &
     &                          'zLanczos_coef', zLanczos_coef,         &
     &                          (/1,1/), (/Ninner,Ninner/),             &
     &                          pioFile = ERR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                 &
     &                          'zLanczos_inv', zLanczos_inv,           &
     &                          (/1,1/), (/Ninner,Ninner/),             &
     &                          piofile = ERR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      CALL pio_netcdf_put_fvar (ng, iTLM, ERR(ng)%name,                 &
     &                          'zLanczos_err', zLanczos_err,           &
     &                          (/1,1/), (/Ninner,Ninner/),             &
     &                          pioFile = ERR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface error variance.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idFsur),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   OCEAN(ng)% tl_zeta(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries error variance.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        IF (ERR(ng)%pioVar(idSbry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iTLM, ERR(ng)%name,                 &
     &                          ERR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          ERR(ng)%pioVar(idSbry(isFsur)),         &
     &                          ERR(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D U-momentum component error variance.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idUbar),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_ubar(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        IF (ERR(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iTLM, ERR(ng)%name,                 &
     &                          ERR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          ERR(ng)%pioVar(idSbry(isUbar)),         &
     &                          ERR(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D V-momentum component error variance.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fwrite2d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idVbar),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % tl_vbar(:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        IF (ERR(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry (ng, iTLM, ERR(ng)%name,                 &
     &                          ERR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          ERR(ng)%pioVar(idSbry(isVbar)),         &
     &                          ERR(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,   &
     &                                                     kout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component error variance.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idUvel),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_u(:,:,:,nout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        IF (ERR(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry (ng, iTLM, ERR(ng)%name,                 &
     &                          ERR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          ERR(ng)%pioVar(idSbry(isUvel)),         &
     &                          ERR(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out 3D V-momentum component error variance.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idVvel),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#   endif
     &                   OCEAN(ng) % tl_v(:,:,:,nout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries error variance.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        IF (ERR(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry (ng, iTLM, ERR(ng)%name,                 &
     &                          ERR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          ERR(ng)%pioVar(idSbry(isVvel)),         &
     &                          ERR(ng)%Rindex,                         &
     &                          ioDesc,                                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,    &
     &                                                  nout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out tracer type variables error variance.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (ERR(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iTLM, ERR(ng)%pioFile,                   &
     &                     ERR(ng)%pioTrc(itrc),                        &
     &                     ERR(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % tl_t(:,:,:,nout,itrc))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        ERR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#   ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries error variance.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          IF (ERR(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry (ng, iTLM, ERR(ng)%name,               &
     &                            ERR(ng)%pioFile,                      &
     &                            Vname(1,ifield),                      &
     &                            ERR(ng)%pioVar(ifield),               &
     &                            ERR(ng)%Rindex,                       &
     &                            ioDesc,                               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,  &
     &                                                    nout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield), ERR(ng)%Rindex
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
!  Write out surface net tracers fluxes error variance. Notice that
!  fluxes have their own fixed time-dimension (of size Nfrec) to allow
!  4DVar adjustments at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp
          IF (ERR(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iTLM, ERR(ng)%pioFile,                 &
     &                       ERR(ng)%pioVar(idTsur(itrc)),              &
     &                       ERR(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % tl_tflux(:,:,:,kout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          ERR(ng)%Rindex
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
!  Write out surface U-momentum stress error variance.  Notice that the
!  stress has its own fixed time-dimension (of size Nfrec) to allow
!  4DVar adjustments at other times in addition to initialization time.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idUsms),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % tl_ustr(:,:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress error variance.
!
      scale=1.0_dp
      IF (ERR(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iTLM, ERR(ng)%pioFile,                     &
     &                   ERR(ng)%pioVar(idVsms),                        &
     &                   ERR(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % tl_vstr(:,:,:,kout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), ERR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize posterior error covariance NetCDF file to disk to allow
!  other processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iTLM, ERR(ng)%name, ERR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, ERR(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, nout, ERR(ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, ERR(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, ERR(ng)%Rindex
#   endif
#  endif
!
  10  FORMAT (/,' WRT_ERROR_PIO - error while writing variable: ',a,    &
     &        /,18x,'into 4DVar error NetCDF file for time record: ',i0)
#  ifdef SOLVE3D
  20  FORMAT (2x,'WRT_ERROR_PIO    - wrote error    fields (Index=', i1,&
#   ifdef NESTING
     &        ',',i1,') into time record = ',i0,2x,i2.2)
#   else
     &        ',',i1,') into time record = ',i0)
#   endif
#  else
  20  FORMAT (2x,'WRT_ERROR_PIO    - wrote error    fields (Index=', i1,&
#   ifdef NESTING
     &        ') into time record = ',i0,2x,i2.2)
#   else
     &        ') into time record = ',i0)
#   endif
#  endif
!
      RETURN
      END SUBROUTINE wrt_error_pio
# endif
#endif
      END MODULE wrt_error_mod

