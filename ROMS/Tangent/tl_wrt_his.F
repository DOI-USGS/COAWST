#include "cppdefs.h"
      MODULE tl_wrt_his_mod
#if defined TANGENT || defined TL_IOMS
!
!git $Id$
!svn $Id: tl_wrt_his.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes fields into output tangent linear history file   !
!  using either the standard NetCDF library or the Parallel-IO (PIO)   !
!  library.                                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS || \
     defined FORCING_SV    || defined STOCHASTIC_OPT || \
     defined HESSIAN_SO    || defined HESSIAN_FSV
      USE mod_forces
# endif
      USE mod_grid
      USE mod_iounits
# ifdef SOLVE3D
      USE mod_mixing
# endif
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
      PUBLIC  :: tl_wrt_his
      PRIVATE :: tl_wrt_his_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: tl_wrt_his_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tl_wrt_his (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
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
      SELECT CASE (TLM(ng)%IOtype)
        CASE (io_nf90)
          CALL tl_wrt_his_nf90 (ng, tile,                               &
# ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
# endif
     &                          LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL tl_wrt_his_pio (ng, tile,                                &
# ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
# endif
     &                         LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TLM(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' TL_WRT_HIS - Illegal output file type, io_type = ',i0,  &
     &        /,14x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE tl_wrt_his
!
!***********************************************************************
      SUBROUTINE tl_wrt_his_nf90 (ng, tile,                             &
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
      integer, intent(in) :: ng, tile
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k
# endif
!
      real(dp) :: scale
      real(r8) :: Tval(1)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tl_wrt_his_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out tangent linear fields.
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
      TLM(ng)%Rindex=TLM(ng)%Rindex+1
      Fcount=TLM(ng)%load
      TLM(ng)%Nrec(Fcount)=TLM(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in the tangent
!  linear file.
!
      IF (LcycleTLM(ng)) THEN
        TLM(ng)%Rindex=MOD(TLM(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtPER(ng)) THEN
        Tval(1)=REAL(TLM(ng)%Rindex,r8)*day2sec
      ELSE
        Tval(1)=time(ng)
      END IF
      CALL netcdf_put_fvar (ng, iTLM, TLM(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), tval,                &
     &                      (/TLM(ng)%Rindex/), (/1/),                  &
     &                      ncid = TLM(ng)%ncid,                        &
     &                      varid = TLM(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp                          ! m2/s2
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idUsms),   &
     &                   TLM(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % tl_ustr(:,:,:,Lfout(ng)))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp                          ! m2/s2
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVsms),   &
     &                   TLM(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % tl_vstr(:,:,:,Lfout(ng)))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp                      ! kinematic flux units
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid,                    &
     &                       TLM(ng)%Vid(idTsur(itrc)),                 &
     &                       TLM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng)% tl_tflux(:,:,:,Lfout(ng),itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))), Lfout(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
# if defined FORCING_SV || defined STOCHASTIC_OPT || \
     defined HESSIAN_SO    || defined HESSIAN_FSV
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idUsms), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % tl_sustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVsms), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % tl_svstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Write out net surface active tracer fluxes.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTsur(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid,                    &
     &                       TLM(ng)%Vid(idTsur(itrc)),                 &
     &                       TLM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tl_stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          TLM(ng)%Rindex
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
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idFsur), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     OCEAN(ng) % tl_zeta(:,:,KOUT),               &
     &                     SetFillVal = .FALSE.)
# else
#  ifdef FORCING_SV
     &                     OCEAN(ng) % f_zeta(:,:))
#  else
     &                     OCEAN(ng) % tl_zeta(:,:,KOUT))
#  endif
# endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# if defined FORWARD_WRITE && defined FORWARD_RHS
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRzet), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_rzeta(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRzet)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iTLM, TLM(ng)%name, TLM(ng)%ncid,    &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         TLM(ng)%Vid(idSbry(isFsur)),             &
     &                         TLM(ng)%Rindex, r2dvar,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        TLM(ng)%Rindex
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
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idUbar), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
# ifdef FORCING_SV
     &                     OCEAN(ng) % f_ubar(:,:))
# else
     &                     OCEAN(ng) % tl_ubar(:,:,KOUT))
# endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRu2d), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_rubar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu2d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRuct), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     COUPLING(ng) % tl_rufrc)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRuct)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idUfx1), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DU_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idUfx2), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DU_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iTLM, TLM(ng)%name, TLM(ng)%ncid,    &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         TLM(ng)%Vid(idSbry(isUbar)),             &
     &                         TLM(ng)%Rindex, u2dvar,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        TLM(ng)%Rindex
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
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVbar), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
# ifdef FORCING_SV
     &                     OCEAN(ng) % f_vbar(:,:))
# else
     &                     OCEAN(ng) % tl_vbar(:,:,KOUT))
# endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRv2d), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_rvbar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv2d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRvct), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     COUPLING(ng) % tl_rvfrc)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRvct)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVfx1), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DV_avg1)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVfx2), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DV_avg2)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iTLM, TLM(ng)%name, TLM(ng)%ncid,    &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         TLM(ng)%Vid(idSbry(isVbar)),             &
     &                         TLM(ng)%Rindex, v2dvar,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idUvel), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
#  ifdef FORCING_SV
     &                     OCEAN(ng) % f_u(:,:,:))
#  else
     &                     OCEAN(ng) % tl_u(:,:,:,NOUT))
#  endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRu3d), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_ru(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu3d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iTLM, TLM(ng)%name, TLM(ng)%ncid,    &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         TLM(ng)%Vid(idSbry(isUvel)),             &
     &                         TLM(ng)%Rindex, u3dvar,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        TLM(ng)%Rindex
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
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVvel), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
#  ifdef FORCING_SV
     &                     OCEAN(ng) % f_v(:,:,:))
#  else
     &                     OCEAN(ng) % tl_v(:,:,:,NOUT))
#  endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idRv3d), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_rv(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv3d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iTLM, TLM(ng)%name, TLM(ng)%ncid,    &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         TLM(ng)%Vid(idSbry(isVvel)),             &
     &                         TLM(ng)%Rindex, v3dvar,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        TLM(ng)%Rindex
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
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Tid(itrc), &
     &                       TLM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
#  ifdef FORCING_SV
     &                       OCEAN(ng) % f_t(:,:,:,itrc))
#  else
     &                       OCEAN(ng) % tl_t(:,:,:,NOUT,itrc))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          TLM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          status=nf_fwrite3d_bry(ng, iTLM, TLM(ng)%name, TLM(ng)%ncid,  &
     &                           Vname(1,idSbry(isTvar(itrc))),         &
     &                           TLM(ng)%Vid(idSbry(isTvar(itrc))),     &
     &                           TLM(ng)%Rindex, r3dvar,                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,   &
     &                                                 Lbout(ng),itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          TLM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idDano), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_rho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined FORWARD_MIXING && \
     (defined BVF_MIXING     || defined GLS_MIXING  || \
      defined LMD_MIXING     || defined MY25_MIXING)
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVvis), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
#   ifdef WEAK_CONSTRAINT
     &                     MIXING(ng) % Akv,                            &
#   else
     &                     MIXING(ng) % tl_Akv,                         &
#   endif
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idTdif), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
#   ifdef WEAK_CONSTRAINT
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
#   else
     &                     MIXING(ng) % tl_Akt(:,:,:,itemp),            &
#   endif
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idSdif), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
#    ifdef WEAK_CONSTRAINT
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
#    else
     &                     MIXING(ng) % tl_Akt(:,:,:,isalt),            &
#    endif
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idMtke), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_tke(:,:,:,NOUT),             &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVmKK), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_Akk)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKK)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idMtls), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_gls(:,:,:,NOUT),             &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVmLS), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_Lscale)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmLS)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#    ifdef GSL_MIXING
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%ncid, TLM(ng)%Vid(idVmKP), &
     &                     TLM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     MIXING(ng) % tl_Akp)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKP)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iTLM, TLM(ng)%name, TLM(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) KOUT, NOUT, TLM(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) KOUT, NOUT, TLM(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) KOUT, TLM(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) KOUT, TLM(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' TL_WRT_HIS_NF90 - error while writing variable: ',a,  &
     &        /,19x,'into tangent NetCDF file for time record: ',i0)
  20  FORMAT (2x,'TL_WRT_HIS_NF90  - wrote history', t40,               &
# ifdef SOLVE3D
#  ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#  endif
# else
#  ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,')   in record = ',i0)
#  endif
# endif
!
      RETURN
      END SUBROUTINE tl_wrt_his_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE tl_wrt_his_pio (ng, tile,                              &
# ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
# endif
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, ifield, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k
# endif
!
      real(dp) :: scale
      real(r8) :: Tval(1)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", tl_wrt_his_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out tangent linear fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      TLM(ng)%Rindex=TLM(ng)%Rindex+1
      Fcount=TLM(ng)%load
      TLM(ng)%Nrec(Fcount)=TLM(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in the tangent
!  linear file.
!
      IF (LcycleTLM(ng)) THEN
        TLM(ng)%Rindex=MOD(TLM(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtPER(ng)) THEN
        Tval(1)=REAL(TLM(ng)%Rindex,r8)*day2sec
      ELSE
        Tval(1)=time(ng)
      END IF
      CALL pio_netcdf_put_fvar (ng, iTLM, TLM(ng)%name,                 &
     &                          TRIM(Vname(1,idtime)), tval,            &
     &                          (/TLM(ng)%Rindex/), (/1/),              &
     &                           pioFile = TLM(ng)%pioFile,             &
     &                           pioVar = TLM(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp                          ! m2/s2
      IF (TLM(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                     &
     &                   TLM(ng)%pioVar(idUsms),                        &
     &                   TLM(ng)%Rindex, ioDesc,                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % tl_ustr(:,:,:,Lfout(ng)))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp                          ! m2/s2
      IF (TLM(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                     &
     &                   TLM(ng)%pioVar(idVsms),                        &
     &                   TLM(ng)%Rindex, ioDesc,                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % tl_vstr(:,:,:,Lfout(ng)))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp                      ! kinematic flux units
          IF (TLM(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                 &
     &                       TLM(ng)%pioVar(idTsur(itrc)),              &
     &                       TLM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng)% tl_tflux(:,:,:,Lfout(ng),itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))), Lfout(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
# if defined FORCING_SV || defined STOCHASTIC_OPT || \
     defined HESSIAN_SO    || defined HESSIAN_FSV
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idUsms),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % tl_sustr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVsms),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % tl_svstr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SOLVE3D
!
!  Write out net surface active tracer fluxes.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTsur(itrc),ng)) THEN
          scale=1.0_dp
          IF (TLM(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
!
          status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                 &
     &                       TLM(ng)%pioVar(idTsur(itrc)),              &
     &                       TLM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tl_stflx(:,:,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          TLM(ng)%Rindex
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
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idFsur),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     OCEAN(ng) % tl_zeta(:,:,KOUT),               &
     &                     SetFillVal = .FALSE.)
# else
#  ifdef FORCING_SV
     &                     OCEAN(ng) % f_zeta(:,:))
#  else
     &                     OCEAN(ng) % tl_zeta(:,:,KOUT))
#  endif
# endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRzet)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRzet),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_rzeta(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRzet)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isFsur)
        IF (TLM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iTLM, TLM(ng)%name,                  &
     &                         TLM(ng)%pioFile,                         &
     &                         Vname(1,ifield),                         &
     &                         TLM(ng)%pioVar(ifield),                  &
     &                         TLM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), TLM(ng)%Rindex
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
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idUbar),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
# ifdef FORCING_SV
     &                     OCEAN(ng) % f_ubar(:,:))
# else
     &                     OCEAN(ng) % tl_ubar(:,:,KOUT))
# endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRu2d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRu2d),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_rubar(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu2d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRuct)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRuct),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     COUPLING(ng) % tl_rufrc)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRuct)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        IF (TLM(ng)%pioVar(idUfx1)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idUfx1),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DU_avg1)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (TLM(ng)%pioVar(idUfx2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idUfx2),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DU_avg2)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isUbar)
        IF (TLM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF

        status=nf_fwrite2d_bry(ng, iTLM, TLM(ng)%name,                  &
     &                         TLM(ng)%pioFile,                         &
     &                         Vname(1,ifield),                         &
     &                         TLM(ng)%pioVar(ifield),                  &
     &                         TLM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), TLM(ng)%Rindex
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
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVbar),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
# ifdef FORCING_SV
     &                     OCEAN(ng) % f_vbar(:,:))
# else
     &                     OCEAN(ng) % tl_vbar(:,:,KOUT))
# endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRv2d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRv2d),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_rvbar(:,:,KOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv2d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRvct)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRvct),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     COUPLING(ng) % tl_rvfrc)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRvct)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        IF (TLM(ng)%pioVar(idVfx1)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVfx1),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DV_avg1)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (TLM(ng)%pioVar(idVfx2)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF

        status=nf_fwrite2d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVfx2),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     COUPLING(ng) % tl_DV_avg2)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isVbar)
        IF (TLM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iTLM, TLM(ng)%name,                  &
     &                         TLM(ng)%pioFile,                         &
     &                         Vname(1,ifield),                         &
     &                         TLM(ng)%pioVar(ifield),                  &
     &                         TLM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idUvel),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
#  ifdef FORCING_SV
     &                     OCEAN(ng) % f_u(:,:,:))
#  else
     &                     OCEAN(ng) % tl_u(:,:,:,NOUT))
#  endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRu3d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRu3d),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_ru(:,:,:,NOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu3d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isUvel)
        IF (TLM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iTLM, TLM(ng)%name,                  &
     &                         TLM(ng)%pioFile,                         &
     &                         Vname(1,ifield),                         &
     &                         TLM(ng)%pioVar(ifield),                  &
     &                         TLM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), TLM(ng)%Rindex
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
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVvel),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
#  ifdef FORCING_SV
     &                     OCEAN(ng) % f_v(:,:,:))
#  else
     &                     OCEAN(ng) % tl_v(:,:,:,NOUT))
#  endif
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (TLM(ng)%pioVar(idRv3d)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idRv3d),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % tl_rv(:,:,:,NOUT))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv3d)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isVvel)
        IF (TLM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iTLM, TLM(ng)%name,                  &
     &                         TLM(ng)%pioFile,                         &
     &                         Vname(1,ifield),                         &
     &                         TLM(ng)%pioVar(ifield),                  &
     &                         TLM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % tl_v_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), TLM(ng)%Rindex
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
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          IF (TLM(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF

          status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                 &
     &                       TLM(ng)%pioTrc(itrc),                      &
     &                       TLM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
#  ifdef FORCING_SV
     &                       OCEAN(ng) % f_t(:,:,:,itrc))
#  else
     &                       OCEAN(ng) % tl_t(:,:,:,NOUT,itrc))
#  endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          TLM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          IF (TLM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry(ng, iTLM, TLM(ng)%name,                &
     &                           TLM(ng)%pioFile,                       &
     &                           Vname(1,ifield),                       &
     &                           TLM(ng)%pioVar(ifield),                &
     &                           TLM(ng)%Rindex, ioDesc,                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,   &
     &                                                 Lbout(ng),itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), TLM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idDano),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % tl_rho)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined FORWARD_MIXING && \
     (defined BVF_MIXING     || defined GLS_MIXING  || \
      defined LMD_MIXING     || defined MY25_MIXING)
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVvis),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
#   ifdef WEAK_CONSTRAINT
     &                     MIXING(ng) % Akv,                            &
#   else
     &                     MIXING(ng) % tl_Akv,                         &
#   endif
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idTdif),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
#   ifdef WEAK_CONSTRAINT
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
#   else
     &                     MIXING(ng) % tl_Akt(:,:,:,itemp),            &
#   endif
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idSdif),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
#    ifdef WEAK_CONSTRAINT
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
#    else
     &                     MIXING(ng) % tl_Akt(:,:,:,isalt),            &
#    endif
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined GLS_MIXING_NOT_YET || defined MY25_MIXING_NOT_YET
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idMtke)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idMtke),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_tke(:,:,:,NOUT),             &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVmKK)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVmKK),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_Akk)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKK)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idMtld)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idMtls),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_gls(:,:,:,NOUT),             &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVmLS)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVmLS),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tl_Lscale)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmLS)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#    ifdef GSL_MIXING
!
        scale=1.0_dp
        IF (TLM(ng)%pioVar(idVmKP)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF

        status=nf_fwrite3d(ng, iTLM, TLM(ng)%pioFile,                   &
     &                     TLM(ng)%pioVar(idVmKP),                      &
     &                     TLM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     MIXING(ng) % tl_Akp)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKP)), TLM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iTLM, TLM(ng)%name, TLM(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) KOUT, NOUT, TLM(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) KOUT, NOUT, TLM(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) KOUT, TLM(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) KOUT, TLM(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' TL_WRT_HIS_PIO - error while writing variable: ',a,   &
     &        /,19x,'into tangent NetCDF file for time record: ',i0)
  20  FORMAT (2x,'TL_WRT_HIS_PIO   - wrote history', t40,               &
# ifdef SOLVE3D
#  ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#  endif
# else
#  ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#  else
     &        'fields (Index=',i1,')   in record = ',i0)
#  endif
# endif
!
      RETURN
      END SUBROUTINE tl_wrt_his_pio
# endif
#endif
      END MODULE tl_wrt_his_mod
