#include "cppdefs.h"
      MODULE ad_wrt_his_mod
#ifdef ADJOINT
!
!git $Id$
!svn $Id: ad_wrt_his.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested adjoint model fields into adjoint     !
!  history NetCDF file.                                                !
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
# ifdef WEAK_CONSTRAINT
      USE mod_fourdvar
# endif
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
      USE omega_mod,           ONLY : scale_omega
# endif
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: ad_wrt_his
      PRIVATE :: ad_wrt_his_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: ad_wrt_his_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ad_wrt_his (ng, tile)
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
      SELECT CASE (ADM(ng)%IOtype)
        CASE (io_nf90)
          CALL ad_wrt_his_nf90 (ng, tile,                               &
# ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
# endif
     &                          LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL ad_wrt_his_pio (ng, tile,                                &
#  ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
#  endif
     &                         LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) ADM(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' AD_WRT_HIS - Illegal output type, io_type = ',i0,       &
     &        /,14x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE ad_wrt_his
!
!***********************************************************************
      SUBROUTINE ad_wrt_his_nf90 (ng, tile,                             &
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
      integer :: Fcount, i, j, gfactor, gtype, status
      integer :: kout
# ifdef WEAK_CONSTRAINT
      integer :: kfout
# endif
# ifdef SOLVE3D
      integer :: itrc, k, nout
# endif
!
      real(dp) :: scale
      real(r8) :: Tval(1)
#ifdef SOLVE3D
!
      real(r8), allocatable :: Wr3d(:,:,:)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_wrt_his_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out adjoint fields.
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
!  Determine time index to write. The "nout" index is updated to the
!  version of "ad_main3d" that updates the "iic" counter at the bottom.
!  Therefore, we need to change the conditional "iic(ng).ne.ntend(ng)"
!  to "iic(ng).gt.ntend(ng)" to get identical solutions.
!
# ifdef SOLVE3D
      kout=kstp(ng)
# else
      kout=kstp(ng)
# endif
# if defined WEAK_CONSTRAINT
      kfout=2
# endif
# ifdef SOLVE3D
      IF (iic(ng).gt.ntend(ng)) THEN
        nout=nnew(ng)
# ifdef AD_OUTPUT_STATE
        LwrtState3d(ng)=.FALSE.
# endif
      ELSE
# ifdef AD_OUTPUT_STATE
        LwrtState3d(ng)=.TRUE.
# endif
        nout=nstp(ng)
      END IF
# endif
!
!  Set time record index.
!
      ADM(ng)%Rindex=ADM(ng)%Rindex+1
      Fcount=ADM(ng)%load
      ADM(ng)%Nrec(Fcount)=ADM(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in the adjoint
!  file.
!
      IF (LcycleADJ(ng)) THEN
        ADM(ng)%Rindex=MOD(ADM(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtTime(ng)) THEN
        IF (LwrtPER(ng)) THEN
          Tval(1)=REAL(ADM(ng)%Rindex,r8)*day2sec
        ELSE
# if defined WEAK_CONSTRAINT && !defined WEAK_NOINTERP
          Tval(1)=ForceTime(ng)
# else
          Tval(1)=time(ng)
# endif
        END IF
        CALL netcdf_put_fvar (ng, iADM, ADM(ng)%name,                   &
     &                        TRIM(Vname(1,idtime)), tval,              &
     &                        (/ADM(ng)%Rindex/), (/1/),                &
     &                        ncid = ADM(ng)%ncid,                      &
     &                        varid = ADM(ng)%Vid(idtime))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp                          ! m2/s2
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idUsms),   &
     &                   ADM(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Lfout(ng)))
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
      status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVsms),   &
     &                   ADM(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Lfout(ng)))
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
!  Write out surface net heat flux. Notice that different tracer fluxes
!  are written at their own fixed time-dimension (of size Nfrec) to
!  allow 4DVAR adjustments at other times in addition to initial time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp                      ! kinematic flux units
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idTsur(itrc)),                 &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng)% ad_tflux(:,:,:,Lfout(ng),itrc))
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
!
!  Write out bathymetry.
!
      IF (Hout(idbath,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idbath), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng)% ad_h,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      IF (Hout(idFsur,ng)) THEN
# ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idFsur),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng)% f_zetaG(:,:,kfout))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
# endif
          scale=1.0_dp
          gtype=gfactor*r2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idFsur),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
# ifdef WET_DRY
     &                         OCEAN(ng)% ad_zeta(:,:,kout),            &
     &                         SetFillVal = .FALSE.)
# else
     &                         OCEAN(ng)% ad_zeta(:,:,kout))
# endif
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idFsur),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
# ifdef WET_DRY
     &                         OCEAN(ng)% ad_zeta_sol,                  &
     &                         SetFillVal = .FALSE.)
# else
     &                         OCEAN(ng)% ad_zeta_sol)
# endif
          ENDIF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
# ifdef WEAK_CONSTRAINT
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,    &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         ADM(ng)%Vid(idSbry(isFsur)),             &
     &                         ADM(ng)%Rindex, r2dvar,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        ADM(ng)%Rindex
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
# ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idUbar),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#  endif
     &                       OCEAN(ng) % f_ubarG(:,:,kfout))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
# endif
          scale=1.0_dp
          gtype=gfactor*u2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idUbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
# endif
     &                         OCEAN(ng) % ad_ubar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idUbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
# endif
     &                         OCEAN(ng) % ad_ubar_sol)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
# ifdef WEAK_CONSTRAINT
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,    &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         ADM(ng)%Vid(idSbry(isUbar)),             &
     &                         ADM(ng)%Rindex, u2dvar,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        ADM(ng)%Rindex
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
# ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idVbar),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#  endif
     &                       OCEAN(ng) % f_vbarG(:,:,kfout))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
# endif
          scale=1.0_dp
          gtype=gfactor*v2dvar
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idVbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         OCEAN(ng) % ad_vbar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idVbar),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, scale,               &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         OCEAN(ng) % ad_vbar_sol)
          END IF
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
# ifdef WEAK_CONSTRAINT
        END IF
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,    &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         ADM(ng)%Vid(idSbry(isVbar)),             &
     &                         ADM(ng)%Rindex, v2dvar,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        ADM(ng)%Rindex
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
#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_dp
          gtype=gfactor*u3dvar
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idUvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#   endif
     &                       OCEAN(ng) % f_uG(:,:,:,kfout))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          scale=1.0_dp
          gtype=gfactor*u3dvar
#  ifdef AD_OUTPUT_STATE
          IF (LwrtState3d(ng)) THEN
#  endif
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idUvel),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#    ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#    endif
     &                         OCEAN(ng) % ad_u(:,:,:,nout))
#  ifdef AD_OUTPUT_STATE
          ELSE
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idUvel),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#    ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#    endif
     &                         OCEAN(ng) % ad_u_sol(:,:,:))
          ENDIF
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,    &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         ADM(ng)%Vid(idSbry(isUvel)),             &
     &                         ADM(ng)%Rindex, u3dvar,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        ADM(ng)%Rindex
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
#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          scale=1.0_dp
          gtype=gfactor*v3dvar
          status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idVvel),                       &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       OCEAN(ng) % f_vG(:,:,:,kfout))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          scale=1.0_dp
          gtype=gfactor*v3dvar
#  ifdef AD_OUTPUT_STATE
          IF (LwrtState3d(ng)) THEN
#  endif
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idVvel),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#   endif
     &                         OCEAN(ng) % ad_v(:,:,:,nout))
#  ifdef AD_OUTPUT_STATE
          ELSE
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Vid(idVvel),                     &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#   endif
     &                         OCEAN(ng) % ad_v_sol(:,:,:))
          END IF
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,    &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         ADM(ng)%Vid(idSbry(isVvel)),             &
     &                         ADM(ng)%Rindex, v3dvar,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % ad_W_sol,                         &
     &                    Wr3d)
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                      &
     &                     ADM(ng)%Vid(idOvel),                         &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
#  ifdef WEAK_CONSTRAINT
          IF (WRTforce(ng)) THEN
            scale=1.0_dp
            gtype=gfactor*r3dvar
            status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                  &
     &                         ADM(ng)%Tid(itrc),                       &
     &                         ADM(ng)%Rindex, gtype,                   &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         OCEAN(ng) % f_tG(:,:,:,kfout,itrc))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            ADM(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
#  endif
            scale=1.0_dp
            gtype=gfactor*r3dvar
#  ifdef AD_OUTPUT_STATE
            IF (LwrtState3d(ng)) THEN
#  endif
              status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                &
     &                           ADM(ng)%Tid(itrc),                     &
     &                           ADM(ng)%Rindex, gtype,                 &
     &                           LBi, UBi, LBj, UBj, 1, N(ng), scale,   &
#   ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#   endif
     &                           OCEAN(ng) % ad_t(:,:,:,nout,itrc))
#  ifdef AD_OUTPUT_STATE
            ELSE
              status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid,                &
     &                           ADM(ng)%Tid(itrc),                     &
     &                           ADM(ng)%Rindex, gtype,                 &
     &                           LBi, UBi, LBj, UBj, 1, N(ng), scale,   &
#   ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#   endif
     &                           OCEAN(ng) % ad_t_sol(:,:,:,itrc))
            END IF
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            ADM(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
#  ifdef WEAK_CONSTRAINT
          END IF
#  endif
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          status=nf_fwrite3d_bry(ng, iADM, ADM(ng)%name, ADM(ng)%ncid,  &
     &                           Vname(1,idSbry(isTvar(itrc))),         &
     &                           ADM(ng)%Vid(idSbry(isTvar(itrc))),     &
     &                           ADM(ng)%Rindex, r3dvar,                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,   &
     &                                                 Lbout(ng),itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          ADM(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idDano), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % ad_rho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVvis), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % ad_Akv,                         &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), ADM(ng)%Rindex
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
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idTdif), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % ad_Akt(:,:,:,itemp),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idSdif), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % ad_Akt(:,:,:,isalt),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifndef ADJUST_STFLUX
!
!  Write out net surface active tracer fluxes.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTsur(itrc),ng)) THEN
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
       defined OPT_OBSERVATIONS
          IF (itrc.eq.itemp) THEN
!!          scale=rho0*Cp
            scale=1.0_dp/(rho0*Cp)
          ELSE
            scale=1.0_dp
          END IF
#   else
          scale=1.0_dp
#   endif
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid,                    &
     &                       ADM(ng)%Vid(idTsur(itrc)),                 &
     &                       ADM(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
# endif
# ifndef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
!!      scale=rho0
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idUsms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % ad_sustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), ADM(ng)%Rindex
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
!!      scale=rho0
#  if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#  else
        scale=1.0_dp
#  endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVsms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % ad_svstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
!!      scale=-rho0
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idUbms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % ad_bustr_sol)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
!!      scale=-rho0
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iADM, ADM(ng)%ncid, ADM(ng)%Vid(idVbms), &
     &                     ADM(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % ad_bvstr_sol)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize adjoint history NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iADM, ADM(ng)%name, ADM(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, ADM(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, nout, ADM(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, ADM(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,20) kout, ADM(ng)%Rindex
#  endif
# endif
!
  10  FORMAT (/,' AD_WRT_HIS_NF90 - error while writing variable: ',a,  &
     &        /,19x,'into adjoint NetCDF file for time record: ',i0)
  20  FORMAT (2x,'AD_WRT_HIS_NF90  - wrote adjoint', t40,               &
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
      END SUBROUTINE ad_wrt_his_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE ad_wrt_his_pio (ng, tile,                              &
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
      integer, intent(in) :: ng, tile
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      integer :: Fcount, i, ifield, j, status
      integer :: kout
#  ifdef WEAK_CONSTRAINT
      integer :: kfout
#  endif
#  ifdef SOLVE3D
      integer :: itrc, k, nout
#  endif
!
      real(dp) :: scale
      real(r8) :: Tval(1)
#ifdef SOLVE3D
!
      real(r8), allocatable :: Wr3d(:,:,:)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ad_wrt_his_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out adjoint fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Determine time index to write. The "nout" index is updated to the
!  version of "ad_main3d" that updates the "iic" counter at the bottom.
!  Therefore, we need to change the conditional "iic(ng).ne.ntend(ng)"
!  to "iic(ng).gt.ntend(ng)" to get identical solutions.
!
#  ifdef SOLVE3D
      kout=kstp(ng)
#  else
      kout=kstp(ng)
#  endif
#  if defined WEAK_CONSTRAINT
      kfout=2
#  endif
#  ifdef SOLVE3D
      IF (iic(ng).gt.ntend(ng)) THEN
        nout=nnew(ng)
#  ifdef AD_OUTPUT_STATE
        LwrtState3d(ng)=.FALSE.
#  endif
      ELSE
#  ifdef AD_OUTPUT_STATE
        LwrtState3d(ng)=.TRUE.
#  endif
        nout=nstp(ng)
      END IF
#  endif
!
!  Set time record index.
!
      ADM(ng)%Rindex=ADM(ng)%Rindex+1
      Fcount=ADM(ng)%load
      ADM(ng)%Nrec(Fcount)=ADM(ng)%Nrec(Fcount)+1
!
!  If requested, set time index to recycle time records in the adjoint
!  file.
!
      IF (LcycleADJ(ng)) THEN
        ADM(ng)%Rindex=MOD(ADM(ng)%Rindex-1,2)+1
      END IF
!
!  Write out model time (s).
!
      IF (LwrtTime(ng)) THEN
        IF (LwrtPER(ng)) THEN
          Tval(1)=REAL(ADM(ng)%Rindex,r8)*day2sec
        ELSE
#  if defined WEAK_CONSTRAINT && !defined WEAK_NOINTERP
          Tval(1)=ForceTime(ng)
#  else
          Tval(1)=time(ng)
#  endif
        END IF
        CALL pio_netcdf_put_fvar (ng, iADM, ADM(ng)%name,               &
     &                            TRIM(Vname(1,idtime)), tval,          &
     &                            (/ADM(ng)%Rindex/), (/1/),            &
     &                            pioFile = ADM(ng)%pioFile,            &
     &                            pioVar = ADM(ng)%pioVar(idtime)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp                          ! m2/s2
      IF (ADM(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                     &
     &                   ADM(ng)%pioVar(idUsms),                        &
     &                   ADM(ng)%Rindex, ioDesc,                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Lfout(ng)))
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
      IF (ADM(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                     &
     &                   ADM(ng)%pioVar(idVsms),                        &
     &                   ADM(ng)%Rindex, ioDesc,                        &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Lfout(ng)))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net heat flux. Notice that different tracer fluxes
!  are written at their own fixed time-dimension (of size Nfrec) to
!  allow 4DVAR adjustments at other times in addition to initial time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp                      ! kinematic flux units
          IF (ADM(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idTsur(itrc)),              &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng)% ad_tflux(:,:,:,Lfout(ng),itrc))
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
#  endif
!
!  Write out bathymetry.
!
      IF (Hout(idbath,ng)) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idbath)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idbath),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     GRID(ng)% ad_h,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      scale=1.0_dp
      IF (ADM(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      IF (Hout(idFsur,ng)) THEN
#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idFsur),                    &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       OCEAN(ng)% f_zetaG(:,:,kfout))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idFsur),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
#  ifdef WET_DRY
     &                         OCEAN(ng)% ad_zeta(:,:,kout),            &
     &                         SetFillVal = .FALSE.)
#  else
     &                         OCEAN(ng)% ad_zeta(:,:,kout))
#  endif
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idFsur),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
#  ifdef WET_DRY
     &                         OCEAN(ng)% ad_zeta_sol,                  &
     &                         SetFillVal = .FALSE.)
#  else
     &                         OCEAN(ng)% ad_zeta_sol)
#  endif
          ENDIF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsur)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idSbry(isFsur))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iADM, ADM(ng)%name,                  &
     &                         ADM(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isFsur)),                 &
     &                         ADM(ng)%pioVar(idSbry(isFsur)),          &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        ADM(ng)%Rindex
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
      IF (ADM(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      IF (Hout(idUbar,ng)) THEN
#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idUbar),                    &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#   endif
     &                       OCEAN(ng) % f_ubarG(:,:,kfout))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idUbar),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#  endif
     &                         OCEAN(ng) % ad_ubar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idUbar),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#  endif
     &                         OCEAN(ng) % ad_ubar_sol)
          END IF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idSbry(isUbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iADM, ADM(ng)%name,                  &
     &                         ADM(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUbar)),                 &
     &                         ADM(ng)%pioVar(idSbry(isUbar)),          &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        ADM(ng)%Rindex
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
      IF (ADM(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      IF (Hout(idVbar,ng)) THEN

#  ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idVbar),                    &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       OCEAN(ng) % f_vbarG(:,:,kfout))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#  endif
          IF (LwrtState2d(ng)) THEN
            status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idVbar),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#  endif
     &                         OCEAN(ng) % ad_vbar(:,:,kout))
          ELSE
            status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idVbar),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, scale,               &
#  ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#  endif
     &                         OCEAN(ng) % ad_vbar_sol)
          END IF
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbar)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#  ifdef WEAK_CONSTRAINT
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idSbry(isVbar))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fwrite2d_bry(ng, iADM, ADM(ng)%name,                  &
     &                         ADM(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVbar)),                 &
     &                         ADM(ng)%pioVar(idSbry(isVbar)),          &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, Nbrec(ng), scale,            &
     &                         BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,    &
     &                                                    Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      scale=1.0_dp
      IF (ADM(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      IF (Hout(idUvel,ng)) THEN
#   ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idUvel),                    &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#    endif
     &                       OCEAN(ng) % f_uG(:,:,:,kfout))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#   endif
#   ifdef AD_OUTPUT_STATE
          IF (LwrtState3d(ng)) THEN
#   endif
            status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idUvel),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#     ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#     endif
     &                         OCEAN(ng) % ad_u(:,:,:,nout))
#   ifdef AD_OUTPUT_STATE
          ELSE
            status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idUvel),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#     ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#     endif
     &                         OCEAN(ng) % ad_u_sol(:,:,:))
          ENDIF
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#   ifdef WEAK_CONSTRAINT
        END IF
#   endif
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idSbry(isUvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iADM, ADM(ng)%name,                  &
     &                         ADM(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isUvel)),                 &
     &                         ADM(ng)%pioVar(idSbry(isUvel)),          &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        ADM(ng)%Rindex
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
      IF (ADM(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      IF (Hout(idVvel,ng)) THEN
#   ifdef WEAK_CONSTRAINT
        IF (WRTforce(ng)) THEN
          status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idVvel),                    &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#    endif
     &                       OCEAN(ng) % f_vG(:,:,:,kfout))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
#   endif
#   ifdef AD_OUTPUT_STATE
          IF (LwrtState3d(ng)) THEN
#   endif
            status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idVvel),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#    ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#    endif
     &                         OCEAN(ng) % ad_v(:,:,:,nout))
#   ifdef AD_OUTPUT_STATE
          ELSE
            status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioVar(idVvel),                  &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#    ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#    endif
     &                         OCEAN(ng) % ad_v_sol(:,:,:))
          END IF
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVvel)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
#   ifdef WEAK_CONSTRAINT
        END IF
#   endif
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idSbry(isVvel))%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fwrite3d_bry(ng, iADM, ADM(ng)%name,                  &
     &                         ADM(ng)%pioFile,                         &
     &                         Vname(1,idSbry(isVvel)),                 &
     &                         ADM(ng)%pioVar(idSbry(isVvel)),          &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng), scale,  &
     &                         BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng)))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % ad_W_sol,                         &
     &                    Wr3d)
!
        IF (ADM(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idOvel),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (ADM(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        IF (Hout(idTvar(itrc),ng)) THEN

#   ifdef WEAK_CONSTRAINT
          IF (WRTforce(ng)) THEN
            status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,               &
     &                         ADM(ng)%pioTrc(itrc),                    &
     &                         ADM(ng)%Rindex, ioDesc,                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#    ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#    endif
     &                         OCEAN(ng) % f_tG(:,:,:,kfout,itrc))
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            ADM(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
#   endif
#   ifdef AD_OUTPUT_STATE
            IF (LwrtState3d(ng)) THEN
#   endif
              status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,             &
     &                           ADM(ng)%pioTrc(itrc),                  &
     &                           ADM(ng)%Rindex, ioDesc,                &
     &                           LBi, UBi, LBj, UBj, 1, N(ng), scale,   &
#    ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#    endif
     &                           OCEAN(ng) % ad_t(:,:,:,nout,itrc))
#   ifdef AD_OUTPUT_STATE
            ELSE
              status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,             &
     &                           ADM(ng)%pioTrc(itrc),                  &
     &                           ADM(ng)%Rindex, ioDesc,                &
     &                           LBi, UBi, LBj, UBj, 1, N(ng), scale,   &
#    ifdef MASKING
     &                           GRID(ng) % rmask,                      &
#    endif
     &                           OCEAN(ng) % ad_t_sol(:,:,:,itrc))
            END IF
#   endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),          &
     &                            ADM(ng)%Rindex
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
#   ifdef WEAK_CONSTRAINT
          END IF
#   endif
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
          IF (ADM(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fwrite3d_bry(ng, iADM, ADM(ng)%name,                &
     &                           ADM(ng)%pioFile,                       &
     &                           Vname(1,ifield),                       &
     &                           ADM(ng)%pioVar(ifield),                &
     &                           ADM(ng)%Rindex, ioDesc,                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale,                                 &
     &                           BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,   &
     &                                                 Lbout(ng),itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,ifield)), ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
!
!  Write out density anomaly.
!
      IF (Hout(idDano,ng)) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idDano),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % ad_rho)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idVvis),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % ad_Akv,                         &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), ADM(ng)%Rindex
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
        IF (ADM(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idTdif),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % ad_Akt(:,:,:,itemp),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), ADM(ng)%Rindex
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
        IF (ADM(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
!
        status=nf_fwrite3d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idSdif),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % ad_Akt(:,:,:,isalt),            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   ifndef ADJUST_STFLUX
!
!  Write out net surface active tracer fluxes.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTsur(itrc),ng)) THEN
#    if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
       defined OPT_OBSERVATIONS
          IF (itrc.eq.itemp) THEN
!!          scale=rho0*Cp
            scale=1.0_dp/(rho0*Cp)
          ELSE
            scale=1.0_dp
          END IF
#    else
          scale=1.0_dp
#    endif
          IF (ADM(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
!
          status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                 &
     &                       ADM(ng)%pioVar(idTsur(itrc)),              &
     &                       ADM(ng)%Rindex, ioDesc,                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       FORCES(ng) % ad_stflx(:,:,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          ADM(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
#  endif
#  ifndef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
!!      scale=rho0
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
        IF (ADM(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idUsms),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     FORCES(ng) % ad_sustr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), ADM(ng)%Rindex
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
!!      scale=rho0
#   if defined AD_SENSITIVITY   || defined I4DVAR_ANA_SENSITIVITY || \
      defined OPT_OBSERVATIONS
        scale=1.0_dp/rho0
#   else
        scale=1.0_dp
#   endif
        IF (ADM(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idVsms),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     FORCES(ng) % ad_svstr)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
!!      scale=-rho0
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idUbms),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % ad_bustr_sol)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
!!      scale=-rho0
        scale=1.0_dp
        IF (ADM(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
!
        status=nf_fwrite2d(ng, iADM, ADM(ng)%pioFile,                   &
     &                     ADM(ng)%pioVar(idVbms),                      &
     &                     ADM(ng)%Rindex, ioDesc,                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % ad_bvstr_sol)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), ADM(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize adjoint history NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iADM, ADM(ng)%name, ADM(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, nout, ADM(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, nout, ADM(ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,20) kout, ADM(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,20) kout, ADM(ng)%Rindex
#   endif
#  endif
!
  10  FORMAT (/,' AD_WRT_HIS_PIO - error while writing variable: ',a,   &
     &        /,18x,'into adjoint NetCDF file for time record: ',i0)
#  ifdef SOLVE3D
  20  FORMAT (2x,'AD_WRT_HIS_PIO   - wrote adjoint', t40,               &
#   ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
#   endif
#  else
  20  FORMAT (2x,'AD_WRT_HIS_PIO   - wrote adjoint', t40,               &
#   ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
#   else
     &        'fields (Index=',i1,')   in record = ',i0)
#   endif
#  endif
!
      RETURN
      END SUBROUTINE ad_wrt_his_pio
# endif
#endif
      END MODULE ad_wrt_his_mod
