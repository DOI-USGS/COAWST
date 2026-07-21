#include "cppdefs.h"
      MODULE wrt_extract_mod

#ifdef GRID_EXTRACT
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes requested model fields into the EXTRACT output   !
!  file using either the standard NetCDF library or the Parallel-IO    !
!  (PIO) library.                                                      !
!                                                                      !
!  If ExtractFlag > 1:   Decimation                                    !
!                                                                      !
!  It decimates the field solution at the prescribed integer factor,   !
!  ExtractFlag . For example, if ExtractFlag=2 (recommended), the      !
!  output fieldsare written at every other point, resulting in coarser !
!  data resolution. This strategy is advantageous in mixed resolution, !
!  split 4D-Var applications where the outer loop background (prior)   !
!  trajectory may be computed at a higher resolution than in the inner !
!  loop minimization to accelerate the calculations. For decimation to !
!  work, the number of parent grid RHO-points (0: Lm+1, 0:Mm+1) must   !
!  be multiples of dec_fator. That is,                                 !
!                                                                      !
!     MOD(Lm+1, ExtractFlag) = 0                                       !
!     MOD(Mm+1, ExtractFlag) = 0                                       !
!                                                                      !
!  If ExtractFlag = 1:   Interpolation                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BBL_MODEL
      USE mod_bbl
# endif
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_grid
      USE mod_extract
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
# endif
      USE mod_stepping
!
# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
      USE bbl_output_mod,      ONLY : bbl_wrt_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,      ONLY : bbl_wrt_pio
#  endif
# endif
      USE extract_field_mod,   ONLY : extract_field
# if defined ICE_MODEL && defined SOLVE3D
      USE ice_output_mod,      ONLY : ice_wrt_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,      ONLY : ice_wrt_pio
#  endif
# endif
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
# if defined SEDIMENT && defined SOLVE3D
      USE sediment_output_mod, ONLY : sediment_wrt_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod, ONLY : sediment_wrt_pio
#  endif
# endif
      USE strings_mod,         ONLY : FoundError
      USE uv_rotate_mod,       ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod,       ONLY : uv_rotate3d
# endif
# if defined WEC_VF && defined SOLVE3D
      USE wec_output_mod,      ONLY : wec_wrt_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,      ONLY : wec_wrt_pio
#  endif
# endif
!
      implicit none
!
      PUBLIC  :: wrt_extract
      PRIVATE :: wrt_extract_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_extract_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_extract (ng, tile)
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
      integer :: iLB, iUB, jLB, jUB
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out decimate fields according to IO type.
!-----------------------------------------------------------------------
!
# ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij                    !  parent boundary bounds
      UBij=BOUNDS(ng)%UBij
# endif
      LBi=BOUNDS(ng)%LBi(tile)                ! parent array bounds
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      iLB=xtr_BOUNDS(ng)%LBi(tile)            ! extract array bounds
      iUB=xtr_BOUNDS(ng)%UBi(tile)
      jLB=xtr_BOUNDS(ng)%LBj(tile)
      jUB=xtr_BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (XTR(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_extract_nf90 (ng, iNLM, tile,                        &
# ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
# endif
     &                           LBi, UBi, LBj, UBj,                    &
     &                           iLB, iUB, jLB, jUB)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_extract_pio (ng, iNLM, tile,                         &
#  ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
#  endif
     &                          LBi, UBi, LBj, UBj,                     &
     &                          iLB, iUB, jLB, jUB)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) XTR(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_EXTRACT - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_extract
!
!***********************************************************************
      SUBROUTINE wrt_extract_nf90 (ng, model, tile,                     &
# ifdef ADJUST_BOUNDARY
     &                             LBij, UBij,                          &
# endif
     &                             LBi, UBi, LBj, UBj,                  &
     &                             iLB, iUB, jLB, jUB)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: iLB, iUB, jLB, jUB
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, ifield, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k
# endif
!
      real(dp) :: scale
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
# ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_extract_nf90"

#  include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out decimate fields.
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
      XTR(ng)%Rindex=XTR(ng)%Rindex+1
      Fcount=XTR(ng)%load
      XTR(ng)%Nrec(Fcount)=XTR(ng)%Nrec(Fcount)+1
!
!  Report.
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, NOUT, XTR(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) KOUT, NOUT, XTR(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, XTR(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) KOUT, XTR(ng)%Rindex
#  endif
# endif
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, XTR(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/XTR(ng)%Rindex/), (/1/),                  &
     &                      ncid = XTR(ng)%ncid,                        &
     &                      varid = XTR(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef WET_DRY
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_dp
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idPwet,               &
     &                   XTR(ng)%Vid(idPwet),                           &
     &                   XTR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   EXTRACT(ng) % pmask,                           &
#  endif
     &                   EXTRACT(ng) % pmask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idPwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idRwet,               &
     &                   XTR(ng)%Vid(idRwet),                           &
     &                   XTR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   EXTRACT(ng) % rmask,                           &
#  endif
     &                   EXTRACT(ng) % rmask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idRwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUwet,               &
     &                   XTR(ng)%Vid(idUwet),                           &
     &                   XTR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   EXTRACT(ng) % umask,                           &
#  endif
     &                   EXTRACT(ng) % umask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVwet,               &
     &                   XTR(ng)%Vid(idVwet),                           &
     &                   XTR(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   EXTRACT(ng) % vmask,                           &
#  endif
     &                   EXTRACT(ng) % vmask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idpthR,             &
     &                     XTR(ng)%Vid(idpthR),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     EXTRACT(ng) % rmask,                         &
#  endif
     &                     EXTRACT(ng) % z_r,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              EXTRACT(ng)%z_v(i,j,k)=0.5_r8*                            &
     &                               (EXTRACT(ng)%z_r(i-1,j,k)+         &
     &                                EXTRACT(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idpthU,             &
     &                     XTR(ng)%Vid(idpthU),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     EXTRACT(ng) % umask,                         &
#  endif
     &                     EXTRACT(ng) % z_v,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              EXTRACT(ng)%z_v(i,j,k)=0.5_r8*(EXTRACT(ng)%z_r(i,j-1,k)+  &
     &                                       EXTRACT(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idpthV,             &
     &                     XTR(ng)%Vid(idpthV),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     EXTRACT(ng) % vmask,                         &
#  endif
     &                     EXTRACT(ng) % z_v,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idpthW,             &
     &                     XTR(ng)%Vid(idpthW),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     EXTRACT(ng) % rmask,                         &
#  endif
     &                     EXTRACT(ng) % z_w,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idFsur,             &
     &                     XTR(ng)%Vid(idFsur),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# if defined FORWARD_WRITE && defined FORWARD_RHS
!
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idRzet,             &
     &                     XTR(ng)%Vid(idRzet),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % rzeta(:,:,KOUT),                 &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRzet)), XTR(ng)%Rindex
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
        status=nf_fwrite2d_bry (ng, model, XTR(ng)%name, XTR(ng)%ncid,  &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          XTR(ng)%Vid(idSbry(isFsur)),            &
     &                          XTR(ng)%Rindex, r2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)),     &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        XTR(ng)%Rindex
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
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUbar,             &
     &                     XTR(ng)%Vid(idUbar),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     OCEAN(ng) % ubar(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idRu2d,             &
     &                     XTR(ng)%Vid(idRu2d),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % rubar(:,:,KOUT),                 &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu2d)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idRuct,             &
     &                     XTR(ng)%Vid(idRuct),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     COUPLING(ng) % rufrc,                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRuct)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        Ur2d=ExtractFlag(ng)*COUPLING(ng)%DU_avg1
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUfx1,             &
     &                     XTR(ng)%Vid(idUfx1),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx1)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        Ur2d=ExtractFlag(ng)*COUPLING(ng)%DU_avg2
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUfx2,             &
     &                     XTR(ng)%Vid(idUfx2),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx2)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, XTR(ng)%name, XTR(ng)%ncid,  &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          XTR(ng)%Vid(idSbry(isUbar)),            &
     &                          XTR(ng)%Rindex, u2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)),     &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        XTR(ng)%Rindex
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
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVbar,             &
     &                     XTR(ng)%Vid(idVbar),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     OCEAN(ng) % vbar(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idRv2d,             &
     &                     XTR(ng)%Vid(idRv2d),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % rvbar(:,:,KOUT),                 &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv2d)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idRvct,             &
     &                     XTR(ng)%Vid(idRvct),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     COUPLING(ng) % rvfrc,                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRvct)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
!
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        Vr2d=ExtractFlag(ng)*COUPLING(ng)%DV_avg1
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVfx1,             &
     &                     XTR(ng)%Vid(idVfx1),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx1)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        Vr2d=ExtractFlag(ng)*COUPLING(ng)%DV_avg2
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVfx2,             &
     &                     XTR(ng)%Vid(idVfx2),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx2)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Vr2d)
#  endif
# endif
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, XTR(ng)%name, XTR(ng)%ncid,  &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          XTR(ng)%Vid(idSbry(isVbar)),            &
     &                          XTR(ng)%Rindex, v2dvar,                 &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)),     &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idu2dE,             &
     &                     XTR(ng)%Vid(idu2dE),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idv2dN,             &
     &                     XTR(ng)%Vid(idv2dN),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idUvel,             &
     &                     XTR(ng)%Vid(idUvel),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % u(:,:,:,NOUT),                   &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idRu3d,             &
     &                     XTR(ng)%Vid(idRu3d),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % ru(:,:,:,NOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu3d)), XTR(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, model, XTR(ng)%name, XTR(ng)%ncid,  &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          XTR(ng)%Vid(idSbry(isUvel)),            &
     &                          XTR(ng)%Rindex, u3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)),        &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        XTR(ng)%Rindex
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
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idVvel,             &
     &                     XTR(ng)%Vid(idVvel),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % v(:,:,:,NOUT),                   &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idRv3d,             &
     &                     XTR(ng)%Vid(idRv3d),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % rv(:,:,:,NOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv3d)), XTR(ng)%Rindex
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
        status=nf_fwrite3d_bry (ng, model, XTR(ng)%name, XTR(ng)%ncid,  &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          XTR(ng)%Vid(idSbry(isVvel)),            &
     &                          XTR(ng)%Rindex, v3dvar,                 &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)),        &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#  ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#  endif
     &                    OCEAN(ng) % u(:,:,:,NOUT),                    &
     &                    OCEAN(ng) % v(:,:,:,NOUT),                    &
     &                    Ur3d, Vr3d)
!
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idu3dE,             &
     &                     XTR(ng)%Vid(idu3dE),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Ur3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idv3dN,             &
     &                     XTR(ng)%Vid(idv3dN),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Vr3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
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
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idOvel,             &
     &                     XTR(ng)%Vid(idOvel),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Wr3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

#  ifdef OMEGA_IMPLICIT
!
!  Write out S-coordinate implicit vertical "omega" momentum component.
!
      IF (Hout(idOvil,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % Wi,                               &
     &                    Wr3d)
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idOvil,             &
     &                     XTR(ng)%Vid(idOvil),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Wr3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvil)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
#  endif
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idWvel,             &
     &                     XTR(ng)%Vid(idWvel),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % wvel,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idTvar(itrc),     &
     &                       XTR(ng)%Tid(itrc),                         &
     &                       XTR(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc),            &
     &                       SetFillVal = .FALSE.,                      &
     &                       ExtractField = ExtractFlag(ng))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          XTR(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          status=nf_fwrite3d_bry (ng, model, XTR(ng)%name, XTR(ng)%ncid,&
     &                            Vname(1,ifield),                      &
     &                            XTR(ng)%Vid(ifield),                  &
     &                            XTR(ng)%Rindex, r3dvar,               &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng),itrc), &
     &                          ExtractField = ExtractFlag(ng))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), XTR(ng)%Rindex
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
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idDano,             &
     &                     XTR(ng)%Vid(idDano),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % rho,                             &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idHsbl,             &
     &                     XTR(ng)%Vid(idHsbl),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % hsbl,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHsbl)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth bottom boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idHbbl,             &
     &                     XTR(ng)%Vid(idHbbl),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % hbbl,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHbbl)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        IF (Hout(idGhat(i),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*w3dvar
          status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idGhat(i),        &
     &                       XTR(ng)%Vid(idGhat(i)),                    &
     &                       XTR(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       MIXING(ng) % ghats(:,:,:,i),               &
     &                       SetFillVal = .FALSE.,                      &
     &                       ExtractField = ExtractFlag(ng))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idGhat(i))), XTR(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idVvis,             &
     &                     XTR(ng)%Vid(idVvis),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), XTR(ng)%Rindex
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
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idTdif,             &
     &                     XTR(ng)%Vid(idTdif),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), XTR(ng)%Rindex
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
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idSdif,             &
     &                     XTR(ng)%Vid(idSdif),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idMtke,             &
     &                     XTR(ng)%Vid(idMtke),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#   ifdef FORWARD_WRITE
!
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idVmKK,             &
     &                     XTR(ng)%Vid(idVmKK),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % Akk,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKK)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idMtls,             &
     &                     XTR(ng)%Vid(idMtls),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#   ifdef FORWARD_WRITE
!
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idVmLS,             &
     &                     XTR(ng)%Vid(idVmLS),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % Lscale,                         &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmLS)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
#   if defined FORWARD_WRITE && defined GLS_MIXING
!
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, XTR(ng)%ncid, idVmKP,             &
     &                     XTR(ng)%Vid(idVmKP),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % Akp,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKP)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idPair,             &
     &                     XTR(ng)%Vid(idPair),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Pair,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES
!
!  Write out surface air temperature.
!
      IF (Hout(idTair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idTair,             &
     &                     XTR(ng)%Vid(idTair),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Tair,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUair,             &
     &                     XTR(ng)%Vid(idUair),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Uwind,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVair,             &
     &                     XTR(ng)%Vid(idVair),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Vwind,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUaiE,ng).and.Hout(idVaiN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#   ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#   endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUaiE,             &
     &                     XTR(ng)%Vid(idUaiE),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVaiN,             &
     &                     XTR(ng)%Vid(idVaiN),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
       deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#  endif
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
#  ifdef SO_SEMI
            scale=1.0_dp
#  else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
#  endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idTsur(itrc),     &
     &                       XTR(ng)%Vid(idTsur(itrc)),                 &
     &                       XTR(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % stflx(:,:,itrc),              &
     &                       SetFillVal = .FALSE.,                      &
     &                       ExtractField = ExtractFlag(ng))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          XTR(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#  if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idLhea,             &
     &                     XTR(ng)%Vid(idLhea),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % lhflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idShea,             &
     &                     XTR(ng)%Vid(idShea),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % shflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idLrad,             &
     &                     XTR(ng)%Vid(idLrad),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % lrflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif

#  ifdef BULK_FLUXES
#   ifdef EMINUSP
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idevap,             &
     &                     XTR(ng)%Vid(idevap),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % evap,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idrain,             &
     &                     XTR(ng)%Vid(idrain),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % rain,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idEmPf,             &
     &                     XTR(ng)%Vid(idEmPf),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % stflux(:,:,isalt),              &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idSrad,             &
     &                     XTR(ng)%Vid(idSrad),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % srflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
# ifdef SO_SEMI
        scale=1.0_dp
# else
        scale=rho0                          ! m2/s2 to Pa
# endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUsms,             &
     &                     XTR(ng)%Vid(idUsms),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % sustr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), XTR(ng)%Rindex
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
# ifdef SO_SEMI
        scale=1.0_dp
# else
        scale=rho0
# endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVsms,             &
     &                     XTR(ng)%Vid(idVsms),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % svstr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idUbms,             &
     &                     XTR(ng)%Vid(idUbms),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % bustr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), XTR(ng)%Rindex
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
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, XTR(ng)%ncid, idVbms,             &
     &                     XTR(ng)%Vid(idVbms),                         &
     &                     XTR(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % bvstr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_nf90 (ng, model, tile,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the  Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_nf90 (ng, model, tile,                               &
    &                    LBi, UBi, LBj, UBj,                            &
    &                    Hout, XTR)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !(defined OUT_NETCDF4 && defined DELAYED_SYNC_NF90)
!
!-----------------------------------------------------------------------
!  Synchronize decimate NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, XTR(ng)%name, XTR(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
  10  FORMAT (2x,'WRT_EXTRACT_NF90 - writing decimate', t42,            &
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
  20  FORMAT (/,' WRT_EXTRACT_NF90 - error while writing variable: ',a, &
     &        /,20x,'into decimate NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_extract_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_extract_pio (ng, model, tile,                      &
#  ifdef ADJUST_BOUNDARY
     &                            LBij, UBij,                           &
#  endif
     &                            LBi, UBi, LBj, UBj,                   &
     &                            iLB, iUB, jLB, jUB)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: iLB, iUB, jLB, jUB
!
!  Local variable declarations.
!
      integer :: Fcount, ifield, status
#  ifdef SOLVE3D
      integer :: i, itrc, j, k
#  endif
!
      real(dp) :: scale
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
#  ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_extract_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesX

#   include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out decimate fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      XTR(ng)%Rindex=XTR(ng)%Rindex+1
      Fcount=XTR(ng)%load
      XTR(ng)%Nrec(Fcount)=XTR(ng)%Nrec(Fcount)+1
!
!  Report.
!
#  ifdef SOLVE3D
#   ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, NOUT, XTR(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,10) KOUT, NOUT, XTR(ng)%Rindex
#   endif
#  else
#   ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, XTR(ng)%Rindex, ng
#   else
      IF (Master) WRITE (stdout,10) KOUT, XTR(ng)%Rindex
#   endif
#  endif
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, model, XTR(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/XTR(ng)%Rindex/), (/1/),              &
     &                           pioFile = XTR(ng)%pioFile,             &
     &                           pioVar = XTR(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef WET_DRY
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_dp
      IF (XTR(ng)%pioVar(idPwet)%dkind.eq.PIO_double) THEN
        ioDesX => ioDesX_dp_p2dvar(ng)
      ELSE
        ioDesX => ioDesX_sp_p2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idPwet,            &
     &                   XTR(ng)%pioVar(idPwet),                        &
     &                   XTR(ng)%Rindex,                                &
     &                   ioDesX,                                        &
     &                   iLB, iUB, jLB, jUB, scale,                     &
#   ifdef MASKING
     &                   EXTRACT(ng) % pmask,                           &
#   endif
     &                   EXTRACT(ng) % pmask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idPwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_dp
      IF (XTR(ng)%pioVar(idRwet)%dkind.eq.PIO_double) THEN
        ioDesX => ioDesX_dp_r2dvar(ng)
      ELSE
        ioDesX => ioDesX_sp_r2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idRwet,            &
     &                   XTR(ng)%pioVar(idRwet),                        &
     &                   XTR(ng)%Rindex,                                &
     &                   ioDesX,                                        &
     &                   iLB, iUB, jLB, jUB, scale,                     &
#   ifdef MASKING
     &                   EXTRACT(ng) % rmask,                           &
#   endif
     &                   EXTRACT(ng) % rmask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idRwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_dp
      IF (XTR(ng)%pioVar(idUwet)%dkind.eq.PIO_double) THEN
        ioDesX => ioDesX_dp_u2dvar(ng)
      ELSE
        ioDesX => ioDesX_sp_u2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUwet,            &
     &                   XTR(ng)%pioVar(idUwet),                        &
     &                   XTR(ng)%Rindex,                                &
     &                   ioDesX,                                        &
     &                   iLB, iUB, jLB, jUB, scale,                     &
#   ifdef MASKING
     &                   EXTRACT(ng) % umask,                           &
#   endif
     &                   EXTRACT(ng) % umask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_dp
      IF (XTR(ng)%pioVar(idVwet)%dkind.eq.PIO_double) THEN
        ioDesX => ioDesX_dp_v2dvar(ng)
      ELSE
        ioDesX => ioDesX_sp_v2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVwet,            &
     &                   XTR(ng)%pioVar(idVwet),                        &
     &                   XTR(ng)%Rindex,                                &
     &                   ioDesX,                                        &
     &                   iLB, iUB, jLB, jUB, scale,                     &
#   ifdef MASKING
     &                   EXTRACT(ng) % vmask,                           &
#   endif
     &                   EXTRACT(ng) % vmask_wet,                       &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVwet)), XTR(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef SOLVE3D
!
!  Write time-varying depths of RHO-points.
!
      IF (Hout(idpthR,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idpthR)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idpthR,          &
     &                     XTR(ng)%pioVar(idpthR),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     iLB, iUB, jLB, jUB, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     EXTRACT(ng) % rmask,                         &
#   endif
     &                     EXTRACT(ng) % z_r,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Hout(idpthU,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              EXTRACT(ng)%z_v(i,j,k)=0.5_r8*                            &
     &                               (EXTRACT(ng)%z_r(i-1,j,k)+         &
     &                                EXTRACT(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        IF (XTR(ng)%pioVar(idpthU)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idpthU,          &
     &                     XTR(ng)%pioVar(idpthU),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     iLB, iUB, jLB, jUB, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     EXTRACT(ng) % umask,                         &
#   endif
     &                     EXTRACT(ng) % z_v,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Hout(idpthV,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              EXTRACT(ng)%z_v(i,j,k)=0.5_r8*                            &
     &                               (EXTRACT(ng)%z_r(i,j-1,k)+         &
     &                                EXTRACT(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        IF (XTR(ng)%pioVar(idpthV)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idpthV,          &
     &                     XTR(ng)%pioVar(idpthV),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     iLB, iUB, jLB, jUB, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     EXTRACT(ng) % vmask,                         &
#   endif
     &                     EXTRACT(ng) % z_v,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Hout(idpthW,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idpthW)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idpthW,          &
     &                     XTR(ng)%pioVar(idpthW),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     iLB, iUB, jLB, jUB, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     EXTRACT(ng) % rmask,                         &
#   endif
     &                     EXTRACT(ng) % z_w,                           &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idFsur,          &
     &                     XTR(ng)%pioVar(idFsur),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRzet)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idRzet,          &
     &                     XTR(ng)%pioVar(idRzet),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % rzeta(:,:,KOUT),                 &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRzet)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, XTR(ng)%name,                &
     &                          XTR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          XTR(ng)%pioVar(idSbry(isFsur)),         &
     &                          XTR(ng)%Rindex,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)),     &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        XTR(ng)%Rindex
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
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUbar,          &
     &                     XTR(ng)%pioVar(idUbar),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % ubar(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  ifdef FORWARD_WRITE
#   ifdef FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRu2d)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idRu2d,          &
     &                     XTR(ng)%pioVar(idRu2d),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     OCEAN(ng) % rubar(:,:,KOUT),                 &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu2d)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
#   ifdef SOLVE3D
#    ifdef FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRuct)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idRuct,          &
     &                     XTR(ng)%pioVar(idRuct),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#     endif
     &                     COUPLING(ng) % rufrc,                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRuct)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
!
        IF (XTR(ng)%pioVar(idUfx1)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        Ur2d=ExtractFlag(ng)*COUPLING(ng)%DU_avg1
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUfx1,          &
     &                     XTR(ng)%pioVar(idUfx1),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx1)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (XTR(ng)%pioVar(idUfx2)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        Ur2d=ExtractFlag(ng)*COUPLING(ng)%DU_avg2
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUfx2,          &
     &                     XTR(ng)%pioVar(idUfx2),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     COUPLING(ng) % DU_avg2,                      &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUfx2)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
#   endif
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, XTR(ng)%name,                &
     &                          XTR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          XTR(ng)%pioVar(idSbry(isUbar)),         &
     &                          XTR(ng)%Rindex,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)),     &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        XTR(ng)%Rindex
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
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVbar,          &
     &                     XTR(ng)%pioVar(idVbar),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % vbar(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#  ifdef FORWARD_WRITE
#   ifdef FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRv2d)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idRv2d,          &
     &                     XTR(ng)%pioVar(idRv2d),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     OCEAN(ng) % rvbar(:,:,KOUT),                 &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv2d)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
#   ifdef SOLVE3D
#    ifdef FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRvct)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idRvct,          &
     &                     XTR(ng)%pioVar(idRvct),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#     endif
     &                     COUPLING(ng) % rvfrc,                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRvct)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
!
        IF (XTR(ng)%pioVar(idVfx1)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        Vr2d=ExtractFlag(ng)*COUPLING(ng)%DV_avg1
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVfx1,          &
     &                     XTR(ng)%pioVar(idVfx1),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx1)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (XTR(ng)%pioVar(idVfx2)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)

        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        Vr2d=ExtractFlag(ng)*COUPLING(ng)%DV_avg2
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVfx2,          &
     &                     XTR(ng)%pioVar(idVfx2),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVfx2)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Vr2d)
#   endif
#  endif
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, model, XTR(ng)%name,                &
     &                          XTR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          XTR(ng)%pioVar(idSbry(isVbar)),         &
     &                          XTR(ng)%Rindex,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)),     &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#  ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#  endif
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idu2dE)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idu2dE,          &
     &                     XTR(ng)%pioVar(idu2dE),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (XTR(ng)%pioVar(idv2dN)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idv2dN,          &
     &                     XTR(ng)%pioVar(idv2dN),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idUvel,          &
     &                     XTR(ng)%pioVar(idUvel),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#   endif
     &                     OCEAN(ng) % u(:,:,:,NOUT),                   &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#   if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRu3d)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idRu3d,          &
     &                     XTR(ng)%pioVar(idRu3d),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#    endif
     &                     OCEAN(ng) % ru(:,:,:,NOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRu3d)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, model, XTR(ng)%name,                &
     &                          XTR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          XTR(ng)%pioVar(idSbry(isUvel)),         &
     &                          XTR(ng)%Rindex,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)),        &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        XTR(ng)%Rindex
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
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idVvel,          &
     &                     XTR(ng)%pioVar(idVvel),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#   endif
     &                     OCEAN(ng) % v(:,:,:,NOUT),                   &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#   if defined FORWARD_WRITE && defined FORWARD_RHS
!
        IF (XTR(ng)%pioVar(idRv3d)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idRv3d,          &
     &                     XTR(ng)%pioVar(idRv3d),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#    endif
     &                     OCEAN(ng) % rv(:,:,:,NOUT),                  &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRv3d)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#   endif
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, model, XTR(ng)%name,                &
     &                          XTR(ng)%pioFile,                        &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          XTR(ng)%pioVar(idSbry(isVvel)),         &
     &                          XTR(ng)%Rindex,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)),        &
     &                          ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#   ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#   endif
     &                    OCEAN(ng) % u(:,:,:,NOUT),                    &
     &                    OCEAN(ng) % v(:,:,:,NOUT),                    &
     &                    Ur3d, Vr3d)
!
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idu3dE)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idu3dE,          &
     &                     XTR(ng)%pioVar(idu3dE),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     Ur3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (XTR(ng)%pioVar(idV3dN)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idv3dN,          &
     &                     XTR(ng)%pioVar(idv3dN),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#   endif
     &                     Vr3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
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
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
!
        IF (XTR(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idOvel,          &
     &                     XTR(ng)%pioVar(idOvel),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Wr3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF

#   ifdef OMEGA_IMPLICIT
!
!  Write out S-coordinate implicit omega vertical velocity (m/s).
!
      IF (Hout(idOvil,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate (Wr3d(LBi:UBi,LBj:UBj,0:N(ng)))
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % Wi,                               &
     &                    Wr3d)
!
        IF (XTR(ng)%pioVar(idOvil)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idOvil,          &
     &                     XTR(ng)%pioVar(idOvil),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     Wr3d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvil)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
#   endif
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idWvel)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idWvel,          &
     &                     XTR(ng)%pioVar(idWvel),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % wvel,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          IF (XTR(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesX => ioDesX_dp_r3dvar(ng)
          ELSE
            ioDesX => ioDesX_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idTvar(itrc),  &
     &                       XTR(ng)%pioTrc(itrc),                      &
     &                       XTR(ng)%Rindex,                            &
     &                       ioDesX,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc),            &
     &                       SetFillVal = .FALSE.,                      &
     &                       ExtractField = ExtractFlag(ng))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          XTR(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          status=nf_fwrite3d_bry (ng, model, XTR(ng)%name,              &
     &                            XTR(ng)%pioFile,                      &
     &                            Vname(1,ifield),                      &
     &                            XTR(ng)%pioVar(ifield),               &
     &                            XTR(ng)%Rindex,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng),itrc), &
     &                            ExtractField = ExtractFlag(ng))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)), XTR(ng)%Rindex
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
        IF (XTR(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idDano,          &
     &                     XTR(ng)%pioVar(idDano),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     OCEAN(ng) % rho,                             &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idHsbl)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idHsbl,          &
     &                     XTR(ng)%pioVar(idHsbl),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % hsbl,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHsbl)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   ifdef LMD_BKPP
!
!  Write out depth bottom boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idHbbl)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idHbbl,          &
     &                     XTR(ng)%pioVar(idHbbl),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % hbbl,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHbbl)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        IF (Hout(idGhat(i),ng)) THEN
          scale=1.0_dp
          IF (XTR(ng)%pioVar(idGhat(i))%dkind.eq.PIO_double) THEN
            ioDesX => ioDesX_dp_w3dvar(ng)
          ELSE
            ioDesX => ioDesX_sp_w3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idGhat(i),     &
     &                       XTR(ng)%pioVar(idGhat(i)),                 &
     &                       XTR(ng)%Rindex,                            &
     &                       ioDesX,                                    &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#    ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#    endif
     &                       MIXING(ng) % ghats(:,:,:,i),               &
     &                       SetFillVal = .FALSE.,                      &
     &                       ExtractField = ExtractFlag(ng))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idGhat(i))), XTR(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idVvis,          &
     &                     XTR(ng)%pioVar(idVvis),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), XTR(ng)%Rindex
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
        IF (XTR(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idTdif,          &
     &                     XTR(ng)%pioVar(idTdif),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), XTR(ng)%Rindex
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
        IF (XTR(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idSdif,          &
     &                     XTR(ng)%pioVar(idSdif),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idMtke)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idMtke,          &
     &                     XTR(ng)%pioVar(idMtke),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#    ifdef FORWARD_WRITE
!
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVmKK)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idVmKK,          &
     &                     XTR(ng)%pioVar(idVmKK),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     MIXING(ng) % Akk,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKK)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idMtls)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idMtls,          &
     &                     XTR(ng)%pioVar(idMtls),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF

#    ifdef FORWARD_WRITE
!
        IF (XTR(ng)%pioVar(idVmLS)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        scale=1.0_dp
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idVmLS,          &
     &                     XTR(ng)%pioVar(idVmLS),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     MIXING(ng) % Lscale,                         &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmLS)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
#    if defined FORWARD_WRITE && defined GLS_MIXING
!
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVmKP)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_w3dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, XTR(ng)%pioFile, idVmKP,          &
     &                     XTR(ng)%pioVar(idVmKP),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#     endif
     &                     MIXING(ng) % Akp,                            &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVmKP)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#    endif
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idPair)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idPair,          &
     &                     XTR(ng)%pioVar(idPair),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % Pair,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined BULK_FLUXES
!
!  Write out surface air temperature.
!
      IF (Hout(idTair,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idTair)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idTair,          &
     &                     XTR(ng)%pioVar(idTair),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % Tair,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idUair)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUair,          &
     &                     XTR(ng)%pioVar(idUair),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % Uwind,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Hout(idVair,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVair)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVair,          &
     &                     XTR(ng)%pioVar(idVair),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % Vwind,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Hout(idUaiE,ng).and.Hout(idVaiN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
          Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
          Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
#    ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#    endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idUaiE)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUaiE,          &
     &                     XTR(ng)%pioVar(idUaiE),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     Ur2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idVaiN)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVaiN,          &
     &                     XTR(ng)%pioVar(idVaiN),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     Vr2d,                                        &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#   endif
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
#   ifdef SO_SEMI
            scale=1.0_dp
#   else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
#   endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          IF (XTR(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesX => ioDesX_dp_r2dvar(ng)
          ELSE
            ioDesX => ioDesX_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idTsur(itrc),  &
     &                       XTR(ng)%pioVar(idTsur(itrc)),              &
     &                       XTR(ng)%Rindex,                            &
     &                       ioDesX,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % stflx(:,:,itrc),              &
     &                       SetFillVal = .FALSE.,                      &
     &                       ExtractField = ExtractFlag(ng))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          XTR(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

#   if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        IF (XTR(ng)%pioVar(idLhea)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idLhea,          &
     &                     XTR(ng)%pioVar(idLhea),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % lhflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        IF (XTR(ng)%pioVar(idShea)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idShea,          &
     &                     XTR(ng)%pioVar(idShea),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % shflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        IF (XTR(ng)%pioVar(idLrad)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idLrad,          &
     &                     XTR(ng)%pioVar(idLrad),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % lrflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif

#   ifdef BULK_FLUXES
#    ifdef EMINUSP
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idevap)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idevap,          &
     &                     XTR(ng)%pioVar(idevap),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % evap,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idrain)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idrain,          &
     &                     XTR(ng)%pioVar(idrain),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % rain,                           &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_dp
        IF (XTR(ng)%pioVar(idEmPf)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idEmPf,          &
     &                     XTR(ng)%pioVar(idEmPf),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % stflux(:,:,isalt),              &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        IF (XTR(ng)%pioVar(idSrad)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_r2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idSrad,          &
     &                     XTR(ng)%pioVar(idSrad),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % srflx,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#  ifdef SO_SEMI
        scale=1.0_dp
#  else
        scale=rho0                          ! m2/s2 to Pa
#  endif
        IF (XTR(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUsms,          &
     &                     XTR(ng)%pioVar(idUsms),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % sustr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), XTR(ng)%Rindex
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
#  ifdef SO_SEMI
        scale=1.0_dp
#  else
        scale=rho0
#  endif
        IF (XTR(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVsms,          &
     &                     XTR(ng)%pioVar(idVsms),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % svstr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        IF (XTR(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_u2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idUbms,          &
     &                     XTR(ng)%pioVar(idUbms),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % bustr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), XTR(ng)%Rindex
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
        scale=-rho0
        IF (XTR(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesX => ioDesX_dp_v2dvar(ng)
        ELSE
          ioDesX => ioDesX_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, XTR(ng)%pioFile, idVbms,          &
     &                     XTR(ng)%pioVar(idVbms),                      &
     &                     XTR(ng)%Rindex,                              &
     &                     ioDesX,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % bvstr,                          &
     &                     SetFillVal = .FALSE.,                        &
     &                     ExtractField = ExtractFlag(ng))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), XTR(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_pio (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_pio (ng, model, tile,                                &
    &                   LBi, UBi, LBj, UBj,                             &
    &                   Hout, XTR)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifndef DELAYED_SYNC_PIO
!
!-----------------------------------------------------------------------
!  Synchronize decimate NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, XTR(ng)%name, XTR(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
  10  FORMAT (2x,'WRT_EXTRACT_PIO  - writing extract', t42,             &
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
  20  FORMAT (/,' WRT_EXTRACT_PIO - error while writing variable: ',a,  &
     &        /,19x,'into decimate NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_extract_pio
# endif
#endif
!
      END MODULE wrt_extract_mod
