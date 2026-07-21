#include "cppdefs.h"
      MODULE wrt_quick_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine writes requested model fields into QUICKSAVE file      !
!  using the standard NetCDF library or the Parallel-IO (PIO) library. !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
      USE mod_stepping
!
#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
      USE bbl_output_mod,      ONLY : bbl_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,      ONLY : bbl_wrt_pio
# endif
#endif
#ifdef SOLVE3D
      USE extract_slice_mod,   ONLY : extract_slice
#endif
#if defined ICE_MODEL && defined SOLVE3D
      USE ice_output_mod,      ONLY : ice_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,      ONLY : ice_wrt_pio
# endif
#endif
      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
#ifdef SOLVE3D
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
      USE omega_mod,           ONLY : scale_omega
#endif
#if defined SEDIMENT && defined SOLVE3D
      USE sediment_output_mod, ONLY : sediment_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod, ONLY : sediment_wrt_pio
# endif
#endif
      USE strings_mod,         ONLY : FoundError
      USE uv_rotate_mod,       ONLY : uv_rotate2d
#ifdef SOLVE3D
      USE uv_rotate_mod,       ONLY : uv_rotate3d
#endif
#ifdef SOLVE3D
      USE vorticity_mod,       ONLY : pvorticity3d, rvorticity3d
#endif
#if defined WEC_VF && defined SOLVE3D
      USE wec_output_mod,      ONLY : wec_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,      ONLY : wec_wrt_pio
# endif
#endif
!
      implicit none
!
      PUBLIC  :: wrt_quick
      PRIVATE :: wrt_quick_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_quick_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_quick (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: IminS, ImaxS, JminS, JmaxS
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
      SELECT CASE (QCK(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_quick_nf90 (ng, iNLM, tile,                          &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS)

#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_quick_pio (ng, iNLM, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
#endif
        CASE DEFAULT
          IF (Master) THEN
            WRITE (stdout,10) QCK(ng)%IOtype
          END IF
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_QUICK - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_quick
!
!***********************************************************************
      SUBROUTINE wrt_quick_nf90 (ng, model, tile,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, gfactor, gtype, status
#ifdef SOLVE3D
      integer :: i, itrc, j, k
#endif
!
      real(dp) :: scale
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
#ifdef SOLVE3D
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_quick_nf90"

#include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out quicksave fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
      gfactor=-1
#else
      gfactor=1
#endif
!
!  Set time record index.
!
      QCK(ng)%Rindex=QCK(ng)%Rindex+1
      Fcount=QCK(ng)%load
      QCK(ng)%Nrec(Fcount)=QCK(ng)%Nrec(Fcount)+1
!
!  Report.
!
#ifdef SOLVE3D
# ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, NOUT, QCK(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,10) KOUT, NOUT, QCK(ng)%Rindex
# endif
#else
# ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, QCK(ng)%Rindex, ng
# else
      IF (Master) WRITE (stdout,10) KOUT, QCK(ng)%Rindex
# endif
#endif
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, QCK(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/QCK(ng)%Rindex/), (/1/),                  &
     &                      ncid = QCK(ng)%ncid,                        &
     &                      varid = QCK(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef WET_DRY
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_dp
      gtype=gfactor*p2dvar
      status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idPwet,               &
     &                   QCK(ng)%Vid(idPwet),                           &
     &                   QCK(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % pmask,                              &
# endif
     &                   GRID(ng) % pmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idPwet)), QCK(ng)%Rindex
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
      status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idRwet,               &
     &                   QCK(ng)%Vid(idRwet),                           &
     &                   QCK(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idRwet)), QCK(ng)%Rindex
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
      status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUwet,               &
     &                   QCK(ng)%Vid(idUwet),                           &
     &                   QCK(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUwet)), QCK(ng)%Rindex
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
      status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVwet,               &
     &                   QCK(ng)%Vid(idVwet),                           &
     &                   QCK(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVwet)), QCK(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
#ifdef SOLVE3D
!
!  Write time-varying depths of RHO-points.
!
      IF (Qout(idpthR,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idpthR,             &
     &                     QCK(ng)%Vid(idpthR),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Qout(idpthU,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u3dvar
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idpthU,             &
     &                     QCK(ng)%Vid(idpthU),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Qout(idpthV,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v3dvar
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idpthV,             &
     &                     QCK(ng)%Vid(idpthV),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Qout(idpthW,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idpthW,             &
     &                     QCK(ng)%Vid(idpthW),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out free-surface (m)
!
      IF (Qout(idFsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idFsur,             &
     &                     QCK(ng)%Vid(idFsur),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#endif
#ifdef WET_DRY
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.)
#else
     &                     OCEAN(ng) % zeta(:,:,KOUT))
#endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Qout(idUbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUbar,             &
     &                     QCK(ng)%Vid(idUbar),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Qout(idVbar,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVbar,             &
     &                     QCK(ng)%Vid(idVbar),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Qout(idu2dE,ng).and.Qout(idv2dN,ng)) THEN
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
#ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
#endif
     &                    OCEAN(ng) % ubar(:,:,KOUT),                   &
     &                    OCEAN(ng) % vbar(:,:,KOUT),                   &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idu2dE,             &
     &                     QCK(ng)%Vid(idu2dE),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idv2dN,             &
     &                     QCK(ng)%Vid(idv2dN),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

#ifdef SOLVE3D
!
!  Write out surface U-momentum component (m/s).
!
      IF (Qout(idUsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUsur,             &
     &                     QCK(ng)%Vid(idUsur),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     OCEAN(ng) % u(:,:,N(ng),NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsur)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum component (m/s).
!
      IF (Qout(idVsur,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVsur,             &
     &                     QCK(ng)%Vid(idVsur),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     OCEAN(ng) % v(:,:,N(ng),NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsur)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Qout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*u3dvar
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idUzsl,             &
     &                     QCK(ng)%Vid(idUzsl),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Qout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*v3dvar
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idVzsl,             &
     &                     QCK(ng)%Vid(idVzsl),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idu3dE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idu3dE,             &
     &                     QCK(ng)%Vid(idu3dE),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idv3dN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idv3dN,             &
     &                     QCK(ng)%Vid(idv3dN),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     OCEAN(ng) % va)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idUsuE,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUsuE,             &
     &                     QCK(ng)%Vid(idUsuE),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     OCEAN(ng) % ua(:,:,N(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsuE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idVsuN,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVsuN,           &
     &                     QCK(ng)%Vid(idVsuN),                       &
     &                     QCK(ng)%Rindex, gtype,                     &
     &                     LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                     &
# endif
     &                     OCEAN(ng) % va(:,:,N(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsuN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Qout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idUzsE,             &
     &                     QCK(ng)%Vid(idUzsE),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Qout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*r3dvar
        CALL extract_slice (ng, model, tile, gtype,                     &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idVzsN,             &
     &                     QCK(ng)%Vid(idVzsN),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (m-1 s-1) at PSI-points at specified
!  constant depth slices.
!
      IF (Qout(idPVzs,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrP,JendT
            DO i=IstrP,IendT
              GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+     &
     &                                     GRID(ng)%z_r(i-1,j  ,k)+     &
     &                                     GRID(ng)%z_r(i  ,j-1,k)+     &
     &                                     GRID(ng)%z_r(i ,j  ,k))
            END DO
          END DO
        END DO
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, NOUT,            &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
# endif
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        CALL extract_slice (ng, model, tile, p3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fr3d,                                       &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % pmask,                           &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idPVzs,             &
     &                     QCK(ng)%Vid(idPVzs),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPVzs)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Fr3d)
        deallocate (Wr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at specified
!  constant depth slices.
!
      IF (Qout(idRVzs,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrP,JendT
            DO i=IstrP,IendT
              GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+     &
     &                                     GRID(ng)%z_r(i-1,j  ,k)+     &
     &                                     GRID(ng)%z_r(i  ,j-1,k)+     &
     &                                     GRID(ng)%z_r(i ,j  ,k))
            END DO
          END DO
        END DO
!
        CALL rvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj, NOUT,                    &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        CALL extract_slice (ng, model, tile, p3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fr3d,                                       &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % pmask,                           &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        gtype=gfactor*p3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idRVzs,             &
     &                     QCK(ng)%Vid(idRVzs),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRVzs)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Fr3d)
        deallocate (Wr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Qout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        gtype=gfactor*w3dvar
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idOvel,             &
     &                     QCK(ng)%Vid(idOvel),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Qout(idWvel,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idWvel,             &
     &                     QCK(ng)%Vid(idWvel),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), QCK(ng)%Rindex
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
        IF (Qout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idTvar(itrc),     &
     &                       QCK(ng)%Tid(itrc),                         &
     &                       QCK(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out surface tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Qout(idsurT(itrc),ng)) THEN
          scale=1.0_dp
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idsurT(itrc),     &
     &                       QCK(ng)%Vid(idsurT(itrc)),                 &
     &                       QCK(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % t(:,:,N(ng),NOUT,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idsurT(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Qout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          scale=1.0_dp
          gtype=gfactor*r3dvar
          CALL extract_slice (ng, model, tile, gtype,                   &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,NOUT,itrc),           &
     &                        GRID(ng) % z_r,                           &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        Zslice, Wr3d)
!
          status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idzslT(itrc),     &
     &                       QCK(ng)%Vid(idzslT(itrc)),                 &
     &                       QCK(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Qout(idDano,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idDano,             &
     &                     QCK(ng)%Vid(idDano),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Qout(idHsbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idHsbl,             &
     &                     QCK(ng)%Vid(idHsbl),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % hsbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHsbl)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth surface boundary layer.
!
      IF (Qout(idHbbl,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idHbbl,             &
     &                     QCK(ng)%Vid(idHbbl),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % hbbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHbbl)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Qout(idVvis,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idVvis,             &
     &                     QCK(ng)%Vid(idVvis),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Qout(idTdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idTdif,             &
     &                     QCK(ng)%Vid(idTdif),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Qout(idSdif,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idSdif,             &
     &                     QCK(ng)%Vid(idSdif),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Qout(idMtke,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idMtke,             &
     &                     QCK(ng)%Vid(idMtke),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Qout(idMtls,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, model, QCK(ng)%ncid, idMtls,             &
     &                     QCK(ng)%Vid(idMtls),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)

        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Qout(idPair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idPair,             &
     &                     QCK(ng)%Vid(idPair),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Pair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
!
!  Write atmosphere boundary layer thickness.
!
      IF (Qout(idPBLh,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idPBLh,             &
     &                     QCK(ng)%Vid(idPBLh),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % PBLh)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPBLh)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface air temperature.
!
      IF (Qout(idTair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idTair,             &
     &                     QCK(ng)%Vid(idTair),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Tair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || \
     defined ATM2OCN_FLUXES || defined SPECTRAL_LIGHT
!
!  Write out surface winds.
!
      IF (Qout(idUair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUair,             &
     &                     QCK(ng)%Vid(idUair),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Uwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Qout(idVair,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVair,             &
     &                     QCK(ng)%Vid(idVair),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % Vwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Qout(idUaiE,ng).and.Qout(idVaiN,ng)) THEN
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
     &                    GRID(ng) % rmask,                             &
#  endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUaiE,             &
     &                     QCK(ng)%Vid(idUaiE),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVaiN,             &
     &                     QCK(ng)%Vid(idVaiN),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
# endif
!
!  Write out surface active tracers fluxes.
!
      DO itrc=1,NAT
        IF (Qout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
# ifdef SO_SEMI
            scale=1.0_dp
# else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
# endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idTsur(itrc),     &
     &                       QCK(ng)%Vid(idTsur(itrc)),                 &
     &                       QCK(ng)%Rindex, gtype,                     &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

# if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Write out latent heat flux.
!
      IF (Qout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idLhea,             &
     &                     QCK(ng)%Vid(idLhea),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % lhflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Qout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idShea,             &
     &                     QCK(ng)%Vid(idShea),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % shflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Qout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idLrad,             &
     &                     QCK(ng)%Vid(idLrad),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % lrflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef EMINUSP
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Qout(idevap,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idevap,             &
     &                     QCK(ng)%Vid(idevap),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % evap)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Qout(idrain,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idrain,             &
     &                     QCK(ng)%Vid(idrain),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % rain)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out E-P (m/s).
!
      IF (Qout(idEmPf,ng)) THEN
        scale=1.0_dp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idEmPf,             &
     &                     QCK(ng)%Vid(idEmPf),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SHORTWAVE
!
!  Write out net shortwave radiation flux.
!
      IF (Qout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idSrad,             &
     &                     QCK(ng)%Vid(idSrad),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % srflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
!
!  Write out surface U-momentum stress.
!
      IF (Qout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUsms,             &
     &                     QCK(ng)%Vid(idUsms),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Qout(idVsms,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVsms,             &
     &                     QCK(ng)%Vid(idVsms),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Qout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idUbms,             &
     &                     QCK(ng)%Vid(idUbms),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask,                            &
#endif
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Qout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, model, QCK(ng)%ncid, idVbms,             &
     &                     QCK(ng)%Vid(idVbms),                         &
     &                     QCK(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#endif
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_nf90 (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_nf90 (ng, model, tile,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_nf90 (ng, model, tile,                               &
    &                    LBi, UBi, LBj, UBj,                            &
    &                    Qout, QCK)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#if !(defined OUT_NETCDF4 && defined DELAYED_SYNC_NF90)
!
!-----------------------------------------------------------------------
!  Synchronize quicksave NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, QCK(ng)%name, QCK(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
  10  FORMAT (2x,'WRT_QUICK_NF90   - writing quicksave', t42,           &
#ifdef SOLVE3D
# ifdef NESTING
     &        'fields (Index=',i1,',',i1,') in record = ',i0,t92,i2.2)
# else
     &        'fields (Index=',i1,',',i1,') in record = ',i0)
# endif
#else
# ifdef NESTING
     &        'fields (Index=',i1,')   in record = ',i0,t92,i2.2)
# else
     &        'fields (Index=',i1,')   in record = ',i0)
# endif
#endif
  20  FORMAT (/,' WRT_QUICK_NF90 - error while writing variable: ',a,   &
     &        /,18x,'into quicksave NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_quick_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_quick_pio (ng, model, tile,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Fcount, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k
# endif
!
      real(dp) :: scale
!
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
# ifdef SOLVE3D
      real(r8), allocatable :: Fr3d(:,:,:)
      real(r8), allocatable :: Wr3d(:,:,:)
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_quick_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc

#  include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out quicksave fields.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      QCK(ng)%Rindex=QCK(ng)%Rindex+1
      Fcount=QCK(ng)%load
      QCK(ng)%Nrec(Fcount)=QCK(ng)%Nrec(Fcount)+1
!
!  Report.
!
# ifdef SOLVE3D
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, NOUT, QCK(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) KOUT, NOUT, QCK(ng)%Rindex
#  endif
# else
#  ifdef NESTING
      IF (Master) WRITE (stdout,10) KOUT, QCK(ng)%Rindex, ng
#  else
      IF (Master) WRITE (stdout,10) KOUT, QCK(ng)%Rindex
#  endif
# endif
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, model, QCK(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/QCK(ng)%Rindex/), (/1/),              &
     &                          pioFile = QCK(ng)%pioFile,              &
     &                          pioVar = QCK(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef WET_DRY
!
!  Write out wet/dry mask at PSI-points.
!
      scale=1.0_dp
      IF (QCK(ng)%pioVar(idPwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_p2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_p2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idPwet,            &
     &                   QCK(ng)%pioVar(idPwet),                        &
     &                   QCK(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % pmask,                              &
#  endif
     &                   GRID(ng) % pmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idPwet)), QCK(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_dp
      IF (QCK(ng)%pioVar(idRwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idRwet,            &
     &                   QCK(ng)%pioVar(idRwet),                        &
     &                   QCK(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idRwet)), QCK(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_dp
      IF (QCK(ng)%pioVar(idUwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUwet,            &
     &                   QCK(ng)%pioVar(idUwet),                        &
     &                   QCK(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUwet)), QCK(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_dp
      IF (QCK(ng)%pioVar(idVwet)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
      status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVwet,            &
     &                   QCK(ng)%pioVar(idVwet),                        &
     &                   QCK(ng)%Rindex,                                &
     &                   ioDesc,                                        &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVwet)), QCK(ng)%Rindex
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
      IF (Qout(idpthR,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idpthR)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idpthR,          &
     &                     QCK(ng)%pioVar(idpthR),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     GRID(ng) % z_r,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthR)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of U-points.
!
      IF (Qout(idpthU,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        IF (QCK(ng)%pioVar(idpthU)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idpthU,          &
     &                     QCK(ng)%pioVar(idpthU),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthU)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of V-points.
!
      IF (Qout(idpthV,ng)) THEN
        scale=1.0_dp
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        IF (QCK(ng)%pioVar(idpthV)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idpthV,          &
     &                     QCK(ng)%pioVar(idpthV),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % z_v,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthV)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write time-varying depths of W-points.
!
      IF (Qout(idpthW,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idpthW)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idpthW,          &
     &                     QCK(ng)%pioVar(idpthW),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     GRID(ng) % z_w,                              &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idpthW)), QCK(ng)%Rindex
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
      IF (Qout(idFsur,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idFsur,          &
     &                     QCK(ng)%pioVar(idFsur),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.)
# else
     &                     OCEAN(ng) % zeta(:,:,KOUT))
# endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idFsur)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Qout(idUbar,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUbar,          &
     &                     QCK(ng)%pioVar(idUbar),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
# endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbar)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Qout(idVbar,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVbar,          &
     &                     QCK(ng)%pioVar(idVbar),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
# endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbar)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Qout(idu2dE,ng).and.Qout(idv2dN,ng)) THEN
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
        IF (QCK(ng)%pioVar(idu2dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idu2dE,          &
     &                     QCK(ng)%pioVar(idu2dE),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Ur2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu2dE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        IF (QCK(ng)%pioVar(idv2dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idv2dN,          &
     &                     QCK(ng)%pioVar(idv2dN),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
# endif
     &                     Vr2d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv2dN)), QCK(ng)%Rindex
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
      IF (Qout(idUvel,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idUvel,          &
     &                     QCK(ng)%pioVar(idUvel),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % u(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUvel)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Qout(idVvel,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idVvel,          &
     &                     QCK(ng)%pioVar(idVvel),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % v(:,:,:,NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvel)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface U-momentum component (m/s).
!
      IF (Qout(idUsur,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUsur,          &
     &                     QCK(ng)%pioVar(idUsur),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     OCEAN(ng) % u(:,:,N(ng),NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsur)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum component (m/s).
!
      IF (Qout(idVsur,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVsur)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVsur,          &
     &                     QCK(ng)%pioVar(idVsur),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     OCEAN(ng) % v(:,:,N(ng),NOUT))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsur)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Qout(idUzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=Jstr-1,Jend+1
            DO i=IstrU-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i-1,j,k)+        &
     &                                    GRID(ng)%z_r(i  ,j,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, u3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % u(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUzsl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_u3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_u3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idUvel,          &
     &                     QCK(ng)%pioVar(idUzsl),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % umask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsl)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D V-momentum component (m/s) at specified constant depth
!  slices.
!
      IF (Qout(idVzsl,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrV-1,Jend+1
            DO i=Istr-1,Iend+1
              GRID(ng)%z_v(i,j,k)=0.5_r8*(GRID(ng)%z_r(i,j-1,k)+        &
     &                                    GRID(ng)%z_r(i,j  ,k))
            END DO
          END DO
        END DO
        CALL extract_slice (ng, model, tile, v3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % v(:,:,:,NOUT),                  &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVzsl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_v3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_v3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idVzsl,          &
     &                     QCK(ng)%pioVar(idVzsl),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % vmask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVzsl)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idu3dE,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idu3dE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idu3dE,          &
     &                     QCK(ng)%pioVar(idu3dE),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % ua)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idu3dE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idv3dN,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idv3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idv3dN,          &
     &                     QCK(ng)%pioVar(idv3dN),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % va)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idv3dN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface Eastward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idUsuE,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUsuE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUsuE,          &
     &                     QCK(ng)%pioVar(idUsuE),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % ua(:,:,N(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsuE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface Northward momentum (m/s) at RHO-points, A-grid.
!
      IF (Qout(idVsuN,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVsuN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVsuN,          &
     &                     QCK(ng)%pioVar(idVsuN),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     OCEAN(ng) % va(:,:,N(ng)))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsuN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D Eastward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Qout(idUzsE,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        CALL extract_slice (ng, model, tile, r3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % ua,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUzsE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idUzsE,          &
     &                     QCK(ng)%pioVar(idUzsE),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUzsE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out 3D Northward momentum (m/s) at RHO-points, A-grid, at
!  specified constant depth slices.
!
      IF (Qout(idVzsN,ng).and.(Nslice.gt.0)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
          Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
        END IF
        CALL extract_slice (ng, model, tile, r3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      OCEAN(ng) % va,                             &
     &                      GRID(ng) % z_r,                             &
# ifdef MASKING
     &                      GRID(ng) % rmask_full,                      &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idv3dN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idVzsN,          &
     &                     QCK(ng)%pioVar(idVzsN),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % rmask_full,                       &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVzsN)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out potential vorticity (1/s) at PSI-points at specified
!  constant depth slices.
!
      IF (Qout(idPVzs,ng)) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrP,JendT
            DO i=IstrP,IendT
              GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+     &
     &                                     GRID(ng)%z_r(i-1,j  ,k)+     &
     &                                     GRID(ng)%z_r(i  ,j-1,k)+     &
     &                                     GRID(ng)%z_r(i ,j  ,k))
            END DO
          END DO
        END DO
!
        CALL pvorticity3d (ng, model, tile,                             &
                           LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS, NOUT,            &
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % umask,   GRID(ng) % vmask,        &
# endif
                           GRID(ng) % f,                                &
                           GRID(ng) % om_u,    GRID(ng) % on_v,         &
                           GRID(ng) % pm,      GRID(ng) % pn,           &
                           GRID(ng) % z_r,     OCEAN(ng) % pden,        &
                           OCEAN(ng) % u,      OCEAN(ng) % v,           &
                           Fr3d)
!
        CALL extract_slice (ng, model, tile, p3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fr3d,                                       &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % pmask,                           &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idPVzs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_p3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_p3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idPVzs,          &
     &                     QCK(ng)%pioVar(idPVzs),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPVzs)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Fr3d)
        deallocate (Wr3d)
      END IF
!
!  Write out relative vorticity (s-1) at PSI-points at specified
!  constant depth slices.
!
      IF (Qout(idRVzs,ng)) THEN
        IF (.not.allocated(Fr3d)) THEN
          allocate ( Fr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Fr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
!
        DO k=1,N(ng)
          DO j=JstrP,JendT
            DO i=IstrP,IendT
              GRID(ng)%z_v(i,j,k)=0.25_r8*(GRID(ng)%z_r(i-1,j-1,k)+     &
     &                                     GRID(ng)%z_r(i-1,j  ,k)+     &
     &                                     GRID(ng)%z_r(i  ,j-1,k)+     &
     &                                     GRID(ng)%z_r(i ,j  ,k))
            END DO
          END DO
        END DO
!
        CALL rvorticity3d (ng, model, tile,                             &
     &                     LBi, UBi, LBj, UBj, NOUT,                    &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
                           GRID(ng) % om_u, GRID(ng) % on_v,            &
                           GRID(ng) % pm,   GRID(ng) % pn,              &
                           OCEAN(ng) % u,   OCEAN(ng) % v,              &
                           Fr3d)
!
        CALL extract_slice (ng, model, tile, p3dvar,                    &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fr3d,                                       &
     &                      GRID(ng) % z_v,                             &
# ifdef MASKING
     &                      GRID(ng) % pmask,                           &
# endif
     &                      Zslice, Wr3d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idRVzs)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesZ_dp_p3dvar(ng)
        ELSE
          ioDesc => ioDesZ_sp_p3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idRVzs,          &
     &                     QCK(ng)%pioVar(idRVzs),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nslice, scale,        &
#  ifdef MASKING
     &                     GRID(ng) % pmask,                            &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idRVzs)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Fr3d)
        deallocate (Wr3d)
      END IF
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Qout(idOvel,ng)) THEN
        IF (.not.allocated(Wr3d)) THEN
          allocate ( Wr3d(LBi:UBi,LBj:UBj,0:N(ng)) )
          Wr3d(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_dp
        CALL scale_omega (ng, tile, LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    OCEAN(ng) % W,                                &
     &                    Wr3d)
!
        IF (QCK(ng)%pioVar(idOvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idOvel,          &
     &                     QCK(ng)%pioVar(idOvel),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     Wr3d)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idOvel)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (Wr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Qout(idWvel,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idWvel)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idWvel,          &
     &                     QCK(ng)%pioVar(idWvel),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % wvel)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idWvel)), QCK(ng)%Rindex
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
        IF (Qout(idTvar(itrc),ng)) THEN
          scale=1.0_dp
          IF (QCK(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idTvar(itrc),  &
     &                       QCK(ng)%pioTrc(itrc),                      &
     &                       QCK(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out surface tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Qout(idsurT(itrc),ng)) THEN
          scale=1.0_dp
          IF (QCK(ng)%pioVar(idsurT(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idsurT(itrc),  &
     &                       QCK(ng)%pioVar(idsurT(itrc)),              &
     &                       QCK(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % t(:,:,N(ng),NOUT,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idsurT(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out tracer type variables at specified constant depth slices.
!
      DO itrc=1,NT(ng)
        IF (Qout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
          IF (.not.allocated(Wr3d)) THEN
            allocate ( Wr3d(LBi:UBi,LBj:UBj,Nslice) )
            Wr3d(LBi:UBi,LBj:UBj,1:Nslice)=0.0_r8
          END IF
          CALL extract_slice (ng, model, tile, r3dvar,                  &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        OCEAN(ng) % t(:,:,:,NOUT,itrc),           &
     &                        GRID(ng) % z_r,                           &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        Zslice, Wr3d)
!
          scale=1.0_dp
          IF (QCK(ng)%pioVar(idzslT(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesZ_dp_r3dvar(ng)
          ELSE
            ioDesc => ioDesZ_sp_r3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idzslT(itrc),  &
     &                       QCK(ng)%pioVar(idzslT(itrc)),              &
     &                       QCK(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, Nslice, scale,      &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       Wr3d)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idzslT(itrc))),            &
     &                          QCK(ng)%Rindex
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
          IF (itrc.eq.NT(ng)) deallocate (Wr3d)
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Qout(idDano,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idDano)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idDano,          &
     &                     QCK(ng)%pioVar(idDano),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % rho)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idDano)), QCK(ng)%Rindex
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
      IF (Qout(idHsbl,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idHsbl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idHsbl,          &
     &                     QCK(ng)%pioVar(idHsbl),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % hsbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHsbl)), QCK(ng)%Rindex
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
      IF (Qout(idHbbl,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idHbbl)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idHbbl,          &
     &                     QCK(ng)%pioVar(idHbbl),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % hbbl)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idHbbl)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Qout(idVvis,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVvis)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idVvis,          &
     &                     QCK(ng)%pioVar(idVvis),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVvis)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Qout(idTdif,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idTdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idTdif,          &
     &                     QCK(ng)%pioVar(idTdif),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTdif)), QCK(ng)%Rindex
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
      IF (Qout(idSdif,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idSdif)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idSdif,          &
     &                     QCK(ng)%pioVar(idSdif),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSdif)), QCK(ng)%Rindex
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
      IF (Qout(idMtke,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idMtke)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idMtke,          &
     &                     QCK(ng)%pioVar(idMtke),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtke)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Qout(idMtls,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idMtls)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_w3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_w3dvar(ng)
        END IF
        status=nf_fwrite3d(ng, model, QCK(ng)%pioFile, idMtls,          &
     &                     QCK(ng)%pioVar(idMtls),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)

        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idMtls)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Qout(idPair,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idPair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idPair,          &
     &                     QCK(ng)%pioVar(idPair),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Pair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
!
!  Write atmosphere boundary layer thickness.
!
      IF (Qout(idPBLh,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idPBLh)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idPBLh,          &
     &                     QCK(ng)%pioVar(idPBLh),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % PBLh)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idPBLh)), QCK(ng)%Rindex
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
      IF (Qout(idTair,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idTair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idTair,          &
     &                     QCK(ng)%pioVar(idTair),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Tair)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTair)), QCK(ng)%Rindex
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
      IF (Qout(idUair,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUair,          &
     &                     QCK(ng)%pioVar(idUair),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Uwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
      IF (Qout(idVair,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVair)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVair,          &
     &                     QCK(ng)%pioVar(idVair),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % Vwind)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVair)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Eastward/Northward surface wind (m/s) at RHO-points.
!
      IF (Qout(idUaiE,ng).and.Qout(idVaiN,ng)) THEN
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
     &                    GRID(ng) % rmask,                             &
#   endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idUaiE)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUaiE,          &
     &                     QCK(ng)%pioVar(idUaiE),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Ur2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUaiE)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idVaiN)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVaiN,          &
     &                     QCK(ng)%pioVar(idVaiN),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     Vr2d)
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVaiN)), QCK(ng)%Rindex
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
        IF (Qout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
#  ifdef SO_SEMI
            scale=1.0_dp
#  else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
#  endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_dp
          END IF
          IF (QCK(ng)%pioVar(idTsur(itrc))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idTsur(itrc),  &
     &                       QCK(ng)%pioVar(idTsur(itrc)),              &
     &                       QCK(ng)%Rindex,                            &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          QCK(ng)%Rindex
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
      IF (Qout(idLhea,ng)) THEN
        scale=rho0*Cp
        IF (QCK(ng)%pioVar(idLhea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idLhea,          &
     &                     QCK(ng)%pioVar(idLhea),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % lhflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLhea)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Qout(idShea,ng)) THEN
        scale=rho0*Cp
        IF (QCK(ng)%pioVar(idShea)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idShea,          &
     &                     QCK(ng)%pioVar(idShea),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % shflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idShea)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out net longwave radiation flux.
!
      IF (Qout(idLrad,ng)) THEN
        scale=rho0*Cp
        IF (QCK(ng)%pioVar(idLrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idLrad,          &
     &                     QCK(ng)%pioVar(idLrad),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % lrflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idLrad)), QCK(ng)%Rindex
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
      IF (Qout(idevap,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idevap)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idevap,          &
     &                     QCK(ng)%pioVar(idevap),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % evap)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idevap)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Qout(idrain,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idrain)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idrain,          &
     &                     QCK(ng)%pioVar(idrain),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     FORCES(ng) % rain)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idrain)), QCK(ng)%Rindex
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
      IF (Qout(idEmPf,ng)) THEN
        scale=1.0_dp
        IF (QCK(ng)%pioVar(idEmPf)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idEmPf,          &
     &                     QCK(ng)%pioVar(idEmPf),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     FORCES(ng) % stflux(:,:,isalt))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idEmPf)), QCK(ng)%Rindex
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
      IF (Qout(idSrad,ng)) THEN
        scale=rho0*Cp
        IF (QCK(ng)%pioVar(idSrad)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idSrad,          &
     &                     QCK(ng)%pioVar(idSrad),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     FORCES(ng) % srflx)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSrad)), QCK(ng)%Rindex
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
      IF (Qout(idUsms,ng)) THEN
        scale=rho0                          ! m2/s2 to Pa
        IF (QCK(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUsms,          &
     &                     QCK(ng)%pioVar(idUsms),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % sustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUsms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Qout(idVsms,ng)) THEN
        scale=rho0
        IF (QCK(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVsms,          &
     &                     QCK(ng)%pioVar(idVsms),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % svstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVsms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Qout(idUbms,ng)) THEN
        scale=-rho0
        IF (QCK(ng)%pioVar(idUbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idUbms,          &
     &                     QCK(ng)%pioVar(idUbms),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     FORCES(ng) % bustr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idUbms)), QCK(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Qout(idVbms,ng)) THEN
        scale=-rho0
        IF (QCK(ng)%pioVar(idVbms)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        status=nf_fwrite2d(ng, model, QCK(ng)%pioFile, idVbms,          &
     &                     QCK(ng)%pioVar(idVbms),                      &
     &                     QCK(ng)%Rindex,                              &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     FORCES(ng) % bvstr)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idVbms)), QCK(ng)%Rindex
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
        CALL bbl_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_wrt_pio (ng, model, tile,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_pio (ng, model, tile,                         &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Write out the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_pio (ng, model, tile,                                &
    &                   LBi, UBi, LBj, UBj,                             &
    &                   Qout, QCK)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifndef DELAYED_SYNC_PIO
!
!-----------------------------------------------------------------------
!  Synchronize quicksave NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, QCK(ng)%name, QCK(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
  10  FORMAT (2x,'WRT_QUICK_PIO    - writing quicksave', t42,           &
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
  20  FORMAT (/,' WRT_QUICK_PIO - error while writing variable: ',a,    &
     &        /,17x,'into quicksave NetCDF file for time record: ',i0)
!
      RETURN
      END SUBROUTINE wrt_quick_pio
#endif
      END MODULE wrt_quick_mod
