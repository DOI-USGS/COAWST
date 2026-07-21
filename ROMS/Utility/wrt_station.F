#include "cppdefs.h"
      MODULE wrt_station_mod
#ifdef STATIONS
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module writes STATION data into output file using the standard !
!  NetCDF library or the Parallel-IO (PIO) library.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BBL_MODEL
      USE mod_bbl
# endif
      USE mod_forces
      USE mod_grid
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
# if defined BBL_MODEL || defined WAVES_OUTPUT
      USE bbl_output_mod,      ONLY : bbl_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,      ONLY : bbl_wrt_station_pio
#  endif
# endif
      USE extract_sta_mod,     ONLY : extract_sta2d
# ifdef SOLVE3D
      USE extract_sta_mod,     ONLY : extract_sta3d
# endif
# ifdef ICE_MODEL
      USE ice_output_mod,      ONLY : ice_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,      ONLY : ice_wrt_station_pio
#  endif
# endif
# ifdef SEDIMENT
      USE sediment_output_mod, ONLY : sediment_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod, ONLY : sediment_wrt_station_pio
#  endif
# endif
      USE uv_rotate_mod,       ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod,       ONLY : uv_rotate3d
# endif
      USE strings_mod,         ONLY : FoundError
# if defined WEC || defined WEC_VF
      USE wec_output_mod,      ONLY : wec_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,      ONLY : wec_wrt_station_pio
#  endif
# endif
!
      implicit none
!
      PUBLIC  :: wrt_station
      PRIVATE :: wrt_station_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_station_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_station (ng, tile)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
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
      SELECT CASE (STA(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_station_nf90 (ng, iNLM, tile,                        &
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_station_pio (ng, iNLM, tile,                         &
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) STA(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_STATION - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_station
!
!***********************************************************************
      SUBROUTINE wrt_station_nf90 (ng, model, tile,                     &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB, NposR, NposW
      integer :: Fcount, i, ifield, k, np, status
!
      real(dp) :: scale

      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
# ifdef SOLVE3D
#  ifdef SEDIMENT
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
#  endif
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
# endif

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
# ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
# endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      STA(ng)%Rindex=STA(ng)%Rindex+1
      Fcount=STA(ng)%Fcount
      STA(ng)%Nrec(Fcount)=STA(ng)%Nrec(Fcount)+1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
# ifdef STATIONS_CGRID
      Cgrid=.TRUE.
# else
      Cgrid=.FALSE.
# endif
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
# ifdef SOLVE3D
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
#  ifdef SEDIMENT
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
#  endif
# endif
      END DO
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, model, STA(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/STA(ng)%Rindex/), (/1/),                  &
     &                      ncid = STA(ng)%ncid,                        &
     &                      varid = STA(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idFsur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idFsur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idFsur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbar, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbar))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbar, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbar))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
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
        CALL extract_sta2d (ng, model, Cgrid, idu2dE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idu2dE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu2dE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idv2dN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idv2dN)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv2dN))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idUvel, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvel, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
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
        CALL extract_sta3d (ng, model, Cgrid, idu3dE, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idu3dE)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu3dE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta3d (ng, model, Cgrid, idv3dN, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idv3dN)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv3dN))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idWvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idWvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idOvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idOvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idOvel))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, ifield, r3dvar,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idTvar(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/N(ng),Nstation(ng),1/),               &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Tid(i))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idDano, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idDano)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idDano))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef NEMURO_SED1
!
!  Write out PON in sediment.
!
      IF (Sout(idPONsed,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idPONsed, r2dvar,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%PONsed,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idPONsed)), psta,            &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idPONsed))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out OPAL in sediment.
!
      IF (Sout(idOPALsed,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idOPALsed, r2dvar,         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%OPALsed,                   &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idOPALsed)), psta,           &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idOPALsed))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out DENIT in sediment.
!
      IF (Sout(idDENITsed,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idDENITsed, r2dvar,        &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%DENITsed,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idDENITsed)), psta,          &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idDENITsed))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out buried PON in sediment.
!
      IF (Sout(idPONbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idPONbur, r2dvar,          &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%PON_burial,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idPONbur)), psta,            &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idPONbur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out buried OPAL in sediment.
!
      IF (Sout(idOPALbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idOPALbur, r2dvar,         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%OPAL_burial,               &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idOPALbur)), psta,           &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idOPALbur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef PRIMARY_PROD
!
!  Write out Net primary productivity.
!
      IF (Sout(idNPP,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idNPP, r2dvar,             &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%Bio_NPP,                   &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idNPP)), psta,               &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idNPP))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHsbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idHsbl)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHsbl))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHbbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idHbbl)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idHbbl))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvis, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVvis)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvis))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idTdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idTdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idSdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idSdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtke, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idMtke)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idMtke))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtls, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idMtls)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idMtls))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Sout(idPair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idPair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Pair,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idPair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idPair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Sout(idUair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Uwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      IF (Sout(idVair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Vwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward surface winds (m/s) at
!  RHO-points
!
      IF (Sout(idUaiE,ng).and.Sout(idVaiN,ng)) THEN
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
        CALL extract_sta2d (ng, model, Cgrid, idUaiE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUaiE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUaiE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idVaiN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVaiN)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVaiN))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#  endif
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef SALINITY
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLhea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idLhea)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idLhea))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idShea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idShea)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idShea))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idLrad)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idLrad))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idSrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idSrad)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSrad))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined EMINUSP && defined BULK_FLUXES
!
!  Write out E-P (m/s).
!
      IF (Sout(idEmPf,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idEmPf, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflux(:,:,isalt),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idEmPf)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idEmPf))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idevap, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idevap)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idevap))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idrain, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idrain)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idrain))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
# endif
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUsms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUsms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVsms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVsms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      IF (Sout(idRwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idRwet, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%rmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idRwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idRwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      IF (Sout(idUwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUwet, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%umask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idUwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      IF (Sout(idVwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVwet, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%vmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, model, STA(ng)%name,                  &
     &                        TRIM(Vname(1,idVwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif

# if defined BBL_MODEL || defined WAVES_OUTPUT
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
      CALL bbl_wrt_station_nf90 (ng, model, tile,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Write out the ice model variables.
!-----------------------------------------------------------------------
!
      CALL ice_wrt_station_nf90 (ng, model, tile,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
      CALL sediment_wrt_station_nf90 (ng, model, tile,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                Sout, STA)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Write out the Waves Effect on Currents station variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_station_nf90 (ng, model, tile,                       &
    &                            LBi, UBi, LBj, UBj,                    &
    &                            Sout, STA)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if !(defined OUT_NETCDF4 && defined DELAYED_SYNC_NF90)
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, STA(ng)%name, STA(ng)%ncid)
# endif
!
      RETURN
      END SUBROUTINE wrt_station_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_station_pio (ng, model, tile,                      &
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB, NposR, NposW
      integer :: Fcount, i, ifield, k, np, status
!
      real(dp) :: scale

      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
#  ifdef SOLVE3D
#   ifdef SEDIMENT
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
#   endif
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
#  endif

      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
#  ifdef SOLVE3D
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
#  endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      STA(ng)%Rindex=STA(ng)%Rindex+1
      Fcount=STA(ng)%Fcount
      STA(ng)%Nrec(Fcount)=STA(ng)%Nrec(Fcount)+1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
#  ifdef SOLVE3D
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
#   ifdef SEDIMENT
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
#   endif
#  endif
      END DO
!
!  Write out model time (s).
!
      CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,                &
     &                          TRIM(Vname(1,idtime)), time(ng:),       &
     &                          (/STA(ng)%Rindex/), (/1/),              &
     &                          pioFile = STA(ng)%pioFile,              &
     &                          pioVar = STA(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idFsur, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idFsur)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idFsur)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUbar, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbar)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbar)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVbar, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbar)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbar)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
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
        CALL extract_sta2d (ng, model, Cgrid, idu2dE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idu2dE)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idu2dE)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idv2dN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idv2dN)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idv2dN)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF

#  ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idUvel, u3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvel, v3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
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
        CALL extract_sta3d (ng, model, Cgrid, idu3dE, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idu3dE)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idu3dE)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta3d (ng, model, Cgrid, idv3dN, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idv3dN)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idv3dN)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out "true" vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idWvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idWvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idWvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idOvel, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idOvel)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idOvel)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, ifield, r3dvar,         &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,            &
     &                              TRIM(Vname(1,idTvar(i))), rsta,     &
     &                              (/1,1,STA(ng)%Rindex/),             &
     &                              (/N(ng),Nstation(ng),1/),           &
     &                              pioFile = STA(ng)%pioFile,          &
     &                              pioVar = STA(ng)%pioTrc(i)%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idDano, r3dvar,           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idDano)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng),Nstation(ng),1/),             &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idDano)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHsbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idHsbl)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idHsbl)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idHbbl, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idHbbl)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idHbbl)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idVvis, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVvis)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVvis)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idTdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idTdif)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idTdif)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idSdif, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idSdif)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idSdif)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtke, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idMtke)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idMtke)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, model, Cgrid, idMtls, w3dvar,           &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idMtls)), rsta,          &
     &                            (/1,1,STA(ng)%Rindex/),               &
     &                            (/N(ng)+1,Nstation(ng),1/),           &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idMtls)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Sout(idPair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idPair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Pair,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idPair)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idPair)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Sout(idUair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Uwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUair)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUair)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      IF (Sout(idVair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVair, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Vwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVair)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVair)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward surface winds (m/s) at
!  RHO-points
!
      IF (Sout(idUaiE,ng).and.Sout(idVaiN,ng)) THEN
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
        CALL extract_sta2d (ng, model, Cgrid, idUaiE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUaiE)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUaiE)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        CALL extract_sta2d (ng, model, Cgrid, idVaiN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVaiN)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVaiN)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
#   endif
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,ifield)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(ifield)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   ifdef SALINITY
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, ifield, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,ifield)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(ifield)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif

#   ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Sout(idLhea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLhea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idLhea)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idLhea)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out sensible heat flux.
!
      IF (Sout(idShea,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idShea, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idShea)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idShea)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Sout(idLrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idLrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idLrad)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idLrad)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#   endif
#   ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Sout(idSrad,ng)) THEN
        scale=rho0*Cp
        CALL extract_sta2d (ng, model, Cgrid, idSrad, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idSrad)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idSrad)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#   if defined EMINUSP && defined BULK_FLUXES
!
!  Write out E-P (m/s).
!
      IF (Sout(idEmPf,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idEmPf, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflux(:,:,isalt),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idEmPf)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idEmPf)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idevap, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idevap)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idevap)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idrain, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idrain)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idrain)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#   endif
#  endif
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idUsms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUsms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUsms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, model, Cgrid, idVsms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVsms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVsms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idUbms, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUbms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUbms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, model, Cgrid, idVbms, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVbms)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVbms)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#  ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      IF (Sout(idRwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idRwet, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%rmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idRwet)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idRwet)%vd)
      END IF
!
!  Write out wet/dry mask at U-points.
!
      IF (Sout(idUwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idUwet, u2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%umask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idUwet)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idUwet)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      IF (Sout(idVwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, model, Cgrid, idVwet, v2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%vmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL pio_netcdf_put_fvar (ng, model, STA(ng)%name,              &
     &                            TRIM(Vname(1,idVwet)), psta,          &
     &                            (/1,STA(ng)%Rindex/),                 &
     &                            (/Nstation(ng),1/),                   &
     &                            pioFile = STA(ng)%pioFile,            &
     &                            pioVar = STA(ng)%pioVar(idVwet)%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  if defined BBL_MODEL || defined WAVES_OUTPUT
!
!-----------------------------------------------------------------------
!  Write out the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
      CALL bbl_wrt_station_pio (ng, model, tile,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Write out the ice model variables.
!-----------------------------------------------------------------------
!
      CALL ice_wrt_station_pio (ng, model, tile,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Write out the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_wrt_station_pio (ng, model, tile,                 &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 Sout, STA)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WEC || defined WEC_VF
!
!-----------------------------------------------------------------------
!  Write out the Waves Effect on Currents station variables.
!-----------------------------------------------------------------------
!
      CALL wec_wrt_station_pio (ng, model, tile,                        &
    &                           LBi, UBi, LBj, UBj,                     &
    &                           Sout, STA)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifndef DELAYED_SYNC_PIO
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, STA(ng)%name, STA(ng)%pioFile)
#  endif
!
      RETURN
      END SUBROUTINE wrt_station_pio
# endif
#endif
      END MODULE wrt_station_mod

