#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE wrt_station (ng)
!
!svn $Id: wrt_station.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out data into stations NetCDF file.             !
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
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
# endif
      USE mod_stepping
!
# ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta2d, extract_sta3d
# else
      USE extract_sta_mod, ONLY : extract_sta2d
# endif
      USE uv_rotate_mod,   ONLY : uv_rotate2d
# ifdef SOLVE3D
      USE uv_rotate_mod,   ONLY : uv_rotate3d
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB, NposR, NposW, LBi, UBi, LBj, UBj
      integer :: Fcount, i, ifield, k, np, status, tile
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
     &  __FILE__
!
      SourceFile=MyFile
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
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
      CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                     &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idFsur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idFsur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idFsur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# if defined SEDIMENT && defined SED_MORPH
!
!  Define time-varying bathymetry.
!
      IF (Sout(idbath,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idbath, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%h,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idbath)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idbath))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbar, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbar, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
# ifdef DISTRIBUTE
        tile=MyRank
# else
        tile=-1
# endif
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

        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idu2dE, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idu2dE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu2dE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL extract_sta2d (ng, iNLM, Cgrid, idv2dN, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idUvel, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvel, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
#  ifdef DISTRIBUTE
        tile=MyRank
#  else
        tile=-1
#  endif
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

        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idu3dE, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idu3dE)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu3dE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL extract_sta3d (ng, iNLM, Cgrid, idv3dN, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
!  Write out vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idWvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idOvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
          CALL extract_sta3d (ng, iNLM, Cgrid, ifield, r3dvar,          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idDano, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idHsbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idHbbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvis, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idTdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idSdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idSdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idSdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
#  if defined MY25_MIXING || defined GLS_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtke, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtls, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idPair, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Pair,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idUair, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Uwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUair)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUair))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

      IF (Sout(idVair,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idVair, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%Vwind,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
      IF (Sout(idUairE,ng).and.Sout(idVairN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
# ifdef DISTRIBUTE
        tile=MyRank
# else
        tile=-1
# endif
        CALL uv_rotate2d (ng, TILE, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
# ifdef MASKING
     &                    GRID(ng) % rmask_full,                        &
# endif
     &                    FORCES(ng) % Uwind,                           &
     &                    FORCES(ng) % Vwind,                           &
     &                    Ur2d, Vr2d)

        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUairE, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUairE)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUairE))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idVairN, r2dvar,           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVairN)), psta,             &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVairN))
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
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idLhea, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lhflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idShea, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%shflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idLrad, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%lrflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idSrad, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%srflx,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idEmPf, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflux(:,:,isalt),        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idEmPf)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idEmPf))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idevap, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idrain, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idUsms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idVsms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbms))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef SOLVE3D
#  ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (Sout(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbrs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbrs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (Sout(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbrs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbrs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (Sout(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbws, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbws)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbws))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Sout(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbws, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbws)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbws))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (Sout(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbcs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbcs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbcs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (Sout(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbcs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbcs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbcs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbot)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbot))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbot)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbot))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (Sout(idUbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbur))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (Sout(idVbvr,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbvr, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbvr)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbvr))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif

#  ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, idfrac(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed_frac(:,:,:,i),      &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idfrac(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idfrac(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, idBmas(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        SEDBED(ng)%bed_mass(:,:,:,NOUT,i),        &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idBmas(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idBmas(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, idSbed(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed(:,:,:,i),           &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idSbed(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idSbed(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif

#  if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (Sout(idBott(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, iNLM, Cgrid, idBott(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBED(ng)%bottom(:,:,i),          &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idBott(i))), rsta,         &
     &                          (/1,STA(ng)%Rindex/),                   &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idBott(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
#  endif
# endif
# ifdef WEC_MELLOR
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Sout(idW2xx,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idW2xx, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % Sxx_bar,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW2xx)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW2xx))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Sout(idW2xy,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idW2xy, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % Sxy_bar,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW2xy)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW2xy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Sout(idW2yy,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idW2yy, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % Syy_bar,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW2yy)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW2yy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out total 2D U-radiation stress.
!
      IF (Sout(idU2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idU2rs, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rustr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU2rs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU2rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out total 2D V-radiation stress.
!
      IF (Sout(idV2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idV2rs, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rvstr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV2rs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV2rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  ifdef SOLVE3D
!
!  Write out 3D radiation stress, Sxx-horizontal component.
!
      IF (Sout(idW3xx,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3xx, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % Sxx,                    &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3xx)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3xx))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D radiation stress, Sxy-horizontal component.
!
      IF (Sout(idW3xy,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3xy, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % Sxy,                    &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3xy)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3xy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D radiation stress, Syy-horizontal component.
!
      IF (Sout(idW3yy,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3yy, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % Syy,                    &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3yy)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3yy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D radiation stress, Szx-vertical component.
!
      IF (Sout(idW3zx,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3zx, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3zx)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3zx))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D radiation stress, Szy-vertical component.
!
      IF (Sout(idW3zy,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3zy, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % Szx,                    &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3zy)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3zy))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
# endif
# ifdef WEC
!
!  Write out 2D U-momentum Stokes drift velocity.
!
      IF (Sout(idU2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idU2Sd, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % ubar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU2Sd)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU2Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 2D V-momentum Stokes drift velocity.
!
      IF (Sout(idV2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idV2Sd, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % vbar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV2Sd)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV2Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out total 2D WEC U-stress.
!
      IF (Sout(idU2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idU2rs, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rustr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU2rs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU2rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out total 2D WEC V-stress.
!
      IF (Sout(idV2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idV2rs, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rvstr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV2rs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV2rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  ifdef SOLVE3D
!
!  Write out 3D U-momentum Stokes drift velocity.
!
      IF (Sout(idU3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idU3Sd, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % u_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU3Sd)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU3Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D V-momentum stokes velocity.
!
      IF (Sout(idV3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idV3Sd, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % v_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV3Sd)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV3Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D Omega-momentum stokes velocity.
!
      IF (Sout(idW3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3Sd, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % W_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3Sd)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3Sd))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D W-momentum stokes velocity
!
      IF (Sout(idW3St,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3St, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % wstvel,                  &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3St)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3St))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D total WEC U-stress.
!
      IF (Sout(idU3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, iNLM, Cgrid, idU3rs, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rustr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU3rs)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU3rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out 3D total WEC V-stress.
!
      IF (Sout(idV3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, iNLM, Cgrid, idV3rs, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rvstr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV3rs)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV3rs))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#  endif
# endif
# ifdef WAVES_HEIGHT
!
!  Write out wind-induced wave height.
!
      IF (Sout(idWamp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWamp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Hwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWamp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWamp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_LENGTH
!
!  Write out wind-induced wave length.
!
      IF (Sout(idWlen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWlen, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWlen)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWlen))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_LENGTHP
!
!  Write out wind-induced peak wave length.
!
      IF (Sout(idWlep,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWlep, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwavep,                 &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWlep)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWlep))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_DIR
!
!  Write out wind-induced wave direction.
!
      IF (Sout(idWdir,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdir, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdir)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdir))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_DIRP
!
!  Write out wind-induced peak wave direction.
!
      IF (Sout(idWdip,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdip, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwavep,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdip)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdip))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Sout(idWptp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWptp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_top,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWptp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWptp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (Sout(idWpbt,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWpbt, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_bot,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWpbt)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWpbt))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
#if defined BBL_MODEL
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (Sout(idWorb,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWorb, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Uwave_rms,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWorb)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWorb))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
!
!  Write out wave dissipation due to bottom friction.
!
      IF (Sout(idWdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdif, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_fric,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdif)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdif))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR
!
!  Write out wave dissipation due to breaking.
!
      IF (Sout(idWdib,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdib, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_break,           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdib)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdib))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to white capping.
!
      IF (Sout(idWdiw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdiw, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_wcap,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdiw)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdiw))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WEC_ROLLER
!
!  Write out wave roller dissipation.
!
      IF (Sout(idWdis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdis, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_roller,          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdis)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdis))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined ROLLER_RENIERS
!
!  Write out wave roller action density.
!
      IF (Sout(idWrol,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWrol, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % rollA,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWrol)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWrol))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined ROLLER_SVENDSEN
!
!  Write out percent wave breaking.
!
      IF (Sout(idWbrk,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWbrk, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_break,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWbrk)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWbrk))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WAVES_DSPR
!
!  Write out waves directional spread.
!
      IF (Sout(idWvds,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWvds, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_ds,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWvds)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvds))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
      IF (Sout(idWvqp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWvqp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Wave_qp,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWvqp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvqp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# if defined WEC_VF
!
!  Write out WEC quasi-static sea level adjustment.
!
      IF (Sout(idWztw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWztw, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % zeta(:,:,KOUT),          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWztw)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWztw))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out WEC quasi-static pressure.
!
      IF (Sout(idWqsp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWqsp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % qsp,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWqsp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWqsp))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
!  Write out WEC Bernoulli head.
!
      IF (Sout(idWbeh,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWbeh, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % bh,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWbeh)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWbeh))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
# ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      IF (Sout(idRwet,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idRwet, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%rmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idUwet, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%umask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
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
        CALL extract_sta2d (ng, iNLM, Cgrid, idVwet, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%vmask_wet,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVwet)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVwet))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, STA(ng)%name, STA(ng)%ncid)

#else
      SUBROUTINE wrt_station
#endif
      RETURN
      END SUBROUTINE wrt_station
