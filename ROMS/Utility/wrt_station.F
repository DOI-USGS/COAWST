#include "cppdefs.h"
#ifdef STATIONS
      SUBROUTINE wrt_station (ng)
!
!svn $Id: wrt_station.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
      USE mod_sediment
# endif
      USE mod_stepping
!
# ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta2d, extract_sta3d
# else
      USE extract_sta_mod, ONLY : extract_sta2d
# endif
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

      integer :: NposB, NposR, NposW, LBi, UBi, LBj, UBj
      integer :: i, ifield, k, np, status, tile

      integer, dimension(3) :: start, total

      real(r8) :: scale

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

      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tSTAindx(ng)=tSTAindx(ng)+1
      NrecSTA(ng)=NrecSTA(ng)+1
!
!  Specify parameters for writing profiles to NetCDF files. These
!  are overridden for some variables below
!
      start(1)=1
      total(1)=N(ng)
      start(2)=1
      total(2)=Nstation(ng)
      start(3)=tSTAindx(ng)
      total(3)=1
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
      IF (Master) THEN
        status=nf90_put_var(ncSTAid(ng), staVid(idtime,ng),             &
     &                      time(ng:), start(3:), total(3:))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idtime)), tSTAindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idFsur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idFsur,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbar, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbar,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbar, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,KOUT),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbar,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idUvel, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUvel,ng),           &
     &                        rsta, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvel, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,NOUT),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVvel,ng),           &
     &                        rsta, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idWvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idWvel,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idOvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idOvel,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, ifield, r3dvar,          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,NOUT,i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          IF (Master) THEN
            status=nf90_put_var(ncSTAid(ng), staTid(i,ng),              &
     &                          rsta, start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(i))), tSTAindx(ng)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idDano, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idDano,ng),           &
     &                        rsta, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Sout(idHsbl,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHsbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hsbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idHsbl,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Sout(idHbbl,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idHbbl, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng)%hbbl,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idHbbl,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvis, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVvis,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idTdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idTdif,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Sout(idSdif,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idSdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,isalt),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idSdif,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#  endif
#  if defined MY25_MIXING || defined GLS_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtke, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idMtke,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_r8
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtls, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,NOUT),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idMtls,ng),           &
     &                        rsta, start, total+(/1,0,0/))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(ifield,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out surface salt flux.
!
      IF (Sout(idTsur(isalt),ng)) THEN
        ifield=idTsur(isalt)
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,isalt),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(ifield,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,ifield)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idLhea,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idShea,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idLrad,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#   ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Sout(idevap,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idevap, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%evap,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idevap,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Sout(idrain,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idrain, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%rain,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idrain,ng),       &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#   endif
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idSrad,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUsms,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVsms,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbms,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbms,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress.
!
      IF (Sout(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbrs,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbrs,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbws,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbws,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbcs,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
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
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbcs,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbot,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbot,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (Sout(idUbur,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idUbur,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (Sout(idVbvr,ng)) THEN
        scale=1.0_r8
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbvr, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        IF (Master) THEN
          status=nf90_put_var(ncSTAid(ng), staVid(idVbvr,ng),           &
     &                        psta, start(2:), total(2:))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), tSTAindx(ng)
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END IF
# endif
# ifdef SEDIMENT
!
!  Write out sediment fraction of each size class in each bed layer.
!
      total(1)=Nbed                       ! Warning, total(1) redefined
      DO i=1,NST      
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, idfrac(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, OCEAN(ng)%bed_frac(:,:,:,i),       &
     &                        NposB, XposB, YposB, ZposB, bsta)
          IF (Master) THEN
            status=nf90_put_var(ncSTAid(ng), staVid(idfrac(i),ng),      &
     &                          bsta, start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tSTAindx(ng)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, idBmas(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        OCEAN(ng)%bed_mass(:,:,:,NOUT,i),         &
     &                        NposB, XposB, YposB, ZposB, bsta)
          IF (Master) THEN
            status=nf90_put_var(ncSTAid(ng), staVid(idBmas(i),ng),      &
     &                          bsta, start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tSTAindx(ng)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP      
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta3d (ng, iNLM, Cgrid, idSbed(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, OCEAN(ng)%bed(:,:,:,i),            &
     &                        NposB, XposB, YposB, ZposB, bsta)
          IF (Master) THEN
            status=nf90_put_var(ncSTAid(ng), staVid(idSbed(i),ng),      &
     &                          bsta, start, total)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tSTAindx(ng)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP      
        IF (Sout(idBott(i),ng)) THEN
          scale=1.0_r8
          CALL extract_sta2d (ng, iNLM, Cgrid, idBott(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, OCEAN(ng)%bottom(:,:,i),           &
     &                        Nstation(ng), Xpos, Ypos, psta)
          IF (Master) THEN
            status=nf90_put_var(ncSTAid(ng), staVid(idBott(i),ng),      &
     &                          bsta, start(2:), total(2:))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tSTAindx(ng)
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          END IF
        END IF
      END DO
# endif
#endif
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        status=nf90_sync(ncSTAid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
  10  FORMAT (/,' WRT_STATION - error while writing variable: ',a,/,    &
     &        15x,'into stations NetCDF file for time record: ',i4)
  20  FORMAT (/,' WRT_STATION - unable to synchronize stations',        &
     &        1x,'NetCDF file to disk.')
#else
      SUBROUTINE wrt_station
#endif
      RETURN
      END SUBROUTINE wrt_station
