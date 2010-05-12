#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE wrt_ini (ng, Tindex)
!
!svn $Id: wrt_ini.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes state variables initial conditions into initial !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Tindex     State variables time index to write.                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
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
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Tindex
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status, varid

      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile='wrt_ini.F'
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
      tINIindx(ng)=tINIindx(ng)+1
      NrecINI(ng)=NrecINI(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, INIname(ng),                      &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tINIindx(ng)/), (/1/),                    &
     &                      ncid = ncINIid(ng),                         &
     &                      varid = iniVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncINIid(ng), iniVid(idFsur,ng),      &
     &                   tINIindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
# ifdef WET_DRY
     &                   OCEAN(ng) % zeta(:,:,Tindex),                  &
     &                   SetFillVal = .FALSE.)
# else
     &                   OCEAN(ng) % zeta(:,:,Tindex))
# endif
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), tINIindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncINIid(ng), iniVid(idUbar,ng),      &
     &                   tINIindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask,                              &
# endif
     &                   OCEAN(ng) % ubar(:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), tINIindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncINIid(ng), iniVid(idVbar,ng),      &
     &                   tINIindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask,                              &
# endif
     &                   OCEAN(ng) % vbar(:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), tINIindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idUvel,ng),      &
     &                   tINIindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask,                              &
#  endif
     &                   OCEAN(ng) % u(:,:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), tINIindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idVvel,ng),      &
     &                   tINIindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   OCEAN(ng) % v(:,:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), tINIindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniTid(itrc,ng),      &
     &                     tINIindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % t(:,:,:,Tindex,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tINIindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO

#  if defined BVF_MIXING  || defined GLS_MIXING || \
      defined MY25_MIXING || defined LMD_MIXING
!
!  If defined, write out vertical viscosity coefficient.
!
      IF (iniVid(idVvis,ng).gt.0) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idVvis,ng),    &
     &                     tINIindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akv)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tINIindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  If defined, write out vertical diffusion coefficient for potential
!  temperature.
!
      IF (iniVid(idTdif,ng).gt.0) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idTdif,ng),    &
     &                     tINIindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#   endif
     &                     MIXING(ng) % Akt(:,:,:,itemp))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tINIindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SALINITY
!
!  If defined, write out vertical diffusion coefficient for salinity.
!
      IF (iniVid(idSdif,ng).gt.0) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idSdif,ng),    &
     &                     tINIindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#    endif
     &                     MIXING(ng) % Akt(:,:,:,isalt))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tINIindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, INIname(ng), ncINIid(ng))
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) THEN
# ifdef SOLVE3D
        WRITE (stdout,20) Nrun, Tindex, Tindex, tINIindx(ng)
# else
        WRITE (stdout,20) Nrun, Tindex, tINIindx(ng)
# endif
      END IF
!
  10  FORMAT (/,' WRT_INI - error while writing variable: ',a,/,11x,    &
     &        'into initial NetCDF file for time record: ',i4)
# ifdef SOLVE3D
  20  FORMAT (6x,'WRT_INI   - wrote initial  fields (Iter=',i4.4,       &
     &           ', Index=',i1,',',i1,', Rec=',i4.4,')')
# else
  20  FORMAT (6x,'WRT_INI   - wrote initial  fields (Iter=',i4.4,       &
     &           ', Index=',i1,', Rec=',i4.4,')')
# endif
      RETURN
      END SUBROUTINE wrt_ini

# if defined ADJUST_BOUNDARY || \
     defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
      SUBROUTINE wrt_frc (ng, Tindex, OutRec)
!
!=======================================================================
!                                                                      !
!  This routine writes surface forcing and/or open boundary background !
!  state into nonlinear model initial conditions NetCDF file.          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Tindex     State variables time index to write.                  !
!     OutRec     NetCDF file unlimited dimension record to write.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#  if defined ADJUST_BOUNDARY
      USE mod_boundary
#  endif
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
      USE mod_grid
#  endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
#  ifdef SOLVE3D
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
#  endif
# endif
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Tindex, OutRec
!
!  Local variable declarations.
!
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      integer :: LBi, UBi, LBj, UBj
#  endif
#  ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#  endif

      integer :: gfactor, gtype, i, itrc, status

      real(r8) :: scale
!
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
#  endif
#  ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#  endif
!
      SourceFile='wrt_ini.F, wrt_frc'
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
      gfactor=-1
#  else
      gfactor=1
#  endif

#  ifdef ADJUST_BOUNDARY
!
!  Write out open boundary fields. Notice that these fields have their
!  own fixed time-dimension (of size Nbrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          iniVid(idSbry(isFsur),ng),              &
     &                          OutRec, r2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          iniVid(idSbry(isUbar),ng),              &
     &                          OutRec, u2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          iniVid(idSbry(isVbar),ng),              &
     &                          OutRec, v2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SOLVE3D
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          iniVid(idSbry(isUvel),ng),              &
     &                          OutRec, u3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          iniVid(idSbry(isVvel),ng),              &
     &                          OutRec, v3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_r8
          status=nf_fwrite3d_bry (ng, iNLM, INIname(ng), ncINIid(ng),   &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            iniVid(idSbry(isTvar(itrc)),ng),      &
     &                            OutRec, r3dvar,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Tindex,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          OutRec
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
      scale=rho0                            ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idUsms,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ustr(:,:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=rho0                            ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idVsms,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % vstr(:,:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          IF (itrc.eq.itemp) THEN
            scale=rho0*Cp                   ! Celsius m/s to W/m2
          ELSE
            scale=1.0_r8
          END IF
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncINIid(ng),                     &
     &                       iniVid(idTsur(itrc),ng),                   &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tflux(:,:,:,Tindex,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, INIname(ng), ncINIid(ng))
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) THEN
        WRITE (stdout,20) Nrun, Tindex, OutRec
      END IF
!
  10  FORMAT (/,' WRT_FRC - error while writing variable: ',a,/,11x,    &
     &        'into initial NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_FRC   - wrote initial  fields (Iter=',i4.4,       &
     &           ', Index=',i1,', Rec=',i4.4,')')

      RETURN
      END SUBROUTINE wrt_frc

      SUBROUTINE wrt_frc_AD (ng, Tindex, OutRec)
!
!=======================================================================
!                                                                      !
!  This routine writes surface forcing and/or open boundary fields     !
!  into initial conditions NetCDF file.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     Tindex     State variables time index to write.                  !
!     OutRec     NetCDF file unlimited dimension record to write.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#  if defined ADJUST_BOUNDARY
      USE mod_boundary
#  endif
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
      USE mod_grid
#  endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
#  ifdef SOLVE3D
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
#  endif
# endif
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Tindex, OutRec
!
!  Local variable declarations.
!
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      integer :: LBi, UBi, LBj, UBj
#  endif
#  ifdef ADJUST_BOUNDARY
      integer :: IorJ, LBij, UBij
#  endif
      integer :: gfactor, gtype, i, itrc, status

      real(r8) :: scale
!
#  if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
#  endif
#  ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
      IorJ=IOBOUNDS(ng)%IorJ
#  endif
!
      SourceFile='wrt_ini.F, wrt_frc_AD'
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
      gfactor=-1
#  else
      gfactor=1
#  endif

#  ifdef ADJUST_BOUNDARY
!
!  Write out open boundary fields. Notice that these fields have their
!  own fixed time-dimension (of size Nbrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          iniVid(idSbry(isFsur),ng),              &
     &                          OutRec, r2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,   &
     &                                                     Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          iniVid(idSbry(isUbar),ng),              &
     &                          OutRec, u2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,   &
     &                                                     Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          iniVid(idSbry(isVbar),ng),              &
     &                          OutRec, v2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,   &
     &                                                     Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#   ifdef SOLVE3D
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          iniVid(idSbry(isUvel),ng),              &
     &                          OutRec, u3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,    &
     &                                                  Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, INIname(ng), ncINIid(ng),     &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          iniVid(idSbry(isVvel),ng),              &
     &                          OutRec, v3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,    &
     &                                                  Tindex))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_r8
          status=nf_fwrite3d_bry (ng, iNLM, INIname(ng), ncINIid(ng),   &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            iniVid(idSbry(isTvar(itrc)),ng),      &
     &                            OutRec, r3dvar,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,  &
     &                                                    Tindex,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          tHSSindx(ng)
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
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idUsms,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, ncINIid(ng), iniVid(idVsms,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif

#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncINIid(ng),                     &
     &                       iniVid(idTsur(itrc),ng),                   &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_tflux(:,:,:,Tindex,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize initial NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, INIname(ng), ncINIid(ng))
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) THEN
        WRITE (stdout,20) Nrun, Tindex, OutRec
      END IF
!
  10  FORMAT (/,' WRT_FRC_AD - error while writing variable: ',a,/,11x, &
     &        'into initial NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_FRC_AD   - wrote initial  fields (Iter=',i4.4,    &
     &           ', Index=',i1,', Rec=',i4.4,')')

      RETURN
      END SUBROUTINE wrt_frc_AD
# endif
#else
      SUBROUTINE wrt_ini
      RETURN
      END SUBROUTINE wrt_ini
#endif
