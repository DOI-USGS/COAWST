#include "cppdefs.h"
#if (defined TANGENT || defined TL_IOMS) && defined FOUR_DVAR
      SUBROUTINE rp_wrt_ini (ng, Tindex, OutRec)
!
!svn $Id: rp_wrt_ini.F 677 2008-08-05 20:17:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine writes representer model initial conditions into       !
!  IRPname initial conditions NetCDF file.                             !
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
# if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
      USE mod_forces
# endif
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
# if defined SEDIMENT_NOT_YET || defined BBL_MODEL_NOT_YET
      USE mod_sediment
# endif
      USE mod_stepping

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Tindex, OutRec
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status, varid

      integer :: nf_fwrite2d
# ifdef SOLVE3D
      integer :: nf_fwrite3d
# endif

      real(r8) :: scale

      character (len=15) :: string
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out tangent linear initial conditions.
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
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iRPM, IRPname(ng),                      &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/OutRec/), (/1/),                          &
     &                      ncid = ncIRPid(ng),                         &
     &                      varid = irpVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iRPM, ncIRPid(ng), irpVid(idFsur,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % tl_zeta(LBi,LBj,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), OutRec
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
      status=nf_fwrite2d(ng, iRPM, ncIRPid(ng), irpVid(idUbar,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % tl_ubar(LBi,LBj,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)), OutRec
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
      status=nf_fwrite2d(ng, iRPM, ncIRPid(ng), irpVid(idVbar,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
# endif
     &                   OCEAN(ng) % tl_vbar(LBi,LBj,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
!!    scale=rho0
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iRPM, ncIRPid(ng), irpVid(idUsms,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % tl_ustr(LBi,LBj,1,Tindex))
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
!!    scale=rho0
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iRPM, ncIRPid(ng), irpVid(idVsms,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   FORCES(ng) % tl_vstr(LBi,LBj,1,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVsms)), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iRPM, ncIRPid(ng), irpVid(idUvel,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % tl_u(LBi,LBj,1,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), OutRec
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
      status=nf_fwrite3d(ng, iRPM, ncIRPid(ng), irpVid(idVvel,ng),      &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask(LBi,LBj),                     &
#  endif
     &                   OCEAN(ng) % tl_v(LBi,LBj,1,Tindex))
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), OutRec
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
        status=nf_fwrite3d(ng, iRPM, ncIRPid(ng), irpTid(itrc,ng),      &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#  endif
     &                     OCEAN(ng) % tl_t(LBi,LBj,1,Tindex,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  ifdef ADJUST_STFLUX
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
!!      scale=rho0*Cp
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iRPM, ncIRPid(ng),                       &
     &                     irpVid(idTsur(itrc),ng),                     &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#   ifdef MASKING
     &                     GRID(ng) % rmask(LBi,LBj),                   &
#   endif
     &                     FORCES(ng) % tl_tflux(LBi,LBj,1,Tindex,itrc))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),              &
     &                        OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent linear initial NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncIRPid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(IRPname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag, 1)
# endif

      IF (Master) THEN
        IF (OutRec.eq.2) THEN
          string='initial  fields'
        END IF
# ifdef SOLVE3D
        WRITE (stdout,30) string, Nrun, Tindex, Tindex, OutRec
# else
        WRITE (stdout,30) string, Nrun, Tindex, OutRec
# endif
      END IF
!
  10  FORMAT (/,' RP_WRT_INI - error while writing variable: ',a,/,14x,    &
     &        'into tangent initial NetCDF file for time record: ',i4)
  20  FORMAT (/,' RP_WRT_INI - unable to synchronize to disk ',            &
     &        'representers model initial file:',/,14x,a)
# ifdef SOLVE3D
  30  FORMAT (3x,'RP_WRT_INI   - wrote ',a,' (Iter=',i4.4,', Index=',i1,   &
     &        ',',i1,', Rec=',i4.4,')')
# else
  30  FORMAT (3x,'RP_WRT_INI   - wrote ',a,' (Iter=',i4.4,', Index=',i1,   &
     &        ', Rec=',i4.4,')')
# endif
      RETURN
      END SUBROUTINE rp_wrt_ini
#else
      SUBROUTINE rp_wrt_ini
      RETURN
      END SUBROUTINE rp_wrt_ini
#endif
