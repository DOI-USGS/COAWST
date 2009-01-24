#include "cppdefs.h"
#if (defined TANGENT || defined TL_IOMS) && defined FOUR_DVAR
      SUBROUTINE tl_wrt_ini (ng, Tindex, OutRec)
!
!svn $Id: tl_wrt_ini.F 677 2008-08-05 20:17:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!======================================================================= 
!                                                                      !
!  This routine writes tangent linear initial conditions into tangent  !
!  linear initial NetCDF file.                                         !
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
      CALL netcdf_put_fvar (ng, iTLM, ITLname(ng),                      &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/OutRec/), (/1/),                          &
     &                      ncid = ncITLid(ng),                         &
     &                      varid = itlVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iTLM, ncITLid(ng), itlVid(idFsur,ng),      &
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
      status=nf_fwrite2d(ng, iTLM, ncITLid(ng), itlVid(idUbar,ng),      &
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
      status=nf_fwrite2d(ng, iTLM, ncITLid(ng), itlVid(idVbar,ng),      &
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
      scale=rho0                            ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iTLM, ncITLid(ng), itlVid(idUsms,ng),      &
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
      scale=rho0                            ! m2/s2 to N/m2 (Pa)
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iTLM, ncITLid(ng), itlVid(idVsms,ng),      &
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
      status=nf_fwrite3d(ng, iTLM, ncITLid(ng), itlVid(idUvel,ng),      &
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
      status=nf_fwrite3d(ng, iTLM, ncITLid(ng), itlVid(idVvel,ng),      &
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
        status=nf_fwrite3d(ng, iTLM, ncITLid(ng), itlTid(itrc,ng),      &
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
        IF (itrc.eq.itemp) THEN
          scale=rho0*Cp                     ! Celsius m/s to W/m2
        ELSE
          scale=1.0_r8
        END IF
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iTLM, ncITLid(ng),                       &
     &                     itlVid(idTsur(itrc),ng),                     &
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

# if defined IS4DVAR || defined IS4DVAR_OLD || defined BACKGROUND
!
!-----------------------------------------------------------------------
!  Write other 4DVAR data to MODname file.
!-----------------------------------------------------------------------

#  if defined IS4DVAR || defined IS4DVAR_OLD
!
!  Write out cost function gradient norm. Notice that it is written
!  into MODname file instead of ITLname file.
!
#   ifdef IS4DVAR_OLD
      IF (Ipass.eq.1) THEN
#   endif
        CALL netcdf_put_fvar (ng, iTLM, MODname(ng), 'cost_gradient',   &
     &                        FOURDVAR(ng)%CostGrad(0:),                &
     &                        (/1,Nrun/), (/NstateVar(ng)+1,1/),        &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
#   ifdef IS4DVAR_OLD
      END IF
#   endif
!
!  Write out model minus observations misfit cost function value to
!  both nonlinear and tangent linear initial file.  Notice that it is
!  written into MODname file instead of ITLname file.
!
#   ifdef IS4DVAR_OLD
      IF (Ipass.eq.1) THEN
#   endif
        CALL netcdf_put_fvar (ng, iTLM, MODname(ng), 'cost_function',   &
     &                        FOURDVAR(ng)%ObsCost(0:),                 &
     &                        (/1,Nrun/), (/NstateVar(ng)+1,1/),        &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
#   ifdef IS4DVAR_OLD
      END IF
#   endif
#  endif

#  ifdef BACKGROUND
!
!  Write out model minus background misfit cost function value to
!  both nonlinear and tangent linear initial file. Notice that it
!  is written into MODname file instead of ITLname file.
!
#   ifdef IS4DVAR_OLD
      IF (Ipass.eq.1) THEN
#   endif
        CALL netcdf_put_fvar (ng, iTLM, MODname(ng), 'back_function',   &
     &                        FOURDVAR(ng)%BackCost(0:),                &
     &                        (/1,Nrun/), (/NstateVar(ng)+1,1/),        &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
#   ifdef IS4DVAR_OLD
      END IF
#   endif
#  endif

#  if defined IS4DVAR || defined IS4DVAR_OLD
!
!  Write out current optimal, normalized cost function minimun.  Notice
!  that it is written into MODname file instead of ITLname file.
!
#   ifdef IS4DVAR_OLD
      IF (Ipass.eq.1) THEN
#   endif
        CALL netcdf_put_fvar (ng, iTLM, MODname(ng), 'Jmin',            &
     &                        Optimality(ng:), (/Nrun/), (/1/),         &
     &                        ncid = ncMODid(ng))
        IF (exit_flag.ne.NoError) RETURN
#   ifdef IS4DVAR_OLD
      END IF
#   endif
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent linear initial NetCDF file to disk to allow other
!  processes to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncITLid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'tangent linear initial', TRIM(ITLname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag, 1)
# endif
      IF (exit_flag.ne.NoError) RETURN

# if defined IS4DVAR || defined IS4DVAR_OLD || defined BACKGROUND
      IF (OutThread) THEN
        status=nf90_sync(ncMODid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) 'model/observation', TRIM(MODname(ng))
          exit_flag=3
          ioerror=status
        END IF
      END IF
#  ifdef DISTRIBUTE
      CALL mp_bcasti (ng, iTLM, exit_flag, 1)
#  endif
# endif
!
      IF (Master) THEN
        IF (OutRec.eq.1) THEN
          string='initial  fields'
        ELSE IF (OutRec.eq.2) THEN
          string='v-increments   '
        ELSE IF (OutRec.eq.3) THEN
          string='v-increments   '
        END IF
# ifdef SOLVE3D
        WRITE (stdout,30) string, Nrun, Tindex, Tindex, OutRec
# else
        WRITE (stdout,30) string, Nrun, Tindex, OutRec
# endif
      END IF
!
  10  FORMAT (/,' TL_WRT_INI - error while writing variable: ',a,/,14x,    &
     &        'into tangent initial NetCDF file for time record: ',i4)
  20  FORMAT (/,' TL_WRT_INI - unable to synchronize to disk ',a,          &
     &        ' file:',/,14x,a)
# ifdef SOLVE3D
  30  FORMAT (3x,'TL_WRT_INI   - wrote ',a,' (Iter=',i4.4,', Index=',i1,   &
     &        ',',i1,', Rec=',i4.4,')')
# else
  30  FORMAT (3x,'TL_WRT_INI   - wrote ',a,' (Iter=',i4.4,', Index=',i1,   &
     &        ', Rec=',i4.4,')')
# endif
      RETURN
      END SUBROUTINE tl_wrt_ini
#else
      SUBROUTINE tl_wrt_ini
      RETURN
      END SUBROUTINE tl_wrt_ini
#endif
