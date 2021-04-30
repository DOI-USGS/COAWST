#include "cppdefs.h"
      MODULE wrt_ini_mod

#ifdef FOUR_DVAR
!
!svn $Id: wrt_ini.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
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
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PUBLIC :: wrt_ini
# if defined ADJUST_BOUNDARY || \
     defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
      PUBLIC :: wrt_frc
      PUBLIC :: wrt_frc_AD
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_ini (ng, Tindex, OutRec)
!***********************************************************************
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
      USE dateclock_mod,   ONLY : time_string
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Tindex

      integer, intent(in), optional :: OutRec
!
!  Local variable declarations.
!
      logical :: SetFillVal
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Fcount, gfactor, gtype, i, itrc, status, varid
!
      real(r8) :: Fmin, Fmax
      real(dp) :: my_time, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set SetFillVal to FALSE. This is essential to ensure identical
!  outer-loop solutions when increments are zero. NOT SURE WHY THIS
!  IS NECESSARY - PERHAPS SOMETHING WRONG WITH THE LOGIC OF
!  nf_fwrite routines.
!-----------------------------------------------------------------------
!
      SetFillVal=.FALSE.
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
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
      IF (PRESENT(OutRec)) THEN
        INI(ng)%Rindex=OutRec
      ELSE
        INI(ng)%Rindex=INI(ng)%Rindex+1
      END IF
      Fcount=INI(ng)%Fcount
      INI(ng)%Nrec(Fcount)=INI(ng)%Nrec(Fcount)+1
!
!  Write out model time (s). Use the "tdays" variable here because of
!  the management of the "time" variable due to nesting.
!
      my_time=tdays(ng)*day2sec

      CALL netcdf_put_fvar (ng, iNLM, INI(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), my_time,             &
     &                      (/INI(ng)%Rindex/), (/1/),                  &
     &                      ncid = INI(ng)%ncid,                        &
     &                      varid = INI(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        RETURN
      ELSE
        IF (Master) THEN
          CALL time_string (my_time, t_code)
          WRITE (Tstring,'(f15.4)') my_time*sec2day
          WRITE (stdout,20) t_code, ng, TRIM(ADJUSTL(Tstring)),         &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex, Tindex
        END IF
      END IF
!
!  Write out free-surface (m)
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idFsur),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   OCEAN(ng) % zeta(:,:,Tindex),                  &
# ifdef WET_DRY
     &                   SetFillVal = .FALSE.,                          &
# else
     &                   SetFillVal = SetFillVal,                       &
# endif
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsur)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUbar),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
# endif
     &                   OCEAN(ng) % ubar(:,:,Tindex),                  &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUbar)),  TRIM(INI(ng)%name), &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVbar),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
# endif
     &                   OCEAN(ng) % vbar(:,:,Tindex),                  &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVbar)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF

# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUvel),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % umask_full,                         &
#  endif
     &                   OCEAN(ng) % u(:,:,:,Tindex),                   &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUvel)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVvel),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % vmask_full,                         &
#  endif
     &                   OCEAN(ng) % v(:,:,:,Tindex),                   &
     &                   SetFillVal = SetFillVal,                       &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvel)), TRIM(INI(ng)%name),  &
     &                      INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Tid(itrc),   &
     &                     INI(ng)%Rindex, gtype,                       &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     OCEAN(ng) % t(:,:,:,Tindex,itrc),            &
     &                     SetFillVal = SetFillVal,                     &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),              &
     &                        TRIM(INI(ng)%name), INI(ng)%Rindex
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
          END IF
        END IF
      END DO

#  if defined BVF_MIXING  || defined GLS_MIXING || \
      defined MY25_MIXING || defined LMD_MIXING
!
!  If defined, write out vertical viscosity coefficient.
!
      IF (INI(ng)%Vid(idVvis).le.0) THEN
        CALL netcdf_inq_varid (ng, iNLM, INI(ng)%name, Vname(1,idVvis), &
     &                         INI(ng)%ncid, INI(ng)%Vid(idVvis))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) RETURN
      END IF
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVvis),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akv,                              &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVvis)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvis)), Fmin, Fmax
        END IF
      END IF
!
!  If defined, write out vertical diffusion coefficient for potential
!  temperature.
!
      IF (INI(ng)%Vid(idTdif).le.0) THEN
        CALL netcdf_inq_varid (ng, iNLM, INI(ng)%name, Vname(1,idTdif), &
     &                         INI(ng)%ncid, INI(ng)%Vid(idTdif))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) RETURN
      END IF
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idTdif),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   MIXING(ng) % Akt(:,:,:,itemp),                 &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTdif)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idTdif)), Fmin, Fmax
        END IF
      END IF

#   ifdef SALINITY
!
!  If defined, write out vertical diffusion coefficient for salinity.
!
      IF (INI(ng)%Vid(idSdif).le.0) THEN
        CALL netcdf_inq_varid (ng, iNLM, INI(ng)%name, Vname(1,idSdif), &
     &                         INI(ng)%ncid, INI(ng)%Vid(idSdif))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) RETURN
      END IF
      scale=1.0_dp
      gtype=gfactor*w3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idSdif),   &
     &                   INI(ng)%Rindex, gtype,                         &
     &                   LBi, UBi, LBj, UBj, 0, N(ng), scale,           &
#    ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#    endif
     &                   MIXING(ng) % Akt(:,:,:,isalt),                 &
     &                   SetFillVal = .FALSE.,                          &
     &                   MinValue = Fmin,                               &
     &                   MaxValue = Fmax)
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idSdif)),                      &
     &                      TRIM(INI(ng)%name), INI(ng)%Rindex
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idSdif)), Fmin, Fmax
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
      CALL netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' WRT_INI - error while writing variable: ',a,/,11x,    &
     &        'into initial NetCDF file for time record: ',i4)
  20  FORMAT ('NLM: WRT_INI - Writing initial state fields,',t75,a,     &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')')
  30  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE wrt_ini

# if defined ADJUST_BOUNDARY || \
     defined ADJUST_WSTRESS  || defined ADJUST_STFLUX
!
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
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
      USE strings_mod,         ONLY : FoundError
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
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc"
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
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          INI(ng)%Vid(idSbry(isFsur)),            &
     &                          OutRec, r2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          INI(ng)%Vid(idSbry(isUbar)),            &
     &                          OutRec, u2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          INI(ng)%Vid(idSbry(isVbar)),            &
     &                          OutRec, v2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          INI(ng)%Vid(idSbry(isUvel)),            &
     &                          OutRec, u3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          INI(ng)%Vid(idSbry(isVvel)),            &
     &                          OutRec, v3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp
          status=nf_fwrite3d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid, &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            INI(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                            OutRec, r3dvar,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ustr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % vstr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
            scale=1.0_dp
          END IF
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid,                    &
     &                       INI(ng)%Vid(idTsur(itrc)),                 &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      CALL netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master) THEN
#  ifdef NESTING
        WRITE (stdout,20) outer, inner, Tindex, OutRec, ng
#  else
        WRITE (stdout,20) outer, inner, Tindex, OutRec
#  endif
      END IF
!
  10  FORMAT (/,' WRT_FRC - error while writing variable: ',a,/,11x,    &
     &        'into initial NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_FRC   - wrote forcing  fields (Outer=',i2.2,      &
#  ifdef NESTING
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,', Grid ',      &
     &        i0,')')
#  else
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE wrt_frc
!
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
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
      USE strings_mod,         ONLY : FoundError
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
!
      real(dp) :: scale
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_frc_AD"
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
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out initial conditions.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          INI(ng)%Vid(idSbry(isFsur)),            &
     &                          OutRec, r2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,   &
     &                                                     Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          INI(ng)%Vid(idSbry(isUbar)),            &
     &                          OutRec, u2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,   &
     &                                                     Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite2d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          INI(ng)%Vid(idSbry(isVbar)),            &
     &                          OutRec, v2dvar,                         &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,   &
     &                                                     Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          INI(ng)%Vid(idSbry(isUvel)),            &
     &                          OutRec, u3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,    &
     &                                                  Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
        scale=1.0_dp
        status=nf_fwrite3d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid,   &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          INI(ng)%Vid(idSbry(isVvel)),            &
     &                          OutRec, v3dvar,                         &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,    &
     &                                                  Tindex))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp
          status=nf_fwrite3d_bry (ng, iNLM, INI(ng)%name, INI(ng)%ncid, &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            INI(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                            OutRec, r3dvar,                       &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,  &
     &                                                    Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*u3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idUsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % umask,                              &
#   endif
     &                   FORCES(ng) % ad_ustr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      scale=1.0_dp
      gtype=gfactor*v3dvar
      status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid, INI(ng)%Vid(idVsms),   &
     &                   OutRec, gtype,                                 &
     &                   LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,       &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                              &
#   endif
     &                   FORCES(ng) % ad_vstr(:,:,:,Tindex))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
          scale=1.0_dp
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, INI(ng)%ncid,                    &
     &                       INI(ng)%Vid(idTsur(itrc)),                 &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,   &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       FORCES(ng) % ad_tflux(:,:,:,Tindex,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
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
      CALL netcdf_sync (ng, iNLM, INI(ng)%name, INI(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master) THEN
# ifdef NESTING
        WRITE (stdout,20) outer, inner, Tindex, OutRec, ng
# else
        WRITE (stdout,20) outer, inner, Tindex, OutRec
# endif
      END IF
!
  10  FORMAT (/,' WRT_FRC_AD - error while writing variable: ',a,/,11x, &
     &        'into initial NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_FRC_AD   - wrote forcing  fields (Outer=',i2.2,   &
#  ifdef NESTING
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,', Grid ',      &
     &        i0,')')
#  else
     &        ', Inner=',i3.3,', Index=',i0,', Rec=',i0,')')
#  endif
!
      RETURN
      END SUBROUTINE wrt_frc_AD
# endif
#endif
      END MODULE wrt_ini_mod
