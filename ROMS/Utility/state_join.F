#include "cppdefs.h"
      SUBROUTINE state_join (ng, model, InpName, OutName, InpStrRec,    &
     &                       InpEndRec, OutStrRec)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and writes ROMS NLM background state solution    !
!  for specified input and output records. Primarily, it is used to    !
!  concatenate the NLM solution from concurrent interval window        !
!  history NetCDF files. The resulting solution is used to linearize   !
!  the tangent linear and adjoint kernels.                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     InpName      Input  history filename (string)                    !
!     InpName      Output history filename (string)                    !
!     InpStrRec    Starting Input  history time record (integer)       !
!     InpEndRec    Ending   Input  history time record (integer)       !
!     OutStrRec    Starting Output history time record (integer)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupling
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
!
      USE dateclock_mod,   ONLY : time_string
      USE nf_fread2d_mod,  ONLY : nf_fread2d
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
#ifdef SOLVE3D
      USE nf_fread3d_mod,  ONLY : nf_fread3d
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
#endif
      USE strings_mod,     ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: InpStrRec, InpEndRec

      integer, intent(inout) :: OutStrRec
!
      character (len=*) :: InpName, OutName
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: InpId, InpRec, OutId, OutRec, Tindex
      integer :: gtype, i, ic, status

      integer :: Vsize(4)
!
      real(r8) :: Fmin, Fmax
      real(dp) :: Fscl, stime
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
!  Process forward background solution.
!-----------------------------------------------------------------------
!
!  Open Input NetCDF file for reading.
!
      CALL netcdf_open (ng, model, InpName, 0, InpId)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,10) TRIM(InpName)
        RETURN
      END IF
!
!  Open Output NetCDF file for writing.
!
      CALL netcdf_open (ng, model, OutName, 1, OutId)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,10) TRIM(InpName)
        RETURN
      END IF
!
!  Inquire about the input variables.
!
      CALL netcdf_inq_var (ng, model, InpName,                          &
     &                     ncid = InpId)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  Scan the variable list and read in needed variables.
!
      OutRec=OutStrRec-1
      Tindex=1
      Fscl=1.0_dp
!
      REC_LOOP : DO InpRec=InpStrRec,InpEndRec
        OutRec=OutRec+1
        VAR_LOOP : DO i=1,n_var
!
!  Time.
!
          CHECK1 : IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            CALL netcdf_get_time (ng, model, InpName,                   &
     &                            TRIM(var_name(i)),                    &
     &                            Rclock%DateNumber, stime,             &
     &                            ncid = InpID,                         &
     &                            start = (/InpRec/),                   &
     &                            total = (/1/))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idtime)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              RETURN
            END IF
!
            CALL netcdf_put_fvar (ng, model, OutName,                   &
     &                            TRIM(var_name(i)), stime,             &
     &                            (/OutRec/), (/1/),                    &
     &                            ncid = OutId,                         &
     &                            varid = var_id(i))
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idtime)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              RETURN
            ELSE
              CALL time_string (stime, t_code)
              WRITE (Tstring,'(f15.4)') stime*sec2day
              IF (Master) THEN
                WRITE (stdout,20) t_code,                               &
     &                            ng, TRIM(ADJUSTL(Tstring)), InpRec,   &
     &                            Tindex, TRIM(InpName),                &
     &                            ng, TRIM(ADJUSTL(Tstring)), OutRec,   &
     &                            Tindex, TRIM(OutName)
              END IF
            END IF
!
!  Free-surface.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            gtype=r2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#endif
     &                        OCEAN(ng) % zeta(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idFsur)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
#ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#endif
     &                         OCEAN(ng) % zeta(:,:,Tindex),            &
     &                         SetFillVal = .FALSE.)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idFsur)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idFsur)), Fmin, Fmax
              END IF
            END IF

#ifdef FORWARD_RHS
!
!  Free-surface equation right-hand-side term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
            gtype=r2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % rmask,                         &
# endif
     &                        OCEAN(ng) % rzeta(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRzet)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % rmask,                        &
# endif
     &                         OCEAN(ng) % rzeta(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRzet)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRzet)), Fmin, Fmax
              END IF
            END IF
#endif
!
!  2D U-momentum component.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            gtype=u2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % umask_full,                    &
#endif
     &                        OCEAN(ng) % ubar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idUbar)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
#ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#endif
     &                         OCEAN(ng) % ubar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idUbar)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idUbar)), Fmin, Fmax
              END IF
            END IF

#ifdef FORWARD_RHS
!
!  2D U-equation right-hand-side term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            gtype=u2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask_full,                    &
# endif
     &                        OCEAN(ng) % rubar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRu2d)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
# endif
     &                         OCEAN(ng) % rubar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRu2d)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRu2d)), Fmin, Fmax
              END IF
            END IF
#endif
!
!  2D V-momentum component.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            gtype=v2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
#endif
     &                        OCEAN(ng) % vbar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idVbar)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
#ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#endif
     &                         OCEAN(ng) % vbar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idVbar)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idVbar)), Fmin, Fmax
              END IF
            END IF

#ifdef FORWARD_RHS
!
!  2D V-equation right-hand-side term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            gtype=v2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
# endif
     &                        OCEAN(ng) % rvbar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRv2d)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         OCEAN(ng) % rvbar(:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRv2d)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRv2d)), Fmin, Fmax
              END IF
            END IF
#endif
#ifdef SOLVE3D
# ifdef FORWARD_RHS
!
!  2D U-equation right-hand-side forcing term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRuct))) THEN
            gtype=u2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask_full,                    &
#  endif
     &                        OCEAN(ng) % rufrc)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRuct)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
#  ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#  endif
     &                         OCEAN(ng) % rufrc)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRuct)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRuct)), Fmin, Fmax
              END IF
            END IF
# endif
!
!  Time-averaged U-flux component for 2D equations.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            gtype=u2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
     &                        COUPLING(ng) % DU_avg1)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idUfx1)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % umask,                        &
# endif
     &                         COUPLING(ng) % DU_avg1)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idUfx1)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idUfx1)), Fmin, Fmax
              END IF
            END IF
!
!  Time-averaged U-flux component for 3D equations.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            gtype=u2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask,                         &
# endif
     &                        COUPLING(ng) % DU_avg2)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idUfx2)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % umask,                        &
# endif
     &                         COUPLING(ng) % DU_avg2)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idUfx2)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idUfx2)), Fmin, Fmax
              END IF
            END IF

# ifdef FORWARD_RHS
!
!  2D V-equation right-hand-side forcing term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRvct))) THEN
            gtype=v2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
#  endif
     &                        OCEAN(ng) % rvfrc)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRvct)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
#  ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#  endif
     &                         OCEAN(ng) % rvfrc)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRvct)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRvct)), Fmin, Fmax
              END IF
            END IF
# endif
!
!  Time-averaged V-flux component for 2D equations.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            gtype=v2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
# endif
     &                        COUPLING(ng) % DV_avg1)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idVfx1)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         COUPLING(ng) % DV_avg1)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idVfx1)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idVfx1)), Fmin, Fmax
              END IF
            END IF
!
!  Time-averaged U-flux component for 3D equations.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            gtype=v2dvar
            status=nf_fread2d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
# endif
     &                        COUPLING(ng) % DV_avg2)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idVfx2)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite2d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, Fscl,                &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         COUPLING(ng) % DV_avg2)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idVfx2)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idVfx2)), Fmin, Fmax
              END IF
            END IF
!
!  3D U-momentum component.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            gtype=u3dvar
            status=nf_fread3d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % umask_full,                    &
# endif
     &                        OCEAN(ng) % u(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idUvel)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite3d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), Fscl,      &
# ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
# endif
     &                         OCEAN(ng) % u(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idUvel)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idUvel)), Fmin, Fmax
              END IF
            END IF

# ifdef FORWARD_RHS
!
!  3D U-momentum right-hand-side term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            gtype=u3dvar
            status=nf_fread3d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % umask_full,                    &
#  endif
     &                        OCEAN(ng) % ru(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRu3d)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite3d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), Fscl,      &
#  ifdef MASKING
     &                         GRID(ng) % umask_full,                   &
#  endif
     &                         OCEAN(ng) % ru(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRu3d)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRu3d)), Fmin, Fmax
              END IF
            END IF
# endif
!
!  3D V-momentum component.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            gtype=v3dvar
            status=nf_fread3d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
# ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
# endif
     &                        OCEAN(ng) % v(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idVvel)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite3d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), Fscl,      &
# ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
# endif
     &                         OCEAN(ng) % v(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idVvel)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idVvel)), Fmin, Fmax
              END IF
            END IF

# ifdef FORWARD_RHS
!
!  3D U-momentum right-hand-side term.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            gtype=v3dvar
            status=nf_fread3d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % vmask_full,                    &
#  endif
     &                        OCEAN(ng) % rv(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idRv3d)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite3d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), Fscl,      &
#  ifdef MASKING
     &                         GRID(ng) % vmask_full,                   &
#  endif
     &                         OCEAN(ng) % rv(:,:,:,Tindex))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idRv3d)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idRv3d)), Fmin, Fmax
              END IF
            END IF
# endif
# ifdef FORWARD_MIXING
!
!  Vertical viscosity.
!
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            gtype=w3dvar
            status=nf_fread3d(ng, model, InpName, InpId,                &
     &                        var_name(i), var_id(i),                   &
     &                        InpRec, gtype, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        Fscl, Fmin, Fmax,                         &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        MIXING(ng) % Akv)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idVvis)), InpRec,        &
     &                            TRIM(InpName)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
!
            status=nf_fwrite3d(ng, model, OutId,                        &
     &                         var_id(i), OutRec, gtype,                &
     &                         LBi, UBi, LBj, UBj, 0, N(ng), Fscl,      &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         MIXING(ng) % Akv)
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idVvis)), OutRec,        &
     &                            TRIM(OutName)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(2,idVvis)), Fmin, Fmax
              END IF
            END IF
# endif
#endif
          END IF CHECK1

#ifdef SOLVE3D
!
!  Tracer type variables.
!
          TRACER1_LOOP : DO itrc=1,NT(ng)
            gtype=r3dvar
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              status=nf_fread3d(ng, model, InpName, InpId,              &
     &                          var_name(i), var_id(i),                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
     &                          OCEAN(ng) % t(:,:,:,Tindex,itrc))
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))),        &
     &                              InpRec, TRIM(InpName)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              END IF
!
              status=nf_fwrite3d(ng, model, OutId,                      &
     &                           var_id(i), OutRec, gtype,              &
     &                           LBi, UBi, LBj, UBj, 1, N(ng), Fscl,    &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           OCEAN(ng) % t(:,:,:,Tindex,itrc))
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))),        &
     &                              OutRec, TRIM(OutName)
                END IF
                exit_flag=3
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
                  WRITE (stdout,50) TRIM(Vname(2,idTvar(itrc))),        &
     &                              Fmin, Fmax
                END IF
              END IF
            END IF
          END DO TRACER1_LOOP
!
!  Tracer type variables vertical diffusion.
!
          TRACER2_LOOP : DO itrc=1,NAT
            gtype=w3dvar
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDiff(itrc)))) THEN
              status=nf_fread3d(ng, model, InpName, InpId,              &
     &                          var_name(i), var_id(i),                 &
     &                          InpRec, gtype, Vsize,                   &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
     &                          MIXING(ng) % Akt(:,:,:,itrc))
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,30) TRIM(Vname(1,idDiff(itrc))),        &
     &                              InpRec, TRIM(InpName)
                END IF
                exit_flag=2
                ioerror=status
                RETURN
              END IF
!
              status=nf_fwrite3d(ng, model, OutId,                      &
     &                           var_id(i), OutRec, gtype,              &
     &                           LBi, UBi, LBj, UBj, 0, N(ng), Fscl,    &
# ifdef MASKING
     &                           GRID(ng) % rmask,                      &
# endif
     &                           MIXING(ng) % Akt(:,:,:,itrc))
              IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
                IF (Master) THEN
                  WRITE (stdout,40) TRIM(Vname(1,idDiff(itrc))),        &
     &                              OutRec, TRIM(OutName)
                END IF
                exit_flag=3
                ioerror=status
                RETURN
              ELSE
                IF (Master) THEN
                  WRITE (stdout,50) TRIM(Vname(2,idDiff(itrc))),        &
     &                              Fmin, Fmax
                END IF
              END IF
            END IF
          END DO TRACER2_LOOP
#endif

        END DO VAR_LOOP

      END DO REC_LOOP
!
!  Update output record value.
!
      OutStrRec=OutRec
!
!  Close input and output files for compleate synchronization.
!
      CALL netcdf_close (ng, iNLM, InpId, InpName, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL netcdf_close (ng, iNLM, OutId, OutName, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (/,' STATE_JOIN - unable to open grid NetCDF file: ',a)
  20  FORMAT ('NLM: STATE_JOIN - Concatenating state fields,',t75,a,    &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', InpRec=',i4.4,          &
     &         ', Index=',i1,', InpFile: ',a,')',                       &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', OutRec=',i4.4,          &
     &         ', Index=',i1,', OutFile: ',a,')')
  30  FORMAT (/,' STATE_JOIN - error while reading variable: ',a,2x,    &
     &        'at time record = ',i0,/,14x,'in input NetCDF file: ',a)
  40  FORMAT (/,' STATE_JOIN - error while writing variable: ',a,2x,    &
     &        'at time record = ',i0,/,14x,'in output NetCDF file: ',a)
  50  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE state_join
