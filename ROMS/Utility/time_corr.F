#include "cppdefs.h"
      MODULE time_corr_mod
#if defined WEAK_CONSTRAINT && defined TIME_CONV
#undef ENDPOINT_TRAPEZOIDAL
!
!git $Id$
!svn $Id: time_corr.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module performs the time correlation. It operates on the       !
!  adjoint state. It uses either the standard NetCDF library or        !
!  the Parallel-IO (PIO) library.                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng              Nested grid number (integer)                     !
!     tile            Domain tile partition (integer)                  !
!     model           Calling model identifier (integer)               !
!     IOtype          I/O NetCDF library to use (integer)              !
!     INPncname       Input adjoint solution NetCDF file name (string) !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
!
      USE nf_fread2d_mod,  ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod,  ONLY : nf_fread3d
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      PUBLIC  :: time_corr
      PRIVATE :: time_corr_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: time_corr_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE time_corr (ng, tile, model, IOtype, INPncname)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, IOtype
!
      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Open existing nonlinear initial conditions NetCDF file and inquire
!  about its contains and define new variables, if needed.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (IOtype)
        CASE (io_nf90)
          CALL time_corr_nf90 (ng, tile, model, INPncname,              &
     &                         LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL time_corr_pio (ng, tile, model, INPncname,               &
     &                        LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' TIME_CORR - Illegal input file type, io_type = ',i0,    &
     &        /,13x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE time_corr
!
!***********************************************************************
      SUBROUTINE time_corr_nf90 (ng, tile, model, INPncname,            &
     &                           LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      integer :: Iinp, Iout, Irec, MyRec, MyType, Nrec
      integer :: INPncid, INPvid
      integer :: i, gtype, status, varid
      integer :: Vsize(4)
      integer :: j, k, itrc
!
      real(r8) :: Fmin, Fmax
!
      real(dp) :: scale
      real(dp) :: inp_time(1)
      real(dp) :: timeI, timeR, fac
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", time_corr_nf90"

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, model, INPncname)
      IF (FoundError(exit_flag, NoError, __LINE__,  MyFile)) RETURN
      Nrec=rec_size
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, INPncname)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!-----------------------------------------------------------------------
!  Read adjoint solution for time convolutions.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.
!
      CALL netcdf_open (ng, model, INPncname, 0, INPncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,10) TRIM(INPncname)
        RETURN
      END IF
!
!  Compute the time convolutions for the model error corrections.
!
      Iinp=1
      scale=1.0_dp
!
!  Read all free-surface records into the "f_zetaS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idZtlf), INPvid)) THEN
        gtype=var_flag(INPvid)*r2dvar
        MyType=gtype
        DO MyRec=1,Nrec-1
          status=nf_fread2d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idZtlf), INPvid,                    &
     &                      MyRec, MyType, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % f_zetaS(:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idZtlf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZtlf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all 2D U-momentum records into the "f_ubarS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idUbtf), INPvid)) THEN
        gtype=var_flag(INPvid)*u2dvar
        MyType=gtype
        DO MyRec=1,Nrec-1
          status=nf_fread2d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idUbtf), INPvid,                    &
     &                      MyRec, MyType, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
# endif
     &                      OCEAN(ng) % f_ubarS(:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idUbtf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUbtf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all 2D V-momentum records into the "f_vbarS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idVbtf), INPvid)) THEN
        gtype=var_flag(INPvid)*v2dvar
        MyType=gtype
        DO MyRec=1,Nrec-1
          status=nf_fread2d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idVbtf), INPvid,                    &
     &                      MyRec, MyType, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
# endif
     &                      OCEAN(ng) % f_vbarS(:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idVbtf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbtf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF

# ifdef SOLVE3D
!
!  Read all 3D U-momentum records into the "f_uS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idUtlf), INPvid)) THEN
        gtype=var_flag(INPvid)*u3dvar
        MyType=gtype
        DO MyRec=1,Nrec-1
          status=nf_fread3d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idUtlf), INPvid,                    &
     &                      MyRec, MyType, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
#  endif
     &                      OCEAN(ng) % f_uS(:,:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idUtlf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUtlf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all 3D V-momentum records into the "f_vS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idVtlf), INPvid)) THEN
        gtype=var_flag(INPvid)*v3dvar
        MyType=gtype
        DO MyRec=1,Nrec-1
          status=nf_fread3d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idVtlf), INPvid,                    &
     &                      MyRec, MyType, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
#  endif
     &                      OCEAN(ng) % f_vS(:,:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idVtlf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVtlf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all the tracers records into the "f_tS" array.
!
      DO itrc=1,NT(ng)
        IF (find_string(var_name, n_var, Vname(1,idTtlf(itrc)),         &
     &                  INPvid)) THEN
          gtype=var_flag(INPvid)*r3dvar
          MyType=gtype
          DO MyRec=1,Nrec-1
            status=nf_fread3d(ng, model, INPncname, INPncid,            &
     &                        Vname(1,idTtlf(itrc)), INPvid,            &
     &                        MyRec, MyType, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % f_tS(:,:,:,MyRec,itrc))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Vname(1,idTtlf(itrc))), MyRec,   &
     &                            TRIM(INPncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END DO
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTtlf(i))),       &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# endif
!
!=======================================================================
!  Perform time convolution and write impulse forcing into output TLF
!  NetCDF file.
!=======================================================================
!
      Iout=0
      REC_LOOP : DO Irec=Nrec-1,1,-1
        Iout=Iout+1
        timeI=dstart*day2sec+(Iout-1)*nADJ(ng)*dt(ng)
!
!  Write out time variable.
!
        IF (find_string(var_name, n_var, Vname(1,idtime), INPvid)) THEN
          CALL netcdf_get_time (ng, model, INPncname, Vname(1,idtime),  &
     &                          Rclock%DateNumber, inp_time,            &
     &                          ncid = INPncid,                         &
     &                          start = (/Irec/), total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          CALL netcdf_put_fvar (ng, model, TLF(ng)%name,                &
     &                          Vname(1,idtime), inp_time,              &
     &                          (/Iout/), (/1/),                        &
     &                          ncid = TLF(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
        END IF
!
!-----------------------------------------------------------------------
!  Free-surface weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            OCEAN(ng) % ad_zeta(i,j,Iinp)=0.0_r8
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idZtlf), INPvid)) THEN
          gtype=var_flag(INPvid)*r2dvar
          MyType=gtype
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isFsur,ng))
# ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
# endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                OCEAN(ng) % ad_zeta(i,j,Iinp)=                          &
     &                              OCEAN(ng) % ad_zeta(i,j,Iinp)+      &
     &                              fac*OCEAN(ng) % f_zetaS(i,j,MyRec)
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          MyType=gtype
          status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idZtlf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % ad_zeta(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idZtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D U-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            OCEAN(ng) % ad_ubar(i,j,Iinp)=0.0_r8
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idUbtf), INPvid)) THEN
          gtype=var_flag(INPvid)*u2dvar
          MyType=gtype
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isUbar,ng))
# ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
# endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO j=JstrR,JendR
             DO i=Istr,IendR
               OCEAN(ng) % ad_ubar(i,j,Iinp)=                           &
     &                             OCEAN(ng) % ad_ubar(i,j,Iinp)+       &
     &                             fac*OCEAN(ng) % f_ubarS(i,j,MyRec)
             END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          MyType=gtype
          status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idUbtf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
# endif
     &                       OCEAN(ng) % ad_ubar(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idUbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D V-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            OCEAN(ng) % ad_vbar(i,j,Iinp)=0.0_r8
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idVbtf), INPvid)) THEN
          gtype=var_flag(INPvid)*v2dvar
          MyType=gtype
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isVbar,ng))
# ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
# endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                OCEAN(ng) % ad_vbar(i,j,Iinp)=                          &
     &                              OCEAN(ng) % ad_vbar(i,j,Iinp)+      &
     &                              fac*OCEAN(ng) % f_vbarS(i,j,MyRec)
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idVbtf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
# endif
     &                       OCEAN(ng) % ad_vbar(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idVbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D U-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              OCEAN(ng) % ad_u(i,j,k,Iinp)=0.0_r8
            END DO
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idUtlf), INPvid)) THEN
          gtype=var_flag(INPvid)*u3dvar
          MyType=gtype
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isUvel,ng))
#  ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
#  endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  OCEAN(ng) % ad_u(i,j,k,Iinp)=                         &
     &                             OCEAN(ng) % ad_u(i,j,k,Iinp)+        &
     &                             fac*OCEAN(ng) % f_uS(i,j,k,MyRec)
                END DO
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite3d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idUtlf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#  endif
     &                       OCEAN(ng) % ad_u(:,:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idUtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D V-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              OCEAN(ng) % ad_v(i,j,k,Iinp)=0.0_r8
            END DO
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idVtlf), INPvid)) THEN
          gtype=var_flag(INPvid)*v3dvar
          MyType=gtype
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isVvel,ng))
#  ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
#  endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  OCEAN(ng) % ad_v(i,j,k,Iinp)=                         &
     &                             OCEAN(ng) % ad_v(i,j,k,Iinp)+        &
     &                             fac*OCEAN(ng) % f_vS(i,j,k,MyRec)
                END DO
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite3d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idVtlf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#  endif
     &                       OCEAN(ng) % ad_v(:,:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idVtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Tracers type variables weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                OCEAN(ng) % ad_t(i,j,k,Iinp,itrc)=0.0_r8
              END DO
            END DO
          END DO
        END DO
!
!  Convolve.
!
        DO itrc=1,NT(ng)
          IF (find_string(var_name, n_var, Vname(1,idTtlf(itrc)),       &
     &                    INPvid)) THEN
            gtype=var_flag(INPvid)*r3dvar
            MyType=gtype
            DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
              timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
              fac=EXP(-ABS(timeI-timeR)/Tdecay(isTvar(itrc),ng))
#  ifdef ENDPOINT_TRAPEZOIDAL
              IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
                fac=0.5_r8*fac
              END IF
#  endif
!
!  Form the weighted sum of the adjoint file records in time.
!
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    OCEAN(ng) % ad_t(i,j,k,Iinp,itrc)=                  &
     &                          OCEAN(ng) % ad_t(i,j,k,Iinp,itrc)+      &
     &                          fac*OCEAN(ng) % f_tS(i,j,k,MyRec,itrc)
                  END DO
                END DO
              END DO
            END DO
!
!  Write out convolved solution.
!
            status=nf_fwrite3d(ng, model, TLF(ng)%ncid,                 &
     &                         TLF(ng)%Tid(itrc),                       &
     &                         Iout, MyType,                            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         OCEAN(ng) % ad_t(:,:,:,Iinp,itrc))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idTtlf(itrc))), Irec,    &
     &                            TRIM(TLF(ng)%name)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTtlf(i))),     &
     &                                    TRIM(INPncname)
            exit_flag=2
            RETURN
          END IF
        END DO
# endif
      END DO REC_LOOP
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, TLF(ng)%name, TLF(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      IF (Master)  WRITE (stdout,50) Nrec-1, TRIM(TLF(ng)%name)
!
  10  FORMAT (/,' TIME_CORR_NF90 - unable to open input NetCDF',        &
     &        ' file: ',a)
  20  FORMAT (/,' TIME_CORR_NF90 - error while reading variable: ',a,   &
     &        2x,'at time record = ',i3,/,18x,'in input NetCDF file:',  &
     &        1x,a)
  30  FORMAT (/,' TIME_CORR_NF90 - error while writing variable: ',a,   &
     &        2x,'at time record = ',i3,/,18x,'into NetCDF file: ',a)
  40  FORMAT (/,' TIME_CORR_NF90 - cannot find state variable: ',a,     &
     &        /,18x,'in input NetCDF file: ',a)
  50  FORMAT (2x,'TIME_CORR_NF90   - wrote convolved adjoint impulses', &
     &        ', records: 001 to ',i3.3,/,21x,'file: ',a)
!
      RETURN
      END SUBROUTINE time_corr_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE time_corr_pio (ng, tile, model, INPncname,             &
     &                          LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      integer :: Iinp, Iout, Irec, MyRec, Nrec
      integer :: i, status, vindex
      integer :: Vsize(4)
      integer :: j, k, itrc
!
      real(r8) :: Fmin, Fmax
!
      real(dp) :: scale
      real(dp) :: inp_time(1)
      real(dp) :: timeI, timeR, fac
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", time_corr_pio"
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
      TYPE (My_VarDesc), pointer :: pioVar

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Inquire about the dimensions and check for consistency.
!
      CALL pio_netcdf_check_dim (ng, model, INPncname)
      IF (FoundError(exit_flag, NoError, __LINE__,  MyFile)) RETURN
      Nrec=rec_size
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, model, INPncname)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!-----------------------------------------------------------------------
!  Read adjoint solution for time convolutions.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.
!
      CALL pio_netcdf_open (ng, model, INPncname, 0, INPpioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,10) TRIM(INPncname)
        RETURN
      END IF
!
!  Compute the time convolutions for the model error corrections.
!
      Iinp=1
      scale=1.0_dp
!
!  Read all free-surface records into the "f_zetaS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idZtlf), vindex)) THEN

        pioVar%vd => var_desc(vindex)
        IF (KIND(OCEAN(ng)%f_zetaS).eq.8) THEN
          pioVar%dkind=PIO_double
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          pioVar%dkind=PIO_real
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF
        pioVar%gtype=r2dvar
!
        DO MyRec=1,Nrec-1
          status=nf_fread2d(ng, model, INPncname, INPpioFile,           &
     &                      Vname(1,idZtlf), pioVar,                    &
     &                      MyRec, ioDesc, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % f_zetaS(:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idZtlf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZtlf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all 2D U-momentum records into the "f_ubarS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idUbtf), vindex)) THEN

        pioVar%vd => var_desc(vindex)
        IF (KIND(OCEAN(ng)%f_ubarS).eq.8) THEN
          pioVar%dkind=PIO_double
          ioDesc => ioDesc_dp_u2dvar(ng)
        ELSE
          pioVar%dkind=PIO_real
          ioDesc => ioDesc_sp_u2dvar(ng)
        END IF
        pioVar%gtype=u2dvar
!
        DO MyRec=1,Nrec-1
          status=nf_fread2d(ng, model, INPncname, INPpioFile,           &
     &                      Vname(1,idUbtf), pioVar,                    &
     &                      MyRec, ioDesc, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
# endif
     &                      OCEAN(ng) % f_ubarS(:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idUbtf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUbtf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all 2D V-momentum records into the "f_vbarS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idVbtf), vindex)) THEN

        pioVar%vd => var_desc(vindex)
        IF (KIND(OCEAN(ng)%f_vbarS).eq.8) THEN
          pioVar%dkind=PIO_double
          ioDesc => ioDesc_dp_v2dvar(ng)
        ELSE
          pioVar%dkind=PIO_real
          ioDesc => ioDesc_sp_v2dvar(ng)
        END IF
        pioVar%gtype=v2dvar
!
        DO MyRec=1,Nrec-1
          status=nf_fread2d(ng, model, INPncname, INPpioFile,           &
     &                      Vname(1,idVbtf), pioVar,                    &
     &                      MyRec, ioDesc, Vsize,                       &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
# endif
     &                      OCEAN(ng) % f_vbarS(:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idVbtf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbtf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF

# ifdef SOLVE3D
!
!  Read all 3D U-momentum records into the "f_uS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idUtlf), vindex)) THEN

        pioVar%vd => var_desc(vindex)
        IF (KIND(OCEAN(ng)%f_uS).eq.8) THEN
          pioVar%dkind=PIO_double
          ioDesc => ioDesc_dp_u3dvar(ng)
        ELSE
          pioVar%dkind=PIO_real
          ioDesc => ioDesc_sp_u3dvar(ng)
        END IF
        pioVar%gtype=u3dvar
!
        DO MyRec=1,Nrec-1
          status=nf_fread3d(ng, model, INPncname, INPpioFile,           &
     &                      Vname(1,idUtlf), pioVar,                    &
     &                      MyRec, ioDesc, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
#  endif
     &                      OCEAN(ng) % f_uS(:,:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idUtlf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUtlf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all 3D V-momentum records into the "f_vS" array.
!
      IF (find_string(var_name, n_var, Vname(1,idVtlf), vindex)) THEN

        pioVar%vd => var_desc(vindex)
        IF (KIND(OCEAN(ng)%f_vS).eq.8) THEN
          pioVar%dkind=PIO_double
          ioDesc => ioDesc_dp_v3dvar(ng)
        ELSE
          pioVar%dkind=PIO_real
          ioDesc => ioDesc_sp_v3dvar(ng)
        END IF
        pioVar%gtype=v3dvar
!
        DO MyRec=1,Nrec-1
          status=nf_fread3d(ng, model, INPncname, INPpioFile,           &
     &                      Vname(1,idVtlf), pioVar,                    &
     &                      MyRec, ioDesc, Vsize,                       &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
#  endif
     &                      OCEAN(ng) % f_vS(:,:,:,MyRec))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idVtlf)), MyRec,           &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          END IF
        END DO
      ELSE
        IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVtlf)),            &
     &                                TRIM(INPncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read all the tracers records into the "f_tS" array.
!
      DO itrc=1,NT(ng)
        IF (find_string(var_name, n_var, Vname(1,idTtlf(itrc)),         &
     &                  vindex)) THEN

          pioVar%vd => var_desc(vindex)
          IF (KIND(OCEAN(ng)%f_tS).eq.8) THEN
            pioVar%dkind=PIO_double
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            pioVar%dkind=PIO_real
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF
          pioVar%gtype=r3dvar
!
          DO MyRec=1,Nrec-1
            status=nf_fread3d(ng, model, INPncname, INPpioFile,         &
     &                        Vname(1,idTtlf(itrc)), pioVar,            &
     &                        MyRec, ioDesc, Vsize,                     &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
     &                        OCEAN(ng) % f_tS(:,:,:,MyRec,itrc))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,20) TRIM(Vname(1,idTtlf(itrc))), MyRec,   &
     &                            TRIM(INPncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            END IF
          END DO
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTtlf(i))),       &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# endif
!
!=======================================================================
!  Perform time convolution and write impulse forcing into output TLF
!  NetCDF file.
!=======================================================================
!
      Iout=0
      REC_LOOP : DO Irec=Nrec-1,1,-1
        Iout=Iout+1
        timeI=dstart*day2sec+(Iout-1)*nADJ(ng)*dt(ng)
!
!  Write out time variable.
!
        IF (find_string(var_name, n_var, Vname(1,idtime), vindex)) THEN
          CALL pio_netcdf_get_time (ng, model, INPncname,               &
     &                              Vname(1,idtime),                    &
     &                              Rclock%DateNumber, inp_time,        &
     &                              pioFile = INPpioFile,               &
     &                              start = (/Irec/),                   &
     &                              total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL pio_netcdf_put_fvar (ng, model, TLF(ng)%name,            &
     &                              Vname(1,idtime), inp_time,          &
     &                              (/Iout/), (/1/),                    &
     &                              pioFile = TLF(ng)%pioFile)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
        END IF
!
!-----------------------------------------------------------------------
!  Free-surface weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            OCEAN(ng) % ad_zeta(i,j,Iinp)=0.0_r8
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idZtlf), vindex)) THEN

          IF (TLF(ng)%pioVar(idZtlf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF

          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isFsur,ng))
# ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
# endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                OCEAN(ng) % ad_zeta(i,j,Iinp)=                          &
     &                              OCEAN(ng) % ad_zeta(i,j,Iinp)+      &
     &                              fac*OCEAN(ng) % f_zetaS(i,j,MyRec)
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idZtlf),                    &
     &                       Iout, ioDesc,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % ad_zeta(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idZtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D U-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            OCEAN(ng) % ad_ubar(i,j,Iinp)=0.0_r8
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idUbtf), vindex)) THEN

          IF (TLF(ng)%pioVar(idUbtf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u2dvar(ng)
          END IF

          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isUbar,ng))
# ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
# endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO j=JstrR,JendR
             DO i=Istr,IendR
               OCEAN(ng) % ad_ubar(i,j,Iinp)=                           &
     &                             OCEAN(ng) % ad_ubar(i,j,Iinp)+       &
     &                             fac*OCEAN(ng) % f_ubarS(i,j,MyRec)
             END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idUbtf),                    &
     &                       Iout, ioDesc,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
# endif
     &                       OCEAN(ng) % ad_ubar(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idUbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  2D V-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            OCEAN(ng) % ad_vbar(i,j,Iinp)=0.0_r8
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idVbtf), vindex)) THEN

          IF (TLF(ng)%pioVar(idVbtf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v2dvar(ng)
          END IF
!
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isVbar,ng))
# ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
# endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                OCEAN(ng) % ad_vbar(i,j,Iinp)=                          &
     &                              OCEAN(ng) % ad_vbar(i,j,Iinp)+      &
     &                              fac*OCEAN(ng) % f_vbarS(i,j,MyRec)
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idVbtf),                    &
     &                       Iout, ioDesc,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
# endif
     &                       OCEAN(ng) % ad_vbar(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idVbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  3D U-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO k=1,N(ng)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              OCEAN(ng) % ad_u(i,j,k,Iinp)=0.0_r8
            END DO
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idUtlf), vindex)) THEN

          IF (TLF(ng)%pioVar(idUtlf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
!
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isUvel,ng))
#  ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
#  endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  OCEAN(ng) % ad_u(i,j,k,Iinp)=                         &
     &                             OCEAN(ng) % ad_u(i,j,k,Iinp)+        &
     &                             fac*OCEAN(ng) % f_uS(i,j,k,MyRec)
                END DO
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idUtlf),                    &
     &                       Iout, ioDesc,                              &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#  endif
     &                       OCEAN(ng) % ad_u(:,:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idUtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  3D V-momentum weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO k=1,N(ng)
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              OCEAN(ng) % ad_v(i,j,k,Iinp)=0.0_r8
            END DO
          END DO
        END DO
!
!  Convolve.
!
        IF (find_string(var_name, n_var, Vname(1,idVtlf), vindex)) THEN

          IF (TLF(ng)%pioVar(idVtlf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
!
          DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
            timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
            fac=EXP(-ABS(timeI-timeR)/Tdecay(isVvel,ng))
#  ifdef ENDPOINT_TRAPEZOIDAL
            IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
              fac=0.5_r8*fac
            END IF
#  endif
!
!  Form the weighted sum of the adjoint file records in time.
!
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  OCEAN(ng) % ad_v(i,j,k,Iinp)=                         &
     &                             OCEAN(ng) % ad_v(i,j,k,Iinp)+        &
     &                             fac*OCEAN(ng) % f_vS(i,j,k,MyRec)
                END DO
              END DO
            END DO
          END DO
!
!  Write out convolved solution.
!
          status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idVtlf),                    &
     &                       Iout, ioDesc,                              &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#  endif
     &                       OCEAN(ng) % ad_v(:,:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idVtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Tracers type variables weak-constraint time convolution.
!-----------------------------------------------------------------------
!
!  Clear adjoint state array.
!
        DO itrc=1,NT(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                OCEAN(ng) % ad_t(i,j,k,Iinp,itrc)=0.0_r8
              END DO
            END DO
          END DO
        END DO
!
!  Convolve.
!
        DO itrc=1,NT(ng)
          IF (find_string(var_name, n_var, Vname(1,idTtlf(itrc)),       &
     &                    vindex)) THEN

            IF (TLF(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r3dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r3dvar(ng)
            END IF
!
            DO MyRec=Nrec-1,1,-1
!
!  Compute weighting factors.
!
              timeR=dstart*day2sec+(Nrec-1-MyRec)*nADJ(ng)*dt(ng)
              fac=EXP(-ABS(timeI-timeR)/Tdecay(isTvar(itrc),ng))
#  ifdef ENDPOINT_TRAPEZOIDAL
              IF ((MyRec.eq.1).or.(MyRec.eq.(Nrec-1))) THEN
                fac=0.5_r8*fac
              END IF
#  endif
!
!  Form the weighted sum of the adjoint file records in time.
!
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    OCEAN(ng) % ad_t(i,j,k,Iinp,itrc)=                  &
     &                          OCEAN(ng) % ad_t(i,j,k,Iinp,itrc)+      &
     &                          fac*OCEAN(ng) % f_tS(i,j,k,MyRec,itrc)
                  END DO
                END DO
              END DO
            END DO
!
!  Write out convolved solution.
!
            status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,              &
     &                         TLF(ng)%pioTrc(itrc),                    &
     &                         Iout, ioDesc,                            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         OCEAN(ng) % ad_t(:,:,:,Iinp,itrc))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(1,idTtlf(itrc))), Irec,    &
     &                            TRIM(TLF(ng)%name)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTtlf(i))),     &
     &                                    TRIM(INPncname)
            exit_flag=2
            RETURN
          END IF
        END DO
# endif
      END DO REC_LOOP
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, TLF(ng)%name, TLF(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      IF (Master)  WRITE (stdout,50) Nrec-1, TRIM(TLF(ng)%name)
!
  10  FORMAT (/,' TIME_CORR_PIO - unable to open input NetCDF',         &
     &        ' file: ',a)
  20  FORMAT (/,' TIME_CORR_PIO - error while reading variable: ',a,    &
     &        2x,'at time record = ',i3,/,17x,'in input NetCDF file:',  &
     &        1x,a)
  30  FORMAT (/,' TIME_CORR_PIO - error while writing variable: ',a,    &
     &        2x,'at time record = ',i3,/,17x,'into NetCDF file: ',a)
  40  FORMAT (/,' TIME_CORR_PIO - cannot find state variable: ',a,      &
     &        /,18x,'in input NetCDF file: ',a)
  50  FORMAT (2x,'TIME_CORR_PIO   - wrote convolved adjoint impulses',  &
     &        ', records: 001 to ',i3.3,/,21x,'file: ',a)
!
      RETURN
      END SUBROUTINE time_corr_pio
# endif
#endif
      END MODULE time_corr_mod

