#include "cppdefs.h"
#if defined WEAK_CONSTRAINT && defined TIME_CONV
#undef ENDPOINT_TRAPEZOIDAL

      SUBROUTINE time_corr (ng, tile, model, INPncname)
!
!svn $Id: time_corr.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine performs the time correlation. It operates on the      !
!  adjoint state.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng              Nested grid number.                              !
!     tile            Domain partition.                                !
!     model           Calling model identifier.                        !
!     INPncname       Input adjoint solution NetCDF file name.         !
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
      USE nf_fread2d_mod, ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod, ONLY : nf_fread3d
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Iinp, Iout, Irec, MyRec, MyType, Nrec
      integer :: INPncid, INPvid
      integer :: i, gtype, status, varid
      integer :: Vsize(4)
      integer :: j, k, itrc

      real(r8) :: Fmin, Fmax, scale
      real(r8) :: inp_time(1)
      real(r8) :: timeI, timeR, fac

# include "set_bounds.h"
!
      SourceFile='time_corr.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, model, INPncname)
      IF (exit_flag.ne.NoError) RETURN
      Nrec=rec_size
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, INPncname)
      IF (exit_flag.ne.NoError) RETURN
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
      IF (exit_flag.ne.NoError) THEN
        WRITE (stdout,10) TRIM(INPncname)
        RETURN
      END IF
!
!  Compute the time convolutions for the model error corrections.
!
      Iinp=1
      scale=1.0_r8
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
            IF (status.ne.nf90_noerr) THEN
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
          CALL netcdf_get_fvar (ng, model, INPncname, Vname(1,idtime),  &
     &                          inp_time,                               &
     &                          ncid = INPncid,                         &
     &                          start = (/Irec/), total = (/1/))
          IF (exit_flag.ne.NoError) RETURN

          CALL netcdf_put_fvar (ng, model, TLF(ng)%name,                &
     &                          Vname(1,idtime), inp_time,              &
     &                          (/Iout/), (/1/),                        &
     &                          ncid = TLF(ng)%ncid)
          IF (exit_flag.ne.NoError) RETURN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
          IF (status.ne.nf90_noerr) THEN
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
            IF (status.ne.nf90_noerr) THEN
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
      CALL netcdf_sync (ng, model, INPncname, TLF(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN

      IF (Master)  WRITE (stdout,50) Nrec-1, TRIM(TLF(ng)%name)
!
  10  FORMAT (/,' TIME_CORR - unable to open input NetCDF file: ',a)
  20  FORMAT (/,' TIME_CORR - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,17x,'in input NetCDF file: ',a)
  30  FORMAT (/,' TIME_CORR - error while writing variable: ',a,2x,    &
     &        'at time record = ',i3,/,17x,'into NetCDF file: ',a)
  40  FORMAT (/,' TIME_CORR - cannot find state variable: ',a,         &
     &        /,12x,'in input NetCDF file: ',a)
  50  FORMAT (4x,'TIME_CORR - wrote convolved adjoint impulses, ',     &
     &           'records: 001 to ',i3.3,/,18x,'file: ',a)

      RETURN
      END SUBROUTINE time_corr
#else
      SUBROUTINE time_corr
      END SUBROUTINE time_corr
#endif
