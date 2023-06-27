#include "cppdefs.h"
      MODULE wrt_impulse_mod
#if defined ADJOINT && defined IMPULSE
!
!git $Id$
!svn $Id: wrt_impulse.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These subroutine read in requested adjoint solution from input      !
!  NetCDF (saved at nADJ time-step intervals) and then writes  to      !
!  output  impulse  forcing  NetCDF  file in ascending time order      !
!  since it is processed by the TL and RP models.                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng              Nested grid number.                              !
!     tile            Domain partition.                                !
!     model           Calling model identifier.                        !
!     INPncname       Input adjoint solution NetCDF file name.         !
!                                                                      !
!  Notice that only momentum is affected by the full time-averaged     !
!  masks.  If applicable, these mask contains information about        !
!  river runoff and time-dependent wetting and drying variations.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
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
      PUBLIC  :: wrt_impulse
      PRIVATE :: wrt_impulse_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_impulse_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_impulse (ng, tile, model, INPncname)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
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
!  Write out time-averaged fields according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (LZE(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_impulse_nf90 (ng, tile, model, INPncname,            &
     &                           LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_impulse_pio (ng, tile, model, INPncname,             &
     &                          LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) LZE(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_IMPULSE - Illegal output file type, io_type = ',i0, &
     &        /,15x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_impulse
!
!***********************************************************************
      SUBROUTINE wrt_impulse_nf90 (ng, tile, model, INPncname,          &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Iinp, Iout, Irec, MyType, Nrec
      integer :: INPncid, INPvid
      integer :: i, gtype, status, varid
      integer :: ibuffer(2), Vsize(4)
# ifdef CHECKSUM
      integer(i8b) :: Fhash
# endif
!
      real(r8) :: Fmin, Fmax
      real(dp) :: scale
      real(dp) :: inp_time(1)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_impulse_nf90"

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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
# ifdef SP4DVAR
      IF (Master)  WRITE (stdout,10) Nrec, TRIM(TLF(ng)%name)
# else
      IF (Master)  WRITE (stdout,10) Nrec-1, TRIM(TLF(ng)%name)
# endif
!
!-----------------------------------------------------------------------
!  Read adjoint solution and convert to impulse forcing.  Then, write
!  impulse forcing into output NetCDF file.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.
!
      CALL netcdf_open (ng, model, INPncname, 0, INPncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,20) TRIM(INPncname)
        RETURN
      END IF
# ifdef SP4DVAR
!
      Iinp=1
      scale=1.0_dp
!
      DO Irec=1,Nrec
        Iout=Irec
# else
!
!  Process each record in input adjoint NetCDF except last. Note that
!  the adjoint records are processed backwards (Nrec-1:1) and written
!  in ascending time order (Iout initialized to 0) since the weak
!  constraint forcing will be read by the TL and RP models. Record
!  Nrec is not processed since it is not needed.
!
      Iinp=1
      Iout=0
      scale=1.0_dp
!
      DO Irec=Nrec-1,1,-1
        Iout=Iout+1
# endif
!
!  Process time.
!
        IF (find_string(var_name, n_var, Vname(1,idtime), INPvid)) THEN
          CALL netcdf_get_time (ng, model, INPncname, Vname(1,idtime),  &
     &                          Rclock%DateNumber, inp_time,            &
     &                          ncid = INPncid,                         &
     &                          start = (/Irec/), total = (/1/))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL netcdf_put_fvar (ng, model, TLF(ng)%name,                &
     &                          Vname(1,idtime), inp_time,              &
     &                          (/Iout/), (/1/),                        &
     &                          ncid = TLF(ng)%ncid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
        END IF
!
!  Process free-surface weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idZtlf), INPvid)) THEN
          gtype=var_flag(INPvid)*r2dvar
          MyType=gtype
          status=nf_fread2d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idZtlf), INPvid,                    &
     &                      Irec, MyType, Vsize,                        &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
# ifdef CHECKSUM
     &                      OCEAN(ng) % ad_zeta(:,:,Iinp),              &
     &                      checksum = Fhash)
# else
     &                      OCEAN(ng) % ad_zeta(:,:,Iinp))
# endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idZtlf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idZtlf)), Irec,            &
# ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
# else
     &                          Fmin, Fmax
# endif
            END IF
          END IF
!
          MyType=gtype
          status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idZtlf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask,                          &
# endif
     &                       OCEAN(ng) % ad_zeta(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idZtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF

# ifndef SOLVE3D
!
!  Process 2D U-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idUbtf), INPvid)) THEN
          gtype=var_flag(INPvid)*u2dvar
          MyType=gtype
          status=nf_fread2d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idUbtf), INPvid,                    &
     &                      Irec, MyType, Vsize,                        &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
#  endif
#  ifdef CHECKSUM
     &                      OCEAN(ng) % ad_ubar(:,:,Iinp),              &
     &                      checksum = Fhash)
#  else
     &                      OCEAN(ng) % ad_ubar(:,:,Iinp))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idUbtf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idUbtf)), Irec,            &
#  ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#  else
     &                          Fmin, Fmax
#  endif
            END IF
          END IF
!
          MyType=gtype
          status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idUbtf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#  endif
     &                       OCEAN(ng) % ad_ubar(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idUbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!  Process 2D V-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idVbtf), INPvid)) THEN
          gtype=var_flag(INPvid)*v2dvar
          MyType=gtype
          status=nf_fread2d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idVbtf), INPvid,                    &
     &                      Irec, MyType, Vsize,                        &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
#  endif
#  ifdef CHECKSUM
     &                      OCEAN(ng) % ad_vbar(:,:,Iinp),              &
     &                      checksum = Fhash)
#  else
     &                      OCEAN(ng) % ad_vbar(:,:,Iinp))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idVbtf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idVbtf)), Irec,            &
#  ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#  else
     &                          Fmin, Fmax
#  endif
            END IF
          END IF
!
          MyType=gtype
          status=nf_fwrite2d(ng, model, TLF(ng)%ncid,                   &
     &                       TLF(ng)%Vid(idVbtf),                       &
     &                       Iout, MyType,                              &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#  endif
     &                       OCEAN(ng) % ad_vbar(:,:,Iinp))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idVbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
# endif
# ifdef SOLVE3D
!
!  Process 3D U-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idUtlf), INPvid)) THEN
          gtype=var_flag(INPvid)*u3dvar
          MyType=gtype
          status=nf_fread3d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idUtlf), INPvid,                    &
     &                      Irec, MyType, Vsize,                        &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
#  endif
#  ifdef CHECKSUM
     &                      OCEAN(ng) % ad_u(:,:,:,Iinp),               &
     &                      checksum = Fhash)
#  else
     &                      OCEAN(ng) % ad_u(:,:,:,Iinp))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idUtlf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idUtlf)), Irec,            &
#  ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#  else
     &                          Fmin, Fmax
#  endif
            END IF
          END IF
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
              WRITE (stdout,60) TRIM(Vname(1,idUtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!  Process 3D V-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idVtlf), INPvid)) THEN
          gtype=var_flag(INPvid)*v3dvar
          MyType=gtype
          status=nf_fread3d(ng, model, INPncname, INPncid,              &
     &                      Vname(1,idVtlf), INPvid,                    &
     &                      Irec, MyType, Vsize,                        &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
#  endif
#  ifdef CHECKSUM
     &                      OCEAN(ng) % ad_v(:,:,:,Iinp),               &
     &                      checksum = Fhash)
#  else
     &                      OCEAN(ng) % ad_v(:,:,:,Iinp))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idVtlf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idVtlf)), Irec,            &
#  ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#  else
     &                          Fmin, Fmax
#  endif
            END IF
          END IF
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
              WRITE (stdout,60) TRIM(Vname(1,idVtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!  Process tracer type variables impulses.
!
        DO i=1,NT(ng)
          IF (find_string(var_name, n_var, Vname(1,idTtlf(i)),          &
     &                    INPvid)) THEN
            gtype=var_flag(INPvid)*r3dvar
            MyType=gtype
            status=nf_fread3d(ng, model, INPncname, INPncid,            &
     &                        Vname(1,idTtlf(i)), INPvid,               &
     &                        Irec, MyType, Vsize,                      &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, Fmin, Fmax,                        &
#  ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#  endif
#  ifdef CHECKSUM
     &                        OCEAN(ng) % ad_t(:,:,:,Iinp,i),           &
     &                        checksum = Fhash)
#  else
     &                        OCEAN(ng) % ad_t(:,:,:,Iinp,i))
#  endif
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idTtlf(i))), Irec,       &
     &                            TRIM(INPncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(1,idTtlf(i))), Irec,       &
#  ifdef CHECKSUM
     &                            Fmin, Fmax, Fhash
#  else
     &                            Fmin, Fmax
#  endif
              END IF
            END IF
!
            status=nf_fwrite3d(ng, model, TLF(ng)%ncid, TLF(ng)%Tid(i), &
     &                         Iout, MyType,                            &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#  ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#  endif
     &                         OCEAN(ng) % ad_t(:,:,:,Iinp,i))
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) TRIM(Vname(1,idTtlf(i))), Irec,       &
     &                            TRIM(TLF(ng)%name)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTtlf(i))),     &
     &                                    TRIM(INPncname)
            exit_flag=2
            RETURN
          END IF
        END DO
# endif
      END DO
!
!  Close input NetCDF file.
!
      CALL netcdf_close (ng, model, INPncid, INPncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,70) TRIM(INPncname)
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, model, INPncname, TLF(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_IMPULSE_NF90 - processing convolved adjoint',     &
     &        ' impulses, records: 1 to ',i0,/,21x,'file: ',a)
  20  FORMAT (/,' WRT_IMPULSE_NF90 - unable to open input NetCDF',      &
     &        ' file: ',a)
  30  FORMAT (/,' WRT_IMPULSE_NF90 - cannot find state variable: ',a,   &
     &        /,20x,'in input NetCDF file: ',a)
  40  FORMAT (/,' WRT_IMPULSE_NF90 - error while reading variable: ',a, &
     &        2x,'at time record = ',i0,                                &
     &        /,20x,'in input NetCDF file: ',a)
# ifdef CHECKSUM
  50  FORMAT (19x,'- ',a,/,22x,'(Rec = ',i0,' Min = ',1p,e15.8,         &
     &        ' Max = ',1p,e15.8,' CheckSum = ',i0,')')
# else
  50  FORMAT (19x,'- ',a,/,22x,'(Rec = ',i0,' Min = ',1p,e15.8,         &
     &        ' Max = ',1p,e15.8,')')
# endif
  60  FORMAT (/,' WRT_IMPULSE_NF90 - error while writing variable: ',a, &
     &        2x,'at time record = ',i0,/,20x,'into NetCDF file: ',a)
  70  FORMAT (/,' WRT_IMPULSE_NF90 - unable to close input NetCDF',     &
     &        ' file: ',a)
!
      RETURN
      END SUBROUTINE wrt_impulse_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_impulse_pio (ng, tile, model, INPncname,           &
     &                            LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
      character (len=*), intent(in) :: INPncname
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: Iinp, Iout, Irec, Nrec
      integer :: INPncid, INPvid
      integer :: i, status
      integer :: ibuffer(2), Vsize(4)
#  ifdef CHECKSUM
      integer(i8b) :: Fhash
#  endif
!
      real(r8) :: Fmin, Fmax
      real(dp) :: scale
      real(dp) :: inp_time(1)
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_impulse_pio"
!
      TYPE (IO_desc_t),  pointer :: ioDesc
      TYPE (File_desc_t) :: pioFile
      TYPE (My_VarDesc)  :: pioVar


#  include "set_bounds.h"
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
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
#  ifdef SP4DVAR
      IF (Master)  WRITE (stdout,10) Nrec, TRIM(TLF(ng)%name)
#  else
      IF (Master)  WRITE (stdout,10) Nrec-1, TRIM(TLF(ng)%name)
#  endif
!
!-----------------------------------------------------------------------
!  Read adjoint solution and convert to impulse forcing.  Then, write
!  impulse forcing into output NetCDF file.
!-----------------------------------------------------------------------
!
!  Open input NetCDF file.
!
      CALL pio_netcdf_open (ng, model, INPncname, 0, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,20) TRIM(INPncname)
        RETURN
      END IF
#  ifdef SP4DVAR
!
      Iinp=1
      scale=1.0_dp
!
      DO Irec=1,Nrec
        Iout=Irec
#  else
!
!  Process each record in input adjoint NetCDF except last. Note that
!  the adjoint records are processed backwards (Nrec-1:1) and written
!  in ascending time order (Iout initialized to 0) since the weak
!  constraint forcing will be read by the TL and RP models. Record
!  Nrec is not processed since it is not needed.
!
      Iinp=1
      Iout=0
      scale=1.0_dp
!
      DO Irec=Nrec-1,1,-1
        Iout=Iout+1
#  endif
!
!  Process time.
!
        IF (find_string(var_name, n_var, Vname(1,idtime), INPvid)) THEN
          CALL pio_netcdf_get_time (ng, model, INPncname,               &
     &                              Vname(1,idtime),                    &
     &                              Rclock%DateNumber, inp_time,        &
     &                              pioFile = pioFile,                  &
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
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
        END IF
!
!  Process free-surface weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idZtlf), INPvid)) THEN

          pioVar%vd=var_desc(INPvid)
          IF (KIND(OCEAN(ng)%ad_zeta).eq.8) THEN
            pioVar%dkind=PIO_double
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            pioVar%dkind=PIO_real
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          pioVar%gtype=r2dvar

          status=nf_fread2d(ng, model, INPncname, pioFile,              &
     &                      Vname(1,idZtlf), pioVar, Irec,              &
     &                      ioDesc, Vsize,                              &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
#  ifdef CHECKSUM
     &                      OCEAN(ng) % ad_zeta(:,:,Iinp),              &
     &                      checksum = Fhash)
#  else
     &                      OCEAN(ng) % ad_zeta(:,:,Iinp))
#  endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idZtlf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idZtlf)), Irec,            &
#  ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#  else
     &                          Fmin, Fmax
#  endif
            END IF
          END IF
!
          IF (TLF(ng)%pioVar(idZtlf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF

          status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idZtlf), Iout,              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % ad_zeta(:,:,Iinp))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idZtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idZtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF

#  ifndef SOLVE3D
!
!  Process 2D U-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idUbtf), INPvid)) THEN

          pioVar%vd=var_desc(INPvid)
          IF (KIND(OCEAN(ng)%ad_ubar).eq.8) THEN
            pioVar%dkind=PIO_double
            ioDesc => ioDesc_dp_u2dvar(ng)
          ELSE
            pioVar%dkind=PIO_real
            ioDesc => ioDesc_sp_u2dvar(ng)
          END IF
          pioVar%gtype=u2dvar

          status=nf_fread2d(ng, model, INPncname, pioFile,              &
     &                      Vname(1,idUbtf), pioVar, Irec,              &
     &                      ioDesc, Vsize,                              &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
#   endif
#   ifdef CHECKSUM
     &                      OCEAN(ng) % ad_ubar(:,:,Iinp),              &
     &                      checksum = Fhash)
#   else
     &                      OCEAN(ng) % ad_ubar(:,:,Iinp))
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idUbtf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idUbtf)), Irec,            &
#   ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#   else
     &                          Fmin, Fmax
#   endif
            END IF
          END IF
!
          IF (TLF(ng)%pioVar(idUbtf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u2dvar(ng)
          END IF

          status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idUbtf), Iout,              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#   endif
     &                       OCEAN(ng) % ad_ubar(:,:,Iinp))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idUbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!  Process 2D V-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idVbtf), INPvid)) THEN

          pioVar%vd=var_desc(INPvid)
          IF (KIND(OCEAN(ng)%ad_vbar).eq.8) THEN
            pioVar%dkind=PIO_double
            ioDesc => ioDesc_dp_v2dvar(ng)
          ELSE
            pioVar%dkind=PIO_real
            ioDesc => ioDesc_sp_v2dvar(ng)
          END IF
          pioVar%gtype=v2dvar

          status=nf_fread2d(ng, model, INPncname, pioFile,              &
     &                      Vname(1,idVbtf), pioVar, Irec,              &
     &                      ioDesc, Vsize,                              &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
#   endif
#   ifdef CHECKSUM
     &                      OCEAN(ng) % ad_vbar(:,:,Iinp),              &
     &                      checksum = Fhash)
#   else
     &                      OCEAN(ng) % ad_vbar(:,:,Iinp))
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idVbtf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idVbtf)), Irec,            &
#   ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#   else
     &                          Fmin, Fmax
#   endif
            END IF
          END IF
!
          IF (TLF(ng)%pioVar(idVbtf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v2dvar(ng)
          END IF

          status=nf_fwrite2d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idVbtf), Iout,              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       OCEAN(ng) % ad_vbar(:,:,Iinp))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idVbtf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVbtf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
!
!  Process 3D U-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idUtlf), INPvid)) THEN

          pioVar%vd=var_desc(INPvid)
          IF (KIND(OCEAN(ng)%ad_u).eq.8) THEN
            pioVar%dkind=PIO_double
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            pioVar%dkind=PIO_real
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF
          pioVar%gtype=u3dvar

          status=nf_fread3d(ng, model, INPncname, pioFile,              &
     &                      Vname(1,idUtlf), pioVar, Irec,              &
     &                      ioDesc, Vsize,                              &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      GRID(ng) % umask_full,                      &
#   endif
#   ifdef CHECKSUM
     &                      OCEAN(ng) % ad_u(:,:,:,Iinp),               &
     &                      checksum = Fhash)
#   else
     &                      OCEAN(ng) % ad_u(:,:,:,Iinp))
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idUtlf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idUtlf)), Irec,            &
#   ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#   else
     &                          Fmin, Fmax
#   endif
            END IF
          END IF
!
          IF (TLF(ng)%pioVar(idUtlf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_u3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_u3dvar(ng)
          END IF

          status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idUtlf), Iout,              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % umask_full,                     &
#   endif
     &                       OCEAN(ng) % ad_u(:,:,:,Iinp))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idUtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!  Process 3D V-momentum weak-constraint impulse forcing.
!
        IF (find_string(var_name, n_var, Vname(1,idVtlf), INPvid)) THEN

          pioVar%vd=var_desc(INPvid)
          IF (KIND(OCEAN(ng)%ad_v).eq.8) THEN
            pioVar%dkind=PIO_double
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            pioVar%dkind=PIO_real
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF
          pioVar%gtype=v3dvar

          status=nf_fread3d(ng, model, INPncname, pioFile,              &
     &                      Vname(1,idVtlf), pioVar, Irec,              &
     &                      ioDesc, Vsize,                              &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      GRID(ng) % vmask_full,                      &
#   endif
#   ifdef CHECKSUM
     &                      OCEAN(ng) % ad_v(:,:,:,Iinp),               &
     &                      checksum = Fhash)
#   else
     &                      OCEAN(ng) % ad_v(:,:,:,Iinp))
#   endif
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(1,idVtlf)), Irec,            &
     &                          TRIM(INPncname)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,50) TRIM(Vname(1,idVtlf)), Irec,            &
#   ifdef CHECKSUM
     &                          Fmin, Fmax, Fhash
#   else
     &                          Fmin, Fmax
#   endif
            END IF
          END IF
!
          IF (TLF(ng)%pioVar(idVtlf)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_v3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_v3dvar(ng)
          END IF

          status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,                &
     &                       TLF(ng)%pioVar(idVtlf), Iout,              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % vmask_full,                     &
#   endif
     &                       OCEAN(ng) % ad_v(:,:,:,Iinp))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,60) TRIM(Vname(1,idVtlf)), Irec,            &
     &                          TRIM(TLF(ng)%name)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        ELSE
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVtlf)),          &
     &                                  TRIM(INPncname)
          exit_flag=2
          RETURN
        END IF
!
!  Process tracer type variables impulses.
!
        DO i=1,NT(ng)
          IF (find_string(var_name, n_var, Vname(1,idTtlf(i)),          &
     &                    INPvid)) THEN

            pioVar%vd=var_desc(INPvid)
            IF (KIND(OCEAN(ng)%ad_t).eq.8) THEN
              pioVar%dkind=PIO_double
              ioDesc => ioDesc_dp_r3dvar(ng)
            ELSE
              pioVar%dkind=PIO_real
              ioDesc => ioDesc_sp_r3dvar(ng)
            END IF
            pioVar%gtype=r3dvar

            status=nf_fread3d(ng, model, INPncname, pioFile,            &
     &                        Vname(1,idTtlf(i)), pioVar,               &
     &                        Irec, ioDesc, Vsize,                      &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, Fmin, Fmax,                        &
#   ifdef MASKING
     &                        GRID(ng) % rmask,                         &
#   endif
#   ifdef CHECKSUM
     &                        OCEAN(ng) % ad_t(:,:,:,Iinp,i),           &
     &                        checksum = Fhash)
#   else
     &                        OCEAN(ng) % ad_t(:,:,:,Iinp,i))
#   endif
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,40) TRIM(Vname(1,idTtlf(i))), Irec,       &
     &                            TRIM(INPncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,50) TRIM(Vname(1,idTtlf(i))), Irec,       &
#   ifdef CHECKSUM
     &                            Fmin, Fmax, Fhash
#   else
     &                            Fmin, Fmax
#   endif
              END IF
            END IF
!
            IF (TLF(ng)%pioTrc(i)%dkind.eq.PIO_double) THEN
              ioDesc => ioDesc_dp_r3dvar(ng)
            ELSE
              ioDesc => ioDesc_sp_r3dvar(ng)
            END IF

            status=nf_fwrite3d(ng, model, TLF(ng)%pioFile,              &
     &                         TLF(ng)%pioTrc(i), Iout,                 &
     &                         ioDesc,                                  &
     &                         LBi, UBi, LBj, UBj, 1, N(ng), scale,     &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         OCEAN(ng) % ad_t(:,:,:,Iinp,i))
            IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,60) TRIM(Vname(1,idTtlf(i))), Irec,       &
     &                            TRIM(TLF(ng)%name)
              END IF
              exit_flag=3
              ioerror=status
              RETURN
            END IF
          ELSE
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTtlf(i))),     &
     &                                    TRIM(INPncname)
            exit_flag=2
            RETURN
          END IF
        END DO
#  endif
      END DO
!
!  Close input NetCDF file.
!
      CALL pio_netcdf_close (ng, model, pioFile, INPncname, .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
        WRITE (stdout,70) TRIM(INPncname)
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize impulse NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, model, INPncname, TLF(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (2x,'WRT_IMPULSE_PIO  - processing convolved adjoint',     &
     &        ' impulses, records: 1 to ',i0,/,21x,'file: ',a)
  20  FORMAT (/,' WRT_IMPULSE_PIO - unable to open input NetCDF',       &
     &        ' file: ',a)
  30  FORMAT (/,' WRT_IMPULSE_PIO - cannot find state variable: ',a,    &
     &        /,20x,'in input NetCDF file: ',a)
  40  FORMAT (/,' WRT_IMPULSE_PIO - error while reading variable: ',a,  &
     &        2x,'at time record = ',i0,                                &
     &        /,20x,'in input NetCDF file: ',a)
#  ifdef CHECKSUM
  50  FORMAT (19x,'- ',a,/,22x,'(Rec = ',i0,' Min = ',1p,e15.8,         &
     &        ' Max = ',1p,e15.8,' CheckSum = ',i0,')')
#  else
  50  FORMAT (19x,'- ',a,/,22x,'(Rec = ',i0,' Min = ',1p,e15.8,         &
     &        ' Max = ',1p,e15.8,')')
#  endif
  60  FORMAT (/,' WRT_IMPULSE_PIO - error while writing variable: ',a,  &
     &        2x,'at time record = ',i0,/,20x,'into NetCDF file: ',a)
  70  FORMAT (/,' WRT_IMPULSE_PIO - unable to close input NetCDF',      &
     &        ' file: ',a)
!
      RETURN
      END SUBROUTINE wrt_impulse_pio
# endif
#endif
      END MODULE wrt_impulse_mod

