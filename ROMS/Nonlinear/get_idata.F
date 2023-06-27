#include "cppdefs.h"
      SUBROUTINE get_idata (ng)
!
!git $Id$
!svn $Id: get_idata.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads input data that needs to be obtained only once.  !
!                                                                      !
!  Currently,  this routine is only executed in serial mode by the     !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef RED_TIDE
      USE mod_biology
#endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
#endif
#ifdef RED_TIDE
      USE mod_ocean
#endif
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
#if defined SSH_TIDES || defined UV_TIDES
      USE mod_tides
#endif
!
#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
      USE def_tides_mod,  ONLY : def_tides
#endif
      USE nf_fread3d_mod, ONLY : nf_fread3d
#ifdef SOLVE3D
      USE nf_fread4d_mod, ONLY : nf_fread4d
#endif
      USE strings_mod,    ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical, save :: recordless = .TRUE.

      logical, dimension(3) :: update =                                 &
     &         (/ .FALSE., .FALSE., .FALSE. /)
!
      integer :: LBi, UBi, LBj, UBj
      integer :: itrc, is

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
      integer :: gtype, status, varid, Vsize(4)
!
      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmin, Fmax, Htime
#endif
      real(r8) :: time_save = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
# if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (IO_Desc_t),  pointer :: ioDesc
      TYPE (My_VarDesc), pointer :: pioVar
# endif
#endif
!
      SourceFile=MyFile
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO is=1,4
        Vsize(is)=0
      END DO
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3, __LINE__, MyFile)
#endif

#if defined SSH_TIDES || defined UV_TIDES
!
!-----------------------------------------------------------------------
!  Tide period, amplitude, phase, and currents.
!-----------------------------------------------------------------------
!
!  Tidal Period.
!
      IF (LprocessTides(ng)) THEN
        IF (iic(ng).eq.0) THEN
          CALL get_ngfld (ng, iNLM, idTper, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), recordless, update(1),           &
     &                    1, MTC, 1, 1, 1, NTC(ng), 1,                  &
     &                    TIDES(ng) % Tperiod)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END IF
#endif

#ifdef SSH_TIDES
!
!  Tidal elevation amplitude and phase. In order to read data as a
!  function of tidal period, we need to reset the model time variables
!  temporarily.
!
      IF (LprocessTides(ng)) THEN
        IF (iic(ng).eq.0) THEN
          time_save=time(ng)
          time(ng)=8640000.0_r8
          tdays(ng)=time(ng)*sec2day
!
          CALL get_2dfld (ng, iNLM, idTzam, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % SSH_Tamp)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL get_2dfld (ng, iNLM, idTzph, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % SSH_Tphase)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          time(ng)=time_save
          tdays(ng)=time(ng)*sec2day
        END IF
      END IF
#endif

#ifdef UV_TIDES
!
!  Tidal currents angle, phase, major and minor ellipse axis.
!
      IF (LprocessTides(ng)) THEN
        IF (iic(ng).eq.0) THEN
          time_save=time(ng)
          time(ng)=8640000.0_r8
          tdays(ng)=time(ng)*sec2day
!
          CALL get_2dfld (ng, iNLM, idTvan, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tangle)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL get_2dfld (ng, iNLM, idTvph, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tphase)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL get_2dfld (ng, iNLM, idTvma, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tmajor)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL get_2dfld (ng, iNLM, idTvmi, TIDE(ng)%ncid,              &
# if defined PIO_LIB && defined DISTRIBUTE
     &                    TIDE(ng)%pioFile,                             &
# endif
     &                    1, TIDE(ng), update(1),                       &
     &                    LBi, UBi, LBj, UBj, MTC, NTC(ng),             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    TIDES(ng) % UV_Tminor)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          time(ng)=time_save
          tdays(ng)=time(ng)*sec2day
        END IF
      END IF
#endif

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
!
!-----------------------------------------------------------------------
!  If detiding and applicable, define additional variable to store
!  time-accumulated tide harmonics variables.  This variable are
!  defined and written into input tide forcing NetCDF file.
!-----------------------------------------------------------------------
!
      CALL def_tides (ng, LdefTIDE(ng))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  If restarting, read in time-accumulated tide harmonics variables.
!-----------------------------------------------------------------------
!
      IF (.not.LdefTIDE(ng).and.(nrrec(ng).ne.0)) THEN
!
!  For consistency, check time of written accumulate harmonics and
!  compare to current time.
!
        SELECT CASE (HAR(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_time (ng, iNLM, HAR(ng)%name,               &
     &                            Vname(1,idtime),                      &
     &                            Rclock%DateNumber, Htime,             &
     &                            ncid = HAR(ng)%ncid)

# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_time (ng, iNLM, HAR(ng)%name,           &
     &                                Vname(1,idtime),                  &
     &                                Rclock%DateNumber, Htime,         &
     &                                pioFile = HAR(ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        IF (time(ng).ne.Htime) THEN
          IF (Master) THEN
            WRITE (stdout,20) tdays(ng), Htime*sec2day
          END IF
          exit_flag=2
          ioerror=0
          RETURN
        END IF
!
!  Number of time-acummulate tide harmonics.
!
        SELECT CASE (HAR(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_ivar (ng, iNLM, HAR(ng)%name,               &
     &                            'Hcount', Hcount(ng),                 &
     &                            ncid = HAR(ng)%ncid)

# if defined PIO_LIB && defined DISTRIBUT
          CASE (io_pio)
            CALL pio_netcdf_get_ivar (ng, iNLM, HAR(ng)%name,           &
     &                                'Hcount', Hcount(ng),             &
     &                                pioFile = HAR(ng)%pioFile)
# endif
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time-accumulated COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idCosW, HAR(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  HAR(ng)%pioFile,                                &
# endif
     &                  1, HAR(ng), recordless, update(1),              &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % CosW_sum)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time-accumulated SIN(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSinW, HAR(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  HAR(ng)%pioFile,                                &
# endif
     &                  1, HAR(ng), recordless, update(1),              &
     &                  1, MTC, 1, 1, 1, NTC(ng), 1,                    &
     &                  TIDES(ng) % SinW_sum)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time-accumulated COS(omega(k)*t)*COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idCos2, HAR(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  HAR(ng)%pioFile,                                &
# endif
     &                  1, HAR(ng), recordless, update(1),              &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) % CosWCosW)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time-accumulated SIN(omega(k)*t)*SIN(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSin2, HAR(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  HAR(ng)%pioFile,                                &
# endif
     &                  1, HAR(ng), recordless, update(1),              &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) %  SinWSinW)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time-accumulated SIN(omega(k)*t)*COS(omega(k)*t) harmonics.
!
        CALL get_ngfld (ng, iNLM, idSWCW, HAR(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  HAR(ng)%pioFile,                                &
# endif
     &                  1, HAR(ng), recordless, update(1),              &
     &                  1, MTC, MTC, 1, 1, NTC(ng), NTC(ng),            &
     &                  TIDES(ng) %  SinWCosW)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Time-accumulated free-surface tide harmonics.
!
        IF (Aout(idFsuD,ng)) THEN
          SELECT CASE (HAR(ng)%IOtype)
            CASE (io_nf90)
              gtype=r3dvar
              status=nf_fread3d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%ncid,                           &
     &                          Vname(1,idFsuH),                        &
     &                          HAR(ng)%Vid(idFsuH),                    &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj, 0, 2*NTC(ng),       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % rmask,                       &
# endif
     &                          TIDES(ng) % zeta_tide)

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)

              pioVar%vd => HAR(ng)%pioVar(idFsuH)%vd
              IF (KIND(TIDES(ng)%zeta_tide).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_r2dhar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_r2dhar(ng)
              END IF
              pioVar%gtype=r2dvar
!
              status=nf_fread3d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%pioFile,                        &
     &                          Vname(1,idFsuH),                        &
     &                          pioVar,                                 &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj, 0, 2*NTC(ng),       &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % rmask,                       &
#  endif
     &                          TIDES(ng) % zeta_tide)

# endif
          END SELECT
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idFsuH)),                  &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idFsuH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated 2D u-momentum tide harmonics.
!
        IF (Aout(idu2dD,ng)) THEN
          SELECT CASE (HAR(ng)%IOtype)
            CASE (io_nf90)
              gtype=u3dvar
              status=nf_fread3d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%ncid,                           &
     &                          Vname(1,idu2dH),                        &
     &                          HAR(ng)%Vid(idu2dH),                    &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj, 0, 2*NTC(ng),       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % umask,                       &
# endif
     &                          TIDES(ng) % ubar_tide)

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)

              pioVar%vd => HAR(ng)%pioVar(idu2dH)%vd
              IF (KIND(TIDES(ng)%ubar_tide).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u2dhar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u2dhar(ng)
              END IF
              pioVar%gtype=u2dvar
!
              status=nf_fread3d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%pioFile,                        &
     &                          Vname(1,idu2dH),                        &
     &                          pioVar,                                 &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj, 0, 2*NTC(ng),       &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % umask,                       &
#  endif
     &                          TIDES(ng) % ubar_tide)
# endif
          END SELECT
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idu2dH)),                  &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idu2dH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated 2D v-momentum tide harmonics.
!
        IF (Aout(idv2dD,ng)) THEN
          SELECT CASE (HAR(ng)%IOtype)
            CASE (io_nf90)
              gtype=v3dvar
              status=nf_fread3d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%ncid,                           &
     &                          Vname(1,idv2dH),                        &
     &                          HAR(ng)%Vid(idv2dH),                    &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj, 0, 2*NTC(ng),       &
     &                          Fscl, Fmin, Fmax,                       &
# ifdef MASKING
     &                          GRID(ng) % vmask,                       &
# endif
     &                          TIDES(ng) % vbar_tide)

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)

              pioVar%vd => HAR(ng)%pioVar(idv2dH)%vd
              IF (KIND(TIDES(ng)%vbar_tide).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v2dhar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v2dhar(ng)
              END IF
              pioVar%gtype=v2dvar
!
              status=nf_fread3d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%pioFile,                        &
     &                          Vname(1,idv2dH),                        &
     &                          pioVar,                                 &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj, 0, 2*NTC(ng),       &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#  endif
     &                          TIDES(ng) % vbar_tide)
# endif
          END SELECT
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idv2dH)),                  &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idv2dH)), Fmin, Fmax
            END IF
          END IF
        END IF

# ifdef SOLVE3D
!
!  Time-accumulated 3D u-momentum tide harmonics.
!
        IF (Aout(idu3dD,ng)) THEN
          SELECT CASE (HAR(ng)%IOtype)
            CASE (io_nf90)
              gtype=u3dvar
              status=nf_fread4d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%ncid,                           &
     &                          Vname(1,idu3dH),                        &
     &                          HAR(ng)%Vid(idu3dH),                    &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          0, 2*NTC(ng),                           &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % umask,                       &
#  endif
     &                          TIDES(ng) % u_tide)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)

              pioVar%vd => HAR(ng)%pioVar(idu3dH)%vd
              IF (KIND(TIDES(ng)%u_tide).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_u3dhar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_u3dhar(ng)
              END IF
              pioVar%gtype=u3dvar
!
              status=nf_fread4d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%pioFile,                        &
     &                          Vname(1,idu3dH),                        &
     &                          pioVar,                                 &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          0, 2*NTC(ng),                           &
     &                          Fscl, Fmin, Fmax,                       &
#   ifdef MASKING
     &                          GRID(ng) % umask,                       &
#   endif
     &                          TIDES(ng) % u_tide)
#  endif
          END SELECT
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idu3dH)),                  &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idu3dH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated 3D v-momentum tide harmonics.
!
        IF (Aout(idv3dD,ng)) THEN
          SELECT CASE (HAR(ng)%IOtype)
            CASE (io_nf90)
              gtype=v3dvar
              status=nf_fread4d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%ncid,                           &
     &                          Vname(1,idv3dH),                        &
     &                          HAR(ng)%Vid(idv3dH),                    &
     &                          0, gtype, Vsize,                        &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          0, 2*NTC(ng),                           &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#  endif
     &                          TIDES(ng) % v_tide)

#  if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)

              pioVar%vd => HAR(ng)%pioVar(idv3dH)%vd
              IF (KIND(TIDES(ng)%v_tide).eq.8) THEN
                pioVar%dkind=PIO_double
                ioDesc => ioDesc_dp_v3dhar(ng)
              ELSE
                pioVar%dkind=PIO_real
                ioDesc => ioDesc_sp_v3dhar(ng)
              END IF
              pioVar%gtype=v3dvar
!
              status=nf_fread4d(ng, iNLM, HAR(ng)%name,                 &
     &                          HAR(ng)%pioFile,                        &
     &                          Vname(1,idv3dH),                        &
     &                          pioVar,                                 &
     &                          0, ioDesc, Vsize,                       &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          0, 2*NTC(ng),                           &
     &                          Fscl, Fmin, Fmax,                       &
#  ifdef MASKING
     &                          GRID(ng) % vmask,                       &
#  endif
     &                          TIDES(ng) % v_tide)
#  endif
          END SELECT
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idv3dH)),                  &
     &                          TRIM(HAR(ng)%name)
            END IF
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idv3dH)), Fmin, Fmax
            END IF
          END IF
        END IF
!
!  Time-accumulated temperature and salinity tide harmonics.
!
        DO itrc=1,NAT
          IF (Aout(idTrcD(itrc),ng)) THEN
            SELECT CASE (HAR(ng)%IOtype)
              CASE (io_nf90)
                gtype=r3dvar
                status=nf_fread4d(ng, iNLM, HAR(ng)%name,               &
     &                            HAR(ng)%ncid,                         &
     &                            Vname(1,idTrcH(itrc)),                &
     &                            HAR(ng)%Vid(idTrcH(itrc)),            &
     &                            0, gtype, Vsize,                      &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            0, 2*NTC(ng),                         &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#  endif
     &                            TIDES(ng) % t_tide(:,:,:,:,itrc))

#  if defined PIO_LIB && defined DISTRIBUTE
              CASE (io_pio)

                pioVar%vd => HAR(ng)%pioVar(idTrcH(itrc))%vd
                IF (KIND(TIDES(ng)%t_tide).eq.8) THEN
                  pioVar%dkind=PIO_double
                  ioDesc => ioDesc_dp_r3dhar(ng)
                ELSE
                  pioVar%dkind=PIO_real
                  ioDesc => ioDesc_sp_r3dhar(ng)
                END IF
                pioVar%gtype=r3dvar
!
                status=nf_fread4d(ng, iNLM, HAR(ng)%name,               &
     &                            HAR(ng)%pioFile,                      &
     &                            Vname(1,idTrcH(itrc)),                &
     &                            pioVar,                               &
     &                            0, ioDesc, Vsize,                     &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            0, 2*NTC(ng),                         &
     &                            Fscl, Fmin, Fmax,                     &
#  ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#  endif
     &                            TIDES(ng) % t_tide(:,:,:,:,itrc))
#  endif
            END SELECT
            IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vname(1,idTrcH(itrc))),          &
     &                            TRIM(HAR(ng)%name)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,30) TRIM(Vname(2,idTrcH(itrc))),          &
     &                            Fmin, Fmax
              END IF
            END IF
          END IF
        END DO
# endif
      END IF
#endif

#ifndef ANA_PSOURCE
!
!-----------------------------------------------------------------------
!  Read in point Sources/Sinks position, direction, special flag, and
!  mass transport nondimensional shape profile.  Point sources are at
!  U- and V-points.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).eq.0).and.                                           &
     &    (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng)))) THEN
        CALL get_ngfld (ng, iNLM, idRxpo, SSF(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  SSF(ng)%pioFile,                                &
# endif
     &                  1, SSF(ng), recordless, update(1),              &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Xsrc)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_ngfld (ng, iNLM, idRepo, SSF(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  SSF(ng)%pioFile,                                &
# endif
     &                  1, SSF(ng), recordless, update(1),              &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Ysrc)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        CALL get_ngfld (ng, iNLM, idRdir, SSF(ng)%ncid,                 &
# if defined PIO_LIB && defined DISTRIBUTE
     &                  SSF(ng)%pioFile,                                &
# endif
     &                  1, SSF(ng), recordless, update(1),              &
     &                  1, Nsrc(ng), 1, 1, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % Dsrc)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SOLVE3D
        CALL get_ngfld (ng, iNLM, idRvsh, SSF(ng)%ncid,                 &
#  if defined PIO_LIB && defined DISTRIBUTE
     &                  SSF(ng)%pioFile,                                &
#  endif
     &                  1, SSF(ng), recordless, update(1),              &
     &                  1, Nsrc(ng), N(ng), 1, 1, Nsrc(ng), N(ng),      &
     &                  SOURCES(ng) % Qshape)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        DO is=1,Nsrc(ng)
          SOURCES(ng)%Isrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Xsrc(is)),Lm(ng)+1))
          SOURCES(ng)%Jsrc(is)=                                         &
     &                MAX(1,MIN(NINT(SOURCES(ng)%Ysrc(is)),Mm(ng)+1))
        END DO
      END IF
#endif

#ifdef RED_TIDE
!
!-----------------------------------------------------------------------
!  Red tides.
!-----------------------------------------------------------------------
!
!  Read initial dynoflagellate bottom cyst concentration.
!
      CALL get_2dfld (ng, iNLM, idCyst, FRCncid(idCyst,ng),             &
# if defined PIO_LIB && defined DISTRIBUTE
     &                FRCpioFile(idCyst,ng),                            &
# endif
     &                nFfiles(ng), FRC(1,ng), update(1),                &
     &                LBi, UBi, LBj, UBj, 1, 1,                         &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                OCEAN(ng) % CystIni)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      IF (.not.Linfo(1,idCyst,ng)) THEN             ! if not gridded,
        OCEAN(ng) % CystIni = Fpoint(1,idCyst,ng)   ! load point data
      END IF
#endif

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3, __LINE__, MyFile)
#endif

#if defined AVERAGES  && defined AVERAGES_DETIDE && \
   (defined SSH_TIDES || defined UV_TIDES)
!
  10  FORMAT (/,' GET_IDATA - error while reading variable: ',a,        &
     &        /,13x,'in input NetCDF file: ',a)
  20  FORMAT (/,' GET_IDATA - incosistent restart and harmonics time:', &
     &        /,13x,f15.4,2x,f15.4)
  30  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
#endif
!
      RETURN
      END SUBROUTINE get_idata
