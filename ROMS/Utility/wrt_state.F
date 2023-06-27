#include "cppdefs.h"
      MODULE wrt_state_mod
#if defined ADJOINT && defined TANGENT
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module writes requested tangent linear or adjoint state        !
!  fields using the standard NetCDF library or the Parallel-IO (PIO)   !
!  library into an output file described by IO structure S. The NetCDF !
!  is opened and closed, and file ID or descriptor is not updated.     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     label        Identification label (string)                       !
!     OutRec       NetCDF file time record (integer)                   !
!     i2d          2D state variables time level index (integer)       !
!     i3d          3D state variables time level index (integer)       !
!     stime        Time variable value (real; seconds)                 !
!     S            File derived type structure (TYPE T_IO)             !
!     CloseFile    Switch to close NetCDF file (logical, OPTIONAL)     !
!                    (defaut = .TRUE.)                                 !
!  On Output:                                                          !
!                                                                      !
!     ad_arrays    ADM state arrays if model = iADM                    !
!              or,                                                     !
!     tl_arrays    TLM state arrays if model = iTLM or iRPM            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef ADJUST_BOUNDARY
      USE mod_boundary
# endif
# ifdef SOLVE3D
      USE mod_coupling
# endif
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      USE mod_forces
# endif
      USE mod_grid
      USE mod_iounits
# ifdef SOLVE3D
      USE mod_mixing
# endif
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE dateclock_mod,       ONLY : time_string
      USE nf_fwrite2d_mod,     ONLY : nf_fwrite2d
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
# endif
# ifdef SOLVE3D
      USE nf_fwrite3d_mod,     ONLY : nf_fwrite3d
#  ifdef ADJUST_BOUNDARY
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
#  endif
# endif
      USE strings_mod,         ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: wrt_state
      PRIVATE :: wrt_state_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: wrt_state_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE wrt_state (ng, tile, model, label,                     &
     &                      OutRec, i2d, i3d, stime, S)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, OutRec, i2d, i3d
!
      real(dp), intent(in) :: stime
!
      character (len=*), intent(in) :: label
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
# ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
# endif
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Write out history fields according to IO type.
!-----------------------------------------------------------------------
!
# ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
# endif
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (HIS(ng)%IOtype)
        CASE (io_nf90)
          CALL wrt_state_nf90 (ng, tile, model, label,                  &
     &                         OutRec, i2d, i3d, stime,                 &
# ifdef ADJUST_BOUNDARY
     &                         LBij, UBij,                              &
     &                         LBi, UBi, LBj, UBj,                      &
# endif
     &                         S)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL wrt_state_pio (ng, tile, model, label,                   &
     &                        OutRec, i2d, i3d, stime,                  &
#  ifdef ADJUST_BOUNDARY
     &                        LBij, UBij,                               &
     &                        LBi, UBi, LBj, UBj,                       &
#  endif
     &                        S)
# endif
        CASE DEFAULT
          IF (Master) THEN
            WRITE (stdout,10) S(ng)%IOtype
          END IF
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' WRT_STATE - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE wrt_state
!
!***********************************************************************
      SUBROUTINE wrt_state_nf90 (ng, tile, model, label,                &
     &                           OutRec, i2d, i3d, stime,               &
# ifdef ADJUST_BOUNDARY
     &                           LBij, UBij,                            &
     &                           LBi, UBi, LBj, UBj,                    &
# endif
     &                           S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, OutRec, i2d, i3d
# ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
# endif
      integer :: LBi, UBi, LBj, UBj
!
      real(dp), intent(in) :: stime
!
      character (len=*), intent(in) :: label
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      integer :: Sstr, Send
      integer :: Fcount, gfactor, gtype, ncid, omode, status
# ifdef SOLVE3D
      integer :: i, itrc, j, k
# endif
!
      real(r8) :: Fmin, Fmax
      real(dp) :: scale
      real(dp) :: Tval(1)
!
      character (len=15) :: Tstring
      character (len=22) :: t_code
      character (len=50) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_state_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out tangent linear or adjoint state fields.
!-----------------------------------------------------------------------
# ifdef PROFILE
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, model, 81, __LINE__, MyFile)
# endif
!
!  Open output NetCDF file.
!
      omode=1                                    ! read and write access
      CALL netcdf_open (ng, model, S(ng)%name, omode, ncid)
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
      S(ng)%Rindex=S(ng)%Rindex+1
      Fcount=S(ng)%Fcount
      S(ng)%Nrec(Fcount)=S(ng)%Nrec(Fcount)+1
!
!  Report information.
!
      IF (Master) THEN
        CALL time_string (stime, t_code)
        IF (model.eq.iNLM) THEN
          string='writing NLM state fields,'
        ELSE IF (model.eq.iTLM) THEN
          string='writing TLM state fields,'
        ELSE IF (model.eq.iADM) THEN
          string='writing ADM state fields,'
        ELSE IF (model.eq.iRPM) THEN
          string='writing RPM state fields,'
        END IF
        Sstr=SCAN(CalledFrom,'/',BACK=.TRUE.)+1
        Send=LEN_TRIM(CalledFrom)
        WRITE (Tstring,'(f15.4)') stime*sec2day
        WRITE (stdout,10) TRIM(label), TRIM(string), t_code,            &
     &                    ng, TRIM(ADJUSTL(Tstring)), TRIM(S(ng)%name), &
     &                    i3d, OutRec, CalledFrom(Sstr:Send)
      END IF
!
!  Write out model time (s).
!
      IF (LwrtPER(ng)) THEN
        Tval(1)=REAL(OutRec,dp)*day2sec
      ELSE
        Tval(1)=stime
      END IF
      CALL netcdf_put_fvar (ng, model, S(ng)%name,                      &
     &                      TRIM(Vname(1,idtime)), tval,                &
     &                      (/OutRec/), (/1/),                          &
     &                      ncid = ncid,                                &
     &                      varid = S(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface.
!
      scale=1.0_dp
      gtype=gfactor*r2dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite2d(ng, model, ncid, S(ng)%Vid(idFsur),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     OCEAN(ng) % tl_zeta(:,:,i2d),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
# else
     &                     OCEAN(ng) % tl_zeta(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
# endif
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite2d(ng, model, ncid, S(ng)%Vid(idFsur),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                            &
# endif
# ifdef WET_DRY
     &                     OCEAN(ng) % ad_zeta(:,:,i2d),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
# else
     &                     OCEAN(ng) % ad_zeta(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
# endif
      END IF
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        IF (model.eq.iTLM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isFsur)),               &
     &                           S(ng)%Vid(idSbry(isFsur)),             &
     &                           OutRec, r2dvar,                        &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isFsur)),               &
     &                           S(ng)%Vid(idSbry(isFsur)),             &
     &                           OutRec, r2dvar,                        &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idSbry(isFsur))),            &
     &                        Fmin, Fmax
          END IF
        END IF
      END IF
# endif
!
!  Write out 2D U-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*u2dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite2d(ng, model, ncid, S(ng)%Vid(idUbar),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     OCEAN(ng) % tl_ubar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)

      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite2d(ng, model, ncid, S(ng)%Vid(idUbar),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                            &
# endif
     &                     OCEAN(ng) % ad_ubar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        IF (model.eq.iTLM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isUbar)),               &
     &                           S(ng)%Vid(idSbry(isUbar)),             &
     &                           OutRec, u2dvar,                        &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isUbar)),               &
     &                           S(ng)%Vid(idSbry(isUbar)),             &
     &                           OutRec, u2dvar,                        &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idSbry(isUbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
      END IF
# endif
!
!  Write out 2D V-momentum component.
!
      scale=1.0_dp
      gtype=gfactor*v2dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite2d(ng, model, ncid, S(ng)%Vid(idVbar),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     OCEAN(ng) % tl_vbar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite2d(ng, model, ncid, S(ng)%Vid(idVbar),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                            &
# endif
     &                     OCEAN(ng) % ad_vbar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        IF (model.eq.iTLM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isVbar)),               &
     &                           S(ng)%Vid(idSbry(isVbar)),             &
     &                           OutRec, v2dvar,                        &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isVbar)),               &
     &                           S(ng)%Vid(idSbry(isVbar)),             &
     &                           OutRec, v2dvar,                        &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idSbry(isVbar))),            &
     &                        Fmin, Fmax
          END IF
        END IF
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idUsms),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % tl_ustr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idUsms),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     FORCES(ng) % ad_ustr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF

      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), TRIM(label),         &
     &                      TRIM(S(ng)%name), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUsms)), Fmin, Fmax
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idVsms),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % tl_vstr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)

      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idVsms),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     FORCES(ng) % ad_vstr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), TRIM(label),         &
     &                      TRIM(S(ng)%name), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVsms)), Fmin, Fmax
        END IF
      END IF
# endif

# ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*u3dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idUvel),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % tl_u(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)

      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idUvel),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % ad_u(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        IF (model.eq.iTLM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isUvel)),               &
     &                           S(ng)%Vid(idSbry(isUvel)),             &
     &                           OutRec, u3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isUvel)),               &
     &                           S(ng)%Vid(idSbry(isUvel)),             &
     &                           OutRec, u3dvar,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idSbry(isUvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Write out 3D V-momentum component (m/s).
!
      scale=1.0_dp
      gtype=gfactor*v3dvar
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idVvel),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % tl_v(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, ncid, S(ng)%Vid(idVvel),          &
     &                     OutRec, gtype,                               &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % ad_v(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        IF (model.eq.iTLM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isVvel)),               &
     &                           S(ng)%Vid(idSbry(isVvel)),             &
     &                           OutRec , v3dvar,                       &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, ncid,           &
     &                           Vname(1,idSbry(isVvel)),               &
     &                           S(ng)%Vid(idSbry(isVvel)),             &
     &                           OutRec , v3dvar,                       &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idSbry(isVvel))),            &
     &                        Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        gtype=gfactor*r3dvar
        IF (model.eq.iTLM) THEN
          status=nf_fwrite3d(ng, model, ncid, S(ng)%Tid(itrc),          &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % tl_t(:,:,:,i3d,itrc),          &
     &                       MinValue = Fmin,                           &
     &                       MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite3d(ng, model, ncid, S(ng)%Tid(itrc),          &
     &                       OutRec, gtype,                             &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       OCEAN(ng) % ad_t(:,:,:,i3d,itrc),          &
     &                       MinValue = Fmin,                           &
     &                       MaxValue = Fmax)
        END IF
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),              &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idTvar(itrc))),              &
     &                        Fmin, Fmax
          END IF
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          IF (model.eq.iTLM) THEN
            status=nf_fwrite3d_bry (ng, model, S(ng)%name, ncid,        &
     &                              Vname(1,idSbry(isTvar(itrc))),      &
     &                              S(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                              OutRec, r3dvar,                     &
     &                              LBij, UBij, 1, N(ng), Nbrec(ng),    &
     &                              scale,                              &
     &                              BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,&
     &                                                Lbout(ng),itrc),  &
     &                              MinValue = Fmin,                    &
     &                              MaxValue = Fmax)
          ELSE IF (model.eq.iADM) THEN
            status=nf_fwrite3d_bry (ng, model, S(ng)%name, ncid,        &
     &                              Vname(1,idSbry(isTvar(itrc))),      &
     &                              S(ng)%Vid(idSbry(isTvar(itrc))),    &
     &                              OutRec, r3dvar,                     &
     &                              LBij, UBij, 1, N(ng), Nbrec(ng),    &
     &                              scale,                              &
     &                              BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,&
     &                                                Lbout(ng),itrc),  &
     &                              MinValue = Fmin,                    &
     &                              MaxValue = Fmax)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          TRIM(label), TRIM(S(ng)%name), OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idSbry(isTvar(itrc)))),    &
     &                          Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp                      ! kinematic flux units
          gtype=gfactor*r3dvar
          IF (model.eq.iTLM) THEN
            status=nf_fwrite3d(ng, iTLM, ncid,                          &
     &                         S(ng)%Vid(idTsur(itrc)),                 &
     &                         OutRec, gtype,                           &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale, &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         FORCES(ng)% tl_tflux(:,:,:,              &
     &                                              Lfout(ng),itrc),    &
     &                         MinValue = Fmin,                         &
     &                         MaxValue = Fmax)
          ELSE IF (model.eq.iADM) THEN
            status=nf_fwrite3d(ng, iTLM, ncid,                          &
     &                         S(ng)%Vid(idTsur(itrc)),                 &
     &                         OutRec, gtype,                           &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale, &
#   ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#   endif
     &                         FORCES(ng)% ad_tflux(:,:,:,              &
     &                                              Lfout(ng),itrc),    &
     &                         MinValue = Fmin,                         &
     &                         MaxValue = Fmax)
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,idTsur(itrc))),            &
     &                          TRIM(label), TRIM(S(ng)%name), Lfout(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,idTsur(itrc))),            &
     &                          Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent NetCDF file to disk to allow other processes
!  to access data immediately after it is written. Then, close NetCDF
!  file.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iTLM, S(ng)%name, ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL netcdf_close (ng, iTLM, ncid, S(ng)%name, .FALSE.)

# ifdef PROFILE
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, model, 81, __LINE__, MyFile)
# endif
!
  10  FORMAT (1x,a,': WRT_STATE_NF90 - ',a,t75,a,                       &
     &        /,24x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Index=',i1,', Rec=',i0,')',                            &
     &        /,24x,'Called from ''',a,'''')
  20  FORMAT (/,' WRT_STATE_NF90 - error while writing variable: ',a,   &
     &        /,18x,'into ',a,' NetCDF file: ',a,                       &
     &        /,18x,'for time record: ',i0)
  30  FORMAT (21x,'- ',a,/,24x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE wrt_state_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE wrt_state_pio (ng, tile, model, label,                 &
     &                          OutRec, i2d, i3d, stime,                &
#  ifdef ADJUST_BOUNDARY
     &                          LBij, UBij,                             &
     &                          LBi, UBi, LBj, UBj,                     &
#  endif
     &                          S)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, OutRec, i2d, i3d
#  ifdef ADJUST_BOUNDARY
      integer, intent(in) :: LBij, UBij
#  endif
      integer :: LBi, UBi, LBj, UBj
!
      real(dp), intent(in) :: stime
!
      character (len=*), intent(in) :: label
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      integer :: Sstr, Send
      integer :: Fcount, ifield, omode, status
#  ifdef SOLVE3D
      integer :: i, itrc, j, k
#  endif
!
      real(r8) :: Fmin, Fmax
      real(dp) :: scale
      real(dp) :: Tval(1)
!
      character (len=15) :: Tstring
      character (len=22) :: t_code
      character (len=50) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", wrt_state_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
      TYPE (file_desc_t) :: pioFile
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out tangent linear or adjoint state fields.
!-----------------------------------------------------------------------
#  ifdef PROFILE
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, model, 81, __LINE__, MyFile)
#  endif
!
!  Open output NetCDF file.
!
      omode=1                                    ! read and write access
      CALL pio_netcdf_open (ng, model, S(ng)%name, omode, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set time record index.
!
      S(ng)%Rindex=S(ng)%Rindex+1
      Fcount=S(ng)%Fcount
      S(ng)%Nrec(Fcount)=S(ng)%Nrec(Fcount)+1
!
!  Report information.
!
      IF (Master) THEN
        CALL time_string (stime, t_code)
        IF (model.eq.iNLM) THEN
          string='writing NLM state fields,'
        ELSE IF (model.eq.iTLM) THEN
          string='writing TLM state fields,'
        ELSE IF (model.eq.iADM) THEN
          string='writing ADM state fields,'
        ELSE IF (model.eq.iRPM) THEN
          string='writing RPM state fields,'
        END IF
        Sstr=SCAN(CalledFrom,'/',BACK=.TRUE.)+1
        Send=LEN_TRIM(CalledFrom)
        WRITE (Tstring,'(f15.4)') stime*sec2day
        WRITE (stdout,10) TRIM(label), TRIM(string), t_code,            &
     &                    ng, TRIM(ADJUSTL(Tstring)), TRIM(S(ng)%name), &
     &                    i3d, OutRec, CalledFrom(Sstr:Send)
      END IF
!
!  Write out model time (s).
!
      IF (LwrtPER(ng)) THEN
        Tval(1)=REAL(OutRec,dp)*day2sec
      ELSE
        Tval(1)=stime
      END IF
      CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idtime)), tval,            &
     &                          (/OutRec/), (/1/),                      &
     &                          pioFile = pioFile,                      &
     &                          pioVar = S(ng)%pioVar(idtime)%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Write out free-surface.
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite2d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idFsur), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
#  ifdef WET_DRY
     &                     OCEAN(ng) % tl_zeta(:,:,i2d),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
#  else
     &                     OCEAN(ng) % tl_zeta(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
#  endif
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite2d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idFsur), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
#  ifdef WET_DRY
     &                     OCEAN(ng) % ad_zeta(:,:,i2d),                &
     &                     SetFillVal = .FALSE.,                        &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
#  else
     &                     OCEAN(ng) % ad_zeta(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
#  endif
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idFsur)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isFsur)
        IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        IF (model.eq.iTLM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % tl_zeta_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % ad_zeta_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.PIO_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)),                    &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Write out 2D U-momentum component (m/s).
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite2d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idUbar), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % tl_ubar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)

      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite2d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idUbar), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                            &
#  endif
     &                     OCEAN(ng) % ad_ubar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUbar)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isUbar)
        IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        IF (model.eq.iTLM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % tl_ubar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % ad_ubar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.PIO_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)),                    &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Write out 2D V-momentum component.
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idVbar)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite2d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idVbar), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % tl_vbar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite2d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idVbar), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     OCEAN(ng) % ad_vbar(:,:,i2d),                &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVbar)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isVbar)
        IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        IF (model.eq.iTLM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % tl_vbar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite2d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, Nbrec(ng), scale,          &
     &                           BOUNDARY(ng) % ad_vbar_obc(LBij:,:,:,  &
     &                                                      Lbout(ng)), &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.PIO_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)),                    &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif

#  ifdef ADJUST_WSTRESS
!
!  Write out surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idUsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idUsms), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     FORCES(ng) % tl_ustr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idUsms), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     FORCES(ng) % ad_ustr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUsms)), TRIM(label),         &
     &                      TRIM(S(ng)%name), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUsms)), Fmin, Fmax
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idVsms)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idVsms), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     FORCES(ng) % tl_vstr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)

      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idVsms), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale,     &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     FORCES(ng) % ad_vstr(:,:,:,Lfout(ng)),       &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVsms)), TRIM(label),         &
     &                      TRIM(S(ng)%name), Lfout(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVsms)), Fmin, Fmax
        END IF
      END IF
#  endif

#  ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idFsur)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idUvel), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     OCEAN(ng) % tl_u(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)

      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idUvel), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % umask,                            &
#   endif
     &                     OCEAN(ng) % ad_u(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idUvel)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isUvel)
        IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        IF (model.eq.iTLM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % tl_u_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % ad_u_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.PIO_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)),                    &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#   endif
!
!  Write out 3D V-momentum component (m/s).
!
      scale=1.0_dp
      IF (S(ng)%pioVar(idVvel)%dkind.eq.PIO_double) THEN
        ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      IF (model.eq.iTLM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idVvel), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     OCEAN(ng) % tl_v(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      ELSE IF (model.eq.iADM) THEN
        status=nf_fwrite3d(ng, model, pioFile,                          &
     &                     S(ng)%pioVar(idVvel), OutRec,                &
     &                     ioDesc,                                      &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % vmask,                            &
#   endif
     &                     OCEAN(ng) % ad_v(:,:,:,i3d),                 &
     &                     MinValue = Fmin,                             &
     &                     MaxValue = Fmax)
      END IF
      IF (status.ne.PIO_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(Vname(1,idVvel)), TRIM(label),         &
     &                      TRIM(S(ng)%name), OutRec
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF

#   ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_dp
        ifield=idSbry(isVvel)
        IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        IF (model.eq.iTLM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite3d_bry(ng, model, S(ng)%name, pioFile,        &
     &                           Vname(1,ifield),                       &
     &                           S(ng)%pioVar(ifield), OutRec,          &
     &                           ioDesc,                                &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng), scale,&
     &                           BOUNDARY(ng) % ad_v_obc(LBij:,:,:,:,   &
     &                                                   Lbout(ng)),    &
     &                           MinValue = Fmin,                       &
     &                           MaxValue = Fmax)
        END IF
        IF (status.ne.PIO_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,ifield)),                    &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#   endif
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        scale=1.0_dp
        IF (S(ng)%pioTrc(itrc)%dkind.eq.PIO_double) THEN
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        IF (model.eq.iTLM) THEN
          status=nf_fwrite3d(ng, model, pioFile,                        &
     &                       S(ng)%pioTrc(itrc), OutRec,                &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       OCEAN(ng) % tl_t(:,:,:,i3d,itrc),          &
     &                       MinValue = Fmin,                           &
     &                       MaxValue = Fmax)
        ELSE IF (model.eq.iADM) THEN
          status=nf_fwrite3d(ng, model, pioFile,                        &
     &                       S(ng)%pioTrc(itrc), OutRec,                &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#   endif
     &                       OCEAN(ng) % ad_t(:,:,:,i3d,itrc),          &
     &                       MinValue = Fmin,                           &
     &                       MaxValue = Fmax)
        END IF
        IF (status.ne.PIO_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(Vname(1,idTvar(itrc))),              &
     &                        TRIM(label), TRIM(S(ng)%name), OutRec
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(2,idTvar(itrc))),              &
     &                        Fmin, Fmax
          END IF
        END IF
      END DO

#   ifdef ADJUST_BOUNDARY
!
!  Write out tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_dp
          ifield=idSbry(isTvar(itrc))
          IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          IF (model.eq.iTLM) THEN
            status=nf_fwrite3d_bry (ng, model, S(ng)%name, pioFile,     &
     &                              Vname(1,ifield),                    &
     &                              S(ng)%pioVar(ifield),               &
     &                              OutRec, ioDesc,                     &
     &                              LBij, UBij, 1, N(ng), Nbrec(ng),    &
     &                              scale,                              &
     &                              BOUNDARY(ng) % tl_t_obc(LBij:,:,:,:,&
     &                                                Lbout(ng),itrc),  &
     &                              MinValue = Fmin,                    &
     &                              MaxValue = Fmax)
          ELSE IF (model.eq.iADM) THEN
            status=nf_fwrite3d_bry (ng, model, S(ng)%name, pioFile,     &
     &                              Vname(1,ifield),                    &
     &                              S(ng)%pioVar(ifield),               &
     &                              OutRec, ioDesc,                     &
     &                              LBij, UBij, 1, N(ng), Nbrec(ng),    &
     &                              scale,                              &
     &                              BOUNDARY(ng) % ad_t_obc(LBij:,:,:,:,&
     &                                                Lbout(ng),itrc),  &
     &                              MinValue = Fmin,                    &
     &                              MaxValue = Fmax)
          END IF
          IF (status.ne.PIO_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)),                  &
     &                          TRIM(label), TRIM(S(ng)%name), OutRec
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#   endif

#   ifdef ADJUST_STFLUX
!
!  Write out surface net tracers fluxes. Notice that fluxes have their
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          scale=1.0_dp                      ! kinematic flux units
          ifield=idTsur(itrc)
          IF (S(ng)%pioVar(ifield)%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          IF (model.eq.iTLM) THEN
            status=nf_fwrite3d(ng, iTLM, pioFile,                       &
     &                         S(ng)%pioVar(ifield),                    &
     &                         OutRec, ioDesc,                          &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale, &
#    ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#    endif
     &                         FORCES(ng)% tl_tflux(:,:,:,              &
     &                                              Lfout(ng),itrc),    &
     &                         MinValue = Fmin,                         &
     &                         MaxValue = Fmax)
          ELSE IF (model.eq.iADM) THEN
            status=nf_fwrite3d(ng, iTLM, pioFile,                       &
     &                         S(ng)%pioVar(ifield),                    &
     &                         OutRec, ioDesc,                          &
     &                         LBi, UBi, LBj, UBj, 1, Nfrec(ng), scale, &
#    ifdef MASKING
     &                         GRID(ng) % rmask,                        &
#    endif
     &                         FORCES(ng)% ad_tflux(:,:,:,              &
     &                                              Lfout(ng),itrc),    &
     &                         MinValue = Fmin,                         &
     &                         MaxValue = Fmax)
          END IF
          IF (status.ne.PIO_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,20) TRIM(Vname(1,ifield)),                  &
     &                          TRIM(label), TRIM(S(ng)%name), Lfout(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Synchronize tangent NetCDF file to disk to allow other processes
!  to access data immediately after it is written. Then, close NetCDF
!  file.
!-----------------------------------------------------------------------
!
      CALL pio_netcdf_sync (ng, iTLM, S(ng)%name, pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL pio_netcdf_close (ng, iTLM, pioFile, S(ng)%name, .FALSE.)

#  ifdef PROFILE
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, model, 81, __LINE__, MyFile)
#  endif
!
  10  FORMAT (1x,a,': WRT_STATE_PIO - ',a,t75,a,                        &
     &        /,23x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Index=',i1,', Rec=',i0,')',                            &
     &        /,23x,'Called from ''',a,'''')
  20  FORMAT (/,' WRT_STATE_PIO - error while writing variable: ',a,    &
     &        /,17x,'into ',a,' NetCDF file: ',a,                       &
     &        /,17x,'for time record: ',i0)
  30  FORMAT (20x,'- ',a,/,23x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE wrt_state_pio
# endif
#endif
      END MODULE wrt_state_mod
