#include "cppdefs.h"
      MODULE state_read_mod
!
!git $Id$
!svn $Id: state_read.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads the ROMS ocean state from requested NetCDF file. !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     tile         Domaing partition (integer)                         !
!     model        Calling model identifier (integer)                  !
!     IOtype       I/O NetCDF library to use (integer)                 !
!     LBi          I-dimension Lower bound (integer)                   !
!     UBi          I-dimension Upper bound (integer)                   !
!     LBj          J-dimension Lower bound (integer)                   !
!     UBj          J-dimension Upper bound (integer)                   !
!     LBij         IJ-dimension Lower bound (integer)                  !
!     UBij         IJ-dimension Upper bound (integer)                  !
!     Lout         State array index to process (integer)              !
!     rec          NetCDF file time record to read (integer)           !
!     nopen        Flag to open NetCDF file (integer)                  !
!                    nopen = 0   => file already opened for access     !
!                    nopen > 0   => open at entry and close at exit    !
!     ncid         NetCDF file ID, if already opened (integer)         !
#if defined PIO_LIB && defined DISTRIBUTE
!     pioFile      PIO file descriptor, TYPE(File_desc_t)              !
!                    pioFile%fh         file handler                   !
!                    pioFile%iosystem   IO system descriptor (struct)  !
#endif
!     ncname       NetCDF filename (string)                            !
#ifdef MASKING
!     rmask        State land/sea mask on RHO-points (real 2D array)   !
!     umask        State land/sea mask on U-points   (real 2D array)   !
!     vmask        State land/sea mask on V-points   (real 2D array)   !
#endif
!                                                                      !
!  On Output:                                                          !
!                                                                      !
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
!     s_t_obc      State tracer OBCs           (:,:,N,4,Nbrec,Lout,NT) !
!     s_u_obc      State 3D U-velocity OBC     (:,:,N,4,Nbrec,Lout)    !
!     s_v_obc      State 3D V-velocity OBC     (:,:,N,4,Nbrec,Lout)    !
# endif
!     s_ubar_obc   State 2D U-velocity OBC     (:,:,4,Nbrec,Lout)      !
!     s_vbar_obc   State 2D V-velocity OBC     (:,:,4,Nbrec,Lout)      !
!     s_zeta_obc   State free-surface OBC      (:,:,4,Nbrec,Lout)      !
#endif
#ifdef ADJUST_WSTRESS
!     s_ustr       State surface wind U-stress (:,:,Nfrec,Lout)        !
!     s_vstr       State surface wind V-stress (:,:,Nfrec,Lout)        !
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
!     s_tflux      State surface tracer flux   (:,:,Nfrec,Lout,NT)     !
# endif
!     s_t          State tracers               (:,:,N,Lout,NT)         !
!     s_u          State 3D U-velocity         (:,:,N,Lout)            !
!     s_v          State 3D V-velocity         (:,:,N,Lout)            !
#else
!     s_ubar       State 2D U-velocity         (:,:,Lout)              !
!     s_vbar       State 2D V-velocity         (:,:,Lout)              !
#endif
!     s_zeta       State free-surface          (:,:,Lout)              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE dateclock_mod,      ONLY : time_string
#ifdef DISTRIBUTE
      USE distribute_mod,     ONLY : mp_bcasti
#endif
#ifdef ADJUST_BOUNDARY
      USE nf_fread2d_bry_mod, ONLY : nf_fread2d_bry
# ifdef SOLVE3D
      USE nf_fread3d_bry_mod, ONLY : nf_fread3d_bry
# endif
#endif
      USE nf_fread2d_mod,     ONLY : nf_fread2d
#ifdef SOLVE3D
      USE nf_fread3d_mod,     ONLY : nf_fread3d
#endif
      USE strings_mod,        ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: state_read
      PRIVATE :: state_read_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: state_read_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE state_read (ng, tile, model, IOtype,                   &
     &                       LBi, UBi, LBj, UBj, LBij, UBij,            &
     &                       Lout, rec,                                 &
     &                       nopen, ncid,                               &
#if defined PIO_LIB && defined DISTRIBUTE
     &                       pioFile,                                   &
#endif
     &                       ncname,                                    &
#ifdef MASKING
     &                       rmask, umask, vmask,                       &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                       s_t_obc, s_u_obc, s_v_obc,                 &
# endif
     &                       s_ubar_obc, s_vbar_obc,                    &
     &                       s_zeta_obc,                                &
#endif
#ifdef ADJUST_WSTRESS
     &                       s_ustr, s_vstr,                            &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                       s_tflux,                                   &
# endif
     &                       s_t, s_u, s_v,                             &
#else
     &                       s_ubar, s_vbar,                            &
#endif
     &                       s_zeta)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, IOtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lout, rec, nopen

      integer, intent(inout) :: ncid

#if defined PIO_LIB && defined DISTRIBUTE
!
      TYPE (File_desc_t), intent(inout) :: pioFile
#endif
!
      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      logical :: Lreport = .FALSE.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in requested state vector
!-----------------------------------------------------------------------
!
      SELECT CASE (IOtype)
        CASE (io_nf90)
          CALL state_read_nf90 (ng, tile, model, Lreport,               &
     &                          LBi, UBi, LBj, UBj, LBij, UBij,         &
     &                          Lout, rec,                              &
     &                          nopen, ncid, ncname,                    &
#ifdef MASKING
     &                          rmask, umask, vmask,                    &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                          s_t_obc, s_u_obc, s_v_obc,              &
# endif
     &                          s_ubar_obc, s_vbar_obc,                 &
     &                          s_zeta_obc,                             &
#endif
#ifdef ADJUST_WSTRESS
     &                          s_ustr, s_vstr,                         &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                          s_tflux,                                &
# endif
     &                          s_t, s_u, s_v,                          &
#else
     &                          s_ubar, s_vbar,                         &
#endif
     &                          s_zeta)

#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL state_read_pio (ng, tile, model, Lreport,                &
     &                         LBi, UBi, LBj, UBj, LBij, UBij,          &
     &                         Lout, rec,                               &
     &                         nopen, pioFile, ncname,                  &
# ifdef MASKING
     &                         rmask, umask, vmask,                     &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                         s_t_obc, s_u_obc, s_v_obc,               &
#  endif
     &                         s_ubar_obc, s_vbar_obc,                  &
     &                         s_zeta_obc,                              &
# endif
# ifdef ADJUST_WSTRESS
     &                         s_ustr, s_vstr,                          &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                         s_tflux,                                 &
#  endif
     &                         s_t, s_u, s_v,                           &
# else
     &                         s_ubar, s_vbar,                          &
# endif
     &                         s_zeta)
#endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' STATE_READ - Illegal inpput file type, io_type = ',i0,  &
     &        /,14x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE state_read
!
!***********************************************************************
      SUBROUTINE state_read_nf90 (ng, tile, model, Lreport,             &
     &                            LBi, UBi, LBj, UBj, LBij, UBij,       &
     &                            Lout, rec,                            &
     &                            nopen, ncid, ncname,                  &
#ifdef MASKING
     &                            rmask, umask, vmask,                  &
#endif
#ifdef ADJUST_BOUNDARY
# ifdef SOLVE3D
     &                            s_t_obc, s_u_obc, s_v_obc,            &
# endif
     &                            s_ubar_obc, s_vbar_obc,               &
     &                            s_zeta_obc,                           &
#endif
#ifdef ADJUST_WSTRESS
     &                            s_ustr, s_vstr,                       &
#endif
#ifdef SOLVE3D
# ifdef ADJUST_STFLUX
     &                            s_tflux,                              &
# endif
     &                            s_t, s_u, s_v,                        &
#else
     &                            s_ubar, s_vbar,                       &
#endif
     &                            s_zeta)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lreport
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lout, rec, nopen

      integer, intent(inout) :: ncid
!
      character (len=*), intent(in) :: ncname
!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#  endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
# else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

#else

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#  endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
# endif
# ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#  endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
# else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
# endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
#endif
!
!  Local variable declarations.
!
      integer :: Sstr, Send
      integer :: i, j, k
      integer :: ifield, itrc
      integer :: gtype, my_ncid, status, varid

      integer, dimension(4) :: Vsize
!
      real(r8) :: Fmin, Fmax
      real(dp) :: stime, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", state_read_nf90"

#include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lout.
!-----------------------------------------------------------------------
#ifdef PROFILE
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, model, 80, __LINE__, MyFile)
#endif
!
!  Determine file and variables ids.
!
      IF ((nopen.gt.0).or.(ncid.eq.-1)) THEN
        CALL netcdf_open (ng, model, ncname, 0, my_ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ncid=my_ncid
      ELSE
        my_ncid=ncid
      END IF
!
      DO i=1,4
        Vsize(i)=0
      END DO

#ifdef SP4DVAR
!
!  Synchronize NetCDF data to disk. It allows processes open for reading
!  data to update the number of records written by other processes.
!  The writer flushes buffers to disk, and the reader makes sure that
!  it is reading from disk rather than from previously cached buffers.
!
      CALL netcdf_sync (ng, model, ncname, my_ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Read in time.
!
      CALL netcdf_get_time (ng, model, ncname, Vname(1,idtime),         &
     &                      Rclock%DateNumber, stime,                   &
     &                      my_ncid, (/rec/), (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Report information.
!
      IF (Master) THEN
        CALL time_string (stime, t_code)
        Sstr=SCAN(CalledFrom,'/',BACK=.TRUE.)+1
        Send=LEN_TRIM(CalledFrom)
        WRITE (Tstring,'(f15.4)') stime*sec2day
        WRITE (stdout,10) 'Reading state fields,', t_code,              &
     &                    ng, Tstring, TRIM(ncname), rec, Lout,         &
     &                    CalledFrom(Sstr:Send)
      END IF
!
!  Read in free-surface.
!
      gtype=r2dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread2d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idFsur), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#ifdef MASKING
     &                  rmask,                                          &
#endif
     &                  s_zeta(:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF

#ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        ifield=idSbry(isFsur)
        gtype=r2dvar
        scale=1.0_dp
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         my_ncid, varid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=nf_fread2d_bry (ng, model, ncname, my_ncid,              &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_zeta_obc(:,:,:,Lout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#endif

#ifndef SOLVE3D
!
!  Read in 2D U-momentum component.
!
      gtype=u2dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUbar),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread2d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idUbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask,                                          &
# endif
     &                  s_ubar(:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 2D V-momentum component.
!
      gtype=v2dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVbar),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread2d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idVbar), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask,                                          &
# endif
     &                  s_vbar(:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF
#endif

#ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        ifield=idSbry(isUbar)
        gtype=u2dvar
        scale=1.0_dp
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         my_ncid, varid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=nf_fread2d_bry (ng, model, ncname, my_ncid,              &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_ubar_obc(:,:,:,Lout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        ifield=idSbry(isVbar)
        gtype=v2dvar
        scale=1.0_dp
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         my_ncid, varid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=nf_fread2d_bry (ng, model, ncname, my_ncid,              &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, Nbrec(ng),                   &
     &                         scale, Fmin, Fmax,                       &
     &                         s_vbar_obc(:,:,:,Lout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#endif

#ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      gtype=u3dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUsms),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread3d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idUsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask,                                          &
# endif
     &                  s_ustr(:,:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUsms)), Fmin, Fmax
        END IF
      END IF
!
      gtype=v3dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVsms),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread3d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idVsms), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask,                                          &
# endif
     &                  s_vstr(:,:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVsms)), Fmin, Fmax
        END IF
      END IF
#endif

#ifdef SOLVE3D
!
!  Read in 3D U-momentum component.
!
      gtype=u3dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread3d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idUvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  umask,                                          &
# endif
     &                  s_u(:,:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        ifield=idSbry(isUvel)
        gtype=u3dvar
        scale=1.0_dp
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         my_ncid, varid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=nf_fread3d_bry (ng, model, ncname, my_ncid,              &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_u_obc(:,:,:,:,Lout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(1,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
# endif
!
!  Read in 3D V-momentum component.
!
      gtype=v3dvar
      scale=1.0_dp
      CALL netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),        &
     &                       my_ncid, varid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      status=nf_fread3d(ng, model, ncname, my_ncid,                     &
     &                  Vname(1,idVvel), varid,                         &
     &                  rec, gtype, Vsize,                              &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  vmask,                                          &
# endif
     &                  s_v(:,:,:,Lout))
      IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        ifield=idSbry(isVvel)
        gtype=v3dvar
        scale=1.0_dp
        CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),      &
     &                         my_ncid, varid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=nf_fread3d_bry (ng, model, ncname, my_ncid,              &
     &                         Vname(1,ifield), varid,                  &
     &                         rec, gtype,                              &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_v_obc(:,:,:,:,Lout))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
# endif
!
!  Read in tracers.
!
      gtype=r3dvar
      scale=1.0_dp
      DO itrc=1,NT(ng)
        CALL netcdf_inq_varid (ng, model, ncname,                       &
     &                         Vname(1,idTvar(itrc)), my_ncid, varid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=nf_fread3d(ng, model, ncname, my_ncid,                   &
     &                    Vname(1,idTvar(itrc)), varid,                 &
     &                    rec, gtype, Vsize,                            &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
# ifdef MASKING
     &                    rmask,                                        &
# endif
     &                    s_t(:,:,:,Lout,itrc))
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
          END IF
        END IF
      END DO

# ifdef ADJUST_BOUNDARY
!
!  Read in tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          ifield=idSbry(isTvar(itrc))
          gtype=r3dvar
          scale=1.0_dp
          CALL netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),    &
     &                           my_ncid, varid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          status=nf_fread3d_bry (ng, model, ncname, my_ncid,            &
     &                           Vname(1,ifield), varid,                &
     &                           rec, gtype,                            &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale, Fmin, Fmax,                     &
     &                           s_t_obc(:,:,:,:,Lout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master.and.Lreport) THEN
              WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
# endif

# ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      gtype=r3dvar
      scale=1.0_dp
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          CALL netcdf_inq_varid (ng, model, ncname,                     &
     &                           Vname(1,idTsur(itrc)), my_ncid, varid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

          status=nf_fread3d(ng, model, ncname, my_ncid,                 &
     &                      Vname(1,idTsur(itrc)), varid,               &
     &                      rec, gtype, Vsize,                          &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#  ifdef MASKING
     &                      rmask,                                      &
#  endif
     &                      s_tflux(:,:,:,Lout,itrc))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idTsur(itrc))), rec,       &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master.and.Lreport) THEN
              WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
# endif
#endif
!
!  Close current file.
!
      IF (nopen.gt.0) THEN
        CALL netcdf_close (ng, model, my_ncid, ncname, .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

#ifdef PROFILE
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, model, 80, __LINE__, MyFile)
#endif
!
  10  FORMAT ('STATE_READ_NF90 - ',a,t75,a,                             &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')',                          &
     &        /,19x,'Called from ''',a,'''')
  20  FORMAT (' STATE_READ_NF90 - unable to open NetCDF file: ',a)
  30  FORMAT (' STATE_READ_NF90 - error while reading variable: ',a,2x, &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)
  40  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE state_read_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE state_read_pio (ng, tile, model, Lreport,              &
     &                           LBi, UBi, LBj, UBj, LBij, UBij,        &
     &                           Lout, rec,                             &
     &                           nopen, pioFile, ncname,                &
# ifdef MASKING
     &                           rmask, umask, vmask,                   &
# endif
# ifdef ADJUST_BOUNDARY
#  ifdef SOLVE3D
     &                           s_t_obc, s_u_obc, s_v_obc,             &
#  endif
     &                           s_ubar_obc, s_vbar_obc,                &
     &                           s_zeta_obc,                            &
# endif
# ifdef ADJUST_WSTRESS
     &                           s_ustr, s_vstr,                        &
# endif
# ifdef SOLVE3D
#  ifdef ADJUST_STFLUX
     &                           s_tflux,                               &
#  endif
     &                           s_t, s_u, s_v,                         &
# else
     &                           s_ubar, s_vbar,                        &
# endif
     &                           s_zeta)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: Lreport
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBij, UBij
      integer, intent(in) :: Lout, rec, nopen
!
      TYPE (File_desc_t), intent(inout) :: pioFile
!
      character (len=*), intent(in) :: ncname
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:,:,:,:,:,:)
      real(r8), intent(inout) :: s_u_obc(LBij:,:,:,:,:)
      real(r8), intent(inout) :: s_v_obc(LBij:,:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_vbar_obc(LBij:,:,:,:)
      real(r8), intent(inout) :: s_zeta_obc(LBij:,:,:,:)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_vstr(LBi:,LBj:,:,:)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:,LBj:,:,:,:)
#   endif
      real(r8), intent(inout) :: s_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: s_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: s_v(LBi:,LBj:,:,:)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: s_vbar(LBi:,LBj:,:)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:,LBj:,:)

# else

#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ADJUST_BOUNDARY
#   ifdef SOLVE3D
      real(r8), intent(inout) :: s_t_obc(LBij:UBij,N(ng),4,             &
     &                                   Nbrec(ng),2,NT(ng))
      real(r8), intent(inout) :: s_u_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_v_obc(LBij:UBij,N(ng),4,Nbrec(ng),2)
#   endif
      real(r8), intent(inout) :: s_ubar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_vbar_obc(LBij:UBij,4,Nbrec(ng),2)
      real(r8), intent(inout) :: s_zeta_obc(LBij:UBij,4,Nbrec(ng),2)
#  endif
#  ifdef ADJUST_WSTRESS
      real(r8), intent(inout) :: s_ustr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
      real(r8), intent(inout) :: s_vstr(LBi:UBi,LBj:UBj,Nfrec(ng),2)
#  endif
#  ifdef SOLVE3D
#   ifdef ADJUST_STFLUX
      real(r8), intent(inout) :: s_tflux(LBi:UBi,LBj:UBj,               &
     &                                   Nfrec(ng),2,NT(ng))
#   endif
      real(r8), intent(inout) :: s_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: s_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: s_v(LBi:UBi,LBj:UBj,N(ng),2)
#  else
      real(r8), intent(inout) :: s_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: s_vbar(LBi:UBi,LBj:UBj,3)
#  endif
      real(r8), intent(inout) :: s_zeta(LBi:UBi,LBj:UBj,3)
# endif
!
!  Local variable declarations.
!
      integer :: Sstr, Send
      integer :: i, j, k
      integer :: ifield, itrc
      integer :: status

      integer, dimension(4) :: Vsize
!
      real(r8) :: Fmin, Fmax
      real(dp) :: stime, scale
!
      character (len=15) :: Tstring
      character (len=22) :: t_code

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", state_read_pio"
!
      TYPE (IO_desc_t), pointer :: my_ioDesc
      TYPE (File_desc_t) :: my_pioFile
      TYPE (My_VarDesc)  :: my_pioVar

# include "set_bounds.h"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in requested model state record. Load data into state array
!  index Lout.
!-----------------------------------------------------------------------
# ifdef PROFILE
!
!  Turn on time wall clock.
!
      CALL wclock_on (ng, model, 80, __LINE__, MyFile)
# endif
!
!  Determine file and variables ids.
!
      IF ((nopen.gt.0).or.(pioFile%fh.eq.-1)) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, my_pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        pioFile=my_pioFile
      ELSE
        my_pioFile=pioFile
      END IF
!
      DO i=1,4
        Vsize(i)=0
      END DO

#ifdef SP4DVAR
!
!  Synchronize NetCDF data to disk. It allows processes open for reading
!  data to update the number of records written by other processes.
!  The writer flushes buffers to disk, and the reader makes sure that
!  it is reading from disk rather than from previously cached buffers.
!
      CALL pio_netcdf_sync (ng, model, ncname, my_pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Read in time.
!
      CALL pio_netcdf_get_time (ng, model, ncname, Vname(1,idtime),     &
     &                          Rclock%DateNumber, stime,               &
     &                          my_pioFile, (/rec/), (/1/))
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Report information.
!
      IF (Master) THEN
        CALL time_string (stime, t_code)
        Sstr=SCAN(CalledFrom,'/',BACK=.TRUE.)+1
        Send=LEN_TRIM(CalledFrom)
        WRITE (Tstring,'(f15.4)') stime*sec2day
        WRITE (stdout,10) 'Reading state fields,', t_code,              &
     &                    ng, Tstring, TRIM(ncname), rec, Lout,         &
     &                    CalledFrom(Sstr:Send)
      END IF
!
!  Read in free-surface.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idFsur),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=r2dvar
      IF (KIND(s_zeta).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_r2dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_r2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idFsur), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
# ifdef MASKING
     &                  rmask,                                          &
# endif
     &                  s_zeta(:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsur)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF

# ifdef ADJUST_BOUNDARY
!
!  Read in free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        ifield=idSbry(isFsur)
        scale=1.0_dp
        CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),  &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        my_pioVar%gtype=r2dobc
        IF (KIND(s_zeta_obc).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_r2dobc(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_r2dobc(ng)
        END IF
!
        status=nf_fread2d_bry(ng, model, ncname, my_pioFile,            &
     &                        Vname(1,ifield), my_pioVar,               &
     &                        rec, my_ioDesc,                           &
     &                        LBij, UBij, Nbrec(ng),                    &
     &                        scale, Fmin, Fmax,                        &
     &                        s_zeta_obc(:,:,:,Lout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
# endif

# ifndef SOLVE3D
!
!  Read in 2D U-momentum component.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idUbar),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=u2dvar
      IF (KIND(s_ubar).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_u2dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_u2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idUbar), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ubar(:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 2D V-momentum component.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idVbar),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=v2dvar
      IF (KIND(s_vbar).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_v2dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_v2dvar(ng)
      END IF
!
      status=nf_fread2d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idVbar), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vbar(:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVbar)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF
# endif

# ifdef ADJUST_BOUNDARY
!
!  Read in 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        ifield=idSbry(isUbar)
        scale=1.0_dp
        CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),  &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        my_pioVar%gtype=u2dobc
        IF (KIND(s_ubar_obc).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_u2dobc(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_u2dobc(ng)
        END IF
!
        status=nf_fread2d_bry(ng, model, ncname, my_pioFile,            &
     &                        Vname(1,ifield), my_pioVar,               &
     &                        rec, my_ioDesc,                           &
     &                        LBij, UBij, Nbrec(ng),                    &
     &                        scale, Fmin, Fmax,                        &
     &                        s_ubar_obc(:,:,:,Lout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
!
!  Read in 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        ifield=idSbry(isVbar)
        scale=1.0_dp
        CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),  &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        my_pioVar%gtype=v2dobc
        IF (KIND(s_vbar_obc).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_v2dobc(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_v2dobc(ng)
        END IF
!
        status=nf_fread2d_bry(ng, model, ncname, my_pioFile,            &
     &                        Vname(1,ifield), my_pioVar,               &
     &                        rec, my_ioDesc,                           &
     &                        LBij, UBij, Nbrec(ng),                    &
     &                        scale, Fmin, Fmax,                        &
     &                        s_vbar_obc(:,:,:,Lout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
# endif

# ifdef ADJUST_WSTRESS
!
!  Read surface momentum stress.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idUsms),    &
     &                       my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=u2dvar
      IF (KIND(s_ustr).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_u2dfrc(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_u2dfrc(ng)
      END IF
!
      status=nf_fread3d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idUsms), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_ustr(:,:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUsms)), Fmin, Fmax
        END IF
      END IF
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idVsms),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=v2dvar
      IF (KIND(s_vstr).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_v2dfrc(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_v2dfrc(ng)
      END IF
!
      status=nf_fread3d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idVsms), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj, 1, Nfrec(ng),               &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_vstr(:,:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVsms)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVsms)), Fmin, Fmax
        END IF
      END IF
# endif

# ifdef SOLVE3D
!
!  Read in 3D U-momentum component.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idUvel),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=u3dvar
      IF (KIND(s_u).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_u3dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_u3dvar(ng)
      END IF
!
      status=nf_fread3d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idUvel), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  umask,                                          &
#  endif
     &                  s_u(:,:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        ifield=idSbry(isUvel)
        scale=1.0_dp
        CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),  &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        my_pioVar%gtype=u3dobc
        IF (KIND(s_u_obc).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_u3dobc(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_u3dobc(ng)
        END IF
!
        status=nf_fread3d_bry (ng, model, ncname, my_pioFile,           &
     &                         Vname(1,ifield), my_pioVar,              &
     &                         rec, my_ioDesc,                          &
     &                         LBij, UBij, 1, N(ng), Nbrec(ng),         &
     &                         scale, Fmin, Fmax,                       &
     &                         s_u_obc(:,:,:,:,Lout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(1,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Read in 3D V-momentum component.
!
      scale=1.0_dp
      CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,idVvel),    &
     &                           my_pioFile, my_pioVar%vd)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      my_pioVar%gtype=v3dvar
      IF (KIND(s_v).eq.8) THEN
        my_pioVar%dkind=PIO_double
        my_ioDesc => ioDesc_dp_v3dvar(ng)
      ELSE
        my_pioVar%dkind=PIO_real
        my_ioDesc => ioDesc_sp_v3dvar(ng)
      END IF
!
      status=nf_fread3d(ng, model, ncname, my_pioFile,                  &
     &                  Vname(1,idVvel), my_pioVar,                     &
     &                  rec, my_ioDesc, Vsize,                          &
     &                  LBi, UBi, LBj, UBj, 1, N(ng),                   &
     &                  scale, Fmin, Fmax,                              &
#  ifdef MASKING
     &                  vmask,                                          &
#  endif
     &                  s_v(:,:,:,Lout))
      IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
        IF (Master) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVvel)), rec, TRIM(ncname)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      ELSE
        IF (Master.and.Lreport) THEN
          WRITE (stdout,40) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF

#  ifdef ADJUST_BOUNDARY
!
!  Read in 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        ifield=idSbry(isVvel)
        scale=1.0_dp
        CALL pio_netcdf_inq_varid (ng, model, ncname, Vname(1,ifield),  &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        my_pioVar%gtype=v3dobc
        IF (KIND(s_v_obc).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_v3dobc(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_v3dobc(ng)
        END IF
!
        status=nf_fread3d_bry(ng, model, ncname, my_pioFile,            &
     &                        Vname(1,ifield), my_pioVar,               &
     &                        rec, my_ioDesc,                           &
     &                        LBij, UBij, 1, N(ng), Nbrec(ng),          &
     &                        scale, Fmin, Fmax,                        &
     &                        s_v_obc(:,:,:,:,Lout))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END IF
#  endif
!
!  Read in tracers.
!
      scale=1.0_dp
      DO itrc=1,NT(ng)
        CALL pio_netcdf_inq_varid (ng, model, ncname,                   &
     &                             Vname(1,idTvar(itrc)),               &
     &                             my_pioFile, my_pioVar%vd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        my_pioVar%gtype=r3dvar
        IF (KIND(s_t).eq.8) THEN
          my_pioVar%dkind=PIO_double
          my_ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          my_pioVar%dkind=PIO_real
          my_ioDesc => ioDesc_sp_r3dvar(ng)
        END IF
!
        status=nf_fread3d(ng, model, ncname, my_pioFile,                &
     &                    Vname(1,idTvar(itrc)), my_pioVar,             &
     &                    rec, my_ioDesc, Vsize,                        &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    scale, Fmin, Fmax,                            &
#  ifdef MASKING
     &                    rmask,                                        &
#  endif
     &                    s_t(:,:,:,Lout,itrc))
        IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
          IF (Master) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), rec,         &
     &                        TRIM(ncname)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        ELSE
          IF (Master.and.Lreport) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
          END IF
        END IF
      END DO

#  ifdef ADJUST_BOUNDARY
!
!  Read in tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          ifield=idSbry(isTvar(itrc))
          scale=1.0_dp
          CALL pio_netcdf_inq_varid (ng, model, ncname,                 &
     &                               Vname(1,ifield),                   &
     &                               my_pioFile, my_pioVar%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          my_pioVar%gtype=r3dobc
          IF (KIND(s_t_obc).eq.8) THEN
            my_pioVar%dkind=PIO_double
            my_ioDesc => ioDesc_dp_r3dobc(ng)
          ELSE
            my_pioVar%dkind=PIO_real
            my_ioDesc => ioDesc_sp_r3dobc(ng)
          END IF
!
          status=nf_fread3d_bry (ng, model, ncname, my_pioFile,         &
     &                           Vname(1,ifield), my_pioVar,            &
     &                           rec, my_ioDesc,                        &
     &                           LBij, UBij, 1, N(ng), Nbrec(ng),       &
     &                           scale, Fmin, Fmax,                     &
     &                           s_t_obc(:,:,:,:,Lout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,ifield)), rec, TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master.and.Lreport) THEN
              WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#  endif

#  ifdef ADJUST_STFLUX
!
!  Read in surface tracers flux.
!
      scale=1.0_dp
      DO itrc=1,NT(ng)
        IF (Lstflux(itrc,ng)) THEN
          CALL pio_netcdf_inq_varid (ng, model, ncname,                 &
     &                               Vname(1,idTsur(itrc)),             &
     &                               my_pioFile, my_pioVar%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          my_pioVar%gtype=r2dvar
          IF (KIND(s_tflux).eq.8) THEN
            my_pioVar%dkind=PIO_double
            my_ioDesc => ioDesc_dp_r2dfrc(ng)
          ELSE
            my_pioVar%dkind=PIO_real
            my_ioDesc => ioDesc_sp_r2dfrc(ng)
          END IF
!
          status=nf_fread3d(ng, model, ncname, my_pioFile,              &
     &                      Vname(1,idTsur(itrc)), my_pioVar,           &
     &                      rec, my_ioDesc, Vsize,                      &
     &                      LBi, UBi, LBj, UBj, 1, Nfrec(ng),           &
     &                      scale, Fmin, Fmax,                          &
#   ifdef MASKING
     &                      rmask,                                      &
#   endif
     &                      s_tflux(:,:,:,Lout,itrc))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,30) TRIM(Vname(1,idTsur(itrc))), rec,       &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          ELSE
            IF (Master.and.Lreport) THEN
              WRITE (stdout,40) TRIM(Vname(2,ifield)), Fmin, Fmax
            END IF
          END IF
        END IF
      END DO
#  endif
# endif
!
!  Close current file.
!
      IF (nopen.gt.0) THEN
        CALL pio_netcdf_close (ng, model, my_pioFile, ncname, .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF

# ifdef PROFILE
!
!  Turn off time wall clock.
!
      CALL wclock_off (ng, model, 80, __LINE__, MyFile)
# endif
!
  10  FORMAT ('STATE_READ_PIO  - ',a,t75,a,                             &
     &        /,19x,'(Grid ',i2.2,', t = ',a,', File: ',a,              &
     &        ', Rec=',i4.4,', Index=',i1,')',                          &
     &        /,19x,'Called from ''',a,'''')
  20  FORMAT (' STATE_READ_PIO  - unable to open NetCDF file: ',a)
  30  FORMAT (' STATE_READ_PIO  - error while reading variable: ',a,2x, &
     &        'at time record = ',i3,/,14x,'in NetCDF file: ',a)
  40  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
!
      RETURN
      END SUBROUTINE state_read_pio
#endif
      END MODULE state_read_mod
