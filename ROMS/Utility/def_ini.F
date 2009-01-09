#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE def_ini (ng)
!
!svn $Id: def_ini.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates a initialization NetCDF file used for          !
!  variational data assimilation.  It defines its dimensions,          !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 24

      logical :: Ldefine = .FALSE.
      logical :: got_var(NV)

      integer :: i, j, nvd4
      integer :: frecdim, status
# ifdef SOLVE3D
      integer :: itrc
# endif

      integer :: DimIDs(29), t4dfrc(4), u4dfrc(4), v4dfrc(4)
      integer :: Vsize(4)

      integer :: def_var

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!=======================================================================
!  Open existing nonlinear model initial conditions and define new
!  new variables.
!=======================================================================
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=INIname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, iNLM, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, iNLM, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!  Check if surface forcing variables have been already defined.
!
      DO i=1,NV
        got_var(i)=.FALSE.
      END DO
!
      DO i=1,n_var
#  ifdef ADJUST_WSTRESS
        IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
          got_var(idUsms)=.TRUE.
          iniVid(idUsms,ng)=var_id(i)
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
          got_var(idVsms)=.TRUE.
          iniVid(idVsms,ng)=var_id(i)
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
        DO itrc=1,NT(ng)
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
            got_var(idTsur(itrc))=.TRUE.
            iniVid(idTsur(itrc),ng)=var_id(i)
          END IF
        END DO
#  endif
      END DO
        
#  ifdef ADJUST_WSTRESS
      IF (.not.got_var(idUsms)) Ldefine=.TRUE.
      IF (.not.got_var(idVsms)) Ldefine=.TRUE.
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
      DO itrc=1,NT(ng)
        IF (.not.got_var(idTsur(itrc))) Ldefine=.TRUE.
      END DO
#  endif
!
!  Put existing file into define mode so new variables can be added.
!
      IF (Ldefine) THEN
        IF (OutThread) THEN
          status=nf90_open(TRIM(ncname),nf90_write,ncINIid(ng))
          IF (status.eq.nf90_noerr) THEN
            status=nf90_redef(ncINIid(ng))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,20) TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#  endif
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (Ldefine.and.OutThread) THEN
        DO i=1,n_dim
          SELECT CASE (TRIM(ADJUSTL(dim_name(i))))
            CASE ('xi_rho')
              DimIDs( 1)=dim_id(i)
            CASE ('xi_u')
              DimIDs( 2)=dim_id(i)
            CASE ('xi_v')
              DimIDs( 3)=dim_id(i)
            CASE ('eta_rho')
              DimIDs( 5)=dim_id(i)
            CASE ('eta_u')
              DimIDs( 6)=dim_id(i)
            CASE ('eta_v')
              DimIDs( 7)=dim_id(i)
#  ifdef SOLVE3D
            CASE ('s_rho')
              DimIDs( 9)=dim_id(i)
            CASE ('s_w')
              DimIDs(10)=dim_id(i)
#  endif
          END SELECT
        END DO
        DimIDs(12)=rec_id
        status=nf90_def_dim(ncINIid(ng),'Nadjust',Nfrec(ng),frecdim)
!
!  Set number of dimensions for output variables.
!
        nvd4=4
!
!  Define dimension vectors for staggered tracer type variables.
!
        t4dfrc(1)=DimIDs( 1)
        t4dfrc(2)=DimIDs( 5)
        t4dfrc(3)=frecdim
        t4dfrc(4)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
        u4dfrc(1)=DimIDs( 2)
        u4dfrc(2)=DimIDs( 6)
        u4dfrc(3)=frecdim
        u4dfrc(4)=DimIDs(12)
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
        v4dfrc(1)=DimIDs( 3)
        v4dfrc(2)=DimIDs( 7)
        v4dfrc(3)=frecdim
        v4dfrc(4)=DimIDs(12)
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,80
            Vinfo(i)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define surface forcing variables.
!-----------------------------------------------------------------------

#  ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
        IF (.not.got_var(idUsms)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncINIid(ng),iniVid(idUsms,ng),NF_FOUT,         &
     &                   nvd4,u4dfrc,Aval,Vinfo,ncname)
        END IF
!
!  Define surface V-momentum stress.
!
        IF (.not.got_var(idVsms)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncINIid(ng),iniVid(idVsms,ng),NF_FOUT,         &
     &                   nvd4,v4dfrc,Aval,Vinfo,ncname)
        END IF
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
!
!  define surface tracer flux. Notice that additional tracer fluxes are
!  written at their own fixed time-dimension (of size Nfrec) to allow
!  4DVAR adjustments at other times in addition to initial time.
!
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTsur(itrc))) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTsur(itrc)))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ncINIid(ng),iniVid(idTsur(itrc),ng),         &
     &                     NF_FOUT,nvd4,t4dfrc,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncINIid(ng))
      END IF
# endif
!
!=======================================================================
!  Open an existing initialization file, check its contents, and
!  prepare for appending data.
!=======================================================================
!
      IF (.not.LdefINI(ng)) THEN
        ncname=INIname(ng)

# if !(defined ADJUST_STFLUX || defined ADJUST_WSTRESS)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open initialization file for read/write.
!
        IF (OutThread) THEN
          status=nf90_open(TRIM(ncname),nf90_write,ncINIid(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#  ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag, 1)
#  endif
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
# endif
!
!  Scan variable list from input NetCDF and activate switches for
!  initialization variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            iniVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            iniVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            iniVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            iniVid(idVbar,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            iniVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            iniVid(idVvel,ng)=var_id(i)
#  if defined BVF_MIXING || defined LMD_MIXING || \
      defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            iniVid(idVvis,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            iniVid(idTdif,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            iniVid(idSdif,ng)=var_id(i)
#  endif
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              iniTid(itrc,ng)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,30) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        tINIindx(ng)=rec_size
      END IF
!
  10  FORMAT (/,' DEF_INI - unable to put in define mode initial',      &
     &        ' NetCDF file: ',a)
  20  FORMAT (/,' DEF_INI - unable to open initial NetCDF file: ',a)
  30  FORMAT (/,' DEF_INI - unable to find variable: ',a,2x,            &
     &        ' in initial NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_ini
#else
      SUBROUTINE def_ini
      RETURN
      END SUBROUTINE def_ini
#endif
