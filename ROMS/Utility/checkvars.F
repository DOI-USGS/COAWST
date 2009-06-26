#include "cppdefs.h"
      SUBROUTINE checkvars (ng, model, ncname, string, Nrec, Nvar,      &
     &                      tvarnam, get_var, have_var)
!
!svn $Id: checkvars.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine checks if needed state variables are available in      !
!  requested NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ncname     NetCDF file name.                                     !
!     string     Identification string.                                !
!     Nvar       Size of logical switches arrays.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Nrec       Number of time records available.                     !
!     tvarnam    Name of time record variable.                         !
!     get_var    Logical switches (T/F), in terms of variable ID,      !
!                  indicating state variables needed by the model.     !
!     have_var   Logical switches (T/F), in terms of variable ID,      !
!                  indicating state variables available in NetCDF      !
!                  file.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Nvar
      integer, intent(inout) :: Nrec

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: string

      character (len=*), intent(inout) :: tvarnam

      logical, dimension(Nvar), intent(out) :: get_var
      logical, dimension(Nvar), intent(out) :: have_var
!
!  Local variable declarations.
!
      integer :: IDmod, i, itrc
!
!-----------------------------------------------------------------------
!  Determine state variables needed and check if they are available in
!  requested NetCDF file.
!-----------------------------------------------------------------------
!
!  Limit model identifier. The profiling is limited to iNLM, iTLM, iRPM,
!  and iADM.
!
      IF ((model.lt.1).or.(model.gt.4)) THEN
        IDmod=iNLM
      ELSE
        IDmod=model
      END IF
!
!  Inquire about the dimensions and check for consistency.
!
      CALL netcdf_check_dim (ng, IDmod, ncname)
      IF (exit_flag.ne.NoError) RETURN
      Nrec=rec_size
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, IDmod, ncname)        
      IF (exit_flag.ne.NoError) RETURN
!
!  Initialize logical switches.
!
      DO i=1,Nvar
        get_var(i)=.FALSE.
        have_var(i)=.FALSE.
      END DO
!
!  Determine state variables to read from input NetCDF file.
!
#ifdef ANA_INITIAL
      IF (nrrec(ng).ne.0) THEN
        get_var(idFsur)=.TRUE.
        get_var(idUbar)=.TRUE.
        get_var(idVbar)=.TRUE.
# ifdef SOLVE3D
        get_var(idUvel)=.TRUE.
        get_var(idVvel)=.TRUE.
        DO itrc=1,NAT
          get_var(idTvar(itrc))=.TRUE.
        END DO
# endif
      END IF
#else
      get_var(idFsur)=.TRUE.
      get_var(idUbar)=.TRUE.
      get_var(idVbar)=.TRUE.
# ifdef SOLVE3D
      get_var(idUvel)=.TRUE.
      get_var(idVvel)=.TRUE.
      DO itrc=1,NAT
        get_var(idTvar(itrc))=.TRUE.
      END DO
# endif
#endif
#ifdef SOLVE3D
# if defined BIOLOGY
#  if defined ANA_BIOLOGY
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NBT
          get_var(idTvar(idbio(itrc)))=.TRUE.
        END DO
      END IF
#  else
      DO itrc=1,NBT
        get_var(idTvar(idbio(itrc)))=.TRUE.
      END DO
#  endif
# endif
# if defined T_PASSIVE
#  if defined ANA_PASSIVE
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NPT
          get_var(idTvar(inert(itrc)))=.TRUE.
        END DO
      END IF
#  else
      DO itrc=1,NPT
        get_var(idTvar(inert(itrc)))=.TRUE.
      END DO
#  endif
# endif
# ifdef SEDIMENT
#  ifdef ANA_SEDIMENT
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NST
          get_var(idTvar(idsed(itrc)))=.TRUE.
          get_var(idfrac(itrc))=.TRUE.
          get_var(idBmas(itrc))=.TRUE.
        END DO
        DO itrc=1,MBEDP
          get_var(idSbed(itrc))=.TRUE.
        END DO
      END IF
#  else
      DO itrc=1,NST
        get_var(idTvar(idsed(itrc)))=.TRUE.
        get_var(idfrac(itrc))=.TRUE.
        get_var(idBmas(itrc))=.TRUE.
      END DO
      DO itrc=1,MBEDP
        get_var(idSbed(itrc))=.TRUE.
      END DO
#  endif
#  ifdef BEDLOAD
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NST
          get_var(idUbld(itrc))=.TRUE.
          get_var(idVbld(itrc))=.TRUE.
        END DO
      END IF
#  endif
#  ifdef SED_MORPH
      IF (nrrec(ng).ne.0) THEN
        get_var(idbath)=.TRUE.
      END IF
#  endif
# endif 
# if defined SEDIMENT || defined BBL_MODEL
      DO itrc=1,MBOTP
        get_var(idBott(itrc))=.TRUE.
      END DO
# endif
#endif
#ifdef PERFECT_RESTART
!
!  Determine perfect restart fields to read from input NetCDF file.
!
      IF ((model.eq.0).and.(nrrec(ng).ne.0)) THEN
# ifdef SOLVE3D
        get_var(idRu3d)=.TRUE.
        get_var(idRv3d)=.TRUE.
# endif
        get_var(idRzet)=.TRUE.
        get_var(idRu2d)=.TRUE.
        get_var(idRv2d)=.TRUE.
# if defined GLS_MIXING || defined MY25_MIXING
        get_var(idMtke)=.TRUE.
        get_var(idMtls)=.TRUE.
        get_var(idVmLS)=.TRUE.
        get_var(idVmKK)=.TRUE.
#  ifdef GLS_MIXING
        get_var(idVmKP)=.TRUE.
#  endif
# endif
      END IF
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Determine aditional 4Dvar state variables to adjust.
!
      IF ((model.eq.iADM).or.(model.eq.iTLM).or.(model.eq.iNLM).or.    &
     &    (model.eq.5).or.(model.eq.6)) THEN
# ifdef ADJUST_STFLUX
        DO itrc=1,NT(ng)
          get_var(idTsur(itrc))=.TRUE.
        END DO
# endif
# ifdef ADJUST_WSTRESS
        get_var(idUsms)=.TRUE.
        get_var(idVsms)=.TRUE.
# endif
      END IF
#endif
!
!  Scan variable list from input NetCDF and activate switches for
!  model state variables.
!
      DO i=1,n_var
        IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
          tvarnam=TRIM(var_name(i))
          have_var(idtime)=.TRUE.
#if defined SEDIMENT && defined SED_MORPH
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbath))) THEN
          have_var(idbath)=.TRUE.
#endif
#ifdef SOLVE3D
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idKver))) THEN
          have_var(idKver)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
          have_var(idUvel)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
          have_var(idRu3d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
          have_var(idVvel)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
          have_var(idRv3d)=.TRUE.
# if defined GLS_MIXING || defined MY25_MIXING || defined LMD_MIXING
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
          have_var(idVvis)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
          have_var(idTdif)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
          have_var(idSdif)=.TRUE.
# endif
# if defined GLS_MIXING || defined MY25_MIXING
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
          have_var(idMtke)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
          have_var(idMtls)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVmLS))) THEN
          have_var(idVmLS)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVmKK))) THEN
          have_var(idVmKK)=.TRUE.
#  ifdef GLS_MIXING
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVmKP))) THEN
          have_var(idVmKP)=.TRUE.
#  endif
# endif
#endif
#ifdef ADJUST_WSTRESS
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
          have_var(idUsms)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
          have_var(idVsms)=.TRUE.
#endif
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
          have_var(idFsur)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
          have_var(idRzet)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
          have_var(idUbar)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
          have_var(idRu2d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
          have_var(idVbar)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
          have_var(idRv2d)=.TRUE.
        ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idKhor))) THEN
          have_var(idKhor)=.TRUE.
        END IF
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
            have_var(idTvar(itrc))=.TRUE.
# ifdef ADJUST_STFLUX
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itrc)))) THEN
            have_var(idTsur(itrc))=.TRUE.
# endif
          END IF
        END DO
# ifdef SEDIMENT
        DO itrc=1,NST
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
            have_var(idfrac(itrc))=.TRUE.
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idBmas(itrc)))) THEN
            have_var(idBmas(itrc))=.TRUE.
#  ifdef BEDLOAD
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idUbld(itrc)))) THEN
            have_var(idUbld(itrc))=.TRUE.
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idVbld(itrc)))) THEN
            have_var(idVbld(itrc))=.TRUE.
#  endif
          END IF
        END DO
        DO itrc=1,MBEDP
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
            have_var(idSbed(itrc))=.TRUE.
          END IF
        END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
        DO itrc=1,MBOTP
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
            have_var(idBott(itrc))=.TRUE.
          END IF
        END DO
# endif
#endif
      END DO
!
!  Check if model state variables are available in input NetCDF file.
!
      IF (.not.have_var(idtime)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(tvarnam),            &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#if defined SEDIMENT && defined SED_MORPH
      IF (.not.have_var(idbath).and.get_var(idbath)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idbath)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#endif
      IF (.not.have_var(idFsur).and.get_var(idFsur)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idFsur)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRzet).and.get_var(idRzet)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRzet)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idUbar).and.get_var(idUbar)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idUbar)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRu2d).and.get_var(idRu2d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRu2d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVbar).and.get_var(idVbar)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVbar)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRv2d).and.get_var(idRv2d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRv2d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#ifdef ADJUST_WSTRESS
      IF ((model.eq.5).or.(model.eq.6)) THEN
        IF (.not.have_var(idUsms).and.get_var(idUsms)) THEN
          IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idUsms)),  &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVsms).and.get_var(idVsms)) THEN
          IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVsms)),  &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
#endif
#ifdef SOLVE3D
      IF (.not.have_var(idUvel).and.get_var(idUvel)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idUvel)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRu3d).and.get_var(idRu3d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRu3d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVvel).and.get_var(idVvel)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVvel)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idRv3d).and.get_var(idRv3d)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idRv3d)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      DO itrc=1,NT(ng)
        IF (.not.have_var(idTvar(itrc)).and.                            &
     &      get_var(idTvar(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idTvar(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# ifdef ADJUST_STFLUX
      IF ((model.eq.5).or.(model.eq.6)) THEN
        DO itrc=1,NT(ng)
          IF (.not.have_var(idTsur(itrc)).and.                          &
     &        get_var(idTsur(itrc))) THEN
            IF (Master) WRITE (stdout,10) string,                       &
     &                                    TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END DO
      END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
      IF (.not.have_var(idMtke).and.get_var(idMtke)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idMtke)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idMtls).and.get_var(idMtls)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idMtls)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVmLS).and.get_var(idVmLS)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVmLS)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
      IF (.not.have_var(idVmKK).and.get_var(idVmKK)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVmKK)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  ifdef GLS_MIXING
      IF (.not.have_var(idVmKP).and.get_var(idVmKP)) THEN
        IF (Master) WRITE (stdout,10) string, TRIM(Vname(1,idVmKP)),    &
     &                                TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
# endif
# ifdef SEDIMENT
      DO itrc=1,NST
        IF (.not.have_var(idfrac(itrc)).and.                            &
     &      get_var(idfrac(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idfrac(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idBmas(itrc)).and.                            &
     &      get_var(idBmas(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idBmas(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#  ifdef BEDLOAD
        IF (.not.have_var(idUbld(itrc)).and.                            &
     &      get_var(idUbld(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idUbld(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVbld(itrc)).and.                            &
     &      get_var(idVbld(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idVbld(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#  endif
      END DO
      DO itrc=1,MBEDP
        IF (.not.have_var(idSbed(itrc)).and.                            &
     &      get_var(idSbed(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idSbed(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
! Check only for the first four required properties, but we looked for
! MBOTP properties above.
!
      DO itrc=1,4
        IF (.not.have_var(idBott(itrc)).and.                            &
     &      get_var(idBott(itrc))) THEN
          IF (Master) WRITE (stdout,10) string,                         &
     &                                  TRIM(Vname(1,idBott(itrc))),    &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
# endif
#endif
!
  10  FORMAT (/,a,'CHECKVARS - unable to find model variable: ',a,      &
     &        /,18x,'in file: ',a)

      RETURN
      END SUBROUTINE checkvars
