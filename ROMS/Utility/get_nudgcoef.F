#include "cppdefs.h"
      MODULE get_nudgcoef_mod

#ifndef ANA_NUDGCOEF
!
!git $Id$
!svn $Id: get_nudgcoef.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module reads grid time scales for nudging to climatology file  !
!  using either the standard NetCDF library or the Parallel-IO (PIO)   !
!  library.                                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_clima
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_2d_mod
# ifdef SOLVE3D
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE mp_exchange_mod, ONLY : mp_exchange4d
#  endif
# endif
      USE nf_fread2d_mod,  ONLY : nf_fread2d
# ifdef SOLVE3D
      USE nf_fread3d_mod,  ONLY : nf_fread3d
# endif
      USE strings_mod,     ONLY : FoundError, find_string
!
      implicit none
!
      PUBLIC  :: get_nudgcoef
      PRIVATE :: get_nudgcoef_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: get_nudgcoef_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE get_nudgcoef (ng, tile, model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Read in nudging coefficients file according to IO type.
!-----------------------------------------------------------------------
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      SELECT CASE (NUD(ng)%IOtype)
        CASE (io_nf90)
          CALL get_nudgcoef_nf90 (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL get_nudgcoef_pio (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) NUD(ng)%IOtype
          exit_flag=2
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' GET_NUDGCOEF - Illegal input file type, io_type = ',i0, &
     &        /,16x,'Check KeyWord ''INP_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE get_nudgcoef
!
!***********************************************************************
      SUBROUTINE get_nudgcoef_nf90 (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
# ifdef SOLVE3D
      logical :: got_generic
      logical :: got_specific(NT(ng))
!
# endif
      integer :: gtype, i, ic, itrc
      integer :: nvatt, nvdim, status, vindex
      integer :: Vsize(4)
# ifdef CHECKSUM
      integer(i8b) :: Fhash
# endif
!
      real(dp) :: Fscl
      real(r8) :: Fmax, Fmin
!
      character (len=40 ) :: tunits
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_nudgcoef_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=NUD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (NUD(ng)%ncid.eq.-1) THEN
        CALL netcdf_open (ng, model, ncname, 0, NUD(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL netcdf_check_dim (ng, model, ncname,                         &
     &                       ncid = NUD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL netcdf_inq_var (ng, model, ncname,                           &
     &                     ncid = NUD(ng)%ncid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
!  Nudging coefficients for 2D momentum.
!
      IF (LnudgeM2CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM2nc),            &
     &                       NUD(ng)%Vid(idM2nc))) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idM2nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF

# ifdef SOLVE3D
!
!  Nudging coefficients for 3D momentum.
!
      IF (LnudgeM3CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM3nc),            &
     &                       NUD(ng)%Vid(idM3nc))) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idM3nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Tracers nudging coefficients.
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        got_generic=find_string(var_name,n_var,Vname(1,idgTnc),         &
     &                          NUD(ng)%Vid(idgTnc))
        DO itrc=1,NT(ng)
          IF (LnudgeTCLM(itrc,ng)) THEN
            vindex=idTnud(itrc)
            got_specific(itrc)=find_string(var_name,n_var,              &
     &                                     Vname(1,vindex),             &
     &                                     NUD(ng)%Vid(vindex))
            IF (.not.(got_generic.or.got_specific(itrc))) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vname(1,idgTnc)),      &
     &                                      TRIM(Vname(1,vindex)),      &
     &                                      TRIM(ncname)
              exit_flag=2
              RETURN
            END IF
          END IF
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Read in nudging coefficients.
!-----------------------------------------------------------------------
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  If appropriate, read nudging coefficients for 2D momentum  (inverse
!  time scales, 1/day).
!
      IF (Master) WRITE (stdout,'(1x)')
!
      IF (LnudgeM2CLM(ng)) THEN
        gtype=r2dvar
        vindex=idM2nc
        Fscl=1.0_dp/day2sec                    ! default units: 1/day
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = NUD(ng)%ncid,                       &
     &                       MyVarName = TRIM(Vname(1,vindex)),         &
     &                       VarID = NUD(ng)%Vid(vindex),               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_dp/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_dp
            END IF
          END IF
        END DO
!
        status=nf_fread2d(ng, model, ncname, NUD(ng)%ncid,              &
     &                    Vname(1,vindex), NUD(ng)%Vid(vindex),         &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
# ifdef CHECKSUM
     &                    CLIMA(ng) % M2nudgcof,                        &
     &                    checksum = Fhash)
# else
     &                    CLIMA(ng) % M2nudgcof)
# endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,vindex))//': '//             &
     &                        TRIM(Vname(1,vindex)),                    &
     &                        ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
            WRITE (stdout,50) Fhash
# endif
          END IF
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            CLIMA(ng) % M2nudgcof)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % M2nudgcof)
# endif
      END IF

# ifdef SOLVE3D
!
!  If appropriate, read nudging coefficients for 3D momentum (inverse
!  time scales, 1/day).
!
      IF (LnudgeM3CLM(ng)) THEN
        gtype=r3dvar
        vindex=idM3nc
        Fscl=1.0_dp/day2sec                    ! default units: 1/day
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = NUD(ng)%ncid,                       &
     &                       MyVarName = TRIM(Vname(1,vindex)),         &
     &                       VarID = NUD(ng)%Vid(vindex),               &
     &                       nVarDim = nvdim,                           &
     &                       nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_dp/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_dp
            END IF
          END IF
        END DO
!
        status=nf_fread3d(ng, model, ncname, NUD(ng)%ncid,              &
     &                    Vname(1,vindex), NUD(ng)%Vid(vindex),         &
     &                    0, gtype, Vsize,                              &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
#  ifdef CHECKSUM
     &                    CLIMA(ng) % M3nudgcof,                        &
     &                    checksum = Fhash)
#  else
     &                    CLIMA(ng) % M3nudgcof)
#  endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,vindex))//': '//             &
     &                        TRIM(Vname(1,vindex)),                    &
     &                        ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
            WRITE (stdout,50) Fhash
#  endif
          END IF
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            CLIMA(ng) % M3nudgcof)
        END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % M3nudgcof)
#  endif
      END IF
!
!  If appropriate, read nudging coefficients for tracer variables
!  (inverse time scales, 1/day). If the nudging coefficients for a
!  specific tracer are available, it will read that NetCDF variable.
!  If NOT AND the generic coefficients are available, it will process
!  those values instead.
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
          gtype=r3dvar
          IF (got_specific(itrc)) THEN
            vindex=idTnud(itrc)                ! specific variable
            Fscl=1.0_dp/day2sec                ! default units: 1/day
          ELSE IF (got_generic) THEN
            vindex=idgTnc                      ! generic variable
            Fscl=1.0_dp/day2sec                ! default units: 1/day
          END IF
!
          CALL netcdf_inq_var (ng, model, ncname,                       &
     &                         ncid = NUD(ng)%ncid,                     &
     &                         MyVarName = TRIM(Vname(1,vindex)),       &
     &                         VarID = NUD(ng)%Vid(vindex),             &
     &                         nVarDim = nvdim,                         &
     &                         nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'units') THEN
              tunits=TRIM(var_Achar(i))
              IF (tunits(1:3).eq.'day') THEN
                Fscl=1.0_dp/day2sec
              ELSE IF (tunits(1:6).eq.'second') THEN
                Fscl=1.0_dp
              END IF
            END IF
          END DO
!
          status=nf_fread3d(ng, model, ncname, NUD(ng)%ncid,            &
     &                      Vname(1,vindex), NUD(ng)%Vid(vindex),       &
     &                      0, gtype, Vsize,                            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
#  ifdef CHECKSUM
     &                      CLIMA(ng) % Tnudgcof(:,:,:,ic),             &
     &                      checksum = Fhash)
#  else
     &                      CLIMA(ng) % Tnudgcof(:,:,:,ic))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(2,vindex))//': '//           &
     &                          TRIM(Vname(1,vindex)),                  &
     &                          ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
              WRITE (stdout,50) Fhash
#  endif
            END IF
          END IF
!
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              CLIMA(ng) % Tnudgcof(:,:,:,ic))
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NTCLM(ng), &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % Tnudgcof)
      END IF
#  endif
# endif
!
  10  FORMAT (/,' GET_NUDGCOEF_NF90 - unable to open nudging NetCDF',   &
     &        ' file: ',a)
  20  FORMAT (/,' GET_NUDGCOEF_NF90 - unable to find nudging',          &
     &        ' variable: ',a,/,21x,'in NetCDF file: ',a)
  30  FORMAT (/,' GET_NUDGCOEF_NF90 - unable to find nudging',          &
     &        ' variable: ',a,                                          &
     &        /,21x,'or generic nudging variable: ',a,                  &
     &        /,21x,'in nudging NetCDF file: ',a)
  40  FORMAT(2x,' GET_NUDGCOEF_NF90 - ',a,/,21x,                        &
     &        '(Grid = ',i2.2,', File: ',a,')',/,21x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
# ifdef CHECKSUM
  50  FORMAT (19x,'(CheckSum = ',i0,')')
# endif
!
      RETURN
      END SUBROUTINE get_nudgcoef_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE get_nudgcoef_pio (ng, tile, model,                     &
     &                             LBi, UBi, LBj, UBj)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
!  Local variable declarations.
!
# ifdef SOLVE3D
      logical :: got_generic
      logical :: got_specific(NT(ng))
!
# endif
      integer :: gtype, i, ic, ifield, itrc
      integer :: nvatt, nvdim, status, vindex
      integer :: Vsize(4)
# ifdef CHECKSUM
      integer(i8b) :: Fhash
# endif
!
      real(dp) :: Fscl
      real(r8) :: Fmax, Fmin
!
      character (len=40 ) :: tunits
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", get_nudgcoef_nf90"
!
      TYPE (IO_Desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Inquire about the contents of grid NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=NUD(ng)%name
!
!  Open grid NetCDF file for reading.
!
      IF (NUD(ng)%pioFile%fh.eq.-1) THEN
        CALL pio_netcdf_open (ng, model, ncname, 0, NUD(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,10) TRIM(ncname)
          RETURN
        END IF
      END IF
!
!  Check grid file dimensions for consitency
!
      CALL pio_netcdf_check_dim (ng, model, ncname,                     &
     &                           pioFile = NUD(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
      CALL pio_netcdf_inq_var (ng, model, ncname,                       &
     &                         pioFile = NUD(ng)%pioFile)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Check if required variables are available.
!-----------------------------------------------------------------------
!
!  Nudging coefficients for 2D momentum.
!
      IF (LnudgeM2CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM2nc),            &
     &                       vindex)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idM2nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        NUD(ng)%pioVar(idM2nc)%vd=var_desc(vindex)
      END IF

# ifdef SOLVE3D
!
!  Nudging coefficients for 3D momentum.
!
      IF (LnudgeM3CLM(ng)) THEN
        IF (.not.find_string(var_name,n_var,Vname(1,idM3nc),            &
     &                       vindex)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idM3nc)),          &
     &                                  TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        NUD(ng)%pioVar(idM3nc)%vd=var_desc(vindex)
      END IF
!
!  Tracers nudging coefficients.
!
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        IF (find_string(var_name,n_var,Vname(1,idgTnc),vindex)) THEN
          got_generic=.TRUE.
          NUD(ng)%pioVar(idgTnc)%vd=var_desc(vindex)
        ELSE
          got_generic=.FALSE.
        END IF
!
        DO itrc=1,NT(ng)
          IF (LnudgeTCLM(itrc,ng)) THEN
            ifield=idTnud(itrc)
            IF (find_string(var_name,n_var,Vname(1,ifield),vindex)) THEN
              got_specific(itrc)=.TRUE.
              NUD(ng)%pioVar(ifield)%vd=var_desc(vindex)
            ELSE
              got_specific(itrc)=.FALSE.
            END IF
!
            IF (.not.(got_generic.or.got_specific(itrc))) THEN
              IF (Master) WRITE (stdout,30) TRIM(Vname(1,idgTnc)),      &
     &                                      TRIM(Vname(1,ifield)),      &
     &                                      TRIM(ncname)
              exit_flag=2
              RETURN
            END IF
          END IF
        END DO
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Read in nudging coefficients.
!-----------------------------------------------------------------------
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!  If appropriate, read nudging coefficients for 2D momentum  (inverse
!  time scales, 1/day).
!
      IF (Master) WRITE (stdout,'(1x)')
!
      IF (LnudgeM2CLM(ng)) THEN
        ifield=idM2nc
        Fscl=1.0_dp/day2sec                    ! default units: 1/day
        NUD(ng)%pioVar(ifield)%gtype=r2dvar
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = NUD(ng)%pioFile,             &
     &                           MyVarName = TRIM(Vname(1,ifield)),     &
     &                           pioVar = NUD(ng)%pioVar(ifield)%vd,    &
     &                           nVarDim = nvdim,                       &
     &                           nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_dp/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_dp
            END IF
          END IF
        END DO
!
        IF (KIND(CLIMA(ng)%M2nudgcof).eq.8) THEN
          NUD(ng)%pioVar(ifield)%dkind=PIO_double
          ioDesc => ioDesc_dp_r2dvar(ng)
        ELSE
          NUD(ng)%pioVar(ifield)%dkind=PIO_real
          ioDesc => ioDesc_sp_r2dvar(ng)
        END IF

        status=nf_fread2d(ng, model, ncname, NUD(ng)%pioFile,           &
     &                    Vname(1,ifield), NUD(ng)%pioVar(ifield),      &
     &                    0, ioDesc, Vsize,                             &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    Fscl, Fmin, Fmax,                             &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
# ifdef CHECKSUM
     &                    CLIMA(ng) % M2nudgcof,                        &
     &                    checksum = Fhash)
# else
     &                    CLIMA(ng) % M2nudgcof)
# endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield))//': '//             &
     &                        TRIM(Vname(1,ifield)),                    &
     &                        ng, TRIM(ncname), Fmin, Fmax
# ifdef CHECKSUM
            WRITE (stdout,50) Fhash
# endif
          END IF
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            CLIMA(ng) % M2nudgcof)
        END IF
# ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % M2nudgcof)
# endif
      END IF

# ifdef SOLVE3D
!
!  If appropriate, read nudging coefficients for 3D momentum (inverse
!  time scales, 1/day).
!
      IF (LnudgeM3CLM(ng)) THEN
        ifield=idM3nc
        Fscl=1.0_dp/day2sec                    ! default units: 1/day
        NUD(ng)%pioVar(ifield)%gtype=r3dvar
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = NUD(ng)%pioFile,             &
     &                           MyVarName = TRIM(Vname(1,ifield)),     &
     &                           pioVar = NUD(ng)%pioVar(ifield)%vd,    &
     &                           nVarDim = nvdim,                       &
     &                           nVarAtt = nvatt)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,nvatt
          IF (TRIM(var_Aname(i)).eq.'units') THEN
            tunits=TRIM(var_Achar(i))
            IF (tunits(1:3).eq.'day') THEN
              Fscl=1.0_dp/day2sec
            ELSE IF (tunits(1:6).eq.'second') THEN
              Fscl=1.0_dp
            END IF
          END IF
        END DO
!
        IF (KIND(CLIMA(ng)%M3nudgcof).eq.8) THEN
          NUD(ng)%pioVar(ifield)%dkind=PIO_double
          ioDesc => ioDesc_dp_r3dvar(ng)
        ELSE
          NUD(ng)%pioVar(ifield)%dkind=PIO_real
          ioDesc => ioDesc_sp_r3dvar(ng)
        END IF

        status=nf_fread3d(ng, model, ncname, NUD(ng)%pioFile,           &
     &                    Vname(1,ifield), NUD(ng)%pioVar(ifield),      &
     &                    0, ioDesc, Vsize,                             &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    Fscl, Fmin, Fmax,                             &
#  ifdef MASKING
     &                    GRID(ng) % rmask,                             &
#  endif
#  ifdef CHECKSUM
     &                    CLIMA(ng) % M3nudgcof,                        &
     &                    checksum = Fhash)
#  else
     &                    CLIMA(ng) % M3nudgcof)
#  endif
        IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
          exit_flag=2
          ioerror=status
          RETURN
        ELSE
          IF (Master) THEN
            WRITE (stdout,40) TRIM(Vname(2,ifield))//': '//             &
     &                        TRIM(Vname(1,ifield)),                    &
     &                        ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
            WRITE (stdout,50) Fhash
#  endif
          END IF
        END IF
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            CLIMA(ng) % M3nudgcof)
        END IF
#  ifdef DISTRIBUTE
        CALL mp_exchange3d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % M3nudgcof)
#  endif
      END IF
!
!  If appropriate, read nudging coefficients for tracer variables
!  (inverse time scales, 1/day). If the nudging coefficients for a
!  specific tracer are available, it will read that NetCDF variable.
!  If NOT AND the generic coefficients are available, it will process
!  those values instead.
!
      ic=0
      DO itrc=1,NT(ng)
        IF (LnudgeTCLM(itrc,ng)) THEN
          ic=ic+1
          IF (got_specific(itrc)) THEN
            ifield=idTnud(itrc)                ! specific variable
            Fscl=1.0_dp/day2sec                ! default units: 1/day
          ELSE IF (got_generic) THEN
            ifield=idgTnc                      ! generic variable
            Fscl=1.0_dp/day2sec                ! default units: 1/day
          END IF
          NUD(ng)%pioVar(ifield)%gtype=r3dvar
!
          CALL pio_netcdf_inq_var (ng, model, ncname,                   &
     &                             pioFile = NUD(ng)%pioFile,           &
     &                             MyVarName = TRIM(Vname(1,ifield)),   &
     &                             pioVar = NUD(ng)%pioVar(ifield)%vd,  &
     &                             nVarDim = nvdim,                     &
     &                             nVarAtt = nvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          DO i=1,nvatt
            IF (TRIM(var_Aname(i)).eq.'units') THEN
              tunits=TRIM(var_Achar(i))
              IF (tunits(1:3).eq.'day') THEN
                Fscl=1.0_dp/day2sec
              ELSE IF (tunits(1:6).eq.'second') THEN
                Fscl=1.0_dp
              END IF
            END IF
          END DO
!
          IF (KIND(CLIMA(ng)%Tnudgcof).eq.8) THEN
            NUD(ng)%pioVar(ifield)%dkind=PIO_double
            ioDesc => ioDesc_dp_r3dvar(ng)
          ELSE
            NUD(ng)%pioVar(ifield)%dkind=PIO_real
            ioDesc => ioDesc_sp_r3dvar(ng)
          END IF

          status=nf_fread3d(ng, model, ncname, NUD(ng)%pioFile,         &
     &                      Vname(1,ifield), NUD(ng)%pioVar(ifield),    &
     &                      0, ioDesc, Vsize,                           &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      Fscl, Fmin, Fmax,                           &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
#  ifdef CHECKSUM
     &                      CLIMA(ng) % Tnudgcof(:,:,:,ic),             &
     &                      checksum = Fhash)
#  else
     &                      CLIMA(ng) % Tnudgcof(:,:,:,ic))
#  endif
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            exit_flag=2
            ioerror=status
            RETURN
          ELSE
            IF (Master) THEN
              WRITE (stdout,40) TRIM(Vname(2,ifield))//': '//           &
     &                          TRIM(Vname(1,ifield)),                  &
     &                          ng, TRIM(ncname), Fmin, Fmax
#  ifdef CHECKSUM
              WRITE (stdout,50) Fhash
#  endif
            END IF
          END IF
!
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              CLIMA(ng) % Tnudgcof(:,:,:,ic))
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (ANY(LnudgeTCLM(:,ng))) THEN
        CALL mp_exchange4d (ng, tile, model, 1,                         &
     &                      LBi, UBi, LBj, UBj, 1, N(ng), 1, NTCLM(ng), &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      CLIMA(ng) % Tnudgcof)
      END IF
#  endif
# endif
!
  10  FORMAT (/,' GET_NUDGCOEF_PIO - unable to open nudging NetCDF',    &
     &        ' file: ',a)
  20  FORMAT (/,' GET_NUDGCOEF_PIO - unable to find nudging',           &
     &        ' variable: ',a,/,20x,'in NetCDF file: ',a)
  30  FORMAT (/,' GET_NUDGCOEF_PIO - unable to find nudging',           &
     &        ' variable: ',a,                                          &
     &        /,20x,'or generic nudging variable: ',a,                  &
     &        /,20x,'in nudging NetCDF file: ',a)
  40  FORMAT(2x,' GET_NUDGCOEF_PIO  - ',a,/,21x,                        &
     &        '(Grid = ',i2.2,', File: ',a,')',/,21x,                   &
     &        '(Min = ', 1p,e15.8,0p,' Max = ',1p,e15.8,0p,')')
# ifdef CHECKSUM
  50  FORMAT (19x,'(CheckSum = ',i0,')')
# endif
!
      RETURN
      END SUBROUTINE get_nudgcoef_pio
# endif
#endif
      END MODULE get_nudgcoef_mod
