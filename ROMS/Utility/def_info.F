#include "cppdefs.h"
      SUBROUTINE def_info (ng, model, ncid, ncname, DimIDs)
!
!svn $Id: def_info.F 807 2016-07-09 02:03:55Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer).                           !
!     model    Calling model identifier (integer).                     !
!     ncid     NetCDF file ID (integer).                               !
!     ncname   NetCDF file name (character).                           !
!     DimIDs   NetCDF dimensions IDs (integer vector):                 !
!                DimIDs( 1) => XI-dimension at RHO-points.             !
!                DimIDs( 2) => XI-dimension at U-points.               !
!                DimIDs( 3) => XI-dimension at V-points.               !
!                DimIDs( 4) => XI-dimension at PSI-points.             !
!                DimIDs( 5) => ETA-dimension at RHO-points.            !
!                DimIDs( 6) => ETA-dimension at U-points.              !
!                DimIDs( 7) => ETA-dimension at V-points.              !
!                DimIDs( 8) => ETA-dimension at PSI-points.            !
!                DimIDs( 9) => S-dimension at RHO-points.              !
!                DimIDs(10) => S-dimension at W-points.                !
!                DimIDs(11) => Number of tracers dimension.            !
!                DimIDs(12) => Unlimited time record dimension.        !
!                DimIDs(13) => Number of stations dimension.           !
!                DimIDs(14) => Boundary dimension.                     !
!                DimIDs(15) => Number of floats dimension.             !
!                DimIDs(16) => Number sediment bed layers dimension.   !
!                DimIDs(17) => Dimension 2D water RHO-points.          !
!                DimIDs(18) => Dimension 2D water U-points.            !
!                DimIDs(19) => Dimension 2D water V-points.            !
!                DimIDs(20) => Dimension 3D water RHO-points.          !
!                DimIDs(21) => Dimension 3D water U-points.            !
!                DimIDs(23) => Dimension 3D water W-points.            !
!                DimIDs(24) => Dimension sediment bed water points.    !
!                DimIDs(25) => Number of EcoSim phytoplankton groups.  !
!                DimIDs(26) => Number of EcoSim bateria groups.        !
!                DimIDs(27) => Number of EcoSim DOM groups.            !
!                DimIDs(28) => Number of EcoSim fecal groups.          !
!                DimIDs(29) => Number of state variables.              !
!                DimIDs(30) => Number of 3D variables time levels (2). !
!                DimIDs(31) => Number of 2D variables time levels (3). !
!                DimIDs(32) => Number of sediment tracers.             !
!                DimIDs(33) => Number of wave directions.              !
!                DimIDs(35) => Number of vegetation properties.        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings
!
      USE def_var_mod, ONLY : def_var
#if !defined PARALLEL_OUT  && defined DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
#endif
      USE strings_mod, ONLY : join_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid
#ifdef INWAVE_MODEL
      integer, intent(in) :: DimIDs(33)
#else
      integer, intent(in) :: DimIDs(32)
#endif
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25

      integer :: brydim, i, ie, is, j, lstr, varid
      integer :: srdim, stadim, status, swdim, trcdim, usrdim
#ifdef SEDIMENT
      integer :: seddim
#endif
#ifdef FOUR_DVAR
      integer :: statedim
#endif
#if defined BIOLOGY && defined ECOSIM
      integer :: bacdim, domdim, fecdim, phydim

      integer :: biodim(2)
#endif
# if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
      integer :: domdim
#endif
#if !defined PARALLEL_OUT  && defined DISTRIBUTE
      integer :: ibuffer(2)
#endif
      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=11 ) :: bryatt
      character (len=11 ) :: frcatt
      character (len=11 ) :: clmatt
      character (len=50 ) :: tiling
      character (len=80 ) :: type
      character (len=120) :: Vinfo(Natt)
#ifdef FOUR_DVAR
      character (len=512) :: state_vector
#endif
#ifdef BIOLOGY
      character (len=512) :: bio_file
#endif
      character (len=1024) :: string
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
#ifdef SEDIMENT
      seddim=DimIDs(32)
#endif
      stadim=DimIDs(13)
      brydim=DimIDs(14)
#ifdef FOUR_DVAR
      statedim=DimIDs(29)
#endif
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
#if defined ECOSIM && defined SOLVE3D
      phydim=DimIDs(25)
      bacdim=DimIDs(26)
      domdim=DimIDs(27)
      fecdim=DimIDs(28)
      biodim(1)=phydim
      biodim(2)=fecdim
#endif
# if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
      domdim=DimIDs(27)
#endif
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(ncid.ne.GST(ng)%ncid)) THEN
        status=def_dim(ng, model, ncid, ncname, 'Nuser', 25, usrdim)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,LEN(Vinfo(1))
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,6
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  Define history global attribute.
!
        IF (LEN_TRIM(date_str).gt.0) THEN
          WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',       &
     &                                      TRIM( version),             &
     &                                      TRIM(date_str)
        ELSE
          WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version',              &
     &                               TRIM(version)
        END IF
!
!  Set tile decomposition global attribute.
!
        WRITE (tiling,10) NtileI(ng), NtileJ(ng)
!
!  Define file name global attribute.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'file',                &
     &                        TRIM(ncname))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

#ifndef DEBUGGING
!
!  Define NetCDF format type.
!
        IF (exit_flag.eq.NoError) THEN
# ifdef HDF5
          status=nf90_put_att(ncid, nf90_global, 'format',              &
     &                        'netCDF-4/HDF5 file')
# else
          status=nf90_put_att(ncid, nf90_global, 'format',              &
     &                        'netCDF-3 64bit offset file')
# endif
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'format', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

#endif
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
        type='CF-1.4'
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'Conventions',         &
     &                        TRIM(type))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'Conventions', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Define file type global attribute.
!
        IF (ncid.eq.ADM(ng)%ncid) THEN
          type='ROMS/TOMS adjoint history file'
        ELSE IF (ncid.eq.AVG(ng)%ncid) THEN
#if defined AD_AVERAGES && defined ADJOINT
          type='ROMS/TOMS adjoint model averages file'
#elif defined RP_AVERAGES && defined TL_IOMS
          type='ROMS/TOMS representer model averages file'
#elif defined TL_AVERAGES && defined TANGENT
          type='ROMS/TOMS tangent linear model averages file'
#else
          type='ROMS/TOMS nonlinear model averages file'
#endif
        ELSE IF (ncid.eq.DIA(ng)%ncid) THEN
          type='ROMS/TOMS diagnostics file'
        ELSE IF (ncid.eq.FLT(ng)%ncid) THEN
          type='ROMS/TOMS floats file'
        ELSE IF (ncid.eq.ERR(ng)%ncid) THEN
          type='ROMS/TOMS posterior analysis error covariance matrix'
        ELSE IF (ncid.eq.GST(ng)%ncid) THEN
          type='ROMS/TOMS GST check pointing restart file'
        ELSE IF (ncid.eq.HSS(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Hessian eigenvectors file'
        ELSE IF (ncid.eq.HIS(ng)%ncid) THEN
          type='ROMS/TOMS history file'
        ELSE IF (ncid.eq.ITL(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear model initial file'
        ELSE IF (ncid.eq.LCZ(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Lanczos vectors file'
        ELSE IF (ncid.eq.LZE(ng)%ncid) THEN
          type='ROMS/TOMS 4DVAR Evolved Lanczos vectors file'
        ELSE IF (ncid.eq.NRM(1,ng)%ncid) THEN
          type='ROMS/TOMS initial conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(2,ng)%ncid) THEN
          type='ROMS/TOMS model error covariance norm file'
        ELSE IF (ncid.eq.NRM(3,ng)%ncid) THEN
         type='ROMS/TOMS boundary conditions error covariance norm file'
        ELSE IF (ncid.eq.NRM(4,ng)%ncid) THEN
          type='ROMS/TOMS surface forcing error covariance norm file'
        ELSE IF (ncid.eq.RST(ng)%ncid) THEN
          type='ROMS/TOMS restart file'
        ELSE IF (ncid.eq.STA(ng)%ncid) THEN
          type='ROMS/TOMS station file'
        ELSE IF (ncid.eq.TLF(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear impulse forcing file'
        ELSE IF (ncid.eq.TLM(ng)%ncid) THEN
          type='ROMS/TOMS tangent linear history file'
        END IF
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'type',                &
     &                        TRIM(type))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'type', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

#ifdef FOUR_DVAR
!
!  Set state vector variables.
!
        is=1
        state_vector=' '
        DO i=1,NstateVar(ng)
          lstr=LEN_TRIM(Vname(1,idSvar(i)))
          ie=is+lstr
          state_vector(is:ie)=TRIM(Vname(1,idSvar(i)))//', '
          is=ie+2
        END DO
#endif
!
!  Define other global attributes to NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'title',               &
     &                        TRIM(title))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'title', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'var_info',            &
     &                        TRIM(varname))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'var_info', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

#ifdef FOUR_DVAR
        IF (exit_flag.eq.NoError) THEN
          lstr=LEN_TRIM(state_vector)-1
          status=nf90_put_att(ncid, nf90_global, 'state_vector',        &
     &                        state_vector(1:lstr))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'state_vector', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef PROPAGATOR
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'gst_file',            &
     &                        TRIM(GST(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'gst_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'rst_file',            &
     &                        TRIM(RST(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'rst_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).gt.0) THEN
              status=nf90_put_att(ncid, nf90_global, 'his_base',        &
     &                            TRIM(HIS(ng)%base))
            ELSE
              status=nf90_put_att(ncid, nf90_global, 'his_file',        &
     &                            TRIM(HIS(ng)%name))
            END IF
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) 'his_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        IF (exit_flag.eq.NoError) THEN
          IF (ndefAVG(ng).gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'avg_base',          &
     &                          TRIM(AVG(ng)%base))
          ELSE
            status=nf90_put_att(ncid, nf90_global, 'avg_file',          &
     &                          TRIM(AVG(ng)%name))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'avg_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif
#ifdef AVERAGES2
        IF (exit_flag.eq.NoError) THEN
          IF (ndefAVG2(ng).gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'avg2_base',         &
     &                          TRIM(AVG2(ng)%base))
          ELSE
            status=nf90_put_att(ncid, nf90_global, 'avg2_file',         &
     &                          TRIM(AVG2(ng)%name))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'avg2_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif
#ifdef HISTORY2
        IF (exit_flag.eq.NoError) THEN
          IF (ndefHIS2(ng).gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'his2_base',         &
     &                          TRIM(HIS2(ng)%base))
          ELSE
            status=nf90_put_att(ncid, nf90_global, 'his2_file',         &
     &                          TRIM(HIS2(ng)%name))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'his2_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef DIAGNOSTICS
        IF (exit_flag.eq.NoError) THEN
          IF (ndefDIA(ng).gt.0) THEN
            status=nf90_put_att(ncid,nf90_global, 'dia_base',           &
     &                          TRIM(DIA(ng)%base))
          ELSE
            status=nf90_put_att(ncid, nf90_global, 'dia_file',          &
     &                          TRIM(DIA(ng)%name))
          END IF
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'dia_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'err_file',            &
     &                        TRIM(ERR(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'err_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef STATIONS
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'sta_file',            &
     &                        TRIM(STA(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'sta_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef FLOATS
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'flt_file',            &
     &                        TRIM(FLT(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'flt_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifndef ANA_GRID
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'grd_file',            &
     &                        TRIM(GRD(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'grd_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef INI_FILE
# ifdef NONLINEAR
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'ini_file',            &
     &                        TRIM(INI(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'ini_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif

# ifdef TANGENT
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid,nf90_global, 'itl_file',             &
     &                        TRIM(ITL(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'itl_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif

# if defined ADJOINT && \
    !(defined AD_SENSITIVITY      || defined FOUR_DVAR        || \
      defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS || \
      defined SENSITIVITY_4DVAR   || defined SO_SEMI          || \
      defined STOCHASTIC_OPT )
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'iad_file',            &
     &                        TRIM(IAD(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'iad_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif
#endif

#if defined IS4DVAR         || defined OPT_OBSERVATIONS || \
    defined WEAK_CONSTRAINT
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'nrm_file',            &
     &                        TRIM(NRM(1,ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'nrm_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef WEAK_CONSTRAINT
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'tlf_file',            &
     &                        TRIM(TLF(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'tlf_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef FOUR_DVAR
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'obs_file',            &
     &                        TRIM(OBS(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'obs_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef FRC_FILE
        IF (exit_flag.eq.NoError) THEN
          DO i=1,nFfiles(ng)
            CALL join_string (FRC(i,ng)%files, FRC(i,ng)%Nfiles,        &
     &                        string, lstr)
            WRITE (frcatt,30) i
            status=nf90_put_att(ncid, nf90_global, frcatt,              &
     &                          string(1:lstr))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) TRIM(frcatt), TRIM(ncname)
              exit_flag=3
              ioerror=status
              EXIT
            END IF
          END DO
        END IF
#endif

        IF (ObcData(ng)) THEN
          DO i=1,nBCfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (BRY(i,ng)%files, BRY(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (bryatt,50) i
              status=nf90_put_att(ncid, nf90_global, bryatt,            &
     &                            string(1:lstr))
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,20) TRIM(bryatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF

#if !(defined ANA_SSH     || defined ANA_M2CLIMA || \
      defined ANA_M3CLIMA || defined ANA_TCLIMA)
        IF (Lclimatology(ng)) THEN
          DO i=1,nCLMfiles(ng)
            IF (exit_flag.eq.NoError) THEN
              CALL join_string (CLM(i,ng)%files, CLM(i,ng)%Nfiles,      &
     &                          string, lstr)
              WRITE (clmatt,40) i
              status=nf90_put_att(ncid, nf90_global, clmatt,            &
     &                            string(1:lstr))
              IF (status.ne.nf90_noerr) THEN
                IF (Master) WRITE (stdout,20) TRIM(clmatt), TRIM(ncname)
                exit_flag=3
                ioerror=status
              END IF
            END IF
          END DO
        END IF
#endif

#ifndef ANA_NUDGCOEF
        IF (Lnudging(ng)) THEN
          IF (exit_flag.eq.NoError) THEN
            status=nf90_put_att(ncid, nf90_global, 'nud_file',          &
     &                        TRIM(NUD(ng)%name))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) 'nud_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
#endif

#ifdef FORWARD_READ
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'fwd_file',            &
     &                        TRIM(FWD(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'fwd_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'ads_file',            &
     &                        TRIM(ADS(ng)%name))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'ads_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#if !defined DEBUGGING && defined DISTRIBUTE
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'script_file',         &
     &                        TRIM(Iname))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'script_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef FOUR_DVAR
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'apar_file',           &
     &                        TRIM(aparnam))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'apar_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef BIOLOGY
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'bpar_file',           &
     &                        TRIM(bparnam))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'bpar_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef FLOATS
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'fpos_file',           &
     &                        TRIM(fposnam))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'fpos_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifdef STATIONS
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'spos_file',           &
     &                        TRIM(sposnam))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'spos_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifndef DEBUGGING
!
!  NLM Lateral boundary conditions.
!
        IF (exit_flag.eq.NoError) THEN
          CALL lbc_putatt (ng, ncid, ncname, 'NLM_LBC', LBC, status)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'NLM_LBC', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

# if defined ADJOINT || defined TANGENT || defined TL_IOMS
!
!  Adjoint-based lateral boundary conditions.
!
        IF (exit_flag.eq.NoError) THEN
          CALL lbc_putatt (ng, ncid, ncname, 'ADM_LBC', ad_LBC, status)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'ADM_LBC', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif
#endif
!
!  SVN repository information.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_url',             &
     &                        TRIM(svn_url))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'svn_url', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

#if !defined DEBUGGING && defined SVN_REV
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'svn_rev',             &
     &                        TRIM(svn_rev))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'svn_rev', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif
#ifndef DEBUGGING
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
# ifdef ROOT_DIR
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'code_dir',            &
     &                        TRIM(Rdir))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'code_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif

# ifdef HEADER_DIR
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_dir',          &
     &                        TRIM(Hdir))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'header_dir', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif

# ifdef ROMS_HEADER
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'header_file',         &
     &                        TRIM(Hfile))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'header_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# endif
#endif
#ifndef DEBUGGING
!
!  Attributes describing platform and compiler
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'os',                  &
     &                        TRIM(my_os))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'os', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'cpu',                 &
     &                        TRIM(my_cpu))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'cpu', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_system',     &
     &                        TRIM(my_fort))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'compiler_system',            &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_command',    &
     &                        TRIM(my_fc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'compiler_command',           &
     &                                    TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'compiler_flags',      &
     &                        TRIM(my_fflags))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'compiler_flags', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Tiling and history attributes.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'tiling',              &
     &                        TRIM(tiling))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'tiling', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF

        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_att(ncid, nf90_global, 'history',             &
     &                        TRIM(history))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'history', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Analytical header files used.
!
        IF (exit_flag.eq.NoError) THEN
          CALL join_string (ANANAME, SIZE(ANANAME), string, lstr)
          IF (lstr.gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'ana_file',          &
     &                          string(1:lstr))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) WRITE (stdout,20) 'ana_file', TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
#endif

#ifdef BIOLOGY
!
!  Biology model header file used.
!
        IF (exit_flag.eq.NoError) THEN
          DO i=1,512
            bio_file(i:i)='-'
          END DO
          status=nf90_put_att(ncid, nf90_global, 'bio_file',            &
     &                        bio_file)
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'bio_file', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif

#ifndef DEBUGGING
!
!  Activated CPP options.
!
        IF (exit_flag.eq.NoError) THEN
          lstr=LEN_TRIM(Coptions)-1
          status=nf90_put_att(ncid, nf90_global, 'CPP_options',         &
     &                        TRIM(Coptions(1:lstr)))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) WRITE (stdout,20) 'CPP_options', TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
#endif
      END IF

#if !defined PARALLEL_OUT  && defined DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
      IF (exit_flag.ne.NoError) RETURN

#ifdef PROPAGATOR
!
!  Avoid writing other information variables if GST check pointing
!  NetCDF file.
!
      IF (ncid.eq.GST(ng)%ncid) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      IF (INT(time_ref).eq.-2) THEN
        Vinfo( 3)='days since 1968-05-23 00:00:00 GMT'
        Vinfo( 4)='gregorian'
      ELSE IF (INT(time_ref).eq.-1) THEN
        Vinfo( 3)='days since 0001-01-01 00:00:00'
        Vinfo( 4)='360_day'
      ELSE IF (INT(time_ref).eq.0) THEN
        Vinfo( 3)='days since 0001-01-01 00:00:00'
        Vinfo( 4)='julian'
      ELSE IF (time_ref.gt.0.0_r8) THEN
        WRITE (Vinfo( 3),'(a,1x,a)') 'days since', TRIM(r_text)
      END IF
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#if defined HDF5 && defined DEFLATE
      Vinfo( 1)='shuffle'
      Vinfo( 2)='NetCDF-4/HDF5 file format shuffle filer flag'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='deflate'
      Vinfo( 2)='NetCDF-4/HDF5 file format deflate filer flag'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='deflate_level'
      Vinfo( 2)='NetCDF-4/HDF5 file format deflate level parameter'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      Vinfo( 1)='ntsAVG'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of time-averaged fields'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nAVG'
      Vinfo( 2)='number of time-steps between time-averaged records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefAVG'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of average files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef AVERAGES2
      Vinfo( 1)='ntsAVG2'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of time-averaged fields'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nAVG2'
      Vinfo( 2)='number of time-steps between time-averaged records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefAVG2'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of average files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef HISTORY2
      Vinfo( 1)='nHIS2'
      Vinfo( 2)='number of time-steps between surface records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefHIS2'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of surface files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef ADJOINT
      Vinfo( 1)='nADJ'
      Vinfo( 2)='number of time-steps between adjoint history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefADJ'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of adjoint files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef TANGENT
      Vinfo( 1)='nTLM'
      Vinfo( 2)='number of time-steps between tangent history records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefTLM'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of tanget files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef ADJUST_BOUNDARY
      Vinfo( 1)='nOBC'
      Vinfo( 2)=                                                        &
     &  'number of time-steps between 4DVAR open boundary adjustment'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      Vinfo( 1)='nSFF'
      Vinfo( 2)=                                                        &
     &  'number of time-steps between 4DVAR surface forcing adjustment'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef PROPAGATOR
      Vinfo( 1)='LmultiGST'
      Vinfo( 2)='Switch to write one GST eigenvector per file'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LrstGST'
      Vinfo( 2)='Switch to restart GST analysis'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='MaxIterGST'
      Vinfo( 2)='maximum number of GST algorithm iterations'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nGST'
      Vinfo( 2)='number GST iterations between check pointing'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='NEV'
      Vinfo( 2)='number Ritz eigenvalues'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='NCV'
      Vinfo( 2)='number Ritz eigenvectors generated'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ritz_tol'
      Vinfo( 2)='relative accuracy of Ritz values'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef DIAGNOSTICS
      Vinfo( 1)='ntsDIA'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of diagnostic fields'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='nDIA'
      Vinfo( 2)='number of time-steps between diagnostic records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ndefDIA'
      Vinfo( 2)=                                                        &
     &   'number of time-steps between the creation of diagnostic files'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef STATIONS
      Vinfo( 1)='nSTA'
      Vinfo( 2)='number of time-steps between stations records'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#ifdef FOUR_DVAR
      Vinfo( 1)='Nouter'
      Vinfo( 2)='number of minimization outer loops'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ninner'
      Vinfo( 2)='number of minimization inner loops'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
      Vinfo( 1)='nl_tnu2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      Vinfo( 1)='ad_tnu2'
      Vinfo( 2)='adjoint model Laplacian mixing coefficient '//         &
     &          'for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      Vinfo( 1)='tl_tnu2'
      Vinfo( 2)='tangent linear model Laplacian mixing coefficient '//  &
     &          'for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined SOLVE3D && defined TS_DIF4
      Vinfo( 1)='nl_tnu4'
      Vinfo( 2)='nonlinear model biharmonic mixing coefficient '//      &
     &          'for tracers'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      Vinfo( 1)='ad_tnu4'
      Vinfo( 2)='adjoint model biharmonic mixing coefficient '//        &
     &          'for tracers'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      Vinfo( 1)='tl_tnu4'
      Vinfo( 2)='tangent linear model biharmonic mixing coefficient '// &
     &          'for tracers'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef UV_VIS2
      Vinfo( 1)='nl_visc2'
      Vinfo( 2)='nonlinear model Laplacian mixing coefficient '//       &
     &          'for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      Vinfo( 1)='ad_visc2'
      Vinfo( 2)='adjoint model Laplacian mixing coefficient '//         &
     &          'for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      Vinfo( 1)='tl_visc2'
      Vinfo( 2)='tangent linear model Laplacian mixing coefficient '//  &
     &          'for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef UV_VIS4
      Vinfo( 1)='nl_visc4'
      Vinfo( 2)='nonlinear model biharmonic mixing coefficient '//      &
     &          'for momentum'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef ADJOINT
      Vinfo( 1)='ad_visc4'
      Vinfo( 2)='adjoint model biharmonic mixing coefficient '//        &
     &          'for momentum'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# if defined TANGENT || defined TL_IOMS
      Vinfo( 1)='tl_visc4'
      Vinfo( 2)='tangent linear model biharmonic mixing coefficient '// &
     &          'for momentum'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
      Vinfo( 1)='tkenu2'
      Vinfo( 2)='harmonic mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef TKE_DIF4
      Vinfo( 1)='tkenu4'
      Vinfo( 2)='biharmonic mixing coefficient for turbulent energy'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ng, model, ncid,varid, NF_TYPE,                    &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif
#if defined UV_VIS2 || defined UV_VIS4
      Vinfo( 1)='LuvSponge'
      Vinfo( 2)='horizontal viscosity sponge activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#if (defined TS_DIF2 || defined TS_DIF4) && defined SOLVE3D
      Vinfo( 1)='LtracerSponge'
      Vinfo( 2)='horizontal diffusivity sponge activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# if defined MY25_MIXING || defined GLS_MIXING
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef FORWARD_MIXING
!
!  Basic state vertical mixing scale used in adjoint-based applications.
!
#  ifdef ADJOINT
      Vinfo( 1)='ad_Akt_fac'
      Vinfo( 2)='adjoint model basic state vertical mixing '//          &
     &          'scale for tracers'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      Vinfo( 1)='tl_Akt_fac'
      Vinfo( 2)='tangent linear model basic state vertical mixing '//   &
     &          'scale for tracers'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ADJOINT
      Vinfo( 1)='ad_Akv_fac'
      Vinfo( 2)='adjoint model basic state vertical mixing '//          &
     &          'scale for momentum'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  if defined TANGENT || defined TL_IOMS
      Vinfo( 1)='tl_Akv_fac'
      Vinfo( 2)='tangent linear model basic state vertical mixing '//   &
     &          'scale for momentum'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
#endif
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo ,ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='minimum value of specific turbulent kinetic energy'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Charnok_alpha'
      Vinfo( 2)='Charnok factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Zos_hsig_alpha'
      Vinfo( 2)='wave amplitude factor for surface roughness'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='sz_alpha'
      Vinfo( 2)='surface flux from wave dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='CrgBan_cw'
      Vinfo( 2)='surface flux due to Craig and Banner wave breaking'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#endif
#ifdef WEC
      Vinfo( 1)='wec_alpha'
      Vinfo( 2)='coefficient for roller/breaking wave dissipation'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Tnudg_SSS'
      Vinfo( 2)='SSS nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# ifdef SSSC_THRESHOLD
      Vinfo( 1)='sss_mismatch_threshold'
      Vinfo( 2)='difference in SSS that triggers SSS correction'
      Vinfo( 3)='PSU'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifndef DEBUGGING
!
!  Open boundary nudging, inverse time scales.
!
      IF (NudgingCoeff(ng)) THEN
        Vinfo( 1)='FSobc_in'
        Vinfo( 2)='free-surface inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='FSobc_out'
        Vinfo( 2)='free-surface outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='M2obc_in'
        Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='M2obc_out'
        Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
        Vinfo( 1)='Tobc_in'
        Vinfo( 2)='tracers inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='Tobc_out'
        Vinfo( 2)='tracers outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 2, tbrydim, Aval, Vinfo, ncname,                 &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='M3obc_in'
        Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='M3obc_out'
        Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
        Vinfo( 3)='second-1'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/brydim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
      END IF
#endif
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#if defined SOLVE3D && defined PROPAGATOR
      Vinfo( 1)='bvf_bak'
      Vinfo( 2)='background Brunt-Vaisala frequency'
      Vinfo( 3)='second-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif

#if defined SOLVE3D     && \
   (!defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR)
      Vinfo( 1)='R0'
      Vinfo( 2)='background density used in linear equation of state'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Tcoef'
      Vinfo( 2)='thermal expansion coefficient'
      Vinfo( 3)='Celsius-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Scoef'
      Vinfo( 2)='Saline contraction coefficient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef SOLVE3D
!
!  Various parameters.
!
# ifdef BODYFORCE
      Vinfo( 1)='levsfrc'
      Vinfo( 2)='shallowest level for body-force stress'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='levbfrc'
      Vinfo( 2)='deepest level for body-force stress'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
! Logical switches to activate horizontal momentum transport
! point Sources/Sinks (like river runoff transport) and mass point
! Sources/Sinks (like volume vertical influx).
!
      Vinfo( 1)='LuvSrc'
      Vinfo( 2)='momentum point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LwSrc'
      Vinfo( 2)='mass point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
!
!  Logical switches indicating which tracer variables are processed
!  during point Sources/Sinks.
!
      Vinfo( 1)='LtracerSrc'
      Vinfo( 2)='tracer point sources and sink activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Logical switches to process climatology fields.
!
      Vinfo( 1)='LsshCLM'
      Vinfo( 2)='sea surface height climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Lm2CLM'
      Vinfo( 2)='2D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      Vinfo( 1)='Lm3CLM'
      Vinfo( 2)='3D momentum climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LtracerCLM'
      Vinfo( 2)='tracer climatology processing switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  Logical switches for nudging of climatology fields.
!
      Vinfo( 1)='LnudgeM2CLM'
      Vinfo( 2)='2D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#ifdef SOLVE3D
!
      Vinfo( 1)='LnudgeM3CLM'
      Vinfo( 2)='3D momentum climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
      Vinfo( 1)='LnudgeTCLM'
      Vinfo( 2)='tracer climatology nudging activation switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
#ifdef FOUR_DVAR
!
!  4DVAR assimilation parameters.
!
# ifdef ADJUST_STFLUX
      Vinfo( 1)='Lstflux'
      Vinfo( 2)='surface tracer fluxes adjustment switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef ADJUST_BOUNDARY
      Vinfo( 1)='Lobc'
      Vinfo( 2)='open boundary conditions adjustment switch'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               2, (/brydim, statedim/), Aval, Vinfo, ncname,      &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifndef IS4DVAR_SENSITIVITY
      Vinfo( 1)='LhessianEV'
      Vinfo( 2)='switch to compute Hessian eigenvectors'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
      Vinfo( 1)='LhotStart'
      Vinfo( 2)='switch for hot start of subsequent outer loops'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif

      Vinfo( 1)='Lprecond'
      Vinfo( 2)='switch for conjugate gradient preconditioning'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Lritz'
      Vinfo( 2)='switch for Ritz limited-memory preconditioning'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
      IF (Lprecond.and.(NritzEV.gt.0)) THEN
        Vinfo( 1)='NritzEV'
        Vinfo( 2)='number of preconditioning eigenpairs to use'
        status=def_var(ng, model, ncid, varid, nf90_int,                &
     &                 1, (/0/), Aval, Vinfo, ncname,                   &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#  endif
# endif
# if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
      Vinfo( 1)='NpostI'
      Vinfo( 2)='number of Lanczos iterations in posterior analysis'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
# ifndef IS4DVAR_SENSITIVITY
      Vinfo( 1)='GradErr'
      Vinfo( 2)='Upper bound on relative error of the gradient'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='HevecErr'
      Vinfo( 2)='Accuracy required for Hessian eigenvectors'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

      Vinfo( 1)='Nmethod'
      Vinfo( 2)='background error covariance normalization method'
      Vinfo( 9)='exact'
      Vinfo(10)='randomization'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Rscheme'
      Vinfo( 2)='Random number generation scheme'
      Vinfo( 9)='intrisic_randon_number'
      Vinfo(10)='Gaussian_distributed_deviates'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Nrandom'
      Vinfo( 2)='number of randomization iterations'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Hgamma'
      Vinfo( 2)='initial conditions error covariance '//                &
     &          'horizontal convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef WEAK_CONSTRAINT
      Vinfo( 1)='HgammaM'
      Vinfo( 2)='model error covariance '//                             &
     &          'horizontal convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef ADJUST_BOUNDARY
      Vinfo( 1)='HgammaB'
      Vinfo( 2)='open boundary conditions error covariance '//          &
     &          'horizontal convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef ADJUST_STFLUX
      Vinfo( 1)='HgammaF'
      Vinfo( 2)='surface forcing error covariance '//                   &
     &          'horizontal convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
      Vinfo( 1)='Vgamma'
      Vinfo( 2)='initial conditions error covariance '//                &
     &          'vertical convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef WEAK_CONSTRAINT
      Vinfo( 1)='VgammaM'
      Vinfo( 2)='model error covariance '//                             &
     &          'vertical convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ADJUST_BOUNDARY
      Vinfo( 1)='VgammaB'
      Vinfo( 2)='open boundary conditions error covariance '//          &
     &          'vertical convolution time-step stability factor'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

      Vinfo( 1)='Hdecay'
      Vinfo( 2)='initial conditions error covariance '//                &
     &          'horizontal decorrelation scale'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/statedim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      Vinfo( 1)='Vdecay'
      Vinfo( 2)='initial conditions error covariance '//                &
     &           'vertical decorrelation scale'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/statedim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

      IF (NSA.eq.2) THEN
        Vinfo( 1)='HdecayM'
        Vinfo( 2)='model error covariance ' //                          &
     &            'horizontal decorrelation scale'
        Vinfo( 3)='meter'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/statedim/), Aval, Vinfo, ncname,            &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
        Vinfo( 1)='VdecayM'
        Vinfo( 2)='model error covariance '//                           &
     &            'vertical decorrelation scale'
        Vinfo( 3)='meter'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/statedim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
# endif
      END IF

# ifdef ADJUST_BOUNDARY
      Vinfo( 1)='HdecayB'
      Vinfo( 2)='open boundary conditions error covariance '//          &
     &          'horizontal decorrelation scale'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               2, (/statedim, brydim/), Aval, Vinfo, ncname,      &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      Vinfo( 1)='VdecayB'
      Vinfo( 2)='open boundary conditions error covariance '//          &
     &           'vertical decorrelation scale'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               2, (/statedim, brydim/), Aval, Vinfo, ncname,      &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef RPM_RELAXATION
      Vinfo( 1)='tl_M2diff'
      Vinfo( 2)='RPM 2D momentum diffusive relaxation coefficient'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#  ifdef SOLVE3D
      Vinfo( 1)='tl_M3diff'
      Vinfo( 2)='RPM 3D momentum diffusive relaxation coefficient'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='tl_Tdiff'
      Vinfo( 2)='RPM tracers diffusive relaxation coefficients'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

# ifdef BALANCE_OPERATOR
#  ifdef ZETA_ELLIPTIC
      Vinfo( 1)='Nbico'
      Vinfo( 2)='number of iterations in SSH elliptic equation'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#  endif

      Vinfo( 1)='Lbalance'
      Vinfo( 2)='switches for state variables included as '//           &
                'constraint in the error covariance balance operator'
      Vinfo( 9)='.FALSE.'
      Vinfo(10)='.TRUE.'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/statedim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LNM_flag'
      Vinfo( 2)='balance operator level of no motion flag'
      Vinfo( 9)='integrate from bottom to surface,'
      Vinfo(10)='integrate from LNM_depth to surface'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='LNM_depth'
      Vinfo( 2)='balance operator level of no motion depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='dTdz_min'
      Vinfo( 2)='minimum dT/dz value used in balaced salinity'
      Vinfo( 3)='Celsius meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ml_depth'
      Vinfo( 2)='mixed layer depth used in balanced salinity'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
!
!  Adjoint sensitivity parameters.
!
      Vinfo( 1)='Lzeta'
      Vinfo( 2)='adjoint sensitivity on free-surface'
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Lubar'
      Vinfo( 2)='adjoint sensitivity on 2D U-momentum'
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Lvbar'
      Vinfo( 2)='adjoint sensitivity on 2D V-momentum'
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      Vinfo( 1)='Luvel'
      Vinfo( 2)='adjoint sensitivity on 3D U-momentum'
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Lvvel'
      Vinfo( 2)='adjoint sensitivity on 3D V-momentum'
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ltracer'
      Vinfo( 2)='adjoint sensitivity on tracer variables'
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KstrS'
      Vinfo( 2)='deepest level for adjoint sensitivity analysis'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='KendS'
      Vinfo( 2)='shallowest level for adjoint sensitivity analysis'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined FORCING_SV || defined SO_SEMI || defined STOCHASTIC_OPT
!
!  Singular Forcing Vectors or Stochastic Optimals state switches.
!
      Vinfo( 1)='Fzeta'
# ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for free-surface'
# else
      Vinfo( 2)='stochastic optimals for free-surface'
# endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifndef SOLVE3D
      Vinfo( 1)='Fubar'
#  ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for 2D U-momentum'
#  else
      Vinfo( 2)='stochastic optimals for 2D U-momentum'
#  endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Fvbar'
#  ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for 2D V-momentum'
#  else
      Vinfo( 2)='stochastic optimals for 2D V-momentum'
#  endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# else

      Vinfo( 1)='Fuvel'
#  ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for 3D U-momentum'
#  else
      Vinfo( 2)='stochastic optimals for 3D U-momentum'
#  endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Fvvel'
#  ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for 3D V-momentum'
#  else
      Vinfo( 2)='stochastic optimals for 3D V-momentum'
#  endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Ftracer'
#  ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for tracer variables'
#  else
      Vinfo( 2)='stochastic optimals for tracer variables'
#  endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

      Vinfo( 1)='Fsustr'
# ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for surface U-stress'
# else
      Vinfo( 2)='stochastic optimals for surface U-stress'
# endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Fsvstr'
# ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for surface V-stress'
# else
      Vinfo( 2)='stochastic optimals for surface V-stress'
# endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      Vinfo( 1)='Fstflx'
#  ifdef FORCING_SV
      Vinfo( 2)='forcing singular vectors for surface tracer flux'
#  else
      Vinfo( 2)='stochastic optimals for surface tracer flux'
#  endif
      Vinfo( 9)='off'
      Vinfo(10)='on'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#ifdef SO_SEMI
!
!  Define Stochatic optimals parameters.
!
# ifndef SO_SEMI_WHITE
      Vinfo( 1)='SO_decay'
      Vinfo( 2)='red noise stochastic optimals time decorrelation'
      Vinfo( 3)='day'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

      Vinfo( 1)='SO_trace'
      Vinfo( 2)='trace of stochastic optimals matrix'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SOsdev_zeta'
      Vinfo( 2)='stochastic optimals scale, free-surface'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifndef SOLVE3D
      Vinfo( 1)='SOsdev_ubar'
      Vinfo( 2)='stochastic optimals scale, 2D U-momentum'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SOsdev_vbar'
      Vinfo( 2)='stochastic optimals scale, 2D V-momentum'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# else

      Vinfo( 1)='SOsdev_uvel'
      Vinfo( 2)='stochastic optimals scale, 3D U-momentum'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SOsdev_vvel'
      Vinfo( 2)='stochastic optimals scale, 3D V-momentum'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SOsdev_tracer'
      Vinfo( 2)='stochastic optimals scale, tracer variables'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif

      Vinfo( 1)='SOsdev_sustr'
      Vinfo( 2)='stochastic optimals scale, surface U-stress'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='SOsdev_svstr'
      Vinfo( 2)='stochastic optimals scale, surface V-stress'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
      Vinfo( 1)='SOsdev_stflx'
      Vinfo( 2)='stochastic optimals scale, surface tracer flux'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/trcdim/), Aval, Vinfo, ncname,                &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# endif
#endif

#if defined BIOLOGY && defined SOLVE3D
# if defined BIO_FENNEL
#  include <fennel_def.h>
# elif defined ESTUARYBGC 
#  include <estuarybgc_def.h>
# elif defined ECOSIM
#  include <ecosim_def.h>
# elif defined NEMURO
#  include <nemuro_def.h>
# elif defined BIO_UMAINE
#  include <umaine_def.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_def.h>
# elif defined NPZD_IRON
#  include <npzd_iron_def.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_def.h>
# elif defined RED_TIDE
#  include <red_tide_def.h>
# endif
#endif

#if defined FLOATS && defined FLOAT_BIOLOGY
# if defined FLOAT_OYSTER
#  include <oyster_floats_def.h>
# endif
#endif

#ifdef SEDIMENT
# include <sediment_def.h>
#endif
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian. Writing characters in
!  parallel I/O is extremely inefficient.  It is better to write
!  this as an integer switch: 0=Cartesian, 1=spherical.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 9)='Cartesian'
      Vinfo(10)='spherical'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
      Vinfo( 1)='Vtransform'
      Vinfo( 2)='vertical terrain-following transformation equation'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Vstretching'
      Vinfo( 2)='vertical terrain-following stretching function'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  SGRID conventions for staggered data on structured grids.
!
      Vinfo( 1)='grid'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_rho, scalar'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
# if defined SEDIMENT && defined SED_MORPH
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: bath depth_c: hc'
# else
      Vinfo(23)='s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc'
# endif
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_w, scalar'
      Vinfo(21)='ocean_s_coordinate'
      IF (Vtransform(ng).eq.1) THEN
        Vinfo(21)='ocean_s_coordinate_g1'
      ELSE IF (Vtransform(ng).eq.2) THEN
        Vinfo(21)='ocean_s_coordinate_g2'
      END IF
# if defined SEDIMENT && defined SED_MORPH
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: bath depth_c: hc'
# else
      Vinfo(23)='s: s_w C: Cs_w eta: zeta depth: h depth_c: hc'
# endif
      vinfo(25)='up'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_r, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/srdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_w, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/swdim/), Aval, Vinfo, ncname,                 &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(14)='user, scalar'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/usrdim/), Aval, Vinfo, ncname,              &
     &               SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#ifdef STATIONS
!
!  Station positions.
!
      IF (ncid.eq.STA(ng)%ncid) THEN
        Vinfo( 1)='Ipos'
        Vinfo( 2)='stations I-direction positions'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/stadim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='Jpos'
        Vinfo( 2)='stations J-direction positions'
        status=def_var(ng, model, ncid, varid, NF_TYPE,                 &
     &                 1, (/stadim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
      END IF
#endif
#ifdef NO_WRITE_GRID
      IF (ncid.eq.STA(ng)%ncid) THEN
#else
      IF (ncid.ne.FLT(ng)%ncid) THEN
#endif
#if !(defined SED_MORPH && defined SEDIMENT)
!
!  Bathymetry.
!
        Vinfo( 1)='h'
        Vinfo( 2)='bathymetry at RHO-points'
        Vinfo( 3)='meter'
        Vinfo(14)='bath, scalar'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
#ifdef ICESHELF
!
!  Ice shelf thickness.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='zice'
          Vinfo( 2)='Ice shelf thickness at RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='zice, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
!
!  Coriolis Parameter.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='f'
          Vinfo( 2)='Coriolis parameter at RHO-points'
          Vinfo( 3)='second-1'
          Vinfo(14)='coriolis, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='pm'
          Vinfo( 2)='curvilinear coordinate metric in XI'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pm, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='pn'
          Vinfo( 2)='curvilinear coordinate metric in ETA'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pn, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_rho'
          Vinfo( 2)='longitude of RHO-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_rho, scalar'
          Vinfo(21)='longitude'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='lat_rho'
          Vinfo( 2)='latitude of RHO-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_rho, scalar'
          Vinfo(21)='latitude'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo,  ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_rho'
          Vinfo( 2)='x-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_rho, scalar'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='y_rho'
          Vinfo( 2)='y-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_rho, scalar'
          IF (ncid.eq.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     1, (/stadim/), Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          ELSE
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Grid coordinates of U-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_u'
          Vinfo( 2)='longitude of U-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_u, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='lat_u'
          Vinfo( 2)='latitude of U-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_u, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_u'
          Vinfo( 2)='x-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_u, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='y_u'
          Vinfo( 2)='y-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_u, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Grid coordinates of V-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_v'
          Vinfo( 2)='longitude of V-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_v, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='lat_v'
          Vinfo( 2)='latitude of V-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_v, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_v'
          Vinfo( 2)='x-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_v, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='y_v'
          Vinfo( 2)='y-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_v, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
!
!  Grid coordinates of PSI-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_psi'
          Vinfo( 2)='longitude of PSI-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_psi, scalar'
          Vinfo(21)='longitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='lat_psi'
          Vinfo( 2)='latitude of PSI-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_psi, scalar'
          Vinfo(21)='latitude'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          Vinfo( 1)='x_psi'
          Vinfo( 2)='x-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_psi, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF

          Vinfo( 1)='y_psi'
          Vinfo( 2)='y-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_psi, scalar'
          IF (ncid.ne.STA(ng)%ncid) THEN
            status=def_var(ng, model, ncid, varid, NF_TYPE,             &
     &                     2, p2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
#ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points.
!
        Vinfo( 1)='angle'
        Vinfo( 2)='angle between XI-axis and EAST'
        Vinfo( 3)='radians'
        Vinfo(14)='angle, scalar'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.STA(ng)%ncid) THEN
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   1, (/stadim/), Aval, Vinfo,  ncname)
          IF (exit_flag.ne.NoError) RETURN
        ELSE
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='mask_rho'
          Vinfo( 2)='mask on RHO-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='mask_u'
          Vinfo( 2)='mask on U-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='mask_v'
          Vinfo( 2)='mask on V-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='mask_psi'
          Vinfo( 2)='mask on psi-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, p2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
!
!  Adjoint sensitivity spatial scope mask at RHO-, U-, and V-points.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
          Vinfo( 1)='scope_rho'
          Vinfo( 2)='adjoint sensitivity spatial scope on RHO-points'
          Vinfo( 9)='inactive'
          Vinfo(10)='active'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='scope_u'
          Vinfo( 2)='adjoint sensitivity spatial scope on U-points'
          Vinfo( 9)='inactive'
          Vinfo(10)='active'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='scope_v'
          Vinfo( 2)='adjoint sensitivity spatial scope on V-points'
          Vinfo( 9)='inactive'
          Vinfo(10)='active'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#endif
#ifdef UV_DRAG_GRID
!
!  Spatially varying bottom friction parameter.
!
        IF (ncid.ne.STA(ng)%ncid) THEN
# ifdef UV_LOGDRAG
          Vinfo( 1)=Vname(1,idZoBL)
          Vinfo( 2)=Vname(2,idZoBL)
          Vinfo( 3)=Vname(3,idZoBL)
          Vinfo(14)=Vname(4,idZoBL)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef UV_LDRAG
          Vinfo( 1)=Vname(1,idragL)
          Vinfo( 2)=Vname(2,idragL)
          Vinfo( 3)=Vname(3,idragL)
          Vinfo(14)=Vname(4,idragL)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef UV_QDRAG
          Vinfo( 1)=Vname(1,idragQ)
          Vinfo( 2)=Vname(2,idragQ)
          Vinfo( 3)=Vname(3,idragQ)
          Vinfo(14)=Vname(4,idragQ)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncid, varid, NF_TYPE,               &
     &                   2, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
# endif
        END IF
#endif
      END IF

  10  FORMAT (i3.3,'x',i3.3)
  20  FORMAT (/,' DEF_INFO - error while creating global attribute: ',  &
     &        a,/,12x,a)
  30  FORMAT ('frc_file_',i2.2)
  40  FORMAT ('clm_file_',i2.2)
  50  FORMAT ('bry_file_',i2.2)

      RETURN
      END SUBROUTINE def_info
