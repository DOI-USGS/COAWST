#include "cppdefs.h"
      FUNCTION def_info (ng,ncid,ncname,DimIDs)
!
!svn $Id: def_info.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     def_info   Error flag (INTEGER_TYPE).                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ncid
      integer, intent(in) :: DimIDs(29)
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 24

      integer :: brydim, i, ie, is, j, lstr, varid
      integer :: srdim, stadim, status, swdim, trcdim, usrdim
#ifdef FOUR_DVAR
      integer :: statedim
#endif
#if defined BIOLOGY && defined ECOSIM
      integer :: bacdim, domdim, fecdim, phydim

      integer :: biodim(2)
#endif

      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, def_var

      real(r8) :: Aval(6)

      character (len=11) :: frcatt
      character (len=50) :: tiling
      character (len=80) :: Vinfo(Natt)
      character (len=1024) :: ana_file
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
!
!  Set dimension for generic user parameters.
!
      IF ((Nuser.gt.0).and.(ncid.ne.ncGSTid(ng))) THEN
        status=nf90_def_dim(ncid,'Nuser',25,usrdim)
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,80
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
!  Define history global attribute.
!
      IF (LEN_TRIM(date_str).gt.0) THEN
        WRITE (history,'(a,1x,a,", ",a)') 'ROMS/TOMS, Version',         &
     &                                    TRIM( version), TRIM(date_str)
      ELSE
        WRITE (history,'(a,1x,a)') 'ROMS/TOMS, Version', TRIM(version)
      END IF
!
!  Define tile decomposition global attribute.
!
      WRITE (tiling,10) NtileI(ng), NtileJ(ng)
  10  FORMAT (i3.3,'x',i3.3)
!
!  Define file type global attribute.
!
      IF (ncid.eq.ncADJid(ng)) THEN
        type='ROMS/TOMS adjoint history file'
      ELSE IF (ncid.eq.ncAVGid(ng)) THEN
#ifdef ADJOINT
        type='ROMS/TOMS adjoint averages file'
#else
        type='ROMS/TOMS averages file'
#endif
      ELSE IF (ncid.eq.ncDIAid(ng)) THEN
        type='ROMS/TOMS diagnostics file'
      ELSE IF (ncid.eq.ncGSTid(ng)) THEN
        type='ROMS/TOMS GST check pointing restart file'
      ELSE IF (ncid.eq.ncHISid(ng)) THEN
        type='ROMS/TOMS history file'
      ELSE IF (ncid.eq.ncNRMid(ng)) THEN
        type='ROMS/TOMS background error covariance norm file'
      ELSE IF (ncid.eq.ncRSTid(ng)) THEN
        type='ROMS/TOMS restart file'
      ELSE IF (ncid.eq.ncSTAid(ng)) THEN
        type='ROMS/TOMS station file'
      ELSE IF (ncid.eq.ncTLMid(ng)) THEN
        type='ROMS/TOMS tangent linear history file'
      ELSE IF (ncid.eq.ncFLTid(ng)) THEN
        type='ROMS/TOMS floats file'
      END IF
      status=nf90_put_att(ncid, nf90_global, 'type',                    &
     &                    TRIM(type))
!
!  Define file climate and forecast metadata convention global
!  attribute.
!
      type='CF-1.0'
      status=nf90_put_att(ncid, nf90_global, 'Conventions',             &
     &                    TRIM(type))
!
!  Define other global attributes to NetCDF file.
!
      status=nf90_put_att(ncid, nf90_global, 'title',                   &
     &                    TRIM(title))
      status=nf90_put_att(ncid, nf90_global, 'var_info',                &
     &                    TRIM(varname))
#ifdef PROPAGATOR
      status=nf90_put_att(ncid, nf90_global, 'gst_file',                &
     &                    TRIM(GSTname(ng)))
#endif
      status=nf90_put_att(ncid, nf90_global, 'rst_file',                &
     &                    TRIM(RSTname(ng)))
      IF (LdefHIS(ng)) THEN
        IF (ndefHIS(ng).gt.0) THEN
          lstr=LEN_TRIM(HISbase(ng))-3
          status=nf90_put_att(ncid, nf90_global, 'his_base',            &
     &                        HISbase(ng)(1:lstr))
        ELSE
          status=nf90_put_att(ncid, nf90_global, 'his_file',            &
     &                        TRIM(HISname(ng)))
        END IF
      END IF
#ifdef AVERAGES
      IF (ndefAVG(ng).gt.0) THEN
        lstr=LEN_TRIM(AVGbase(ng))-3
        status=nf90_put_att(ncid, nf90_global, 'avg_base',              &
     &                      AVGbase(ng)(1:lstr))
      ELSE
        status=nf90_put_att(ncid, nf90_global, 'avg_file',              &
     &                      TRIM(AVGname(ng)))
      END IF
#endif
#ifdef DIAGNOSTICS
      IF (ndefDIA(ng).gt.0) THEN
        lstr=LEN_TRIM(DIAbase(ng))-3
        status=nf90_put_att(ncid,nf90_global, 'dia_base',               &
     &                      DIAbase(ng)(1:lstr))
      ELSE
        status=nf90_put_att(ncid, nf90_global, 'dia_file',              &
     &                      TRIM(DIAname(ng)))
      END IF
#endif
#ifdef STATIONS
      status=nf90_put_att(ncid, nf90_global, 'sta_file',                &
     &                    TRIM(STAname(ng)))
#endif
#ifdef FLOATS
      status=nf90_put_att(ncid, nf90_global, 'flt_file',                &
     &                    TRIM(FLTname(ng)))
#endif
#ifndef ANA_GRID
      status=nf90_put_att(ncid, nf90_global, 'grd_file',                &
     &                    TRIM(GRDname(ng)))
#endif
#ifdef INI_FILE
# ifdef NONLINEAR
      status=nf90_put_att(ncid, nf90_global, 'ini_file',                &
     &                    TRIM(INIname(ng)))
# endif
# ifdef TANGENT
      status=nf90_put_att(ncid,nf90_global,'itl_file',                  &
     &                    TRIM(ITLname(ng)))
# endif
# if defined ADJOINT && \
    !(defined AD_SENSITIVITY  || defined FOUR_DVAR        || \
      defined OBS_SENSITIVITY || defined OPT_OBSERVATIONS || \
      defined SO_SEMI)
      status=nf90_put_att(ncid, nf90_global, 'iad_file',                &
     &                    TRIM(IADname(ng)))
# endif
#endif
#if defined IS4DVAR || defined IS4DVAR_OLD || defined IOM || \
    defined OPT_OBSERVATIONS || defined WEAK_CONSTRAINT
      status=nf90_put_att(ncid, nf90_global, 'nrm_file',                &
     &                    TRIM(NRMname(ng)))
#endif
#if defined WEAK_CONSTRAINT || defined IOM
      status=nf90_put_att(ncid, nf90_global, 'tlf_file',                &
     &                    TRIM(TLFname(ng)))
#endif
#ifdef FOUR_DVAR
      status=nf90_put_att(ncid, nf90_global, 'obs_file',                &
     &                    TRIM(OBSname(ng)))
#endif
#ifdef FRC_FILE
      DO i=1,nFfiles(ng)
        WRITE (frcatt,20) i
  20    FORMAT ('frc_file_',i2.2)
        status=nf90_put_att(ncid, nf90_global, frcatt,                  &
     &                      TRIM(FRCname(i,ng)))
      END DO
#endif
#ifdef OBC_DATA
      status=nf90_put_att(ncid, nf90_global, 'bry_file',                &
     &                    TRIM(BRYname(ng)))
#endif
#if (defined  ZCLIMATOLOGY && !defined ANA_SSH)     || \
    (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
    (defined  TCLIMATOLOGY && !defined ANA_TCLIMA)  || \
    (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
      status=nf90_put_att(ncid, nf90_global, 'clm_file',                &
     &                    TRIM(CLMname(ng)))
#endif
#ifdef FORWARD_READ
      status=nf90_put_att(ncid, nf90_global, 'fwd_file',                &
     &                    TRIM(FWDname(ng)))
#endif
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
      status=nf90_put_att(ncid, nf90_global, 'ads_file',                &
     &                    TRIM(ADSname(ng)))
#endif
#if !defined DEBUGGING && defined DISTRIBUTE
      status=nf90_put_att(ncid, nf90_global, 'script_file',             &
     &                    TRIM(Iname))
#endif
#if defined ASSIMILATION || defined NUDGING
      status=nf90_put_att(ncid, nf90_global, 'apar_file',               &
     &                    TRIM(aparnam))
#endif
#ifdef BIOLOGY
      status=nf90_put_att(ncid, nf90_global, 'bpar_file',               &
     &                    TRIM(bparnam))
#endif
#ifdef FLOATS
      status=nf90_put_att(ncid, nf90_global, 'fpos_file',               &
     &                    TRIM(fposnam))
#endif
#ifdef STATIONS
      status=nf90_put_att(ncid, nf90_global, 'spos_file',               &
     &                    TRIM(sposnam))
#endif
!
!  SVN repository information.
!
      status=nf90_put_att(ncid,nf90_global,'svn_url',                   &
     &                    TRIM(svn_url))
#if !defined DEBUGGING && defined SVN_REV
      status=nf90_put_att(ncid,nf90_global,'svn_rev',                   &
     &                    TRIM(svn_rev))
#endif
!
!  Local root directory, cpp header directory and file, and analytical
!  directory
!
#ifndef DEBUGGING
# ifdef ROOT_DIR
      status=nf90_put_att(ncid,nf90_global,'code_dir',                  &
     &                    TRIM(Rdir))
# endif
# ifdef HEADER_DIR
      status=nf90_put_att(ncid,nf90_global,'header_dir',                &
     &                    TRIM(Hdir))
# endif
# ifdef ROMS_HEADER
      status=nf90_put_att(ncid,nf90_global,'header_file',               &
     &                    TRIM(Hfile))
# endif
#endif
#ifndef DEBUGGING
!
!  Attributes describing platform and compiler
!
      status=nf90_put_att(ncid,nf90_global,'os',                        &
     &                    TRIM(my_os))
      status=nf90_put_att(ncid,nf90_global,'cpu',                       &
     &                    TRIM(my_cpu))
      status=nf90_put_att(ncid,nf90_global,'compiler_system',           &
     &                    TRIM(my_fort))
      status=nf90_put_att(ncid,nf90_global,'compiler_command',          &
     &                    TRIM(my_fc))
      status=nf90_put_att(ncid,nf90_global,'compiler_flags',            &
     &                    TRIM(my_fflags))
!
!  Tiling and history attributes.
!
      status=nf90_put_att(ncid,nf90_global,'tiling',                    &
     &                    TRIM(tiling))
      status=nf90_put_att(ncid,nf90_global,'history',                   &
     &                    TRIM(history))
#endif
!
!  Analytical header files used.
!
      DO i=1,1024
        ana_file(i:i)=' '
      END DO
      is=1
      DO i=1,37
        lstr=LEN_TRIM(ANANAME(i))
        IF (lstr.gt.0) THEN
          ie=is+lstr-1
          ana_file(is:ie)=TRIM(ANANAME(i))
          is=ie+1
          ana_file(is:is)=','
          is=is+2
        END IF
      END DO
      lstr=LEN_TRIM(ana_file)-1
      IF (lstr.gt.0) THEN
        status=nf90_put_att(ncid,nf90_global,'ana_file',                &
     &                      ana_file(1:lstr))
      END IF
#ifndef DEBUGGING
!
!  Activated CPP options.
!
      status=nf90_put_att(ncid,nf90_global,'CPP_options',               &
     &                    TRIM(Coptions))
#endif
#ifdef PROPAGATOR
!
!  Avoid writing other information variables if GST check pointing
!  NetCDF file.
!
      IF (ncid.eq.ncGSTid(ng)) RETURN
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
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

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
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='nHIS'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndefHIS'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of history files'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='nRST'
      Vinfo( 2)='number of time-steps between restart records'
      IF (LcycleRST(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

#ifdef AVERAGES
      Vinfo( 1)='ntsAVG'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of time-averaged fields'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='nAVG'
      Vinfo( 2)='number of time-steps between time-averaged records'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndefAVG'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of average files'
        status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#ifdef ADJOINT
      Vinfo( 1)='nADJ'
      Vinfo( 2)='number of time-steps between adjoint history records'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndefADJ'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of adjoint files'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#ifdef TANGENT
      Vinfo( 1)='nTLM'
      Vinfo( 2)='number of time-steps between tangent history records'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndefTLM'
      Vinfo( 2)=                                                        &
     &    'number of time-steps between the creation of tanget files'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      Vinfo( 1)='nSFF'
      Vinfo( 2)=                                                        &
     &  'number of time-steps between 4DVAR surface forcing adjustment'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#ifdef PROPAGATOR
      Vinfo( 1)='LrstGST'
      Vinfo( 2)='Switch to restart GST analysis'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='MaxIterGST'
      Vinfo( 2)='maximum number of GST algorithm iterations'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='nGST'
      Vinfo( 2)='number GST iterations between check pointing'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Ritz_tol'
      Vinfo( 2)='relative accuracy of Ritz values'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif

#ifdef DIAGNOSTICS
      Vinfo( 1)='ntsDIA'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of diagnostic fields'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='nDIA'
      Vinfo( 2)='number of time-steps between diagnostic records'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndefDIA'
      Vinfo( 2)=                                                        &
     &   'number of time-steps between the creation of diagnostic files'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#ifdef STATIONS
      Vinfo( 1)='nSTA'
      Vinfo( 2)='number of time-steps between stations records'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#ifdef FOUR_DVAR
      Vinfo( 1)='Nouter'
      Vinfo( 2)='number of minimization outer loops'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Ninner'
      Vinfo( 2)='number of minimization inner loops'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
      Vinfo( 1)='tnu2'
      Vinfo( 2)='Laplacian mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#endif

#if defined SOLVE3D && defined TS_DIF4
      Vinfo( 1)='tnu4'
      Vinfo( 2)='biharmonic mixing coefficient for tracers'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#endif

#ifdef UV_VIS2
      Vinfo( 1)='visc2'
      Vinfo( 2)='Laplacian mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif
#ifdef UV_VIS4
      Vinfo( 1)='visc4'
      Vinfo( 2)='biharmonic mixing coefficient for momentum'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
      Vinfo( 1)='tkenu2'
      Vinfo( 2)='harmonic mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
# endif

# ifdef TKE_DIF4
      Vinfo( 1)='tkenu4'
      Vinfo( 2)='biharmonic mixing coefficient for turbulent energy'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
# endif
#endif
#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)

      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

# if defined MY25_MIXING || defined GLS_MIXING
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
# endif
#endif
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

#ifdef SOLVE3D
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif
#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='minimum value of specific turbulent kinetic energy'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

#ifdef SOLVE3D
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#endif
#ifdef NUDGING
!
!  Nudging inverse time scales used in data assimilation.
!
      Vinfo( 1)='Znudass'
      Vinfo( 2)='free-surface nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='M2nudass'
      Vinfo( 2)='2D momentum nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='M3nudass'
      Vinfo( 2)='3D momentum nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tnudass'
      Vinfo( 2)='Tracers nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
# endif
#endif
#ifdef NUDGING_COFF
!
!  Open boundary nudging, inverse time scales.
!
      Vinfo( 1)='FSobc_in'
      Vinfo( 2)='free-surface inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='FSobc_out'
      Vinfo( 2)='free-surface outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M2obc_in'
      Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M2obc_out'
      Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='Tobc_in'
      Vinfo( 2)='tracers inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,2,tbrydim,Aval,Vinfo,ncname)

      Vinfo( 1)='Tobc_out'
      Vinfo( 2)='tracers outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,2,tbrydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M3obc_in'
      Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M3obc_out'
      Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)
# endif
#endif
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

#if defined SOLVE3D && (defined FOUR_DVAR || defined PROPAGATOR)
      Vinfo( 1)='bvf_bak'
      Vinfo( 2)='background Brunt-Vaisala frequency'
      Vinfo( 3)='second-2'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif

#if defined SOLVE3D     && \
   (!defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR)
      Vinfo( 1)='R0'
      Vinfo( 2)='background density used in linear equation of state'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tcoef'
      Vinfo( 2)='thermal expansion coefficient'
      Vinfo( 3)='Celsius-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Scoef'
      Vinfo( 2)='Saline contraction coefficient'
      Vinfo( 3)='PSU-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#endif
#ifdef SOLVE3D
!
!  Various parameters.
!
# ifdef BODYFORCE
      Vinfo( 1)='levsfrc'
      Vinfo( 2)='shallowest level for body-force stress'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='levbfrc'
      Vinfo( 2)='deepest level for body-force stress'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
# endif
#endif
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='slipperiness parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

#ifdef FOUR_DVAR
!
!  Incremental 4DVAR parameters.
!
      Vinfo( 1)='GradErr'
      Vinfo( 2)='Upper bound on relative error of the gradient'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='HevecErr'
      Vinfo( 2)='Accuracy required for Hessian eigenvectors'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Nmethod'
      Vinfo( 2)='background error covariance normalization method'
      Vinfo( 9)='exact'
      Vinfo(10)='approximated, randomization'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Rscheme'
      Vinfo( 2)='Random number generation scheme'
      Vinfo( 9)='intrisic F90 routine randon_number'
      Vinfo(10)='Gaussian distributed deviates'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Nrandom'
      Vinfo( 2)='number of randomization iterations'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Hgamma'
      Vinfo( 2)='horizontal convolution time-step stability factor'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='Vgamma'
      Vinfo( 2)='vertical convolution time-step stability factor'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
# endif

      Vinfo( 1)='Hdecay'
      Vinfo( 2)='background covariance horizontal decorrelation scale'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,statedim,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='Vdecay'
      Vinfo( 2)='background covariance vertical decorrelation scale'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,statedim,Aval,Vinfo,ncname)
# endif

# ifdef RPM_RELAXATION
      Vinfo( 1)='tl_M2diff'
      Vinfo( 2)='RPM 2D momentum diffusive relaxation coefficient'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

#  ifdef SOLVE3D
      Vinfo( 1)='tl_M3diff'
      Vinfo( 2)='RPM 3D momentum diffusive relaxation coefficient'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='tl_Tdiff'
      Vinfo( 2)='RPM tracers diffusive relaxation coefficients'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#  endif
# endif
#endif

#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
!
!  Adjoint sensitivity parameters.
!
      Vinfo( 1)='Lzeta'
      Vinfo( 2)='adjoint sensitivity on free-surface'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Lubar'
      Vinfo( 2)='adjoint sensitivity on 2D U-momentum'
      Vinfo( 7)='on'
      Vinfo( 8)='off'

      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)
      Vinfo( 1)='Lvbar'
      Vinfo( 2)='adjoint sensitivity on 2D V-momentum'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='Luvel'
      Vinfo( 2)='adjoint sensitivity on 3D U-momentum'
      Vinfo( 7)='on'
      Vinfo( 8)='off'

      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)
      Vinfo( 1)='Lvvel'
      Vinfo( 2)='adjoint sensitivity on 3D V-momentum'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Ltracer'
      Vinfo( 2)='adjoint sensitivity on tracer variables'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,trcdim,Aval,Vinfo,ncname)

      Vinfo( 1)='KstrS'
      Vinfo( 2)='deepest level for adjoint sensitivity analysis'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='KendS'
      Vinfo( 2)='shallowest level for adjoint sensitivity analysis'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)
# endif
# ifdef SO_SEMI
!
!  Define Stochatic optimals parameters.
!
#  ifndef SO_SEMI_WHITE
      Vinfo( 1)='SO_decay'
      Vinfo( 2)='red noise stochastic optimals time decorrelation'
      Vinfo( 3)='day'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#  endif

      Vinfo( 1)='SO_trace'
      Vinfo( 2)='trace of stochastic optimals matrix'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Lso_Ustr'
      Vinfo( 2)='Stochastic optimats on surface u-stress'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Lso_Vstr'
      Vinfo( 2)='Stochastic optimats on surface v-stress'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

#  ifdef SOLVE3D
      Vinfo( 1)='Lso_tracer'
      Vinfo( 2)='Stochastic optimats on surface flux of tracer'
      Vinfo( 7)='on'
      Vinfo( 8)='off'
      status=def_var(ncid,varid,nf90_char,1,trcdim,Aval,Vinfo,ncname)
#  endif

      Vinfo( 1)='SOsdev_Ustr'
      Vinfo( 2)='stochastic optimals surface u-stress scale'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='SOsdev_Vstr'
      Vinfo( 2)='stochastic optimals surface V-stress scale'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

#  ifdef SOLVE3D
      Vinfo( 1)='SOsdev_stflx'
      Vinfo( 2)='stochastic optimals surface trace flux scales'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#  endif
# endif
#endif

#if defined BIOLOGY && defined SOLVE3D
# if defined BIO_FASHAM
!
!  Define Fasham biological model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation by seawater'
      Vinfo( 3)='meter-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='AttChl'
      Vinfo( 2)='light attenuation by chlorophyll'
      Vinfo( 3)='meter-2 milligram_Chl-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PARfrac'
      Vinfo( 2)='photosynthetically available radiation fraction'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Vp0'
      Vinfo( 2)='Eppley temperature-limited growth parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='I_thNH4'
      Vinfo( 2)='radiation threshold for nitrification'
      Vinfo( 3)='watt meter-2'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='D_p5NH4'
      Vinfo( 2)='half-saturation radiation for nitrification'
      Vinfo( 3)='watt meter-2'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='NitriR'
      Vinfo( 2)='nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_NH4'
      Vinfo( 2)='inverse half-saturation for phytoplankton NH4 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_Phy'
      Vinfo( 2)='zooplankton half-saturation constant for ingestion'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Chl2C_m'
      Vinfo( 2)='maximum chlorophyll to carbon ratio'
      Vinfo( 3)='milligram_chl milligram_carbon-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ChlMin'
      Vinfo( 2)='minimum chlorophyll threshold'
      Vinfo( 3)='milligram_chl meter-3'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyCN'
      Vinfo( 2)='phytoplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyIP'
      Vinfo( 2)='phytoplankton, NH4 inhibition parameter'
      Vinfo( 3)='millimole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton, initial slope of P-I curve'
      Vinfo( 3)='milligram_C milligram_Chl-1 watt-1 meter2 day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMin'
      Vinfo( 2)='minimum phytoplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMR'
      Vinfo( 2)='phytoplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooAE_N'
      Vinfo( 2)='zooplankton mitrogen assimilation efficiency'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooBM'
      Vinfo( 2)='rate for zooplankton basal metabolism'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooCN'
      Vinfo( 2)='zooplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooER'
      Vinfo( 2)='zooplankton specific excretion rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton maximum growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMin'
      Vinfo( 2)='minimum zooplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMR'
      Vinfo( 2)='zooplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='LDeRRN'
      Vinfo( 2)='rate of large detritus nitrogen re-mineralization'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='LDeRRC'
      Vinfo( 2)='rate of large detritus carbon re-mineralization'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='CoagR'
      Vinfo( 2)='coagulation rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='SDeRRN'
      Vinfo( 2)='remineralization rate for small detritus nitrogen'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='SDeRRC'
      Vinfo( 2)='remineralization rate for small detritus carbon'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wPhy'
      Vinfo( 2)='vertical sinking velocity for phytoplankton'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wLDet'
      Vinfo( 2)='vertical sinking velocity for large detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wSDet'
      Vinfo( 2)='vertical sinking velocity for small detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='pCO2air'
      Vinfo( 2)='partial pressure of CO2 in the air'
      Vinfo( 3)='parts per million by volume'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

# elif defined NPZD_FRANKS
!
!  Define NPZD (Franks et al., 1986) biological model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_ext'
      Vinfo( 2)='light extinction coefficient'
      Vinfo( 3)='meter-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_Phy'
      Vinfo( 2)='phytoplankton saturation coefficient'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Vm_NO3'
      Vinfo( 2)='nitrate uptake rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMR'
      Vinfo( 2)='phytoplankton senescence/mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton maximum growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMR'
      Vinfo( 2)='zooplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMD'
      Vinfo( 2)='zooplankton death bits rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooGA'
      Vinfo( 2)='zooplankton grazing inefficiency'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooEC'
      Vinfo( 2)='zooplankton excreted fraction'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='DetRR'
      Vinfo( 2)='detritus remineralization rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wDet'
      Vinfo( 2)='detrital sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

# elif defined NPZD_POWELL
!
!  Define NPZD (Powell et al., 2006) biological model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation due to sea water'
      Vinfo( 3)='meter-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='AttPhy'
      Vinfo( 2)='light attenuation due to phytoplankton'
      Vinfo( 3)='meter2/millimole_N'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton initial slopw of the P-I cureve'
      Vinfo( 3)='m2/W'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Vm_NO3'
      Vinfo( 2)='nitrate uptake rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMRD'
      Vinfo( 2)='phytoplankton mortality rate to the detritus pool'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMRN'
      Vinfo( 2)='phytoplankton mortality rate to the Nitrogen pool'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Ivlev'
      Vinfo( 2)='Ivlev constant for zooplankton grazing'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton grazing rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooEED'
      Vinfo( 2)='zooplankton excretion efficiency to detritus pool'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooEEN'
      Vinfo( 2)='zooplankton excretion efficiency to Nitrogen pool'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMRD'
      Vinfo( 2)='zooplankton mortality rate to the detritus pool'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMRN'
      Vinfo( 2)='zootoplankton mortality rate to the Nitrogen pool'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='DetRR'
      Vinfo( 2)='detritus remineralization rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wPhy'
      Vinfo( 2)='phytoplankton sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wDet'
      Vinfo( 2)='detrital sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

# elif defined ECOSIM
!
!  Define EcoSim bio-optical model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='maximum number of iterations to achieve convergence'
      status=def_var(ncid,varid,nf90_int,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='RtUVR_flag'
      Vinfo( 2)='switch to calculate CDOC UV photolysis.'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='NFIX_flag'
      Vinfo( 2)='switch to calculate temperature based N fixation'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Regen_flag'
      Vinfo( 2)='switch to calculate fecal matter regeneration'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='HsNO3'
      Vinfo( 2)='half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='micromole_NO3 liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsNH4'
      Vinfo( 2)='half-saturation for phytoplankton NH4 uptake'
      Vinfo( 3)='micromole_NH4 liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsSiO'
      Vinfo( 2)='half-saturation for phytoplankton SiO uptake'
      Vinfo( 3)='micromole_SiO liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsPO4'
      Vinfo( 2)='half-saturation for phytoplankton PO4 uptake'
      Vinfo( 3)='micromole_PO4 liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsFe'
      Vinfo( 2)='half-saturation for phytoplankton Fe uptake'
      Vinfo( 3)='micromole_Fe liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='GtALG_max'
      Vinfo( 2)='maximum phytoplankton 24 hour growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyTbase'
      Vinfo( 2)=                                                        &
     &         'phytoplankton temperature base for exponential response'
      Vinfo( 3)='Celsius'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyTfac'
      Vinfo( 2)='phytoplankton exponential temperature factor'
      Vinfo( 3)='Celsius-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='BET_'
      Vinfo( 2)='nitrate uptake inhibition for NH4'
      Vinfo( 3)='micromole-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='maxC2nALG'
      Vinfo( 2)='maximum phytoplankton C:N ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='minC2nALG'
      Vinfo( 2)='balanced phytoplankton C:N ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2nALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:N ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='maxC2SiALG'
      Vinfo( 2)='maximum phytoplankton C:Si ratio'
      Vinfo( 3)='micromole_C micromole_Si-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='minC2SiALG'
      Vinfo( 2)='balanced phytoplankton C:Si ratio'
      Vinfo( 3)='micromole_C micromole_Si-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2SiALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:Si ratio'
      Vinfo( 3)='micromole_C micromole_Si-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='maxC2pALG'
      Vinfo( 2)='maximum phytoplankton C:P ratio'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='minC2pALG'
      Vinfo( 2)='balanced phytoplankton C:P ratio'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2pALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:P ratio'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='maxC2FeALG'
      Vinfo( 2)='maximum phytoplankton C:Fe ratio'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='minC2FeALG'
      Vinfo( 2)='balanced phytoplankton C:Fe ratio'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2FeALGminABS'
      Vinfo( 2)='absolute minimum phytoplankton C:Fe ratio'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='qu_yld'
      Vinfo( 2)='maximum quantum yield'
      Vinfo( 3)='micromole_C micromole_quanta-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='E0_comp'
      Vinfo( 2)='compensation light level'
      Vinfo( 3)='micromole_quanta'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='E0_inhib'
      Vinfo( 2)='light level for onset of photoinhibition'
      Vinfo( 3)='micromole_quanta'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='inhib_fac'
      Vinfo( 2)='exponential decay factor for light limited growth'
      Vinfo( 3)='micromole_quanta-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2CHL_max'
      Vinfo( 2)='maximum lighted limited C:Chl ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxC2Cl'
      Vinfo( 2)='rate of change in light limited C:CHL ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1 micromole_quanta-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_C2Cl'
      Vinfo( 2)='minimum lighted limited C:Chl ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxC2Cn'
      Vinfo( 2)='rate of change in nutient limited C:Chl ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1 micromole_N micromole_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_C2Cn'
      Vinfo( 2)='minimum nutrient limited C:CHL ratio'
      Vinfo( 3)='microgram_C microgram_Chl-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxPacEff'
      Vinfo( 2)='rate of change in package effect'
      Vinfo( 3)='microgram_Chl microgram_C-1' 
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_PacEff'
      Vinfo( 2)='maximum package effect'
      Vinfo( 3)='microgram_Chl microgram_C-1' 
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxChlB'
      Vinfo( 2)='rate of change in the Chl_b:CHL_a ratio'
      Vinfo( 3)=                                                        &
     & 'microgram_Chl_b microgram_Chl_a-1 microgram_Chl_a microgrma_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_ChlB'
      Vinfo( 2)='maximum Chl_b:Chl_a ratio'
      Vinfo( 3)='microgram_Chl_b microgram_Chl_a-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxChlC'
      Vinfo( 2)='rate of change in the Chl_c:Chl_a ratio'
      Vinfo( 3)=                                                        &
     & 'microgram_Chl_c microgram_Chl_a-1 microgram_Chl_a microgrma_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_ChlC'
      Vinfo( 2)='maximum Chl_c:Chl_a ratio'
      Vinfo( 3)='microgram_Chl_c microgram_Chl_a-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxPSC'
      Vinfo( 2)='rate of change in the PSC:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &   'microgram_PSC microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_PSC'
      Vinfo( 2)='maximum PSC:Chl_a ratio'
      Vinfo( 3)='microgram_PSC microgram_Chl_a-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxPPC'
      Vinfo( 2)='rate of change in the PPC:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &   'microgram_PPC microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_PPC'
      Vinfo( 2)='maximum PPC:Chl_a ratio'
      Vinfo( 3)='microgram_Chl_c microgram_Chl_a-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxLPUb'
      Vinfo( 2)='rate of change in the LPUb:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &  'microgram_LPUb microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_LPUb'
      Vinfo( 2)='Maximum LPUb:Chl_a ratio'
      Vinfo( 3)='microgram_HPUb microgram_Chl_a-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='mxHPUb'
      Vinfo( 2)='rate of change in the HPUb:Chl_a ratio'
      Vinfo( 3)=                                                        &
     &  'microgram_HPUb microgram_Chl_a-1 microgram_Chl_a microgram_C-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='b_HPUb'
      Vinfo( 2)='maximum HPUb:Chl_a ratio'
      Vinfo( 3)='microgram_HPUb microgram_Chl_a-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='FecDOC'
      Vinfo( 2)='proportion of grazing stress apportioned to DOM'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='FecPEL'
      Vinfo( 2)='proportion of grazing stress apportioned to fecal'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='FecCYC'
      Vinfo( 2)='proportion of grazing stress that is recycled'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='ExALG'
      Vinfo( 2)='proportion of daily production lost to excretion'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='WS'
      Vinfo( 2)='phytoplankton sinking speed'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsGRZ'
      Vinfo( 2)='phytoplankton grazing parameter'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='MinRefuge'
      Vinfo( 2)='refuge phytoplankton population'
      Vinfo( 3)='micromole_C liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='RefugeDep'
      Vinfo( 2)='maximum refuge phytoplankton depth'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='Norm_Vol'
      Vinfo( 2)='normalized volume factor'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='norm_Surf'
      Vinfo( 2)='normalized surface area factor'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsDOP'
      Vinfo( 2)='half-saturation constant for DOP uptake'
      Vinfo( 3)='micromole_DOP liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2pALKPHOS'
      Vinfo( 2)='C:P ratio where DOP uptake begins'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsDON'
      Vinfo( 2)='Half Saturation Constant for DON uptake'
      Vinfo( 3)='micromole_DON liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2nNupDON'
      Vinfo( 2)='C:N ratio where DON uptake begins'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,phydim,Aval,Vinfo,ncname)

      Vinfo( 1)='HsDOC_ba'
      Vinfo( 2)='half-saturation constant for bacteria DOC uptake'
      Vinfo( 3)='micromole_DOC liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,bacdim,Aval,Vinfo,ncname)

      Vinfo( 1)='GtBAC_max'
      Vinfo( 2)='maximum 24 hour bacterial growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,bacdim,Aval,Vinfo,ncname)

      Vinfo( 1)='BacTbase'
      Vinfo( 2)='bacteria temperature base for exponential response'
      Vinfo( 3)='Celsius'
      status=def_var(ncid,varid,NF_TYPE,1,bacdim,Aval,Vinfo,ncname)

      Vinfo( 1)='BacTfac'
      Vinfo( 2)='bacteria exponential temperature factor'
      Vinfo( 3)='Celsius'
      status=def_var(ncid,varid,NF_TYPE,1,bacdim,Aval,Vinfo,ncname)

      Vinfo( 1)='C2nBAC'
      Vinfo( 2)='carbon to nitrogen ratio of bacteria'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='C2pBAC'
      Vinfo( 2)='carbon to phosphorus ratio of bacteria'
      Vinfo( 3)='micromole_C micromole_P-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='C2FeBAC'
      Vinfo( 2)='carbon to iron ratio of bacteria'
      Vinfo( 3)='micromole_C micromole_Fe-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='BacDOC'
      Vinfo( 2)=                                                        &
     &        'proportion of bacteria grazing stress apportioned to DOM'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='BacPEL'
      Vinfo( 2)=                                                        &
     &      'proportion of bacteria grazing stress apportioned to fecal'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='BacCYC'
      Vinfo( 2)='proportion of bacteria grazing stress recycled'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ExBAC_c'
      Vinfo( 2)=                                                        &
     &      'bacterial recalcitrant C excretion as proportion of uptake'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ExBacC2N'
      Vinfo( 2)='bacterial recalcitrant excretion carbon:nitrogen ratio'
      Vinfo( 3)='micromole_C micromole_N-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Bac_Ceff'
      Vinfo( 2)='bacterial gross growth carbon efficiency'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='RtNIT'
      Vinfo( 2)='maximum bacterial nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='HsNIT'
      Vinfo( 2)='half-saturation constant for bacterial nitrification'
      Vinfo( 3)='micromole_NH4 liter-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='cDOCfrac_c'
      Vinfo( 2)=                                                        &
     &  'colored fraction of DOC from phytoplankton and bacteria losses'
      status=def_var(ncid,varid,NF_TYPE,1,domdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RtUVR_DIC'
      Vinfo( 2)='UV degradation of DOC into DIC at 410 nanometer'
      Vinfo( 3)='micromole meter-1 liter-1 hour-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='RtUVR_DOC'
      Vinfo( 2)=                                                        &
     &  'UV degradation of DOC to colorless labile DOC at 410 nanometer'
      Vinfo( 3)='micromole meter-1 liter-1 hour-1'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='WF'
      Vinfo( 2)='fecal sinking flux'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegTbase'
      Vinfo( 2)=                                                        &
     &    'fecal regeneration temperature base for exponential response'
      Vinfo( 3)='Celsius'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegTfac'
      Vinfo( 2)='fecal regeneration exponential temperature factor'
      Vinfo( 3)='Celsius-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegCR'
      Vinfo( 2)='fecal carbon regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegNR'
      Vinfo( 2)='fecal nitrogen regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegSR'
      Vinfo( 2)='fecal silica regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegPR'
      Vinfo( 2)='fecal phosphorus regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

      Vinfo( 1)='RegFR'
      Vinfo( 2)='fecal iron regeneration rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,fecdim,Aval,Vinfo,ncname)

# endif
#endif
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 7)='spherical'
      Vinfo( 8)='Cartesian'
      status=def_var(ncid,varid,nf90_char,1,0,Aval,Vinfo,ncname)
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)

      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,1,0,Aval,Vinfo,ncname)
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='s_rho'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_rho, scalar'
      Vinfo(21)='ocean_s_coordinate'
#if defined SEDIMENT && defined SED_MORPH
      Vinfo(23)=                                                        &
     &'s: s_rho eta: zeta depth: bath a: theta_s b: theta_b depth_c: hc'
#else
      Vinfo(23)=                                                        &
     &   's: s_rho eta: zeta depth: h a: theta_s b: theta_b depth_c: hc'
#endif
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ncid,varid,NF_TYPE,1,srdim,Aval,Vinfo,ncname)
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='s_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='s_w, scalar'
      Vinfo(21)='ocean_s_coordinate'
#if defined SEDIMENT && defined SED_MORPH
      Vinfo(23)=                                                        &
     &  's: s_w eta: zeta depth: bath a: theta_s b: theta_b depth_c: hc'
#else
      Vinfo(23)=                                                        &
     &     's: s_w eta: zeta depth: h a: theta_s b: theta_b depth_c: hc'
#endif
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ncid,varid,NF_TYPE,1,swdim,Aval,Vinfo,ncname)
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
      status=def_var(ncid,varid,NF_TYPE,1,srdim,Aval,Vinfo,ncname)
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
      status=def_var(ncid,varid,NF_TYPE,1,swdim,Aval,Vinfo,ncname)
#endif
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(14)='user, scalar'
        status=def_var(ncid,varid,NF_TYPE,1,usrdim,Aval,Vinfo,ncname)
      END IF
#ifdef STATIONS
!
!  Station positions.
!
      IF (ncid.eq.ncSTAid(ng)) THEN
        Vinfo( 1)='Ipos'
        Vinfo( 2)='stations I-direction positions'
        status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,          &
     &                 ncname)

        Vinfo( 1)='Jpos'
        Vinfo( 2)='stations J-direction positions'
        status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,          &
     &                 ncname)
      END IF
#endif
#ifdef NO_WRITE_GRID
      IF (ncid.eq.ncSTAid(ng)) THEN
#else
      IF (ncid.ne.ncFLTid(ng)) THEN
#endif
# if !(defined SED_MORPH && defined SEDIMENT)
!
!  Bathymetry.
!
        Vinfo( 1)='h'
        Vinfo( 2)='bathymetry at RHO-points'
        Vinfo( 3)='meter'
        Vinfo(14)='bath, scalar'
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        IF (ncid.eq.ncSTAid(ng)) THEN
          status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,        &
     &                   ncname)
        ELSE
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
# endif
!
!  Coriolis Parameter.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          Vinfo( 1)='f'
          Vinfo( 2)='Coriolis parameter at RHO-points'
          Vinfo( 3)='second-1'
          Vinfo(14)='coriolis, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          Vinfo( 1)='pm'
          Vinfo( 2)='curvilinear coordinate metric in XI'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pm, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='pn'
          Vinfo( 2)='curvilinear coordinate metric in ETA'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pn, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
!
!  Grid coordinates of RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_rho'
          Vinfo( 2)='longitude of RHO-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_rho, scalar'
          IF (ncid.eq.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF

          Vinfo( 1)='lat_rho'
          Vinfo( 2)='latitude of RHO-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_rho, scalar'
          IF (ncid.eq.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        ELSE
          Vinfo( 1)='x_rho'
          Vinfo( 2)='x-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_rho, scalar'
          IF (ncid.eq.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='y_rho'
          Vinfo( 2)='y-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_rho, scalar'
          IF (ncid.eq.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
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
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='lat_u'
          Vinfo( 2)='latitude of U-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_u, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        ELSE
          Vinfo( 1)='x_u'
          Vinfo( 2)='x-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_u, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='y_u'
          Vinfo( 2)='y-locations of U-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_u, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,      &
     &                     ncname)
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
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='lat_v'
          Vinfo( 2)='latitude of V-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_v, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        ELSE
          Vinfo( 1)='x_v'
          Vinfo( 2)='x-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_v, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='y_v'
          Vinfo( 2)='y-locations of V-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_v, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,      &
     &                     ncname)
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
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,p2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='lat_psi'
          Vinfo( 2)='latitude of PSI-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_psi, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,p2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        ELSE
          Vinfo( 1)='x_psi'
          Vinfo( 2)='x-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_psi, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,p2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='y_psi'
          Vinfo( 2)='y-locations of PSI-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_psi, scalar'
          IF (ncid.ne.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,2,p2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        END IF
# ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points.
!
        IF (spherical) THEN
          Vinfo( 1)='angle'
          Vinfo( 2)='angle between XI-axis and EAST'
          Vinfo( 3)='radians'
          Vinfo(14)='angle, scalar'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          IF (ncid.eq.ncSTAid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        END IF
# endif
# ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          Vinfo( 1)='mask_rho'
          Vinfo( 2)='mask on RHO-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='mask_u'
          Vinfo( 2)='mask on U-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='mask_v'
          Vinfo( 2)='mask on V-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='mask_psi'
          Vinfo( 2)='mask on psi-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,p2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
# endif
# if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
     defined OPT_OBSERVATIONS || defined SO_SEMI
!
!  Adjoint sensitivity spatial scope mask at RHO-, U-, and V-points.
!
        IF (ncid.ne.ncSTAid(ng)) THEN
          Vinfo( 1)='scope_rho'
          Vinfo( 2)='adjoint sensitivity spatial scope on RHO-points'
          Vinfo( 9)='inactive'
          Vinfo(10)='active'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='scope_u'
          Vinfo( 2)='adjoint sensitivity spatial scope on U-points'
          Vinfo( 9)='inactive'
          Vinfo(10)='active'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='scope_v'
          Vinfo( 2)='adjoint sensitivity spatial scope on V-points'
          Vinfo( 9)='inactive'
          Vinfo(10)='active'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
# endif
      END IF
      def_info=status
      RETURN
      end function def_info
