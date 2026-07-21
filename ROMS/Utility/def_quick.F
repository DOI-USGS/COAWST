#include "cppdefs.h"
      MODULE def_quick_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module creates output QUICKSAVE file using either the standard !
!  NetCDF library or the Parallel-IO (PIO) library.  It defines its    !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
#ifdef SEDIMENT
      USE mod_sediment
#endif
#ifdef VEGETATION
      USE mod_vegetation
#endif
!
#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
      USE bbl_output_mod,        ONLY : bbl_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,        ONLY : bbl_def_pio
# endif
#endif
      USE def_dim_mod,           ONLY : def_dim
      USE def_info_mod,          ONLY : def_info
      USE def_var_mod,           ONLY : def_var
#if defined ICE_MODEL && defined SOLVE3D
      USE ice_output_mod,        ONLY : ice_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,        ONLY : ice_def_pio
# endif
#endif
#if defined SEDIMENT && defined SOLVE3D
      USE sediment_output_mod,   ONLY : sediment_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod,   ONLY : sediment_def_pio
# endif
#endif
      USE strings_mod,           ONLY : FoundError
#if defined VEGETATION && defined SOLVE3D
      USE vegetation_output_mod, ONLY : vegetation_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE vegetation_output_mod, ONLY : vegetation_def_pio
# endif
#endif
#if defined WEC_VF && defined SOLVE3D
      USE wec_output_mod,        ONLY : wec_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,        ONLY : wec_def_pio
# endif
#endif
      USE wrt_info_mod,          ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_quick
      PRIVATE :: def_quick_nf90
#if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_quick_pio
#endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_quick (ng, ldef)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Create a new history file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (QCK(ng)%IOtype)
        CASE (io_nf90)
          CALL def_quick_nf90 (ng, iNLM, ldef)

#if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_quick_pio (ng, iNLM, ldef)
#endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) QCK(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_QUICK - Illegal output file type, io_type = ',i0,   &
     &        /,13x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_quick
!
!***********************************************************************
      SUBROUTINE def_quick_nf90 (ng, model, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, ifield, itrc, nvd3, nvd4, varid
      integer :: recdim, status
#ifdef ADJUST_BOUNDARY
      integer :: IorJdim
#endif
      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)

#ifdef SOLVE3D
# ifdef SEDIMENT
      integer :: b3dgrd(4)
# endif
# if defined VEGETATION 
      integer :: sapgrd(4)
# endif
      integer :: p3dgrd(4),  t3dgrd(4),  u3dgrd(4),  v3dgrd(4)
      integer :: p3dzgrd(4), t3dzgrd(4), u3dzgrd(4), v3dzgrd(4)
      integer :: w3dgrd(4)
#endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_quick_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=QCK(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new quicksave file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), QCK(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xi_rho',       &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xi_u',         &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xi_v',         &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xi_psi',       &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'eta_rho',      &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'eta_u',        &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'eta_v',        &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'eta_psi',      &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'IorJ',         &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xy_rho',       &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xy_u',         &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xy_v',         &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef SOLVE3D
# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xyz_rho',      &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xyz_u',        &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xyz_v',        &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xyz_w',        &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'N',            &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        IF (Nslice.gt.0) THEN
          status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'z_slice',    &
     &                   Nslice, DimIDs(34))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SEDIMENT
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'NST',          &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'xybed',        &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

# if defined VEGETATION
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nveg',         &
     &                 NVEG, DimIDs(35))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ECOSIM
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nbands',       &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif

        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef FOUR_DVAR
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef MULTI_SCALE_B
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'Nscale',       &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, QCK(ng)%ncid, ncname, 'obc_adjust',   &
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

        status=def_dim(ng, model, QCK(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
#if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#else
        nvd3=3
        nvd4=4
#endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
# endif
#else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
!
        t3dzgrd(1)=DimIDs( 1)
        t3dzgrd(2)=DimIDs( 5)
        t3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        t3dzgrd(4)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
#ifdef SOLVE3D
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
!
        p3dzgrd(1)=DimIDs( 4)
        p3dzgrd(2)=DimIDs( 8)
        p3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        p3dzgrd(4)=DimIDs(12)
#endif
#if defined VEGETATION && defined SOLVE3D
!
!  Define dimensions for the submerged aquatic plant (SAP) variables.
!
        sapgrd(1)=DimIDs( 1)
        sapgrd(2)=DimIDs( 5)
        sapgrd(3)=DimIDs(35)
        sapgrd(4)=DimIDs(12)
#endif 
!
!  Define dimension vectors for staggered u-momentum type variables.
!
#if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
# endif
#else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
!
        u3dzgrd(1)=DimIDs( 2)
        u3dzgrd(2)=DimIDs( 6)
        u3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        u3dzgrd(4)=DimIDs(12)
# endif
#endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
#if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
# endif
#else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
!
        v3dzgrd(1)=DimIDs( 3)
        v3dzgrd(2)=DimIDs( 7)
        v3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        v3dzgrd(4)=DimIDs(12)
# endif
#endif
#ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
# else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
# endif
# ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#  if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#  else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#  endif
# endif
#endif
!
!  Initialize unlimited time record dimension.
!
        QCK(ng)%Rindex=0
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, QCK(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idtime),    &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#ifdef WET_DRY
!
!  Define wet/dry mask on PSI-points.
!
        Vinfo( 1)=Vname(1,idPwet)
        Vinfo( 2)=Vname(2,idPwet)
        Vinfo( 3)=Vname(3,idPwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idPwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(21)=Vname(6,idPwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idPwet,ng),r8)
        status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idPwet),    &
     &                 NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define wet/dry mask on RHO-points.
!
        Vinfo( 1)=Vname(1,idRwet)
        Vinfo( 2)=Vname(2,idRwet)
        Vinfo( 3)=Vname(3,idRwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idRwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(21)=Vname(6,idRwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRwet,ng),r8)
        status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idRwet),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define wet/dry mask on U-points.
!
        Vinfo( 1)=Vname(1,idUwet)
        Vinfo( 2)=Vname(2,idUwet)
        Vinfo( 3)=Vname(3,idUwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idUwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(21)=Vname(6,idUwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUwet,ng),r8)
        status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUwet),    &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define wet/dry mask on V-points.
!
        Vinfo( 1)=Vname(1,idVwet)
        Vinfo( 2)=Vname(2,idVwet)
        Vinfo( 3)=Vname(3,idVwet)
        Vinfo(14)=Vname(4,idVwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(21)=Vname(6,idVwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVwet,ng),r8)
        status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVwet),    &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef SOLVE3D
!
!  Define time-varying depth of RHO-points.
!
        IF (Qout(idpthR,ng)) THEN
          Vinfo( 1)=Vname(1,idpthR)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthR))
          Vinfo( 3)=Vname(3,idpthR)
          Vinfo(14)=Vname(4,idpthR)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idpthR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthR,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idpthR),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Qout(idpthU,ng)) THEN
          Vinfo( 1)=Vname(1,idpthU)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthU))
          Vinfo( 3)=Vname(3,idpthU)
          Vinfo(14)=Vname(4,idpthU)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idpthU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthU,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idpthU),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Qout(idpthV,ng)) THEN
          Vinfo( 1)=Vname(1,idpthV)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthV))
          Vinfo( 3)=Vname(3,idpthV)
          Vinfo(14)=Vname(4,idpthV)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idpthV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthV,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idpthV),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Qout(idpthW,ng)) THEN
          Vinfo( 1)=Vname(1,idpthW)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthW))
          Vinfo( 3)=Vname(3,idpthW)
          Vinfo(14)=Vname(4,idpthW)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idpthW)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthW,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idpthW),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#endif
!
!  Define free-surface.
!
        IF (Qout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idFsur),  &
#ifdef WET_DRY
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
#else
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
#endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D U-momentum component.
!
        IF (Qout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUbar),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D V-momentum component.
!
        IF (Qout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVbar),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Qout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idu2dE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Qout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idv2dN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        IF (Qout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUvel),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component.
!
        IF (Qout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVvel),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface U-momentum component.
!
        IF (Qout(idUsur,ng)) THEN
          Vinfo( 1)=Vname(1,idUsur)
          Vinfo( 2)=Vname(2,idUsur)
          Vinfo( 3)=Vname(3,idUsur)
          Vinfo(14)=Vname(4,idUsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsur,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUsur),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface V-momentum component.
!
        IF (Qout(idVsur,ng)) THEN
          Vinfo( 1)=Vname(1,idVsur)
          Vinfo( 2)=Vname(2,idVsur)
          Vinfo( 3)=Vname(3,idVsur)
          Vinfo(14)=Vname(4,idVsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsur,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVsur),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D U-momentum component at specified constant depth slices.
!
        IF (Qout(idUzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsl)
          Vinfo( 2)=Vname(2,idUzsl)
          Vinfo( 3)=Vname(3,idUzsl)
          Vinfo(14)=Vname(4,idUzsl)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsl,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUzsl),  &
     &                   NF_FOUT, nvd4, u3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component at specified constant depth slices.
!
        IF (Qout(idVzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsl)
          Vinfo( 2)=Vname(2,idVzsl)
          Vinfo( 3)=Vname(3,idVzsl)
          Vinfo(14)=Vname(4,idVzsl)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idVzsl,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVzsl),  &
     &                   NF_FOUT, nvd4, v3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Qout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idu3dE),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Qout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idv3dN),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface Eastward momentum component at RHO-points.
!
        IF (Qout(idUsuE,ng)) THEN
          Vinfo( 1)=Vname(1,idUsuE)
          Vinfo( 2)=Vname(2,idUsuE)
          Vinfo( 3)=Vname(3,idUsuE)
          Vinfo(14)=Vname(4,idUsuE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idUsuE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsuE,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUsuE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface Northward momentum component at RHO-points.
!
        IF (Qout(idVsuN,ng)) THEN
          Vinfo( 1)=Vname(1,idVsuN)
          Vinfo( 2)=Vname(2,idVsuN)
          Vinfo( 3)=Vname(3,idVsuN)
          Vinfo(14)=Vname(4,idVsuN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idVsuN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsuN,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVsuN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-point for specified
!  constant depth slices.
!
        IF (Qout(idUzsE,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsE)
          Vinfo( 2)=Vname(2,idUzsE)
          Vinfo( 3)=Vname(3,idUzsE)
          Vinfo(14)=Vname(4,idUzsE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idUzsE)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsE,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUzsE),  &
     &                   NF_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points for specified
!  constant depth slices.
!
        IF (Qout(idVzsN,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsN)
          Vinfo( 2)=Vname(2,idVzsN)
          Vinfo( 3)=Vname(3,idVzsN)
          Vinfo(14)=Vname(4,idVzsN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idVzsN)
          Vinfo(22)='coordinates, zslice'
          Aval(5)=REAL(Iinfo(1,idVzsN,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVzsN),  &
     &                   NF_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model potential vorticity at PSI-points for specified constant
!  depth slices.
!
        IF (Qout(idPVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idPVzs)
          Vinfo( 2)=Vname(2,idPVzs)
          Vinfo( 3)=Vname(3,idPVzs)
          Vinfo(14)=Vname(4,idPVzs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
# endif
          Vinfo(21)=Vname(6,idPVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idPVzs,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idPVzs),  &
     &                   NF_FOUT, nvd3, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model relative vorticity at PSI-points for specified constant
!  depth slices.
!
        IF (Qout(idRVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idRVzs)
          Vinfo( 2)=Vname(2,idRVzs)
          Vinfo( 3)=Vname(3,idRVzs)
          Vinfo(14)=Vname(4,idRVzs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
# endif
          Vinfo(21)=Vname(6,idRVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idRVzs,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idRVzs),  &
     &                   NF_FOUT, nvd3, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Qout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idWvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Qout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idOvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Qout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTvar(itrc),ng),r8)
            status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Tid(itrc),  &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define surface tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Qout(idsurT(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idsurT(itrc))
            Vinfo( 2)=Vname(2,idsurT(itrc))
            Vinfo( 3)=Vname(3,idsurT(itrc))
            Vinfo(14)=Vname(4,idsurT(itrc))
            Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(21)=Vname(6,idsurT(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, model, QCK(ng)%ncid,                     &
     &                     QCK(ng)%Vid(idsurT(itrc)),                   &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define tracer type variables at specified depth slices.
!
        DO itrc=1,NT(ng)
          IF (Qout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
            Vinfo( 1)=Vname(1,idzslT(itrc))
            Vinfo( 2)=Vname(2,idzslT(itrc))
            Vinfo( 3)=Vname(3,idzslT(itrc))
            Vinfo(14)=Vname(4,idzslT(itrc))
            Vinfo(16)=Vname(1,idtime)
# ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
# endif
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(21)=Vname(6,idzslT(itrc))
            Vinfo(22)='coordinates, z_slice'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, QCK(ng)%ncid,                     &
     &                     QCK(ng)%Vid(idzslT(itrc)),                   &
     &                     NF_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Qout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idDano),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Qout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idHsbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idHsbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Qout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idHbbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idHbbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define vertical viscosity coefficient.
!
        IF (Qout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVvis),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Qout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idTdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Qout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idSdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Qout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idMtke)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idMtke),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Qout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idMtls)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idMtls),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Qout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idPair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif 
# if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
!
!  Define atmosphere boundary layer thickness.
!
        IF (Qout(idPBLh,ng)) THEN
          Vinfo( 1)=Vname(1,idPBLh)
          Vinfo( 2)=Vname(2,idPBLh)
          Vinfo( 3)=Vname(3,idPBLh)
          Vinfo(14)=Vname(4,idPBLh)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idPBLh)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPBLh,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idPBLh),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Qout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Qout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface winds at RHO-points.
!
        IF (Qout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          Vinfo( 2)=Vname(2,idUaiE)
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUaiE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Qout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          Vinfo( 2)=Vname(2,idVaiN)
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVaiN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Qout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
# endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            status=def_var(ng, model, QCK(ng)%ncid,                     &
     &                     QCK(ng)%Vid(idTsur(itrc)), NF_FOUT,          &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

# if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Define latent heat flux.
!
        IF (Qout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Qout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Qout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif

# if defined BULK_FLUXES
!
!  Define atmospheric air temperature.
!
        IF (Qout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idTair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Qout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idevap),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Qout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idrain),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
!
!  Define E-P flux.
!
        IF (Qout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idEmPf)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idEmPf),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SHORTWAVE
!
!  Define net shortwave radiation flux.
!
        IF (Qout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idSrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
#endif
!
!  Define surface U-momentum stress.
!
        IF (Qout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUsms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Qout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVsms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Qout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#endif
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idUbms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Qout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#endif
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, QCK(ng)%ncid, QCK(ng)%Vid(idVbms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_nf90 (ng, model, ldef, Qout, QCK,                  &
     &                     t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_nf90 (ng, model, ldef, Qout, QCK,                  &
     &                     t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_nf90 (ng, model, ldef, Qout, QCK,             &
     &                          t2dgrd, u2dgrd, v2dgrd,                 &
     &                          b3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined VEGETATION && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the submerged aquatic vegetation model variables.
!-----------------------------------------------------------------------
!
        CALL vegetation_def_nf90 (ng, model, ldef, Qout, QCK,           &
     &                            sapgrd, t2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
        CALL wec_def_nf90 (ng, model, ldef, Qout, QCK,                  &
     &                     t2dgrd, u2dgrd, v2dgrd,                      &
     &                     t3dgrd, u3dgrd, v3dgrd, w3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, QCK(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, QCK(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing quicksave file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=QCK(ng)%name
!
!  Open quicksave file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, QCK(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,60) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = QCK(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = QCK(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  quicksave variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            QCK(ng)%Vid(idtime)=var_id(i)
#if defined WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            QCK(ng)%Vid(idPwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            QCK(ng)%Vid(idRwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            QCK(ng)%Vid(idUwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            QCK(ng)%Vid(idVwet)=var_id(i)
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            QCK(ng)%Vid(idpthR)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            QCK(ng)%Vid(idpthU)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            QCK(ng)%Vid(idpthV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            QCK(ng)%Vid(idpthW)=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            QCK(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            QCK(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            QCK(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            QCK(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            QCK(ng)%Vid(idv2dN)=var_id(i)
#ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            QCK(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            QCK(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsur))) THEN
            got_var(idUsur)=.TRUE.
            QCK(ng)%Vid(idUsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsur))) THEN
            got_var(idVsur)=.TRUE.
            QCK(ng)%Vid(idVsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsl))) THEN
            got_var(idUzsl)=.TRUE.
            QCK(ng)%Vid(idUzsl)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsl))) THEN
            got_var(idVzsl)=.TRUE.
            QCK(ng)%Vid(idVzsl)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            QCK(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            QCK(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsuE))) THEN
            got_var(idUsuE)=.TRUE.
            QCK(ng)%Vid(idUsuE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsuN))) THEN
            got_var(idVsuN)=.TRUE.
            QCK(ng)%Vid(idVsuN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsE))) THEN
            got_var(idUzsE)=.TRUE.
            QCK(ng)%Vid(idUzsE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsN))) THEN
            got_var(idVzsN)=.TRUE.
            QCK(ng)%Vid(idVzsN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPVzs))) THEN
            got_var(idPVzs)=.TRUE.
            QCK(ng)%Vid(idPVzs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRVzs))) THEN
            got_var(idRVzs)=.TRUE.
            QCK(ng)%Vid(idRVzs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            QCK(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            QCK(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            QCK(ng)%Vid(idDano)=var_id(i)
# ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            QCK(ng)%Vid(idHsbl)=var_id(i)
# endif
# ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            QCK(ng)%Vid(idHbbl)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            QCK(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            QCK(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            QCK(ng)%Vid(idSdif)=var_id(i)
# if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            QCK(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            QCK(ng)%Vid(idMtls)=var_id(i)
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            QCK(ng)%Vid(idPair)=var_id(i)
# endif
# if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPBLh))) THEN
            got_var(idPBLh)=.TRUE.
            QCK(ng)%Vid(idPBLh)=var_id(i)
# endif
# if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            QCK(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            QCK(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            QCK(ng)%Vid(idUaiE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            QCK(ng)%Vid(idVaiN)=var_id(i)
# endif
# if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            QCK(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            QCK(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            QCK(ng)%Vid(idLrad)=var_id(i)
# endif
# ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            QCK(ng)%Vid(idTair)=var_id(i)
#  ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            QCK(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            QCK(ng)%Vid(idrain)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            QCK(ng)%Vid(idEmPf)=var_id(i)
# ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            QCK(ng)%Vid(idSrad)=var_id(i)
# endif
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            QCK(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            QCK(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            QCK(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            QCK(ng)%Vid(idVbms)=var_id(i)
          END IF
#ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              QCK(ng)%Tid(itrc)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idsurT(itrc)))) THEN
              got_var(idsurT(itrc))=.TRUE.
              QCK(ng)%Vid(idsurT(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              QCK(ng)%Vid(idTsur(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idzslT(itrc)))) THEN
              got_var(idzslT(itrc))=.TRUE.
              QCK(ng)%Vid(idzslT(itrc))=var_id(i)
            END IF
          END DO
#endif
        END DO
!
!  Check if quicksave variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#if defined WET_DRY
        IF (.not.got_var(idPwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
        IF (.not.got_var(idpthR).and.Qout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Qout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Qout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Qout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idFsur).and.Qout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Qout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Qout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Qout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Qout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Qout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Qout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsur).and.Qout(idUsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsur).and.Qout(idVsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsl).and.Qout(idUzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsl).and.Qout(idVzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Qout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Qout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsuE).and.Qout(idUsuE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsuE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsuN).and.Qout(idVsuN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsuN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsE).and.Qout(idUzsE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsN).and.Qout(idVzsN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPVzs).and.Qout(idPVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRVzs).and.Qout(idRVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Qout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Qout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Qout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Qout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Qout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idVvis).and.Qout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Qout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SALINITY
        IF (.not.got_var(idSdif).and.Qout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Qout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Qout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Qout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
        IF (.not.got_var(idPBLh).and.Qout(idPBLh,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPBLh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Qout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Qout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Qout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Qout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined BULK_FLUXES || defined FRC_COUPLING
        IF (.not.got_var(idLhea).and.Qout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Qout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Qout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef BULK_FLUXES
        IF (.not.got_var(idTair).and.Qout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef EMINUSP
        IF (.not.got_var(idevap).and.Qout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Qout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idEmPf).and.Qout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Qout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
        IF (.not.got_var(idUsms).and.Qout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Qout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Qout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Qout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Qout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idsurT(itrc)).and.Qout(idsurT(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idsurT(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Qout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idzslT(itrc)).and.Qout(idzslT(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idzslT(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#endif

#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for quicksave variables. Get variable IDs.
!
        CALL bbl_def_nf90 (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined ICE_MODEL && defined SOLVE3D
!
!  Scan sea-ice variables from input NetCDF and activate switches for
!  quicksave variables. Get variable IDs.
!
        CALL ice_def_nf90 (ng, model, ldef, Qout, QCK)
#endif

#if defined SEDIMENT && defined SOLVE3D
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for quicksave variables. Get variable IDs.
!
        CALL sediment_def_nf90 (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined VEGETATION && defined SOLVE3D
!
!  Scan submerged aquatic vegetation model variables from input NetCDF
!  and activate switches for quicksave variables. Ger variable IDs.
!
        CALL vegetation_def_nf90 (ng, model, ldef, Hout, HIS)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined WEC_VF && defined SOLVE3D
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for quicksave variables. Get variable IDs.
!
        CALL wec_def_nf90 (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefQCK(ng).gt.0) THEN
          QCK(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefQCK(ng)*((ntstart(ng)-1)/ndefQCK(ng)))/   &
     &                   nQCK(ng)
        ELSE
          QCK(ng)%Rindex=(ntstart(ng)-1)/nQCK(ng)
        END IF
        QCK(ng)%Rindex=MIN(QCK(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_QUICK_NF90   - creating quicksave file,',t56,     &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_QUICK_NF90   - inquiring quicksave file,',t56,    &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_QUICK_NF90 - unable to create quicksave NetCDF',  &
     &        ' file:', 1x,a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_QUICK_NF90 - unable to open quicksave NetCDF',    &
     &        ' file: ',a)
  70  FORMAT (/,' DEF_QUICK_NF90 - unable to find variable: ',a,2x,     &
     &        ' in quicksave NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_quick_nf90

#if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_quick_pio (ng, model, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, j, ifield, itrc, nvd3, nvd4, varid
      integer :: recdim, status
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim
# endif
      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
#  if defined VEGETATION 
      integer :: sapgrd(4)
#  endif
      integer :: p3dgrd(4),  t3dgrd(4),  u3dgrd(4),  v3dgrd(4)
      integer :: p3dzgrd(4), t3dzgrd(4), u3dzgrd(4), v3dzgrd(4)
      integer :: w3dgrd(4)
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_quick_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=QCK(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new quicksave file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname), QCK(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xi_rho',    &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xi_u',      &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xi_v',      &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xi_psi',    &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'eta_rho',   &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'eta_u',     &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'eta_v',     &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'eta_psi',   &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'IorJ',      &
     &                 IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xy_rho',    &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xy_u',      &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xy_v',      &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xyz_rho',   &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xyz_u',     &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xyz_v',     &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xyz_w',     &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'N',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 's_rho',     &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 's_w',       &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        IF (Nslice.gt.0) THEN
          status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'z_slice', &
     &                   Nslice, DimIDs(34))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'tracer',    &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'NST',       &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nbed',      &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'xybed',     &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  if defined VEGETATION
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nveg',      &
     &                 NVEG, DimIDs(35))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nbands',    &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nphy',      &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nbac',      &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Ndom',      &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nfec',      &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'boundary',  &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nstate',    &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef MULTI_SCALE_B
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'Nscale',    &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, QCK(ng)%pioFile, ncname, 'obc_adjust',&
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, model, QCK(ng)%pioFile, ncname,              &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=4
# endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
#  endif
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
!
        t3dzgrd(1)=DimIDs( 1)
        t3dzgrd(2)=DimIDs( 5)
        t3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        t3dzgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
# ifdef SOLVE3D
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
!
        p3dzgrd(1)=DimIDs( 4)
        p3dzgrd(2)=DimIDs( 8)
        p3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        p3dzgrd(4)=DimIDs(12)
# endif
# if defined VEGETATION && defined SOLVE3D
!
!  Define dimensions for the submerged aquatic plant (SAP) variables.
!
        sapgrd(1)=DimIDs( 1)
        sapgrd(2)=DimIDs( 5)
        sapgrd(3)=DimIDs(35)
        sapgrd(4)=DimIDs(12)
# endif 
!
!  Define dimension vectors for staggered u-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
#  endif
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
!
        u3dzgrd(1)=DimIDs( 1)
        u3dzgrd(2)=DimIDs( 5)
        u3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        u3dzgrd(4)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
#  endif
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
!
        v3dzgrd(1)=DimIDs( 1)
        v3dzgrd(2)=DimIDs( 5)
        v3dzgrd(3)=DimIDs(34)          ! selected constant depth slices
        v3dzgrd(4)=DimIDs(12)
#  endif
# endif
# ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#  else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#  endif
#  ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#   if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#   else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#   endif
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        QCK(ng)%Rindex=0
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, QCK(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        QCK(ng)%pioVar(idtime)%dkind=PIO_TOUT
        QCK(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, model, QCK(ng)%pioFile,                      &
     &                 QCK(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef WET_DRY
!
!  Define wet/dry mask on PSI-points.
!
        Vinfo( 1)=Vname(1,idPwet)
        Vinfo( 2)=Vname(2,idPwet)
        Vinfo( 3)=Vname(3,idPwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idPwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(21)=Vname(6,idPwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idPwet,ng),r8)
        QCK(ng)%pioVar(idPwet)%dkind=PIO_FOUT
        QCK(ng)%pioVar(idPwet)%gtype=p2dvar
!
        status=def_var(ng, model, QCK(ng)%pioFile,                      &
     &                 QCK(ng)%pioVar(idPwet)%vd,                       &
     &                 PIO_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define wet/dry mask on RHO-points.
!
        Vinfo( 1)=Vname(1,idRwet)
        Vinfo( 2)=Vname(2,idRwet)
        Vinfo( 3)=Vname(3,idRwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idRwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(21)=Vname(6,idRwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idRwet,ng),r8)
        QCK(ng)%pioVar(idRwet)%dkind=PIO_FOUT
        QCK(ng)%pioVar(idRwet)%gtype=r2dvar
!
        status=def_var(ng, model, QCK(ng)%pioFile,                      &
     &                 QCK(ng)%pioVar(idRwet)%vd,                       &
     &                 PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define wet/dry mask on U-points.
!
        Vinfo( 1)=Vname(1,idUwet)
        Vinfo( 2)=Vname(2,idUwet)
        Vinfo( 3)=Vname(3,idUwet)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(14)=Vname(4,idUwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(21)=Vname(6,idUwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUwet,ng),r8)
        QCK(ng)%pioVar(idUwet)%dkind=PIO_FOUT
        QCK(ng)%pioVar(idUwet)%gtype=u2dvar
!
        status=def_var(ng, model, QCK(ng)%pioFile,                      &
     &                 QCK(ng)%pioVar(idUwet)%vd,                       &
     &                 PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define wet/dry mask on V-points.
!
        Vinfo( 1)=Vname(1,idVwet)
        Vinfo( 2)=Vname(2,idVwet)
        Vinfo( 3)=Vname(3,idVwet)
        Vinfo(14)=Vname(4,idVwet)
        Vinfo(16)=Vname(1,idtime)
        Vinfo( 9)='land'
        Vinfo(10)='water'
        Vinfo(21)=Vname(6,idVwet)
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVwet,ng),r8)
        QCK(ng)%pioVar(idVwet)%dkind=PIO_FOUT
        QCK(ng)%pioVar(idVwet)%gtype=v2dvar
!
        status=def_var(ng, model, QCK(ng)%pioFile,                      &
     &                 QCK(ng)%pioVar(idVwet)%vd,                       &
     &                 PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef SOLVE3D
!
!  Define time-varying depth of RHO-points.
!
        IF (Qout(idpthR,ng)) THEN
          Vinfo( 1)=Vname(1,idpthR)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthR))
          Vinfo( 3)=Vname(3,idpthR)
          Vinfo(14)=Vname(4,idpthR)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idpthR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthR,ng),r8)
          QCK(ng)%pioVar(idpthR)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idpthR)%gtype=r3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idpthR)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Qout(idpthU,ng)) THEN
          Vinfo( 1)=Vname(1,idpthU)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthU))
          Vinfo( 3)=Vname(3,idpthU)
          Vinfo(14)=Vname(4,idpthU)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idpthU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthU,ng),r8)
          QCK(ng)%pioVar(idpthU)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idpthU)%gtype=u3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idpthU)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Qout(idpthV,ng)) THEN
          Vinfo( 1)=Vname(1,idpthV)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthV))
          Vinfo( 3)=Vname(3,idpthV)
          Vinfo(14)=Vname(4,idpthV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idpthV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthV,ng),r8)
          QCK(ng)%pioVar(idpthV)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idpthV)%gtype=v3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idpthV)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Qout(idpthW,ng)) THEN
          Vinfo( 1)=Vname(1,idpthW)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthW))
          Vinfo( 3)=Vname(3,idpthW)
          Vinfo(14)=Vname(4,idpthW)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idpthW)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthW,ng),r8)
          QCK(ng)%pioVar(idpthW)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idpthW)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idpthW)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define free-surface.
!
        IF (Qout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          QCK(ng)%pioVar(idFsur)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idFsur)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idFsur)%vd,                     &
# ifdef WET_DRY
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
# else
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D U-momentum component.
!
        IF (Qout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          QCK(ng)%pioVar(idUbar)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUbar)%gtype=u2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUbar)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D V-momentum component.
!
        IF (Qout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          QCK(ng)%pioVar(idVbar)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVbar)%gtype=v2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVbar)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Qout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          QCK(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idu2dE)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idu2dE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Qout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          QCK(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idv2dN)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idv2dN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        IF (Qout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          QCK(ng)%pioVar(idUvel)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUvel)%gtype=u3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUvel)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component.
!
        IF (Qout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          QCK(ng)%pioVar(idVvel)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVvel)%gtype=v3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVvel)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface U-momentum component.
!
        IF (Qout(idUsur,ng)) THEN
          Vinfo( 1)=Vname(1,idUsur)
          Vinfo( 2)=Vname(2,idUsur)
          Vinfo( 3)=Vname(3,idUsur)
          Vinfo(14)=Vname(4,idUsur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsur,ng),r8)
          QCK(ng)%pioVar(idUsur)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUsur)%gtype=u2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUsur)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface V-momentum component.
!
        IF (Qout(idVsur,ng)) THEN
          Vinfo( 1)=Vname(1,idVsur)
          Vinfo( 2)=Vname(2,idVsur)
          Vinfo( 3)=Vname(3,idVsur)
          Vinfo(14)=Vname(4,idVsur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsur,ng),r8)
          QCK(ng)%pioVar(idVsur)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVsur)%gtype=v2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVsur)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D U-momentum component at specified constant depth slices.
!
        IF (Qout(idUzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsl)
          Vinfo( 2)=Vname(2,idUzsl)
          Vinfo( 3)=Vname(3,idUzsl)
          Vinfo(14)=Vname(4,idUzsl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsl,ng),r8)
          QCK(ng)%pioVar(idUzsl)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUzsl)%gtype=u3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUzsl)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D V-momentum component at specified constant depth slices.
!
        IF (Qout(idVzsl,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsl)
          Vinfo( 2)=Vname(2,idVzsl)
          Vinfo( 3)=Vname(3,idVzsl)
          Vinfo(14)=Vname(4,idVzsl)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVzsl)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idVzsl,ng),r8)
          QCK(ng)%pioVar(idVzsl)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVzsl)%gtype=v3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVzsl)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Qout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          QCK(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idu3dE)%gtype=r3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idu3dE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Qout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          QCK(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idv3dN)%gtype=r3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idv3dN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface Eastward momentum component at RHO-points.
!
        IF (Qout(idUsuE,ng)) THEN
          Vinfo( 1)=Vname(1,idUsuE)
          Vinfo( 2)=Vname(2,idUsuE)
          Vinfo( 3)=Vname(3,idUsuE)
          Vinfo(14)=Vname(4,idUsuE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUsuE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsuE,ng),r8)
          QCK(ng)%pioVar(idUsuE)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUsuE)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUsuE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model surface Northward momentum component at RHO-points.
!
        IF (Qout(idVsuN,ng)) THEN
          Vinfo( 1)=Vname(1,idVsuN)
          Vinfo( 2)=Vname(2,idVsuN)
          Vinfo( 3)=Vname(3,idVsuN)
          Vinfo(14)=Vname(4,idVsuN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVsuN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsuN,ng),r8)
          QCK(ng)%pioVar(idVsuN)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVsuN)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVsuN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points for selected
!  constant depth slices.
!
        IF (Qout(idUzsE,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idUzsE)
          Vinfo( 2)=Vname(2,idUzsE)
          Vinfo( 3)=Vname(3,idUzsE)
          Vinfo(14)=Vname(4,idUzsE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUzsE)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idUzsE,ng),r8)
          QCK(ng)%pioVar(idUzsE)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUzsE)%gtype=r3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUzsE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points for selected
!  constant depth slices.
!
        IF (Qout(idVzsN,ng).and.(Nslice.gt.0)) THEN
          Vinfo( 1)=Vname(1,idVzsN)
          Vinfo( 2)=Vname(2,idVzsN)
          Vinfo( 3)=Vname(3,idVzsN)
          Vinfo(14)=Vname(4,idVzsN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVzsN)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idVzsN,ng),r8)
          QCK(ng)%pioVar(idVzsN)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVzsN)%gtype=r3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVzsN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model potential vorticity at PSI-points for selected constant
!  depth slices.
!
        IF (Qout(idPVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idPVzs)
          Vinfo( 2)=Vname(2,idPVzs)
          Vinfo( 3)=Vname(3,idPVzs)
          Vinfo(14)=Vname(4,idPVzs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idPVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idPVzs,ng),r8)
          QCK(ng)%pioVar(idPVzs)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idPVzs)%gtype=p3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idPVzs)%vd,                     &
     &                   PIO_FOUT, nvd3, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define model relative vorticity at PSI-points for selected constant
!  depth slices.
!
        IF (Qout(idRVzs,ng)) THEN
          Vinfo( 1)=Vname(1,idRVzs)
          Vinfo( 2)=Vname(2,idRVzs)
          Vinfo( 3)=Vname(3,idRVzs)
          Vinfo(14)=Vname(4,idRVzs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idRVzs)
          Vinfo(22)='coordinates, z_slice'
          Aval(5)=REAL(Iinfo(1,idRVzs,ng),r8)
          QCK(ng)%pioVar(idRVzs)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idRVzs)%gtype=p3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idRVzs)%vd,                     &
     &                   PIO_FOUT, nvd3, p3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the S-direction.
!
        IF (Qout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          QCK(ng)%pioVar(idWvel)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idWvel)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idWvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Qout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          QCK(ng)%pioVar(idOvel)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idOvel)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idOvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Qout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            QCK(ng)%pioTrc(itrc)%dkind=PIO_FOUT
            QCK(ng)%pioTrc(itrc)%gtype=r3dvar
!
            status=def_var(ng, model, QCK(ng)%pioFile,                  &
     &                     QCK(ng)%pioTrc(itrc)%vd,                     &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define surface tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Qout(idsurT(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idsurT(itrc))
            Vinfo( 2)=Vname(2,idsurT(itrc))
            Vinfo( 3)=Vname(3,idsurT(itrc))
            Vinfo(14)=Vname(4,idsurT(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idsurT(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            QCK(ng)%pioVar(idsurT(itrc))%dkind=PIO_FOUT
            QCK(ng)%pioVar(idsurT(itrc))%gtype=r2dvar
!
            status=def_var(ng, model, QCK(ng)%pioFile,                  &
     &                     QCK(ng)%pioVar(idsurT(itrc))%vd,             &
     &                     PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define surface tracer type variables at specified depth slices.
!
        DO itrc=1,NT(ng)
          IF (Qout(idzslT(itrc),ng).and.(Nslice.gt.0)) THEN
            Vinfo( 1)=Vname(1,idzslT(itrc))
            Vinfo( 2)=Vname(2,idzslT(itrc))
            Vinfo( 3)=Vname(3,idzslT(itrc))
            Vinfo(14)=Vname(4,idzslT(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idzslT(itrc))
            Vinfo(22)='coordinates, z_slice'
            Aval(5)=REAL(r3dvar,r8)
            QCK(ng)%pioVar(idzslT(itrc))%dkind=PIO_FOUT
            QCK(ng)%pioVar(idzslT(itrc))%gtype=r2dvar
!
            status=def_var(ng, model, QCK(ng)%pioFile,                  &
     &                     QCK(ng)%pioVar(idzslT(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, t3dzgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Qout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          QCK(ng)%pioVar(idDano)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idDano)%gtype=r3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idDano)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Qout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idHsbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          QCK(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idHsbl)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idHsbl)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Qout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idHbbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          QCK(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idHbbl)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idHbbl)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define vertical viscosity coefficient.
!
        IF (Qout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          QCK(ng)%pioVar(idVvis)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVvis)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVvis)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Qout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          QCK(ng)%pioVar(idTdif)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idTdif)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idTdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Qout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          QCK(ng)%pioVar(idSdif)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idSdif)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idSdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Qout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idMtke)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          QCK(ng)%pioVar(idMtke)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idMtke)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idMtke)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Qout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idMtls)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          QCK(ng)%pioVar(idMtls)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idMtls)%gtype=w3dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idMtls)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Qout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          QCK(ng)%pioVar(idPair)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idPair)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idPair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
!
!  Define atmosphere boundary layer thickness.
!
        IF (Qout(idPBLh,ng)) THEN
          Vinfo( 1)=Vname(1,idPBLh)
          Vinfo( 2)=Vname(2,idPBLh)
          Vinfo( 3)=Vname(3,idPBLh)
          Vinfo(14)=Vname(4,idPBLh)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idPBLh)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPBLh,ng),r8)
          QCK(ng)%pioVar(idPBLh)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idPBLh)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idPBLh)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Qout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          QCK(ng)%pioVar(idUair)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUair)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Qout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          QCK(ng)%pioVar(idVair)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVair)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface winds at RHO-points.
!
        IF (Qout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          Vinfo( 2)=Vname(2,idUaiE)
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          QCK(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUaiE)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUaiE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Qout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          Vinfo( 2)=Vname(2,idVaiN)
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          QCK(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVaiN)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVaiN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Qout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            Vinfo( 2)=Vname(2,idTsur(itrc))
            Vinfo( 3)=Vname(3,idTsur(itrc))
            IF (itrc.eq.itemp) THEN
              Vinfo(11)='upward flux, cooling'
              Vinfo(12)='downward flux, heating'
            ELSE IF (itrc.eq.isalt) THEN
              Vinfo(11)='upward flux, freshening (net precipitation)'
              Vinfo(12)='downward flux, salting (net evaporation)'
            END IF
            Vinfo(14)=Vname(4,idTsur(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            QCK(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
            QCK(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
!
            status=def_var(ng, model, QCK(ng)%pioFile,                  &
     &                     QCK(ng)%pioVar(idTsur(itrc))%vd,             &
     &                     PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#  if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Define latent heat flux.
!
        IF (Qout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          QCK(ng)%pioVar(idLhea)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idLhea)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idLhea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Qout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          QCK(ng)%pioVar(idShea)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idShea)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idShea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Qout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          QCK(ng)%pioVar(idLrad)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idLrad)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idLrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined BULK_FLUXES
!
!  Define atmospheric air temperature.
!
        IF (Qout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          QCK(ng)%pioVar(idTair)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idTair)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idTair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Qout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          QCK(ng)%pioVar(idevap)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idevap)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idevap)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Qout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          QCK(ng)%pioVar(idrain)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idrain)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idrain)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
!
!  Define E-P flux.
!
        IF (Qout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idEmPf)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          QCK(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idEmPf)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idEmPf)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SHORTWAVE
!
!  Define net shortwave radiation flux.
!
        IF (Qout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          QCK(ng)%pioVar(idSrad)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idSrad)%gtype=r2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idSrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
!
!  Define surface U-momentum stress.
!
        IF (Qout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          QCK(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUsms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Qout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          QCK(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVsms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Qout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          QCK(ng)%pioVar(idUbms)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idUbms)%gtype=u2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idUbms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Qout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          QCK(ng)%pioVar(idVbms)%dkind=PIO_FOUT
          QCK(ng)%pioVar(idVbms)%gtype=v2dvar
!
          status=def_var(ng, model, QCK(ng)%pioFile,                    &
     &                   QCK(ng)%pioVar(idVbms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_pio (ng, model, ldef, Qout, QCK,                   &
     &                    t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_pio (ng, model, ldef, Qout, QCK,                   &
     &                    t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_pio (ng, model, ldef, Qout, QCK,              &
     &                         t2dgrd, u2dgrd, v2dgrd,                  &
     &                         b3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined VEGETATION && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the submerged aquatic vegetation model variables.
!-----------------------------------------------------------------------
!
        CALL vegetation_def_pio (ng, model, ldef, Qout, QCK,            &
     &                           sapgrd, t2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
        CALL wec_def_pio (ng, model, ldef, Qout, QCK,                   &
     &                    t2dgrd, u2dgrd, v2dgrd,                       &
     &                    t3dgrd, u3dgrd, v3dgrd, w3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, QCK(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, QCK(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing quicksave file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=QCK(ng)%name
!
!  Open quicksave file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1, QCK(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,60) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = QCK(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = QCK(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  quicksave variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            QCK(ng)%pioVar(idtime)%vd=var_desc(i)
            QCK(ng)%pioVar(idtime)%dkind=PIO_TOUT
            QCK(ng)%pioVar(idtime)%gtype=0
# if defined WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            QCK(ng)%pioVar(idPwet)%vd=var_desc(i)
            QCK(ng)%pioVar(idPwet)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idPwet)%gtype=p2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            QCK(ng)%pioVar(idRwet)%vd=var_desc(i)
            QCK(ng)%pioVar(idRwet)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idRwet)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            QCK(ng)%pioVar(idUwet)%vd=var_desc(i)
            QCK(ng)%pioVar(idUwet)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUwet)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            QCK(ng)%pioVar(idVwet)%vd=var_desc(i)
            QCK(ng)%pioVar(idVwet)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVwet)%gtype=v2dvar
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            QCK(ng)%pioVar(idpthR)%vd=var_desc(i)
            QCK(ng)%pioVar(idpthR)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idpthR)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            QCK(ng)%pioVar(idpthU)%vd=var_desc(i)
            QCK(ng)%pioVar(idpthU)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idpthU)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            QCK(ng)%pioVar(idpthV)%vd=var_desc(i)
            QCK(ng)%pioVar(idpthV)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idpthV)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            QCK(ng)%pioVar(idpthW)%vd=var_desc(i)
            QCK(ng)%pioVar(idpthW)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idpthW)%gtype=w3dvar
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            QCK(ng)%pioVar(idFsur)%vd=var_desc(i)
            QCK(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            QCK(ng)%pioVar(idUbar)%vd=var_desc(i)
            QCK(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            QCK(ng)%pioVar(idVbar)%vd=var_desc(i)
            QCK(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVbar)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            QCK(ng)%pioVar(idu2dE)%vd=var_desc(i)
            QCK(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idu2dE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            QCK(ng)%pioVar(idv2dN)%vd=var_desc(i)
            QCK(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idv2dN)%gtype=r2dvar
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            QCK(ng)%pioVar(idUvel)%vd=var_desc(i)
            QCK(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            QCK(ng)%pioVar(idVvel)%vd=var_desc(i)
            QCK(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVvel)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsur))) THEN
            got_var(idUsur)=.TRUE.
            QCK(ng)%pioVar(idUsur)%vd=var_desc(i)
            QCK(ng)%pioVar(idUsur)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUsur)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsur))) THEN
            got_var(idVsur)=.TRUE.
            QCK(ng)%pioVar(idVsur)%vd=var_desc(i)
            QCK(ng)%pioVar(idVsur)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVsur)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsl))) THEN
            got_var(idUzsl)=.TRUE.
            QCK(ng)%pioVar(idUzsl)%vd=var_desc(i)
            QCK(ng)%pioVar(idUzsl)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUzsl)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsl))) THEN
            got_var(idVzsl)=.TRUE.
            QCK(ng)%pioVar(idVzsl)%vd=var_desc(i)
            QCK(ng)%pioVar(idVzsl)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVzsl)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            QCK(ng)%pioVar(idu3dE)%vd=var_desc(i)
            QCK(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idu3dE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            QCK(ng)%pioVar(idv3dN)%vd=var_desc(i)
            QCK(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idv3dN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsuE))) THEN
            got_var(idUsuE)=.TRUE.
            QCK(ng)%pioVar(idUsuE)%vd=var_desc(i)
            QCK(ng)%pioVar(idUsuE)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUsuE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsuN))) THEN
            got_var(idVsuN)=.TRUE.
            QCK(ng)%pioVar(idVsuN)%vd=var_desc(i)
            QCK(ng)%pioVar(idVsuN)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVsuN)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUzsE))) THEN
            got_var(idUzsE)=.TRUE.
            QCK(ng)%pioVar(idUzsE)%vd=var_desc(i)
            QCK(ng)%pioVar(idUzsE)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUzsE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVzsN))) THEN
            got_var(idVzsN)=.TRUE.
            QCK(ng)%pioVar(idVzsN)%vd=var_desc(i)
            QCK(ng)%pioVar(idVzsN)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVzsN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPVzs))) THEN
            got_var(idPVzs)=.TRUE.
            QCK(ng)%pioVar(idPVzs)%vd=var_desc(i)
            QCK(ng)%pioVar(idPVzs)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idPVzs)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRVzs))) THEN
            got_var(idRVzs)=.TRUE.
            QCK(ng)%pioVar(idRVzs)%vd=var_desc(i)
            QCK(ng)%pioVar(idRVzs)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idRVzs)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            QCK(ng)%pioVar(idWvel)%vd=var_desc(i)
            QCK(ng)%pioVar(idWvel)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idWvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            QCK(ng)%pioVar(idOvel)%vd=var_desc(i)
            QCK(ng)%pioVar(idOvel)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idOvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            QCK(ng)%pioVar(idDano)%vd=var_desc(i)
            QCK(ng)%pioVar(idDano)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idDano)%gtype=r3dvar
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            QCK(ng)%pioVar(idHsbl)%vd=var_desc(i)
            QCK(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idHsbl)%gtype=r2dvar
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            QCK(ng)%pioVar(idHbbl)%vd=var_desc(i)
            QCK(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idHbbl)%gtype=r2dvar
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            QCK(ng)%pioVar(idVvis)%vd=var_desc(i)
            QCK(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVvis)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            QCK(ng)%pioVar(idTdif)%vd=var_desc(i)
            QCK(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idTdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            QCK(ng)%pioVar(idSdif)%vd=var_desc(i)
            QCK(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idSdif)%gtype=w3dvar
#  if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            QCK(ng)%pioVar(idMtke)%vd=var_desc(i)
            QCK(ng)%pioVar(idMtke)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idMtke)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            QCK(ng)%pioVar(idMtls)%vd=var_desc(i)
            QCK(ng)%pioVar(idMtls)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idMtls)%gtype=w3dvar
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            QCK(ng)%pioVar(idPair)%vd=var_desc(i)
            QCK(ng)%pioVar(idPair)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idPair)%gtype=r2dvar
#  endif
#  if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
            got_var(idPBLh)=.TRUE.
            QCK(ng)%pioVar(idPBLh)%vd=var_desc(i)
            QCK(ng)%pioVar(idPBLh)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idPBLh)%gtype=r2dvar
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            QCK(ng)%pioVar(idUair)%vd=var_desc(i)
            QCK(ng)%pioVar(idUair)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            QCK(ng)%pioVar(idVair)%vd=var_desc(i)
            QCK(ng)%pioVar(idVair)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUair)=.TRUE.
            QCK(ng)%pioVar(idUaiE)%vd=var_desc(i)
            QCK(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUaiE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVair)=.TRUE.
            QCK(ng)%pioVar(idVaiN)%vd=var_desc(i)
            QCK(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVaiN)%gtype=r2dvar
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            QCK(ng)%pioVar(idLhea)%vd=var_desc(i)
            QCK(ng)%pioVar(idLhea)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idLhea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            QCK(ng)%pioVar(idShea)%vd=var_desc(i)
            QCK(ng)%pioVar(idShea)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idShea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            QCK(ng)%pioVar(idLrad)%vd=var_desc(i)
            QCK(ng)%pioVar(idLrad)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idLrad)%gtype=r2dvar
#  endif
#  ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            QCK(ng)%pioVar(idTair)%vd=var_desc(i)
            QCK(ng)%pioVar(idTair)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idTair)%gtype=r2dvar
#   ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            QCK(ng)%pioVar(idevap)%vd=var_desc(i)
            QCK(ng)%pioVar(idevap)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idevap)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            QCK(ng)%pioVar(idrain)%vd=var_desc(i)
            QCK(ng)%pioVar(idrain)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idrain)%gtype=r2dvar
#   endif
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            QCK(ng)%pioVar(idEmPf)%vd=var_desc(i)
            QCK(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idEmPf)%gtype=r2dvar
#  ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            QCK(ng)%pioVar(idSrad)%vd=var_desc(i)
            QCK(ng)%pioVar(idSrad)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idSrad)%gtype=r2dvar
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            QCK(ng)%pioVar(idUsms)%vd=var_desc(i)
            QCK(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            QCK(ng)%pioVar(idVsms)%vd=var_desc(i)
            QCK(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVsms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            QCK(ng)%pioVar(idUbms)%vd=var_desc(i)
            QCK(ng)%pioVar(idUbms)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idUbms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            QCK(ng)%pioVar(idVbms)%vd=var_desc(i)
            QCK(ng)%pioVar(idVbms)%dkind=PIO_FOUT
            QCK(ng)%pioVar(idVbms)%gtype=v2dvar
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              QCK(ng)%pioTrc(itrc)%vd=var_desc(i)
              QCK(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              QCK(ng)%pioTrc(itrc)%gtype=r3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idsurT(itrc)))) THEN
              got_var(idsurT(itrc))=.TRUE.
              QCK(ng)%pioVar(idsurT(itrc))%vd=var_desc(i)
              QCK(ng)%pioVar(idsurT(itrc))%dkind=PIO_FOUT
              QCK(ng)%pioVar(idsurT(itrc))%gtype=r2dvar
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              QCK(ng)%pioVar(idTsur(itrc))%vd=var_desc(i)
              QCK(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              QCK(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idzslT(itrc)))) THEN
              got_var(idzslT(itrc))=.TRUE.
              QCK(ng)%pioVar(idzslT(itrc))%vd=var_desc(i)
              QCK(ng)%pioVar(idzslT(itrc))%dkind=PIO_FOUT
              QCK(ng)%pioVar(idzslT(itrc))%gtype=r3dvar
            END IF
          END DO
# endif
        END DO
!
!  Check if quicksave variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined WET_DRY
        IF (.not.got_var(idPwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwet)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVwet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idpthR).and.Qout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Qout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Qout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Qout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idFsur).and.Qout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Qout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Qout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Qout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Qout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Qout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Qout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsur).and.Qout(idUsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsur).and.Qout(idVsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsl).and.Qout(idUzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsl).and.Qout(idVzsl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Qout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Qout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUsuE).and.Qout(idUsuE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsuE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsuN).and.Qout(idVsuN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsuN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUzsE).and.Qout(idUzsE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUzsE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVzsN).and.Qout(idVzsN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVzsN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPVzs).and.Qout(idPVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRVzs).and.Qout(idRVzs,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRVzs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Qout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Qout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Qout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Qout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Qout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idVvis).and.Qout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Qout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SALINITY
        IF (.not.got_var(idSdif).and.Qout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Qout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Qout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Qout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES && defined FRC_COUPLING && defined PBLH
        IF (.not.got_var(idPBLh).and.Qout(idPBLh,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPBLh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Qout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Qout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Qout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Qout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
        IF (.not.got_var(idLhea).and.Qout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Qout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Qout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (.not.got_var(idTair).and.Qout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef EMINUSP
        IF (.not.got_var(idevap).and.Qout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Qout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
        IF (.not.got_var(idEmPf).and.Qout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Qout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Qout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Qout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Qout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Qout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Qout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idsurT(itrc)).and.Qout(idsurT(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idsurT(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Qout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idzslT(itrc)).and.Qout(idzslT(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idzslT(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for quicksave variables. Get variable
!  IDs.
!
        CALL bbl_def_pio (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!  Scan sea-ice variables from input NetCDF and activate switches for
!  quicksave variables. Get variable IDs.
!
        CALL ice_def_pio (ng, model, ldef, Qout, QCK)
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for quicksave variables. Get variable IDs.
!
        CALL sediment_def_pio (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined VEGETATION && defined SOLVE3D
!
!  Scan submerged aquatic vegetation model variables from input NetCDF
!  and activate switches for quicksave variables. Ger variable IDs.
!
        CALL vegetation_def_pio (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for quicksave variables. Get variable IDs.
!
        CALL wec_def_pio (ng, model, ldef, Qout, QCK)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefQCK(ng).gt.0) THEN
          QCK(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefQCK(ng)*((ntstart(ng)-1)/ndefQCK(ng)))/   &
     &                   nQCK(ng)
        ELSE
          QCK(ng)%Rindex=(ntstart(ng)-1)/nQCK(ng)
        END IF
        QCK(ng)%Rindex=MIN(QCK(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_QUICK_PIO    - creating quicksave file,',t56,     &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_QUICK_PIO    - inquiring quicksave file,',t56,    &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_QUICK_PIO - unable to create quicksave NetCDF',   &
     &        ' file:', 1x,a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_QUICK_PIO - unable to open quicksave NetCDF',     &
     &        ' file: ',a)
  70  FORMAT (/,' DEF_QUICK_PIO - unable to find variable: ',a,2x,      &
     &        ' in quicksave NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_quick_pio
#endif

      END MODULE def_quick_mod
