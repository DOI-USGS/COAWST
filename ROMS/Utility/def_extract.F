#include "cppdefs.h"
      MODULE def_extract_mod

#ifdef GRID_EXTRACT
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module creates output extract file using either the standard   !
!  NetCDF library or the Parallel-IO (PIO) library.  It defines its    !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!  If ExtractFlag > 1:   Decimation                                    !
!                                                                      !
!  It decimates the field solution at the prescribed integer factor,   !
!  ExtractFlag . For example, if ExtractFlag=2 (recommended), the      !
!  output fieldsare written at every other point, resulting in coarser !
!  data resolution. This strategy is advantageous in mixed resolution, !
!  split 4D-Var applications where the outer loop background (prior)   !
!  trajectory may be computed at a higher resolution than in the inner !
!  loop minimization to accelerate the calculations. For decimation to !
!  work, the number of parent grid RHO-points (0: Lm+1, 0:Mm+1) must   !
!  be multiples of dec_fator. That is,                                 !
!                                                                      !
!     MOD(Lm+1, ExtractFlag) = 0                                       !
!     MOD(Mm+1, ExtractFlag) = 0                                       !
!                                                                      !
!  If ExtractFlag < 2:   Interpolation                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef BIOLOGY
      USE mod_biology
# endif
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
!
# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
      USE bbl_output_mod,      ONLY : bbl_def_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE bbl_output_mod,      ONLY : bbl_def_pio
#  endif
# endif
      USE def_dim_mod,         ONLY : def_dim
      USE def_info_mod,        ONLY : def_info
      USE def_var_mod,         ONLY : def_var
# if defined ICE_MODEL && defined SOLVE3D
      USE ice_output_mod,      ONLY : ice_def_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE ice_output_mod,      ONLY : ice_def_pio
#  endif
# endif
# if defined SEDIMENT && defined SOLVE3D
      USE sediment_output_mod, ONLY : sediment_def_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE sediment_output_mod, ONLY : sediment_def_pio
#  endif
# endif
      USE strings_mod,         ONLY : FoundError
# if defined WEC_VF && defined SOLVE3D
      USE wec_output_mod,      ONLY : wec_def_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      USE wec_output_mod,      ONLY : wec_def_pio
#  endif
# endif
      USE wrt_info_mod,        ONLY : wrt_info
!
      implicit none
!
      PUBLIC  :: def_extract
      PRIVATE :: def_extract_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_extract_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_extract (ng, ldef)
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
!  Create a new decimate file according to IO type.
!-----------------------------------------------------------------------
!
      SELECT CASE (XTR(ng)%IOtype)
        CASE (io_nf90)
          CALL def_extract_nf90 (ng, iNLM, ldef)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_extract_pio (ng, iNLM, ldef)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) XTR(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_EXTRACT - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_extract
!
!***********************************************************************
      SUBROUTINE def_extract_nf90 (ng, model, ldef)
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
# ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
# endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
# ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
# endif

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: p2grd(3), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#  ifdef ADJUST_BOUNDARY
      integer :: t3dobc(5)
#  endif
# endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_extract_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=XTR(ng)%name
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
!  Create a new decimate file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), XTR(ng)%ncid)
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
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xi_rho',       &
     &                 xtr_IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xi_u',         &
     &                 xtr_IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xi_v',         &
     &                 xtr_IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xi_psi',       &
     &                 xtr_IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'eta_rho',      &
     &                 xtr_IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'eta_u',        &
     &                 xtr_IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'eta_v',        &
     &                 xtr_IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'eta_psi',      &
     &                 xtr_IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'IorJ',         &
     &                 xtr_IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xy_rho',       &
     &                 xtr_IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xy_u',         &
     &                 xtr_IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xy_v',         &
     &                 xtr_IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xyz_rho',      &
     &                 xtr_IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xyz_u',        &
     &                 xtr_IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xyz_v',        &
     &                 xtr_IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xyz_w',        &
     &                 xtr_IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'N',            &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'NST',          &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'xybed',        &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nbands',       &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef MULTI_SCALE_B
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'Nscale',       &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, XTR(ng)%ncid, ncname, 'obc_adjust',   &
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, model, XTR(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
# ifdef ADJUST_BOUNDARY
        brecdim=DimIDs(31)
# endif
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
#  endif
# endif
# ifdef ADJUST_BOUNDARY
        t2dobc(1)=IorJdim
        t2dobc(2)=DimIDs(14)
        t2dobc(3)=brecdim
        t2dobc(4)=DimIDs(12)
#  ifdef SOLVE3D
        t3dobc(1)=IorJdim
        t3dobc(2)=DimIDs( 9)
        t3dobc(3)=DimIDs(14)
        t3dobc(4)=brecdim
        t3dobc(5)=DimIDs(12)
#  endif
# endif
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        p2grd(1)=DimIDs( 4)
        p2grd(2)=DimIDs( 8)
        p2grd(3)=DimIDs(12)
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
        XTR(ng)%Rindex=0
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
        CALL def_info (ng, model, XTR(ng)%ncid, ncname, DimIDs)
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
        status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idtime),    &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef SO_SEMI
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='real Ritz eigenvalues'
        status=def_var(ng, model, XTR(ng)%ncid, varid, NF_TYPE,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ng, model, XTR(ng)%ncid, varid, NF_TYPE,         &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

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
        status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idPwet),    &
     &                 NF_FOUT, nvd3, p2grd, Aval, Vinfo, ncname,       &
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
        status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRwet),    &
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
        status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUwet),    &
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
        status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVwet),    &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname,      &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
# ifdef SOLVE3D
!
!  Define time-varying depth of RHO-points.
!
        IF (Hout(idpthR,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idpthR),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Hout(idpthU,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idpthU),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Hout(idpthV,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idpthV),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Hout(idpthW,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idpthW),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idFsur),  &
# ifdef WET_DRY
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
# else
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined FORWARD_WRITE && defined FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRzet)
          Vinfo( 2)=Vname(2,idRzet)
          Vinfo( 3)=Vname(3,idRzet)
          Vinfo(14)=Vname(4,idRzet)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idRzet)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRzet,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRzet),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
        END IF

# ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idFsur)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(ifield),  &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D U-momentum component.
!
        IF (Hout(idUbar,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUbar),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRu2d)
          Vinfo( 2)=Vname(2,idRu2d)
          Vinfo( 3)=Vname(3,idRu2d)
          Vinfo(14)=Vname(4,idRu2d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idRu2d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRu2d,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRu2d),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRuct)
          Vinfo( 2)=Vname(2,idRuct)
          Vinfo( 3)=Vname(3,idRuct)
          Vinfo(14)=Vname(4,idRuct)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idRuct)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRuct,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRuct),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
!
          Vinfo( 1)=Vname(1,idUfx1)
          Vinfo( 2)=Vname(2,idUfx1)
          Vinfo( 3)=Vname(3,idUfx1)
          Vinfo(14)=Vname(4,idUfx1)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUfx1)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUfx1,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUfx1),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idUfx2)
          Vinfo( 2)=Vname(2,idUfx2)
          Vinfo( 3)=Vname(3,idUfx2)
          Vinfo(14)=Vname(4,idUfx2)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUfx2)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUfx2,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUfx2),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
        END IF

# ifdef ADJUST_BOUNDARY
!
!  Define 2D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbar)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(ifield),  &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D V-momentum component.
!
        IF (Hout(idVbar,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVbar),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRv2d)
          Vinfo( 2)=Vname(2,idRv2d)
          Vinfo( 3)=Vname(3,idRv2d)
          Vinfo(14)=Vname(4,idRv2d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idRv2d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRv2d,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRv2d),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRvct)
          Vinfo( 2)=Vname(2,idRvct)
          Vinfo( 3)=Vname(3,idRvct)
          Vinfo(14)=Vname(4,idRvct)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idRvct)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRvct,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRvct),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
!
          Vinfo( 1)=Vname(1,idVfx1)
          Vinfo( 2)=Vname(2,idVfx1)
          Vinfo( 3)=Vname(3,idVfx1)
          Vinfo(14)=Vname(4,idVfx1)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVfx1)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVfx1,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVfx1),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idVfx2)
          Vinfo( 2)=Vname(2,idVfx2)
          Vinfo( 3)=Vname(3,idVfx2)
          Vinfo(14)=Vname(4,idVfx2)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVfx2)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVfx2,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVfx2),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif
        END IF

# ifdef ADJUST_BOUNDARY
!
!  Define 2D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbar)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(ifield),  &
     &                   NF_FOUT, 4, t2dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Hout(idu2dE,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idu2dE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Hout(idv2dN,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idv2dN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        IF (Hout(idUvel,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUvel),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRu3d)
          Vinfo( 2)=Vname(2,idRu3d)
          Vinfo( 3)=Vname(3,idRu3d)
          Vinfo(14)=Vname(4,idRu3d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idRu3d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRu3d,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRu3d),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
        END IF

#  ifdef ADJUST_BOUNDARY
!
!  Define 3D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUvel)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(ifield),  &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 3D V-momentum component.
!
        IF (Hout(idVvel,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVvel),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRv3d)
          Vinfo( 2)=Vname(2,idRv3d)
          Vinfo( 3)=Vname(3,idRv3d)
          Vinfo(14)=Vname(4,idRv3d)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idRv3d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRv3d,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idRv3d),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
        END IF

#  ifdef ADJUST_BOUNDARY
!
!  Define 3D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvel)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(ifield),  &
     &                   NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,       &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Hout(idu3dE,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idu3dE),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Hout(idv3dN,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idv3dN),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the Z-direction.
!
        IF (Hout(idWvel,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idWvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idOvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef OMEGA_IMPLICIT
!
!  Define S-coordinate implicit vertical "omega" momentum component.
!
        IF (Hout(idOvil,ng)) THEN
          Vinfo( 1)=Vname(1,idOvil)
          Vinfo( 2)=Vname(2,idOvil)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvil)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idOvil)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvil,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idOvil),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
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
            Aval(5)=REAL(Iinfo(1,idTvar(itrc),ng),r8)
            status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Tid(itrc),  &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#  ifdef ADJUST_BOUNDARY
!
!  Define tracer type variables open boundaries.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),60) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#   endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(ifield),&
     &                     NF_FOUT, 5, t3dobc, Aval, Vinfo, ncname,     &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idDano),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idHsbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idHbbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Define out KPP nonlocal transport.
!
        DO itrc=1,NAT
          IF (Hout(idGhat(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idGhat(itrc))
            Vinfo( 2)=Vname(2,idGhat(itrc))
            Vinfo( 3)=Vname(3,idGhat(itrc))
            Vinfo(14)=Vname(4,idGhat(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idGhat(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idGhat(itrc),ng),r8)
            status=def_var(ng, model, XTR(ng)%ncid,                     &
     &                     XTR(ng)%Vid(idGhat(itrc)), NF_FOUT,          &
     &                     nvd4, w3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVvis),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idTdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idSdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idMtke),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef FORWARD_WRITE
!
          Vinfo( 1)=Vname(1,idVmKK)
          Vinfo( 2)=Vname(2,idVmKK)
          Vinfo( 3)=Vname(3,idVmKK)
          Vinfo(14)=Vname(4,idVmKK)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idVmKK)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVmKK),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idMtls),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef FORWARD_WRITE
!
          Vinfo( 1)=Vname(1,idVmLS)
          Vinfo( 2)=Vname(2,idVmLS)
          Vinfo( 3)=Vname(3,idVmLS)
          Vinfo(14)=Vname(4,idVmLS)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idVmLS)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVmLS),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   if defined FORWARD_WRITE && defined GLS_MIXING
!
          Vinfo( 1)=Vname(1,idVmKP)
          Vinfo( 2)=Vname(2,idVmKP)
          Vinfo( 3)=Vname(3,idVmKP)
          Vinfo(14)=Vname(4,idVmKP)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idVmKP)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKP,ng),r8)
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVmKP),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Hout(idPair,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idPair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Hout(idUair,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Hout(idVair,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface winds at RHO-points.
!
        IF (Hout(idUaiE,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUaiE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Hout(idVaiN,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVaiN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
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
            status=def_var(ng, model, XTR(ng)%ncid,                     &
     &                     XTR(ng)%Vid(idTsur(itrc)), NF_FOUT,          &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#  if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined BULK_FLUXES
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idTair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idevap),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idrain),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
!
!  Define E-P flux.
!
        IF (Hout(idEmPf,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idEmPf),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SHORTWAVE
!
!  Define net shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idSrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUsms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVsms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idUbms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
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
          status=def_var(ng, model, XTR(ng)%ncid, XTR(ng)%Vid(idVbms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_nf90 (ng, model, ldef, Hout, XTR,                  &
     &                     t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_nf90 (ng, model, ldef, Hout, XTR,                  &
     &                     t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_nf90 (ng, model, ldef, Hout, XTR,             &
     &                          t2dgrd, u2dgrd, v2dgrd,                 &
     &                          b3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
        CALL wec_def_nf90 (ng, model, ldef, Hout, XTR,                  &
     &                     t2dgrd, u2dgrd, v2dgrd,                      &
     &                     t3dgrd, u3dgrd, v3dgrd, w3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, XTR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, XTR(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing decimate file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=XTR(ng)%name
!
!  Open decimate file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, XTR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,60) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = XTR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = XTR(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  decimate variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            XTR(ng)%Vid(idtime)=var_id(i)
# if defined WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            XTR(ng)%Vid(idPwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            XTR(ng)%Vid(idRwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            XTR(ng)%Vid(idUwet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            XTR(ng)%Vid(idVwet)=var_id(i)
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            XTR(ng)%Vid(idpthR)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            XTR(ng)%Vid(idpthU)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            XTR(ng)%Vid(idpthV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            XTR(ng)%Vid(idpthW)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            XTR(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            XTR(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            XTR(ng)%Vid(idVbar)=var_id(i)
# ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            XTR(ng)%Vid(idSbry(isFsur))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            XTR(ng)%Vid(idSbry(isUbar))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            XTR(ng)%Vid(idSbry(isVbar))=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            XTR(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            XTR(ng)%Vid(idv2dN)=var_id(i)
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            XTR(ng)%Vid(idRzet)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            XTR(ng)%Vid(idRu2d)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            XTR(ng)%Vid(idRv2d)=var_id(i)
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRuct))) THEN
            got_var(idRuct)=.TRUE.
            XTR(ng)%Vid(idRuct)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRvct))) THEN
            got_var(idRvct)=.TRUE.
            XTR(ng)%Vid(idRvct)=var_id(i)
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            XTR(ng)%Vid(idUfx1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            XTR(ng)%Vid(idUfx2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            XTR(ng)%Vid(idVfx1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            XTR(ng)%Vid(idVfx2)=var_id(i)
#   ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            XTR(ng)%Vid(idRu3d)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            XTR(ng)%Vid(idRv3d)=var_id(i)
#   endif
#  endif
# endif
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            XTR(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            XTR(ng)%Vid(idVvel)=var_id(i)
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            XTR(ng)%Vid(idSbry(isUvel))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            XTR(ng)%Vid(idSbry(isVvel))=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            XTR(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            XTR(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            XTR(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            XTR(ng)%Vid(idOvel)=var_id(i)
#  ifdef OMEGA_IMPLICIT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvil))) THEN
            got_var(idOvil)=.TRUE.
            XTR(ng)%Vid(idOvil)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            XTR(ng)%Vid(idDano)=var_id(i)
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            XTR(ng)%Vid(idHsbl)=var_id(i)
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            XTR(ng)%Vid(idHbbl)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            XTR(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            XTR(ng)%Vid(idTdif)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            XTR(ng)%Vid(idSdif)=var_id(i)
#  if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            XTR(ng)%Vid(idMtke)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            XTR(ng)%Vid(idMtls)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            XTR(ng)%Vid(idPair)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            XTR(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            XTR(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            XTR(ng)%Vid(idUaiE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            XTR(ng)%Vid(idVaiN)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            XTR(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            XTR(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            XTR(ng)%Vid(idLrad)=var_id(i)
#  endif
#  ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            XTR(ng)%Vid(idTair)=var_id(i)
#   ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            XTR(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            XTR(ng)%Vid(idrain)=var_id(i)
#   endif
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            XTR(ng)%Vid(idEmPf)=var_id(i)
#  ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            XTR(ng)%Vid(idSrad)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            XTR(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            XTR(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            XTR(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            XTR(ng)%Vid(idVbms)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              XTR(ng)%Tid(itrc)=var_id(i)
#  ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              XTR(ng)%Vid(idSbry(isTvar(itrc)))=var_id(i)
#  endif
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              XTR(ng)%Vid(idTsur(itrc))=var_id(i)
#  if defined FORWARD_WRITE && defined LMD_NONLOCAL
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idGhat(itrc)))) THEN
              got_var(idGhat(itrc))=.TRUE.
              XTR(ng)%Vid(idGhat(itrc))=var_id(i)
#  endif
            END IF
          END DO
# endif
        END DO
!
!  Check if decimate variables are available in input NetCDF file.
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
        IF (.not.got_var(idpthR).and.Hout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Hout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Hout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Hout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idu2dE).and.Hout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Hout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef FORWARD_WRITE
#  ifdef FORWARD_RHS
        IF (.not.got_var(idRzet).and.Hout(idRzet,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRzet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRu2d).and.Hout(idRu2d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRu2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRv2d).and.Hout(idRv2d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRv2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef SOLVE3D
#   ifdef FORWARD_RHS
        IF (.not.got_var(idRuct).and.Hout(idRuct,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRuct)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRcvt).and.Hout(idRcvt,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRcvt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idUfx1).and.Hout(idUfx1,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2).and.Hout(idUfx2,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1).and.Hout(idVfx1,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2).and.Hout(idVfx2,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef FORWARD_RHS
        IF (.not.got_var(idRu3d).and.Hout(idRu3d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRu3d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRv3d).and.Hout(idRv3d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRv3d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
# endif
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idu3dE).and.Hout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Hout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef OMEGA_IMPLICIT
        IF (.not.got_var(idOvil).and.Hout(idOvil,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvil)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SALINITY
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Hout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Hout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Hout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef EMINUSP
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
        IF (.not.got_var(idEmPf).and.Hout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,70)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  if defined FORWARD_WRITE && defined LMD_NONLOCAL
          IF (.not.got_var(idGhat(itrc)).and.Hout(idGhat(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idGhat(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#  endif
        END DO
# endif

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for decimate variables. Get variable
!  IDs.
!
        CALL bbl_def_nf90 (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!  Scan sea-ice variables from input NetCDF and activate switches for
!  decimate variables. Get variable IDs.
!
        CALL ice_def_nf90 (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for decimate variables. Get variable IDs.
!
        CALL sediment_def_nf90 (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for decimate variables. Get variable IDs.
!
        CALL wec_def_nf90 (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefXTR(ng).gt.0) THEN
          XTR(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefXTR(ng)*((ntstart(ng)-1)/ndefXTR(ng)))/   &
     &                   nXTR(ng)
        ELSE
          XTR(ng)%Rindex=(ntstart(ng)-1)/nXTR(ng)
        END IF
        XTR(ng)%Rindex=MIN(XTR(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_EXTRACT_NF90 - creating decimate file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_EXTRACT_NF90 - inquiring decimate file,',t56,     &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_EXTRACT_NF90 - unable to create decimate NetCDF', &
     &        ' file: ',a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_EXTRACT_NF90 - unable to open decimate NetCDF',   &
     &        ' file: ',a)
  70  FORMAT (/,' DEF_EXTRACT_NF90 - unable to find variable: ',a,2x,   &
     &        ' in decimate NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_extract_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_extract_pio (ng, model, ldef)
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
      integer :: i, j, ifield, itrc, nvd3, nvd4
      integer :: recdim, status
#  ifdef ADJUST_BOUNDARY
      integer :: IorJdim, brecdim
#  endif
      integer :: DimIDs(nDimID)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)
#  ifdef ADJUST_BOUNDARY
      integer :: t2dobc(4)
#  endif

#  ifdef SOLVE3D
#   ifdef SEDIMENT
      integer :: b3dgrd(4)
#   endif
      integer :: p2grd(3), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#   ifdef ADJUST_BOUNDARY
      integer :: t3dobc(5)
#   endif
#  endif
!
      real(r8) :: Aval(6)
!
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_extract_pio"
!
      TYPE (Var_desc_t) :: varDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=XTR(ng)%name
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
!  Create a new decimate file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname), XTR(ng)%pioFile)
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
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xi_rho',    &
     &                 xtr_IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xi_u',      &
     &                 xtr_IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xi_v',      &
     &                 xtr_IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xi_psi',    &
     &                 xtr_IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'eta_rho',   &
     &                 xtr_IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'eta_u',     &
     &                 xtr_IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'eta_v',     &
     &                 xtr_IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'eta_psi',   &
     &                 xtr_IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'IorJ',      &
     &                 xtr_IOBOUNDS(ng)%IorJ, IorJdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xy_rho',    &
     &                 xtr_IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xy_u',      &
     &                 xtr_IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xy_v',      &
     &                 xtr_IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xyz_rho',   &
     &                 xtr_IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xyz_u',     &
     &                 xtr_IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xyz_v',     &
     &                 xtr_IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xyz_w',     &
     &                 xtr_IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'N',         &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 's_rho',     &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 's_w',       &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'tracer',    &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SEDIMENT
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'NST',       &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nbed',      &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'xybed',     &
     &                 xtr_IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif

#   ifdef ECOSIM
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nbands',    &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nphy',      &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nbac',      &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Ndom',      &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nfec',      &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'boundary',  &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FOUR_DVAR
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nstate',    &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef MULTI_SCALE_B
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'Nscale',    &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef ADJUST_BOUNDARY
        status=def_dim(ng, model, XTR(ng)%pioFile, ncname, 'obc_adjust',&
     &                 Nbrec(ng), DimIDs(31))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, model, XTR(ng)%pioFile, ncname,              &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
#  ifdef ADJUST_BOUNDARY
        brecdim=DimIDs(31)
#  endif
!
!  Set number of dimensions for output variables.
!
#  if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#  else
        nvd3=3
        nvd4=4
#  endif
!
!  Define dimension vectors for staggered tracer type variables.
!
#  if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(12)
#   endif
#  else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(12)
#   endif
#  endif
#  ifdef ADJUST_BOUNDARY
        t2dobc(1)=IorJdim
        t2dobc(2)=DimIDs(14)
        t2dobc(3)=brecdim
        t2dobc(4)=DimIDs(12)
#   ifdef SOLVE3D
        t3dobc(1)=IorJdim
        t3dobc(2)=DimIDs( 9)
        t3dobc(3)=DimIDs(14)
        t3dobc(4)=brecdim
        t3dobc(5)=DimIDs(12)
#   endif
#  endif
!
!  Define dimension vectors for staggered type variables at PSI-points.
!
        p2grd(1)=DimIDs( 4)
        p2grd(2)=DimIDs( 8)
        p2grd(3)=DimIDs(12)
!
!  Define dimension vectors for staggered u-momentum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(12)
#   endif
#  else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(12)
#   endif
#  endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(12)
#   endif
#  else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(12)
#   endif
#  endif
#  ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momentum type variables.
!
#   if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(12)
#   else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(12)
#   endif
#   ifdef SEDIMENT
!
!  Define dimension vector for sediment bed layer type variables.
!
#    if defined WRITE_WATER && defined MASKING
        b3dgrd(1)=DimIDs(24)
        b3dgrd(2)=DimIDs(12)
#    else
        b3dgrd(1)=DimIDs( 1)
        b3dgrd(2)=DimIDs( 5)
        b3dgrd(3)=DimIDs(16)
        b3dgrd(4)=DimIDs(12)
#    endif
#   endif
#  endif
!
!  Initialize unlimited time record dimension.
!
        XTR(ng)%Rindex=0
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
        CALL def_info (ng, model, XTR(ng)%pioFile, ncname, DimIDs)
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
        XTR(ng)%pioVar(idtime)%dkind=PIO_TOUT
        XTR(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, model, XTR(ng)%pioFile,                      &
     &                 XTR(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SO_SEMI
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='real Ritz eigenvalues'
        status=def_var(ng, model, XTR(ng)%pioFile, varDesc, PIO_TYPE,   &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ng, model, XTR(ng)%pioFile, varDesc, PIO_TYPE,   &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef WET_DRY
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
        XTR(ng)%pioVar(idPwet)%dkind=PIO_FOUT
        XTR(ng)%pioVar(idPwet)%gtype=p2dvar
!
        status=def_var(ng, model, XTR(ng)%pioFile,                      &
     &                 XTR(ng)%pioVar(idPwet)%vd,                       &
     &                 PIO_FOUT, nvd3, p2grd, Aval, Vinfo, ncname,      &
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
        XTR(ng)%pioVar(idRwet)%dkind=PIO_FOUT
        XTR(ng)%pioVar(idRwet)%gtype=r2dvar
!
        status=def_var(ng, model, XTR(ng)%pioFile,                      &
     &                 XTR(ng)%pioVar(idRwet)%vd,                       &
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
        XTR(ng)%pioVar(idUwet)%dkind=PIO_FOUT
        XTR(ng)%pioVar(idUwet)%gtype=u2dvar
!
        status=def_var(ng, model, XTR(ng)%pioFile,                      &
     &                 XTR(ng)%pioVar(idUwet)%vd,                       &
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
        XTR(ng)%pioVar(idVwet)%dkind=PIO_FOUT
        XTR(ng)%pioVar(idVwet)%gtype=v2dvar
!
        status=def_var(ng, model, XTR(ng)%pioFile,                      &
     &                 XTR(ng)%pioVar(idVwet)%vd,                       &
     &                 PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname,     &
     &                 SetFillVal = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef SOLVE3D
!
!  Define time-varying depth of RHO-points.
!
        IF (Hout(idpthR,ng)) THEN
          Vinfo( 1)=Vname(1,idpthR)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthR))
          Vinfo( 3)=Vname(3,idpthR)
          Vinfo(14)=Vname(4,idpthR)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idpthR)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthR,ng),r8)
          XTR(ng)%pioVar(idpthR)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idpthR)%gtype=r3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idpthR)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of U-points.
!
        IF (Hout(idpthU,ng)) THEN
          Vinfo( 1)=Vname(1,idpthU)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthU))
          Vinfo( 3)=Vname(3,idpthU)
          Vinfo(14)=Vname(4,idpthU)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idpthU)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthU,ng),r8)
          XTR(ng)%pioVar(idpthU)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idpthU)%gtype=u3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idpthU)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of V-points.
!
        IF (Hout(idpthV,ng)) THEN
          Vinfo( 1)=Vname(1,idpthV)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthV))
          Vinfo( 3)=Vname(3,idpthV)
          Vinfo(14)=Vname(4,idpthV)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idpthV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthV,ng),r8)
          XTR(ng)%pioVar(idpthV)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idpthV)%gtype=v3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idpthV)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define time-varying depth of W-points.
!
        IF (Hout(idpthW,ng)) THEN
          Vinfo( 1)=Vname(1,idpthW)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idpthW))
          Vinfo( 3)=Vname(3,idpthW)
          Vinfo(14)=Vname(4,idpthW)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idpthW)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idpthW,ng),r8)
          XTR(ng)%pioVar(idpthW)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idpthW)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idpthW)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          Vinfo( 2)=Vname(2,idFsur)
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          XTR(ng)%pioVar(idFsur)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idFsur)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idFsur)%vd,                     &
#  ifdef WET_DRY
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
#  else
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
#  endif
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined FORWARD_WRITE && defined FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRzet)
          Vinfo( 2)=Vname(2,idRzet)
          Vinfo( 3)=Vname(3,idRzet)
          Vinfo(14)=Vname(4,idRzet)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idRzet)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRzet,ng),r8)
          XTR(ng)%pioVar(idRzet)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRzet)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRzet)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
        END IF

#  ifdef ADJUST_BOUNDARY
!
!  Define free-surface open boundaries.
!
        IF (ANY(Lobc(:,isFsur,ng))) THEN
          ifield=idSbry(isFsur)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idFsur)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          XTR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          XTR(ng)%pioVar(ifield)%gtype=r2dobc
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 4, t2dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D U-momentum component.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          Vinfo( 2)=Vname(2,idUbar)
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          XTR(ng)%pioVar(idUbar)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUbar)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUbar)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FORWARD_WRITE
#   ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRu2d)
          Vinfo( 2)=Vname(2,idRu2d)
          Vinfo( 3)=Vname(3,idRu2d)
          Vinfo(14)=Vname(4,idRu2d)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idRu2d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRu2d,ng),r8)
          XTR(ng)%pioVar(idRu2d)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRu2d)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRu2d)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef SOLVE3D
#    ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRuct)
          Vinfo( 2)=Vname(2,idRuct)
          Vinfo( 3)=Vname(3,idRuct)
          Vinfo(14)=Vname(4,idRuct)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#     endif
          Vinfo(21)=Vname(6,idRuct)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRuct,ng),r8)
          XTR(ng)%pioVar(idRuct)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRuct)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRuct)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
!
          Vinfo( 1)=Vname(1,idUfx1)
          Vinfo( 2)=Vname(2,idUfx1)
          Vinfo( 3)=Vname(3,idUfx1)
          Vinfo(14)=Vname(4,idUfx1)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idUfx1)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUfx1,ng),r8)
          XTR(ng)%pioVar(idUfx1)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUfx1)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUfx1)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idUfx2)
          Vinfo( 2)=Vname(2,idUfx2)
          Vinfo( 3)=Vname(3,idUfx2)
          Vinfo(14)=Vname(4,idUfx2)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idUfx2)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUfx2,ng),r8)
          XTR(ng)%pioVar(idUfx2)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUfx2)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUfx2)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
        END IF

#  ifdef ADJUST_BOUNDARY
!
!  Define 2D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUbar,ng))) THEN
          ifield=idSbry(isUbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUbar)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          XTR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          XTR(ng)%pioVar(ifield)%gtype=u2dobc
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 4, t2dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D V-momentum component.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          Vinfo( 2)=Vname(2,idVbar)
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          XTR(ng)%pioVar(idVbar)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVbar)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVbar)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FORWARD_WRITE
#   ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRv2d)
          Vinfo( 2)=Vname(2,idRv2d)
          Vinfo( 3)=Vname(3,idRv2d)
          Vinfo(14)=Vname(4,idRv2d)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idRv2d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRv2d,ng),r8)
          XTR(ng)%pioVar(idRv2d)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRv2d)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRv2d)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#   ifdef SOLVE3D
#    ifdef FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRvct)
          Vinfo( 2)=Vname(2,idRvct)
          Vinfo( 3)=Vname(3,idRvct)
          Vinfo(14)=Vname(4,idRvct)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#     endif
          Vinfo(21)=Vname(6,idRvct)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRvct,ng),r8)
          XTR(ng)%pioVar(idRvct)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRvct)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRvct)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
!
          Vinfo( 1)=Vname(1,idVfx1)
          Vinfo( 2)=Vname(2,idVfx1)
          Vinfo( 3)=Vname(3,idVfx1)
          Vinfo(14)=Vname(4,idVfx1)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idVfx1)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVfx1,ng),r8)
          XTR(ng)%pioVar(idVfx1)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVfx1)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVfx1)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          Vinfo( 1)=Vname(1,idVfx2)
          Vinfo( 2)=Vname(2,idVfx2)
          Vinfo( 3)=Vname(3,idVfx2)
          Vinfo(14)=Vname(4,idVfx2)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idVfx2)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVfx2,ng),r8)
          XTR(ng)%pioVar(idVfx2)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVfx2)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVfx2)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
        END IF

#  ifdef ADJUST_BOUNDARY
!
!  Define 2D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVbar,ng))) THEN
          ifield=idSbry(isVbar)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVbar)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          XTR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          XTR(ng)%pioVar(ifield)%gtype=v2dobc
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 4, t2dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Hout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          Vinfo( 2)=Vname(2,idu2dE)
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          XTR(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idu2dE)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idu2dE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Hout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          Vinfo( 2)=Vname(2,idv2dN)
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          XTR(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idv2dN)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idv2dN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define 3D U-momentum component.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          Vinfo( 2)=Vname(2,idUvel)
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          XTR(ng)%pioVar(idUvel)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUvel)%gtype=u3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUvel)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined FORWARD_WRITE && defined FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRu3d)
          Vinfo( 2)=Vname(2,idRu3d)
          Vinfo( 3)=Vname(3,idRu3d)
          Vinfo(14)=Vname(4,idRu3d)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idRu3d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRu3d,ng),r8)
          XTR(ng)%pioVar(idRu3d)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRu3d)%gtype=u3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRu3d)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
        END IF

#   ifdef ADJUST_BOUNDARY
!
!  Define 3D U-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isUvel,ng))) THEN
          ifield=idSbry(isUvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idUvel)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          XTR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          XTR(ng)%pioVar(ifield)%gtype=u3dobc
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 5, t3dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define 3D V-momentum component.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          Vinfo( 2)=Vname(2,idVvel)
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          XTR(ng)%pioVar(idVvel)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVvel)%gtype=v3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVvel)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined FORWARD_WRITE && defined FORWARD_RHS
!
          Vinfo( 1)=Vname(1,idRv3d)
          Vinfo( 2)=Vname(2,idRv3d)
          Vinfo( 3)=Vname(3,idRv3d)
          Vinfo(14)=Vname(4,idRv3d)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idRv3d)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRv3d,ng),r8)
          XTR(ng)%pioVar(idRv3d)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idRv3d)%gtype=v3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idRv3d)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
        END IF

#   ifdef ADJUST_BOUNDARY
!
!  Define 3D V-momentum component open boundaries.
!
        IF (ANY(Lobc(:,isVvel,ng))) THEN
          ifield=idSbry(isVvel)
          Vinfo( 1)=Vname(1,ifield)
          Vinfo( 2)=Vname(2,ifield)
          Vinfo( 3)=Vname(3,ifield)
          Vinfo(14)=Vname(4,ifield)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvel)
          Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
          XTR(ng)%pioVar(ifield)%dkind=PIO_FOUT
          XTR(ng)%pioVar(ifield)%gtype=v3dobc
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(ifield)%vd,                     &
     &                   PIO_FOUT, 5, t3dobc, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Hout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          Vinfo( 2)=Vname(2,idu3dE)
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          XTR(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idu3dE)%gtype=r3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idu3dE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Hout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          Vinfo( 2)=Vname(2,idv3dN)
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          XTR(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idv3dN)%gtype=r3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idv3dN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D momentum component in the Z-direction.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          Vinfo( 2)=Vname(2,idWvel)
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          XTR(ng)%pioVar(idWvel)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idWvel)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idWvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          Vinfo( 2)=Vname(2,idOvel)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          XTR(ng)%pioVar(idOvel)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idOvel)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idOvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef OMEGA_IMPLICIT
!
!  Define S-coordinate implicit vertical "omega" momentum component.
!
        IF (Hout(idOvil,ng)) THEN
          Vinfo( 1)=Vname(1,idOvil)
          Vinfo( 2)=Vname(2,idOvil)
          Vinfo( 3)='meter second-1'
          Vinfo(14)=Vname(4,idOvil)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idOvil)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvil,ng),r8)
          XTR(ng)%pioVar(idOvil)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idOvil)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idOvil)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#   endif
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTvar(itrc),ng),r8)
            XTR(ng)%pioTrc(itrc)%dkind=PIO_FOUT
            XTR(ng)%pioTrc(itrc)%gtype=r3dvar
!
            status=def_var(ng, model, XTR(ng)%pioFile,                  &
     &                     XTR(ng)%pioTrc(itrc)%vd,                     &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#   ifdef ADJUST_BOUNDARY
!
!  Define tracer type variables open boundaries.
!
        DO itrc=1,NT(ng)
          IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
            ifield=idSbry(isTvar(itrc))
            Vinfo( 1)=Vname(1,ifield)
            Vinfo( 2)=Vname(2,ifield)
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),60) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#    endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Aval(5)=REAL(Iinfo(1,ifield,ng),r8)
            XTR(ng)%pioVar(ifield)%dkind=PIO_FOUT
            XTR(ng)%pioVar(ifield)%gtype=r3dobc
!
            status=def_var(ng, model, XTR(ng)%pioFile,                  &
     &                     XTR(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, 5, t3dobc, Aval, Vinfo, ncname,    &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          Vinfo( 2)=Vname(2,idDano)
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          XTR(ng)%pioVar(idDano)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idDano)%gtype=r3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idDano)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          Vinfo( 2)=Vname(2,idHsbl)
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idHsbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          XTR(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idHsbl)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idHsbl)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          Vinfo( 2)=Vname(2,idHbbl)
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idHbbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          XTR(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idHbbl)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idHbbl)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Define out KPP nonlocal transport.
!
        DO itrc=1,NAT
          IF (Hout(idGhat(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idGhat(itrc))
            Vinfo( 2)=Vname(2,idGhat(itrc))
            Vinfo( 3)=Vname(3,idGhat(itrc))
            Vinfo(14)=Vname(4,idGhat(itrc))
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,idGhat(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idGhat(itrc),ng),r8)
            XTR(ng)%pioVar(idGhat(itrc))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idGhat(itrc))%gtype=w3dvar
!
            status=def_var(ng, model, XTR(ng)%pioFile,                  &
     &                     XTR(ng)%pioVar(idGhat(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          Vinfo( 2)=Vname(2,idVvis)
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          XTR(ng)%pioVar(idVvis)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVvis)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVvis)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          Vinfo( 2)=Vname(2,idTdif)
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          XTR(ng)%pioVar(idTdif)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idTdif)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idTdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          Vinfo( 2)=Vname(2,idSdif)
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          XTR(ng)%pioVar(idSdif)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idSdif)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idSdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
!
!  Define turbulent kinetic energy.
!
        IF (Hout(idMtke,ng)) THEN
          Vinfo( 1)=Vname(1,idMtke)
          Vinfo( 2)=Vname(2,idMtke)
          Vinfo( 3)=Vname(3,idMtke)
          Vinfo(14)=Vname(4,idMtke)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idMtke)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtke,ng),r8)
          XTR(ng)%pioVar(idMtke)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idMtke)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idMtke)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    ifdef FORWARD_WRITE
!
          Vinfo( 1)=Vname(1,idVmKK)
          Vinfo( 2)=Vname(2,idVmKK)
          Vinfo( 3)=Vname(3,idVmKK)
          Vinfo(14)=Vname(4,idVmKK)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idVmKK)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKK,ng),r8)
          XTR(ng)%pioVar(idVmKK)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVmKK)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVmKK)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
        END IF
!
!  Define turbulent kinetic energy time length scale.
!
        IF (Hout(idMtls,ng)) THEN
          Vinfo( 1)=Vname(1,idMtls)
          Vinfo( 2)=Vname(2,idMtls)
          Vinfo( 3)=Vname(3,idMtls)
          Vinfo(14)=Vname(4,idMtls)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idMtls)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idMtls,ng),r8)
          XTR(ng)%pioVar(idMtls)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idMtls)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idMtls)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    ifdef FORWARD_WRITE
!
          Vinfo( 1)=Vname(1,idVmLS)
          Vinfo( 2)=Vname(2,idVmLS)
          Vinfo( 3)=Vname(3,idVmLS)
          Vinfo(14)=Vname(4,idVmLS)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idVmLS)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmLS,ng),r8)
          XTR(ng)%pioVar(idVmLS)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVmLS)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVmLS)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#    if defined FORWARD_WRITE && defined GLS_MIXING
!
          Vinfo( 1)=Vname(1,idVmKP)
          Vinfo( 2)=Vname(2,idVmKP)
          Vinfo( 3)=Vname(3,idVmKP)
          Vinfo(14)=Vname(4,idVmKP)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idVmKP)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVmKP,ng),r8)
          XTR(ng)%pioVar(idVmKP)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVmKP)%gtype=w3dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVmKP)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Hout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          Vinfo( 2)=Vname(2,idPair)
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          XTR(ng)%pioVar(idPair)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idPair)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idPair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Hout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          Vinfo( 2)=Vname(2,idUair)
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          XTR(ng)%pioVar(idUair)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUair)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Hout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          Vinfo( 2)=Vname(2,idVair)
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          XTR(ng)%pioVar(idVair)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVair)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface winds at RHO-points.
!
        IF (Hout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          Vinfo( 2)=Vname(2,idUaiE)
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          XTR(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUaiE)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUaiE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Hout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          Vinfo( 2)=Vname(2,idVaiN)
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          XTR(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVaiN)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVaiN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define surface active tracer fluxes.
!
        DO itrc=1,NAT
          IF (Hout(idTsur(itrc),ng)) THEN
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
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTsur(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            XTR(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
!
            status=def_var(ng, model, XTR(ng)%pioFile,                  &
     &                     XTR(ng)%pioVar(idTsur(itrc))%vd,             &
     &                     PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#   if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          Vinfo( 2)=Vname(2,idLhea)
          Vinfo( 3)=Vname(3,idLhea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLhea)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idLhea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          XTR(ng)%pioVar(idLhea)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idLhea)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idLhea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          Vinfo( 2)=Vname(2,idShea)
          Vinfo( 3)=Vname(3,idShea)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idShea)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idShea)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          XTR(ng)%pioVar(idShea)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idShea)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idShea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          Vinfo( 2)=Vname(2,idLrad)
          Vinfo( 3)=Vname(3,idLrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idLrad)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idLrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          XTR(ng)%pioVar(idLrad)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idLrad)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idLrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   if defined BULK_FLUXES
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          Vinfo( 2)=Vname(2,idTair)
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          XTR(ng)%pioVar(idTair)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idTair)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idTair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#    ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          Vinfo( 2)=Vname(2,idevap)
          Vinfo( 3)=Vname(3,idevap)
          Vinfo(11)='downward flux, freshening (condensation)'
          Vinfo(12)='upward flux, salting (evaporation)'
          Vinfo(14)=Vname(4,idevap)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idevap)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          XTR(ng)%pioVar(idevap)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idevap)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idevap)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          Vinfo( 2)=Vname(2,idrain)
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(11)='upward flux, salting (NOT POSSIBLE)'
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          XTR(ng)%pioVar(idrain)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idrain)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idrain)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#    endif
#   endif
!
!  Define E-P flux.
!
        IF (Hout(idEmPf,ng)) THEN
          Vinfo( 1)=Vname(1,idEmPf)
          Vinfo( 2)=Vname(2,idEmPf)
          Vinfo( 3)=Vname(3,idEmPf)
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idEmPf)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idEmPf)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idEmPf,ng),r8)
          XTR(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idEmPf)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idEmPf)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SHORTWAVE
!
!  Define net shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          Vinfo( 2)=Vname(2,idSrad)
          Vinfo( 3)=Vname(3,idSrad)
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idSrad)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idSrad)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          XTR(ng)%pioVar(idSrad)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idSrad)%gtype=r2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idSrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          Vinfo( 2)=Vname(2,idUsms)
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          XTR(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUsms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          Vinfo( 2)=Vname(2,idVsms)
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          XTR(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVsms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          Vinfo( 2)=Vname(2,idUbms)
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          XTR(ng)%pioVar(idUbms)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idUbms)%gtype=u2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idUbms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          Vinfo( 2)=Vname(2,idVbms)
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          XTR(ng)%pioVar(idVbms)%dkind=PIO_FOUT
          XTR(ng)%pioVar(idVbms)%gtype=v2dvar
!
          status=def_var(ng, model, XTR(ng)%pioFile,                    &
     &                   XTR(ng)%pioVar(idVbms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_pio (ng, model, ldef, Hout, XTR,                   &
     &                    t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_pio (ng, model, ldef, Hout, XTR,                   &
     &                    t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_pio (ng, model, ldef, Hout, XTR,              &
     &                         t2dgrd, u2dgrd, v2dgrd,                  &
     &                         b3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WEC_VF && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the Waves Effect on Currents variables.
!-----------------------------------------------------------------------
!
        CALL wec_def_pio (ng, model, ldef, Hout, XTR,                   &
     &                    t2dgrd, u2dgrd, v2dgrd,                       &
     &                    t3dgrd, u3dgrd, v3dgrd, w3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, XTR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, XTR(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing decimate file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=XTR(ng)%name
!
!  Open decimate file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1, XTR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,60) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = XTR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = XTR(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  decimate variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            XTR(ng)%pioVar(idtime)%vd=var_desc(i)
            XTR(ng)%pioVar(idtime)%dkind=PIO_TOUT
            XTR(ng)%pioVar(idtime)%gtype=0
#  if defined WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPwet))) THEN
            got_var(idPwet)=.TRUE.
            XTR(ng)%pioVar(idPwet)%vd=var_desc(i)
            XTR(ng)%pioVar(idPwet)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idPwet)%gtype=p2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            XTR(ng)%pioVar(idRwet)%vd=var_desc(i)
            XTR(ng)%pioVar(idRwet)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRwet)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwet))) THEN
            got_var(idUwet)=.TRUE.
            XTR(ng)%pioVar(idUwet)%vd=var_desc(i)
            XTR(ng)%pioVar(idUwet)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUwet)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwet))) THEN
            got_var(idVwet)=.TRUE.
            XTR(ng)%pioVar(idVwet)%vd=var_desc(i)
            XTR(ng)%pioVar(idVwet)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVwet)%gtype=v2dvar
#  endif
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthR))) THEN
            got_var(idpthR)=.TRUE.
            XTR(ng)%pioVar(idpthR)%vd=var_desc(i)
            XTR(ng)%pioVar(idpthR)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idpthR)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthU))) THEN
            got_var(idpthU)=.TRUE.
            XTR(ng)%pioVar(idpthU)%vd=var_desc(i)
            XTR(ng)%pioVar(idpthU)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idpthU)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthV))) THEN
            got_var(idpthV)=.TRUE.
            XTR(ng)%pioVar(idpthV)%vd=var_desc(i)
            XTR(ng)%pioVar(idpthV)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idpthV)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idpthW))) THEN
            got_var(idpthW)=.TRUE.
            XTR(ng)%pioVar(idpthW)%vd=var_desc(i)
            XTR(ng)%pioVar(idpthW)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idpthW)%gtype=w3dvar
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            XTR(ng)%pioVar(idFsur)%vd=var_desc(i)
            XTR(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            XTR(ng)%pioVar(idUbar)%vd=var_desc(i)
            XTR(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            XTR(ng)%pioVar(idVbar)%vd=var_desc(i)
            XTR(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVbar)%gtype=v2dvar
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isFsur)))) THEN
            got_var(idSbry(isFsur))=.TRUE.
            XTR(ng)%pioVar(idSbry(isFsur))%vd=var_desc(i)
            XTR(ng)%pioVar(idSbry(isFsur))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSbry(isFsur))%gtype=r2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUbar)))) THEN
            got_var(idSbry(isUbar))=.TRUE.
            XTR(ng)%pioVar(idSbry(isUbar))%vd=var_desc(i)
            XTR(ng)%pioVar(idSbry(isUbar))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSbry(isUbar))%gtype=u2dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVbar)))) THEN
            got_var(idSbry(isVbar))=.TRUE.
            XTR(ng)%pioVar(idSbry(isVbar))%vd=var_desc(i)
            XTR(ng)%pioVar(idSbry(isVbar))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSbry(isVbar))%gtype=v2dobc
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            XTR(ng)%pioVar(idu2dE)%vd=var_desc(i)
            XTR(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idu2dE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            XTR(ng)%pioVar(idv2dN)%vd=var_desc(i)
            XTR(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idv2dN)%gtype=r2dvar
#  ifdef FORWARD_WRITE
#   ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRzet))) THEN
            got_var(idRzet)=.TRUE.
            XTR(ng)%pioVar(idRzet)%vd=var_desc(i)
            XTR(ng)%pioVar(idRzet)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRzet)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            got_var(idRu2d)=.TRUE.
            XTR(ng)%pioVar(idRu2d)%vd=var_desc(i)
            XTR(ng)%pioVar(idRu2d)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRu2d)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            got_var(idRv2d)=.TRUE.
            XTR(ng)%pioVar(idRv2d)%vd=var_desc(i)
            XTR(ng)%pioVar(idRv2d)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRv2d)%gtype=v2dvar
#   endif
#   ifdef SOLVE3D
#    ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRuct))) THEN
            got_var(idRuct)=.TRUE.
            XTR(ng)%pioVar(idRuct)%vd=var_desc(i)
            XTR(ng)%pioVar(idRuct)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRuct)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRvct))) THEN
            got_var(idRvct)=.TRUE.
            XTR(ng)%pioVar(idRvct)%vd=var_desc(i)
            XTR(ng)%pioVar(idRvct)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRvct)%gtype=v2dvar
#    endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            XTR(ng)%pioVar(idUfx1)%vd=var_desc(i)
            XTR(ng)%pioVar(idUfx1)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUfx1)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            XTR(ng)%pioVar(idUfx2)%vd=var_desc(i)
            XTR(ng)%pioVar(idUfx2)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUfx2)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            XTR(ng)%pioVar(idVfx1)%vd=var_desc(i)
            XTR(ng)%pioVar(idVfx1)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVfx1)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            XTR(ng)%pioVar(idVfx2)%vd=var_desc(i)
            XTR(ng)%pioVar(idVfx2)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVfx2)%gtype=v2dvar
#    ifdef FORWARD_RHS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            got_var(idRu3d)=.TRUE.
            XTR(ng)%pioVar(idRu3d)%vd=var_desc(i)
            XTR(ng)%pioVar(idRu3d)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRu3d)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            got_var(idRv3d)=.TRUE.
            XTR(ng)%pioVar(idRv3d)%vd=var_desc(i)
            XTR(ng)%pioVar(idRv3d)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idRv3d)%gtype=v3dvar
#    endif
#   endif
#  endif
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            XTR(ng)%pioVar(idUvel)%vd=var_desc(i)
            XTR(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            XTR(ng)%pioVar(idVvel)%vd=var_desc(i)
            XTR(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVvel)%gtype=v3dvar
#   ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isUvel)))) THEN
            got_var(idSbry(isUvel))=.TRUE.
            XTR(ng)%pioVar(idSbry(isUvel))%vd=var_desc(i)
            XTR(ng)%pioVar(idSbry(isUvel))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSbry(isUvel))%gtype=u3dobc
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idSbry(isVvel)))) THEN
            got_var(idSbry(isVvel))=.TRUE.
            XTR(ng)%pioVar(idSbry(isVvel))%vd=var_desc(i)
            XTR(ng)%pioVar(idSbry(isVvel))%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSbry(isVvel))%gtype=v3dobc
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            XTR(ng)%pioVar(idu3dE)%vd=var_desc(i)
            XTR(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idu3dE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            XTR(ng)%pioVar(idv3dN)%vd=var_desc(i)
            XTR(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idv3dN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            XTR(ng)%pioVar(idWvel)%vd=var_desc(i)
            XTR(ng)%pioVar(idWvel)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idWvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            XTR(ng)%pioVar(idOvel)%vd=var_desc(i)
            XTR(ng)%pioVar(idOvel)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idOvel)%gtype=w3dvar
#   ifdef OMEGA_IMPLICIT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvil))) THEN
            got_var(idOvil)=.TRUE.
            XTR(ng)%pioVar(idOvil)%vd=var_desc(i)
            XTR(ng)%pioVar(idOvil)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idOvil)%gtype=w3dvar
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            XTR(ng)%pioVar(idDano)%vd=var_desc(i)
            XTR(ng)%pioVar(idDano)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idDano)%gtype=r3dvar
#   ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            XTR(ng)%pioVar(idHsbl)%vd=var_desc(i)
            XTR(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idHsbl)%gtype=r2dvar
#   endif
#   ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            XTR(ng)%pioVar(idHbbl)%vd=var_desc(i)
            XTR(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idHbbl)%gtype=r2dvar
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            XTR(ng)%pioVar(idVvis)%vd=var_desc(i)
            XTR(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVvis)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            XTR(ng)%pioVar(idTdif)%vd=var_desc(i)
            XTR(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idTdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            XTR(ng)%pioVar(idSdif)%vd=var_desc(i)
            XTR(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSdif)%gtype=w3dvar
#   if defined GLS_MIXING || defined MY25_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtke))) THEN
            got_var(idMtke)=.TRUE.
            XTR(ng)%pioVar(idMtke)%vd=var_desc(i)
            XTR(ng)%pioVar(idMtke)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idMtke)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idMtls))) THEN
            got_var(idMtls)=.TRUE.
            XTR(ng)%pioVar(idMtls)%vd=var_desc(i)
            XTR(ng)%pioVar(idMtls)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idMtls)%gtype=w3dvar
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            XTR(ng)%pioVar(idPair)%vd=var_desc(i)
            XTR(ng)%pioVar(idPair)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idPair)%gtype=r2dvar
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            XTR(ng)%pioVar(idUair)%vd=var_desc(i)
            XTR(ng)%pioVar(idUair)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            XTR(ng)%pioVar(idVair)%vd=var_desc(i)
            XTR(ng)%pioVar(idVair)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUair)=.TRUE.
            XTR(ng)%pioVar(idUaiE)%vd=var_desc(i)
            XTR(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUaiE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVair)=.TRUE.
            XTR(ng)%pioVar(idVaiN)%vd=var_desc(i)
            XTR(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVaiN)%gtype=r2dvar
#   endif
#   if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            XTR(ng)%pioVar(idLhea)%vd=var_desc(i)
            XTR(ng)%pioVar(idLhea)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idLhea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            XTR(ng)%pioVar(idShea)%vd=var_desc(i)
            XTR(ng)%pioVar(idShea)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idShea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            XTR(ng)%pioVar(idLrad)%vd=var_desc(i)
            XTR(ng)%pioVar(idLrad)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idLrad)%gtype=r2dvar
#   endif
#   ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            XTR(ng)%pioVar(idTair)%vd=var_desc(i)
            XTR(ng)%pioVar(idTair)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idTair)%gtype=r2dvar
#    ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            XTR(ng)%pioVar(idevap)%vd=var_desc(i)
            XTR(ng)%pioVar(idevap)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idevap)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            XTR(ng)%pioVar(idrain)%vd=var_desc(i)
            XTR(ng)%pioVar(idrain)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idrain)%gtype=r2dvar
#    endif
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idEmPf))) THEN
            got_var(idEmPf)=.TRUE.
            XTR(ng)%pioVar(idEmPf)%vd=var_desc(i)
            XTR(ng)%pioVar(idEmPf)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idEmPf)%gtype=r2dvar
#   ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            XTR(ng)%pioVar(idSrad)%vd=var_desc(i)
            XTR(ng)%pioVar(idSrad)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idSrad)%gtype=r2dvar
#   endif
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            XTR(ng)%pioVar(idUsms)%vd=var_desc(i)
            XTR(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            XTR(ng)%pioVar(idVsms)%vd=var_desc(i)
            XTR(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVsms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            XTR(ng)%pioVar(idUbms)%vd=var_desc(i)
            XTR(ng)%pioVar(idUbms)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idUbms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            XTR(ng)%pioVar(idVbms)%vd=var_desc(i)
            XTR(ng)%pioVar(idVbms)%dkind=PIO_FOUT
            XTR(ng)%pioVar(idVbms)%gtype=v2dvar
          END IF
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              XTR(ng)%pioTrc(itrc)%vd=var_desc(i)
              XTR(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              XTR(ng)%pioTrc(itrc)%gtype=r3dvar
#   ifdef ADJUST_BOUNDARY
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idSbry(isTvar(itrc))))) THEN
              got_var(idSbry(isTvar(itrc)))=.TRUE.
              XTR(ng)%pioVar(idSbry(isTvar(itrc)))%vd=var_desc(i)
              XTR(ng)%pioVar(idSbry(isTvar(itrc)))%dkind=PIO_FOUT
              XTR(ng)%pioVar(idSbry(isTvar(itrc)))%gtype=r3dobc
#   endif
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              XTR(ng)%pioVar(idTsur(itrc))%vd=var_desc(i)
              XTR(ng)%pioVar(idTsur(itrc))%dkind=PIO_FOUT
              XTR(ng)%pioVar(idTsur(itrc))%gtype=r2dvar
#   if defined FORWARD_WRITE && defined LMD_NONLOCAL
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idGhat(itrc)))) THEN
              got_var(idGhat(itrc))=.TRUE.
              XTR(ng)%pioVar(idGhat(itrc))%vd=var_desc(i)
              XTR(ng)%pioVar(idGhat(itrc))%dkind=PIO_FOUT
              XTR(ng)%pioVar(idGhat(itrc))%gtype=w3dvar
#   endif
            END IF
          END DO
#  endif
        END DO
!
!  Check if decimate variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined WET_DRY
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
#  endif
#  ifdef SOLVE3D
        IF (.not.got_var(idpthR).and.Hout(idpthR,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthR)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthU).and.Hout(idpthU,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthU)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthV).and.Hout(idpthV,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idpthW).and.Hout(idpthW,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idpthW)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isFsur))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isFsur))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isUbar))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isUbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVbar))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isVbar))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idu2dE).and.Hout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Hout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef FORWARD_WRITE
#   ifdef FORWARD_RHS
        IF (.not.got_var(idRzet).and.Hout(idRzet,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRzet)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRu2d).and.Hout(idRu2d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRu2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRv2d).and.Hout(idRv2d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRv2d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef SOLVE3D
#    ifdef FORWARD_RHS
        IF (.not.got_var(idRuct).and.Hout(idRuct,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRuct)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRcvt).and.Hout(idRcvt,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRcvt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
        IF (.not.got_var(idUfx1).and.Hout(idUfx1,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2).and.Hout(idUfx2,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1).and.Hout(idVfx1,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2).and.Hout(idVfx2,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef FORWARD_RHS
        IF (.not.got_var(idRu3d).and.Hout(idRu3d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRu3d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idRv3d).and.Hout(idRv3d,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idRv3d)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#  endif
#  ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef ADJUST_BOUNDARY
        IF (.not.got_var(idSbry(isUvel))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isUvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSbry(isVvel))) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSbry(isVvel))),  &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idu3dE).and.Hout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Hout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef OMEGA_IMPLICIT
        IF (.not.got_var(idOvil).and.Hout(idOvil,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idOvil)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idVvis).and.Hout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Hout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SALINITY
        IF (.not.got_var(idSdif).and.Hout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.got_var(idMtke).and.Hout(idMtke,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtke)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idMtls).and.Hout(idMtls,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idMtls)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Hout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Hout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Hout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Hout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Hout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined FRC_COUPLING
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef BULK_FLUXES
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef EMINUSP
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
        IF (.not.got_var(idEmPf).and.Hout(idEmPf,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idEmPf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   ifdef ADJUST_BOUNDARY
          IF (.not.got_var(idSbry(isTvar(itrc)))) THEN
            IF (Master) WRITE (stdout,70)                               &
     &                        TRIM(Vname(1,idSbry(isTvar(itrc)))),      &
     &                        TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   endif
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   if defined FORWARD_WRITE && defined LMD_NONLOCAL
          IF (.not.got_var(idGhat(itrc)).and.Hout(idGhat(itrc),ng)) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idGhat(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
#   endif
        END DO
#  endif

#  if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for decimate variables. Get variable
!  IDs.
!
        CALL bbl_def_pio (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ICE_MODEL && defined SOLVE3D
!
!  Scan sea-ice variables from input NetCDF and activate switches for
!  decimate variables. Get variable IDs.
!
        CALL ice_def_pio (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SEDIMENT && defined SOLVE3D
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for decimate variables. Get variable IDs.
!
        CALL sediment_def_pio (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined WEC || defined WEC_VF
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for decimate variables. Get variable IDs.
!
        CALL wec_def_pio (ng, model, ldef, Hout, XTR)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefXTR(ng).gt.0) THEN
          XTR(ng)%Rindex=((ntstart(ng)-1)-                              &
     &                    ndefXTR(ng)*((ntstart(ng)-1)/ndefXTR(ng)))/   &
     &                   nXTR(ng)
        ELSE
          XTR(ng)%Rindex=(ntstart(ng)-1)/nXTR(ng)
        END IF
        XTR(ng)%Rindex=MIN(XTR(ng)%Rindex,rec_size)
      END IF QUERY
!
  10  FORMAT (2x,'DEF_EXTRACT_PIO  - creating decimate file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_EXTRACT_PIO  - inquiring decimate file,',t56,     &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_EXTRACT_PIO - unable to create decimate NetCDF',  &
     &        ' file: ',a)
  40  FORMAT ('time dependent',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' DEF_EXTRACT_PIO - unable to open decimate NetCDF',    &
     &        ' file: ', a)
  70  FORMAT (/,' DEF_EXTRACT_PIO - unable to find variable: ',a,2x,    &
     &        ' in decimate NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_extract_pio
# endif
#endif
!
      END MODULE def_extract_mod
