#include "cppdefs.h"
      MODULE def_avg_mod
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module creates output TIME-AVERAGED file using either the      !
!  standard NetCDF library or the Parallel-IO (PIO) library.  It       !
!  defines its dimensions, attributes, and variables.                  !
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
      USE mod_netcdf
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
      PUBLIC  :: def_avg
      PRIVATE :: def_avg_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PRIVATE :: def_avg_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE def_avg (ng, ldef)
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
      integer :: model
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
!-----------------------------------------------------------------------
!  Create a new time-averaged file according to IO type.
!-----------------------------------------------------------------------
!
# ifdef ADJOINT
      model=iADM
# else
      model=iNLM
# endif
!
      SELECT CASE (AVG(ng)%IOtype)
        CASE (io_nf90)
          CALL def_avg_nf90 (ng, model, ldef)

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL def_avg_pio (ng, model, ldef)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) AVG(ng)%IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' DEF_AVG - Illegal output file type, io_type = ',i0,     &
     &        /,11x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE def_avg
!
!***********************************************************************
      SUBROUTINE def_avg_nf90 (ng, model, ldef)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, ifield, itrc, j, nvd3, nvd4
# if defined WRITE_WATER && defined MASKING
      integer :: xy_pdim, xyz_pdim
# endif
      integer :: recdim, status

      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)

!#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!      integer :: r3degrd(4)
!#endif
#if defined WEC_VF && defined SOLVE3D
      integer :: a3degrd(4)
      integer :: u3degrd(4)
      integer :: v3degrd(4)
      integer :: t3degrd(4)
#endif

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: p3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif
!
      real(r8) :: Aval(6)
!
# ifdef ADJOINT
      character (len=21)     :: Prefix
# else
      character (len=13)     :: Prefix
# endif
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_avg_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=AVG(ng)%name
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
!  Create a new averages file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), AVG(ng)%ncid)
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
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_rho',       &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_u',         &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_v',         &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_psi',       &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_rho',      &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_u',        &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_v',        &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_psi',      &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_rho',       &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_u',         &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_v',         &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_psi',       &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

#if defined INWAVE_MODEL
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'ND',           &
     &                 ND, DimIDs(37))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'NDp1',         &
     &                 ND+1, DimIDs(38))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

#if defined UV_BANIHASHEMI
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'MSCs',         &
     &                 MSCs, DimIDs(39))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_rho',      &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_u',        &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_v',        &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_w',        &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_psi',      &
     &                 IOBOUNDS(ng)%xy_psi*N(ng), xyz_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'NST',          &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xybed',        &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif
#  if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
#  ifdef ECOSIM
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nbands',       &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
# endif

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# ifdef MULTI_SCALE_B
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nscale',       &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

        status=def_dim(ng, model, AVG(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
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
#  endif
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
!
!  Define dimension vectors for staggered variables at PSI-points.
!
# if defined WRITE_WATER && defined MASKING
        p2dgrd(1)=xy_pdim
        p2dgrd(2)=DimIDs(12)
#  ifdef SOLVE3D
        p3dgrd(1)=xyz_pdim
        p3dgrd(2)=DimIDs(12)
#  endif
# else
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
#  ifdef SOLVE3D
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
#  endif
# endif
!#if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!        r3degrd(1)=DimIDs( 1)
!        r3degrd(2)=DimIDs( 5)
!        r3degrd(3)=DimIDs(39)
!        r3degrd(4)=DimIDs(12)
!#endif
# if defined WEC_VF && defined SOLVE3D
!
!  Define dimension vector for energy wave layer type variables.
!
        a3degrd(1)=DimIDs( 1)
        a3degrd(2)=DimIDs( 5)
        a3degrd(3)=DimIDs(37)
        a3degrd(4)=DimIDs(12)

        u3degrd(1)=DimIDs( 2)
        u3degrd(2)=DimIDs( 6)
        u3degrd(3)=DimIDs(37)
        u3degrd(4)=DimIDs(12)

        v3degrd(1)=DimIDs( 3)
        v3degrd(2)=DimIDs( 7)
        v3degrd(3)=DimIDs(37)
        v3degrd(4)=DimIDs(12)

        t3degrd(1)=DimIDs( 1)
        t3degrd(2)=DimIDs( 5)
        t3degrd(3)=DimIDs(38)
        t3degrd(4)=DimIDs(12)
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
        AVG(ng)%Rindex=0
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
!  Set long name prefix string.
!
# ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
# else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
# endif
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, AVG(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idtime),    &
     &                 NF_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface.
!
        IF (Aout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if !defined WET_DRY && (defined WRITE_WATER && defined MASKING)
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idFsur),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided free-surface.
!
        IF (Aout(idFsuD,ng)) THEN
          Vinfo( 1)=Vname(1,idFsuD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsuD))
          Vinfo( 3)=Vname(3,idFsuD)
          Vinfo(14)=Vname(4,idFsuD)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idFsuD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsuD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idFsuD),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D momentum in the XI-direction.
!
        IF (Aout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbar),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 2D momentum in the XI-direction.
!
        IF (Aout(idu2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dD))
          Vinfo( 3)=Vname(3,idu2dD)
          Vinfo(14)=Vname(4,idu2dD)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idu2dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu2dD),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Aout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbar),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 2D momentum in the ETA-direction.
!
        IF (Aout(idv2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dD))
          Vinfo( 3)=Vname(3,idv2dD)
          Vinfo(14)=Vname(4,idv2dD)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idv2dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv2dD),  &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Aout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dE))
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu2dE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Aout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dN))
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv2dN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
!
!  Define time-averaged mass fluxes for 3D momentum coupling.
!
        IF (Aout(idUfx1,ng)) THEN
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
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUfx1),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idUfx2,ng)) THEN
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
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUfx2),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVfx1,ng)) THEN
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
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVfx1),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVfx2,ng)) THEN
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
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVfx2),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Aout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUvel),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 3D momentum component in the XI-direction.
!
        IF (Aout(idu3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dD))
          Vinfo( 3)=Vname(3,idu3dD)
          Vinfo(14)=Vname(4,idu3dD)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idu3dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu3dD),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Aout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVvel),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 3D momentum component in the ETA-direction.
!
        IF (Aout(idv3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dD))
          Vinfo( 3)=Vname(3,idv3dD)
          Vinfo(14)=Vname(4,idv3dD)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idv3dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv3dD),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Aout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dE))
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu3dE),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Aout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dN))
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv3dN),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Aout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOvel))
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idOvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define "true" vertical momentum component.
!
        IF (Aout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvel))
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Tid(itrc),  &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided temperature and salinity.
!
        DO itrc=1,NAT
          IF (Aout(idTrcD(itrc),ng)) THEN
            ifield=idTrcD(itrc)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(ifield), NF_FOUT,                &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#  endif
!
!  Define density anomaly.
!
        IF (Aout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idDano),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  ifdef NEMURO_SED1
!
!  Define PON in sediment.
!
        IF (Aout(idPONsed,ng)) THEN
          Vinfo( 1)=Vname(1,idPONsed)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPONsed))
          Vinfo( 3)=Vname(3,idPONsed)
          Vinfo(14)=Vname(4,idPONsed)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPONsed,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idPONsed),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define OPAL in sediment.
!
        IF (Aout(idOPALsed,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALsed)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOPALsed))
          Vinfo( 3)=Vname(3,idOPALsed)
          Vinfo(14)=Vname(4,idOPALsed)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOPALsed,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,AVG(ng)%Vid(idOPALsed),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define DENIT in sediment.
!
        IF (Aout(idDENITsed,ng)) THEN
          Vinfo( 1)=Vname(1,idDENITsed)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDENITsed))
          Vinfo( 3)=Vname(3,idDENITsed)
          Vinfo(14)=Vname(4,idDENITsed)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDENITsed,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idDENITsed),                       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define buried PON in sediment.
!
        IF (Aout(idPONbur,ng)) THEN
          Vinfo( 1)=Vname(1,idPONbur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPONbur))
          Vinfo( 3)=Vname(3,idPONbur)
          Vinfo(14)=Vname(4,idPONbur)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPONbur,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,AVG(ng)%Vid(idPONbur), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define buried OPAL in sediment.
!
        IF (Aout(idOPALbur,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALbur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOPALbur))
          Vinfo( 3)=Vname(3,idOPALbur)
          Vinfo(14)=Vname(4,idOPALbur)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOPALbur,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idOPALbur),                        &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef PRIMARY_PROD
!
!  Define Net primary production.
!
        IF (Aout(idNPP,ng)) THEN
          Vinfo( 1)=Vname(1,idNPP)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idNPP))
          Vinfo( 3)=Vname(3,idNPP)
          Vinfo(14)=Vname(4,idNPP)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idNPP,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idNPP),                        &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Aout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHsbl))
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idHsbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHsbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Aout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHbbl))
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idHbbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHbbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
!
!  Define 2D potential vorticity.
!
        IF (Aout(id2dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dPV))
          Vinfo( 3)=Vname(3,id2dPV)
          Vinfo(14)=Vname(4,id2dPV)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
# endif
          Vinfo(21)=Vname(6,id2dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id2dPV),  &
     &                   NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D relative vorticity.
!
        IF (Aout(id2dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dRV))
          Vinfo( 3)=Vname(3,id2dRV)
          Vinfo(14)=Vname(4,id2dRV)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
# endif
          Vinfo(21)=Vname(6,id2dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id2dRV),  &
     &                   NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D potential vorticity.
!
        IF (Aout(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dPV))
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#  endif
          Vinfo(21)=Vname(6,id3dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id3dPV),  &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D relative vorticity.
!
        IF (Aout(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dRV))
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#  endif
          Vinfo(21)=Vname(6,id3dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(id3dRV),  &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
# endif
!
!  Define quadratic <zeta*zeta> term.
!
        IF (Aout(idZZav,ng)) THEN
          Vinfo( 1)=Vname(1,idZZav)
          Vinfo( 2)=TRIM(Vname(2,idZZav))
          Vinfo( 3)=Vname(3,idZZav)
          Vinfo(14)=Vname(4,idZZav)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idZZav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idZZav),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <ubar*ubar> term.
!
        IF (Aout(idU2av,ng)) THEN
          Vinfo( 1)=Vname(1,idU2av)
          Vinfo( 2)=TRIM(Vname(2,idU2av))
          Vinfo( 3)=Vname(3,idU2av)
          Vinfo(14)=Vname(4,idU2av)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idU2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU2av),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <vbar*vbar> term.
!
        IF (Aout(idV2av,ng)) THEN
          Vinfo( 1)=Vname(1,idV2av)
          Vinfo( 2)=TRIM(Vname(2,idV2av))
          Vinfo( 3)=Vname(3,idV2av)
          Vinfo(14)=Vname(4,idV2av)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idV2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV2av),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define u-volume flux.
!
        IF (Aout(idHUav,ng)) THEN
          Vinfo( 1)=Vname(1,idHUav)
          Vinfo( 2)=TRIM(Vname(2,idHUav))
          Vinfo( 3)=Vname(3,idHUav)
          Vinfo(14)=Vname(4,idHUav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idHUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHUav),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define v-volume flux.
!
        IF (Aout(idHVav,ng)) THEN
          Vinfo( 1)=Vname(1,idHVav)
          Vinfo( 2)=TRIM(Vname(2,idHVav))
          Vinfo( 3)=Vname(3,idHVav)
          Vinfo(14)=Vname(4,idHVav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idHVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHVav),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*u> term.
!
        IF (Aout(idUUav,ng)) THEN
          Vinfo( 1)=Vname(1,idUUav)
          Vinfo( 2)=TRIM(Vname(2,idUUav))
          Vinfo( 3)=Vname(3,idUUav)
          Vinfo(14)=Vname(4,idUUav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUUav),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*v> term.
!
        IF (Aout(idUVav,ng)) THEN
          Vinfo( 1)=Vname(1,idUVav)
          Vinfo( 2)=TRIM(Vname(2,idUVav))
          Vinfo( 3)=Vname(3,idUVav)
          Vinfo(14)=Vname(4,idUVav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUVav),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <v*v> term.
!
        IF (Aout(idVVav,ng)) THEN
          Vinfo( 1)=Vname(1,idVVav)
          Vinfo( 2)=TRIM(Vname(2,idVVav))
          Vinfo( 3)=Vname(3,idVVav)
          Vinfo(14)=Vname(4,idVVav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVVav),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <t*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTTav(itrc)))
            Vinfo( 3)=Vname(3,idTTav(itrc))
            Vinfo(14)=Vname(4,idTTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(21)=Vname(6,idTTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTTav(itrc)), NF_FOUT,          &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define active tracers volume fluxes.
!
        DO itrc=1,NT(ng)
          IF (Aout(iHUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHUTav(itrc)))
            Vinfo( 3)=Vname(3,iHUTav(itrc))
            Vinfo(14)=Vname(4,iHUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#  endif
            Vinfo(21)=Vname(6,iHUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(iHUTav(itrc)), NF_FOUT,          &
     &                     nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (Aout(iHVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHVTav(itrc)))
            Vinfo( 3)=Vname(3,iHVTav(itrc))
            Vinfo(14)=Vname(4,iHVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#  endif
            Vinfo(21)=Vname(6,iHVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(iHVTav(itrc)), NF_FOUT,          &
     &                     nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define quadratic <u*t> and <v*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idUTav(itrc)))
            Vinfo( 3)=Vname(3,idUTav(itrc))
            Vinfo(14)=Vname(4,idUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#  endif
            Vinfo(21)=Vname(6,idUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUTav(itrc)), NF_FOUT,          &
     &                     nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (Aout(idVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idVTav(itrc)))
            Vinfo( 3)=Vname(3,idVTav(itrc))
            Vinfo(14)=Vname(4,idVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#  endif
            Vinfo(21)=Vname(6,idVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVTav(itrc)), NF_FOUT,          &
     &                     nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!
!  Define vertical viscosity coefficient.
!
        IF (Aout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvis))
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVvis),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Aout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTdif))
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idTdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
#  ifdef TS_VAR
!
!  Define vertical diffusion coefficient for salinity.
!
!        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idRdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idRdif))
          Vinfo( 3)=Vname(3,idRdif)
          Vinfo(14)=Vname(4,idRdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idRdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Aout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPair))
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idPair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES
!
!  Define surface air temperature.
!
        IF (Aout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTair))
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idTair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Aout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUair))
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVair))
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVair),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface wind at RHO-points.
!
        IF (Aout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUaiE))
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUaiE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVaiN))
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVaiN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define surface net heat flux.
!
        IF (Aout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idTsur(itemp))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTsur(itemp)), NF_FOUT,           &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SALINITY
!
!  Define surface net salt flux.
!
        IF (Aout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTsur(isalt))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTsur(isalt)), NF_FOUT,           &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#  if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Define latent heat flux.
!
        IF (Aout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLhea))
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
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Aout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idShea))
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
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Aout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLrad))
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
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Aout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idevap))
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
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idevap),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Aout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idrain),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
#  ifdef SHORTWAVE
!
!  Define net shortwave radiation flux.
!
        IF (Aout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSrad))
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
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
# endif
!
!  Define surface u-momentum stress.
!
        IF (Aout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUsms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Aout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVsms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Aout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Aout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_nf90 (ng, model, ldef, Aout, AVG,                  &
     &                     t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_nf90 (ng, model, ldef, Aout, AVG,                  &
     &                     t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_nf90 (ng, model, ldef, Aout, AVG,             &
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
        CALL wec_def_nf90 (ng, model, ldef, Aout, AVG,                  &
     &                     t2dgrd, u2dgrd, v2dgrd,                      &
     &                     t3dgrd, u3dgrd, v3dgrd, w3dgrd,              &
     &                     a3degrd, u3degrd, v3degrd, t3degrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, AVG(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVG(ng)%name
!
!  Open averages file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = AVG(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  average variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            AVG(ng)%Vid(idtime)=var_id(i)
# ifdef WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            AVG(ng)%Vid(idRwet)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            AVG(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            AVG(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            AVG(ng)%Vid(idVbar)=var_id(i)
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsuD))) THEN
            got_var(idFsuD)=.TRUE.
            AVG(ng)%Vid(idFsuD)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dD))) THEN
            got_var(idu2dD)=.TRUE.
            AVG(ng)%Vid(idu2dD)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dD))) THEN
            got_var(idv2dD)=.TRUE.
            AVG(ng)%Vid(idv2dD)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            AVG(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            AVG(ng)%Vid(idv2dN)=var_id(i)
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            AVG(ng)%Vid(idUfx1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            AVG(ng)%Vid(idUfx2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            AVG(ng)%Vid(idVfx1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            AVG(ng)%Vid(idVfx2)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            AVG(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            AVG(ng)%Vid(idVvel)=var_id(i)
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dD))) THEN
            got_var(idu3dD)=.TRUE.
            AVG(ng)%Vid(idu3dD)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dD))) THEN
            got_var(idv3dD)=.TRUE.
            AVG(ng)%Vid(idv3dD)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            AVG(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            AVG(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            AVG(ng)%Vid(idOvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            AVG(ng)%Vid(idWvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            AVG(ng)%Vid(idDano)=var_id(i)
#  ifdef NEMURO_SED1
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONsed))) THEN
            got_var(idPONsed)=.TRUE.
            AVG(ng)%Vid(idPONsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALsed))) THEN
            got_var(idOPALsed)=.TRUE.
            AVG(ng)%Vid(idOPALsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDENITsed))) THEN
            got_var(idDENITsed)=.TRUE.
            AVG(ng)%Vid(idDENITsed)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONbur))) THEN
            got_var(idPONbur)=.TRUE.
            AVG(ng)%Vid(idPONbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALbur))) THEN
            got_var(idOPALbur)=.TRUE.
            AVG(ng)%Vid(idOPALbur)=var_id(i)
#  endif
#  ifdef PRIMARY_PROD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNPP))) THEN
            got_var(idNPP)=.TRUE.
            AVG(ng)%Vid(idNPP)=var_id(i)
#  endif
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            AVG(ng)%Vid(idHsbl)=var_id(i)
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            AVG(ng)%Vid(idHbbl)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dPV))) THEN
            got_var(id2dPV)=.TRUE.
            AVG(ng)%Vid(id2dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dRV))) THEN
            got_var(id2dRV)=.TRUE.
            AVG(ng)%Vid(id2dRV)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            AVG(ng)%Vid(id3dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            AVG(ng)%Vid(id3dRV)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZZav))) THEN
            got_var(idZZav)=.TRUE.
            AVG(ng)%Vid(idZZav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2av))) THEN
            got_var(idU2av)=.TRUE.
            AVG(ng)%Vid(idU2av)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2av))) THEN
            got_var(idV2av)=.TRUE.
            AVG(ng)%Vid(idV2av)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHUav))) THEN
            got_var(idHUav)=.TRUE.
            AVG(ng)%Vid(idHUav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHVav))) THEN
            got_var(idHVav)=.TRUE.
            AVG(ng)%Vid(idHVav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUUav))) THEN
            got_var(idUUav)=.TRUE.
            AVG(ng)%Vid(idUUav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVav))) THEN
            got_var(idUVav)=.TRUE.
            AVG(ng)%Vid(idUVav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVVav))) THEN
            got_var(idVVav)=.TRUE.
            AVG(ng)%Vid(idVVav)=var_id(i)
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            AVG(ng)%Vid(idVvis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            AVG(ng)%Vid(idTdif)=var_id(i)
#   ifdef SALINITY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            AVG(ng)%Vid(idSdif)=var_id(i)
#   endif
#  endif
#  ifdef TS_VAR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRdif))) THEN
            got_var(idRdif)=.TRUE.
            AVG(ng)%Vid(idRdif)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            AVG(ng)%Vid(idPair)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            AVG(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            AVG(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            AVG(ng)%Vid(idUaiE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            AVG(ng)%Vid(idVaiN)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            AVG(ng)%Vid(idTsur(itemp))=var_id(i)
#  ifdef SALINITY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            AVG(ng)%Vid(idTsur(isalt))=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            AVG(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            AVG(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            AVG(ng)%Vid(idLrad)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            AVG(ng)%Vid(idTair)=var_id(i)
#   ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            AVG(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            AVG(ng)%Vid(idrain)=var_id(i)
#   endif
#  endif
#  ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            AVG(ng)%Vid(idSrad)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            AVG(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            AVG(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            AVG(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            AVG(ng)%Vid(idVbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            AVG(ng)%Vid(idUbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            AVG(ng)%Vid(idVbrs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            AVG(ng)%Vid(idUbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            AVG(ng)%Vid(idVbws)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            AVG(ng)%Vid(idUbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            AVG(ng)%Vid(idVbcs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            AVG(ng)%Vid(idUVwc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            AVG(ng)%Vid(idUbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            AVG(ng)%Vid(idVbot)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            AVG(ng)%Vid(idUbur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            AVG(ng)%Vid(idVbvr)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             AVG(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
!         DO itrc=1,NAT
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
              got_var(iHUTav(itrc))=.TRUE.
              AVG(ng)%Vid(iHUTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,iHVTav(itrc)))) THEN
              got_var(iHVTav(itrc))=.TRUE.
              AVG(ng)%Vid(iHVTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUTav(itrc)))) THEN
              got_var(idUTav(itrc))=.TRUE.
              AVG(ng)%Vid(idUTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVTav(itrc)))) THEN
              got_var(idVTav(itrc))=.TRUE.
              AVG(ng)%Vid(idVTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTTav(itrc)))) THEN
              got_var(idTTav(itrc))=.TRUE.
              AVG(ng)%Vid(idTTav(itrc))=var_id(i)
            END IF
          END DO
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTrcD(itrc)))) THEN
              got_var(idTrcD(itrc))=.TRUE.
              AVG(ng)%Vid(idTrcD(itrc))=var_id(i)
            END IF
          END DO
#  endif
# endif
        END DO
!
!  Check if averages variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Aout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Aout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Aout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idFsuD).and.Aout(idFsuD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsuD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dD).and.Aout(idu2dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dD).and.Aout(idv2dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idu2dE).and.Aout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Aout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
        IF (.not.got_var(idUfx1).and.Aout(idUfx1,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2).and.Aout(idUfx2,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1).and.Aout(idVfx1,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2).and.Aout(idVfx2,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idUvel).and.Aout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Aout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idu3dD).and.Aout(idu3dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dD).and.Aout(idv3dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idu3dE).and.Aout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Aout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Aout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Aout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Aout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef NEMURO_SED1
        IF (.not.got_var(idPONsed).and.Aout(idPONsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPONsed)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALsed).and.Aout(idOPALsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOPALsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDENITsed).and.Aout(idDENITsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDENITsed)),      &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPONbur).and.Aout(idPONbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPONbur)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALbur).and.Aout(idOPALbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOPALbur)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (.not.got_var(idNPP).and.Aout(idNPP,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idNPP)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Aout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Aout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(id2dPV).and.Aout(id2dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dRV).and.Aout(id2dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(id3dPV).and.Aout(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Aout(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idZZav).and.Aout(idZZav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZZav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2av).and.Aout(idU2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2av).and.Aout(idV2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idHUav).and.Aout(idHUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHVav).and.Aout(idHVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUUav).and.Aout(idUUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav).and.Aout(idUVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav).and.Aout(idVVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (.not.got_var(idVvis).and.Aout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Aout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SALINITY
        IF (.not.got_var(idSdif).and.Aout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef TS_VAR
        IF (.not.got_var(idRdif).and.Aout(idRdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idRdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Aout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Aout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Aout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Aout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Aout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idTsur(itemp)).and.Aout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SALINITY
        IF (.not.got_var(idTsur(isalt)).and.Aout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined FRC_COUPLING
        IF (.not.got_var(idLhea).and.Aout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Aout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Aout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Aout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef EMINUSP
        IF (.not.got_var(idevap).and.Aout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Aout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Aout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Aout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Aout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Aout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Aout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Aout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
!       DO itrc=1,NAT
        DO itrc=1,NT(ng)
          IF (.not.got_var(iHUTav(itrc)).and.Aout(iHUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc)).and.Aout(iHVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc)).and.Aout(idUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc)).and.Aout(idVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idTTav(itrc)).and.Aout(idTTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        DO itrc=1,NAT
          IF (.not.got_var(idTrcD(itrc)).and.Aout(idTrcD(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTrcD(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
# endif

# if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for time-averaged variables. Get
!  variable IDs.
!
        CALL bbl_def_nf90 (ng, model, ldef, Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined ICE_MODEL && defined SOLVE3D
!
!  Scan sea-ice variables from input NetCDF and activate switches for
!  time-averaged variables. Get variable IDs.
!
        CALL ice_def_nf90 (ng, model, ldef, Aout, AVG)
# endif

# if defined SEDIMENT && defined SOLVE3D
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for time-averaged variables. Get variable IDs.
!
        CALL sediment_def_nf90 (ng, model, ldef, Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

# if defined WEC_VF && defined SOLVE3D
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for time-averaged variables. Get variable IDs.
!
        CALL wec_def_nf90 (ng, model, ldef, Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nAVG(ng)) THEN
          IF (ndefAVG(ng).gt.0) THEN
            AVG(ng)%Rindex=((ntstart(ng)-1)-                            &
     &                      ndefAVG(ng)*((ntstart(ng)-1)/ndefAVG(ng)))/ &
     &                     nAVG(ng)
          ELSE
            AVG(ng)%Rindex=(ntstart(ng)-1)/nAVG(ng)
          END IF
        ELSE
          AVG(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial average time. Notice that the value is offset by half
!  nAVG*dt so there is not a special case when computing its value

# ifdef ADJOINT
!  in "ad_set_avg".
!
        AVGtime(ng)=time(ng)+0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
# else
!  in "set_avg".
!
      IF (ntsAVG(ng).eq.1) THEN
        AVGtime(ng)=time(ng)-0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      ELSE
        AVGtime(ng)=time(ng)+REAL(ntsAVG(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      END IF
# endif
!
  10  FORMAT (2x,'DEF_AVG_NF90     - creating average file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_AVG_NF90     - inquiring average file,',t56,      &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_AVG_NF90 - unable to create averages NetCDF',     &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG_NF90 - unable to open averages NetCDF',       &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_AVG_NF90 - unable to find variable: ',a,2x,       &
     &        ' in averages NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_avg_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE def_avg_pio (ng, model, ldef)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer :: i, ifield, itrc, j, nvd3, nvd4
#  if defined WRITE_WATER && defined MASKING
      integer :: xy_pdim, xyz_pdim
#  endif
      integer :: recdim, status

      integer :: DimIDs(nDimID)
      integer :: p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)

#  ifdef SOLVE3D
#   ifdef SEDIMENT
      integer :: b3dgrd(4)
#   endif
      integer :: p3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#  endif
!
      real(r8) :: Aval(6)
!
#  ifdef ADJOINT
      character (len=21)     :: Prefix
#  else
      character (len=13)     :: Prefix
#  endif
      character (len=256)    :: ncname
      character (len=MaxLen) :: Vinfo(Natt)

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", def_avg_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=AVG(ng)%name
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
!  Create a new averages file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL pio_netcdf_create (ng, model, TRIM(ncname),                &
     &                          AVG(ng)%pioFile)
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
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_rho',    &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_u',      &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_v',      &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xi_psi',    &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_rho',   &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_u',     &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_v',     &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'eta_psi',   &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xy_rho',    &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xy_u',      &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xy_v',      &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xy_psi',    &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef INWAVE_MODEL
        status=def_dim(ng, iNLM, AVG(ng)%pioFile, ncname, 'ND',         &
     &                 IOBOUNDS(ng)%ND, DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, iNLM, AVG(ng)%pioFile, ncname, 'NDp1',       &
     &                 IOBOUNDS(ng)%ND+1, DimIDs(37))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef SOLVE3D
#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xyz_rho',   &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xyz_u',     &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xyz_v',     &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xyz_w',     &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xyz_psi',   &
     &                 IOBOUNDS(ng)%xy_psi*N(ng), xyz_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 's_rho',     &
     &                 N(ng), DimIDs( 9))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 's_w',       &
     &                 N(ng)+1, DimIDs(10))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'tracer',    &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#   ifdef SEDIMENT
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'NST',       &
     &                 NST, DimIDs(32))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nbed',      &
     &                 Nbed, DimIDs(16))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#    if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'xybed',     &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#    endif
#   endif

#  if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#   ifdef ECOSIM
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nbands',    &
     &                 NBands, DimIDs(33))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nphy',      &
     &                 Nphy, DimIDs(25))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nbac',      &
     &                 Nbac, DimIDs(26))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Ndom',      &
     &                 Ndom, DimIDs(27))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nfec',      &
     &                 Nfec, DimIDs(28))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#   endif
#  endif

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'boundary',  &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#  ifdef FOUR_DVAR
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nstate',    &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  ifdef MULTI_SCALE_B
        status=def_dim(ng, model, AVG(ng)%pioFile, ncname, 'Nscale',    &
     &                 Nscale(ng), DimIDs(36))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

        status=def_dim(ng, model, AVG(ng)%pioFile, ncname,              &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 PIO_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

        recdim=DimIDs(12)
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
!
!  Define dimension vectors for staggered variables at PSI-points.
!
#  if defined WRITE_WATER && defined MASKING
        p2dgrd(1)=xy_pdim
        p2dgrd(2)=DimIDs(12)
#   ifdef SOLVE3D
        p3dgrd(1)=xyz_pdim
        p3dgrd(2)=DimIDs(12)
#   endif
#  else
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
#   ifdef SOLVE3D
        p3dgrd(1)=DimIDs( 4)
        p3dgrd(2)=DimIDs( 8)
        p3dgrd(3)=DimIDs( 9)
        p3dgrd(4)=DimIDs(12)
#   endif
#  endif
#  if defined WEC_VF && defined SOLVE3D
!
!  Define dimension vector for energy wave layer type variables.
!
        a3degrd(1)=DimIDs( 1)
        a3degrd(2)=DimIDs( 5)
        a3degrd(3)=DimIDs(37)
        a3degrd(4)=DimIDs(12)

        u3degrd(1)=DimIDs( 2)
        u3degrd(2)=DimIDs( 6)
        u3degrd(3)=DimIDs(37)
        u3degrd(4)=DimIDs(12)

        v3degrd(1)=DimIDs( 3)
        v3degrd(2)=DimIDs( 7)
        v3degrd(3)=DimIDs(37)
        v3degrd(4)=DimIDs(12)

        t3degrd(1)=DimIDs( 1)
        t3degrd(2)=DimIDs( 5)
        t3degrd(3)=DimIDs(38)
        t3degrd(4)=DimIDs(12)
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
        AVG(ng)%Rindex=0
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
!  Set long name prefix string.
!
#  ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
#  else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
#  endif
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, AVG(ng)%pioFile, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,a)') 'averaged ', TRIM(Vname(2,idtime))
        WRITE (Vinfo( 3),'(a,a)') 'seconds since ', TRIM(Rclock%string)
        Vinfo( 4)=TRIM(Rclock%calendar)
        Vinfo(14)=Vname(4,idtime)
        Vinfo(21)=Vname(6,idtime)
        AVG(ng)%pioVar(idtime)%dkind=PIO_TOUT
        AVG(ng)%pioVar(idtime)%gtype=0
!
        status=def_var(ng, model, AVG(ng)%pioFile,                      &
     &                 AVG(ng)%pioVar(idtime)%vd,                       &
     &                 PIO_TOUT, 1, (/recdim/), Aval, Vinfo, ncname,    &
     &                 SetParAccess = .TRUE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Define free-surface.
!
        IF (Aout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
#  if !defined WET_DRY && (defined WRITE_WATER && defined MASKING)
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idFsur)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          AVG(ng)%pioVar(idFsur)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idFsur)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idFsur)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided free-surface.
!
        IF (Aout(idFsuD,ng)) THEN
          Vinfo( 1)=Vname(1,idFsuD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsuD))
          Vinfo( 3)=Vname(3,idFsuD)
          Vinfo(14)=Vname(4,idFsuD)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idFsuD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsuD,ng),r8)
          AVG(ng)%pioVar(idFsuD)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idFsuD)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idFsuD)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D momentum in the XI-direction.
!
        IF (Aout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          AVG(ng)%pioVar(idUbar)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUbar)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUbar)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 2D momentum in the XI-direction.
!
        IF (Aout(idu2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dD))
          Vinfo( 3)=Vname(3,idu2dD)
          Vinfo(14)=Vname(4,idu2dD)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idu2dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dD,ng),r8)
          AVG(ng)%pioVar(idu2dD)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idu2dD)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idu2dD)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Aout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVbar)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          AVG(ng)%pioVar(idVbar)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVbar)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVbar)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 2D momentum in the ETA-direction.
!
        IF (Aout(idv2dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dD))
          Vinfo( 3)=Vname(3,idv2dD)
          Vinfo(14)=Vname(4,idv2dD)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idv2dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dD,ng),r8)
          AVG(ng)%pioVar(idv2dD)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idv2dD)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idv2dD)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Aout(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dE))
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idu2dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          AVG(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idu2dE)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idu2dE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Aout(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dN))
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idv2dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          AVG(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idv2dN)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idv2dN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
#   ifdef FORWARD_WRITE
!
!  Define time-averaged mass fluxes for 3D momentum coupling.
!
        IF (Aout(idUfx1,ng)) THEN
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
          Aval(5)=REAL(u2dvar,r8)
          AVG(ng)%pioVar(idUfx1)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUfx1)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUfx1)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idUfx2,ng)) THEN
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
          Aval(5)=REAL(u2dvar,r8)
          AVG(ng)%pioVar(idUfx2)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUfx2)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUfx2)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVfx1,ng)) THEN
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
          Aval(5)=REAL(v2dvar,r8)
          AVG(ng)%pioVar(idVfx1)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVfx1)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVfx1)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVfx2,ng)) THEN
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
          Aval(5)=REAL(v2dvar,r8)
          AVG(ng)%pioVar(idVfx2)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVfx2)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVfx2)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Aout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          AVG(ng)%pioVar(idUvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUvel)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUvel)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 3D momentum component in the XI-direction.
!
        IF (Aout(idu3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dD))
          Vinfo( 3)=Vname(3,idu3dD)
          Vinfo(14)=Vname(4,idu3dD)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(21)=Vname(6,idu3dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dD,ng),r8)
          AVG(ng)%pioVar(idu3dD)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idu3dD)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idu3dD)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Aout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          AVG(ng)%pioVar(idVvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVvel)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVvel)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 3D momentum component in the ETA-direction.
!
        IF (Aout(idv3dD,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dD)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dD))
          Vinfo( 3)=Vname(3,idv3dD)
          Vinfo(14)=Vname(4,idv3dD)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(21)=Vname(6,idv3dD)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dD,ng),r8)
          AVG(ng)%pioVar(idv3dD)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idv3dD)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idv3dD)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Aout(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dE))
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idu3dE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          AVG(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idu3dE)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idu3dE)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Aout(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dN))
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idv3dN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          AVG(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idv3dN)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idv3dN)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Aout(idOvel,ng)) THEN
          Vinfo( 1)=Vname(1,idOvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOvel))
          Vinfo( 3)=Vname(3,idOvel)
          Vinfo(14)=Vname(4,idOvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idOvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          AVG(ng)%pioVar(idOvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idOvel)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idOvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define "true" vertical momentum component.
!
        IF (Aout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvel))
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idWvel)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          AVG(ng)%pioVar(idWvel)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idWvel)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idWvel)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#   ifdef SEDIMENT
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#   endif
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTvar(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            AVG(ng)%pioTrc(itrc)%dkind=PIO_FOUT
            AVG(ng)%pioTrc(itrc)%gtype=r3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioTrc(itrc)%vd,                     &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided temperature and salinity.
!
        DO itrc=1,NAT
          IF (Aout(idTrcD(itrc),ng)) THEN
            ifield=idTrcD(itrc)
            Vinfo( 1)=Vname(1,ifield)
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,ifield))
            Vinfo( 3)=Vname(3,ifield)
            Vinfo(14)=Vname(4,ifield)
            Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#    endif
            Vinfo(21)=Vname(6,ifield)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            AVG(ng)%pioVar(ifield)%dkind=PIO_FOUT
            AVG(ng)%pioVar(ifield)%gtype=r3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(ifield)%vd,                   &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
#   endif
!
!  Define density anomaly.
!
        IF (Aout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idDano)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          AVG(ng)%pioVar(idDano)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idDano)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idDano)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  ifdef NEMURO_SED1
!
!  Define PON in sediment.
!
        IF (Aout(idPONsed,ng)) THEN
          Vinfo( 1)=Vname(1,idPONsed)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPONsed))
          Vinfo( 3)=Vname(3,idPONsed)
          Vinfo(14)=Vname(4,idPONsed)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPONsed,ng),r8)
          AVG(ng)%pioVar(idPONsed)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idPONsed)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idPONsed)%vd,                   &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define OPAL in sediment.
!
        IF (Aout(idOPALsed,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALsed)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOPALsed))
          Vinfo( 3)=Vname(3,idOPALsed)
          Vinfo(14)=Vname(4,idOPALsed)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOPALsed,ng),r8)
          AVG(ng)%pioVar(idOPALsed)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idOPALsed)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idOPALsed)%vd,                  &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define DENIT in sediment.
!
        IF (Aout(idDENITsed,ng)) THEN
          Vinfo( 1)=Vname(1,idDENITsed)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDENITsed))
          Vinfo( 3)=Vname(3,idDENITsed)
          Vinfo(14)=Vname(4,idDENITsed)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDENITsed,ng),r8)
          AVG(ng)%pioVar(idDENITsed)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idDENITsed)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idDENITsed)%vd,                 &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define buried PON in sediment.
!
        IF (Aout(idPONbur,ng)) THEN
          Vinfo( 1)=Vname(1,idPONbur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPONbur))
          Vinfo( 3)=Vname(3,idPONbur)
          Vinfo(14)=Vname(4,idPONbur)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPONbur,ng),r8)
          AVG(ng)%pioVar(idPONbur)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idPONbur)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idPONbur)%vd,                   &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define buried OPAL in sediment.
!
        IF (Aout(idOPALbur,ng)) THEN
          Vinfo( 1)=Vname(1,idOPALbur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idOPALbur))
          Vinfo( 3)=Vname(3,idOPALbur)
          Vinfo(14)=Vname(4,idOPALbur)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOPALbur,ng),r8)
          AVG(ng)%pioVar(idOPALbur)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idOPALbur)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idOPALbur)%vd,                  &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
#  ifdef PRIMARY_PROD
!
!  Define Net primary production.
!
        IF (Aout(idNPP,ng)) THEN
          Vinfo( 1)=Vname(1,idNPP)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idNPP))
          Vinfo( 3)=Vname(3,idNPP)
          Vinfo(14)=Vname(4,idNPP)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idNPP,ng),r8)
          AVG(ng)%pioVar(idNPP)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idNPP)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idNPP)%vd,                      &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif

#   ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Aout(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHsbl))
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idHsbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          AVG(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idHsbl)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idHsbl)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Aout(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHbbl))
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idHbbl)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          AVG(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idHbbl)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idHbbl)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
!
!  Define 2D potential vorticity.
!
        IF (Aout(id2dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dPV))
          Vinfo( 3)=Vname(3,id2dPV)
          Vinfo(14)=Vname(4,id2dPV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#  endif
          Vinfo(21)=Vname(6,id2dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          AVG(ng)%pioVar(id2dPV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id2dPV)%gtype=p2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id2dPV)%vd,                     &
     &                   PIO_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 2D relative vorticity.
!
        IF (Aout(id2dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dRV))
          Vinfo( 3)=Vname(3,id2dRV)
          Vinfo(14)=Vname(4,id2dRV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#  endif
          Vinfo(21)=Vname(6,id2dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          AVG(ng)%pioVar(id2dRV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id2dRV)%gtype=p2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id2dRV)%vd,                     &
     &                   PIO_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define 3D potential vorticity.
!
        IF (Aout(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dPV))
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#   endif
          Vinfo(21)=Vname(6,id3dPV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          AVG(ng)%pioVar(id3dPV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id3dPV)%gtype=p3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id3dPV)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define 3D relative vorticity.
!
        IF (Aout(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dRV))
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#   endif
          Vinfo(21)=Vname(6,id3dRV)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          AVG(ng)%pioVar(id3dRV)%dkind=PIO_FOUT
          AVG(ng)%pioVar(id3dRV)%gtype=p3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(id3dRV)%vd,                     &
     &                   PIO_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#  endif
!
!  Define quadratic <zeta*zeta> term.
!
        IF (Aout(idZZav,ng)) THEN
          Vinfo( 1)=Vname(1,idZZav)
          Vinfo( 2)=TRIM(Vname(2,idZZav))
          Vinfo( 3)=Vname(3,idZZav)
          Vinfo(14)=Vname(4,idZZav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idZZav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          AVG(ng)%pioVar(idZZav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idZZav)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idZZav)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <ubar*ubar> term.
!
        IF (Aout(idU2av,ng)) THEN
          Vinfo( 1)=Vname(1,idU2av)
          Vinfo( 2)=TRIM(Vname(2,idU2av))
          Vinfo( 3)=Vname(3,idU2av)
          Vinfo(14)=Vname(4,idU2av)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idU2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          AVG(ng)%pioVar(idU2av)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idU2av)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idU2av)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <vbar*vbar> term.
!
        IF (Aout(idV2av,ng)) THEN
          Vinfo( 1)=Vname(1,idV2av)
          Vinfo( 2)=TRIM(Vname(2,idV2av))
          Vinfo( 3)=Vname(3,idV2av)
          Vinfo(14)=Vname(4,idV2av)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idV2av)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          AVG(ng)%pioVar(idV2av)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idV2av)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idV2av)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  ifdef SOLVE3D
!
!  Define u-volume flux.
!
        IF (Aout(idHUav,ng)) THEN
          Vinfo( 1)=Vname(1,idHUav)
          Vinfo( 2)=TRIM(Vname(2,idHUav))
          Vinfo( 3)=Vname(3,idHUav)
          Vinfo(14)=Vname(4,idHUav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idHUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          AVG(ng)%pioVar(idHUav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idHUav)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idHUav)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define v-volume flux.
!
        IF (Aout(idHVav,ng)) THEN
          Vinfo( 1)=Vname(1,idHVav)
          Vinfo( 2)=TRIM(Vname(2,idHVav))
          Vinfo( 3)=Vname(3,idHVav)
          Vinfo(14)=Vname(4,idHVav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idHVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          AVG(ng)%pioVar(idHVav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idHVav)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idHVav)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*u> term.
!
        IF (Aout(idUUav,ng)) THEN
          Vinfo( 1)=Vname(1,idUUav)
          Vinfo( 2)=TRIM(Vname(2,idUUav))
          Vinfo( 3)=Vname(3,idUUav)
          Vinfo(14)=Vname(4,idUUav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(21)=Vname(6,idUUav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          AVG(ng)%pioVar(idUUav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUUav)%gtype=u3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUUav)%vd,                     &
     &                   PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <u*v> term.
!
        IF (Aout(idUVav,ng)) THEN
          Vinfo( 1)=Vname(1,idUVav)
          Vinfo( 2)=TRIM(Vname(2,idUVav))
          Vinfo( 3)=Vname(3,idUVav)
          Vinfo(14)=Vname(4,idUVav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idUVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          AVG(ng)%pioVar(idUVav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUVav)%gtype=r3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUVav)%vd,                     &
     &                   PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <v*v> term.
!
        IF (Aout(idVVav,ng)) THEN
          Vinfo( 1)=Vname(1,idVVav)
          Vinfo( 2)=TRIM(Vname(2,idVVav))
          Vinfo( 3)=Vname(3,idVVav)
          Vinfo(14)=Vname(4,idVVav)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(21)=Vname(6,idVVav)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          AVG(ng)%pioVar(idVVav)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVVav)%gtype=v3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVVav)%vd,                     &
     &                   PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define quadratic <t*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idTTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idTTav(itrc)))
            Vinfo( 3)=Vname(3,idTTav(itrc))
            Vinfo(14)=Vname(4,idTTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#   endif
            Vinfo(21)=Vname(6,idTTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            AVG(ng)%pioVar(idTTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTTav(itrc))%gtype=r3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idTTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define active tracers volume fluxes.
!
        DO itrc=1,NT(ng)
          IF (Aout(iHUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHUTav(itrc)))
            Vinfo( 3)=Vname(3,iHUTav(itrc))
            Vinfo(14)=Vname(4,iHUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(21)=Vname(6,iHUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            AVG(ng)%pioVar(iHUTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(iHUTav(itrc))%gtype=u3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(iHUTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (Aout(iHVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,iHVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,iHVTav(itrc)))
            Vinfo( 3)=Vname(3,iHVTav(itrc))
            Vinfo(14)=Vname(4,iHVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(21)=Vname(6,iHVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            AVG(ng)%pioVar(iHVTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(iHVTav(itrc))%gtype=v3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(iHVTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define quadratic <u*t> and <v*t> terms.
!
        DO itrc=1,NT(ng)
          IF (Aout(idUTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idUTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idUTav(itrc)))
            Vinfo( 3)=Vname(3,idUTav(itrc))
            Vinfo(14)=Vname(4,idUTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#   endif
            Vinfo(21)=Vname(6,idUTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            AVG(ng)%pioVar(idUTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUTav(itrc))%gtype=u3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idUTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          IF (Aout(idVTav(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idVTav(itrc))
            Vinfo( 2)=TRIM(Vname(2,idVTav(itrc)))
            Vinfo( 3)=Vname(3,idVTav(itrc))
            Vinfo(14)=Vname(4,idVTav(itrc))
            Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#   endif
            Vinfo(21)=Vname(6,idVTav(itrc))
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            AVG(ng)%pioVar(idVTav(itrc))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVTav(itrc))%gtype=v3dvar
!
            status=def_var(ng, model, AVG(ng)%pioFile,                  &
     &                     AVG(ng)%pioVar(idVTav(itrc))%vd,             &
     &                     PIO_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

#   if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
!
!  Define vertical viscosity coefficient.
!
        IF (Aout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvis))
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idVvis)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          AVG(ng)%pioVar(idVvis)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVvis)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVvis)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Aout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTdif))
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idTdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          AVG(ng)%pioVar(idTdif)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTdif)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#    ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(21)=Vname(6,idSdif)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          AVG(ng)%pioVar(idSdif)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idSdif)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idSdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#    endif
#   endif
#  ifdef TS_VAR
!
!  Define vertical diffusion coefficient for salinity.
!
!        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idRdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idRdif))
          Vinfo( 3)=Vname(3,idRdif)
          Vinfo(14)=Vname(4,idRdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRdif,ng),r8)
          AVG(ng)%pioVar(idRdif)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idRdif)%gtype=w3dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idRdif)%vd,                     &
     &                   PIO_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,   &
     &                   SetFillVal = .FALSE.)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!        END IF
#  endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Aout(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPair))
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idPair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          AVG(ng)%pioVar(idPair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idPair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idPair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined BULK_FLUXES
!
!  Define surface air temperature.
!
        IF (Aout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTair))
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          AVG(ng)%pioVar(idTair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Aout(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUair))
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idUair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          AVG(ng)%pioVar(idUair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVair))
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idVair)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          AVG(ng)%pioVar(idVair)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVair)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVair)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define Eastward/Northward surface wind at RHO-points.
!
        IF (Aout(idUaiE,ng)) THEN
          Vinfo( 1)=Vname(1,idUaiE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUaiE))
          Vinfo( 3)=Vname(3,idUaiE)
          Vinfo(14)=Vname(4,idUaiE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idUaiE)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUaiE,ng),r8)
          AVG(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUaiE)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUaiE)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        IF (Aout(idVaiN,ng)) THEN
          Vinfo( 1)=Vname(1,idVaiN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVaiN))
          Vinfo( 3)=Vname(3,idVaiN)
          Vinfo(14)=Vname(4,idVaiN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)=Vname(6,idVaiN)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVaiN,ng),r8)
          AVG(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVaiN)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVaiN)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
!
!  Define surface net heat flux.
!
        IF (Aout(idTsur(itemp),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(itemp))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(itemp)))
          Vinfo( 3)=Vname(3,idTsur(itemp))
          Vinfo(11)='upward flux, cooling'
          Vinfo(12)='downward flux, heating'
          Vinfo(14)=Vname(4,idTsur(itemp))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(21)=Vname(6,idTsur(itemp))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          AVG(ng)%pioVar(idTsur(itemp))%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTsur(itemp))%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTsur(itemp))%vd,              &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#   ifdef SALINITY
!
!  Define surface net salt flux.
!
        IF (Aout(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(21)=Vname(6,idTsur(isalt))
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          AVG(ng)%pioVar(idTsur(isalt))%dkind=PIO_FOUT
          AVG(ng)%pioVar(idTsur(isalt))%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idTsur(isalt))%vd,              &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   if defined BULK_FLUXES || defined FRC_COUPLING
!
!  Define latent heat flux.
!
        IF (Aout(idLhea,ng)) THEN
          Vinfo( 1)=Vname(1,idLhea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLhea))
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
          AVG(ng)%pioVar(idLhea)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idLhea)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idLhea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Aout(idShea,ng)) THEN
          Vinfo( 1)=Vname(1,idShea)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idShea))
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
          AVG(ng)%pioVar(idShea)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idShea)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idShea)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define net longwave radiation flux.
!
        IF (Aout(idLrad,ng)) THEN
          Vinfo( 1)=Vname(1,idLrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLrad))
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
          AVG(ng)%pioVar(idLrad)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idLrad)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idLrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif

#   ifdef BULK_FLUXES
#    ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Aout(idevap,ng)) THEN
          Vinfo( 1)=Vname(1,idevap)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idevap))
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
          AVG(ng)%pioVar(idevap)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idevap)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idevap)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Aout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#     if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#     endif
          Vinfo(21)=Vname(6,idrain)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          AVG(ng)%pioVar(idrain)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idrain)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idrain)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#    endif
#   endif
#   ifdef SHORTWAVE
!
!  Define net shortwave radiation flux.
!
        IF (Aout(idSrad,ng)) THEN
          Vinfo( 1)=Vname(1,idSrad)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSrad))
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
          AVG(ng)%pioVar(idSrad)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idSrad)%gtype=r2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idSrad)%vd,                     &
     &                   PIO_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
#   endif
#  endif
!
!  Define surface u-momentum stress.
!
        IF (Aout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          AVG(ng)%pioVar(idUsms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUsms)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUsms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Aout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVsms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          AVG(ng)%pioVar(idVsms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVsms)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVsms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Aout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(21)=Vname(6,idUbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          AVG(ng)%pioVar(idUbms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idUbms)%gtype=u2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idUbms)%vd,                     &
     &                   PIO_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Aout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(21)=Vname(6,idVbms)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          AVG(ng)%pioVar(idVbms)%dkind=PIO_FOUT
          AVG(ng)%pioVar(idVbms)%gtype=v2dvar
!
          status=def_var(ng, model, AVG(ng)%pioFile,                    &
     &                   AVG(ng)%pioVar(idVbms)%vd,                     &
     &                   PIO_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF

#  if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the bottom boundary layer model or waves variables.
!-----------------------------------------------------------------------
!
        CALL bbl_def_pio (ng, model, ldef, Aout, AVG,                   &
     &                    t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ICE_MODEL && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sea-ice model variables.
!-----------------------------------------------------------------------
!
        CALL ice_def_pio (ng, model, ldef, Aout, AVG,                   &
     &                    t2dgrd, u2dgrd, v2dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined SEDIMENT && defined SOLVE3D
!
!-----------------------------------------------------------------------
!  Define the sediment model variables.
!-----------------------------------------------------------------------
!
        CALL sediment_def_pio (ng, model, ldef, Aout, AVG,              &
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
        CALL wec_def_pio (ng, model, ldef, Aout, AVG,                   &
     &                    t2dgrd, u2dgrd, v2dgrd,                       &
     &                    t3dgrd, u3dgrd, v3dgrd, w3dgrd)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL pio_netcdf_enddef (ng, model, ncname, AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, AVG(ng)%pioFile, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVG(ng)%name
!
!  Open averages file for read/write.
!
        CALL pio_netcdf_open (ng, model, ncname, 1, AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL pio_netcdf_check_dim (ng, model, ncname,                   &
     &                             pioFile = AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about the variables.
!
        CALL pio_netcdf_inq_var (ng, model, ncname,                     &
     &                           pioFile = AVG(ng)%pioFile)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  average variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            AVG(ng)%pioVar(idtime)%vd=var_desc(i)
            AVG(ng)%pioVar(idtime)%dkind=PIO_TOUT
            AVG(ng)%pioVar(idtime)%gtype=0
#  ifdef WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            AVG(ng)%pioVar(idRwet)%vd=var_desc(i)
            AVG(ng)%pioVar(idRwet)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idRwet)%gtype=r2dvar
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            AVG(ng)%pioVar(idFsur)%vd=var_desc(i)
            AVG(ng)%pioVar(idFsur)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idFsur)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            AVG(ng)%pioVar(idUbar)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbar)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbar)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            AVG(ng)%pioVar(idVbar)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbar)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbar)%gtype=v2dvar
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsuD))) THEN
            got_var(idFsuD)=.TRUE.
            AVG(ng)%pioVar(idFsuD)%vd=var_desc(i)
            AVG(ng)%pioVar(idFsuD)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idFsuD)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dD))) THEN
            got_var(idu2dD)=.TRUE.
            AVG(ng)%pioVar(idu2dD)%vd=var_desc(i)
            AVG(ng)%pioVar(idu2dD)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idu2dD)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dD))) THEN
            got_var(idv2dD)=.TRUE.
            AVG(ng)%pioVar(idv2dD)%vd=var_desc(i)
            AVG(ng)%pioVar(idv2dD)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idv2dD)%gtype=v2dvar
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            AVG(ng)%pioVar(idu2dE)%vd=var_desc(i)
            AVG(ng)%pioVar(idu2dE)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idu2dE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            AVG(ng)%pioVar(idv2dN)%vd=var_desc(i)
            AVG(ng)%pioVar(idv2dN)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idv2dN)%gtype=r2dvar
#  ifdef SOLVE3D
#   ifdef FORWARD_WRITE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            AVG(ng)%pioVar(idUfx1)%vd=var_desc(i)
            AVG(ng)%pioVar(idUfx1)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUfx1)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            AVG(ng)%pioVar(idUfx2)%vd=var_desc(i)
            AVG(ng)%pioVar(idUfx2)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUfx2)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            AVG(ng)%pioVar(idVfx1)%vd=var_desc(i)
            AVG(ng)%pioVar(idVfx1)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVfx1)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            AVG(ng)%pioVar(idVfx2)%vd=var_desc(i)
            AVG(ng)%pioVar(idVfx2)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVfx2)%gtype=v2dvar
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            AVG(ng)%pioVar(idUvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idUvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUvel)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            AVG(ng)%pioVar(idVvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idVvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVvel)%gtype=v3dvar
#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dD))) THEN
            got_var(idu3dD)=.TRUE.
            AVG(ng)%pioVar(idu3dD)%vd=var_desc(i)
            AVG(ng)%pioVar(idu3dD)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idu3dD)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dD))) THEN
            got_var(idv3dD)=.TRUE.
            AVG(ng)%pioVar(idv3dD)%vd=var_desc(i)
            AVG(ng)%pioVar(idv3dD)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idv3dD)%gtype=v3dvar
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            AVG(ng)%pioVar(idu3dE)%vd=var_desc(i)
            AVG(ng)%pioVar(idu3dE)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idu3dE)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            AVG(ng)%pioVar(idv3dN)%vd=var_desc(i)
            AVG(ng)%pioVar(idv3dN)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idv3dN)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            AVG(ng)%pioVar(idOvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idOvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idOvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            AVG(ng)%pioVar(idWvel)%vd=var_desc(i)
            AVG(ng)%pioVar(idWvel)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idWvel)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            AVG(ng)%pioVar(idDano)%vd=var_desc(i)
            AVG(ng)%pioVar(idDano)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idDano)%gtype=r3dvar
#  ifdef NEMURO_SED1
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONsed))) THEN
            got_var(idPONsed)=.TRUE.
            AVG(ng)%pioVar(idPONsed)%vd=var_desc(i)
            AVG(ng)%pioVar(idPONsed)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idPONsed)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALsed))) THEN
            got_var(idOPALsed)=.TRUE.
            AVG(ng)%pioVar(idOPALsed)%vd=var_desc(i)
            AVG(ng)%pioVar(idOPALsed)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idOPALsed)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDENITsed))) THEN
            got_var(idDENITsed)=.TRUE.
            AVG(ng)%pioVar(idDENITsed)%vd=var_desc(i)
            AVG(ng)%pioVar(idDENITsed)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idDENITsed)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPONbur))) THEN
            got_var(idPONbur)=.TRUE.
            AVG(ng)%pioVar(idPONbur)%vd=var_desc(i)
            AVG(ng)%pioVar(idPONbur)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idPONbur)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOPALbur))) THEN
            got_var(idOPALbur)=.TRUE.
            AVG(ng)%pioVar(idOPALbur)%vd=var_desc(i)
            AVG(ng)%pioVar(idOPALbur)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idOPALbur)%gtype=r3dvar
#  endif
#  ifdef PRIMARY_PROD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idNPP))) THEN
            got_var(idNPP)=.TRUE.
            AVG(ng)%pioVar(idNPP)%vd=var_desc(i)
            AVG(ng)%pioVar(idNPP)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idNPP)%gtype=r3dvar
#  endif
#   ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            AVG(ng)%pioVar(idHsbl)%vd=var_desc(i)
            AVG(ng)%pioVar(idHsbl)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHsbl)%gtype=r2dvar
#   endif
#   ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            AVG(ng)%pioVar(idHbbl)%vd=var_desc(i)
            AVG(ng)%pioVar(idHbbl)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHbbl)%gtype=r2dvar
#   endif
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dPV))) THEN
            got_var(id2dPV)=.TRUE.
            AVG(ng)%pioVar(id2dPV)%vd=var_desc(i)
            AVG(ng)%pioVar(id2dPV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id2dPV)%gtype=p2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dRV))) THEN
            got_var(id2dRV)=.TRUE.
            AVG(ng)%pioVar(id2dRV)%vd=var_desc(i)
            AVG(ng)%pioVar(id2dRV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id2dRV)%gtype=p2dvar
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            AVG(ng)%pioVar(id3dPV)%vd=var_desc(i)
            AVG(ng)%pioVar(id3dPV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id3dPV)%gtype=p3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            AVG(ng)%pioVar(id3dRV)%vd=var_desc(i)
            AVG(ng)%pioVar(id3dRV)%dkind=PIO_FOUT
            AVG(ng)%pioVar(id3dRV)%gtype=p3dvar
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZZav))) THEN
            got_var(idZZav)=.TRUE.
            AVG(ng)%pioVar(idZZav)%vd=var_desc(i)
            AVG(ng)%pioVar(idZZav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idZZav)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2av))) THEN
            got_var(idU2av)=.TRUE.
            AVG(ng)%pioVar(idU2av)%vd=var_desc(i)
            AVG(ng)%pioVar(idU2av)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idU2av)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2av))) THEN
            got_var(idV2av)=.TRUE.
            AVG(ng)%pioVar(idV2av)%vd=var_desc(i)
            AVG(ng)%pioVar(idV2av)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idV2av)%gtype=v2dvar
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHUav))) THEN
            got_var(idHUav)=.TRUE.
            AVG(ng)%pioVar(idHUav)%vd=var_desc(i)
            AVG(ng)%pioVar(idHUav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHUav)%gtype=u3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHVav))) THEN
            got_var(idHVav)=.TRUE.
            AVG(ng)%pioVar(idHVav)%vd=var_desc(i)
            AVG(ng)%pioVar(idHVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHVav)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUUav))) THEN
            got_var(idUUav)=.TRUE.
            AVG(ng)%pioVar(idUUav)%vd=var_desc(i)
            AVG(ng)%pioVar(idHVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idHVav)%gtype=v3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVav))) THEN
            got_var(idUVav)=.TRUE.
            AVG(ng)%pioVar(idUVav)%vd=var_desc(i)
            AVG(ng)%pioVar(idUVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUVav)%gtype=r3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVVav))) THEN
            got_var(idVVav)=.TRUE.
            AVG(ng)%pioVar(idVVav)%vd=var_desc(i)
            AVG(ng)%pioVar(idVVav)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVVav)%gtype=v3dvar
#   if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            AVG(ng)%pioVar(idVvis)%vd=var_desc(i)
            AVG(ng)%pioVar(idVvis)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVvis)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            AVG(ng)%pioVar(idTdif)%vd=var_desc(i)
            AVG(ng)%pioVar(idTdif)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTdif)%gtype=w3dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            AVG(ng)%pioVar(idSdif)%vd=var_desc(i)
            AVG(ng)%pioVar(idSdif)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idSdif)%gtype=w3dvar
#   endif
#  ifdef TS_VAR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRdif))) THEN
            got_var(idRdif)=.TRUE.
            AVG(ng)%pioVar(idRdif)%vd=var_desc(i)
            AVG(ng)%pioVar(idRdif)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idRdif)%gtype=w3dvar
#  endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            AVG(ng)%pioVar(idPair)%vd=var_desc(i)
            AVG(ng)%pioVar(idPair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idPair)%gtype=r2dvar
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            AVG(ng)%pioVar(idUair)%vd=var_desc(i)
            AVG(ng)%pioVar(idUair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            AVG(ng)%pioVar(idVair)%vd=var_desc(i)
            AVG(ng)%pioVar(idVair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVair)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUaiE))) THEN
            got_var(idUaiE)=.TRUE.
            AVG(ng)%pioVar(idUaiE)%vd=var_desc(i)
            AVG(ng)%pioVar(idUaiE)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUaiE)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVaiN))) THEN
            got_var(idVaiN)=.TRUE.
            AVG(ng)%pioVar(idVaiN)%vd=var_desc(i)
            AVG(ng)%pioVar(idVaiN)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVaiN)%gtype=r2dvar
#   endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            AVG(ng)%pioVar(idTsur(itemp))%vd=var_desc(i)
            AVG(ng)%pioVar(idTsur(itemp))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTsur(itemp))%gtype=r2dvar
#   ifdef SALINITY
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            AVG(ng)%pioVar(idTsur(isalt))%vd=var_desc(i)
            AVG(ng)%pioVar(idTsur(isalt))%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTsur(isalt))%gtype=r2dvar
#   endif
#   if defined BULK_FLUXES || defined FRC_COUPLING
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            AVG(ng)%pioVar(idLhea)%vd=var_desc(i)
            AVG(ng)%pioVar(idLhea)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idLhea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            AVG(ng)%pioVar(idShea)%vd=var_desc(i)
            AVG(ng)%pioVar(idShea)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idShea)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            AVG(ng)%pioVar(idLrad)%vd=var_desc(i)
            AVG(ng)%pioVar(idLrad)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idLrad)%gtype=r2dvar
#   endif
#   ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            AVG(ng)%pioVar(idTair)%vd=var_desc(i)
            AVG(ng)%pioVar(idTair)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idTair)%gtype=r2dvar
#    ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            AVG(ng)%pioVar(idevap)%vd=var_desc(i)
            AVG(ng)%pioVar(idevap)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idevap)%gtype=r2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            AVG(ng)%pioVar(idrain)%vd=var_desc(i)
            AVG(ng)%pioVar(idrain)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idrain)%gtype=r2dvar
#    endif
#   endif
#   ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            AVG(ng)%pioVar(idSrad)%vd=var_desc(i)
            AVG(ng)%pioVar(idSrad)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idSrad)%gtype=r2dvar
#   endif
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            AVG(ng)%pioVar(idUsms)%vd=var_desc(i)
            AVG(ng)%pioVar(idUsms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUsms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            AVG(ng)%pioVar(idVsms)%vd=var_desc(i)
            AVG(ng)%pioVar(idVsms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVsms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            AVG(ng)%pioVar(idUbms)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbms)%gtype=u2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            AVG(ng)%pioVar(idVbms)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbms)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbms)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbrs))) THEN
            got_var(idUbrs)=.TRUE.
            AVG(ng)%pioVar(idUbrs)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbrs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbrs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbrs))) THEN
            got_var(idVbrs)=.TRUE.
            AVG(ng)%pioVar(idVbrs)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbrs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbrs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbws))) THEN
            got_var(idUbws)=.TRUE.
            AVG(ng)%pioVar(idUbws)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbws)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbws)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbws))) THEN
            got_var(idVbws)=.TRUE.
            AVG(ng)%pioVar(idVbws)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbws)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbws)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbcs))) THEN
            got_var(idUbcs)=.TRUE.
            AVG(ng)%pioVar(idUbcs)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbcs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbcs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbcs))) THEN
            got_var(idVbcs)=.TRUE.
            AVG(ng)%pioVar(idVbcs)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbcs)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbcs)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVwc))) THEN
            got_var(idUVwc)=.TRUE.
            AVG(ng)%pioVar(idUVwc)%vd=var_desc(i)
            AVG(ng)%pioVar(idUVwc)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUVwc)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbot))) THEN
            got_var(idUbot)=.TRUE.
            AVG(ng)%pioVar(idUbot)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbot)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbot)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbot))) THEN
            got_var(idVbot)=.TRUE.
            AVG(ng)%pioVar(idVbot)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbot)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbot)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbur))) THEN
            got_var(idUbur)=.TRUE.
            AVG(ng)%pioVar(idUbur)%vd=var_desc(i)
            AVG(ng)%pioVar(idUbur)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idUbur)%gtype=v2dvar
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbvr))) THEN
            got_var(idVbvr)=.TRUE.
            AVG(ng)%pioVar(idVbvr)%vd=var_desc(i)
            AVG(ng)%pioVar(idVbvr)%dkind=PIO_FOUT
            AVG(ng)%pioVar(idVbvr)%gtype=v2dvar
          END IF
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              AVG(ng)%pioTrc(itrc)%vd=var_desc(i)
              AVG(ng)%pioTrc(itrc)%dkind=PIO_FOUT
              AVG(ng)%pioTrc(itrc)%gtype=r3dvar
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
              got_var(iHUTav(itrc))=.TRUE.
              AVG(ng)%pioVar(iHUTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(iHUTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(iHUTav(itrc))%gtype=u3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,iHVTav(itrc)))) THEN
              got_var(iHVTav(itrc))=.TRUE.
              AVG(ng)%pioVar(iHVTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(iHVTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(iHVTav(itrc))%gtype=v3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUTav(itrc)))) THEN
              got_var(idUTav(itrc))=.TRUE.
              AVG(ng)%pioVar(idUTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idUTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idUTav(itrc))%gtype=u3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVTav(itrc)))) THEN
              got_var(idVTav(itrc))=.TRUE.
              AVG(ng)%pioVar(idVTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idVTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idVTav(itrc))%gtype=v3dvar
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTTav(itrc)))) THEN
              got_var(idTTav(itrc))=.TRUE.
              AVG(ng)%pioVar(idTTav(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idTTav(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idTTav(itrc))%gtype=r3dvar
            END IF
          END DO
#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTrcD(itrc)))) THEN
              got_var(idTrcD(itrc))=.TRUE.
              AVG(ng)%pioVar(idTrcD(itrc))%vd=var_desc(i)
              AVG(ng)%pioVar(idTrcD(itrc))%dkind=PIO_FOUT
              AVG(ng)%pioVar(idTrcD(itrc))%gtype=r3dvar
            END IF
          END DO
#   endif
#  endif
        END DO
!
!  Check if averages variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Aout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Aout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Aout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idFsuD).and.Aout(idFsuD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsuD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dD).and.Aout(idu2dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dD).and.Aout(idv2dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idu2dE).and.Aout(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Aout(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
#   ifdef FORWARD_WRITE
        IF (.not.got_var(idUfx1).and.Aout(idUfx1,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2).and.Aout(idUfx2,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1).and.Aout(idVfx1,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2).and.Aout(idVfx2,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idUvel).and.Aout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Aout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idu3dD).and.Aout(idu3dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dD).and.Aout(idv3dD,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idu3dE).and.Aout(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Aout(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOvel).and.Aout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Aout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Aout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef NEMURO_SED1
        IF (.not.got_var(idPONsed).and.Aout(idPONsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPONsed)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALsed).and.Aout(idOPALsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOPALsed)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDENITsed).and.Aout(idDENITsed,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDENITsed)),      &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPONbur).and.Aout(idPONbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPONbur)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idOPALbur).and.Aout(idOPALbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOPALbur)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef PRIMARY_PROD
        IF (.not.got_var(idNPP).and.Aout(idNPP,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idNPP)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#   ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Aout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Aout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
        IF (.not.got_var(id2dPV).and.Aout(id2dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dRV).and.Aout(id2dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(id3dPV).and.Aout(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Aout(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idZZav).and.Aout(idZZav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZZav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2av).and.Aout(idU2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2av).and.Aout(idV2av,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idHUav).and.Aout(idHUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHVav).and.Aout(idHVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUUav).and.Aout(idUUav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav).and.Aout(idUVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav).and.Aout(idVVav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (.not.got_var(idVvis).and.Aout(idVvis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTdif).and.Aout(idTdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef SALINITY
        IF (.not.got_var(idSdif).and.Aout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#  ifdef TS_VAR
        IF (.not.got_var(idRdif).and.Aout(idRdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idRdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#   if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Aout(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Aout(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Aout(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUaiE).and.Aout(idUaiE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUaiE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVaiN).and.Aout(idVaiN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVaiN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idTsur(itemp)).and.Aout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef SALINITY
        IF (.not.got_var(idTsur(isalt)).and.Aout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   if defined BULK_FLUXES || defined FRC_COUPLING
        IF (.not.got_var(idLhea).and.Aout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Aout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Aout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#   ifdef BULK_FLUXES
        IF (.not.got_var(idTair).and.Aout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef EMINUSP
        IF (.not.got_var(idevap).and.Aout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Aout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#   ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Aout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
        IF (.not.got_var(idUsms).and.Aout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Aout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Aout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Aout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Aout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(iHUTav(itrc)).and.Aout(iHUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc)).and.Aout(iHVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc)).and.Aout(idUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc)).and.Aout(idVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idTTav(itrc)).and.Aout(idTTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        DO itrc=1,NAT
          IF (.not.got_var(idTrcD(itrc)).and.Aout(idTrcD(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTrcD(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#   endif
#  endif

#  if (defined BBL_MODEL || defined WAVES_OUTPUT) && defined SOLVE3D
!
!  Scan bottom boundary layer model and waves variables from input
!  NetCDF and activate switches for time-averaged variables. Get
!  variable IDs.
!
        CALL bbl_def_pio (ng, model, ldef, Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif

#  if defined ICE_MODEL && defined SOLVE3D
!
!  Scan sea-ice variables from input NetCDF and activate switches for
!  time-averaged variables. Get variable IDs.
!
        CALL ice_def_pio (ng, model, ldef, Aout, AVG)
#  endif

# if defined SEDIMENT && defined SOLVE3D
!
!  Scan sediment model variables from input NetCDF and activate
!  switches for time-averaged variables. Get variable IDs.
!
        CALL sediment_def_pio (ng, model, ldef, Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif

#  if defined WEC_VF && defined SOLVE3D
!
!  Scan Waves Effect on Currents variables from input NetCDF and
!  activate switches for time-averaged variables. Get variable IDs.
!
        CALL wec_def_pio (ng, model, ldef, Aout, AVG)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#  endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nAVG(ng)) THEN
          IF (ndefAVG(ng).gt.0) THEN
            AVG(ng)%Rindex=((ntstart(ng)-1)-                            &
     &                      ndefAVG(ng)*((ntstart(ng)-1)/ndefAVG(ng)))/ &
     &                     nAVG(ng)
          ELSE
            AVG(ng)%Rindex=(ntstart(ng)-1)/nAVG(ng)
          END IF
        ELSE
          AVG(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial average time. Notice that the value is offset by half
!  nAVG*dt so there is not a special case when computing its value

#  ifdef ADJOINT
!  in "ad_set_avg".
!
        AVGtime(ng)=time(ng)+0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
#  else
!  in "set_avg".
!
      IF (ntsAVG(ng).eq.1) THEN
        AVGtime(ng)=time(ng)-0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      ELSE
        AVGtime(ng)=time(ng)+REAL(ntsAVG(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      END IF
#  endif
!
  10  FORMAT (2x,'DEF_AVG_PIO      - creating average file,',t56,       &
     &        'Grid ',i2.2,': ',a)
  20  FORMAT (2x,'DEF_AVG_PIO      - inquiring average file',t56,       &
     &        'Grid ',i2.2,': ',a)
  30  FORMAT (/,' DEF_AVG_PIO - unable to create averages NetCDF',      &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG_PIO - unable to open averages NetCDF'         &
     &        ' file: ',a)
  60  FORMAT (/,' DEF_AVG_PIO - unable to find variable: ',a,2x,        &
     &        ' in averages NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE def_avg_pio
# endif
#endif
      END MODULE def_avg_mod
