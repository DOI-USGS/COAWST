#include "cppdefs.h"
#ifdef AVERAGES2
      SUBROUTINE def_avg2 (ng, ldef)
!
!svn $Id: def_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates averages NetCDF file, it defines its           !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      logical :: got_var(NV)

      integer, parameter :: Natt = 25

      integer :: i, ifield, itrc, j, model, nvd3, nvd4
# if defined WRITE_WATER && defined MASKING
      integer :: xy_pdim, xyz_pdim
# endif
      integer :: recdim, status

      integer :: DimIDs(32), p2dgrd(3), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_dim

      real(r8) :: Aval(6)

# ifdef ADJOINT
      character (len=21) :: Prefix
# else
      character (len=13) :: Prefix
# endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile=__FILE__
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
      ncname=AVG2(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) ng, TRIM(ncname)
        ELSE
          WRITE (stdout,20) ng, TRIM(ncname)
        END IF
      END IF
# ifdef ADJOINT
      model=iADM
# else
      model=iNLM
# endif

!
!=======================================================================
!  Create a new averages file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, model, TRIM(ncname), AVG2(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) THEN
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
        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xi_rho',      &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xi_u',        &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xi_v',        &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xi_psi',      &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'eta_rho',     &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'eta_u',       &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'eta_v',       &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'eta_psi',     &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xy_rho',      &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xy_u',        &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xy_v',        &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xy_psi',      &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xyz_rho',     &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xyz_u',       &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xyz_v',       &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'xyz_w',       &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
#  endif

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'tracer',      &
     &                 NT(ng), DimIDs(11))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# endif

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'boundary',    &
     &                 4, DimIDs(14))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, AVG2(ng)%ncid, ncname, 'Nstate',      &
     &                 NstateVar(ng), DimIDs(29))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
# endif

        status=def_dim(ng, model, AVG2(ng)%ncid, ncname,                &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN

        recdim=DimIDs(12)
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=3
# endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(12)
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered u-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(12)
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered v-momentum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(12)
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(12)
# endif
!
!  Define dimension vectors for staggered variables at PSI-points.
!
# if defined WRITE_WATER && defined MASKING
        p2dgrd(1)=xy_pdim
        p2dgrd(2)=DimIDs(12)
# else
        p2dgrd(1)=DimIDs( 4)
        p2dgrd(2)=DimIDs( 8)
        p2dgrd(3)=DimIDs(12)
# endif
!
!  Initialize unlimited time record dimension.
!
        AVG2(ng)%Rindex=0
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
        Prefix='time-averaged adjoint'
# else
        Prefix='time-averaged'
# endif
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, AVG2(ng)%ncid, ncname, DimIDs)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
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
        status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idtime),  &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Define free-surface.
!
        IF (Aout2(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if !defined WET_DRY && (defined WRITE_WATER && defined MASKING)
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idFsur),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Aout2(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbar))
          Vinfo( 3)=Vname(3,idUbar)
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idUbar),&
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Aout2(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbar))
          Vinfo( 3)=Vname(3,idVbar)
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idVbar),&
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D Eastward momentum component at RHO-points.
!
        IF (Aout2(idu2dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu2dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu2dE))
          Vinfo( 3)=Vname(3,idu2dE)
          Vinfo(14)=Vname(4,idu2dE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)='barotropic_eastward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idu2dE), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D Northward momentum component at RHO-points.
!
        IF (Aout2(idv2dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv2dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv2dN))
          Vinfo( 3)=Vname(3,idv2dN)
          Vinfo(14)=Vname(4,idv2dN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(21)='barotropic_northward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idv2dN), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Aout2(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUvel))
          Vinfo( 3)=Vname(3,idUvel)
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idUvel),&
     &                   NF_FOUT, nvd4, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Aout2(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvel))
          Vinfo( 3)=Vname(3,idVvel)
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idVvel),&
     &                   NF_FOUT, nvd4, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 3D Eastward momentum component at RHO-points.
!
        IF (Aout2(idu3dE,ng)) THEN
          Vinfo( 1)=Vname(1,idu3dE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idu3dE))
          Vinfo( 3)=Vname(3,idu3dE)
          Vinfo(14)=Vname(4,idu3dE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)='eastward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idu3dE), &
     &                   NF_FOUT, nvd4, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 3D Northward momentum component at RHO-points.
!
        IF (Aout2(idv3dN,ng)) THEN
          Vinfo( 1)=Vname(1,idv3dN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idv3dN))
          Vinfo( 3)=Vname(3,idv3dN)
          Vinfo(14)=Vname(4,idv3dN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(21)='northward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idv3dN), &
     &                   NF_FOUT, nvd4, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Tid(itrc),&
     &                     NF_FOUT, nvd4, t2dgrd, Aval, Vinfo, ncname)
            IF (FoundError(exit_flag, NoError, __LINE__,                &
     &                     __FILE__)) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Aout2(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idDano),&
     &                   NF_FOUT, nvd4, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Aout2(idHsbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHsbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHsbl))
          Vinfo( 3)=Vname(3,idHsbl)
          Vinfo(14)=Vname(4,idHsbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idHsbl),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Aout2(idHbbl,ng)) THEN
          Vinfo( 1)=Vname(1,idHbbl)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHbbl))
          Vinfo( 3)=Vname(3,idHbbl)
          Vinfo(14)=Vname(4,idHbbl)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idHbbl),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#  endif
# endif
!
!  Define 2D potential vorticity.
!
        IF (Aout2(id2dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dPV))
          Vinfo( 3)=Vname(3,id2dPV)
          Vinfo(14)=Vname(4,id2dPV)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(id2dPV),   &
     &                   NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D relative vorticity.
!
        IF (Aout2(id2dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id2dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id2dRV))
          Vinfo( 3)=Vname(3,id2dRV)
          Vinfo(14)=Vname(4,id2dRV)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(id2dRV),   &
     &                 NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
# ifdef SOLVE3D
!
!  Define 3D potential vorticity.
!
        IF (Aout2(id3dPV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dPV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dPV))
          Vinfo( 3)=Vname(3,id3dPV)
          Vinfo(14)=Vname(4,id3dPV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(id3dPV), &
     &                   NF_FOUT, nvd4, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 3D relative vorticity.
!
        IF (Aout2(id3dRV,ng)) THEN
          Vinfo( 1)=Vname(1,id3dRV)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,id3dRV))
          Vinfo( 3)=Vname(3,id3dRV)
          Vinfo(14)=Vname(4,id3dRV)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_psi'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(id3dRV),   &
     &                   NF_FOUT, nvd4, p2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
#  if defined RUNOFF
!
!  Define surface fresh water flux
!
        IF (Aout2(idRunoff,ng)) THEN
          Vinfo( 1)=Vname(1,idRunoff)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idRunoff))
          Vinfo( 3)=Vname(3,idRunoff)
          Vinfo(14)=Vname(4,idRunoff)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idRunoff,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idRunoff), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Define surface air pressure.
!
        IF (Aout2(idPair,ng)) THEN
          Vinfo( 1)=Vname(1,idPair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idPair))
          Vinfo( 3)=Vname(3,idPair)
          Vinfo(14)=Vname(4,idPair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idPair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
!
!  Define surface winds.
!
        IF (Aout2(idUair,ng)) THEN
          Vinfo( 1)=Vname(1,idUair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUair))
          Vinfo( 3)=Vname(3,idUair)
          Vinfo(14)=Vname(4,idUair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idUair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF

        IF (Aout2(idVair,ng)) THEN
          Vinfo( 1)=Vname(1,idVair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVair))
          Vinfo( 3)=Vname(3,idVair)
          Vinfo(14)=Vname(4,idVair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVair,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idVair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D Eastward surface wind at RHO-points.
!
        IF (Aout2(idUairE,ng)) THEN
          Vinfo( 1)=Vname(1,idUairE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUairE))
          Vinfo( 3)=Vname(3,idUairE)
          Vinfo(14)=Vname(4,idUairE)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUairE,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idUairE),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D Northward surface wind at RHO-points.
!
        IF (Aout2(idVairN,ng)) THEN
          Vinfo( 1)=Vname(1,idVairN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVairN))
          Vinfo( 3)=Vname(3,idVairN)
          Vinfo(14)=Vname(4,idVairN)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVairN,ng),r8)
          status=def_var(ng, iNLM, AVG2(ng)%ncid, AVG2(ng)%Vid(idVairN),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF

#  endif
!
!  Define surface net heat flux.
!
        IF (Aout2(idTsur(itemp),ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idTsur(itemp)), NF_FOUT,          &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Aout2(idTsur(isalt),ng)) THEN
          Vinfo( 1)=Vname(1,idTsur(isalt))
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                          &
     &                                 TRIM(Vname(2,idTsur(isalt)))
          Vinfo( 3)=Vname(3,idTsur(isalt))
          Vinfo(11)='upward flux, freshening (net precipitation)'
          Vinfo(12)='downward flux, salting (net evaporation)'
          Vinfo(14)=Vname(4,idTsur(isalt))
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idTsur(isalt)), NF_FOUT,          &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF

#  ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        IF (Aout2(idLhea,ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Aout2(idShea,ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Aout2(idLrad,ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#   ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Aout2(idevap,ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idevap),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Aout2(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idrain),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#   endif
#  endif
#  ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Aout2(idSrad,ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idSrad),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#  endif
# endif
!
!  Define surface u-momentum stress.
!
        IF (Aout2(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idUsms),&
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Aout2(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idVsms),&
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Aout2(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idUbms),&
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Aout2(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idVbms),&
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
# ifdef SSSFLX
        IF (Aout2(idSSSf,ng)) THEN
          Vinfo( 1)=Vname(1,idSSSf)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSSSf))
          Vinfo( 3)=Vname(3,idSSSf)
          Vinfo(14)=Vname(4,idSSSf)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idSSSf),&
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
# endif
# ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        IF (Aout2(idUice,ng)) THEN
          Vinfo( 1)=Vname(1,idUice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUice))
          Vinfo( 3)=Vname(3,idUice)
          Vinfo(14)=Vname(4,idUice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idUice), NF_FOUT,                 &
     &                   nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Aout2(idVice,ng)) THEN
          Vinfo( 1)=Vname(1,idVice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVice))
          Vinfo( 3)=Vname(3,idVice)
          Vinfo(14)=Vname(4,idVice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idVice), NF_FOUT,                 &
     &                   nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D ice momentum to the East.
!
        IF (Aout2(idUiceE,ng)) THEN
          Vinfo( 1)=Vname(1,idUiceE)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUiceE))
          Vinfo( 3)=Vname(3,idUiceE)
          Vinfo(14)=Vname(4,idUiceE)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUiceE,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idUiceE), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define 2D ice momentum to the North.
!
        IF (Aout2(idViceN,ng)) THEN
          Vinfo( 1)=Vname(1,idViceN)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idViceN))
          Vinfo( 3)=Vname(3,idViceN)
          Vinfo(14)=Vname(4,idViceN)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idViceN,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idViceN), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice concentration.
!
        IF (Aout2(idAice,ng)) THEN
          Vinfo( 1)=Vname(1,idAice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idAice))
          Vinfo( 3)=Vname(3,idAice)
          Vinfo(14)=Vname(4,idAice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idAice), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice average thickness.
!
        IF (Aout2(idHice,ng)) THEN
          Vinfo( 1)=Vname(1,idHice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHice))
          Vinfo( 3)=Vname(3,idHice)
          Vinfo(14)=Vname(4,idHice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idHice), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice/snow surface temperature.
!
        IF (Aout2(idTice,ng)) THEN
          Vinfo( 1)=Vname(1,idTice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTice))
          Vinfo( 3)=Vname(3,idTice)
          Vinfo(14)=Vname(4,idTice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idTice), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define snow thickness.
!
        IF (Aout2(idHsno,ng)) THEN
          Vinfo( 1)=Vname(1,idHsno)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHsno))
          Vinfo( 3)=Vname(3,idHsno)
          Vinfo(14)=Vname(4,idHsno)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsno,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idHsno), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice age.
!
        IF (Aout2(idAgeice,ng)) THEN
          Vinfo( 1)=Vname(1,idAgeice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idAgeice))
          Vinfo( 3)=Vname(3,idAgeice)
          Vinfo(14)=Vname(4,idAgeice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAgeice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idAgeice), NF_FOUT,               &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice-ocean mass flux
!
        IF (Aout2(idIomflx,ng)) THEN
        Vinfo( 1)=Vname(1,idIomflx)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idIomflx))
          Vinfo( 3)=Vname(3,idIomflx)
          Vinfo(14)=Vname(4,idIomflx)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIomflx,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idIomflx), NF_FOUT,               &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define internal ice temperature.
!
        IF (Aout2(idTimid,ng)) THEN
          Vinfo( 1)=Vname(1,idTimid)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTimid))
          Vinfo( 3)=Vname(3,idTimid)
          Vinfo(14)=Vname(4,idTimid)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTimid,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idTimid), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define top ice temperature.
!
        IF (Aout2(idT2ice,ng)) THEN
          Vinfo( 1)=Vname(1,idT2ice)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idT2ice))
          Vinfo( 3)=Vname(3,idT2ice)
          Vinfo(14)=Vname(4,idT2ice)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idT2ice,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idT2ice), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define internal ice stress component 11.
!
        IF (Aout2(idSig11,ng)) THEN
          Vinfo( 1)=Vname(1,idSig11)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSig11))
          Vinfo( 3)=Vname(3,idSig11)
          Vinfo(14)=Vname(4,idSig11)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig11,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idSig11), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define internal ice stress component 12.
!
        IF (Aout2(idSig11,ng)) THEN
          Vinfo( 1)=Vname(1,idSig12)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSig12))
          Vinfo( 3)=Vname(3,idSig12)
          Vinfo(14)=Vname(4,idSig12)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig12,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idSig12), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define internal ice stress component 22.
!
        IF (Aout2(idSig22,ng)) THEN
          Vinfo( 1)=Vname(1,idSig22)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSig22))
          Vinfo( 3)=Vname(3,idSig22)
          Vinfo(14)=Vname(4,idSig22)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSig22,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idSig22), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice-water friction velocity.
!
        IF (Aout2(idTauiw,ng)) THEN
          Vinfo( 1)=Vname(1,idTauiw)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTauiw))
          Vinfo( 3)=Vname(3,idTauiw)
          Vinfo(14)=Vname(4,idTauiw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTauiw,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idTauiw), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice-water momentum transfer coefficient.
!
        IF (Aout2(idChuiw,ng)) THEN
          Vinfo( 1)=Vname(1,idChuiw)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idChuiw))
          Vinfo( 3)=Vname(3,idChuiw)
          Vinfo(14)=Vname(4,idChuiw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idChuiw,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idChuiw), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define temperature of molecular sublayer under ice.
!
        IF (Aout2(idT0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idT0mk)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idT0mk))
          Vinfo( 3)=Vname(3,idT0mk)
          Vinfo(14)=Vname(4,idT0mk)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idT0mk,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idT0mk), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define salinity of molecular sublayer under ice.
!
        IF (Aout2(idS0mk,ng)) THEN
          Vinfo( 1)=Vname(1,idS0mk)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idS0mk))
          Vinfo( 3)=Vname(3,idS0mk)
          Vinfo(14)=Vname(4,idS0mk)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idS0mk,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idS0mk), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define frazil ice growth rate.
!
        IF (Aout2(idWfr,ng)) THEN
          Vinfo( 1)=Vname(1,idWfr)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWfr))
          Vinfo( 3)=Vname(3,idWfr)
          Vinfo(14)=Vname(4,idWfr)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWfr,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idWfr), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice growth/melt rate.
!
        IF (Aout2(idWai,ng)) THEN
          Vinfo( 1)=Vname(1,idWai)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWai))
          Vinfo( 3)=Vname(3,idWai)
          Vinfo(14)=Vname(4,idWai)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWfr,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idWai), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice growth/melt rate.
!
        IF (Aout2(idWao,ng)) THEN
          Vinfo( 1)=Vname(1,idWao)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWao))
          Vinfo( 3)=Vname(3,idWao)
          Vinfo(14)=Vname(4,idWao)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWao,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idWao), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice growth/melt rate.
!
        IF (Aout2(idWio,ng)) THEN
          Vinfo( 1)=Vname(1,idWio)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWio))
          Vinfo( 3)=Vname(3,idWio)
          Vinfo(14)=Vname(4,idWio)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWio,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idWio), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice melt runoff rate.
!
        IF (Aout2(idWro,ng)) THEN
          Vinfo( 1)=Vname(1,idWro)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWro))
          Vinfo( 3)=Vname(3,idWro)
          Vinfo(14)=Vname(4,idWro)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWro,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idWro), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice divergence rate.
!
        IF (Aout2(idWdiv,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiv)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdiv))
          Vinfo( 3)=Vname(3,idWdiv)
          Vinfo(14)=Vname(4,idWdiv)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiv,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idWdiv), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
# endif
#ifdef ICE_DIAGS
!
!  Define ocean albedo
!
        IF (Aout2(idAlbe,ng)) THEN
          Vinfo( 1)=Vname(1,idAlbe)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idAlbe))
          Vinfo( 3)=Vname(3,idAlbe)
          Vinfo(14)=Vname(4,idAlbe)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAlbe,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idAlbe),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice albedo
!
        IF (Aout2(idAlbI,ng)) THEN
          Vinfo( 1)=Vname(1,idAlbI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idAlbI))
          Vinfo( 3)=Vname(3,idAlbI)
          Vinfo(14)=Vname(4,idAlbI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idAlbI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idAlbI),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define downward longwave radiation
!
        IF (Aout2(idLWdn,ng)) THEN
          Vinfo( 1)=Vname(1,idLWdn)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLWdn))
          Vinfo( 3)=Vname(3,idLWdn)
          Vinfo(14)=Vname(4,idLWdn)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLWdn,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idLWdn),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define downward shortwave radiation
!
        IF (Aout2(idSWdn,ng)) THEN
          Vinfo( 1)=Vname(1,idSWdn)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSWdn))
          Vinfo( 3)=Vname(3,idSWdn)
          Vinfo(14)=Vname(4,idSWdn)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSWdn,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idSWdn),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define latent heat over ice
!
        IF (Aout2(idLatI,ng)) THEN
          Vinfo( 1)=Vname(1,idLatI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLatI))
          Vinfo( 3)=Vname(3,idLatI)
          Vinfo(14)=Vname(4,idLatI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLatI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idLatI),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define sensible heat over ice
!
        IF (Aout2(idSenI,ng)) THEN
          Vinfo( 1)=Vname(1,idSenI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSenI))
          Vinfo( 3)=Vname(3,idSenI)
          Vinfo(14)=Vname(4,idSenI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSenI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idSenI),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define upward longwave over ice
!
        IF (Aout2(idLWupI,ng)) THEN
          Vinfo( 1)=Vname(1,idLWupI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idLWupI))
          Vinfo( 3)=Vname(3,idLWupI)
          Vinfo(14)=Vname(4,idLWupI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLWupI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                     AVG2(ng)%Vid(idLWupI),                       &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define upward shortwave over ice
!
        IF (Aout2(idSWupI,ng)) THEN
          Vinfo( 1)=Vname(1,idSWupI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSWupI))
          Vinfo( 3)=Vname(3,idSWupI)
          Vinfo(14)=Vname(4,idSWupI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSWupI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                     AVG2(ng)%Vid(idSWupI),                       &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define shortwave into ice
!
        IF (Aout2(idSWinI,ng)) THEN
          Vinfo( 1)=Vname(1,idSWinI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSWinI))
          Vinfo( 3)=Vname(3,idSWinI)
          Vinfo(14)=Vname(4,idSWinI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSWinI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                     AVG2(ng)%Vid(idSWinI),                       &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define shortwave through ice
!
        IF (Aout2(idSWthruI,ng)) THEN
          Vinfo( 1)=Vname(1,idSWthruI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSWthruI))
          Vinfo( 3)=Vname(3,idSWthruI)
          Vinfo(14)=Vname(4,idSWthruI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSWthruI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                   AVG2(ng)%Vid(idSWthruI),                       &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define salt flux under ice
!
        IF (Aout2(idSaltI,ng)) THEN
          Vinfo( 1)=Vname(1,idSaltI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSaltI))
          Vinfo( 3)=Vname(3,idSaltI)
          Vinfo(14)=Vname(4,idSaltI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSaltI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                     AVG2(ng)%Vid(idSaltI),                       &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define salt flux through leads
!
        IF (Aout2(idSaltL,ng)) THEN
          Vinfo( 1)=Vname(1,idSaltL)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSaltL))
          Vinfo( 3)=Vname(3,idSaltL)
          Vinfo(14)=Vname(4,idSaltL)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSaltL,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid,                      &
     &                     AVG2(ng)%Vid(idSaltL),                       &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define heat flux through leads
!
        IF (Aout2(idQao,ng)) THEN
          Vinfo( 1)=Vname(1,idQao)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idQao))
          Vinfo( 3)=Vname(3,idQao)
          Vinfo(14)=Vname(4,idQao)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idQao,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idQao), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice-air heat flux
!
        IF (Aout2(idQai,ng)) THEN
          Vinfo( 1)=Vname(1,idQai)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idQai))
          Vinfo( 3)=Vname(3,idQai)
          Vinfo(14)=Vname(4,idQai)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idQai,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idQai), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice-ocean heat flux
!
        IF (Aout2(idQio,ng)) THEN
          Vinfo( 1)=Vname(1,idQio)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idQio))
          Vinfo( 3)=Vname(3,idQio)
          Vinfo(14)=Vname(4,idQio)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idQio,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idQio), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define ice heat flux
!
        IF (Aout2(idQi2,ng)) THEN
          Vinfo( 1)=Vname(1,idQi2)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idQi2))
          Vinfo( 3)=Vname(3,idQi2)
          Vinfo(14)=Vname(4,idQi2)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idQi2,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idQi2), &
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
!
!  Define snow-ice conversion
!
        IF (Aout2(idSnoI,ng)) THEN
          Vinfo( 1)=Vname(1,idSnoI)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSnoI))
          Vinfo( 3)=Vname(3,idSnoI)
          Vinfo(14)=Vname(4,idSnoI)
          Vinfo(16)=Vname(1,idtime)
#if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSnoI,ng),r8)
          status=def_var(ng, model, AVG2(ng)%ncid, AVG2(ng)%Vid(idSnoI),&
     &                     NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (FoundError(exit_flag, NoError, __LINE__,                  &
     &                   __FILE__)) RETURN
        END IF
#endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, AVG2(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, AVG2(ng)%ncid, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVG2(ng)%name
!
!  Open averages file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, AVG2(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname,                       &
     &                         ncid = AVG2(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname,                         &
     &                       ncid = AVG2(ng)%ncid)
        IF (FoundError(exit_flag, NoError, __LINE__,                    &
     &                 __FILE__)) RETURN
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
            AVG2(ng)%Vid(idtime)=var_id(i)
# ifdef WET_DRY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRwet))) THEN
            got_var(idRwet)=.TRUE.
            AVG2(ng)%Vid(idRwet)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            AVG2(ng)%Vid(idFsur)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            AVG2(ng)%Vid(idUbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            AVG2(ng)%Vid(idVbar)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dE))) THEN
            got_var(idu2dE)=.TRUE.
            AVG2(ng)%Vid(idu2dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dN))) THEN
            got_var(idv2dN)=.TRUE.
            AVG2(ng)%Vid(idv2dN)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            AVG2(ng)%Vid(idUvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            AVG2(ng)%Vid(idVvel)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dE))) THEN
            got_var(idu3dE)=.TRUE.
            AVG2(ng)%Vid(idu3dE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dN))) THEN
            got_var(idv3dN)=.TRUE.
            AVG2(ng)%Vid(idv3dN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            AVG2(ng)%Vid(idDano)=var_id(i)
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            AVG2(ng)%Vid(idHsbl)=var_id(i)
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            AVG2(ng)%Vid(idHbbl)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dPV))) THEN
            got_var(id2dPV)=.TRUE.
            AVG2(ng)%Vid(id2dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id2dRV))) THEN
            got_var(id2dRV)=.TRUE.
            AVG2(ng)%Vid(id2dRV)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dPV))) THEN
            got_var(id3dPV)=.TRUE.
            AVG2(ng)%Vid(id3dPV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,id3dRV))) THEN
            got_var(id3dRV)=.TRUE.
            AVG2(ng)%Vid(id3dRV)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSSSf))) THEN
            got_var(idSSSf)=.TRUE.
            AVG2(ng)%Vid(idSSSf)=var_id(i)
#  ifdef ICE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.true.
            AVG2(ng)%Vid(idUice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.true.
            AVG2(ng)%Vid(idVice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiceE))) THEN
            got_var(idUiceE)=.true.
            AVG2(ng)%Vid(idUiceE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViceN))) THEN
            got_var(idViceN)=.true.
            AVG2(ng)%Vid(idViceN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.true.
            AVG2(ng)%Vid(idAice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.true.
            AVG2(ng)%Vid(idHice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.true.
            AVG2(ng)%Vid(idTice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.true.
            AVG2(ng)%Vid(idHsno)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.true.
            AVG2(ng)%Vid(idAgeice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
            got_var(idIomflx)=.true.
            AVG2(ng)%Vid(idIomflx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.true.
            AVG2(ng)%Vid(idTimid)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT2ice))) THEN
            got_var(idT2ice)=.true.
            AVG2(ng)%Vid(idT2ice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.true.
            AVG2(ng)%Vid(idSig11)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.true.
            AVG2(ng)%Vid(idSig12)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.true.
            AVG2(ng)%Vid(idSig22)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.true.
            AVG2(ng)%Vid(idTauiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.true.
            AVG2(ng)%Vid(idChuiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.true.
            AVG2(ng)%Vid(idT0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.true.
            AVG2(ng)%Vid(idS0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
            got_var(idWfr)=.true.
            AVG2(ng)%Vid(idWfr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
            got_var(idWai)=.true.
            AVG2(ng)%Vid(idWai)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
            got_var(idWao)=.true.
            AVG2(ng)%Vid(idWao)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
            got_var(idWio)=.true.
            AVG2(ng)%Vid(idWio)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
            got_var(idWro)=.true.
            AVG2(ng)%Vid(idWro)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiv))) THEN
            got_var(idWdiv)=.true.
            AVG2(ng)%Vid(idWdiv)=var_id(i)
#  endif
#  ifdef ICE_DIAGS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAlbe))) THEN
            got_var(idAlbe)=.true.
            AVG2(ng)%Vid(idAlbe)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAlbI))) THEN
            got_var(idAlbI)=.true.
            AVG2(ng)%Vid(idAlbI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLWdn))) THEN
            got_var(idLWdn)=.true.
            AVG2(ng)%Vid(idLWdn)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSWdn))) THEN
            got_var(idSWdn)=.true.
            AVG2(ng)%Vid(idSWdn)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLatI))) THEN
            got_var(idLatI)=.true.
            AVG2(ng)%Vid(idLatI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSenI))) THEN
            got_var(idSenI)=.true.
            AVG2(ng)%Vid(idSenI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLWupI))) THEN
            got_var(idLWupI)=.true.
            AVG2(ng)%Vid(idLWupI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSWupI))) THEN
            got_var(idSWupI)=.true.
            AVG2(ng)%Vid(idSWupI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSWinI))) THEN
            got_var(idSWinI)=.true.
            AVG2(ng)%Vid(idSWinI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSWthruI))) THEN
            got_var(idSWthruI)=.true.
            AVG2(ng)%Vid(idSWthruI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSaltI))) THEN
            got_var(idSaltI)=.true.
            AVG2(ng)%Vid(idSaltI)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSaltL))) THEN
            got_var(idSaltL)=.true.
            AVG2(ng)%Vid(idSaltL)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idQao))) THEN
            got_var(idQao)=.true.
            AVG2(ng)%Vid(idQao)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idQai))) THEN
            got_var(idQai)=.true.
            AVG2(ng)%Vid(idQai)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idQio))) THEN
            got_var(idQio)=.true.
            AVG2(ng)%Vid(idQio)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idQi2))) THEN
            got_var(idQi2)=.true.
            AVG2(ng)%Vid(idQi2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSnoI))) THEN
            got_var(idSnoI)=.true.
            AVG2(ng)%Vid(idSnoI)=var_id(i)
#  endif
#  if defined RUNOFF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRunoff))) THEN
            got_var(idRunoff)=.TRUE.
            AVG2(ng)%Vid(idRunoff)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idPair))) THEN
            got_var(idPair)=.TRUE.
            AVG2(ng)%Vid(idPair)=var_id(i)
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUair))) THEN
            got_var(idUair)=.TRUE.
            AVG2(ng)%Vid(idUair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVair))) THEN
            got_var(idVair)=.TRUE.
            AVG2(ng)%Vid(idVair)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUairE))) THEN
            got_var(idUairE)=.TRUE.
            AVG2(ng)%Vid(idUairE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVairN))) THEN
            got_var(idVairN)=.TRUE.
            AVG2(ng)%Vid(idVairN)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            AVG2(ng)%Vid(idTsur(itemp))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            AVG2(ng)%Vid(idTsur(isalt))=var_id(i)
#  ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            AVG2(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            AVG2(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            AVG2(ng)%Vid(idLrad)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            AVG2(ng)%Vid(idTair)=var_id(i)
#   ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            AVG2(ng)%Vid(idevap)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            AVG2(ng)%Vid(idrain)=var_id(i)
#   endif
#  endif
#  ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            AVG2(ng)%Vid(idSrad)=var_id(i)
#  endif
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            AVG2(ng)%Vid(idUsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            AVG2(ng)%Vid(idVsms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            AVG2(ng)%Vid(idUbms)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            AVG2(ng)%Vid(idVbms)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             AVG2(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
              got_var(iHUTav(itrc))=.TRUE.
              AVG2(ng)%Vid(iHUTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,iHVTav(itrc)))) THEN
              got_var(iHVTav(itrc))=.TRUE.
              AVG2(ng)%Vid(iHVTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUTav(itrc)))) THEN
              got_var(idUTav(itrc))=.TRUE.
              AVG2(ng)%Vid(idUTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVTav(itrc)))) THEN
              got_var(idVTav(itrc))=.TRUE.
              AVG2(ng)%Vid(idVTav(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTTav(itrc)))) THEN
              got_var(idTTav(itrc))=.TRUE.
              AVG2(ng)%Vid(idTTav(itrc))=var_id(i)
            END IF
          END DO
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
# ifdef WET_DRY
        IF (.not.got_var(idFsur).and.Aout2(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idUbar).and.Aout2(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Aout2(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dE).and.Aout2(idu2dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dN).and.Aout2(idv2dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Aout2(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Aout2(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu3dE).and.Aout2(idu3dE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dE)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dN).and.Aout2(idv3dN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dN)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Aout2(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Aout2(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Aout2(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(id2dPV).and.Aout2(id2dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id2dRV).and.Aout2(id2dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id2dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(id3dPV).and.Aout2(id3dPV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dPV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(id3dRV).and.Aout2(id3dRV,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,id3dRV)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSSSf).and.Aout2(idSSSf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSSSf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined RUNOFF
        IF (.not.got_var(idRunoff).and.Aout2(idRunoff,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idRunoff)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
        IF (.not.got_var(idPair).and.Aout2(idPair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ICE_MODEL
        IF (.not.got_var(idUice).and.Aout2(idUice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Aout2(idVice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUiceE).and.Aout2(idUiceE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUiceE)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViceN).and.Aout2(idViceN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idViceN)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Aout2(idAice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Aout2(idHice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Aout2(idTice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Aout2(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAgeice).and.Aout2(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Aout2(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Aout2(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT2ice).and.Aout2(idT2ice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idT2ice)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Aout2(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Aout2(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Aout2(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Aout2(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Aout2(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Aout2(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Aout2(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Aout2(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Aout2(idWai,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Aout2(idWao,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Aout2(idWio,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Aout2(idWro,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiv).and.Aout2(idWdiv,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdiv)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef ICE_DIAGS
        IF (.not.got_var(idAlbe).and.Aout2(idAlbe,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAlbe)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAlbI).and.Aout2(idAlbI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAlbI)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLWdn).and.Aout2(idLWdn,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLWdn)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWdn).and.Aout2(idSWdn,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSWdn)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLatI).and.Aout2(idLatI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLatI)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSenI).and.Aout2(idSenI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSenI)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLWupI).and.Aout2(idLWupI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLWupI)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWupI).and.Aout2(idSWupI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSWupI)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWinI).and.Aout2(idSWinI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSWinI)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSWthruI).and.Aout2(idSWthruI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSWthruI)),       &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSaltI).and.Aout2(idSaltI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSaltI)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSaltL).and.Aout2(idSaltL,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSaltL)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idQao).and.Aout2(idQao,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idQao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idQai).and.Aout2(idQai,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idQai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idQio).and.Aout2(idQio,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idQio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idQi2).and.Aout2(idQi2,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idQi2)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSnoI).and.Aout2(idSnoI,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSnoI)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
        IF (.not.got_var(idUair).and.Aout2(idUair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVair).and.Aout2(idVair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUairE).and.Aout2(idUairE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUairE)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVairN).and.Aout2(idVairN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVairN)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idTsur(itemp)).and.Aout2(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Aout2(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Aout2(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Aout2(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Aout2(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Aout2(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef EMINUSP
        IF (.not.got_var(idevap).and.Aout2(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Aout2(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
#  ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Aout2(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
        IF (.not.got_var(idUsms).and.Aout2(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Aout2(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Aout2(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Aout2(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Aout2(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO itrc=1,NAT
          IF (.not.got_var(iHUTav(itrc)).and.Aout2(iHUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc)).and.Aout2(iHVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc)).and.Aout2(idUTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc)).and.Aout2(idVTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idTTav(itrc)).and.Aout2(idTTav(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (nRST(ng).eq.nAVG2(ng)) THEN
          IF (ndefAVG2(ng).gt.0) THEN
            AVG2(ng)%Rindex=((ntstart(ng)-1)-                           &
     &                      ndefAVG2(ng)*((ntstart(ng)-1)/              &
     &                      ndefAVG2(ng)))/nAVG2(ng)
          ELSE
            AVG2(ng)%Rindex=(ntstart(ng)-1)/nAVG2(ng)
          END IF
        ELSE
          AVG2(ng)%Rindex=rec_size
        END IF
      END IF QUERY
!
!  Set initial average time. Notice that the value is offset by half
!  nAVG2*dt so there is not a special case when computing its value

# ifdef ADJOINT
!  in "ad_set_avg".
!
      AVG2time(ng)=time(ng)+0.5_r8*REAL(nAVG2(ng),r8)*dt(ng)
# else
!  in "set_avg".
!
      IF (ntsAVG2(ng).eq.1) THEN
        AVG2time(ng)=time(ng)-0.5_r8*REAL(nAVG2(ng),r8)*dt(ng)
      ELSE
        AVG2time(ng)=time(ng)+REAL(ntsAVG2(ng),r8)*dt(ng)-              &
     &              0.5_r8*REAL(nAVG2(ng),r8)*dt(ng)
      END IF
# endif
!
  10  FORMAT (6x,'DEF_AVG2     - creating  average', t43,               &
     &        ' file, Grid ',i2.2,': ', a)
  20  FORMAT (6x,'DEF_AVG2     - inquiring average', t43,               &
     &        ' file, Grid ',i2.2,': ', a)
  30  FORMAT (/,' DEF_AVG2 - unable to create averages NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG2 - unable to open averages NetCDF file: ',a)
  60  FORMAT (/,' DEF_AVG2 - unable to find variable: ',a,2x,           &
     &        ' in averages NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_avg2
#else
      SUBROUTINE def_avg2 (ng, ldef)
      integer, intent(in) :: ng
      logical, intent(in) :: ldef
      RETURN
      END SUBROUTINE def_avg2
#endif
