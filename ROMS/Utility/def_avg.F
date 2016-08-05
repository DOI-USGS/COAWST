#include "cppdefs.h"
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      SUBROUTINE def_avg (ng, ldef)
!
!svn $Id: def_avg.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
# ifdef FILTERED
      USE mod_filter, ONLY: nfil
# endif
# ifdef SEDIMENT
      USE mod_sediment
# endif
!
      USE def_var_mod, ONLY : def_var
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

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: p3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      real(r8) :: Aval(6)

# ifdef ADJOINT
      character (len=21) :: Prefix
# else
      character (len=13) :: Prefix
# endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_avg.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=AVG(ng)%name
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
        CALL netcdf_create (ng, model, TRIM(ncname), AVG(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
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
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_u',         &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_v',         &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xi_psi',       &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_rho',      &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_u',        &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_v',        &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'eta_psi',      &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_rho',       &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_u',         &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_v',         &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xy_psi',       &
     &                 IOBOUNDS(ng)%xy_psi, xy_pdim)
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_rho',      &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_u',        &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_v',        &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_w',        &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xyz_psi',      &
     &                 IOBOUNDS(ng)%xy_psi*N(ng), xyz_pdim)
        IF (exit_flag.ne.NoError) RETURN
#  endif

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, AVG(ng)%ncid, ncname, 'NST',           &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'xybed',        &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
#  if defined VARIABLE_CDOM && defined SPECTRAL_LIGHT
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef ECOSIM
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, AVG(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, model, AVG(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

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
        Prefix='time-averaged adjoint'
# else
        Prefix='time-averaged'
# endif
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, model, AVG(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        WRITE (Vinfo( 2),'(a,1x,a)') 'averaged', TRIM(Vname(2,idtime))
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idtime),    &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)

        IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idFsur),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsuD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idFsuD),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbar),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu2dD),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbar),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv2dD),  &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(21)='barotropic_eastward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu2dE,ng),r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idu2dE),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(21)='barotropic_northward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv2dN,ng),r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idv2dN),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUfx1),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUfx2),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVfx1),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVfx2),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUvel),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idu3dD),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVvel),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dD,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idv3dD),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(21)='eastward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idu3dE,ng),r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idu3dE),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(21)='northward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idv3dN,ng),r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idv3dN),   &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idOvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(21)='upward_sea_water_velocity'
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWvel),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Tid(itrc),  &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO

# ifdef BEST_NPZ
#ifdef STATIONARY
!-------------------------
!3D Stationary Variable Array
!--------------------------
        DO itrc=1,NTS(ng)
          IF (Aout(idTSvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTSvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTSvar(itrc)))
            Vinfo( 3)=Vname(3,idTSvar(itrc))
            Vinfo(14)=Vname(4,idTSvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model,AVG(ng)%ncid,AVG(ng)%TSid(itrc),   &
     &                     NF_FOUT,nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#endif
#ifdef STATIONARY2
!-------------------------
!2D Stationary Array
!--------------------------
        DO itrc=1,NTS2(ng)
          IF (Aout(idTS2var(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTS2var(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTS2var(itrc)))
            Vinfo( 3)=Vname(3,idTS2var(itrc))
            Vinfo(14)=Vname(4,idTS2var(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, model,AVG(ng)%ncid,AVG(ng)%TS2id(itrc),  &
     &                     NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#endif
#ifdef BIOFLUX
!-------------------------
!bioflux Array
!--------------------------
!       DO itrc=1,NT(ng)
!         DO itrc2=1,NT(ng)
            IF (Aout(idTBFvar(iBF),ng)) THEN
              Vinfo( 1)=Vname(1,idTBFvar(iBF))
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTBFvar(iBF)))
              Vinfo( 3)=Vname(3,idTBFvar(iBF))
              Vinfo(14)=Vname(4,idTBFvar(iBF))
              Vinfo(16)=Vname(1,idtime)
              Vinfo(22)='coordinates'
              Aval(5)=REAL(r2dvar,r8)
              status=def_var(ng, model,AVG(ng)%ncid,AVG(ng)%TBFid(iBF),   &
     &                     NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
            END IF
!        END DO
!        END DO
#endif
#ifdef PROD3
!-------------------------------
!3D Production variables
!-------------------------------
        DO itrc=1,NPT3(ng)
          IF (Aout(idPT3var(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idPT3var(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idPT3var(itrc)))
            Vinfo( 3)=Vname(3,idPT3var(itrc))
            Vinfo(14)=Vname(4,idPT3var(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model,AVG(ng)%ncid,AVG(ng)%PT3id(itrc),  &
     &                    NF_FOUT,nvd4,t3dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#endif
#ifdef PROD2
!-------------------------------
!2D Production variables
!-------------------------------
        DO itrc=1,NPT2(ng)
          IF (Aout(idPT2var(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idPT2var(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idPT2var(itrc)))
            Vinfo( 3)=Vname(3,idPT2var(itrc))
            Vinfo(14)=Vname(4,idPT2var(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model,AVG(ng)%ncid,AVG(ng)%PT2id(itrc),  &
     &                     NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#endif
!------------------
!benthic bio variables
!-------------------
#  if defined BENTHIC
        DO itrc=1,NBeT(ng)
          IF (Aout(idBeTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idBeTvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idBeTvar(itrc)))

            Vinfo( 1)=Vname(1,idBeTvar(itrc))
            Vinfo( 2)=Vname(2,idBeTvar(itrc))
            Vinfo(3)=Vname(3,idBeTvar(itrc))
            Vinfo(14)=Vname(4,idBeTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, model,ncavgid(ng),avgBid(itrc,ng),  &
     &                     NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#  endif
!------------------
!ice bio variables
!-------------------
# if defined ICE_BIO
#   ifdef CLIM_ICE_1D
        DO itrc=1,NIceT(ng)
          IF (Aout(idIceBvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idIceBvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                  TRIM(Vname(2,idIceBvar(itrc)))
            Vinfo( 2)=Vname(2,idIceBvar(itrc))
            Vinfo(3)=Vname(3,idIceBvar(itrc))
            Vinfo(14)=Vname(4,idIceBvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Aval(5)=REAL(r2dvar,r8)
            status=def_var(ng, model,ncavgid(ng),avgIceBid(itrc,ng),    &
     &                     NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
#elif defined BERING_10K
!
!  Define ice algae.
!
        IF (Aout(idIcePhL,ng)) THEN
          Vinfo( 1)=Vname(1,idIcePhL)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idIcePhL))
          Vinfo( 3)=Vname(3,idIcePhL)
          Vinfo(14)=Vname(4,idIcePhL)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIcePhL,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idIcePhL), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!  Define ice nitrate.
!
        IF (Aout(idIceNO3,ng)) THEN
          Vinfo( 1)=Vname(1,idIceNO3)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idIceNO3))
          Vinfo( 3)=Vname(3,idIceNO3)
          Vinfo(14)=Vname(4,idIceNO3)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIceNO3,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idIceNO3), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!  Define ice ammonium.
!
        IF (Aout(idIceNH4,ng)) THEN
          Vinfo( 1)=Vname(1,idIceNH4)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idIceNH4))
          Vinfo( 3)=Vname(3,idIceNH4)
          Vinfo(14)=Vname(4,idIceNH4)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIceNH4,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idIceNH4), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!  Define ice log.
!
        IF (Aout(idIceLog,ng)) THEN
          Vinfo( 1)=Vname(1,idIceLog)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idIceLog))
          Vinfo( 3)=Vname(3,idIceLog)
          Vinfo(14)=Vname(4,idIceLog)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idIceLog,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idIceLog), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
#endif

# endif
# ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
          IF (Aout(idTSvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTSvar(itrc))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idTSvar(itrc)))
            Vinfo( 3)=Vname(3,idTSvar(itrc))
            Vinfo(14)=Vname(4,idTSvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%TSid(itrc), &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif

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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(ifield), NF_FOUT,                &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idDano),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          IF (exit_flag.ne.NoError) RETURN
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
          IF (exit_flag.ne.NoError) RETURN
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
          IF (exit_flag.ne.NoError) RETURN
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
          IF (exit_flag.ne.NoError) RETURN
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
          IF (exit_flag.ne.NoError) RETURN
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
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHsbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHbbl),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id2dPV),   &
     &                   NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p2dvar,r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id2dRV),   &
     &                 NF_FOUT, nvd3, p2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id3dPV),   &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(p3dvar,r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(id3dRV),   &
     &                   NF_FOUT, nvd4, p3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idZZav),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU2av),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV2av),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHUav),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idHVav),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUUav),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUVav),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVVav),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idTTav(itrc)), NF_FOUT,          &
     &                     nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(iHUTav(itrc)), NF_FOUT,          &
     &                     nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(iHVTav(itrc)), NF_FOUT,          &
     &                     nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(u3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUTav(itrc)), NF_FOUT,          &
     &                     nvd4, u3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(v3dvar,r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVTav(itrc)), NF_FOUT,          &
     &                     nvd4, v3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVvis),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idTdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Aout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSdif),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  if defined RUNOFF
!
!  Define surface air pressure.
!
        IF (Aout(idRunoff,ng)) THEN
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
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idRunoff), &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idPair,ng),r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idPair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUair,ng),r8)
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

        IF (Aout(idVair,ng)) THEN
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
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVair),   &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Eastward surface wind at RHO-points.
!
        IF (Aout(idUairE,ng)) THEN
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
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idUairE),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Northward surface wind at RHO-points.
!
        IF (Aout(idVairN,ng)) THEN
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
          status=def_var(ng, iNLM, AVG(ng)%ncid, AVG(ng)%Vid(idVairN),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTsur(itemp)), NF_FOUT,           &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
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
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTsur(isalt)), NF_FOUT,           &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  ifdef BULK_FLUXES
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idLhea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idShea),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idLrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idevap,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idevap),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idrain),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#   endif
#  endif
#  ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSrad),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUsms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVsms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbms),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbms),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifdef BBL
!
!  Define bottom u-momentum stress at rho points.
!
        IF (Aout(idUbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbrs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbrs))
          Vinfo( 3)=Vname(3,idUbrs)
          Vinfo(14)=Vname(4,idUbrs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbrs,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbrs),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom v-momentum stress at rho points.
!
        IF (Aout(idVbrs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbrs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbrs))
          Vinfo( 3)=Vname(3,idVbrs)
          Vinfo(14)=Vname(4,idVbrs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbrs,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbrs),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom u-wave stress at rho points.
!
        IF (Aout(idUbws,ng)) THEN
          Vinfo( 1)=Vname(1,idUbws)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbws))
          Vinfo( 3)=Vname(3,idUbws)
          Vinfo(14)=Vname(4,idUbws)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbws,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbws),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom v-wave stress at rho points.
!
        IF (Aout(idVbws,ng)) THEN
          Vinfo( 1)=Vname(1,idVbws)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbws))
          Vinfo( 3)=Vname(3,idVbws)
          Vinfo(14)=Vname(4,idVbws)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbws,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbws),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom max u wave-curr stress at rho points.
!
        IF (Aout(idUbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idUbcs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbcs))
          Vinfo( 3)=Vname(3,idUbcs)
          Vinfo(14)=Vname(4,idUbcs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbcs,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbcs),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom max v wave-curr stress at rho points.
!
        IF (Aout(idVbcs,ng)) THEN
          Vinfo( 1)=Vname(1,idVbcs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbcs))
          Vinfo( 3)=Vname(3,idVbcs)
          Vinfo(14)=Vname(4,idVbcs)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbcs,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbcs),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom max uv wave-curr stress magnitude at rho points.
!
        IF (Aout(idUVwc,ng)) THEN
          Vinfo( 1)=Vname(1,idUVwc)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUVwc))
          Vinfo( 3)=Vname(3,idUVwc)
          Vinfo(14)=Vname(4,idUVwc)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUVwc,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUVwc),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

!
!  Define bed uwave orbital vel at rho points.
!
        IF (Aout(idUbot,ng)) THEN
          Vinfo( 1)=Vname(1,idUbot)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbot))
          Vinfo( 3)=Vname(3,idUbot)
          Vinfo(14)=Vname(4,idUbot)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbot,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbot),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bed vwave orbital vel at rho points.
!
        IF (Aout(idVbot,ng)) THEN
          Vinfo( 1)=Vname(1,idVbot)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbot))
          Vinfo( 3)=Vname(3,idVbot)
          Vinfo(14)=Vname(4,idVbot)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbot,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbot),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom u-vel above bed at rho points.
!
        IF (Aout(idUbur,ng)) THEN
          Vinfo( 1)=Vname(1,idUbur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbur))
          Vinfo( 3)=Vname(3,idUbur)
          Vinfo(14)=Vname(4,idUbur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbur,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUbur),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom u-vel above bed at rho points.
!
        IF (Aout(idVbvr,ng)) THEN
          Vinfo( 1)=Vname(1,idVbvr)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbvr))
          Vinfo( 3)=Vname(3,idVbvr)
          Vinfo(14)=Vname(4,idVbvr)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbvr,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVbvr),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif

# if defined SEDIMENT && defined BEDLOAD
        DO i=1,NST
!
!  Define Bedload U-direction.
!
          IF (Aout(idUbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idUbld(i))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idUbld(i)))
            Vinfo( 3)=Vname(3,idUbld(i))
            Vinfo(14)=Vname(4,idUbld(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idUbld(i)), NF_FOUT,             &
     &                     nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Define Bedload V-direction.
!
          IF (Aout(idVbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idVbld(i))
            WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idVbld(i)))
            Vinfo( 3)=Vname(3,idVbld(i))
            Vinfo(14)=Vname(4,idVbld(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
            status=def_var(ng, model, AVG(ng)%ncid,                     &
     &                     AVG(ng)%Vid(idVbld(i)), NF_FOUT,             &
     &                     nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
# ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idSSSf),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          Aval(5)=REAL(r2dvar,r8)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        IF (Aout(idUice,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idUice), NF_FOUT,                  &
     &                   nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Aout(idVice,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idVice), NF_FOUT,                  &
     &                   nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum to the East.
!
        IF (Aout(idUiceE,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idUiceE), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D ice momentum to the North.
!
        IF (Aout(idViceN,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idViceN), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice concentration.
!
        IF (Aout(idAice,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idAice), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice average thickness.
!
        IF (Aout(idHice,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idHice), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice/snow surface temperature.
!
        IF (Aout(idTice,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTice), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define snow thickness.
!
        IF (Aout(idHsno,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idHsno), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef MELT_PONDS
!
!  Define surface water (on ice) fraction.
!
        IF (Aout(idApond,ng)) THEN
          Vinfo( 1)=Vname(1,idApond)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idApond))
          Vinfo( 3)=Vname(3,idApond)
          Vinfo(14)=Vname(4,idApond)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idApond,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idApond), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface water (on ice) thickness.
!
        IF (Aout(idHpond,ng)) THEN
          Vinfo( 1)=Vname(1,idHpond)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idHpond))
          Vinfo( 3)=Vname(3,idHpond)
          Vinfo(14)=Vname(4,idHpond)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idHpond,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idHpond), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
!
!  Define ice age.
!
        IF (Aout(idAgeice,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idAgeice), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-ocean mass flux
!
        IF (Aout(idIomflx,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idIomflx), NF_FOUT,                &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice temperature.
!
        IF (Aout(idTimid,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTimid), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress component 11.
!
        IF (Aout(idSig11,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idSig11), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress component 12.
!
        IF (Aout(idSig11,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idSig12), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define internal ice stress component 22.
!
        IF (Aout(idSig22,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idSig22), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water friction velocity.
!
        IF (Aout(idTauiw,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idTauiw), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice-water momentum transfer coefficient.
!
        IF (Aout(idChuiw,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idChuiw), NF_FOUT,                 &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define temperature of molecular sublayer under ice.
!
        IF (Aout(idT0mk,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idT0mk), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define salinity of molecular sublayer under ice.
!
        IF (Aout(idS0mk,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idS0mk), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define frazil ice growth rate.
!
        IF (Aout(idWfr,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idWfr), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate.
!
        IF (Aout(idWai,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idWai), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate.
!
        IF (Aout(idWao,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idWao), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice growth/melt rate.
!
        IF (Aout(idWio,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idWio), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice melt runoff rate.
!
        IF (Aout(idWro,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idWro), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define ice divergence rate.
!
        IF (Aout(idWdiv,ng)) THEN
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
          status=def_var(ng, model, AVG(ng)%ncid,                       &
     &                   AVG(ng)%Vid(idWdiv), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WEC_MELLOR
!
!  Define 2D radiation stress Sxx.
!
        IF (Aout(idW2xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xx)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW2xx))
          Vinfo( 3)=Vname(3,idW2xx)
          Vinfo(14)=Vname(4,idW2xx)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW2xx),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress Sxy.
!
        IF (Aout(idW2xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xy)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW2xy))
          Vinfo( 3)=Vname(3,idW2xy)
          Vinfo(14)=Vname(4,idW2xy)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW2xy),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress Syy.
!
        IF (Aout(idW2yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2yy)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW2yy))
          Vinfo( 3)=Vname(3,idW2yy)
          Vinfo(14)=Vname(4,idW2yy)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW2yy),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WEC
!
!  Define 2D total u-wec stress.
!
        IF (Aout(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU2rs))
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU2rs),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D total v-wec stress.
!
        IF (Aout(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV2rs))
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV2rs),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D u-stokes velocity.
!
        IF (Aout(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU2Sd))
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU2Sd),  &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D v-stokes velocity.
!
        IF (Aout(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV2Sd))
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV2Sd),  &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WEC_VF
!
!  Define 2D quasi-static seal level adjustment.
!
        IF (Aout(idWztw,ng)) THEN
          Vinfo( 1)=Vname(1,idWztw)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWztw))
          Vinfo( 3)=Vname(3,idWztw)
          Vinfo(14)=Vname(4,idWztw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWztw,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWztw),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D quasi-static pressure.
!
        IF (Aout(idWqsp,ng)) THEN
          Vinfo( 1)=Vname(1,idWqsp)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWqsp))
          Vinfo( 3)=Vname(3,idWqsp)
          Vinfo(14)=Vname(4,idWqsp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWqsp,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWqsp),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D Bernoulli head.
!
        IF (Aout(idWbeh,ng)) THEN
          Vinfo( 1)=Vname(1,idWbeh)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWbeh))
          Vinfo( 3)=Vname(3,idWbeh)
          Vinfo(14)=Vname(4,idWbeh)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbeh,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWbeh),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WEC_MELLOR
#  ifdef SOLVE3D
!
!  Define 3D radiation stress Sxx.
!
        IF (Aout(idW3xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xx)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3xx))
          Vinfo( 3)=Vname(3,idW3xx)
          Vinfo(14)=Vname(4,idW3xx)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3xx),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Sxy.
!
        IF (Aout(idW3xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xy)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3xy))
          Vinfo( 3)=Vname(3,idW3xy)
          Vinfo(14)=Vname(4,idW3xy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3xy),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Syy.
!
        IF (Aout(idW3yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3yy)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3yy))
          Vinfo( 3)=Vname(3,idW3yy)
          Vinfo(14)=Vname(4,idW3yy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3yy),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Szx.
!
        IF (Aout(idW3zx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zx)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3zx))
          Vinfo( 3)=Vname(3,idW3zx)
          Vinfo(14)=Vname(4,idW3zx)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3zx),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Szy.
!
        IF (Aout(idW3zy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zy)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3zy))
          Vinfo( 3)=Vname(3,idW3zy)
          Vinfo(14)=Vname(4,idW3zy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3zy),  &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
# ifdef WEC
#  ifdef SOLVE3D
!
!  Define 3D total u-wec stress.
!
        IF (Aout(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU3rs))
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU3rs),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D total v-wec stress.
!
        IF (Aout(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV3rs))
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV3rs),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D u-stokes velocity.
!
        IF (Aout(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idU3Sd))
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idU3Sd),  &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D v-stokes velocity.
!
        IF (Aout(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idV3Sd))
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idV3Sd),  &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D Omega-stokes velocity.
!
        IF (Aout(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3Sd))
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3Sd,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3Sd),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D w-stokes velocity.
!
        IF (Aout(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idW3St))
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3St,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idW3St),  &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
# ifdef WAVES_HEIGHT
        IF (Aout(idWamp,ng)) THEN
          Vinfo( 1)=Vname(1,idWamp)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWamp))
          Vinfo( 3)=Vname(3,idWamp)
          Vinfo(14)=Vname(4,idWamp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWamp,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWamp),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (Aout(idWlen,ng)) THEN
          Vinfo( 1)=Vname(1,idWlen)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWlen))
          Vinfo( 3)=Vname(3,idWlen)
          Vinfo(14)=Vname(4,idWlen)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlen,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWlen),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_LENGTHP
        IF (Aout(idWlep,ng)) THEN
          Vinfo( 1)=Vname(1,idWlep)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWlep))
          Vinfo( 3)=Vname(3,idWlep)
          Vinfo(14)=Vname(4,idWlep)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWlep,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWlep),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_DIR
        IF (Aout(idWdir,ng)) THEN
          Vinfo( 1)=Vname(1,idWdir)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdir))
          Vinfo( 3)=Vname(3,idWdir)
          Vinfo(14)=Vname(4,idWdir)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdir,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWdir),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (Aout(idWptp,ng)) THEN
          Vinfo( 1)=Vname(1,idWptp)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWptp))
          Vinfo( 3)=Vname(3,idWptp)
          Vinfo(14)=Vname(4,idWptp)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWptp,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWptp),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (Aout(idWpbt,ng)) THEN
          Vinfo( 1)=Vname(1,idWpbt)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWpbt))
          Vinfo( 3)=Vname(3,idWpbt)
          Vinfo(14)=Vname(4,idWpbt)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWpbt,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWpbt),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WAVES_OCEAN
        IF (Aout(idWorb,ng)) THEN
          Vinfo( 1)=Vname(1,idWorb)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWorb))
          Vinfo( 3)=Vname(3,idWorb)
          Vinfo(14)=Vname(4,idWorb)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWorb,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWorb),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (Aout(idWdif,ng)) THEN
          Vinfo( 1)=Vname(1,idWdif)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdif))
          Vinfo( 3)=Vname(3,idWdif)
          Vinfo(14)=Vname(4,idWdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdif,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWdif),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (Aout(idWdib,ng)) THEN
          Vinfo( 1)=Vname(1,idWdib)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdib))
          Vinfo( 3)=Vname(3,idWdib)
          Vinfo(14)=Vname(4,idWdib)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdib,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWdib),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (Aout(idWdiw,ng)) THEN
          Vinfo( 1)=Vname(1,idWdiw)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdiw))
          Vinfo( 3)=Vname(3,idWdiw)
          Vinfo(14)=Vname(4,idWdiw)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdiw,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWdiw),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef ROLLER_SVENDSEN
        IF (Aout(idWbrk,ng)) THEN
          Vinfo( 1)=Vname(1,idWbrk)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWbrk))
          Vinfo( 3)=Vname(3,idWbrk)
          Vinfo(14)=Vname(4,idWbrk)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWbrk,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWbrk),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef WEC_ROLLER
        IF (Aout(idWdis,ng)) THEN
          Vinfo( 1)=Vname(1,idWdis)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWdis))
          Vinfo( 3)=Vname(3,idWdis)
          Vinfo(14)=Vname(4,idWdis)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWdis,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWdis),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (Aout(idWrol,ng)) THEN
          Vinfo( 1)=Vname(1,idWrol)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWrol))
          Vinfo( 3)=Vname(3,idWrol)
          Vinfo(14)=Vname(4,idWrol)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWrol,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idWrol),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef UV_KIRBY
        IF (Aout(idUwav,ng)) THEN
          Vinfo( 1)=Vname(1,idUwav)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUwav))
          Vinfo( 3)=Vname(3,idUwav)
          Vinfo(14)=Vname(4,idUwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUwav,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idUwav),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
        IF (Aout(idVwav,ng)) THEN
          Vinfo( 1)=Vname(1,idVwav)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVwav))
          Vinfo( 3)=Vname(3,idVwav)
          Vinfo(14)=Vname(4,idVwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVwav,ng),r8)
          status=def_var(ng, model, AVG(ng)%ncid, AVG(ng)%Vid(idVwav),  &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, AVG(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, AVG(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
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
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, model, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open averages file for read/write.
!
        CALL netcdf_open (ng, model, ncname, 1, AVG(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
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
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            AVG(ng)%Vid(idSdif)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSSSf))) THEN
            got_var(idSSSf)=.TRUE.
            AVG(ng)%Vid(idSSSf)=var_id(i)
#  ifdef ICE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.true.
            AVG(ng)%Vid(idUice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.true.
            AVG(ng)%Vid(idVice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUiceE))) THEN
            got_var(idUiceE)=.true.
            AVG(ng)%Vid(idUiceE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idViceN))) THEN
            got_var(idViceN)=.true.
            AVG(ng)%Vid(idViceN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.true.
            AVG(ng)%Vid(idAice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.true.
            AVG(ng)%Vid(idHice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.true.
            AVG(ng)%Vid(idTice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.true.
            AVG(ng)%Vid(idHsno)=var_id(i)
#ifdef MELT_PONDS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idApond))) THEN
            got_var(idApond)=.true.
            AVG(ng)%Vid(idApond)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHpond))) THEN
            got_var(idHpond)=.true.
            AVG(ng)%Vid(idHpond)=var_id(i)
#endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.true.
            AVG(ng)%Vid(idAgeice)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
            got_var(idIomflx)=.true.
            AVG(ng)%Vid(idIomflx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.true.
            AVG(ng)%Vid(idTimid)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.true.
            AVG(ng)%Vid(idSig11)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.true.
            AVG(ng)%Vid(idSig12)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.true.
            AVG(ng)%Vid(idSig22)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.true.
            AVG(ng)%Vid(idTauiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.true.
            AVG(ng)%Vid(idChuiw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.true.
            AVG(ng)%Vid(idT0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.true.
            AVG(ng)%Vid(idS0mk)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
            got_var(idWfr)=.true.
            AVG(ng)%Vid(idWfr)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
            got_var(idWai)=.true.
            AVG(ng)%Vid(idWai)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
            got_var(idWao)=.true.
            AVG(ng)%Vid(idWao)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
            got_var(idWio)=.true.
            AVG(ng)%Vid(idWio)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
            got_var(idWro)=.true.
            AVG(ng)%Vid(idWro)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiv))) THEN
            got_var(idWdiv)=.true.
            AVG(ng)%Vid(idWdiv)=var_id(i)
#  endif
#  if defined RUNOFF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idRunoff))) THEN
            got_var(idRunoff)=.TRUE.
            AVG(ng)%Vid(idRunoff)=var_id(i)
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
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUairE))) THEN
            got_var(idUairE)=.TRUE.
            AVG(ng)%Vid(idUairE)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVairN))) THEN
            got_var(idVairN)=.TRUE.
            AVG(ng)%Vid(idVairN)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            AVG(ng)%Vid(idTsur(itemp))=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            AVG(ng)%Vid(idTsur(isalt))=var_id(i)
#  ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            AVG(ng)%Vid(idLhea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            AVG(ng)%Vid(idShea)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            AVG(ng)%Vid(idLrad)=var_id(i)
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
# ifdef WEC_MELLOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xx))) THEN
            got_var(idW2xx)=.TRUE.
            AVG(ng)%Vid(idW2xx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xy))) THEN
            got_var(idW2xy)=.TRUE.
            AVG(ng)%Vid(idW2xy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2yy))) THEN
            got_var(idW2yy)=.TRUE.
            AVG(ng)%Vid(idW2yy)=var_id(i)
# endif
# ifdef WEC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            AVG(ng)%Vid(idU2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            AVG(ng)%Vid(idV2rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            AVG(ng)%Vid(idU2Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            AVG(ng)%Vid(idV2Sd)=var_id(i)
# endif
# ifdef WEC_MELLOR
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xx))) THEN
            got_var(idW3xx)=.TRUE.
            AVG(ng)%Vid(idW3xx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xy))) THEN
            got_var(idW3xy)=.TRUE.
            AVG(ng)%Vid(idW3xy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3yy))) THEN
            got_var(idW3yy)=.TRUE.
            AVG(ng)%Vid(idW3yy)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zx))) THEN
            got_var(idW3zx)=.TRUE.
            AVG(ng)%Vid(idW3zx)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zy))) THEN
            got_var(idW3zy)=.TRUE.
            AVG(ng)%Vid(idW3zy)=var_id(i)
#  endif
# endif
# ifdef WEC
#  ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            AVG(ng)%Vid(idU3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            AVG(ng)%Vid(idV3rs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            AVG(ng)%Vid(idU3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            AVG(ng)%Vid(idV3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            AVG(ng)%Vid(idW3Sd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            AVG(ng)%Vid(idW3St)=var_id(i)
#  endif
# endif
# ifdef WEC_VF
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWztw))) THEN
            got_var(idWztw)=.TRUE.
            AVG(ng)%Vid(idWztw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWqsp))) THEN
            got_var(idWqsp)=.TRUE.
            AVG(ng)%Vid(idWqsp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbeh))) THEN
            got_var(idWbeh)=.TRUE.
            AVG(ng)%Vid(idWbeh)=var_id(i)
# endif
# ifdef WAVES_HEIGHT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWamp))) THEN
            got_var(idWamp)=.TRUE.
            AVG(ng)%Vid(idWamp)=var_id(i)
# endif
# ifdef WAVES_LENGTH
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlen))) THEN
            got_var(idWlen)=.TRUE.
            AVG(ng)%Vid(idWlen)=var_id(i)
# endif
# ifdef WAVES_LENGTHP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWlep))) THEN
            got_var(idWlep)=.TRUE.
            AVG(ng)%Vid(idWlep)=var_id(i)
# endif
# ifdef WAVES_DIR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdir))) THEN
            got_var(idWdir)=.TRUE.
            AVG(ng)%Vid(idWdir)=var_id(i)
# endif
# ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWptp))) THEN
            got_var(idWptp)=.TRUE.
            AVG(ng)%Vid(idWptp)=var_id(i)
# endif
# ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWpbt))) THEN
            got_var(idWpbt)=.TRUE.
            AVG(ng)%Vid(idWpbt)=var_id(i)
# endif
# ifdef WAVES_OCEAN
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWorb))) THEN
            got_var(idWorb)=.TRUE.
            AVG(ng)%Vid(idWorb)=var_id(i)
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdif))) THEN
            got_var(idWdif)=.TRUE.
            AVG(ng)%Vid(idWdif)=var_id(i)
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdib))) THEN
            got_var(idWdib)=.TRUE.
            AVG(ng)%Vid(idWdib)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdiw))) THEN
            got_var(idWdiw)=.TRUE.
            AVG(ng)%Vid(idWdiw)=var_id(i)
# endif
# ifdef ROLLER_SVENDSEN
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWbrk))) THEN
            got_var(idWbrk)=.TRUE.
            AVG(ng)%Vid(idWbrk)=var_id(i)
# endif
# ifdef WEC_ROLLER
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWdis))) THEN
            got_var(idWdis)=.TRUE.
            AVG(ng)%Vid(idWdis)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWrol))) THEN
            got_var(idWrol)=.TRUE.
            AVG(ng)%Vid(idWrol)=var_id(i)
# endif
# ifdef UV_KIRBY
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUwav))) THEN
            got_var(idUwav)=.TRUE.
            AVG(ng)%Vid(idUwav)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVwav))) THEN
            got_var(idVwav)=.TRUE.
            AVG(ng)%Vid(idVwav)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             AVG(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
          DO itrc=1,NAT
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
#  if defined SEDIMENT && defined BEDLOAD
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              AVG(ng)%Vid(idUbld(itrc))=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              AVG(ng)%Vid(idVbld(itrc))=var_id(i)
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
        IF (.not.got_var(idSdif).and.Aout(idSdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idSSSf).and.Aout(idSSSf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSSSf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined RUNOFF
        IF (.not.got_var(idRunoff).and.Aout(idRunoff,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idRunoff)),        &
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
#  ifdef ICE_MODEL
        IF (.not.got_var(idUice).and.Aout(idUice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Aout(idVice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUiceE).and.Aout(idUiceE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUiceE)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idViceN).and.Aout(idViceN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idViceN)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Aout(idAice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Aout(idHice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Aout(idTice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Aout(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef MELT_PONDS
        IF (.not.got_var(idApond).and.Aout(idApond,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idApond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHpond).and.Aout(idHpond,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHpond)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
        IF (.not.got_var(idAgeice).and.Aout(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Aout(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Aout(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Aout(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Aout(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Aout(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Aout(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Aout(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Aout(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Aout(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Aout(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Aout(idWai,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Aout(idWao,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Aout(idWio,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Aout(idWro,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiv).and.Aout(idWdiv,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdiv)),          &
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
        IF (.not.got_var(idUairE).and.Aout(idUairE,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUairE)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVairN).and.Aout(idVairN,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVairN)),         &
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
        IF (.not.got_var(idTsur(isalt)).and.Aout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef BULK_FLUXES
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
# ifdef BBL
        IF (.not.got_var(idUbrs).and.Aout(idUbrs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbrs).and.Aout(idVbrs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbrs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbws).and.Aout(idUbws,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbws).and.Aout(idVbws,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbws)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbcs).and.Aout(idUbcs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbcs).and.Aout(idVbcs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbcs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVwc).and.Aout(idUVwc,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUVwc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        IF (.not.got_var(idUbot).and.Aout(idUbot,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbot).and.Aout(idVbot,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbot)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbur).and.Aout(idUbur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbvr).and.Aout(idVbvr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbvr)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WEC_MELLOR
        IF (.not.got_var(idW2xx).and.Aout(idW2xx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2xy).and.Aout(idW2xy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2yy).and.Aout(idW2yy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WEC
        IF (.not.got_var(idU2rs).and.Aout(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Aout(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2Sd).and.Aout(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Aout(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WEC_MELLOR
#  ifdef SOLVE3D
        IF (.not.got_var(idW3xx).and.Aout(idW3xx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3xy).and.Aout(idW3xy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3yy).and.Aout(idW3yy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zx).and.Aout(idW3zx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3zx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zy).and.Aout(idW3zy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3zy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef WEC
#  ifdef SOLVE3D
        IF (.not.got_var(idU3rs).and.Aout(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Aout(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Aout(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Aout(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.Aout(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.Aout(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef WEC_VF
        IF (.not.got_var(idWztw).and.Aout(idWztw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWztw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWqsp).and.Aout(idWqsp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWqsp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWbeh).and.Aout(idWbeh,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWbeh)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_HEIGHT
        IF (.not.got_var(idWamp).and.Aout(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_LENGTH
        IF (.not.got_var(idWamp).and.Aout(idWamp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWamp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
            AVG(ng)%Vid(idWlen)=var_id(i)
# endif
# ifdef WAVES_LENGTHP
        IF (.not.got_var(idWlep).and.Aout(idWlep,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWlep)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_DIR
        IF (.not.got_var(idWdir).and.Aout(idWdir,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdir)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_TOP_PERIOD
        IF (.not.got_var(idWptp).and.Aout(idWptp,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWptp)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_BOT_PERIOD
        IF (.not.got_var(idWpbt).and.Aout(idWpbt,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWpbt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WAVES_OCEAN
        IF (.not.got_var(idWorb).and.Aout(idWorb,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWorb)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
        IF (.not.got_var(idWdif).and.Aout(idWdif,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
        IF (.not.got_var(idWdib).and.Aout(idWdib,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdib)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWdiw).and.Aout(idWdiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdiw)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef ROLLER_SVENDSEN
        IF (.not.got_var(idWbrk).and.Aout(idWbrk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWbrk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef WEC_ROLLER
        IF (.not.got_var(idWdis).and.Aout(idWdis,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWdis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWrol).and.Aout(idWrol,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWrol)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef UV_KIRBY
        IF (.not.got_var(idUwav).and.Aout(idUwav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVwav).and.Aout(idVwav,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVwav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
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
#  if defined SEDIMENT && defined BEDLOAD
        DO i=1,NST
          IF (.not.got_var(idUbld(i)).and.Aout(idUbld(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.Aout(idVbld(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
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
# elif defined FILTERED
      AVGtime(ng)=time(ng)-0.5_r8*REAL(nfil,r8)*dt(ng)
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
  10  FORMAT (6x,'DEF_AVG   - creating average file, Grid ',i2.2,': ',  &
     &        a)
  20  FORMAT (6x,'DEF_AVG   - inquiring average file, Grid ',i2.2,': ', &
     &        a)
  30  FORMAT (/,' DEF_AVG - unable to create averages NetCDF file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_AVG - unable to open averages NetCDF file: ',a)
  60  FORMAT (/,' DEF_AVG - unable to find variable: ',a,2x,            &
     &        ' in averages NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_avg
#else
      SUBROUTINE def_avg
      RETURN
      END SUBROUTINE def_avg
#endif
