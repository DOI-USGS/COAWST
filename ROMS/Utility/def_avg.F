#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE def_avg (ng, ldef)
!
!svn $Id: def_avg.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
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

      integer :: i, itrc, j, model, nvd3, nvd4
      integer :: recdim, status

      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_dim

# ifdef SOLVE3D
#  ifdef SEDIMENT
      integer :: b3dgrd(4)
#  endif
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      real(r8) :: Aval(6)

# ifdef ADJOINT
      character (len=21) :: Prefix
# else
      character (len=13) :: Prefix
# endif
      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname
!
      SourceFile='def_avg.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=AVGname(ng)
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
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
        CALL netcdf_create (ng, model, TRIM(ncname), ncAVGid(ng))
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
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_rho',        &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_u',          &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_v',          &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xi_psi',        &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_rho',       &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_u',         &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_v',         &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'eta_psi',       &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xy_rho',        &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xy_u',          &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xy_v',          &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xyz_rho',       &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xyz_u',         &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xyz_v',         &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xyz_w',         &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
#  endif

        status=def_dim(ng, model, ncAVGid(ng), ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, ncAVGid(ng), ncname, 'NST',            &
     &                 NST, DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'xybed',         &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, model, ncAVGid(ng), ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, model, ncAVGid(ng), ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, model, ncAVGid(ng), ncname,                  &
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
!  Define dimension vectors for staggered u-momemtum type variables.
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
!  Define dimension vectors for staggered v-momemtum type variables.
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
!  Define dimension vector for staggered w-momemtum type variables.
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
        tAVGindx(ng)=0
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
        CALL def_info (ng, model, ncAVGid(ng), ncname, DimIDs)
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idtime,ng),       &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)

        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idFsur))
          Vinfo( 3)=Vname(3,idFsur)
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idFsur,ng),     &
# ifdef WET_DRY
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname,    &
     &                   SetFillVal = .FALSE.)
# else
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
# endif
          IF (exit_flag.ne.NoError) RETURN
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided free-surface.
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idFsuD,ng),       &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUbar,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 2D momentum in the XI-direction.
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idu2dD,ng),       &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVbar,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 2D momentum in the ETA-direction.
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idv2dD,ng),       &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
!
!  Define time-averaged mass fluxes for 3D momentum coupling.
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idUfx1,ng),       &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idUfx2,ng),       &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idVfx1,ng),       &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idVfx2,ng),       &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUvel,ng),     &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 3D momentum component in the XI-direction.
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idu3dD,ng),       &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVvel,ng),     &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Define detided 3D momentum component in the ETA-direction.
!
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
        status=def_var(ng, model, ncAVGid(ng), avgVid(idv3dD,ng),       &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
!
!  Define S-coordinate vertical "omega" momentum component.
!
        IF (Hout(idOvel,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idOvel,ng),     &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define "true" vertical momentum component.
!
        IF (Hout(idWvel,ng)) THEN
          Vinfo( 1)=Vname(1,idWvel)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idWvel))
          Vinfo( 3)=Vname(3,idWvel)
          Vinfo(14)=Vname(4,idWvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idWvel,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idWvel,ng),     &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
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
            status=def_var(ng, model, ncAVGid(ng), avgTid(itrc,ng),     &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idDano,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        IF (Hout(idHsbl,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idHsbl,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        IF (Hout(idHbbl,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idHbbl,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
!
!  Define quadratic <zeta*zeta> term.
!
        Vinfo( 1)=Vname(1,idZZav)
        Vinfo( 2)=TRIM(Vname(2,idZZav))
        Vinfo( 3)=Vname(3,idZZav)
        Vinfo(14)=Vname(4,idZZav)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r2dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idZZav,ng),       &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define quadratic <ubar*ubar> term.
!
        Vinfo( 1)=Vname(1,idU2av)
        Vinfo( 2)=TRIM(Vname(2,idU2av))
        Vinfo( 3)=Vname(3,idU2av)
        Vinfo(14)=Vname(4,idU2av)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u2dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idU2av,ng),       &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define quadratic <vbar*vbar> term.
!
        Vinfo( 1)=Vname(1,idV2av)
        Vinfo( 2)=TRIM(Vname(2,idV2av))
        Vinfo( 3)=Vname(3,idV2av)
        Vinfo(14)=Vname(4,idV2av)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v2dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idV2av,ng),       &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
!
!  Define u-volume flux.
!
        Vinfo( 1)=Vname(1,idHUav)
        Vinfo( 2)=TRIM(Vname(2,idHUav))
        Vinfo( 3)=Vname(3,idHUav)
        Vinfo(14)=Vname(4,idHUav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u3dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idHUav,ng),       &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define v-volume flux.
!
        Vinfo( 1)=Vname(1,idHVav)
        Vinfo( 2)=TRIM(Vname(2,idHVav))
        Vinfo( 3)=Vname(3,idHVav)
        Vinfo(14)=Vname(4,idHVav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v3dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idHVav,ng),       &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define quadratic <u*u> term.
!
        Vinfo( 1)=Vname(1,idUUav)
        Vinfo( 2)=TRIM(Vname(2,idUUav))
        Vinfo( 3)=Vname(3,idUUav)
        Vinfo(14)=Vname(4,idUUav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(u3dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idUUav,ng),       &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define quadratic <u*v> term.
!
        Vinfo( 1)=Vname(1,idUVav)
        Vinfo( 2)=TRIM(Vname(2,idUVav))
        Vinfo( 3)=Vname(3,idUVav)
        Vinfo(14)=Vname(4,idUVav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(r3dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idUVav,ng),       &
     &                 NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define quadratic <v*v> term.
!
        Vinfo( 1)=Vname(1,idVVav)
        Vinfo( 2)=TRIM(Vname(2,idVVav))
        Vinfo( 3)=Vname(3,idVVav)
        Vinfo(14)=Vname(4,idVVav)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(v3dvar,r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idVVav,ng),       &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define quadratic <t*t> terms.
!
        DO itrc=1,NAT
          Vinfo( 1)=Vname(1,idTTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,idTTav(itrc)))
          Vinfo( 3)=Vname(3,idTTav(itrc))
          Vinfo(14)=Vname(4,idTTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTTav(itrc),ng), NF_FOUT,              &
     &                   nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Define active tracers volume fluxes.
!
        DO itrc=1,NAT
          Vinfo( 1)=Vname(1,iHUTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,iHUTav(itrc)))
          Vinfo( 3)=Vname(3,iHUTav(itrc))
          Vinfo(14)=Vname(4,iHUTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(iHUTav(itrc),ng), NF_FOUT,              &
     &                   nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
          Vinfo( 1)=Vname(1,iHVTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,iHVTav(itrc)))
          Vinfo( 3)=Vname(3,iHVTav(itrc))
          Vinfo(14)=Vname(4,iHVTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(iHVTav(itrc),ng), NF_FOUT,              &
     &                   nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Define quadratic <u*t> and <v*t> terms.
!
        DO itrc=1,NAT
          Vinfo( 1)=Vname(1,idUTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,idUTav(itrc)))
          Vinfo( 3)=Vname(3,idUTav(itrc))
          Vinfo(14)=Vname(4,idUTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idUTav(itrc),ng), NF_FOUT,              &
     &                   nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
!
          Vinfo( 1)=Vname(1,idVTav(itrc))
          Vinfo( 2)=TRIM(Vname(2,idVTav(itrc)))
          Vinfo( 3)=Vname(3,idVTav(itrc))
          Vinfo(14)=Vname(4,idVTav(itrc))
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idVTav(itrc),ng), NF_FOUT,              &
     &                   nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#  endif
#  ifdef AVERAGES_AKV
!
!  Define vertical viscosity coefficient.
!
        Vinfo( 1)=Vname(1,idVvis)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVvis))
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idVvis,ng),       &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_AKT
!
!  Define vertical diffusion coefficient for potential temperature.
!
        Vinfo( 1)=Vname(1,idTdif)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTdif))
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idTdif,ng),       &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_AKS
!
!  Define vertical diffusion coefficient for salinity.
!
        Vinfo( 1)=Vname(1,idSdif)
        WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSdif))
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
        status=def_var(ng, model, ncAVGid(ng), avgVid(idSdif,ng),       &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Define surface net heat flux.
!
        IF (Hout(idTsur(itemp),ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTsur(itemp),ng), NF_FOUT,             &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface net salt flux.
!
        IF (Hout(idTsur(isalt),ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTsur(isalt),ng), NF_FOUT,             &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

#   if defined BULK_FLUXES || defined ATM2OCN_FLUXES
!
!  Define latent heat flux.
!
        IF (Hout(idLhea,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idLhea,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define sensible heat flux.
!
        IF (Hout(idShea,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idShea,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define longwave radiation flux.
!
        IF (Hout(idLrad,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idLrad,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define atmospheric air temperature.
!
        IF (Hout(idTair,ng)) THEN
          Vinfo( 1)=Vname(1,idTair)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idTair))
          Vinfo( 3)=Vname(3,idTair)
          Vinfo(14)=Vname(4,idTair)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTair,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idTair,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#   ifdef EMINUSP
!
!  Define evaporation rate.
!
        IF (Hout(idevap,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng), avgVid(idevap,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define precipitation rate.
!
        IF (Hout(idrain,ng)) THEN
          Vinfo( 1)=Vname(1,idrain)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idrain))
          Vinfo( 2)='averaged '//TRIM(Vname(2,idrain))
          Vinfo( 3)=Vname(3,idrain)
          Vinfo(12)='downward flux, freshening (precipitation)'
          Vinfo(14)=Vname(4,idrain)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idrain,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idrain,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#    endif
#   endif
#   ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        IF (Hout(idSrad,ng)) THEN
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
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idSrad,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Define surface u-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUsms,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface v-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVsms,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom u-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idUbms,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom v-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idVbms,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# if defined SEDIMENT && defined BEDLOAD
        DO i=1,NST
!
!  Define Bedload U-direction.
!
          IF (Hout(idUbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idUbld(i))
            Vinfo( 2)=Vname(2,idUbld(i))
            Vinfo( 3)=Vname(3,idUbld(i))
            Vinfo(14)=Vname(4,idUbld(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_u'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idUbld(i),ng),r8)
            status=def_var(ng, model, ncAVGid(ng),                      &
     &                     avgVid(idUbld(i),ng), NF_FOUT,               &
     &                     nvd3, u2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
!
!  Define Bedload V-direction.
!
          IF (Hout(idVbld(i),ng)) THEN
            Vinfo( 1)=Vname(1,idVbld(i))
            Vinfo( 2)=Vname(2,idVbld(i))
            Vinfo( 3)=Vname(3,idVbld(i))
            Vinfo(14)=Vname(4,idVbld(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_v'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idVbld(i),ng),r8)
            status=def_var(ng, model, ncAVGid(ng),                      &
     &                     avgVid(idVbld(i),ng), NF_FOUT,               &
     &                     nvd3, v2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
# endif
# ifdef WRT_SSSFLX
        IF (Hout(idSSSf,ng)) THEN
          Vinfo( 1)=Vname(1,idSSSf)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSSSf))
          Vinfo( 3)=Vname(3,idSSSf)
          Vinfo(14)=Vname(4,idSSSf)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSSSf,ng),r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idSSSf,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
# ifdef ICE_MODEL      
! 
!  Define 2D ice momentum in the XI-direction.
! 
        IF (Hout(idUice,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idUice,ng), NF_FOUT,                    &
     &                   nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
! 
!  Define 2D ice momentum in the ETA-direction.
!
        IF (Hout(idVice,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idVice,ng), NF_FOUT,                    &
     &                   nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF  
! 
!  Define ice concentration.
!  
        IF (Hout(idAice,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idAice,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!       
!  Define ice average thickness.
!  
        IF (Hout(idHice,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idHice,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice/snow surface temperature.
!  
        IF (Hout(idTice,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTice,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define snow thickness.
!  
        IF (Hout(idHsno,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idHsno,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define surface water (on ice) thickness.
! 
        IF (Hout(idSfwat,ng)) THEN
          Vinfo( 1)=Vname(1,idSfwat)
          WRITE (Vinfo( 2),'(a,1x,a)') Prefix, TRIM(Vname(2,idSfwat))
          Vinfo( 3)=Vname(3,idSfwat)
          Vinfo(14)=Vname(4,idSfwat)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif   
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSfwat,ng),r8)
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idSfwat,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice age.
! 
        IF (Hout(idAgeice,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idAgeice,ng), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice-ocean mass flux
!  
        IF (Hout(idIomflx,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idIomflx,ng), NF_FOUT,                  &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define internal ice temperature.
!  
        IF (Hout(idTimid,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTimid,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define internal ice stress component 11.
!  
        IF (Hout(idSig11,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idSig11,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define internal ice stress component 12.
!  
        IF (Hout(idSig11,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idSig12,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define internal ice stress component 22.
!  
        IF (Hout(idSig22,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idSig22,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice-water friction velocity.
!  
        IF (Hout(idTauiw,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idTauiw,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice-water momentum transfer coefficient.
!  
        IF (Hout(idChuiw,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idChuiw,ng), NF_FOUT,                   &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define temperature of molecular sublayer under ice.
!  
        IF (Hout(idT0mk,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idT0mk,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define salinity of molecular sublayer under ice.
!  
        IF (Hout(idS0mk,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idS0mk,ng), NF_FOUT,                    &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define frazil ice growth rate.
!  
        IF (Hout(idWfr,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idWfr,ng), NF_FOUT,                     &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice growth/melt rate.
!  
        IF (Hout(idWai,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idWai,ng), NF_FOUT,                     &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice growth/melt rate.
!  
        IF (Hout(idWao,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idWao,ng), NF_FOUT,                     &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice growth/melt rate.
!  
        IF (Hout(idWio,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idWio,ng), NF_FOUT,                     &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
!       
!  Define ice melt runoff rate.
!  
        IF (Hout(idWro,ng)) THEN
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
          status=def_var(ng, model, ncAVGid(ng),                        &
     &                   avgVid(idWro,ng), NF_FOUT,                     &
     &                   nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF           
# endif 
# ifdef AVERAGES_UV_KIRBY
!
!  Define coupling U-current given by Kirby and Chen (1989).
!
        IF (Hout(idUwav,ng)) THEN
          Vinfo( 1)=Vname(1,idUwav)
          Vinfo( 2)=Vname(2,idUwav)
          Vinfo( 3)=Vname(3,idUwav)
          Vinfo(14)=Vname(4,idUwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUwav,ng),r8)
          status=def_var(ng,model,ncAVGid(ng),avgVid(idUwav,ng),        &
     &                   NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF
!
!  Define coupling V-current given by Kirby and Chen (1989).
!
        IF (Hout(idVwav,ng)) THEN
          Vinfo( 1)=Vname(1,idVwav)
          Vinfo( 2)=Vname(2,idVwav)
          Vinfo( 3)=Vname(3,idVwav)
          Vinfo(14)=Vname(4,idVwav)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_r'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVwav,ng),r8)
          status=def_var(ng,model,ncAVGid(ng),avgVid(idVwav,ng),        &
     &                   NF_FOUT,nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF

# endif
# ifdef AVERAGES_WEC
#  ifdef WEC_MELLOR
!
!  Define 2D radiation stress Sxx.
!
        IF (Hout(idW2xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xx)
          Vinfo( 2)=Vname(2,idW2xx)
          Vinfo( 3)=Vname(3,idW2xx)
          Vinfo(14)=Vname(4,idW2xx)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW2xx,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress Sxy.
!
        IF (Hout(idW2xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xy)
          Vinfo( 2)=Vname(2,idW2xy)
          Vinfo( 3)=Vname(3,idW2xy)
          Vinfo(14)=Vname(4,idW2xy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW2xy,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D radiation stress Syy.
!
        IF (Hout(idW2yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW2yy)
          Vinfo( 2)=Vname(2,idW2yy)
          Vinfo( 3)=Vname(3,idW2yy)
          Vinfo(14)=Vname(4,idW2yy)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW2yy,ng),     &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D total u-radiation stress.
!
        IF (Hout(idU2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU2rs)
          Vinfo( 2)=Vname(2,idU2rs)
          Vinfo( 3)=Vname(3,idU2rs)
          Vinfo(14)=Vname(4,idU2rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idU2rs,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D total v-radiation stress.
!
        IF (Hout(idV2rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV2rs)
          Vinfo( 2)=Vname(2,idV2rs)
          Vinfo( 3)=Vname(3,idV2rs)
          Vinfo(14)=Vname(4,idV2rs)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idV2rs,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
!
!  Define 2D u-stokes velocity.
!
        IF (Hout(idU2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU2Sd)
          Vinfo( 2)=Vname(2,idU2Sd)
          Vinfo( 3)=Vname(3,idU2Sd)
          Vinfo(14)=Vname(4,idU2Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU2Sd,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idU2Sd,ng),     &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D v-stokes velocity.
!
        IF (Hout(idV2Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV2Sd)
          Vinfo( 2)=Vname(2,idV2Sd)
          Vinfo( 3)=Vname(3,idV2Sd)
          Vinfo(14)=Vname(4,idV2Sd)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV2Sd,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idV2Sd,ng),     &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef SOLVE3D
#   ifdef WEC_MELLOR
!
!  Define 3D radiation stress Sxx.
!
        IF (Hout(idW3xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xx)
          Vinfo( 2)=Vname(2,idW3xx)
          Vinfo( 3)=Vname(3,idW3xx)
          Vinfo(14)=Vname(4,idW3xx)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3xx,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Sxy.
!
        IF (Hout(idW3xy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xy)
          Vinfo( 2)=Vname(2,idW3xy)
          Vinfo( 3)=Vname(3,idW3xy)
          Vinfo(14)=Vname(4,idW3xy)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3xy,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Syy.
!
        IF (Hout(idW3yy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3yy)
          Vinfo( 2)=Vname(2,idW3yy)
          Vinfo( 3)=Vname(3,idW3yy)
          Vinfo(14)=Vname(4,idW3yy)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3yy,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Szx.
!
        IF (Hout(idW3zx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zx)
          Vinfo( 2)=Vname(2,idW3zx)
          Vinfo( 3)=Vname(3,idW3zx)
          Vinfo(14)=Vname(4,idW3zx)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3zx,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D radiation stress Szy.
!
        IF (Hout(idW3zy,ng)) THEN
          Vinfo( 1)=Vname(1,idW3zy)
          Vinfo( 2)=Vname(2,idW3zy)
          Vinfo( 3)=Vname(3,idW3zy)
          Vinfo(14)=Vname(4,idW3zy)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3zy,ng),     &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D total u-radiation stress.
!
        IF (Hout(idU3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idU3rs)
          Vinfo( 2)=Vname(2,idU3rs)
          Vinfo( 3)=Vname(3,idU3rs)
          Vinfo(14)=Vname(4,idU3rs)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idU3rs,ng),     &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D total v-radiation stress.
!
        IF (Hout(idV3rs,ng)) THEN
          Vinfo( 1)=Vname(1,idV3rs)
          Vinfo( 2)=Vname(2,idV3rs)
          Vinfo( 3)=Vname(3,idV3rs)
          Vinfo(14)=Vname(4,idV3rs)
          Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#    endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idV3rs,ng),     &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#   endif
!
!  Define 3D u-stokes velocity.
!
        IF (Hout(idU3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idU3Sd)
          Vinfo( 2)=Vname(2,idU3Sd)
          Vinfo( 3)=Vname(3,idU3Sd)
          Vinfo(14)=Vname(4,idU3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idU3Sd,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idU3Sd,ng),     &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D v-stokes velocity.
!
        IF (Hout(idV3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idV3Sd)
          Vinfo( 2)=Vname(2,idV3Sd)
          Vinfo( 3)=Vname(3,idV3Sd)
          Vinfo(14)=Vname(4,idV3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idV3Sd,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idV3Sd,ng),     &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D omega-stokes velocity.
!
        IF (Hout(idW3Sd,ng)) THEN
          Vinfo( 1)=Vname(1,idW3Sd)
          Vinfo( 2)=Vname(2,idW3Sd)
          Vinfo( 3)=Vname(3,idW3Sd)
          Vinfo(14)=Vname(4,idW3Sd)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3Sd,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3Sd,ng),     &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D w-stokes velocity test
!
        IF (Hout(idW3St,ng)) THEN
          Vinfo( 1)=Vname(1,idW3St)
          Vinfo( 2)=Vname(2,idW3St)
          Vinfo( 3)=Vname(3,idW3St)
          Vinfo(14)=Vname(4,idW3St)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idW3St,ng),r8)
          status=def_var(ng, model, ncAVGid(ng), avgVid(idW3St,ng),     &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, model, ncname, ncAVGid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, model, ncAVGid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=AVGname(ng)
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
        CALL netcdf_open (ng, model, ncname, 1, ncAVGid(ng))
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
            avgVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            avgVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            avgVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            avgVid(idVbar,ng)=var_id(i)
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsuD))) THEN
            got_var(idFsuD)=.TRUE.
            avgVid(idFsuD,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu2dD))) THEN
            got_var(idu2dD)=.TRUE.
            avgVid(idu2dD,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv2dD))) THEN
            got_var(idv2dD)=.TRUE.
            avgVid(idv2dD,ng)=var_id(i)
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            avgVid(idUfx1,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            avgVid(idUfx2,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            avgVid(idVfx1,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            avgVid(idVfx2,ng)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            avgVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            avgVid(idVvel,ng)=var_id(i)
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idu3dD))) THEN
            got_var(idu3dD)=.TRUE.
            avgVid(idu3dD,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idv3dD))) THEN
            got_var(idv3dD)=.TRUE.
            avgVid(idv3dD,ng)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            avgVid(idOvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            avgVid(idWvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            avgVid(idDano,ng)=var_id(i)
#  ifdef LMD_SKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            avgVid(idHsbl,ng)=var_id(i)
#  endif
#  ifdef LMD_BKPP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            avgVid(idHbbl,ng)=var_id(i)
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idZZav))) THEN
            got_var(idZZav)=.TRUE.
            avgVid(idZZav,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2av))) THEN
            got_var(idU2av)=.TRUE.
            avgVid(idU2av,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2av))) THEN
            got_var(idV2av)=.TRUE.
            avgVid(idV2av,ng)=var_id(i)
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHUav))) THEN
            got_var(idHUav)=.TRUE.
            avgVid(idHUav,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHVav))) THEN
            got_var(idHVav)=.TRUE.
            avgVid(idHVav,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUUav))) THEN
            got_var(idUUav)=.TRUE.
            avgVid(idUUav,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUVav))) THEN
            got_var(idUVav)=.TRUE.
            avgVid(idUVav,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVVav))) THEN
            got_var(idVVav)=.TRUE.
            avgVid(idVVav,ng)=var_id(i)
#  endif
#  ifdef AVERAGES_AKV
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            avgVid(idVvis,ng)=var_id(i)
#  endif
#  ifdef AVERAGES_AKT
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            avgVid(idTdif,ng)=var_id(i)
#  endif
#  ifdef AVERAGES_AKS
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            avgVid(idSdif,ng)=var_id(i)
#  endif
#  ifdef WRT_SSSFLX
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSSSf))) THEN
            got_var(idSSSf)=.TRUE.
            avgVid(idSSSf,ng)=var_id(i)
#  endif
#  ifdef ICE_MODEL
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUice))) THEN
            got_var(idUice)=.true.
            avgVid(idUice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVice))) THEN
            got_var(idVice)=.true.
            avgVid(idVice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAice))) THEN
            got_var(idAice)=.true.
            avgVid(idAice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHice))) THEN
            got_var(idHice)=.true.
            avgVid(idHice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTice))) THEN
            got_var(idTice)=.true.
            avgVid(idTice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsno))) THEN
            got_var(idHsno)=.true.
            avgVid(idHsno,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSfwat))) THEN
            got_var(idSfwat)=.true.
            avgVid(idSfwat,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idAgeice))) THEN
            got_var(idAgeice)=.true.
            avgVid(idAgeice,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idIomflx))) THEN
            got_var(idIomflx)=.true.
            avgVid(idIomflx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTimid))) THEN
            got_var(idTimid)=.true.
            avgVid(idTimid,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig11))) THEN
            got_var(idSig11)=.true.
            avgVid(idSig11,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig12))) THEN
            got_var(idSig12)=.true.
            avgVid(idSig12,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSig22))) THEN
            got_var(idSig22)=.true.
            avgVid(idSig22,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTauiw))) THEN
            got_var(idTauiw)=.true.
            avgVid(idTauiw,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idChuiw))) THEN
            got_var(idChuiw)=.true.
            avgVid(idChuiw,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idT0mk))) THEN
            got_var(idT0mk)=.true.
            avgVid(idT0mk,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idS0mk))) THEN
            got_var(idS0mk)=.true.
            avgVid(idS0mk,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWfr))) THEN
            got_var(idWfr)=.true.
            avgVid(idWfr,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWai))) THEN
            got_var(idWai)=.true.
            avgVid(idWai,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWao))) THEN
            got_var(idWao)=.true.
            avgVid(idWao,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWio))) THEN
            got_var(idWio)=.true.
            avgVid(idWio,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idWro))) THEN
            got_var(idWro)=.true.
            avgVid(idWro,ng)=var_id(i)
#  endif
#  ifdef AVERAGES_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            avgVid(idTsur(itemp),ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.                                &
     &             TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            avgVid(idTsur(isalt),ng)=var_id(i)
#   ifdef BULK_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            avgVid(idLhea,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            avgVid(idShea,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            avgVid(idLrad,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            avgVid(idTair,ng)=var_id(i)
#    ifdef EMINUSP
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            avgVid(idevap,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            avgVid(idrain,ng)=var_id(i)
#    endif
#   endif
#   ifdef SHORTWAVE
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            avgVid(idSrad,ng)=var_id(i)
#   endif
#  endif
# endif

# ifdef AVERAGES_FLUXES
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            avgVid(idUsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            avgVid(idVsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            avgVid(idUbms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            avgVid(idVbms,ng)=var_id(i)
# endif
# ifdef AVERAGES_WEC
#  ifdef WEC_MELLOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xx))) THEN
            got_var(idW2xx)=.TRUE.
            avgVid(idW2xx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2xy))) THEN
            got_var(idW2xy)=.TRUE.
            avgVid(idW2xy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW2yy))) THEN
            got_var(idW2yy)=.TRUE.
            avgVid(idW2yy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            avgVid(idU2rs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            avgVid(idV2rs,ng)=var_id(i)
#  endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            avgVid(idU2Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            avgVid(idV2Sd,ng)=var_id(i)
#  ifdef SOLVE3D
#   ifdef WEC_MELLOR
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xx))) THEN
            got_var(idW3xx)=.TRUE.
            avgVid(idW3xx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3xy))) THEN
            got_var(idW3xy)=.TRUE.
            avgVid(idW3xy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3yy))) THEN
            got_var(idW3yy)=.TRUE.
            avgVid(idW3yy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zx))) THEN
            got_var(idW3zx)=.TRUE.
            avgVid(idW3zx,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3zy))) THEN
            got_var(idW3zy)=.TRUE.
            avgVid(idW3zy,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            avgVid(idU3rs,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            avgVid(idV3rs,ng)=var_id(i)
#   endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            avgVid(idU3Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            avgVid(idV3Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3Sd))) THEN
            got_var(idW3Sd)=.TRUE.
            avgVid(idW3Sd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idW3St))) THEN
            got_var(idW3St)=.TRUE.
            avgVid(idW3St,ng)=var_id(i)
#  endif
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             avgTid(itrc,ng)=var_id(i)
            END IF
          END DO
#  ifdef AVERAGES_QUADRATIC
          DO itrc=1,NAT
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
              got_var(iHUTav(itrc))=.TRUE.
              avgVid(iHUTav(itrc),ng)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,iHVTav(itrc)))) THEN
              got_var(iHVTav(itrc))=.TRUE.
              avgVid(iHVTav(itrc),ng)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idUTav(itrc)))) THEN
              got_var(idUTav(itrc))=.TRUE.
              avgVid(idUTav(itrc),ng)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVTav(itrc)))) THEN
              got_var(idVTav(itrc))=.TRUE.
              avgVid(idVTav(itrc),ng)=var_id(i)
            END IF
          END DO
#  endif
#  if defined SEDIMENT && defined BEDLOAD
          DO itrc=1,NST
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              avgVid(idUbld(itrc),ng)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              avgVid(idVbld(itrc),ng)=var_id(i)
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
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idFsuD)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFsuD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dD)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu2dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dD)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv2dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
        IF (.not.got_var(idUfx1)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVfx1)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVfx2)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idu3dD)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idu3dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dD)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idv3dD)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idOvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHbbl)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
        IF (.not.got_var(idZZav)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idZZav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2av)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2av)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2av)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
        IF (.not.got_var(idHUav)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHVav)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUUav)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUUav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVVav)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKV
        IF (.not.got_var(idVvis)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVvis)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKT
        IF (.not.got_var(idTdif)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKS
        IF (.not.got_var(idSdif)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSdif)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef WRT_SSSFLX
        IF (.not.got_var(idSSSf).and.Hout(idSSSf,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSSSf)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# ifdef ICE_MODEL
        IF (.not.got_var(idUice).and.Hout(idUice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVice).and.Hout(idVice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAice).and.Hout(idAice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHice).and.Hout(idHice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTice).and.Hout(idTice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTice)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHsno).and.Hout(idHsno,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idHsno)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSfwat).and.Hout(idSfwat,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSfwat)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idAgeice).and.Hout(idAgeice,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idAgeice)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idIomflx).and.Hout(idIomflx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idIomflx)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTimid).and.Hout(idTimid,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTimid)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig11).and.Hout(idSig11,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig11)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig12).and.Hout(idSig12,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig12)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idSig22).and.Hout(idSig22,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSig22)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTauiw).and.Hout(idTauiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTauiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idChuiw).and.Hout(idChuiw,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idChuiw)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idT0mk).and.Hout(idT0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idT0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idS0mk).and.Hout(idS0mk,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idS0mk)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWfr).and.Hout(idWfr,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWfr)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWai).and.Hout(idWai,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWai)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWao).and.Hout(idWao,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWao)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWio).and.Hout(idWio,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWio)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWro).and.Hout(idWro,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idWro)),           &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#  ifdef AVERAGES_FLUXES
        IF (.not.got_var(idTsur(itemp)).and.Hout(idTsur(itemp),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Hout(idTsur(isalt),ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))),   &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLhea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idShea)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idLrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTair)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef EMINUSP
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idevap)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idrain)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#   ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idSrad)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef AVERAGES_WEC
#  ifdef WEC_MELLOR
        IF (.not.got_var(idW2xx).and.Hout(idW2xx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2xy).and.Hout(idW2xy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2yy).and.Hout(idW2yy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW2yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Hout(idU2rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Hout(idV2rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idU2Sd).and.Hout(idU2Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Hout(idV2Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV2Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
#   ifdef WEC_MELLOR
        IF (.not.got_var(idW3xx).and.Hout(idW3xx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3xx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3xy).and.Hout(idW3xy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3xy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3yy).and.Hout(idW3yy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3yy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zx).and.Hout(idW3zx,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3zx)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zy).and.Hout(idW3zy,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3zy)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3rs).and.Hout(idU3rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Hout(idV3rs,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV3rs)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
        IF (.not.got_var(idU3Sd).and.Hout(idU3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idU3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Hout(idV3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idV3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3Sd).and.Hout(idW3Sd,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3Sd)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3St).and.Hout(idW3St,ng)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idW3St)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  ifdef AVERAGES_QUADRATIC
        DO itrc=1,NAT
          IF (.not.got_var(iHUTav(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  if defined SEDIMENT && defined BEDLOAD
        DO i=1,NST
          IF (.not.got_var(idUbld(i)).and.Hout(idUbld(i),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idUbld(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.Hout(idVbld(i),ng)) THEN
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
            tAVGindx(ng)=((ntstart(ng)-1)-                              &
     &                    ndefAVG(ng)*((ntstart(ng)-1)/ndefAVG(ng)))/   &
     &                   nAVG(ng)
          ELSE
            tAVGindx(ng)=(ntstart(ng)-1)/nAVG(ng)
          END IF
        ELSE
          tAVGindx(ng)=rec_size
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
  10  FORMAT (6x,'DEF_AVG   - creating average file: ',a)
  20  FORMAT (6x,'DEF_AVG   - inquiring average file: ',a)
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
