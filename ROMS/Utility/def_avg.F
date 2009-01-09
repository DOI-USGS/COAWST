#include "cppdefs.h"
#ifdef AVERAGES
      SUBROUTINE def_avg (ng,ldef)
!
!svn $Id: def_avg.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
# ifdef SEDIMENT
      USE mod_sediment
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
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
      integer, parameter :: Natt = 24

      logical :: got_var(NV)

      integer :: i, itrc, j, nrec, nvd, nvd3, nvd4
      integer :: recdim, status

      integer :: DimIDs(29), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

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
!
!=======================================================================
!  Create a new averages NetCDF file.
!=======================================================================
!
      IF (ldef.and.OutThread) THEN
        status=nf90_create(TRIM(ncname),nf90_clobber,ncAVGid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iNLM, ncAVGid(ng), 1)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
        status=nf90_def_dim(ncAVGid(ng),'xi_rho'  ,Lm(ng)+2,            &
     &                      DimIDs( 1))
        status=nf90_def_dim(ncAVGid(ng),'xi_u'    ,Lm(ng)+1,            &
     &                      DimIDs( 2))
        status=nf90_def_dim(ncAVGid(ng),'xi_v'    ,Lm(ng)+2,            &
     &                      DimIDs( 3))
        status=nf90_def_dim(ncAVGid(ng),'xi_psi'  ,Lm(ng)+1,            &
     &                      DimIDs( 4))
        status=nf90_def_dim(ncAVGid(ng),'eta_rho' ,Mm(ng)+2,            &
     &                      DimIDs( 5))
        status=nf90_def_dim(ncAVGid(ng),'eta_u'   ,Mm(ng)+2,            &
     &                      DimIDs( 6))
        status=nf90_def_dim(ncAVGid(ng),'eta_v'   ,Mm(ng)+1,            &
     &                      DimIDs( 7))
        status=nf90_def_dim(ncAVGid(ng),'eta_psi' ,Mm(ng)+1,            &
     &                      DimIDs( 8))
# if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncAVGid(ng),'xy_rho'  ,Nxyr(ng),            &
     &                      DimIDs(17))
        status=nf90_def_dim(ncAVGid(ng),'xy_u'    ,Nxyu(ng),            &
     &                      DimIDs(18))
        status=nf90_def_dim(ncAVGid(ng),'xy_v'    ,Nxyv(ng),            &
     &                      DimIDs(19))
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncAVGid(ng),'xyz_rho' ,Nxyr(ng)*N(ng),      &
     &                      DimIDs(20))
        status=nf90_def_dim(ncAVGid(ng),'xyz_u'   ,Nxyu(ng)*N(ng),      &
     &                      DimIDs(21))
        status=nf90_def_dim(ncAVGid(ng),'xyz_v'   ,Nxyv(ng)*N(ng),      &
     &                      DimIDs(22))
        status=nf90_def_dim(ncAVGid(ng),'xyz_w'   ,Nxyr(ng)*(N(ng)+1),  &
     &                      DimIDs(23))
#  endif
        status=nf90_def_dim(ncAVGid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncAVGid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncAVGid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
#  ifdef SEDIMENT
        status=nf90_def_dim(ncAVGid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#   if defined WRITE_WATER && defined MASKING
        status=nf90_def_dim(ncAVGid(ng),'xybed'   ,Nxyr(ng)*Nbed,       &
     &                      DimIDs(24))
#   endif
#  endif
#  ifdef ECOSIM
        status=nf90_def_dim(ncAVGid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncAVGid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncAVGid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncAVGid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
#  endif
# endif
        status=nf90_def_dim(ncAVGid(ng),'boundary',4,                   &
     &                      DimIDs(14))
# ifdef FOUR_DVAR
        status=nf90_def_dim(ncAVGid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
# endif
        status=nf90_def_dim(ncAVGid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
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
          DO j=1,80
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
        status=def_info(ng,ncAVGid(ng),ncname,DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
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
        status=def_var(ncAVGid(ng),avgVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idFsur,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idFsuD,ng),NF_FOUT,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idUbar,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idu2dD,ng),NF_FOUT,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idVbar,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idv2dD,ng),NF_FOUT,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idUfx1,ng),NF_FOUT,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idUfx2,ng),NF_FOUT,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idVfx1,ng),NF_FOUT,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idVfx2,ng),NF_FOUT,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idUvel,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idu3dD,ng),NF_FOUT,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idVvel,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idv3dD,ng),NF_FOUT,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idOvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idWvel,ng),NF_FOUT,         &
     &                   nvd4,w3dgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncAVGid(ng),avgTid(itrc,ng),NF_FOUT,         &
     &                     nvd4,t3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idDano,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idHsbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idHbbl,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idZZav,ng),NF_FOUT,           &
     &                 nvd3,t2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idU2av,ng),NF_FOUT,           &
     &                 nvd3,u2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idV2av,ng),NF_FOUT,           &
     &                 nvd3,v2dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idHUav,ng),NF_FOUT,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idHVav,ng),NF_FOUT,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idUUav,ng),NF_FOUT,           &
     &                 nvd4,u3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idUVav,ng),NF_FOUT,           &
     &                 nvd4,t3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idVVav,ng),NF_FOUT,           &
     &                 nvd4,v3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idTTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(iHUTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(iHVTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idUTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idVTav(itrc),ng),NF_FOUT,   &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idVvis,ng),NF_FOUT,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idTdif,ng),NF_FOUT,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
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
        status=def_var(ncAVGid(ng),avgVid(idSdif,ng),NF_FOUT,           &
     &                 nvd4,w3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idTsur(itemp),ng),NF_FOUT,  &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
          status=def_var(ncAVGid(ng),avgVid(idTsur(isalt),ng),NF_FOUT,  &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
        END IF

#   ifdef BULK_FLUXES
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
          status=def_var(ncAVGid(ng),avgVid(idLhea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idShea,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idLrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idTair,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idevap,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idrain,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idSrad,ng),NF_FOUT,         &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idUsms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idVsms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idUbms,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idVbms,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncAVGid(ng),avgVid(idUbld(i),ng),NF_FOUT,    &
     &                     nvd3,u2dgrd,Aval,Vinfo,ncname)
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
            status=def_var(ncAVGid(ng),avgVid(idVbld(i),ng),NF_FOUT,    &
     &                     nvd3,v2dgrd,Aval,Vinfo,ncname)
          END IF
        END DO
# endif
# ifdef AVERAGES_NEARSHORE
!
!  Define 2D radiation stress Sxx.
!
        IF (Hout(idW2xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW2xx)
          Vinfo( 2)=Vname(2,idW2xx)
          Vinfo( 3)=Vname(3,idW2xx)
          Vinfo(14)=Vname(4,idW2xx)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW2xx,ng),NF_FOUT,           &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW2xy,ng),NF_FOUT,           &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r2dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW2yy,ng),NF_FOUT,           &
     &                   nvd3,t2dgrd,Aval,Vinfo,ncname)
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
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u2dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idU2rs,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
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
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v2dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idV2rs,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
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
          status=def_var(ncAVGid(ng),avgVid(idU2Sd,ng),NF_FOUT,         &
     &                   nvd3,u2dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idV2Sd,ng),NF_FOUT,         &
     &                   nvd3,v2dgrd,Aval,Vinfo,ncname)
        END IF
#  ifdef SOLVE3D
!
!  Define 3D radiation stress Sxx.
!
        IF (Hout(idW3xx,ng)) THEN
          Vinfo( 1)=Vname(1,idW3xx)
          Vinfo( 2)=Vname(2,idW3xx)
          Vinfo( 3)=Vname(3,idW3xx)
          Vinfo(14)=Vname(4,idW3xx)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW3xx,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW3xy,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW3yy,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW3zx,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idW3zy,ng),NF_FOUT,         &
     &                   nvd4,t3dgrd,Aval,Vinfo,ncname)
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
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(u3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idU3rs,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
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
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(v3dvar,r8)
          status=def_var(ncAVGid(ng),avgVid(idV3rs,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
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
          status=def_var(ncAVGid(ng),avgVid(idU3Sd,ng),NF_FOUT,         &
     &                   nvd4,u3dgrd,Aval,Vinfo,ncname)
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
          status=def_var(ncAVGid(ng),avgVid(idV3Sd,ng),NF_FOUT,         &
     &                   nvd4,v3dgrd,Aval,Vinfo,ncname)
        END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncAVGid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iNLM, ncAVGid(ng), OutThread, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing averages file, check its contents, and prepare
!  for appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of averages NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open averages file for read/write.
!
        status=nf90_open(TRIM(ncname), nf90_write, ncAVGid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,50) TRIM(ncname)
          exit_flag=3
          ioerror=status
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
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idtime)),    &
     &                            avgVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idFsur)),    &
     &                            avgVid(idFsur,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUbar)),    &
     &                            avgVid(idUbar,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVbar)),    &
     &                            avgVid(idVbar,ng))
          END IF
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsuD))) THEN
            got_var(idFsuD)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idFsuD)),    &
     &                            avgVid(idFsuD,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idu2dD))) THEN
            got_var(idu2dD)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idu2dD)),    &
     &                            avgVid(idu2dD,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idv2dD))) THEN
            got_var(idv2dD)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idv2dD)),    &
     &                            avgVid(idv2dD,ng))
          END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx1))) THEN
            got_var(idUfx1)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUfx1)),    &
     &                            avgVid(idUfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUfx2))) THEN
            got_var(idUfx2)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUfx2)),    &
     &                            avgVid(idUfx2,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVfx1))) THEN
            got_var(idVfx1)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVfx1)),    &
     &                            avgVid(idVfx1,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVfx2))) THEN
            got_var(idVfx2)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVfx2)),    &
     &                            avgVid(idVfx2,ng))
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUvel)),    &
     &                            avgVid(idUvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVvel)),    &
     &                            avgVid(idVvel,ng))
          END IF
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idu3dD))) THEN
            got_var(idu3dD)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idu3dD)),    &
     &                            avgVid(idu3dD,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idv3dD))) THEN
            got_var(idv3dD)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idv3dD)),    &
     &                            avgVid(idv3dD,ng))
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOvel))) THEN
            got_var(idOvel)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idOvel)),    &
     &                            avgVid(idOvel,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idWvel))) THEN
            got_var(idWvel)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idWvel)),    &
     &                            avgVid(idWvel,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             status=nf90_inq_varid(ncAVGid(ng),                         &
     &                             TRIM(Vname(1,idTvar(itrc))),         &
     &                             avgTid(itrc,ng))
            END IF
          END DO
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idDano)),    &
     &                            avgVid(idDano,ng))
          END IF
#  ifdef LMD_SKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHsbl))) THEN
            got_var(idHsbl)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idHsbl)),    &
     &                            avgVid(idHsbl,ng))
          END IF
#  endif
#  ifdef LMD_BKPP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHbbl))) THEN
            got_var(idHbbl)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idHbbl)),    &
     &                            avgVid(idHbbl,ng))
          END IF
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idZZav))) THEN
            got_var(idZZav)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idZZav)),    &
     &                            avgVid(idZZav,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU2av))) THEN
            got_var(idU2av)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idU2av)),    &
     &                            avgVid(idU2av,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV2av))) THEN
            got_var(idV2av)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idV2av)),    &
     &                            avgVid(idV2av,ng))
          END IF
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHUav))) THEN
            got_var(idHUav)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idHUav)),    &
     &                            avgVid(idHUav,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idHVav))) THEN
            got_var(idHVav)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idHVav)),    &
     &                            avgVid(idHVav,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUUav))) THEN
            got_var(idUUav)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUUav)),    &
     &                            avgVid(idUUav,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUVav))) THEN
            got_var(idUVav)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUVav)),    &
     &                            avgVid(idUVav,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVVav))) THEN
            got_var(idVVav)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVVav)),    &
     &                            avgVid(idVVav,ng))
          END IF
          DO itrc=1,NAT
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,iHUTav(itrc)))) THEN
             got_var(iHUTav(itrc))=.TRUE.
             status=nf90_inq_varid(ncAVGid(ng),                         &
     &                             TRIM(Vname(1,iHUTav(itrc))),         &
     &                             avgVid(iHUTav(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,iHVTav(itrc)))) THEN
             got_var(iHVTav(itrc))=.TRUE.
             status=nf90_inq_varid(ncAVGid(ng),                         &
     &                             TRIM(Vname(1,iHVTav(itrc))),         &
     &                             avgVid(iHVTav(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUTav(itrc)))) THEN
             got_var(idUTav(itrc))=.TRUE.
             status=nf90_inq_varid(ncAVGid(ng),                         &
     &                             TRIM(Vname(1,idUTav(itrc))),         &
     &                             avgVid(idUTav(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVTav(itrc)))) THEN
             got_var(idVTav(itrc))=.TRUE.
             status=nf90_inq_varid(ncAVGid(ng),                         &
     &                             TRIM(Vname(1,idVTav(itrc))),         &
     &                             avgVid(idVTav(itrc),ng))
            END IF
          END DO
#  endif
#  ifdef AVERAGES_AKV
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVvis)),    &
     &                            avgVid(idVvis,ng))
          END IF
#  endif
#  ifdef AVERAGES_AKT
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
           status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idTdif)),     &
     &                           avgVid(idTdif,ng))
          END IF
#  endif
#  ifdef AVERAGES_AKS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idSdif)),    &
     &                            avgVid(idSdif,ng))
          END IF
#  endif
#  ifdef AVERAGES_FLUXES
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itemp)))) THEN
            got_var(idTsur(itemp))=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),                          &
     &                            TRIM(Vname(1,idTsur(itemp))),         &
     &                            avgVid(idTsur(itemp),ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(isalt)))) THEN
            got_var(idTsur(isalt))=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),                          &
     &                            TRIM(Vname(1,idTsur(isalt))),         &
     &                            avgVid(idTsur(isalt),ng))
          END IF
#   ifdef BULK_FLUXES
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idLhea))) THEN
            got_var(idLhea)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idLhea)),    &
     &                            avgVid(idLhea,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idShea))) THEN
            got_var(idShea)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idShea)),    &
     &                            avgVid(idShea,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idLrad))) THEN
            got_var(idLrad)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idLrad)),    &
     &                            avgVid(idLrad,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTair))) THEN
            got_var(idTair)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idTair)),    &
     &                            avgVid(idTair,ng))
          END IF
#    ifdef EMINUSP
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idevap))) THEN
            got_var(idevap)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idevap)),    &
     &                            avgVid(idevap,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idrain))) THEN
            got_var(idrain)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idrain)),    &
     &                            avgVid(idrain,ng))
          END IF
#    endif
#   endif
#   ifdef SHORTWAVE
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSrad))) THEN
            got_var(idSrad)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idSrad)),    &
     &                            avgVid(idSrad,ng))
          END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUsms)),    &
     &                            avgVid(idUsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVsms)),    &
     &                            avgVid(idVsms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idUbms)),    &
     &                            avgVid(idUbms,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idVbms)),    &
     &                            avgVid(idVbms,ng))
          END IF
# endif
# if defined SEDIMENT && defined BEDLOAD
          DO itrc=1,NST
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
              got_var(idUbld(itrc))=.true.
              status=nf90_inq_varid(ncAVGid(ng),                        &
     &                              TRIM(Vname(1,idUbld(itrc))),        &
     &                              avgVid(idUbld(itrc),ng))
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbld(itrc)))) THEN
              got_var(idVbld(itrc))=.true.
              status=nf90_inq_varid(ncAVGid(ng),                        &
     &                              TRIM(Vname(1,idVbld(itrc))),        &
     &                              avgVid(idVbld(itrc),ng))
            END IF
          END DO
# endif
# ifdef AVERAGES_NEARSHORE
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW2xx))) THEN
            got_var(idW2xx)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW2xx)),    &
     &                            avgVid(idW2xx,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW2xy))) THEN
            got_var(idW2xy)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW2xy)),    &
     &                            avgVid(idW2xy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW2yy))) THEN
            got_var(idW2yy)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW2yy)),    &
     &                            avgVid(idW2yy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU2rs))) THEN
            got_var(idU2rs)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idU2rs)),    &
     &                            avgVid(idU2rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV2rs))) THEN
            got_var(idV2rs)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idV2rs)),    &
     &                            avgVid(idV2rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU2Sd))) THEN
            got_var(idU2Sd)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idU2Sd)),    &
     &                            avgVid(idU2Sd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV2Sd))) THEN
            got_var(idV2Sd)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idV2Sd)),    &
     &                            avgVid(idV2Sd,ng))
          END IF
#  ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3xx))) THEN
            got_var(idW3xx)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW3xx)),    &
     &                            avgVid(idW3xx,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3xy))) THEN
            got_var(idW3xy)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW3xy)),    &
     &                            avgVid(idW3xy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3yy))) THEN
            got_var(idW3yy)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW3yy)),    &
     &                            avgVid(idW3yy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3zx))) THEN
            got_var(idW3zx)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW3zx)),    &
     &                            avgVid(idW3zx,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idW3zy))) THEN
            got_var(idW3zy)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idW3zy)),    &
     &                            avgVid(idW3zy,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU3rs))) THEN
            got_var(idU3rs)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idU3rs)),    &
     &                            avgVid(idU3rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV3rs))) THEN
            got_var(idV3rs)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idV3rs)),    &
     &                            avgVid(idV3rs,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idU3Sd))) THEN
            got_var(idU3Sd)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idU3Sd)),    &
     &                            avgVid(idU3Sd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idV3Sd))) THEN
            got_var(idV3Sd)=.TRUE.
            status=nf90_inq_varid(ncAVGid(ng),TRIM(Vname(1,idV3Sd)),    &
     &                            avgVid(idV3Sd,ng))
          END IF
#  endif
# endif
        END DO
!
!  Check if averages variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idFsuD)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idFsuD)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idu2dD)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idu2dD)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv2dD)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idv2dD)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_WRITE
        IF (.not.got_var(idUfx1)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUfx1)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUfx2)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUfx2)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx1)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVfx1)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVfx2)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVfx2)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
        IF (.not.got_var(idu3dD)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idu3dD)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idv3dD)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idv3dD)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
        IF (.not.got_var(idOvel).and.Hout(idOvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idOvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idWvel).and.Hout(idWvel,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idWvel)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef LMD_SKPP
        IF (.not.got_var(idHsbl).and.Hout(idHsbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHsbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (.not.got_var(idHbbl).and.Hout(idHbbl,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHbbl)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
# endif
# ifdef AVERAGES_QUADRATIC
        IF (.not.got_var(idZZav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idZZav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2av)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU2av)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2av)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV2av)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_QUADRATIC
        IF (.not.got_var(idHUav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHUav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idHVav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idHVav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUUav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUUav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUVav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUVav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVVav)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVVav)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NAT
          IF (.not.got_var(iHUTav(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,iHUTav(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(iHVTav(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,iHVTav(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idUTav(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,idUTav(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVTav(itrc))) THEN
            WRITE (stdout,60) TRIM(Vname(1,idVTav(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
#  endif
#  ifdef AVERAGES_AKV
        IF (.not.got_var(idVvis)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVvis)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKT
        IF (.not.got_var(idTdif)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_AKS
        IF (.not.got_var(idSdif)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSdif)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif
#  ifdef AVERAGES_FLUXES
        IF (.not.got_var(idTsur(itemp)).and.Hout(idTsur(itemp),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(itemp))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTsur(isalt)).and.Hout(idTsur(isalt),ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTsur(isalt))), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   ifdef BULK_FLUXES
        IF (.not.got_var(idLhea).and.Hout(idLhea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLhea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idShea).and.Hout(idShea,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idShea)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idLrad).and.Hout(idLrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idLrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idTair).and.Hout(idTair,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idTair)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    ifdef EMINUSP
        IF (.not.got_var(idevap).and.Hout(idevap,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idevap)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idrain).and.Hout(idrain,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idrain)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#    endif
#   endif
#   ifdef SHORTWAVE
        IF (.not.got_var(idSrad).and.Hout(idSrad,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idSrad)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVsms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idUbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idVbms)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# if defined SEDIMENT && defined BEDLOAD
        DO i=1,NST
          IF (.not.got_var(idUbld(i)).and.Hout(idUbld(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idUbld(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idVbld(i)).and.Hout(idVbld(i),ng)) THEN
            WRITE (stdout,60) TRIM(Vname(1,idVbld(i))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
# ifdef AVERAGES_NEARSHORE
        IF (.not.got_var(idW2xx).and.Hout(idW2xx,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW2xx)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2xy).and.Hout(idW2xy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW2xy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW2yy).and.Hout(idW2yy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW2yy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2rs).and.Hout(idU2rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU2rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2rs).and.Hout(idV2rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV2rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU2Sd).and.Hout(idU2Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU2Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV2Sd).and.Hout(idV2Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV2Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idW3xx).and.Hout(idW3xx,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3xx)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3xy).and.Hout(idW3xy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3xy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3yy).and.Hout(idW3yy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3yy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zx).and.Hout(idW3zx,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3zx)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idW3zy).and.Hout(idW3zy,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idW3zy)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3rs).and.Hout(idU3rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU3rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3rs).and.Hout(idV3rs,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV3rs)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idU3Sd).and.Hout(idU3Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idU3Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idV3Sd).and.Hout(idV3Sd,ng)) THEN
          WRITE (stdout,60) TRIM(Vname(1,idV3Sd)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
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
          tAVGindx(ng)=tsize
        END IF
      END IF
!
!  Set initial average time.
!
      IF (ntsAVG(ng).eq.1) THEN
        AVGtime(ng)=time(ng)-0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      ELSE
        AVGtime(ng)=time(ng)+REAL(ntsAVG(ng),r8)*dt(ng)-                &
     &              0.5_r8*REAL(nAVG(ng),r8)*dt(ng)
      END IF
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
