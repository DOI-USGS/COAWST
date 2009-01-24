#include "cppdefs.h"
#ifdef ADJOINT
      SUBROUTINE ad_def_his (ng, ldef)
!
!svn $Id: ad_def_his.F 677 2008-08-05 20:17:30Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates adjoint history NetCDF file, it defines its    !
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
# if defined SEDIMENT_NOT_YET || defined BBL_MODEL_NOT_YET
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
      logical :: got_var(NV)

      integer, parameter :: Natt = 24

      integer :: i, j, itrc, nvd3, nvd4
      integer :: recdim, status, varid
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
      integer :: frecdim
# endif
# if defined IS4DVAR && defined LANCZOS
      integer :: MinnerDim, NinnerDim, NouterDim
      integer :: vardim(2)
# endif
      integer :: DimIDs(31), t2dgrd(3), u2dgrd(3), v2dgrd(3)
      integer :: Vsize(4)

      integer :: def_dim, def_var

# ifdef SOLVE3D
      integer :: t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
#  ifdef ADJUST_STFLUX
      integer :: t3dfrc(4)
#  endif
# endif
# ifdef ADJUST_WSTRESS
      integer :: u3dfrc(4), v3dfrc(4)
# endif

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=ADJname(ng)
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
!  Create a new adjoint history file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        IF (OutThread) THEN
          status=nf90_create(TRIM(ncname), nf90_clobber, ncADJid(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag, 1)
# endif
        IF (exit_flag.ne.NoError) RETURN
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, ncADJid(ng), 1)
# endif
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xi_rho',         &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xi_u',           &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xi_v',           &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xi_psi',         &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'eta_rho',        &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'eta_u',          &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'eta_v',          &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'eta_psi',        &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

# if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xy_rho',         &
     &                 IOBOUNDS(ng)%xy_rho, DimIDs(17))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xy_u',           &
     &                 IOBOUNDS(ng)%xy_u, DimIDs(18))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xy_v',           &
     &                 IOBOUNDS(ng)%xy_v, DimIDs(19))
        IF (exit_flag.ne.NoError) RETURN
# endif

# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xyz_rho',        &
     &                 IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xyz_u',          &
     &                 IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xyz_v',          &
     &                 IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xyz_w',          &
     &                 IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
        IF (exit_flag.ne.NoError) RETURN
#  endif

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'N',              &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 's_rho',          &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 's_w',            &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'tracer',         &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT_NOT_YET
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nbed',           &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN

#   if defined WRITE_WATER && defined MASKING
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'xybed',          &
     &                 IOBOUNDS(ng)%xy_rho*Nbed, DimIDs(24))
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nphy',           &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nbac',           &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Ndom',           &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nfec',           &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'boundary',       &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nstate',         &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nadjust',        &
     &                 Nfrec(ng), DimIDs(30))
        IF (exit_flag.ne.NoError) RETURN
# endif

# if defined IS4DVAR && defined LANCZOS
        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Ninner',         &
     &                 Ninner, NinnerDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Minner',         &
     &                 Ninner+1, MinnerDim)
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iADM, ncADJid(ng), ncname, 'Nouter',         &
     &                 Nouter, NouterDim)
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iADM, ncADJid(ng), ncname,                   &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
# if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
        frecdim=DimIDs(30)
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
#  ifdef ADJUST_STFLUX
        t3dfrc(1)=DimIDs( 1)
        t3dfrc(2)=DimIDs( 5)
        t3dfrc(3)=frecdim
        t3dfrc(4)=DimIDs(12)
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
#  ifdef ADJUST_WSTRESS
        u3dfrc(1)=DimIDs( 2)
        u3dfrc(2)=DimIDs( 6)
        u3dfrc(3)=frecdim
        u3dfrc(4)=DimIDs(12)
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
#  ifdef ADJUST_WSTRESS
        v3dfrc(1)=DimIDs( 3)
        v3dfrc(2)=DimIDs( 7)
        v3dfrc(3)=frecdim
        v3dfrc(4)=DimIDs(12)
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
# endif
!
!  Initialize unlimited time record dimension.
!
        tADJindx(ng)=0
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
        CALL def_info (ng, iADM, ncADJid(ng), ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define time-varying variables.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iADM, ncADJid(ng), adjVid(idtime,ng),        &
     &                 NF_TYPE, 1, recdim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef PROPAGATOR
!
!  Define Ritz eigenvalues and Ritz eigenvectors Euclidean norm.
!
        Vinfo( 1)='Ritz_rvalue'
        Vinfo( 2)='Ritz eigenvalues'
        status=def_var(ng, iADM, ncADJid(ng), varid,                    &
     &                 NF_TYPE, 1, recdim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

#  if defined AFT_EIGENMODES
        Vinfo( 1)='Ritz_ivalue'
        Vinfo( 2)='imaginary Ritz eigenvalues'
        status=def_var(ng, iADM, ncADJid(ng), varid,                    &
     &                 NF_TYPE, 1, recdim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif

        Vinfo( 1)='Ritz_norm'
        Vinfo( 2)='Ritz eigenvectors Euclidean norm'
        status=def_var(ng, iADM, ncADJid(ng), varid,                    &
     &                 NF_TYPE, 1, recdim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined IS4DVAR && defined LANCZOS
!
!  Define Lanczos algorithm coefficients which can be used to
!  compute the sensitivity of the observations to the 4DVAR
!  data assimilation system.
!
        Vinfo( 1)='cg_beta'
        Vinfo( 2)='conjugate gradient beta coefficient'
        vardim(1)=MinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iADM, ncADJid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='cg_delta'
        Vinfo( 2)='Lanczos algorithm delta coefficient'
        vardim(1)=NinnerDim
        vardim(2)=NouterDim
        status=def_var(ng, iADM, ncADJid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

        Vinfo( 1)='cg_zv'
        Vinfo( 2)='Lanczos recurrence eigenvectors'
        vardim(1)=NinnerDim
        vardim(2)=NinnerDim
        status=def_var(ng, iADM, ncADJid(ng), varid, NF_FRST,           &
     &                 2, vardim, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.  Notice that the stress has its
!  own fixed time-dimension (of size Nfrec) to allow 4DVAR adjustments
!  at other times in addition to initialization time.
!
        Vinfo( 1)=Vname(1,idUsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
        Vinfo( 3)=Vname(3,idUsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
        status=def_var(ng, iADM, ncADJid(ng), adjVid(idUsms,ng),        &
     &                 NF_FOUT, nvd4, u3dfrc, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface V-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
        Vinfo( 3)=Vname(3,idVsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
        status=def_var(ng, iADM, ncADJid(ng), adjVid(idVsms,ng),        &
     &                 NF_FOUT, nvd4, v3dfrc, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# if defined ADJUST_STFLUX && defined SOLVE3D
!
!  Define surface net heat flux. Notice that different tracer fluxes
!  are written at their own fixed time-dimension (of size Nfrec) to
!  allow 4DVAR adjustments at other times in addition to initial time.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTsur(itrc))
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
          Vinfo( 3)=Vname(3,idTsur(itrc))
          IF (itrc.eq.itemp) THEN
            Vinfo(11)='upward flux, cooling'
            Vinfo(12)='downward flux, heating'
          ELSE IF (itrc.eq.isalt) THEN
            Vinfo(11)='upward flux, freshening (net precipitation)'
            Vinfo(12)='downward flux, salting (net evaporation)'
          END IF
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
          status=def_var(ng, iADM, ncADJid(ng),                         &
     &                   adjVid(idTsur(itrc),ng), NF_FOUT,              &
     &                   nvd4, t3dfrc, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
!
!  Define bathymetry.
!
        Vinfo( 1)=Vname(1,idbath)
        WRITE (Vinfo( 2),40) TRIM(Vname(2,idbath))
        Vinfo( 3)='meter-1'
        Vinfo(14)=Vname(4,idbath)
        Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idbath,ng),r8)
        status=def_var(ng, iADM, ncADJid(ng), adjVid(idbath,ng),        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define free-surface.
!
        IF (Hout(idFsur,ng)) THEN
          Vinfo( 1)=Vname(1,idFsur)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idFsur))
          Vinfo( 3)='meter-1'
          Vinfo(14)=Vname(4,idFsur)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idFsur,ng),      &
     &                   NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the XI-direction.
!
        IF (Hout(idUbar,ng)) THEN
          Vinfo( 1)=Vname(1,idUbar)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbar))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idUbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idUbar,ng),      &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 2D momentum in the ETA-direction.
!
        IF (Hout(idVbar,ng)) THEN
          Vinfo( 1)=Vname(1,idVbar)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbar))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idVbar)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idVbar,ng),      &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        IF (Hout(idUvel,ng)) THEN
          Vinfo( 1)=Vname(1,idUvel)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUvel))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idUvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idUvel,ng),      &
     &                   NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define 3D momentum component in the ETA-direction.
!
        IF (Hout(idVvel,ng)) THEN
          Vinfo( 1)=Vname(1,idVvel)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvel))
          Vinfo( 3)='second meter-1'
          Vinfo(14)=Vname(4,idVvel)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idVvel,ng),      &
     &                   NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTvar(itrc)))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
#  ifdef SEDIMENT_NOT_YET
            DO i=1,NST
              IF (itrc.eq.idsed(i)) THEN
                WRITE (Vinfo(19),50) 1000.0_r8*Sd50(i,ng)
              END IF
            END DO
#  endif
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(r3dvar,r8)
            status=def_var(ng, iADM, ncADJid(ng), adjTid(itrc,ng),      &
     &                     NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
!
!  Define density anomaly.
!
        IF (Hout(idDano,ng)) THEN
          Vinfo( 1)=Vname(1,idDano)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idDano))
          Vinfo( 3)=Vname(3,idDano)
          Vinfo(14)=Vname(4,idDano)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idDano,ng),      &
     &                   NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical viscosity coefficient.
!
        IF (Hout(idVvis,ng)) THEN
          Vinfo( 1)=Vname(1,idVvis)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVvis))
          Vinfo( 3)=Vname(3,idVvis)
          Vinfo(14)=Vname(4,idVvis)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idVvis,ng),      &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define vertical diffusion coefficient for potential temperature.
!
        IF (Hout(idTdif,ng)) THEN
          Vinfo( 1)=Vname(1,idTdif)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idTdif))
          Vinfo( 3)=Vname(3,idTdif)
          Vinfo(14)=Vname(4,idTdif)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idTdif,ng),      &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  ifdef SALINITY
!
!  Define vertical diffusion coefficient for salinity.
!
        IF (Hout(idSdif,ng)) THEN
          Vinfo( 1)=Vname(1,idSdif)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idSdif))
          Vinfo( 3)=Vname(3,idSdif)
          Vinfo(14)=Vname(4,idSdif)
          Vinfo(16)=Vname(1,idtime)
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idSdif,ng),      &
     &                   NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
#  endif
#  ifndef ADJUST_STFLUX
!
!  Define surface tracer fluxes.
!
        DO itrc=1,NT(ng)
          IF (Hout(idTsur(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTsur(itrc))
            WRITE (Vinfo( 2),40) TRIM(Vname(2,idTsur(itrc)))
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
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idTsur(itrc),ng),r8)
            status=def_var(ng, iADM, ncADJid(ng),                       &
     &                     adjVid(idTsur(itrc),ng), NF_FOUT,            &
     &                     nvd3, t2dgrd, Aval, Vinfo, ncname)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END DO
#  endif
# endif
# ifndef ADJUST_WSTRESS
!
!  Define surface U-momentum stress.
!
        IF (Hout(idUsms,ng)) THEN
          Vinfo( 1)=Vname(1,idUsms)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUsms))
          Vinfo( 3)=Vname(3,idUsms)
          Vinfo(14)=Vname(4,idUsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idUsms,ng),      &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define surface V-momentum stress.
!
        IF (Hout(idVsms,ng)) THEN
          Vinfo( 1)=Vname(1,idVsms)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVsms))
          Vinfo( 3)=Vname(3,idVsms)
          Vinfo(14)=Vname(4,idVsms)
          Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
#  endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idVsms,ng),      &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
# endif
!
!  Define bottom U-momentum stress.
!
        IF (Hout(idUbms,ng)) THEN
          Vinfo( 1)=Vname(1,idUbms)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idUbms))
          Vinfo( 3)=Vname(3,idUbms)
          Vinfo(14)=Vname(4,idUbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_u'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idUbms,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idUbms,ng),      &
     &                   NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!  Define bottom V-momentum stress.
!
        IF (Hout(idVbms,ng)) THEN
          Vinfo( 1)=Vname(1,idVbms)
          WRITE (Vinfo( 2),40) TRIM(Vname(2,idVbms))
          Vinfo( 3)=Vname(3,idVbms)
          Vinfo(14)=Vname(4,idVbms)
          Vinfo(16)=Vname(1,idtime)
# if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_v'
# endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(Iinfo(1,idVbms,ng),r8)
          status=def_var(ng, iADM, ncADJid(ng), adjVid(idVbms,ng),      &
     &                   NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        IF (OutThread) THEN
          status=nf90_enddef(ncADJid(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,60) TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag, 1)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iADM, ncADJid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing adjoint file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=ADJname(ng)
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iADM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iADM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open adjoint file for read/write.
!
        IF (OutThread) THEN
          status=nf90_open(TRIM(ncname), nf90_write, ncADJid(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,70) TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, exit_flag, 1)
# endif
     IF (exit_flag.ne.NoError) RETURN
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iADM, ncADJid(ng), 1)
# endif
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  adjoint variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            adjVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFsur))) THEN
            got_var(idFsur)=.TRUE.
            adjVid(idFsur,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbar))) THEN
            got_var(idUbar)=.TRUE.
            adjVid(idUbar,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbar))) THEN
            got_var(idVbar)=.TRUE.
            adjVid(idVbar,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
            got_var(idUvel)=.TRUE.
            adjVid(idUvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
            got_var(idVvel)=.TRUE.
            adjVid(idVvel,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            adjVid(idDano,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
            got_var(idVvis)=.TRUE.
            adjVid(idVvis,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
            got_var(idTdif)=.TRUE.
            adjVid(idTdif,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
            got_var(idSdif)=.TRUE.
            adjVid(idSdif,ng)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUsms))) THEN
            got_var(idUsms)=.TRUE.
            adjVid(idUsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVsms))) THEN
            got_var(idVsms)=.TRUE.
            adjVid(idVsms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUbms))) THEN
            got_var(idUbms)=.TRUE.
            adjVid(idUbms,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVbms))) THEN
            got_var(idVbms)=.TRUE.
            adjVid(idVbms,ng)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              adjTid(itrc,ng)=var_id(i)
            ELSE IF (TRIM(var_name(i)).eq.                              &
     &               TRIM(Vname(1,idTsur(itrc)))) THEN
              got_var(idTsur(itrc))=.TRUE.
              adjVid(idTsur(itemp),ng)=var_id(i)
            END IF
          END DO
#endif
        END DO
!
!  Check if adjoint variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur).and.Hout(idFsur,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar).and.Hout(idUbar,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar).and.Hout(idVbar,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idUvel).and.Hout(idUvel,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel).and.Hout(idVvel,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano).and.Hout(idDano,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idUsms).and.Hout(idUsms,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVsms).and.Hout(idVsms,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVsms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbms).and.Hout(idUbms,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idUbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbms).and.Hout(idVbms,ng)) THEN
          IF (Master) WRITE (stdout,80) TRIM(Vname(1,idVbms)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Hout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,80) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
          IF (.not.got_var(idTsur(itrc)).and.Hout(idTsur(itrc),ng)) THEN
            IF (Master) WRITE (stdout,80) TRIM(Vname(1,idTsur(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
        IF (ndefADJ(ng).gt.0) THEN
          tADJindx(ng)=((ntstart(ng)-1)-                                &
     &                  ndefADJ(ng)*((ntstart(ng)-1)/ndefADJ(ng)))/     &
     &                 nADJ(ng)
        ELSE
          tADJindx(ng)=(ntstart(ng)-1)/nADJ(ng)
        END IF
        tADJindx(ng)=MIN(tADJindx(ng),rec_size)
      END IF QUERY
!
  10  FORMAT (3x,'AD_DEF_HIS   - creating adjoint  file: ',a)
  20  FORMAT (3x,'AD_DEF_HIS   - inquiring adjoint file: ',a)
  30  FORMAT (/,' AD_DEF_HIS - unable to create adjoint NetCDF file: ', &
     &        a)
  40  FORMAT ('adjoint',1x,a)
  50  FORMAT (1pe11.4,1x,'millimeter')
  60  FORMAT (/,' AD_DEF_HIS - unable to leave definition mode for',    &
     &        ' file: ', a)
  70  FORMAT (/,' AD_DEF_HIS - unable to open adjoint NetCDF file: ',a)
  80  FORMAT (/,' AD_DEF_HIS - unable to find variable: ',a,2x,         &
     &        ' in adjoint NetCDF file: ',a)

      RETURN
      END SUBROUTINE ad_def_his
#else
      SUBROUTINE ad_def_his
      RETURN
      END SUBROUTINE ad_def_his
#endif
