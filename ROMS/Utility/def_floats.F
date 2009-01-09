#include "cppdefs.h"
#ifdef FLOATS
      SUBROUTINE def_floats (ng,ldef)
!
!svn $Id: def_floats.F 733 2008-09-07 01:56:45Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates FLOATS NetCDF file, it defines dimensions,     !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_collect, mp_bcasti
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

      logical :: got_var(-6:NV)

      integer :: fltdim, i, itrc, j, l, nrec, nvd
      integer :: recdim, status

      integer :: DimIDs(29), fgrd(2), start(2), total(2)
      integer :: Vsize(4)

      integer :: def_info, def_var, wrt_info

# ifdef DISTRIBUTE
      integer :: Npts

      real(r8), parameter :: Fspv = 0.0_r8
      real(r8), dimension(Nfloats(ng)*NFV(ng)) :: Fwrk
# endif

      real(r8) :: Aval(6), Tinp(Nfloats(ng))

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: dimnam, fname, ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=FLTname(ng)
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
!  Create a new floats data file.
!=======================================================================
!
!  Create floats NetCDF file.
!
      IF (exit_flag.ne.NoError) RETURN
      IF (ldef.and.OutThread) THEN
        status=nf90_create(TRIM(ncname), nf90_clobber, ncFLTid(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef DISTRIBUTE
      IF (ldef) THEN
        CALL mp_bcasti (ng, iNLM, ncFLTid(ng), 1)
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
      IF (ldef.and.OutThread) THEN
# ifdef SOLVE3D
        status=nf90_def_dim(ncFLTid(ng),'s_rho'   ,N(ng),               &
     &                      DimIDs( 9))
        status=nf90_def_dim(ncFLTid(ng),'s_w'     ,N(ng)+1,             &
     &                      DimIDs(10))
        status=nf90_def_dim(ncFLTid(ng),'tracer'  ,NT(ng),              &
     &                      DimIDs(11))
#  ifdef SEDIMENT
        status=nf90_def_dim(ncFLTid(ng),'Nbed'    ,Nbed,                &
     &                      DimIDs(16))
#  endif
#  ifdef ECOSIM
        status=nf90_def_dim(ncFLTid(ng),'Nphy'    ,Nphy,                &
     &                      DimIDs(25))
        status=nf90_def_dim(ncFLTid(ng),'Nbac'    ,Nbac,                &
     &                      DimIDs(26))
        status=nf90_def_dim(ncFLTid(ng),'Ndom'    ,Ndom,                &
     &                      DimIDs(27))
        status=nf90_def_dim(ncFLTid(ng),'Nfec'    ,Nfec,                &
     &                      DimIDs(28))
#  endif
# endif
        status=nf90_def_dim(ncFLTid(ng),'drifter' ,Nfloats(ng),         &
     &                      DimIDs(15))
        status=nf90_def_dim(ncFLTid(ng),'boundary',4,                   &
     &                      DimIDs(14))
# ifdef FOUR_DVAR
        status=nf90_def_dim(ncFLTid(ng),'Nstate  ',NstateVar(ng),       &
     &                      DimIDs(29))
# endif
        status=nf90_def_dim(ncFLTid(ng),TRIM(ADJUSTL(Vname(5,idtime))), &
     &                      nf90_unlimited,DimIDs(12))
        recdim=DimIDs(12)
        fltdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        tFLTindx(ng)=0
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
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        status=def_info(ng,ncFLTid(ng),ncname,DimIDs)
        if (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
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
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ncFLTid(ng),fltVid(idtime,ng),NF_TYPE,           &
     &                 1,recdim,Aval,Vinfo,ncname)
!
!  Define floats X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Lm(ng)+1,r8)
        Vinfo(14)='Xgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncFLTid(ng),fltVid(idXgrd,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
!
!  Define floats Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Mm(ng)+1,r8)
        Vinfo(14)='Ygrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncFLTid(ng),fltVid(idYgrd,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
# ifdef SOLVE3D
!
!  Define floats Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid floats locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(N(ng),r8)
        Vinfo(14)='Zgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncFLTid(ng),fltVid(idZgrd,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
# endif
!
!  Define floats (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of floats trajectories'
          Vinfo( 3)='degree_east'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lon, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(2)=-180.0_r8
          Aval(3)=180.0_r8
          Aval(4)=spval
          status=def_var(ncFLTid(ng),fltVid(idglon,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of floats trajectories'
          Vinfo( 3)='degree_north'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lat, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(2)=-90.0_r8
          Aval(3)=90.0_r8
          Aval(4)=spval
          status=def_var(ncFLTid(ng),fltVid(idglat,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(4)=spval
          status=def_var(ncFLTid(ng),fltVid(idglon,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
          Vinfo( 1)='y'
          Vinfo( 2)='y-location of floats trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(4)=spval
          status=def_var(ncFLTid(ng),fltVid(idglat,ng),NF_FOUT,         &
     &                   2,fgrd,Aval,Vinfo,ncname)
        END IF
# ifdef SOLVE3D
!
!  Define floats depths.
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of floats trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncFLTid(ng),fltVid(iddpth,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
        Vinfo(17)='missing_value'
        Aval(4)=spval
        status=def_var(ncFLTid(ng),fltVid(idDano,ng),NF_FOUT,           &
     &                 2,fgrd,Aval,Vinfo,ncname)
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
          Vinfo(17)='missing_value'
          Aval(4)=spval
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),20) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
          status=def_var(ncFLTid(ng),fltTid(itrc,ng),NF_FOUT,           &
     &                   2,fgrd,Aval,Vinfo,ncname)
        END DO
# endif
!
!  Initialize unlimited time record dimension.
!
        tFLTindx(ng)=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        status=nf90_enddef(ncFLTid(ng))
      END IF
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
      IF (ldef) THEN
        status=wrt_info(ng, iNLM, ncFLTid(ng), Master, ncname)
        if (exit_flag.ne.NoError) RETURN
      END IF
!
!=======================================================================
!  Open an existing floats file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      IF (.not.ldef.and.OutThread) THEN
!
!  Inquire about the contents of floats NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
        ncname=FLTname(ng)
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, nrec, nvd, Vsize)
        if (exit_flag.ne.NoError) RETURN
!
!  Open floats file for read/write.
!
        status=nf90_open(TRIM(ncname),nf90_write,ncFLTid(ng))
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
!  Inquire size of the "drifter" dimension.
!
        status=nf90_inq_dimid(ncFLTid(ng),'drifter',fltdim)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,60) 'drifter', TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf90_inquire_dimension(ncFLTid(ng),fltdim,dimnam,        &
     &                                Nfloats(ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,70) 'drifter', TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
!
!  Scan variable list from input NetCDF and activate switches for
!  float variables. Get variable IDs.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            status=nf90_inq_varid(ncFLTid(ng),TRIM(Vname(1,idtime)),    &
     &                            fltVid(idtime,ng))
          END IF
          IF (TRIM(varnam(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            status=nf90_inq_varid(ncFLTid(ng),'Xgrid',fltVid(idXgrd,ng))
          END IF
          IF (TRIM(varnam(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            status=nf90_inq_varid(ncFLTid(ng),'Ygrid',fltVid(idYgrd,ng))
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            status=nf90_inq_varid(ncFLTid(ng),'Zgrid',fltVid(idZgrd,ng))
          END IF
# endif
          IF (spherical) THEN
            IF (TRIM(varnam(i)).eq.'lon') THEN
              got_var(idglon)=.TRUE.
              status=nf90_inq_varid(ncFLTid(ng),'lon',fltVid(idglon,ng))
            END IF
            IF (TRIM(varnam(i)).eq.'lat') THEN
              got_var(idglat)=.TRUE.
              status=nf90_inq_varid(ncFLTid(ng),'lat',fltVid(idglat,ng))
            END IF
          ELSE
            IF (TRIM(varnam(i)).eq.'x') THEN
              got_var(idglon)=.TRUE.
              status=nf90_inq_varid(ncFLTid(ng),'x',fltVid(idglon,ng))
            END IF
            IF (TRIM(varnam(i)).eq.'y') THEN
              got_var(idglat)=.TRUE.
              status=nf90_inq_varid(ncFLTid(ng),'y',fltVid(idglat,ng))
            END IF
          END IF
# ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            status=nf90_inq_varid(ncFLTid(ng),'depth',fltVid(iddpth,ng))
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            status=nf90_inq_varid(ncFLTid(ng),TRIM(Vname(1,idDano)),    &
     &                            fltVid(idDano,ng))
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
             got_var(idTvar(itrc))=.TRUE.
             status=nf90_inq_varid(ncFLTid(ng),                         &
     &                             TRIM(Vname(1,idTvar(itrc))),         &
     &                             fltTid(itrc,ng))
            END IF
          END DO
# endif
        END DO
!
!  Check if floats variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          WRITE (stdout,80) TRIM(Vname(1,idtime)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idXgrd)) THEN
          WRITE (stdout,80) 'Xgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idYgrd)) THEN
          WRITE (stdout,80) 'Ygrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idZgrd)) THEN
          WRITE (stdout,80) 'Zgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idglon)) THEN
          IF (spherical) THEN
            WRITE (stdout,80) 'lon', TRIM(ncname)
          ELSE
            WRITE (stdout,80) 'x', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idglat)) THEN
          IF (spherical) THEN
            WRITE (stdout,80) 'lat', TRIM(ncname)
          ELSE
            WRITE (stdout,80) 'y', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(iddpth)) THEN
          WRITE (stdout,80) 'depth', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          WRITE (stdout,80) TRIM(Vname(1,idDano)), TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            WRITE (stdout,80) TRIM(Vname(1,idTvar(itrc))), TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize floats positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up floats time record.
!
        IF (frrec(ng).lt.0) THEN
          tFLTindx(ng)=tsize
        ELSE
         tFLTindx(ng)=ABS(frrec(ng))
!!       tFLTindx(ng)=1+(ntstart(ng)-1)/nFLT(ng)
        END IF
        NrecFLT=tFLTindx(ng)
!
!  Read in floats nondimentional horizontal positions.
!
        start(1)=1
        total(1)=Nfloats(ng)
        start(2)=tFLTindx(ng)
        total(2)=1
        status=nf90_get_var(ncFLTid(ng),fltVid(idXgrd,ng),              &
     &                      Tinp,start,total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,90) 'Xgrid', tFLTindx(ng), TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            FLT(ng)%bounded(l)=.FALSE.
          ELSE
            FLT(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              FLT(ng)%track(ixgrd,i,l)=Tinp(l)
              FLT(ng)%track(ixrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
        status=nf90_get_var(ncFLTid(ng),fltVid(idYgrd,ng),              &
     &                      Tinp,start,total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,90) 'Ygrid', tFLTindx(ng), TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            FLT(ng)%bounded(l)=.FALSE.
          ELSE
            FLT(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              FLT(ng)%track(iygrd,i,l)=Tinp(l)
              FLT(ng)%track(iyrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
# ifdef SOLVE3D
        status=nf90_get_var(ncFLTid(ng),fltVid(idZgrd,ng),              &
     &                      Tinp,start,total)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,90) 'Zgrid', tFLTindx(ng), TRIM(ncname)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        DO l=1,Nfloats(ng)
          IF ((Tinp(l).gt.REAL(N(ng),r8)).or.(Tinp(l).lt.0.0_r8)) THEN
            FLT(ng)%bounded(l)=.FALSE.
          ELSE
            FLT(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              FLT(ng)%track(izgrd,i,l)=Tinp(l)
              FLT(ng)%track(izrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
# endif
      END IF

# ifdef DISTRIBUTE
!
!  Broadcast restart float information to all nodes.
!
      IF (.not.ldef) THEN

        Npts=NFV(ng)*(NFT+1)*Nfloats(ng)

        CALL mp_collect (ng, iNLM, Npts, Fspv, FLT(ng)%track)
!
!  Collect the bounded status switch.
!
        Fwrk=Fspv
        DO l=1,Nfloats(ng)
          IF (FLT(ng)%bounded(l)) THEN
            Fwrk(l)=1.0_r8
          END IF
        END DO
        CALL mp_collect (ng, iNLM, Nfloats(ng), Fspv, Fwrk)
        DO l=1,Nfloats(ng)
          IF (Fwrk(l).ne.Fspv) THEN
            FLT(ng)%bounded(l)=.TRUE.
          ELSE
            FLT(ng)%bounded(l)=.FALSE.
          END IF
        END DO
      END IF
# endif
!
  10  FORMAT (6x,'DEF_FLOATS  - creating floats file: ',a)
  20  FORMAT (6x,'DEF_FLOATS  - inquiring history file: ',a)
  30  FORMAT (/,' DEF_FLOATS - unable to create floats NetCDF',         &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_FLOATS - unable to open floats NetCDF file: ',a)
  60  FORMAT (/,' DEF_FLOATS - error while inquiring dimension ID ',    &
     &        ' for: ',a,2x,/,13x,'in floats NetCDF file: ',a)
  70  FORMAT (/,' DEF_FLOATS - error while inquiring size of',          &
     &        ' dimension: ',a,2x,/,13x,'in floats NetCDF file: ',a)
  80  FORMAT (/,' DEF_FLOATS - unable to find variable: ',a,2x,         &
     &        ' in floats NetCDF file: ',a)
  90  FORMAT (/,' DEF_FLOATS - error while reading variable: ',a,2x,    &
     &        ' at time record = ', i6.6,/,13x,                         &
     &        ' in floats NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_floats
#else
      SUBROUTINE def_floats
      RETURN
      END SUBROUTINE def_floats
#endif
