#include "cppdefs.h"
      SUBROUTINE read_PhyPar (model, inp, out, Lwrite)
!
!svn $Id: read_phypar.F 793 2016-05-10 00:08:25Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports physical model input parameters.     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
#ifdef MODEL_COUPLING
      USE mod_coupler
#endif
#ifdef COAWST_COUPLING
      USE mct_coupler_params
#endif
#ifdef INWAVE_MODEL
      USE mod_inwave_params
      USE inwave_iounits
#endif
#if defined FOUR_DVAR || defined VERIFICATION
      USE mod_fourdvar
#endif
      USE mod_iounits
      USE mod_ncparam
#ifdef NESTING
      USE mod_nesting
#endif
      USE mod_netcdf
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
#if defined VEGETATION
      USE mod_vegetation
#endif
      USE mod_stepping
#ifdef PROPAGATOR
      USE mod_storage
#endif
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: got_Ngrids, got_NestLayers
      logical :: find_file, obc_data
      logical :: Lvalue(1)
#if defined SOLVE3D && defined SEDIMENT
      logical :: LreadNCS = .FALSE.
      logical :: LreadNNS = .FALSE.
#endif
#if defined SOLVE3D && defined T_PASSIVE
      logical, allocatable :: Linert(:,:)
#endif
#if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
      logical, allocatable :: Lbottom(:,:)
#endif
      logical, allocatable :: Ltracer(:,:)
#if defined AD_SENSITIVITY    || defined IS4DVAR_SENSITIVITY || \
    defined FORCING_SV        || defined OPT_OBSERVATIONS    || \
    defined SENSITIVITY_4DVAR || defined SO_SEMI             || \
    defined STOCHASTIC_OPT
      logical, allocatable :: Ladsen(:)
#endif
#if defined DIAGNOSTICS_UV
      logical, allocatable :: Ldiag(:)
#endif

      integer :: Lstr, Npts, Nval, i, itrc, ivar, k, ng, nl, status
      integer :: ifield, ifile, igrid, itracer, nline, max_Ffiles
      integer :: ibcfile, iclmfile
      integer :: Cdim, Clen, Rdim
      integer :: Ivalue(1)

      integer :: decode_line, load_i, load_l, load_lbc, load_r
      integer :: load_s1d, load_s2d

      integer, allocatable :: Nfiles(:)
      integer, allocatable :: Ncount(:,:)
      integer, allocatable :: NBCcount(:,:)
      integer, allocatable :: NCLMcount(:,:)

      real(r8), allocatable :: Rtracer(:,:)
      real(r8), allocatable :: tracer(:,:)

      real(r8) :: Rvalue(1)

      real(r8), dimension(200) :: Rval

      character (len=1  ), parameter :: blank = ' '
      character (len=19 ) :: ref_att
      character (len=40 ) :: KeyWord, text
      character (len=50 ) :: label
      character (len=256) :: fname, line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      ifile=1                            ! multiple file counter
      ibcfile=1                          ! multiple BC file counter
      iclmfile=1                         ! multiple CLM file counter
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      nline=0                            ! LBC multi-line counter
      DO i=1,LEN(label)
        label(i:i)=blank
      END DO
      got_Ngrids=.FALSE.
      got_NestLayers=.FALSE.
      Cdim=SIZE(Cval,1)
      Clen=LEN(Cval(1))
      Rdim=SIZE(Rval,1)
!
!-----------------------------------------------------------------------
!  Read in physical model parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('TITLE')
              IF (Nval.eq.1) THEN
                title=TRIM(ADJUSTL(Cval(Nval)))
              ELSE
                WRITE(title,'(a,1x,a)') TRIM(ADJUSTL(title)),           &
     &                                  TRIM(ADJUSTL(Cval(Nval)))
              END IF
            CASE ('MyAppCPP')
              DO i=1,LEN(MyAppCPP)
                MyAppCPP(i:i)=blank
              END DO
              MyAppCPP=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('VARNAME')
              DO i=1,LEN(varname)
                varname(i:i)=blank
              END DO
              varname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('Ngrids')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Ngrids=Ivalue(1)
              IF (Ngrids.le.0) THEN
                IF (Master) WRITE (out,290) 'Ngrids', Ngrids,           &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              got_Ngrids=.TRUE.
              CALL allocate_param       ! Start allocating variables in
              CALL allocate_parallel    ! modules that solely depend on
              CALL allocate_iounits     ! the number of nested grids
              CALL allocate_stepping
#if defined INWAVE_MODEL
              CALL allocate_inwave_params
              CALL allocate_inwave_iounits
#endif
#if defined FOUR_DVAR || defined VERIFICATION
              CALL allocate_fourdvar
#endif
#if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
              IF (.not.allocated(Lbottom)) THEN
                allocate ( Lbottom(MBOTP,Ngrids) )
              END IF
#endif
#if defined DIAGNOSTICS_UV
              IF (.not.allocated(Ldiag)) THEN
                allocate ( Ldiag(Ngrids) )
              END IF
#endif
              IF (.not.allocated(Nfiles)) THEN
                allocate ( Nfiles(Ngrids) )
                Nfiles(1:Ngrids)=0
              END IF
            CASE ('NestLayers')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NestLayers=Ivalue(1)
              IF (NestLayers.lt.1) THEN
                IF (Master) WRITE (out,290) 'NestLayers', NestLayers,   &
     &            'must be greater or equal than one.'
                exit_flag=5
                RETURN
              END IF
#ifndef NESTING
              IF (NestLayers.gt.1) THEN
                IF (Master) WRITE (out,290) 'NestLayers', NestLayers,   &
     &            'must be equal to one in non-nesting applications.'
                exit_flag=5
                RETURN
              END IF
#endif
              got_NestLayers=.TRUE.
              IF (.not.allocated(GridsInLayer)) THEN
                allocate ( GridsInLayer(NestLayers) )
#ifndef NESTING
                GridsInLayer(1:NestLayers)=1
#endif
              END IF
              IF (.not.allocated(GridNumber)) THEN
                allocate ( GridNumber(Ngrids,NestLayers) )
                GridNumber(1:Ngrids,1:NestLayers)=0        ! Important
              END IF
            CASE ('GridsInLayer')
              IF (.not.got_NestLayers) THEN
                IF (Master) WRITE (out,320) 'NestLayers',               &
     &            'Add "NestLayers" keyword before GridsInLayer.'
                exit_flag=5
                RETURN
              END IF
              Npts=load_i(Nval, Rval, NestLayers, GridsInLayer)
              ng=0
              DO nl=1,NestLayers
                DO i=1,GridsInLayer(nl)
                  ng=ng+1                  ! order of grids are very in
                  GridNumber(i,nl)=ng      ! nesting applications. See
                END DO                     ! WikiROMS for details.
              END DO
            CASE ('Lm')
              IF (.not.got_Ngrids) THEN
                IF (Master) WRITE (out,320) 'Ngrids',                   &
     &            'Add "Ngrids" keyword before grid dimension (Lm, Mm).'
                exit_flag=5
                RETURN
              END IF
              Npts=load_i(Nval, Rval, Ngrids, Lm)
              DO ng=1,Ngrids
                IF (Lm(ng).le.0) THEN
                  IF (Master) WRITE (out,300) 'Lm', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('Mm')
              Npts=load_i(Nval, Rval, Ngrids, Mm)
              DO ng=1,Ngrids
                IF (Mm(ng).le.0) THEN
                  IF (Master) WRITE (out,300) 'Mm', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('N')
              Npts=load_i(Nval, Rval, Ngrids, N)
              DO ng=1,Ngrids
                IF (N(ng).lt.0) THEN
                  IF (Master) WRITE (out,300) 'N', ng,                  &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
#ifdef INWAVE_MODEL
            CASE ('ND')
              Npts=load_i(Nval, Rval, Ngrids, ND)
              DO ng=1,Ngrids
                IF (ND(ng).lt.0) THEN
                  IF (Master) WRITE (out,300) 'ND', ng,                 &
     &              'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
#endif
#if defined ESTUARYBGC && defined SPECTRAL_LIGHT
            CASE ('NBAND') 
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NBAND=Ivalue(1)
              IF (NBAND.le.0) THEN
                IF (Master) WRITE (out,290) 'NBAND = ', NBAND,          &
     &            'must be greater than zero.'
                exit_flag=5
              RETURN
            END IF
#endif
#if defined SEDIMENT && defined SOLVE3D
            CASE ('Nbed')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nbed=Ivalue(1)
              IF (Nbed.le.0) THEN
                IF (Master) WRITE (out,290) 'Nbed = ', Nbed,            &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
#endif
#ifdef SOLVE3D
            CASE ('NAT')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NAT=Ivalue(1)
              IF ((NAT.lt.1).or.(NAT.gt.2)) THEN
                IF (Master) WRITE (out,290) 'NAT = ', NAT,              &
     &            'make sure that NAT is either 1 or 2.'
                exit_flag=5
                RETURN
              END IF
# ifdef SALINITY
              IF (NAT.ne.2) THEN
                IF (Master) WRITE (out,290) 'NAT = ', NAT,              &
     &            'make sure that NAT is equal to 2.'
                exit_flag=5
                RETURN
              END IF
# endif
#endif
#if defined T_PASSIVE && defined SOLVE3D
            CASE ('NPT')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NPT=Ivalue(1)
              IF (NPT.le.0) THEN
                IF (Master) WRITE (out,290) 'NPT = ', NPT,              &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
# ifdef AGE_MEAN
              IF (MOD(NPT,2).ne.0) THEN
                IF (Master) WRITE (out,290) 'NPT = ', NPT,              &
     &            'must be an even number when Mean Age is activated.'
                exit_flag=5
                RETURN
              END IF
# endif
#endif
#if defined SEDIMENT && defined SOLVE3D
            CASE ('NCS')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NCS=Ivalue(1)
              IF (NCS.lt.0) THEN
                IF (Master) WRITE (out,290) 'NCS = ', NCS,              &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              LreadNCS=.TRUE.
              IF (LreadNNS.and.((NCS+NNS).le.0)) THEN
                IF (Master) WRITE (out,290) 'NST = ', NCS+NNS,          &
     &            'either NCS or NNS must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              NST=NST+NCS
            CASE ('NNS')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NNS=Ivalue(1)
              IF (NNS.lt.0) THEN
                IF (Master) WRITE (out,290) 'NNS = ',                   &
     &            'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              LreadNNS=.TRUE.
              IF (LreadNCS.and.((NCS+NNS).le.0)) THEN
                IF (Master) WRITE (out,290) 'NST = ', NCS+NNS,          &
     &            'either NCS or NNS must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
              NST=NST+NNS
#endif
            CASE ('NtileI')
              Npts=load_i(Nval, Rval, Ngrids, NtileI)
              NtileX(1:Ngrids)=NtileI(1:Ngrids)
            CASE ('NtileJ')
              Npts=load_i(Nval, Rval, Ngrids, NtileJ)
              NtileE(1:Ngrids)=NtileJ(1:Ngrids)
#ifdef BIOLOGY
              CALL initialize_biology
#endif
#if defined SEDIMENT || defined BBL_MODEL
              CALL initialize_sediment
#endif
#ifdef VEGETATION
              CALL initialize_vegetation
#endif
              CALL initialize_param    ! Continue allocating/initalizing
              CALL allocate_scalars    ! variables since the application
              CALL initialize_scalars  ! number of nested grids and
              CALL allocate_ncparam    ! domain parameters are known
              CALL initialize_ncparam
#if defined AD_SENSITIVITY    || defined IS4DVAR_SENSITIVITY || \
    defined FORCING_SV        || defined OPT_OBSERVATIONS    || \
    defined SENSITIVITY_4DVAR || defined SO_SEMI             || \
    defined STOCHASTIC_OPT
              IF (.not.allocated(Ladsen)) THEN
                allocate (Ladsen(MT*Ngrids))
              END IF
#endif
              IF (.not.allocated(Ltracer)) THEN
                allocate (Ltracer(NAT+NPT,Ngrids))
              END IF
#if defined SOLVE3D && defined T_PASSIVE
              IF (.not.allocated(Linert)) THEN
                allocate (Linert(NPT,Ngrids))
              END IF
              IF (MAXVAL(inert).eq.0) THEN
                IF (Master) WRITE (out,280) 'inert'
                exit_flag=5
                RETURN
              END IF
#endif
              IF (.not.allocated(Rtracer)) THEN
                allocate (Rtracer(NAT+NPT,Ngrids))
              END IF
              IF (.not.allocated(tracer)) THEN
                allocate (tracer(MT,Ngrids))
              END IF
            CASE ('LBC(isFsur)')
              Npts=load_lbc(Nval, Cval, line, nline, isFsur, igrid,     &
     &                      0, 0, Vname(1,idFsur), LBC)
            CASE ('LBC(isUbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isUbar, igrid,     &
     &                      0, 0, Vname(1,idUbar), LBC)
            CASE ('LBC(isVbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isVbar, igrid,     &
     &                      0, 0, Vname(1,idVbar), LBC)
#ifdef WEC
            CASE ('LBC(isU2Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isU2Sd, igrid,     &
     &                      0, 0, Vname(1,idU2Sd), LBC)
            CASE ('LBC(isV2Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isV2Sd, igrid,     &
     &                      0, 0, Vname(1,idV2Sd), LBC)
#endif
#ifdef SOLVE3D
            CASE ('LBC(isUvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isUvel, igrid,     &
     &                      0, 0, Vname(1,idUvel), LBC)
            CASE ('LBC(isVvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isVvel, igrid,     &
     &                      0, 0, Vname(1,idVvel), LBC)
# ifdef WEC
            CASE ('LBC(isU3Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isU3Sd, igrid,     &
     &                      0, 0, Vname(1,idU3Sd), LBC)
            CASE ('LBC(isV3Sd)')
              Npts=load_lbc(Nval, Cval, line, nline, isV3Sd, igrid,     &
     &                      0, 0, Vname(1,idV3Sd), LBC)
# endif
# ifdef INWAVE_MODEL
            CASE ('LBC(isAC3d)')
              Npts=load_lbc(Nval, Cval, line, nline, isAC3d, igrid,     &
     &                      0, 0, Vname(1,idACen), LBC)
            CASE ('LBC(isCT3d)')
              Npts=load_lbc(Nval, Cval, line, nline, isCT3d, igrid,     &
     &                      0, 0, Vname(1,idACct), LBC)
            CASE ('LBC(isCX3d)')
              Npts=load_lbc(Nval, Cval, line, nline, isCX3d, igrid,     &
     &                      0, 0, Vname(1,idACcx), LBC)
            CASE ('LBC(isCY3d)')
              Npts=load_lbc(Nval, Cval, line, nline, isCY3d, igrid,     &
     &                      0, 0, Vname(1,idACcy), LBC)
# endif
# if defined GLS_MIXING || defined MY25_MIXING
            CASE ('LBC(isMtke)')
              Npts=load_lbc(Nval, Cval, line, nline, isMtke, igrid,     &
     &                      0, 0, Vname(1,idMtke), LBC)
# endif
# ifdef ICE_MODEL
            CASE ('LBC(isUice)')
              Npts=load_lbc(Nval, Cval, line, nline, isUice, igrid,     &
     &                      0, 0, Vname(1,idUice), LBC)
            CASE ('LBC(isVice)')
              Npts=load_lbc(Nval, Cval, line, nline, isVice, igrid,     &
     &                      0, 0, Vname(1,idVice), LBC)
            CASE ('LBC(isAice)')
              Npts=load_lbc(Nval, Cval, line, nline, isAice, igrid,     &
     &                      0, 0, Vname(1,idAice), LBC)
            CASE ('LBC(isHice)')
              Npts=load_lbc(Nval, Cval, line, nline, isHice, igrid,     &
     &                      0, 0, Vname(1,idHice), LBC)
            CASE ('LBC(isTice)')
              Npts=load_lbc(Nval, Cval, line, nline, isTice, igrid,     &
     &                      0, 0, Vname(1,idTice), LBC)
            CASE ('LBC(isHsno)')
              Npts=load_lbc(Nval, Cval, line, nline, isHsno, igrid,     &
     &                      0, 0, Vname(1,idHsno), LBC)
            CASE ('LBC(isApond)')
              Npts=load_lbc(Nval, Cval, line, nline, isApond, igrid,    &
     &                      0, 0, Vname(1,idApond), LBC)
            CASE ('LBC(isHpond)')
              Npts=load_lbc(Nval, Cval, line, nline, isHpond, igrid,    &
     &                      0, 0, Vname(1,idHpond), LBC)
            CASE ('LBC(isSig11)')
              Npts=load_lbc(Nval, Cval, line, nline, isSig11, igrid,    &
     &                      0, 0, Vname(1,idSig11), LBC)
            CASE ('LBC(isSig12)')
              Npts=load_lbc(Nval, Cval, line, nline, isSig12, igrid,    &
     &                      0, 0, Vname(1,idSig12), LBC)
            CASE ('LBC(isSig22)')
              Npts=load_lbc(Nval, Cval, line, nline, isSig22, igrid,    &
     &                      0, 0, Vname(1,idSig22), LBC)
# endif
            CASE ('LBC(isTvar)')
              IF (itracer.lt.(NAT+NPT)) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(itracer)
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      1, NAT+NPT, Vname(1,idTvar(itracer)), LBC)
#endif
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_LBC(isFsur)')
              Npts=load_lbc(Nval, Cval, line, nline, isFsur, igrid,     &
     &                      0, 0, Vname(1,idFsur), ad_LBC)
            CASE ('ad_LBC(isUbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isUbar, igrid,     &
     &                      0, 0, Vname(1,idUbar), ad_LBC)
            CASE ('ad_LBC(isVbar)')
              Npts=load_lbc(Nval, Cval, line, nline, isVbar, igrid,     &
     &                      0, 0, Vname(1,idVbar), ad_LBC)
# ifdef SOLVE3D
            CASE ('ad_LBC(isUvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isUvel, igrid,     &
     &                      0, 0, Vname(1,idUvel), ad_LBC)
            CASE ('ad_LBC(isVvel)')
              Npts=load_lbc(Nval, Cval, line, nline, isVvel, igrid,     &
     &                      0, 0, Vname(1,idVvel), ad_LBC)
#  if defined GLS_MIXING || defined MY25_MIXING
            CASE ('ad_LBC(isMtke)')
              Npts=load_lbc(Nval, Cval, line, nline, isMtke, igrid,     &
     &                      0, 0, Vname(1,idMtke), ad_LBC)
#  endif
            CASE ('ad_LBC(isTvar)')
              IF (itracer.lt.(NAT+NPT)) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(itracer)
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      1, NAT+NPT, Vname(1,idTvar(itracer)),       &
     &                      ad_LBC)
# endif
#endif
            CASE ('VolCons(west)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(iwest,:))
            CASE ('VolCons(east)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(ieast,:))
            CASE ('VolCons(south)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(isouth,:))
            CASE ('VolCons(north)')
              Npts=load_l(Nval, Cval, Ngrids, VolCons(inorth,:))
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
            CASE ('ad_VolCons(west)')
              Npts=load_l(Nval, Cval, Ngrids, ad_VolCons(iwest,:))
              DO ng=1,Ngrids
                tl_VolCons(iwest,ng)=ad_VolCons(iwest,ng)
              END DO
            CASE ('ad_VolCons(east)')
              Npts=load_l(Nval, Cval, Ngrids, ad_VolCons(ieast,:))
              DO ng=1,Ngrids
                tl_VolCons(ieast,ng)=ad_VolCons(ieast,ng)
              END DO
            CASE ('ad_VolCons(south)')
              Npts=load_l(Nval, Cval, Ngrids, ad_VolCons(isouth,:))
              DO ng=1,Ngrids
                tl_VolCons(isouth,ng)=ad_VolCons(isouth,ng)
              END DO
            CASE ('ad_VolCons(north)')
              Npts=load_l(Nval, Cval, Ngrids, ad_VolCons(inorth,:))
              DO ng=1,Ngrids
                tl_VolCons(inorth,ng)=ad_VolCons(inorth,ng)
              END DO
#endif
            CASE ('NTIMES')
              Npts=load_i(Nval, Rval, Ngrids, ntimes)
            CASE ('DT')
              Npts=load_r(Nval, Rval, Ngrids, dt)
#ifdef MODEL_COUPLING
              IF (.not.allocated(CoupleSteps)) THEN
                allocate (CoupleSteps(Nmodels,Ngrids))
                CoupleSteps=0
              END IF
              DO ng=1,Ngrids
                DO i=1,Nmodels
                  CoupleSteps(i,ng)=MAX(1,                              &
     &                                  INT(TimeInterval(Iocean,i)/     &
     &                                  dt(ng)))
                END DO
              END DO
#endif
            CASE ('NDTFAST')
              Npts=load_i(Nval, Rval, Ngrids, ndtfast)
            CASE ('ERstr')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ERstr=Ivalue(1)
            CASE ('ERend')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              ERend=Ivalue(1)
            CASE ('Nouter')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nouter=Ivalue(1)
            CASE ('Ninner')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Ninner=Ivalue(1)
            CASE ('Nintervals')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nintervals=Ivalue(1)
#ifdef PROPAGATOR
            CASE ('NEV')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NEV=Ivalue(1)
            CASE ('NCV')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              NCV=Ivalue(1)
# if defined FT_EIGENMMODES || defined AFT_EIGENMODES
              IF (NCV.lt.(2*NEV+1)) THEN
                IF (Master) WRITE (out,260) 'NCV = ', NCV,              &
     &            'Must be greater than or equal to 2*NEV+1'
                exit_flag=5
                RETURN
              END IF
# else
              IF (NCV.lt.NEV) THEN
                IF (Master) WRITE (out,260) 'NCV = ', NCV,              &
     &            'Must be greater than NEV'
                exit_flag=5
                RETURN
              END IF
# endif
#endif
            CASE ('NRREC')
              Npts=load_i(Nval, Rval, Ngrids, nrrec)
              DO ng=1,Ngrids
                IF (nrrec(ng).lt.0) THEN
                  LastRec(ng)=.TRUE.
                ELSE
                  LastRec(ng)=.FALSE.
                END IF
              END DO
            CASE ('LcycleRST')
              Npts=load_l(Nval, Cval, Ngrids, LcycleRST)
            CASE ('NRST')
              Npts=load_i(Nval, Rval, Ngrids, nRST)
            CASE ('NSTA')
              Npts=load_i(Nval, Rval, Ngrids, nSTA)
            CASE ('NFLT')
              Npts=load_i(Nval, Rval, Ngrids, nFLT)
            CASE ('NINFO')
              Npts=load_i(Nval, Rval, Ngrids, ninfo)
              DO ng=1,Ngrids
                IF (ninfo(ng).le.0) THEN
                  WRITE (text,'(a,i2.2,a)') 'ninfo(', ng, ') = '
                  IF (Master) WRITE (out,260) TRIM(text), ninfo(ng),    &
     &               'must be greater than zero.'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('LDEFOUT')
              Npts=load_l(Nval, Cval, Ngrids, ldefout)
            CASE ('NHIS')
              Npts=load_i(Nval, Rval, Ngrids, nHIS)
            CASE ('NDEFHIS')
              Npts=load_i(Nval, Rval, Ngrids, ndefHIS)
            CASE ('NHIS2')
              Npts=load_i(Nval, Rval, Ngrids, nHIS2)
            CASE ('NDEFHIS2')
              Npts=load_i(Nval, Rval, Ngrids, ndefHIS2)
            CASE ('NTSAVG')
              Npts=load_i(Nval, Rval, Ngrids, ntsAVG)
#ifdef ADJOINT
              DO ng=1,Ngrids
                IF (ntsAVG(ng).eq.1) ntsAVG(ng)=ntimes(ng)
              END DO
#endif
            CASE ('NAVG')
              Npts=load_i(Nval, Rval, Ngrids, nAVG)
            CASE ('NDEFAVG')
              Npts=load_i(Nval, Rval, Ngrids, ndefAVG)
            CASE ('NTSAVG2')
              Npts=load_i(Nval, Rval, Ngrids, ntsAVG2)
            CASE ('NAVG2')
              Npts=load_i(Nval, Rval, Ngrids, nAVG2)
            CASE ('NDEFAVG2')
              Npts=load_i(Nval, Rval, Ngrids, ndefAVG2)
            CASE ('NTSDIA')
              Npts=load_i(Nval, Rval, Ngrids, ntsDIA)
            CASE ('NDIA')
              Npts=load_i(Nval, Rval, Ngrids, nDIA)
            CASE ('NDEFDIA')
              Npts=load_i(Nval, Rval, Ngrids, ndefDIA)
            CASE ('LcycleTLM')
              Npts=load_l(Nval, Cval, Ngrids, LcycleTLM)
            CASE ('NTLM')
              Npts=load_i(Nval, Rval, Ngrids, nTLM)
            CASE ('NDEFTLM')
              Npts=load_i(Nval, Rval, Ngrids, ndefTLM)
            CASE ('LcycleADJ')
              Npts=load_l(Nval, Cval, Ngrids, LcycleADJ)
            CASE ('NADJ')
              Npts=load_i(Nval, Rval, Ngrids, nADJ)
            CASE ('NDEFADJ')
              Npts=load_i(Nval, Rval, Ngrids, ndefADJ)
            CASE ('NOBC')
              Npts=load_i(Nval, Rval, Ngrids, nOBC)
#ifdef ADJUST_BOUNDARY
              DO ng=1,Ngrids
                Nbrec(ng)=1+ntimes(ng)/nOBC(ng)
              END DO
              allocate ( OBC_time(MAXVAL(Nbrec),Ngrids) )
#endif
            CASE ('NSFF')
              Npts=load_i(Nval, Rval, Ngrids, nSFF)
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
              DO ng=1,Ngrids
                Nfrec(ng)=1+ntimes(ng)/nSFF(ng)
              END DO
              allocate ( SF_time(MAXVAL(Nfrec),Ngrids) )
#endif
            CASE ('LmultiGST')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LmultiGST=Lvalue(1)
            CASE ('LrstGST')
              Npts=load_l(Nval, Cval, 1, Lvalue)
              LrstGST=Lvalue(1)
            CASE ('MaxIterGST')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              MaxIterGST=Ivalue(1)
            CASE ('NGST')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              nGST=Ivalue(1)
#ifdef PROPAGATOR
            CASE ('Ritz_tol')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              Ritz_tol=Rvalue(1)
#endif
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  nl_tnu2(itrc,ng)=Rtracer(itrc,ng)
                END DO
#ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  nl_tnu2(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
#endif
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  nl_tnu4(itrc,ng)=Rtracer(itrc,ng)
                END DO
#ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  nl_tnu4(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
#endif
              END DO
            CASE ('ad_TNU2')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_tnu2(itrc,ng)=Rtracer(itrc,ng)
                  tl_tnu2(itrc,ng)=Rtracer(itrc,ng)
                END DO
#ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  ad_tnu2(itrc,ng)=Rtracer(NAT+i,ng)
                  tl_tnu2(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
#endif
              END DO
            CASE ('ad_TNU4')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_tnu4(itrc,ng)=Rtracer(itrc,ng)
                  tl_tnu4(itrc,ng)=Rtracer(itrc,ng)
                END DO
#ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  ad_tnu4(itrc,ng)=Rtracer(NAT+i,ng)
                  tl_tnu4(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
#endif
              END DO
            CASE ('VISC2')
              Npts=load_r(Nval, Rval, Ngrids, nl_visc2)
            CASE ('VISC4')
              Npts=load_r(Nval, Rval, Ngrids, nl_visc4)
            CASE ('ad_VISC2')
              Npts=load_r(Nval, Rval, Ngrids, ad_visc2)
              DO ng=1,Ngrids
                tl_visc2(ng)=ad_visc2(ng)
              END DO
            CASE ('ad_VISC4')
              Npts=load_r(Nval, Rval, Ngrids, ad_visc4)
              DO ng=1,Ngrids
                tl_visc4(ng)=ad_visc4(ng)
              END DO
            CASE ('LuvSponge')
              Npts=load_l(Nval, Cval, Ngrids, LuvSponge)
#ifdef SOLVE3D
            CASE ('LtracerSponge')
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LtracerSponge(itrc,ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerSponge(itrc,ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
#endif
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Akt_bak(itrc,ng)=Rtracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Akt_bak(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('ad_AKT_fac')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  ad_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                  tl_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  ad_Akt_fac(itrc,ng)=Rtracer(NAT+i,ng)
                  tl_Akt_fac(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('AKV_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
            CASE ('ad_AKV_fac')
              Npts=load_r(Nval, Rval, Ngrids, ad_Akv_fac)
              DO ng=1,Ngrids
                tl_Akv_fac(ng)=ad_AKv_fac(ng)
              END DO
            CASE ('AKK_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akk_bak)
            CASE ('AKP_BAK')
              Npts=load_r(Nval, Rval, Ngrids, Akp_bak)
            CASE ('TKENU2')
              Npts=load_r(Nval, Rval, Ngrids, tkenu2)
            CASE ('TKENU4')
              Npts=load_r(Nval, Rval, Ngrids, tkenu4)
            CASE ('GLS_P')
              Npts=load_r(Nval, Rval, Ngrids, gls_p)
            CASE ('GLS_M')
              Npts=load_r(Nval, Rval, Ngrids, gls_m)
            CASE ('GLS_N')
              Npts=load_r(Nval, Rval, Ngrids, gls_n)
            CASE ('GLS_Kmin')
              Npts=load_r(Nval, Rval, Ngrids, gls_Kmin)
            CASE ('GLS_Pmin')
              Npts=load_r(Nval, Rval, Ngrids, gls_Pmin)
            CASE ('GLS_CMU0')
              Npts=load_r(Nval, Rval, Ngrids, gls_cmu0)
            CASE ('GLS_C1')
              Npts=load_r(Nval, Rval, Ngrids, gls_c1)
            CASE ('GLS_C2')
              Npts=load_r(Nval, Rval, Ngrids, gls_c2)
            CASE ('GLS_C3M')
              Npts=load_r(Nval, Rval, Ngrids, gls_c3m)
            CASE ('GLS_C3P')
              Npts=load_r(Nval, Rval, Ngrids, gls_c3p)
            CASE ('GLS_SIGK')
              Npts=load_r(Nval, Rval, Ngrids, gls_sigk)
            CASE ('GLS_SIGP')
              Npts=load_r(Nval, Rval, Ngrids, gls_sigp)
            CASE ('CHARNOK_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, charnok_alpha)
            CASE ('ZOS_HSIG_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, zos_hsig_alpha)
            CASE ('SZ_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, sz_alpha)
            CASE ('CRGBAN_CW')
              Npts=load_r(Nval, Rval, Ngrids, crgban_cw)
            CASE ('WEC_ALPHA')
              Npts=load_r(Nval, Rval, Ngrids, wec_alpha)
            CASE ('RDRG')
              Npts=load_r(Nval, Rval, Ngrids, rdrg)
            CASE ('RDRG2')
              Npts=load_r(Nval, Rval, Ngrids, rdrg2)
            CASE ('Zob')
              Npts=load_r(Nval, Rval, Ngrids, Zob)
            CASE ('Zos')
              Npts=load_r(Nval, Rval, Ngrids, Zos)
#ifdef BULK_FLUXES
            CASE ('BLK_ZQ')
              Npts=load_r(Nval, Rval, Ngrids, blk_ZQ)
            CASE ('BLK_ZT')
              Npts=load_r(Nval, Rval, Ngrids, blk_ZT)
            CASE ('BLK_ZW')
              Npts=load_r(Nval, Rval, Ngrids, blk_ZW)
#endif
            CASE ('DCRIT')
              Npts=load_r(Nval, Rval, Ngrids, Dcrit)
            CASE ('WTYPE')
              Npts=load_i(Nval, Rval, Ngrids, lmd_Jwt)
            CASE ('LEVSFRC')
              Npts=load_i(Nval, Rval, Ngrids, levsfrc)
            CASE ('LEVBFRC')
              Npts=load_i(Nval, Rval, Ngrids, levbfrc)
            CASE ('Vtransform')
              Npts=load_i(Nval, Rval, Ngrids, Vtransform)
              DO ng=1,Ngrids
                IF ((Vtransform(ng).lt.0).or.                           &
     &              (Vtransform(ng).gt.2)) THEN
                  IF (Master) WRITE (out,260) 'Vtransform = ',          &
     &                                        Vtransform(ng),           &
     &                                        'Must be either 1 or 2'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('Vstretching')
              Npts=load_i(Nval, Rval, Ngrids, Vstretching)
              DO ng=1,Ngrids
                IF ((Vstretching(ng).lt.0).or.                          &
     &              (Vstretching(ng).gt.4)) THEN
                  IF (Master) WRITE (out,260) 'Vstretching = ',         &
     &                                        Vstretching(ng),          &
     &                                        'Must between 1 and 4'
                  exit_flag=5
                  RETURN
                END IF
              END DO
            CASE ('THETA_S')
              Npts=load_r(Nval, Rval, Ngrids, theta_s)
            CASE ('THETA_B')
              Npts=load_r(Nval, Rval, Ngrids, theta_b)
            CASE ('TCLINE')
              Npts=load_r(Nval, Rval, Ngrids, Tcline)
              DO ng=1,Ngrids
                hc(ng)=Tcline(ng)
              END DO
            CASE ('RHO0')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              rho0=Rvalue(1)
            CASE ('BVF_BAK')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              bvf_bak=Rvalue(1)
            CASE ('DSTART')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              dstart=Rvalue(1)
            CASE ('TIDE_START')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              tide_start=Rvalue(1)
            CASE ('TIME_REF')
              Npts=load_r(Nval, Rval, 1, Rvalue)
              time_ref=Rvalue(1)
              r_text=ref_att(time_ref,r_date)
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Tnudg(itrc,ng)=Rtracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Tnudg(itrc,ng)=Rtracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('ZNUDG')
              Npts=load_r(Nval, Rval, Ngrids, Znudg)
            CASE ('M2NUDG')
              Npts=load_r(Nval, Rval, Ngrids, M2nudg)
            CASE ('M3NUDG')
              Npts=load_r(Nval, Rval, Ngrids, M3nudg)
            CASE ('TNUDG_SSS')
              Npts=load_r(Nval, Rval, Ngrids, Tnudg_SSS)
            CASE ('SSS_MISMATCH_THRESHOLD')
              Npts=load_r(Nval, Rval, Ngrids, sss_mismatch_threshold)
            CASE ('OBCFAC')
              Npts=load_r(Nval, Rval, Ngrids, obcfac)
            CASE ('R0')
              Npts=load_r(Nval, Rval, Ngrids, R0)
              DO ng=1,Ngrids
                IF (R0(ng).lt.100.0_r8) R0(ng)=R0(ng)+1000.0_r8
              END DO
            CASE ('T0')
              Npts=load_r(Nval, Rval, Ngrids, T0)
            CASE ('S0')
              Npts=load_r(Nval, Rval, Ngrids, S0)
            CASE ('TCOEF')
              Npts=load_r(Nval, Rval, Ngrids, Tcoef)
              DO ng=1,Ngrids
                Tcoef(ng)=ABS(Tcoef(ng))
              END DO
            CASE ('SCOEF')
              Npts=load_r(Nval, Rval, Ngrids, Scoef)
              DO ng=1,Ngrids
                Scoef(ng)=ABS(Scoef(ng))
              END DO
            CASE ('GAMMA2')
              Npts=load_r(Nval, Rval, Ngrids, gamma2)
            CASE ('LuvSrc')
              Npts=load_l(Nval, Cval, Ngrids, LuvSrc)
            CASE ('LwSrc')
              Npts=load_l(Nval, Cval, Ngrids, LwSrc)
#ifdef SOLVE3D
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LtracerSrc(itrc,ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerSrc(itrc,ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
#endif
            CASE ('LsshCLM')
              Npts=load_l(Nval, Cval, Ngrids, LsshCLM)
            CASE ('Lm2CLM')
              Npts=load_l(Nval, Cval, Ngrids, Lm2CLM)
#ifdef SOLVE3D
            CASE ('Lm3CLM')
              Npts=load_l(Nval, Cval, Ngrids, Lm3CLM)
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LtracerCLM(itrc,ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  LtracerCLM(itrc,ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
#endif
            CASE ('LnudgeM2CLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeM2CLM)
#ifdef SOLVE3D
            CASE ('LnudgeM3CLM')
              Npts=load_l(Nval, Cval, Ngrids, LnudgeM3CLM)
            CASE ('LnudgeTCLM')
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  LnudgeTCLM(itrc,ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  LnudgeTCLM(itrc,ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
#endif
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
            CASE ('DstrS')
              Npts=load_r(Nval, Rval, Ngrids, DstrS)
            CASE ('DendS')
              Npts=load_r(Nval, Rval, Ngrids, DendS)
            CASE ('KstrS')
              Npts=load_i(Nval, Rval, Ngrids, KstrS)
            CASE ('KendS')
              Npts=load_i(Nval, Rval, Ngrids, KendS)
            CASE ('Lstate(isFsur)')
              IF (isFsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'isFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Lstate(isFsur)=Ladsen(ng)
              END DO
            CASE ('Lstate(isUbar)')
              IF (isUbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'isUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Lstate(isUbar)=Ladsen(ng)
              END DO
            CASE ('Lstate(isVbar)')
              IF (isVbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'isVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Lstate(isVbar)=Ladsen(ng)
              END DO
# ifdef SOLVE3D
            CASE ('Lstate(isUvel)')
              IF (isUvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'isUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Lstate(isUvel)=Ladsen(ng)
              END DO
            CASE ('Lstate(isVvel)')
              IF (isVvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'isVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Lstate(isVvel)=Ladsen(ng)
              END DO
            CASE ('Lstate(isTvar)')
              IF (MAXVAL(isTvar).eq.0) THEN
                IF (Master) WRITE (out,280) 'isTvar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, MT*Ngrids, Ladsen)
              k=0
              DO ng=1,Ngrids
                DO itrc=1,MT
                  k=k+1
                  i=isTvar(itrc)
                  SCALARS(ng)%Lstate(i)=Ladsen(k)
                END DO
              END DO
# endif
#endif
#if defined FORCING_SV     || defined SO_SEMI || \
    defined STOCHASTIC_OPT
            CASE ('Fstate(isFsur)')
              IF (isFsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'isFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isFsur)=Ladsen(ng)
              END DO
            CASE ('Fstate(isUbar)')
              IF (isUbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'isUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isUbar)=Ladsen(ng)
              END DO
            CASE ('Fstate(isVbar)')
              IF (isVbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'isVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isVbar)=Ladsen(ng)
              END DO
# ifdef SOLVE3D
            CASE ('Fstate(isUvel)')
              IF (isUvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'isUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isUvel)=Ladsen(ng)
              END DO
            CASE ('Fstate(isVvel)')
              IF (isVvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'isVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isVvel)=Ladsen(ng)
              END DO
            CASE ('Fstate(isTvar)')
              IF (MAXVAL(isTvar).eq.0) THEN
                IF (Master) WRITE (out,280) 'isTvar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, MT*Ngrids, Ladsen)
              k=0
              DO ng=1,Ngrids
                DO itrc=1,MT
                  k=k+1
                  i=isTvar(itrc)
                  SCALARS(ng)%Fstate(i)=Ladsen(k)
                END DO
              END DO
# endif
            CASE ('Fstate(isUstr)')
              IF (isUstr.eq.0) THEN
                IF (Master) WRITE (out,280) 'isUstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isUstr)=Ladsen(ng)
              END DO
            CASE ('Fstate(isVstr)')
              IF (isUstr.eq.0) THEN
                IF (Master) WRITE (out,280) 'isVstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ladsen)
              DO ng=1,Ngrids
                SCALARS(ng)%Fstate(isVstr)=Ladsen(ng)
              END DO
# ifdef SOLVE3D
            CASE ('Fstate(isTsur)')
              IF (MAXVAL(isTsur).eq.0) THEN
                IF (Master) WRITE (out,280) 'isTsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, MT*Ngrids, Ladsen)
              k=0
              DO ng=1,Ngrids
                DO itrc=1,MT
                  k=k+1
                  i=isTsur(itrc)
                  SCALARS(ng)%Fstate(i)=Ladsen(k)
                END DO
              END DO
# endif
#endif
#if defined SO_SEMI        || \
   (defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE)
            CASE ('SO_decay')
              Npts=load_r(Nval, Rval, Ngrids, SO_decay)
            CASE ('SO_sdev(isFsur)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isFsur,1))
            CASE ('SO_sdev(isUbar)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isUbar,1))
            CASE ('SO_sdev(isVbar)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isVbar,1))
# ifdef SOLVE3D
            CASE ('SO_sdev(isUvel)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isUvel,1))
            CASE ('SO_sdev(isVvel)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isVvel,1))
            CASE ('SO_sdev(isTvar)')
              Npts=load_r(Nval, Rval, MT*Ngrids, tracer)
              k=0
              DO ng=1,Ngrids
                DO itrc=1,MT
                  k=k+1
                  i=isTvar(itrc)
                  SO_sdev(i,ng)=tracer(k,ng)
                END DO
              END DO
# endif
            CASE ('SO_sdev(isUstr)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isUstr,1))
            CASE ('SO_sdev(isVstr)')
              Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isVstr,1))
# ifdef SOLVE3D
            CASE ('SO_sdev(isTsur)')
              Npts=load_r(Nval, Rval, MT*Ngrids, tracer)
              k=0
              DO ng=1,Ngrids
                DO itrc=1,MT
                  k=k+1
                  i=isTsur(itrc)
                  SO_sdev(i,ng)=tracer(k,ng)
                END DO
              END DO
# endif
#endif
            CASE ('Hout(idUvel)')
              IF (idUvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUvel,:))
            CASE ('Hout(idVvel)')
              IF (idVvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVvel,:))
            CASE ('Hout(idWvel)')
              IF (idWvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWvel,:))
            CASE ('Hout(idOvel)')
              IF (idOvel.eq.0) THEN
                IF (Master) WRITE (out,280) 'idOvel'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idOvel,:))
            CASE ('Hout(idUbar)')
              IF (idUbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbar,:))
            CASE ('Hout(idVbar)')
              IF (idVbar.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbar,:))
            CASE ('Hout(idFsur)')
              IF (idFsur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idFsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idFsur,:))
#if defined SEDIMENT && defined SED_MORPH
            CASE ('Hout(idBath)')
              IF (idbath.eq.0) THEN
                IF (Master) WRITE (out,280) 'idbath'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idBath,:))
#endif
            CASE ('Hout(idu2dE)')
              IF (idu2dE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu2dE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idu2dE,:))
            CASE ('Hout(idv2dN)')
              IF (idv2dN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv2dN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idv2dN,:))
            CASE ('Hout(idu3dE)')
              IF (idu3dE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu3dE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idu3dE,:))
            CASE ('Hout(idv3dN)')
              IF (idv3dN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv3dN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idv3dN,:))
            CASE ('Hout(idTvar)')
              IF (MAXVAL(idTvar).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTvar'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Hout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idpthR)')
              IF (idpthR.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthR'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idpthR,:))
            CASE ('Hout(idpthU)')
              IF (idpthU.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthU'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idpthU,:))
            CASE ('Hout(idpthV)')
              IF (idpthV.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idpthV,:))
            CASE ('Hout(idpthW)')
              IF (idpthW.eq.0) THEN
                IF (Master) WRITE (out,280) 'idpthW'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idpthW,:))
            CASE ('Hout(idUsms)')
              IF (idUsms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUsms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUsms,:))
            CASE ('Hout(idVsms)')
              IF (idVsms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVsms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVsms,:))
            CASE ('Hout(idUbms)')
              IF (idUbms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbms,:))
            CASE ('Hout(idVbms)')
              IF (idVbms.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbms'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbms,:))
#ifdef BBL_MODEL
            CASE ('Hout(idUbrs)')
              IF (idUbrs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbrs,:))
            CASE ('Hout(idVbrs)')
              IF (idVbrs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbrs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbrs,:))
            CASE ('Hout(idUbws)')
              IF (idUbws.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbws'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbws,:))
            CASE ('Hout(idVbws)')
              IF (idVbws.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbws'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbws,:))
            CASE ('Hout(idUbcs)')
              IF (idUbcs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbcs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbcs,:))
            CASE ('Hout(idVbcs)')
              IF (idVbcs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbcs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbcs,:))
            CASE ('Hout(idUVwc)')
              IF (idUVwc.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUVwc'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUVwc,:))
            CASE ('Hout(idUbot)')
              IF (idUbot.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbot,:))
            CASE ('Hout(idVbot)')
              IF (idVbot.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbot'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbot,:))
            CASE ('Hout(idUbur)')
              IF (idUbur.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUbur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUbur,:))
            CASE ('Hout(idVbvr)')
              IF (idVbvr.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVbvr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVbvr,:))
#endif
#ifdef WEC_MELLOR
            CASE ('Hout(idW2xx)')
              IF (idW2xx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW2xx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xx,:))
            CASE ('Hout(idW2xy)')
              IF (idW2xy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW2xy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xy,:))
            CASE ('Hout(idW2yy)')
              IF (idW2yy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW2yy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW2yy,:))
            CASE ('Hout(idW3xx)')
              IF (idW3xx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3xx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xx,:))
            CASE ('Hout(idW3xy)')
              IF (idW3xy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3xy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xy,:))
            CASE ('Hout(idW3yy)')
              IF (idW3yy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3yy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3yy,:))
            CASE ('Hout(idW3zx)')
              IF (idW3zx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3zx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zx,:))
            CASE ('Hout(idW3zy)')
              IF (idW3zy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3zy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zy,:))
#endif
#ifdef WEC_VF
            CASE ('Hout(idWztw)')
              IF (idWztw.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWztw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWztw,:))
            CASE ('Hout(idWqsp)')
              IF (idWqsp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWqsp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWqsp,:))
            CASE ('Hout(idWbeh)')
              IF (idWbeh.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWbeh'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWbeh,:))
#endif
#ifdef WEC
            CASE ('Hout(idU2rs)')
              IF (idU2rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU2rs,:))
            CASE ('Hout(idV2rs)')
              IF (idV2rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV2rs,:))
            CASE ('Hout(idU3rs)')
              IF (idU3rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU3rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU3rs,:))
            CASE ('Hout(idV3rs)')
              IF (idV3rs.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV3rs'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV3rs,:))
            CASE ('Hout(idU2Sd)')
              IF (idU2Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU2Sd,:))
            CASE ('Hout(idV2Sd)')
              IF (idV2Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV2Sd,:))
            CASE ('Hout(idU3Sd)')
              IF (idU3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idU3Sd,:))
            CASE ('Hout(idV3Sd)')
              IF (idV3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idV3Sd,:))
            CASE ('Hout(idW3Sd)')
              IF (idW3Sd.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3Sd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3Sd,:))
            CASE ('Hout(idW3St)')
              IF (idW3St.eq.0) THEN
                IF (Master) WRITE (out,280) 'idW3St'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idW3St,:))
#endif
#ifdef WAVES_HEIGHT
            CASE ('Hout(idWamp)')
              IF (idWamp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWamp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWamp,:))
#endif
#ifdef WAVES_LENGTH

            CASE ('Hout(idWlen)')
              IF (idWlen.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWlen'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWlen,:))
#endif
#ifdef WAVES_LENGTHP
            CASE ('Hout(idWlep)')
              IF (idWlep.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWlep'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWlep,:))
#endif
#ifdef WAVES_DIR
            CASE ('Hout(idWdir)')
              IF (idWdir.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdir'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdir,:))
#endif
#ifdef WAVES_TOP_PERIOD
            CASE ('Hout(idWptp)')
              IF (idWptp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWptp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWptp,:))
#endif
#ifdef WAVES_BOT_PERIOD
            CASE ('Hout(idWpbt)')
              IF (idWpbt.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWpbt'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWpbt,:))
#endif
#ifdef WAVES_OCEAN
            CASE ('Hout(idWorb)')
              IF (idWorb.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWorb'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWorb,:))
#endif
#if defined ROLLER_SVENDSEN
            CASE ('Hout(idWbrk)')
              IF (idWbrk.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWbrk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWbrk,:))
#endif
#if defined UV_KIRBY
            CASE ('Hout(idUwav)')
              IF (idUwav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUwav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUwav,:))
            CASE ('Hout(idVwav)')
              IF (idVwav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVwav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVwav,:))
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
            CASE ('Hout(idWdif)')
              IF (idWdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdif,:))
#endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR || \
    defined WDISS_INWAVE
            CASE ('Hout(idWdib)')
              IF (idWdib.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdib'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdib,:))
            CASE ('Hout(idWdiw)')
              IF (idWdiw.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdiw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdiw,:))
#endif
#if defined WEC_ROLLER
            CASE ('Hout(idWdis)')
              IF (idWdis.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWdis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWdis,:))
            CASE ('Hout(idWrol)')
              IF (idWrol.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWrol'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWrol,:))
#endif
#if defined WAVES_DSPR
            CASE ('Hout(idWvds)')
              IF (idWvds.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWvds'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWvds,:))
            CASE ('Hout(idWvqp)')
              IF (idWvqp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idWvqp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idWvqp,:))
#endif
#if defined INWAVE_MODEL
            CASE ('Hout(idACen)')
              IF (idACen.eq.0) THEN
                IF (Master) WRITE (out,280) 'idACen'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idACen,:))
            CASE ('Hout(idACct)')
              IF (idACct.eq.0) THEN
                IF (Master) WRITE (out,280) 'idACct'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idACct,:))
            CASE ('Hout(idACcx)')
              IF (idACcx.eq.0) THEN
                IF (Master) WRITE (out,280) 'idACcx'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idACcx,:))
            CASE ('Hout(idACcy)')
              IF (idACcy.eq.0) THEN
                IF (Master) WRITE (out,280) 'idACcy'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idACcy,:))
            CASE ('Hout(idACtp)')
              IF (idACtp.eq.0) THEN
                IF (Master) WRITE (out,280) 'idACtp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idACtp,:))
#endif
#ifdef SOLVE3D
# if defined RUNOFF
            CASE ('Hout(idRunoff)')
              IF (idRunoff.eq.0) THEN
                IF (Master) WRITE (out,280) 'idRunoff'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idRunoff,:))
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
  || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
            CASE ('Hout(idPair)')
              IF (idPair.eq.0) THEN
                IF (Master) WRITE (out,280) 'idPair'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idPair,:))
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM2OCN_FLUXES
            CASE ('Hout(idUair)')
              IF (idUair.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUair'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUair,:))
            CASE ('Hout(idVair)')
              IF (idVair.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVair'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVair,:))
            CASE ('Hout(idUairE)')
              IF (idUairE.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUairE'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idUairE,:))
            CASE ('Hout(idVairN)')
              IF (idVairN.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVairN'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVairN,:))
# endif
#endif
            CASE ('Hout(idTsur)')
              IF (idTsur(itemp).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTsur'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Hout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idLhea)')
              IF (idLhea.eq.0) THEN
                IF (Master) WRITE (out,280) 'idLhea'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idLhea,:))
            CASE ('Hout(idShea)')
              IF (idShea.eq.0) THEN
                IF (Master) WRITE (out,280) 'idShea'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idShea,:))
            CASE ('Hout(idLrad)')
              IF (idLrad.eq.0) THEN
                IF (Master) WRITE (out,280) 'idLrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idLrad,:))
            CASE ('Hout(idSrad)')
              IF (idSrad.eq.0) THEN
                IF (Master) WRITE (out,280) 'idSrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,:))
            CASE ('Hout(idEmPf)')
              IF (idEmPf.eq.0) THEN
                IF (Master) WRITE (out,280) 'idEmPf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idEmPf,:))
            CASE ('Hout(idevap)')
              IF (idevap.eq.0) THEN
                IF (Master) WRITE (out,280) 'idevap'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idevap,:))
            CASE ('Hout(idrain)')
              IF (idrain.eq.0) THEN
                IF (Master) WRITE (out,280) 'idrain'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idrain,:))
            CASE ('Hout(idDano)')
              IF (idDano.eq.0) THEN
                IF (Master) WRITE (out,280) 'idDano'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idDano,:))
            CASE ('Hout(idVvis)')
              IF (idVvis.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idVvis,:))
            CASE ('Hout(idTdif)')
              IF (idTdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idTdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idTdif,:))
            CASE ('Hout(idSdif)')
              IF (idSdif.eq.0) THEN
                IF (Master) WRITE (out,280) 'idSdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idSdif,:))
            CASE ('Hout(idHsbl)')
              IF (idHsbl.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHsbl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHsbl,:))
            CASE ('Hout(idHbbl)')
              IF (idHbbl.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHbbl'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idHbbl,:))
            CASE ('Hout(idMtke)')
              IF (idMtke.eq.0) THEN
                IF (Master) WRITE (out,280) 'idMtke'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idMtke,:))
            CASE ('Hout(idMtls)')
              IF (idMtls.eq.0) THEN
                IF (Master) WRITE (out,280) 'idMtls'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idMtls,:))
#if defined SOLVE3D && defined T_PASSIVE
            CASE ('Hout(inert)')
              Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=idTvar(inert(i))
                  Hout(itrc,ng)=Linert(i,ng)
                END DO
              END DO
#endif
#if defined HISTORY2
            CASE ('Hout2(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUvel,:))
            CASE ('Hout2(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVvel,:))
            CASE ('Hout2(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idWvel,:))
            CASE ('Hout2(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idOvel,:))
            CASE ('Hout2(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUbar,:))
            CASE ('Hout2(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVbar,:))
            CASE ('Hout2(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idFsur,:))
            CASE ('Hout2(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idu2dE,:))
            CASE ('Hout2(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idv2dN,:))
            CASE ('Hout2(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idu3dE,:))
            CASE ('Hout2(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idv3dN,:))
            CASE ('Hout2(idTvar)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout2(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUsms,:))
            CASE ('Hout2(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVsms,:))
            CASE ('Hout2(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUbms,:))
            CASE ('Hout2(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVbms,:))
            CASE ('Hout2(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW2xx,:))
            CASE ('Hout2(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW2xy,:))
            CASE ('Hout2(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW2yy,:))
            CASE ('Hout2(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idU2rs,:))
            CASE ('Hout2(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idV2rs,:))
            CASE ('Hout2(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idU2Sd,:))
            CASE ('Hout2(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idV2Sd,:))
            CASE ('Hout2(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW3xx,:))
            CASE ('Hout2(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW3xy,:))
            CASE ('Hout2(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW3yy,:))
            CASE ('Hout2(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW3zx,:))
            CASE ('Hout2(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idW3zy,:))
            CASE ('Hout2(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idU3rs,:))
            CASE ('Hout2(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idV3rs,:))
            CASE ('Hout2(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idU3Sd,:))
            CASE ('Hout2(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idV3Sd,:))
# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
            CASE ('Hout2(idPair)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idPair,:))
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
            CASE ('Hout2(idUair)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUair,:))
            CASE ('Hout2(idVair)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVair,:))
            CASE ('Hout2(idUairE)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUairE,:))
            CASE ('Hout2(idVairN)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVairN,:))
#  endif
# endif
            CASE ('Hout2(idTsur)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout2(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idLhea,:))
            CASE ('Hout2(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idShea,:))
            CASE ('Hout2(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idLrad,:))
            CASE ('Hout2(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idSrad,:))
            CASE ('Hout2(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idevap,:))
            CASE ('Hout2(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idrain,:))
            CASE ('Hout2(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idDano,:))
            CASE ('Hout2(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVvis,:))
            CASE ('Hout2(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idTdif,:))
            CASE ('Hout2(idSdif)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idSdif,:))
            CASE ('Hout2(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHsbl,:))
            CASE ('Hout2(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHbbl,:))
            CASE ('Hout2(id2dRV)')
              IF (id2dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(id2dRV,:))
            CASE ('Hout2(id3dRV)')
              IF (id3dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(id3dRV,:))
            CASE ('Hout2(id2dPV)')
              IF (id2dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(id2dPV,:))
            CASE ('Hout2(id3dPV)')
              IF (id3dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(id3dPV,:))
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
            CASE ('Hout2(idFsuD)')
              IF (idFsuD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idFsuD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idFsuD,:))
            CASE ('Hout2(idu2dD)')
              IF (idu2dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu2dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idu2dD,:))
            CASE ('Hout2(idv2dD)')
              IF (idv2dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv2dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idv2dD,:))
#  ifdef SOLVE3D
            CASE ('Hout2(idu3dD)')
              IF (idu3dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu3dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idu3dD,:))
            CASE ('Hout2(idv3dD)')
              IF (idv3dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv3dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idv3dD,:))
            CASE ('Hout2(idTrcD)')
              IF (MAXVAL(idTrcD).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTrcD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTrcD(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
#  endif
# endif
# ifdef SOLVE3D
            CASE ('Hout2(idHUav)')
              IF (idHUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHUav,:))
            CASE ('Hout2(idHVav)')
              IF (idHVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idHVav,:))
            CASE ('Hout2(idUUav)')
              IF (idUUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUUav,:))
            CASE ('Hout2(idUVav)')
              IF (idUVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idUVav,:))
            CASE ('Hout2(idVVav)')
              IF (idVVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idVVav,:))
# endif
            CASE ('Hout2(idU2av)')
              IF (idU2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idU2av,:))
            CASE ('Hout2(idV2av)')
              IF (idV2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idV2av,:))
            CASE ('Hout2(idZZav)')
              IF (idZZav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idZZav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout2(idZZav,:))
# ifdef SOLVE3D
            CASE ('Hout2(idTTav)')
              IF (MAXVAL(idTTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTTav(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout2(idUTav)')
              IF (MAXVAL(idUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idUTav(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout2(idVTav)')
              IF (MAXVAL(idVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idVTav(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout2(iHUTav)')
              IF (MAXVAL(iHUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=iHUTav(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Hout2(iHVTav)')
              IF (MAXVAL(iHVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=iHVTav(itrc)
                  Hout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
# endif
# if defined SOLVE3D && defined T_PASSIVE
            CASE ('Hout2(inert)')
              Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=idTvar(inert(i))
                  Hout2(itrc,ng)=Linert(i,ng)
                END DO
              END DO
# endif
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
            CASE ('Aout(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUvel,:))
            CASE ('Aout(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVvel,:))
            CASE ('Aout(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWvel,:))
            CASE ('Aout(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idOvel,:))
            CASE ('Aout(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbar,:))
            CASE ('Aout(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbar,:))
            CASE ('Aout(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idFsur,:))
            CASE ('Aout(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idu2dE,:))
            CASE ('Aout(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idv2dN,:))
            CASE ('Aout(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idu3dE,:))
            CASE ('Aout(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idv3dN,:))
            CASE ('Aout(idTvar)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUsms,:))
            CASE ('Aout(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVsms,:))
            CASE ('Aout(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbms,:))
            CASE ('Aout(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbms,:))
# ifdef BBL_MODEL
            CASE ('Aout(idUbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbrs,:))
            CASE ('Aout(idVbrs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbrs,:))
            CASE ('Aout(idUbws)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbws,:))
            CASE ('Aout(idVbws)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbws,:))
            CASE ('Aout(idUbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbcs,:))
            CASE ('Aout(idVbcs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbcs,:))
            CASE ('Aout(idUVwc)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUVwc,:))
            CASE ('Aout(idUbot)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbot,:))
            CASE ('Aout(idVbot)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbot,:))
            CASE ('Aout(idUbur)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUbur,:))
            CASE ('Aout(idVbvr)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVbvr,:))
# endif
# ifdef WEC_MELLOR
            CASE ('Aout(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW2xx,:))
            CASE ('Aout(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW2xy,:))
            CASE ('Aout(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW2yy,:))
# endif
# ifdef WEC
            CASE ('Aout(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU2rs,:))
            CASE ('Aout(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV2rs,:))
            CASE ('Aout(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU2Sd,:))
            CASE ('Aout(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV2Sd,:))
# endif
# ifdef WEC_MELLOR
#  ifdef SOLVE3D
            CASE ('Aout(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3xx,:))
            CASE ('Aout(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3xy,:))
            CASE ('Aout(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3yy,:))
            CASE ('Aout(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3zx,:))
            CASE ('Aout(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3zy,:))
#  endif
# endif
# ifdef WEC
#  ifdef SOLVE3D
            CASE ('Aout(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU3rs,:))
            CASE ('Aout(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV3rs,:))
            CASE ('Aout(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU3Sd,:))
            CASE ('Aout(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV3Sd,:))
            CASE ('Aout(idW3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3Sd,:))
            CASE ('Aout(idW3St)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idW3St,:))
#  endif
            CASE ('Aout(idWztw)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWztw,:))
            CASE ('Aout(idWqsp)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWqsp,:))
            CASE ('Aout(idWbeh)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWbeh,:))
# endif
# ifdef WAVES_HEIGHT
            CASE ('Aout(idWamp)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWamp,:))
# endif
# ifdef WAVES_LENGTH
            CASE ('Aout(idWlen)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWlen,:))
# endif
# ifdef WAVES_LENGTHP
            CASE ('Aout(idWlep)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWlep,:))
# endif
# ifdef WAVES_DIR
            CASE ('Aout(idWdir)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdir,:))
# endif
# ifdef WAVES_TOP_PERIOD
            CASE ('Aout(idWptp)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWptp,:))
# endif
# ifdef WAVES_BOT_PERIOD
            CASE ('Aout(idWpbt)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWpbt,:))
# endif
# ifdef WAVES_OCEAN
            CASE ('Aout(idWorb)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWorb,:))
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
            CASE ('Aout(idWdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdif,:))
# endif
# if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
     defined WDISS_THORGUZA || defined WDISS_CHURTHOR
            CASE ('Aout(idWdib)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdib,:))
            CASE ('Aout(idWdiw)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdiw,:))
# endif
# ifdef ROLLER_SVENDSEN
            CASE ('Aout(idWbrk)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWbrk,:))
# endif
# ifdef WEC_ROLLER
            CASE ('Aout(idWdis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWdis,:))
# endif
# ifdef ROLLER_RENIERS
            CASE ('Aout(idWrol)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWrol,:))
# endif
# ifdef WAVES_DSPR
            CASE ('Aout(idWvds)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWvds,:))
            CASE ('Aout(idWvqp)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idWvqp,:))
# endif
# ifdef UV_KIRBY
            CASE ('Aout(idUwav)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUwav,:))
            CASE ('Aout(idVwav)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVwav,:))
# endif
# ifdef SOLVE3D
#  if defined RUNOFF
            CASE ('Aout(idRunoff)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idRunoff,:))
#  endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
  || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
            CASE ('Aout(idPair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idPair,:))
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
            CASE ('Aout(idUair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUair,:))
            CASE ('Aout(idVair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVair,:))
            CASE ('Aout(idUairE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUairE,:))
            CASE ('Aout(idVairN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVairN,:))
#  endif
# endif
            CASE ('Aout(idTsur)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idLhea,:))
            CASE ('Aout(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idShea,:))
            CASE ('Aout(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idLrad,:))
            CASE ('Aout(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSrad,:))
            CASE ('Aout(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idevap,:))
            CASE ('Aout(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idrain,:))
            CASE ('Aout(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idDano,:))
            CASE ('Aout(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVvis,:))
            CASE ('Aout(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idTdif,:))
            CASE ('Aout(idSdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idSdif,:))
            CASE ('Aout(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHsbl,:))
            CASE ('Aout(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHbbl,:))
            CASE ('Aout(id2dRV)')
              IF (id2dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id2dRV,:))
            CASE ('Aout(id3dRV)')
              IF (id3dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id3dRV,:))
            CASE ('Aout(id2dPV)')
              IF (id2dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id2dPV,:))
            CASE ('Aout(id3dPV)')
              IF (id3dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(id3dPV,:))
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
            CASE ('Aout(idFsuD)')
              IF (idFsuD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idFsuD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idFsuD,:))
            CASE ('Aout(idu2dD)')
              IF (idu2dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu2dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idu2dD,:))
            CASE ('Aout(idv2dD)')
              IF (idv2dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv2dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idv2dD,:))
#  ifdef SOLVE3D
            CASE ('Aout(idu3dD)')
              IF (idu3dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu3dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idu3dD,:))
            CASE ('Aout(idv3dD)')
              IF (idv3dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv3dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idv3dD,:))
            CASE ('Aout(idTrcD)')
              IF (MAXVAL(idTrcD).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTrcD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTrcD(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
#  endif
# endif
# ifdef SOLVE3D
            CASE ('Aout(idHUav)')
              IF (idHUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHUav,:))
            CASE ('Aout(idHVav)')
              IF (idHVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idHVav,:))
            CASE ('Aout(idUUav)')
              IF (idUUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUUav,:))
            CASE ('Aout(idUVav)')
              IF (idUVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idUVav,:))
            CASE ('Aout(idVVav)')
              IF (idVVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idVVav,:))
# endif
            CASE ('Aout(idU2av)')
              IF (idU2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idU2av,:))
            CASE ('Aout(idV2av)')
              IF (idV2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idV2av,:))
            CASE ('Aout(idZZav)')
              IF (idZZav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idZZav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout(idZZav,:))
# ifdef SOLVE3D
            CASE ('Aout(idTTav)')
              IF (MAXVAL(idTTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT+NPT
                  i=idTTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idUTav)')
              IF (MAXVAL(idUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT+NPT
                  i=idUTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idVTav)')
              IF (MAXVAL(idVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT+NPT
                  i=idVTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHUTav)')
              IF (MAXVAL(iHUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT+NPT
                  i=iHUTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iHVTav)')
              IF (MAXVAL(iHVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT+NPT
                  i=iHVTav(itrc)
                  Aout(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
# endif
# if defined SOLVE3D && defined T_PASSIVE
            CASE ('Aout(inert)')
              Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=idTvar(inert(i))
                  Aout(itrc,ng)=Linert(i,ng)
                END DO
              END DO
# endif
#endif
#if defined AVERAGES2
            CASE ('Aout2(idUvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUvel,:))
            CASE ('Aout2(idVvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVvel,:))
            CASE ('Aout2(idWvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idWvel,:))
            CASE ('Aout2(idOvel)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idOvel,:))
            CASE ('Aout2(idUbar)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUbar,:))
            CASE ('Aout2(idVbar)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVbar,:))
            CASE ('Aout2(idFsur)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idFsur,:))
            CASE ('Aout2(idu2dE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idu2dE,:))
            CASE ('Aout2(idv2dN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idv2dN,:))
            CASE ('Aout2(idu3dE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idu3dE,:))
            CASE ('Aout2(idv3dN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idv3dN,:))
            CASE ('Aout2(idTvar)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTvar(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout2(idUsms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUsms,:))
            CASE ('Aout2(idVsms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVsms,:))
            CASE ('Aout2(idUbms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUbms,:))
            CASE ('Aout2(idVbms)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVbms,:))
            CASE ('Aout2(idW2xx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW2xx,:))
            CASE ('Aout2(idW2xy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW2xy,:))
            CASE ('Aout2(idW2yy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW2yy,:))
            CASE ('Aout2(idU2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idU2rs,:))
            CASE ('Aout2(idV2rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idV2rs,:))
            CASE ('Aout2(idU2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idU2Sd,:))
            CASE ('Aout2(idV2Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idV2Sd,:))
            CASE ('Aout2(idW3xx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW3xx,:))
            CASE ('Aout2(idW3xy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW3xy,:))
            CASE ('Aout2(idW3yy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW3yy,:))
            CASE ('Aout2(idW3zx)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW3zx,:))
            CASE ('Aout2(idW3zy)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idW3zy,:))
            CASE ('Aout2(idU3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idU3rs,:))
            CASE ('Aout2(idV3rs)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idV3rs,:))
            CASE ('Aout2(idU3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idU3Sd,:))
            CASE ('Aout2(idV3Sd)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idV3Sd,:))
# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
            CASE ('Aout2(idPair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idPair,:))
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
            CASE ('Aout2(idUair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUair,:))
            CASE ('Aout2(idVair)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVair,:))
            CASE ('Aout2(idUairE)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUairE,:))
            CASE ('Aout2(idVairN)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVairN,:))
#  endif
# endif
            CASE ('Aout2(idTsur)')
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTsur(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout2(idLhea)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idLhea,:))
            CASE ('Aout2(idShea)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idShea,:))
            CASE ('Aout2(idLrad)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idLrad,:))
            CASE ('Aout2(idSrad)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idSrad,:))
            CASE ('Aout2(idevap)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idevap,:))
            CASE ('Aout2(idrain)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idrain,:))
            CASE ('Aout2(idDano)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idDano,:))
            CASE ('Aout2(idVvis)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVvis,:))
            CASE ('Aout2(idTdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idTdif,:))
            CASE ('Aout2(idSdif)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idSdif,:))
            CASE ('Aout2(idHsbl)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHsbl,:))
            CASE ('Aout2(idHbbl)')
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHbbl,:))
            CASE ('Aout2(id2dRV)')
              IF (id2dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(id2dRV,:))
            CASE ('Aout2(id3dRV)')
              IF (id3dRV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dRV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(id3dRV,:))
            CASE ('Aout2(id2dPV)')
              IF (id2dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id2dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(id2dPV,:))
            CASE ('Aout2(id3dPV)')
              IF (id3dPV.eq.0) THEN
                IF (Master) WRITE (out,280) 'id3dPV'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(id3dPV,:))
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
            CASE ('Aout2(idFsuD)')
              IF (idFsuD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idFsuD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idFsuD,:))
            CASE ('Aout2(idu2dD)')
              IF (idu2dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu2dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idu2dD,:))
            CASE ('Aout2(idv2dD)')
              IF (idv2dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv2dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idv2dD,:))
#  ifdef SOLVE3D
            CASE ('Aout2(idu3dD)')
              IF (idu3dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idu3dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idu3dD,:))
            CASE ('Aout2(idv3dD)')
              IF (idv3dD.eq.0) THEN
                IF (Master) WRITE (out,280) 'idv3dD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idv3dD,:))
            CASE ('Aout2(idTrcD)')
              IF (MAXVAL(idTrcD).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTrcD'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTrcD(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
#  endif
# endif
# ifdef SOLVE3D
            CASE ('Aout2(idHUav)')
              IF (idHUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHUav,:))
            CASE ('Aout2(idHVav)')
              IF (idHVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idHVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idHVav,:))
            CASE ('Aout2(idUUav)')
              IF (idUUav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUUav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUUav,:))
            CASE ('Aout2(idUVav)')
              IF (idUVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idUVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idUVav,:))
            CASE ('Aout2(idVVav)')
              IF (idVVav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idVVav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idVVav,:))
# endif
            CASE ('Aout2(idU2av)')
              IF (idU2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idU2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idU2av,:))
            CASE ('Aout2(idV2av)')
              IF (idV2av.eq.0) THEN
                IF (Master) WRITE (out,280) 'idV2av'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idV2av,:))
            CASE ('Aout2(idZZav)')
              IF (idZZav.eq.0) THEN
                IF (Master) WRITE (out,280) 'idZZav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Aout2(idZZav,:))
# ifdef SOLVE3D
            CASE ('Aout2(idTTav)')
              IF (MAXVAL(idTTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idTTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idTTav(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout2(idUTav)')
              IF (MAXVAL(idUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idUTav(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout2(idVTav)')
              IF (MAXVAL(idVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'idVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=idVTav(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout2(iHUTav)')
              IF (MAXVAL(iHUTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHUTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=iHUTav(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
            CASE ('Aout2(iHVTav)')
              IF (MAXVAL(iHVTav).eq.0) THEN
                IF (Master) WRITE (out,280) 'iHVTav'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  i=iHVTav(itrc)
                  Aout2(i,ng)=Ltracer(itrc,ng)
                END DO
              END DO
# endif
# if defined SOLVE3D && defined T_PASSIVE
            CASE ('Aout2(inert)')
              Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
              DO ng=1,Ngrids
                DO i=1,NPT
                  itrc=idTvar(inert(i))
                  Aout2(itrc,ng)=Linert(i,ng)
                END DO
              END DO
# endif
#endif
#ifdef DIAGNOSTICS_UV
            CASE ('Dout(M2rate)')
              IF (M2rate.le.0) THEN
                IF (Master) WRITE (out,280) 'M2rate'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2rate),ng)=Ldiag(ng)
                Dout(idDv2d(M2rate),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2pgrd)')
              IF (M2pgrd.le.0) THEN
                IF (Master) WRITE (out,280) 'M2pgrd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2pgrd),ng)=Ldiag(ng)
                Dout(idDv2d(M2pgrd),ng)=Ldiag(ng)
              END DO
# ifdef UV_COR
            CASE ('Dout(M2fcor)')
              IF (M2fcor.le.0) THEN
                IF (Master) WRITE (out,280) 'M2fcor'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2fcor),ng)=Ldiag(ng)
                Dout(idDv2d(M2fcor),ng)=Ldiag(ng)
              END DO
# endif
# ifdef UV_ADV
            CASE ('Dout(M2hadv)')
              IF (M2hadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M2hadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2hadv),ng)=Ldiag(ng)
                Dout(idDv2d(M2hadv),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2xadv)')
              IF (M2xadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M2xadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2xadv),ng)=Ldiag(ng)
                Dout(idDv2d(M2xadv),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2yadv)')
              IF (M2yadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M2yadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2yadv),ng)=Ldiag(ng)
                Dout(idDv2d(M2yadv),ng)=Ldiag(ng)
              END DO
# endif
# ifdef WEC_MELLOR
            CASE ('Dout(M2hrad)')
              IF (M2hrad.le.0) THEN
                IF (Master) WRITE (out,280) 'M2hrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2hrad),ng)=Ldiag(ng)
                Dout(idDv2d(M2hrad),ng)=Ldiag(ng)
              END DO
# endif
# ifdef WEC_VF
            CASE ('Dout(M2hjvf)')
              IF (M2hjvf.le.0) THEN
                IF (Master) WRITE (out,280) 'M2hjvf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2hjvf),ng)=Ldiag(ng)
                Dout(idDv2d(M2hjvf),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2kvrf)')
              IF (M2kvrf.le.0) THEN
                IF (Master) WRITE (out,280) 'M2kvrf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2kvrf),ng)=Ldiag(ng)
                Dout(idDv2d(M2kvrf),ng)=Ldiag(ng)
              END DO
#  ifdef UV_COR
            CASE ('Dout(M2fsco)')
              IF (M2fsco.le.0) THEN
                IF (Master) WRITE (out,280) 'M2fsco'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2fsco),ng)=Ldiag(ng)
                Dout(idDv2d(M2fsco),ng)=Ldiag(ng)
              END DO
#  endif
#  ifdef BOTTOM_STREAMING
            CASE ('Dout(M2bstm)')
              IF (M2bstm.le.0) THEN
                IF (Master) WRITE (out,280) 'M2bstm'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2bstm),ng)=Ldiag(ng)
                Dout(idDv2d(M2bstm),ng)=Ldiag(ng)
              END DO
#  endif
#  ifdef SURFACE_STREAMING
            CASE ('Dout(M2sstm)')
              IF (M2sstm.le.0) THEN
                IF (Master) WRITE (out,280) 'M2sstm'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2sstm),ng)=Ldiag(ng)
                Dout(idDv2d(M2sstm),ng)=Ldiag(ng)
              END DO
#  endif
            CASE ('Dout(M2wrol)')
              IF (M2wrol.le.0) THEN
                IF (Master) WRITE (out,280) 'M2wrol'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2wrol),ng)=Ldiag(ng)
                Dout(idDv2d(M2wrol),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2wbrk)')
              IF (M2wbrk.le.0) THEN
                IF (Master) WRITE (out,280) 'M2wbrk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2wbrk),ng)=Ldiag(ng)
                Dout(idDv2d(M2wbrk),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2zeta)')
              IF (M2zeta.le.0) THEN
                IF (Master) WRITE (out,280) 'M2zeta'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2zeta),ng)=Ldiag(ng)
                Dout(idDv2d(M2zeta),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2zetw)')
              IF (M2zetw.le.0) THEN
                IF (Master) WRITE (out,280) 'M2zetw'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2zetw),ng)=Ldiag(ng)
                Dout(idDv2d(M2zetw),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2zqsp)')
              IF (M2zqsp.le.0) THEN
                IF (Master) WRITE (out,280) 'M2zqsp'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2zqsp),ng)=Ldiag(ng)
                Dout(idDv2d(M2zqsp),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2zbeh)')
              IF (M2zbeh.le.0) THEN
                IF (Master) WRITE (out,280) 'M2zbeh'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2zbeh),ng)=Ldiag(ng)
                Dout(idDv2d(M2zbeh),ng)=Ldiag(ng)
              END DO
# endif
# if defined UV_VIS2 || defined UV_VIS4
            CASE ('Dout(M2hvis)')
              IF (M2hvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M2hvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2hvis),ng)=Ldiag(ng)
                Dout(idDv2d(M2hvis),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2xvis)')
              IF (M2xvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M2xvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2xvis),ng)=Ldiag(ng)
                Dout(idDv2d(M2xvis),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2yvis)')
              IF (M2yvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M2yvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2yvis),ng)=Ldiag(ng)
                Dout(idDv2d(M2yvis),ng)=Ldiag(ng)
              END DO
# endif
            CASE ('Dout(M2sstr)')
              IF (M2sstr.le.0) THEN
                IF (Master) WRITE (out,280) 'M2sstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2sstr),ng)=Ldiag(ng)
                Dout(idDv2d(M2sstr),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2bstr)')
              IF (M2bstr.le.0) THEN
                IF (Master) WRITE (out,280) 'M2bstr'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2bstr),ng)=Ldiag(ng)
                Dout(idDv2d(M2bstr),ng)=Ldiag(ng)
              END DO
# ifdef SOLVE3D
            CASE ('Dout(M3rate)')
              IF (M3rate.le.0) THEN
                IF (Master) WRITE (out,280) 'M3rate'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3rate),ng)=Ldiag(ng)
                Dout(idDv3d(M3rate),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3pgrd)')
              IF (M3pgrd.le.0) THEN
                IF (Master) WRITE (out,280) 'M3pgrd'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3pgrd),ng)=Ldiag(ng)
                Dout(idDv3d(M3pgrd),ng)=Ldiag(ng)
              END DO
#  ifdef UV_COR
            CASE ('Dout(M3fcor)')
              IF (M3fcor.le.0) THEN
                IF (Master) WRITE (out,280) 'M3fcor'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3fcor),ng)=Ldiag(ng)
                Dout(idDv3d(M3fcor),ng)=Ldiag(ng)
              END DO
#  endif
#  ifdef UV_ADV
            CASE ('Dout(M3hadv)')
              IF (M3hadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M3hadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3hadv),ng)=Ldiag(ng)
                Dout(idDv3d(M3hadv),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3xadv)')
              IF (M3xadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M3xadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3xadv),ng)=Ldiag(ng)
                Dout(idDv3d(M3xadv),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3yadv)')
              IF (M3yadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M3yadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3yadv),ng)=Ldiag(ng)
                Dout(idDv3d(M3yadv),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3vadv)')
              IF (M3vadv.le.0) THEN
                IF (Master) WRITE (out,280) 'M3vadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3vadv),ng)=Ldiag(ng)
                Dout(idDv3d(M3vadv),ng)=Ldiag(ng)
              END DO
#  endif
#  ifdef WEC_MELLOR
            CASE ('Dout(M3hrad)')
              IF (M3hrad.le.0) THEN
                IF (Master) WRITE (out,280) 'M3hrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3hrad),ng)=Ldiag(ng)
                Dout(idDv3d(M3hrad),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3vrad)')
              IF (M3vrad.le.0) THEN
                IF (Master) WRITE (out,280) 'M3vrad'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3vrad),ng)=Ldiag(ng)
                Dout(idDv3d(M3vrad),ng)=Ldiag(ng)
              END DO
#  endif
#  ifdef WEC_VF
            CASE ('Dout(M3vjvf)')
              IF (M3vjvf.le.0) THEN
                IF (Master) WRITE (out,280) 'M3vjvf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3vjvf),ng)=Ldiag(ng)
                Dout(idDv3d(M3vjvf),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3hjvf)')
              IF (M3hjvf.le.0) THEN
                IF (Master) WRITE (out,280) 'M3hjvf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3hjvf),ng)=Ldiag(ng)
                Dout(idDv3d(M3hjvf),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3kvrf)')
              IF (M3kvrf.le.0) THEN
                IF (Master) WRITE (out,280) 'M3kvrf'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3kvrf),ng)=Ldiag(ng)
                Dout(idDv3d(M3kvrf),ng)=Ldiag(ng)
              END DO
#   ifdef UV_COR
            CASE ('Dout(M3fsco)')
              IF (M3fsco.le.0) THEN
                IF (Master) WRITE (out,280) 'M3fsco'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3fsco),ng)=Ldiag(ng)
                Dout(idDv3d(M3fsco),ng)=Ldiag(ng)
              END DO
#   endif
#   ifdef BOTTOM_STREAMING
            CASE ('Dout(M3bstm)')
              IF (M3bstm.le.0) THEN
                IF (Master) WRITE (out,280) 'M3bstm'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3bstm),ng)=Ldiag(ng)
                Dout(idDv3d(M3bstm),ng)=Ldiag(ng)
              END DO
#   endif
#   ifdef SURFACE_STREAMING
            CASE ('Dout(M3sstm)')
              IF (M3sstm.le.0) THEN
                IF (Master) WRITE (out,280) 'M3sstm'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3sstm),ng)=Ldiag(ng)
                Dout(idDv3d(M3sstm),ng)=Ldiag(ng)
              END DO
#   endif
            CASE ('Dout(M3wrol)')
              IF (M3wrol.le.0) THEN
                IF (Master) WRITE (out,280) 'M3wrol'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3wrol),ng)=Ldiag(ng)
                Dout(idDv3d(M3wrol),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3wbrk)')
              IF (M3wbrk.le.0) THEN
                IF (Master) WRITE (out,280) 'M3wbrk'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3wbrk),ng)=Ldiag(ng)
                Dout(idDv3d(M3wbrk),ng)=Ldiag(ng)
              END DO
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
            CASE ('Dout(M3hvis)')
              IF (M3hvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M3hvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3hvis),ng)=Ldiag(ng)
                Dout(idDv3d(M3hvis),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3xvis)')
              IF (M3xvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M3xvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3xvis),ng)=Ldiag(ng)
                Dout(idDv3d(M3xvis),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M3yvis)')
              IF (M3yvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M3yvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3yvis),ng)=Ldiag(ng)
                Dout(idDv3d(M3yvis),ng)=Ldiag(ng)
              END DO
#  endif
            CASE ('Dout(M3vvis)')
              IF (M3vvis.le.0) THEN
                IF (Master) WRITE (out,280) 'M3vvis'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3vvis),ng)=Ldiag(ng)
                Dout(idDv3d(M3vvis),ng)=Ldiag(ng)
              END DO
#  if defined VEGETATION && defined VEG_DRAG
            CASE ('Dout(M3fveg)')
              IF (M3fveg.le.0) THEN
                IF (Master) WRITE (out,280) 'M3fveg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu3d(M3fveg),ng)=Ldiag(ng)
                Dout(idDv3d(M3fveg),ng)=Ldiag(ng)
              END DO
            CASE ('Dout(M2fveg)')
              IF (M2fveg.le.0) THEN
                IF (Master) WRITE (out,280) 'M2fveg'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Ldiag)
              DO ng=1,Ngrids
                Dout(idDu2d(M2fveg),ng)=Ldiag(ng)
                Dout(idDv2d(M2fveg),ng)=Ldiag(ng)
              END DO
#  endif
# endif
#endif
#if defined DIAGNOSTICS_TS && defined SOLVE3D
            CASE ('Dout(iTrate)')
              IF (iTrate.le.0) THEN
                IF (Master) WRITE (out,280) 'iTrate'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTrate),ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTrate),ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('Dout(iThadv)')
              IF (iThadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iThadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iThadv),ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iThadv),ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('Dout(iTxadv)')
              IF (iTxadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iTxadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTxadv),ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('Dout(iTyadv)')
              IF (iTyadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iTyadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTyadv),ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
            CASE ('Dout(iTvadv)')
              IF (iTvadv.le.0) THEN
                IF (Master) WRITE (out,280) 'iTvadv'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTvadv),ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
# if defined TS_DIF2 || defined TS_DIF4
            CASE ('Dout(iThdif)')
              IF (iThdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iThdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iThdif),ng)=Ltracer(itrc,ng)
                END DO
#  ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iThdif),ng)=Ltracer(NAT+i,ng)
                END DO
#  endif
              END DO
            CASE ('Dout(iTxdif)')
              IF (iTxdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTxdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltracer(itrc,ng)
                END DO
#  ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTxdif),ng)=Ltracer(NAT+i,ng)
                END DO
#  endif
              END DO
            CASE ('Dout(iTydif)')
              IF (iTydif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTydif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTydif),ng)=Ltracer(itrc,ng)
                END DO
#  ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTydif),ng)=Ltracer(NAT+i,ng)
                END DO
#  endif
              END DO
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
            CASE ('Dout(iTsdif)')
              IF (iTsdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTsdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltracer(itrc,ng)
                END DO
#   ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTsdif),ng)=Ltracer(NAT+i,ng)
                END DO
#   endif
              END DO
#  endif
# endif
            CASE ('Dout(iTvdif)')
              IF (iTvdif.le.0) THEN
                IF (Master) WRITE (out,280) 'iTvdif'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
              DO ng=1,Ngrids
                DO itrc=1,NAT
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltracer(itrc,ng)
                END DO
# ifdef T_PASSIVE
                DO i=1,NPT
                  itrc=inert(i)
                  Dout(idDtrc(itrc,iTvdif),ng)=Ltracer(NAT+i,ng)
                END DO
# endif
              END DO
#endif
            CASE ('NUSER')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              Nuser=Ivalue(1)
            CASE ('USER')
              Npts=load_r(Nval, Rval, MAX(1,Nuser), user)
            CASE ('NC_SHUFFLE')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              shuffle=Ivalue(1)
            CASE ('NC_DEFLATE')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              deflate=Ivalue(1)
            CASE ('NC_DLEVEL')
              Npts=load_i(Nval, Rval, 1, Ivalue)
              deflate_level=Ivalue(1)
            CASE ('DAINAME')
              label='DAI - Data Assimilation Initial/Restart fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, DAI)
            CASE ('GSTNAME')
              label='GST - generalized stability theory analysis'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, GST)
            CASE ('RSTNAME')
              label='RST - restart fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, RST)
            CASE ('HISNAME')
              label='HIS - nonlinear model history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, HIS)
            CASE ('HIS2NAME')
              label='HIS2 - nonlinear model surface fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, HIS2)
            CASE ('TLMNAME')
              label='TLM - tangent linear model history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, TLM)
            CASE ('TLFNAME')
              label='TLF - tangent linear model impulse forcing'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, TLF)
            CASE ('ADJNAME')
              label='ADM - adjoint model history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, ADM)
            CASE ('AVGNAME')
              label='AVG - time-averaged history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, AVG)
            CASE ('AVG2NAME')
              label='AVG2 - time-averaged history fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, AVG2)
            CASE ('DIANAME')
              label='DIA - time-averaged diagnostics fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, DIA)
            CASE ('STANAME')
              label='STA - stations time-series'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, STA)
            CASE ('FLTNAME')
              label='FLT - Lagragian particles trajectories'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, FLT)
            CASE ('GRDNAME')
              label='GRD - application grid'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, GRD)
            CASE ('ININAME')
              label='INI - nonlinear model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, INI)
            CASE ('IRPNAME')
              label='IRP - representer model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, IRP)
            CASE ('ITLNAME')
              label='ITL - tangent linear model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, ITL)
            CASE ('IADNAME')
              label='IAD - adjoint model initial conditions'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, IAD)
            CASE ('FWDNAME')
              label='FWD - basic state forward trajectory'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, FWD)
            CASE ('ADSNAME')
              label='ADS - adjoint sensitivity functional'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, ADS)
            CASE ('NGCNAME')
              DO i=1,LEN(NGCname)
                NGCname(i:i)=blank
              END DO
              NGCname=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('NBCFILES')
              Npts=load_i(Nval, Rval, Ngrids, nBCfiles)
              DO ng=1,Ngrids
                IF (nBCfiles(ng).le.0) THEN 
                  IF (Master) WRITE (out,260) 'NBCFILES', nBCfiles(ng), &
     &                            'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nBCfiles)
              allocate ( BRY(max_Ffiles,Ngrids) )
              allocate ( BRYids(max_Ffiles,Ngrids) )
              allocate ( NBCcount(max_Ffiles,Ngrids) )
              BRYids(1:max_Ffiles,1:Ngrids)=-1
              NBCcount(1:max_Ffiles,1:Ngrids)=0
            CASE ('BRYNAME')
              label='BRY - lateral open boundary conditions'
              Npts=load_s2d(Nval, Cval, Cdim, line, label, ibcfile,     &
     &                      igrid, nBCfiles, NBCcount, nBCfiles, BRY)
            CASE ('NCLMFILES')
              Npts=load_i(Nval, Rval, Ngrids, nCLMfiles)
              DO ng=1,Ngrids
                IF (nCLMfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NCLMFILES',              &
     &               nCLMfiles(ng),                                     &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nCLMfiles)
              allocate ( CLM(max_Ffiles,Ngrids) )
              allocate ( CLMids(max_Ffiles,Ngrids) )
              allocate ( NCLMcount(max_Ffiles,Ngrids) )
              CLMids(1:max_Ffiles,1:Ngrids)=-1
              NCLMcount(1:max_Ffiles,1:Ngrids)=0
            CASE ('CLMNAME')
              label='CLM - climatology fields'
              Npts=load_s2d(Nval, Cval, Cdim, line, label, iclmfile,    &
     &                      igrid, nCLMfiles, NCLMcount, nCLMfiles, CLM)
            CASE ('NUDNAME')
              label='NUD - nudging coefficients'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, NUD)
            CASE ('SSFNAME')
              label='SSF - Sources/Sinks forcing fields'
              Npts=load_s1d(Nval, Cval, Cdim, line, label, igrid,       &
     &                      Nfiles, SSF)
            CASE ('NFFILES')
              Npts=load_i(Nval, Rval, Ngrids, nFfiles)
              DO ng=1,Ngrids
                IF (nFfiles(ng).le.0) THEN
                  IF (Master) WRITE (out,260) 'NFFILES', nFfiles(ng),   &
     &              'Must be equal or greater than one.'
                  exit_flag=4
                  RETURN
                END IF
              END DO
              max_Ffiles=MAXVAL(nFfiles)
              allocate ( FRC(max_Ffiles,Ngrids) )
              allocate ( FRCids(max_Ffiles,Ngrids) )
              allocate ( Ncount(max_Ffiles,Ngrids) )
              FRCids(1:max_Ffiles,1:Ngrids)=-1
              Ncount(1:max_Ffiles,1:Ngrids)=0
            CASE ('FRCNAME')
              label='FRC - forcing fields'
              Npts=load_s2d(Nval, Cval, Cdim, line, label, ifile,       &
     &                      igrid, nFfiles, Ncount, max_Ffiles, FRC)
            CASE ('APARNAM')
              DO i=1,LEN(aparnam)
                aparnam(i:i)=blank
              END DO
              aparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('SPOSNAM')
              DO i=1,LEN(sposnam)
                sposnam(i:i)=blank
              END DO
              sposnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('FPOSNAM')
              DO i=1,LEN(fposnam)
                fposnam(i:i)=blank
              END DO
              fposnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('BPARNAM')
              DO i=1,LEN(bparnam)
                bparnam(i:i)=blank
              END DO
              bparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('IPARNAM')
              DO i=1,LEN(iparnam)
                iparnam(i:i)=blank
              END DO
              iparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('SPARNAM')
              DO i=1,LEN(sparnam)
                sparnam(i:i)=blank
              END DO
              sparnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('VEGNAM')
              DO i=1,LEN(vegnam)
                vegnam(i:i)=blank
              END DO
              vegnam=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('USRNAME')
              DO i=1,LEN(USRname)
                USRname(i:i)=blank
              END DO
              USRname=TRIM(ADJUSTL(Cval(Nval)))
#ifdef INWAVE_MODEL
            CASE ('IWININAME')
              DO i=1,LEN(IWININAME(1))
                IWININAME(1)(i:i)=blank
              END DO
              IWININAME(1)=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('IWBNDNAME')
              DO i=1,LEN(IWBNDNAME(1))
                IWBNDNAME(1)(i:i)=blank
              END DO
              IWBNDNAME(1)=TRIM(ADJUSTL(Cval(Nval)))
            CASE ('IWSWNNAME')
              DO i=1,LEN(IWSWNNAME(1))
                IWSWNNAME(1)(i:i)=blank
              END DO
              IWSWNNAME(1)=TRIM(ADJUSTL(Cval(Nval)))
#endif
          END SELECT
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Process input parameters.
!-----------------------------------------------------------------------
!
!  Check if nesting parameters "NestLayers", "GridsInLayer", and
!  "GridNumber" have been assigned.  The code below is necessary
!  for compatability with old "ocean.in" input scripts.
!
#ifndef NESTING
      IF (.not.got_NestLayers) THEN
        NestLayers=1
        IF (.not.allocated(GridsInLayer)) THEN
          allocate ( GridsInLayer(NestLayers) )
        END IF
        IF (.not.allocated(GridNumber)) THEN
          allocate ( GridNumber(Ngrids,NestLayers) )
        END IF
      END IF
      GridsInLayer=1              ! In case that users set illegal
      GridNumber=1                ! values in non-nesting applications
#else
      IF (.not.got_NestLayers) THEN
        IF (Master) THEN
          WRITE (out,320) 'NestLayers',                                 &
     &            'Add "NestLayers" keyword after "Ngrids".'
          exit_flag=5
          RETURN
        END IF
      END IF
      IF (.not.allocated(GridsInLayer)) THEN
        IF (Master) THEN
          WRITE (out,320) 'GridsInLayer',                               &
     &            'Add "GridsInLayer" keyword after "NestLayers".'
          exit_flag=5
          RETURN
        END IF
      END IF
#endif
!
!  Check if both point sources methodologies are activated.  Only one
!  method is allowed for a particular grid.  Otherwise, the point
!  source will be applies twice.
!
      DO ng=1,Ngrids
       IF (LuvSrc(ng).and.LwSrc(ng)) THEN
         IF (Master) THEN
           WRITE (out,260) 'LuvSrc', LuvSrc(ng),                        &
     &           'Because LwSrc  is also ON; only one method is legal.'
           WRITE (out,260) 'LwSrc', LwSrc(ng),                          &
     &           'Because LuvSrc is also ON; only one method is legal.'
           exit_flag=4
           RETURN
         END IF
       END IF
     END DO

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
!
!  Allocate and initialize BLK I/O structure with the same values as
!  HIS.
!
      DO ng=1,Ngrids
        Nfiles(ng)=MAX(1,HIS(ng)%Nfiles)
        allocate ( BLK(ng)%Nrec(Nfiles(ng)) )
        allocate ( BLK(ng)%time_min(Nfiles(ng)) )
        allocate ( BLK(ng)%time_max(Nfiles(ng)) )
        allocate ( BLK(ng)%Vid(NV) )
        allocate ( BLK(ng)%Tid(MT) )
        allocate ( BLK(ng)%files(Nfiles(ng)) )
        BLK(ng)%Nfiles=Nfiles(ng)
        BLK(ng)%Fcount=1
        BLK(ng)%Rindex=0
        BLK(ng)%ncid=-1
        BLK(ng)%Vid(1:NV)=-1
        BLK(ng)%Tid(1:MT)=-1
        BLK(ng)%Nrec(1:Nfiles(ng))=0
        BLK(ng)%time_min(1:Nfiles(ng))=0.0_r8
        BLK(ng)%time_max(1:Nfiles(ng))=0.0_r8
        BLK(ng)%label='BLK - nonlinear model bulk fluxes'
      END DO
!
!  Make sure that logical output switches are activated for wind
!  stress and surface active tracers fluxes when using fluxes from
!  the nonlinear model via (bulk_flux).
!
      DO ng=1,Ngrids
        Hout(idUsms,ng)=.TRUE.
        Hout(idVsms,ng)=.TRUE.
# ifdef SOLVE3D
        Hout(idTsur(itemp),ng)=.TRUE.
#  ifdef EMINUSP
        Hout(idEmPf,ng)=.TRUE.
#  endif
# endif
      END DO
#endif

#if defined FORCING_SV     || defined SO_SEMI || \
    defined STOCHASTIC_OPT
!
!  Make sure that logical output switches are activated for the
!  Forcing Singular Vectors or Stochastic Optimals variables.
!
      DO ng=1,Ngrids
        IF (SCALARS(ng)%Fstate(isFsur)) Hout(idFsur,ng)=.TRUE.
# ifndef SOLVE3D
        IF (SCALARS(ng)%Fstate(isUbar)) Hout(idUbar,ng)=.TRUE.
        IF (SCALARS(ng)%Fstate(isVbar)) Hout(idVbar,ng)=.TRUE.
# else
        IF (SCALARS(ng)%Fstate(isUvel)) Hout(idUvel,ng)=.TRUE.
        IF (SCALARS(ng)%Fstate(isVvel)) Hout(idVvel,ng)=.TRUE.
        DO itrc=1,NT(ng)
          IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
            Hout(idTvar(itrc),ng)=.TRUE.
          END IF
          IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
            Hout(idTsur(itrc),ng)=.TRUE.
          END IF
        END DO
# endif
        IF (SCALARS(ng)%Fstate(isUstr)) Hout(idUsms,ng)=.TRUE.
        IF (SCALARS(ng)%Fstate(isVstr)) Hout(idVsms,ng)=.TRUE.
      END DO
#endif
!
!  Make sure that both component switches are activated when processing
!  (Eastward,Northward) momentum components at RHO-points.
!
      DO ng=1,Ngrids
        IF (.not.Hout(idu2dE,ng).and.Hout(idv2dN,ng)) THEN
          Hout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv2dN,ng).and.Hout(idu2dE,ng)) THEN
          Hout(idv2dN,ng)=.TRUE.
        END IF
#ifdef SOLVE3D
        IF (.not.Hout(idu3dE,ng).and.Hout(idv3dN,ng)) THEN
          Hout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Hout(idv3dN,ng).and.Hout(idu3dE,ng)) THEN
          Hout(idv3dN,ng)=.TRUE.
        END IF
#endif
#ifdef AVERAGES
        IF (.not.Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          Aout(idu2dE,ng)=.TRUE.
        END IF
        IF (.not.Aout(idv2dN,ng).and.Aout(idu2dE,ng)) THEN
          Aout(idv2dN,ng)=.TRUE.
        END IF
# ifdef SOLVE3D
        IF (.not.Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          Aout(idu3dE,ng)=.TRUE.
        END IF
        IF (.not.Aout(idv3dN,ng).and.Aout(idu3dE,ng)) THEN
          Aout(idv3dN,ng)=.TRUE.
        END IF
# endif
#endif
      END DO
!
!  Set various parameters.
!
      DO ng=1,Ngrids
!
!  Set switch to craeate NetCDF file.
!
#ifdef SOLVE3D
        DO i=1,NV
          IF (Hout(i,ng)) LdefHIS(ng)=.TRUE.
        END DO
#else
        IF (Hout(idFsur,ng).or.Hout(idUbar,ng).or.Hout(idVbar,ng)) THEN
          LdefHIS(ng)=.TRUE.
        END IF
#endif
!
!  Set switch to process climatology file.
!
#ifndef ANA_SSH
        IF (LsshCLM(ng)) CLM_FILE(ng)=.TRUE.
#endif
#ifndef ANA_M2CLIMA
        IF (Lm2CLM(ng)) CLM_FILE(ng)=.TRUE.
#endif
#ifdef SOLVE3D
# ifndef ANA_M3CLIMA
        IF (Lm3CLM(ng)) CLM_FILE(ng)=.TRUE.
# endif
# ifndef ANA_TCLIMA
        IF (ANY(LtracerCLM(:,ng))) CLM_FILE(ng)=.TRUE.
# endif
# if defined CLIMA_TS_MIX && (defined TS_DIF2 || defined TS_DIF4)
        CLM_FILE(ng)=.TRUE.
# endif
#endif

#if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  If weak constraint, disallow recycling of the adjoint model.
!
        LcycleADJ(ng)=.FALSE.
#endif
#if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  If computing time-averaged detided fields, make sure the the regular
!  time-averaged switches are activated for such variables.
!
        IF (.not.Aout(idFsur,ng).and.Aout(idFsuD,ng)) THEN
          Aout(idFsur,ng)=.TRUE.
        END IF
        IF (.not.Aout(idUbar,ng).and.Aout(idu2dD,ng)) THEN
          Aout(idUbar,ng)=.TRUE.
        END IF
        IF (.not.Aout(idVbar,ng).and.Aout(idv2dD,ng)) THEN
          Aout(idVbar,ng)=.TRUE.
        END IF
# ifdef SOLVE3D
        IF (.not.Aout(idUvel,ng).and.Aout(idu3dD,ng)) THEN
          Aout(idUvel,ng)=.TRUE.
        END IF
        IF (.not.Aout(idVvel,ng).and.Aout(idv3dD,ng)) THEN
          Aout(idVvel,ng)=.TRUE.
        END IF
        DO itrc=1,NAT
          IF (.not.Aout(idTvar(itrc),ng).and.                           &
     &        Aout(idTrcD(itrc),ng)) THEN
            Aout(idTvar(itrc),ng)=.TRUE.
          END IF
        END DO
# endif
#endif
#if defined IS4DVAR
!
!  If strong constraint, write only final adjoint solution since only
!  we are estimating initial conditions.
!
        nADJ(ng)=ntimes(ng)
!
!  Insure that restart file is written only at the end.  In sequential
!  data assimilation the restart file is used as the first guess for
!  the next assimilation cycle.
!
        nRST(ng)=ntimes(ng)
#endif
#if defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE
!
!  For red noise stochastic optimals, nADJ must equal the number of
!  trapezoidal intervals.
!
        nADJ(ng)=ntimes(ng)/Nintervals
#endif
#if defined FOUR_DVAR || defined IMPULSE
!
!  Set size of additonal dimension for error covariance normalization
!  and standard deviation factors.
!
# if defined CORRELATION || defined WEAK_CONSTRAINT
        IF (nADJ(ng).lt.ntimes(ng)) THEN
          NSA=2
        ELSE
          NSA=1
        END IF
# else
        NSA=1
# endif
#endif
#ifdef WEAK_CONSTRAINT
!
!  If weak constraint assimilation, set tangent linear number of
!  time-steps between writing of fields the same as the basic
!  state.  Disallow recycling.
!
        nTLM(ng)=nHIS(ng)
        LcycleTLM(ng)=.FALSE.
#endif
#if defined TIME_CONV && defined WEAK_CONSTRAINT
!
!  Set number of full state records needed for time convolutions.
!
        NrecTC(ng)=(ntimes(ng)/nADJ(ng))+1
#endif
#if defined FOUR_DVAR
!
!  If variational data assimilation, disallow creation of multiple
!  output history files.  Otherwise, the forward basic state will
!  not processed correctly.
!
        ndefHIS(ng)=0
        ndefTLM(ng)=0
        ndefADJ(ng)=0
#endif
        IF (((nrrec(ng).eq.0).and.(nAVG(ng).gt.ntimes(ng))).or.         &
     &      (nAVG(ng).eq.0)) THEN
          LdefAVG(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nAVG2(ng).gt.ntimes(ng))).or.        &
     &      (nAVG2(ng).eq.0)) THEN
          LdefAVG2(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nDIA(ng).gt.ntimes(ng))).or.         &
     &      (nDIA(ng).eq.0)) THEN
          LdefDIA(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nFLT(ng).gt.ntimes(ng))).or.         &
     &      (nFLT(ng).eq.0)) THEN
          LdefFLT(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nHIS(ng).gt.ntimes(ng))).or.         &
     &      (nHIS(ng).eq.0)) THEN
          LdefHIS(ng)=.FALSE.
        END IF
        IF (((nrrec(ng).eq.0).and.(nRST(ng).gt.ntimes(ng))).or.         &
     &      (nRST(ng).eq.0)) THEN
          LdefRST(ng)=.FALSE.
        END  IF
        IF (((nrrec(ng).eq.0).and.(nSTA(ng).gt.ntimes(ng))).or.         &
     &      (nSTA(ng).eq.0)) THEN
          LdefSTA(ng)=.FALSE.
        END IF
!
!  Determine switch to process boundary NetCDF file.
!
        ObcData(ng)=.FALSE.
#ifndef ANA_FSOBC
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isFsur,ng)%acquire)
# if defined ADJOINT || defined TANGENT || defined TL_IOMS
        ObcData(ng)=ObcData(ng).or.ANY(ad_LBC(:,isFsur,ng)%acquire)
# endif
#endif
#ifndef ANA_M2OBC
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isUbar,ng)%acquire)        &
     &                         .or.ANY(LBC(:,isVbar,ng)%acquire)
# if defined ADJOINT || defined TANGENT || defined TL_IOMS
        ObcData(ng)=ObcData(ng).or.ANY(ad_LBC(:,isUbar,ng)%acquire)     &
     &                         .or.ANY(ad_LBC(:,isVbar,ng)%acquire)
# endif
#endif
#ifdef SOLVE3D
# ifndef ANA_M3OBC
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isUvel,ng)%acquire)        &
     &                         .or.ANY(LBC(:,isVvel,ng)%acquire)
#  if defined ADJOINT || defined TANGENT || defined TL_IOMS
        ObcData(ng)=ObcData(ng).or.ANY(ad_LBC(:,isUvel,ng)%acquire)     &
     &                         .or.ANY(ad_LBC(:,isVvel,ng)%acquire)
#  endif
# endif
# ifndef ANA_TOBC
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isTvar(:),ng)%acquire)
#  if defined ADJOINT || defined TANGENT || defined TL_IOMS
        ObcData(ng)=ObcData(ng).or.ANY(ad_LBC(:,isTvar(:),ng)%acquire)
#  endif
# endif
# ifdef INWAVE_MODEL
        ObcData(ng)=ObcData(ng).or.ANY(LBC(:,isAC3d,ng)%acquire)
# endif
#endif
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
#ifdef DISTRIBUTE
        WRITE (out,60) TRIM(title), TRIM(my_os), TRIM(my_cpu),          &
     &                 TRIM(my_fort), TRIM(my_fc), TRIM(my_fflags),     &
     &                 TRIM(Iname), TRIM(svn_url), TRIM(svn_rev),       &
     &                 TRIM(Rdir), TRIM(Hdir), TRIM(Hfile), TRIM(Adir)
#else
        WRITE (out,60) TRIM(title), TRIM(my_os), TRIM(my_cpu),          &
     &                 TRIM(my_fort), TRIM(my_fc), TRIM(my_fflags),     &
     &                 TRIM(svn_url), TRIM(svn_rev), TRIM(Rdir),        &
     &                 TRIM(Hdir), TRIM(Hfile), TRIM(Adir)
#endif

        DO ng=1,Ngrids
!
!  Report grid size and domain decomposition.  Check for correct tile
!  decomposition.
!
#ifdef DISTRIBUTE
          WRITE (out,70) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
     &                   NtileI(ng), NtileJ(ng)
          IF ((NtileI(ng)*NtileJ(ng)).ne.numthreads) THEN
            WRITE (out,80) ng
            exit_flag=6
            RETURN
          END IF
#else
          WRITE (out,90) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
     &                   NtileI(ng), NtileJ(ng)
          IF (NtileI(ng)*NtileJ(ng).le.0) THEN
            WRITE (out,100) ng
            exit_flag=6
            RETURN
          END IF
          IF (MOD(NtileI(ng)*NtileJ(ng),numthreads).ne.0) THEN
            WRITE (out,100) ng
            exit_flag=6
            RETURN
          END IF
#endif
!
!  Report physical parameters.
!
          WRITE (out,110) ng
          WRITE (out,120) ntimes(ng), 'ntimes',                         &
     &          'Number of timesteps for 3-D equations.'
          WRITE (out,140) dt(ng), 'dt',                                 &
     &          'Timestep size (s) for 3-D equations.'
          WRITE (out,130) ndtfast(ng), 'ndtfast',                       &
     &          'Number of timesteps for 2-D equations between',        &
     &          'each 3D timestep.'
          WRITE (out,120) ERstr, 'ERstr',                               &
     &          'Starting ensemble/perturbation run number.'
          WRITE (out,120) ERend, 'ERend',                               &
     &          'Ending ensemble/perturbation run number.'
#ifdef FOUR_DVAR
          WRITE (out,120) Nouter, 'Nouter',                             &
     &          'Maximun number of 4DVAR outer loop iterations.'
#endif
#if defined IS4DVAR           || defined IS4DVAR_SENSITIVITY || \
    defined SENSITIVITY_4DVAR || \
    defined TL_W4DPSAS        || defined TL_W4DVAR           || \
    defined W4DPSAS           || defined W4DVAR              || \
    defined HESSIAN_SV        || defined HESSIAN_SO          || \
    defined HESSIAN_FSV
          WRITE (out,120) Ninner, 'Ninner',                             &
     &          'Maximun number of 4DVAR inner loop iterations.'
#endif
#ifdef STOCHASTIC_OPT
          WRITE (out,120) Nintervals, 'Nintervals',                     &
     &          'Number of stochastic optimals timestep intervals.'
#endif
#ifdef PROPAGATOR
          WRITE (out,120) NEV, 'NEV',                                   &
     &          'Number of Lanczos/Arnoldi eigenvalues to compute.'
          WRITE (out,120) NCV, 'NCV',                                   &
     &          'Number of Lanczos/Arnoldi eigenvectors to compute.'
#endif
          WRITE (out,120) nrrec(ng), 'nrrec',                           &
     &          'Number of restart records to read from disk.'
          WRITE (out,170) LcycleRST(ng), 'LcycleRST',                   &
     &          'Switch to recycle time-records in restart file.'
          WRITE (out,130) nRST(ng), 'nRST',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into restart fields.'
          WRITE (out,130) ninfo(ng), 'ninfo',                           &
     &          'Number of timesteps between print of information',     &
     &          'to standard output.'
#ifdef STATIONS
          WRITE (out,130) nSTA(ng), 'nSTA',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'the stations file.'
#endif
#ifdef FLOATS
          WRITE (out,130) nFLT(ng), 'nFLT',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into floats file.'
#endif
          WRITE (out,170) ldefout(ng), 'ldefout',                       &
     &          'Switch to create a new output NetCDF file(s).'
          WRITE (out,130) nHIS(ng), 'nHIS',                             &
     &          'Number of timesteps between the writing fields',       &
     &          'into history file.'
          IF (ndefHIS(ng).gt.0) THEN
            WRITE (out,130) ndefHIS(ng), 'ndefHIS',                     &
     &            'Number of timesteps between creation of new',        &
     &            'history files.'
          END IF
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          WRITE (out,130) ntsAVG(ng), 'ntsAVG',                         &
     &          'Starting timestep for the accumulation of output',     &
     &          'time-averaged data.'
          WRITE (out,130) nAVG(ng), 'nAVG',                             &
     &          'Number of timesteps between the writing of',           &
     &          'time-averaged data into averages file.'
          IF (ndefAVG(ng).gt.0) THEN
            WRITE (out,130) ndefAVG(ng), 'ndefAVG',                     &
     &            'Number of timesteps between creation of new',        &
     &            'time-averaged file.'
          END IF
#endif
#if defined AVERAGES2
          WRITE (out,130) ntsAVG2(ng), 'ntsAVG2',                       &
     &          'Starting timestep for the accumulation of output',     &
     &          'secondary time-averaged data.'
          WRITE (out,130) nAVG2(ng), 'nAVG2',                           &
     &          'Number of timesteps between the writing of',           &
     &          'secondary time-averaged data into averages file.'
          IF (ndefAVG2(ng).gt.0) THEN
            WRITE (out,130) ndefAVG2(ng), 'ndefAVG2',                   &
     &            'Number of timesteps between creation of new',        &
     &            'secondary time-averaged file.'
          END IF
#endif
#if defined HISTORY2
          WRITE (out,130) nHIS2(ng), 'nHIS2',                           &
     &          'Number of timesteps between the writing of',           &
     &          'surface file.'
          IF (ndefHIS2(ng).gt.0) THEN
            WRITE (out,130) ndefHIS2(ng), 'ndefHIS2',                   &
     &            'Number of timesteps between creation of new',        &
     &            'surface file.'
          END IF
#endif
#ifdef DIAGNOSTICS
          WRITE (out,130) ntsDIA(ng), 'ntsDIA',                         &
     &          'Starting timestep for the accumulation of output',     &
     &          'time-averaged diagnostics data.'
          WRITE (out,130) nDIA(ng), 'nDIA',                             &
     &          'Number of timesteps between the writing of',           &
     &          'time-averaged data into diagnostics file.'
          IF (ndefDIA(ng).gt.0) THEN
            WRITE (out,130) ndefDIA(ng), 'ndefDIA',                     &
     &            'Number of timesteps between creation of new',        &
     &            'diagnostic file.'
          END IF
#endif
#ifdef TANGENT
          WRITE (out,170) LcycleTLM(ng), 'LcycleTLM',                   &
     &          'Switch to recycle time-records in tangent file.'
          WRITE (out,130) nTLM(ng), 'nTLM',                             &
     &          'Number of timesteps between the writing of',           &
     &          'data into tangent file.'
          IF (ndefTLM(ng).gt.0) THEN
            WRITE (out,130) ndefTLM(ng), 'ndefTLM',                     &
     &            'Number of timesteps between creation of new',        &
     &            'tangent file.'
          END IF
#endif
#ifdef ADJOINT
          WRITE (out,170) LcycleADJ(ng), 'LcycleADJ',                   &
     &          'Switch to recycle time-records in adjoint file.'
          WRITE (out,130) nADJ(ng), 'nADJ',                             &
     &          'Number of timesteps between the writing of',           &
     &          'data into adjoint file.'
          IF (ndefADJ(ng).gt.0) THEN
            WRITE (out,130) ndefADJ(ng), 'ndefADJ',                     &
     &            'Number of timesteps between creation of new',        &
     &            'adjoint file.'
          END IF
#endif
#ifdef ADJUST_BOUNDARY
          WRITE (out,130) nOBC(ng), 'nOBC',                             &
     &          'Number of timesteps between 4DVAR adjustment of',      &
     &          'open boundaries.'
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          WRITE (out,130) nSFF(ng), 'nSFF',                             &
     &          'Number of timesteps between 4DVAR adjustment of',      &
     &          'surface forcing fields.'
#endif
#ifdef PROPAGATOR
          WRITE (out,170) LmultiGST, 'LmultiGST',                       &
     &          'Switch to write one GST eigenvector per file.'
          WRITE (out,170) LrstGST, 'LrstGST',                           &
     &          'Switch to restart GST analysis.'
          WRITE (out,120) MaxIterGST, 'MaxIterGST',                     &
     &          'Maximum number of GST algorithm iterations.'
          WRITE (out,130) nGST, 'nGST',                                 &
     &          'Number of GST iterations between storing check',       &
     &          'pointing data into NetCDF file.'
          WRITE (out,210) Ritz_tol, 'Ritz_tol',                         &
     &          'Relative accuracy of Ritz values computed in the',     &
     &          'GST analysis.'
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
          DO i=1,NAT+NPT
            itrc=i
#  ifdef T_PASSIVE
            IF (i.gt.NAT) itrc=inert(i-NAT)
#  endif
            WRITE (out,190) nl_tnu2(itrc,ng), 'nl_tnu2', itrc,          &
     &            'NLM Horizontal, harmonic mixing coefficient',        &
     &            '(m2/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
#  ifdef ADJOINT
            WRITE (out,190) ad_tnu2(itrc,ng), 'ad_tnu2', itrc,          &
     &            'ADM Horizontal, harmonic mixing coefficient',        &
     &            '(m2/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
#  endif
#  if defined TANGENT || defined TL_IOMS
            WRITE (out,190) tl_tnu2(itrc,ng), 'tl_tnu2', itrc,          &
     &            'TLM Horizontal, harmonic mixing coefficient',        &
     &            '(m2/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
#  endif
          END DO
# endif
# ifdef TS_DIF4
          DO i=1,NAT+NPT
            itrc=i
#  ifdef T_PASSIVE
            IF (i.gt.NAT) itrc=inert(i-NAT)
#  endif
            WRITE (out,190) nl_tnu4(itrc,ng), 'nl_tnu4', itrc,          &
     &            'NLM Horizontal, biharmonic mixing coefficient',      &
     &            '(m4/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
#  ifdef ADJOINT
            WRITE (out,190) ad_tnu4(itrc,ng), 'ad_tnu4', itrc,          &
     &            'ADM Horizontal, biharmonic mixing coefficient',      &
     &            '(m4/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
#  endif
#  if defined TANGENT || defined TL_IOMS
            WRITE (out,190) tl_tnu4(itrc,ng), 'tl_tnu4', itrc,          &
     &            'TLM Horizontal, biharmonic mixing coefficient',      &
     &            '(m4/s) for tracer ', itrc,                           &
     &            TRIM(Vname(1,idTvar(itrc)))
#  endif
          END DO
# endif
#endif
#ifdef UV_VIS2
          WRITE (out,210) nl_visc2(ng), 'nl_visc2',                     &
     &          'NLM Horizontal, harmonic mixing coefficient',          &
     &          '(m2/s) for momentum.'
# ifdef ADJOINT
          WRITE (out,210) ad_visc2(ng), 'ad_visc2',                     &
     &          'ADM Horizontal, harmonic mixing coefficient',          &
     &          '(m2/s) for momentum.'
# endif
# if defined TANGENT || defined TL_IOMS
          WRITE (out,210) tl_visc2(ng), 'tl_visc2',                     &
     &          'TLM Horizontal, harmonic mixing coefficient',          &
     &          '(m2/s) for momentum.'
# endif
#endif
#ifdef UV_VIS4
          WRITE (out,210) nl_visc4(ng), 'nl_visc4',                     &
     &          'NLM Horizontal, biharmonic mixing coefficient',        &
     &          '(m4/s) for momentum.'
# ifdef ADJOINT
          WRITE (out,210) ad_visc4(ng), 'ad_visc4',                     &
     &          'ADM Horizontal, biharmonic mixing coefficient',        &
     &          '(m4/s) for momentum.'
# endif
# if defined TANGENT || defined TL_IOMS
          WRITE (out,210) tl_visc4(ng), 'tl_visc4',                     &
     &          'TLM Horizontal, biharmonic mixing coefficient',        &
     &          '(m4/s) for momentum.'
# endif
#endif
          IF (LuvSponge(ng)) THEN
            WRITE (out,170) LuvSponge(ng), 'LuvSponge',                 &
     &          'Turning ON  sponge on horizontal momentum.'
          ELSE
            WRITE (out,170) LuvSponge(ng), 'LuvSponge',                 &
     &          'Turning OFF sponge on horizontal momentum.'
          END IF
#ifdef SOLVE3D
          DO i=1,NAT
            IF (LtracerSponge(i,ng)) THEN
              WRITE (out,185) LtracerSponge(i,ng), 'LtracerSponge', i,  &
     &            'Turning ON  sponge on tracer ', i,                   &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LtracerSponge(i,ng), 'LtracerSponge', i,  &
     &            'Turning OFF sponge on tracer ', i,                   &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
# ifdef T_PASSIVE
          DO itrc=1,NPT
            i=inert(itrc)
            IF (LtracerSponge(i,ng)) THEN
              WRITE (out,185) LtracerSponge(i,ng), 'LtracerSponge', i,  &
     &            'Turning ON  sponge on tracer ', i,                   &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LtracerSponge(i,ng), 'LtracerSponge', i,  &
     &            'Turning OFF sponge on tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
# endif
          DO i=1,NAT+NPT
            itrc=i
# ifdef T_PASSIVE
            IF (i.gt.NAT) itrc=inert(i-NAT)
# endif
            WRITE (out,190) Akt_bak(itrc,ng), 'Akt_bak', itrc,          &
     &            'Background vertical mixing coefficient (m2/s)',      &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Akv_bak(ng), 'Akv_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for momentum.'
# if defined MY25_MIXING || defined GLS_MIXING
          WRITE (out,210) Akk_bak(ng), 'Akk_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for turbulent energy.'
          WRITE (out,210) Akp_bak(ng), 'Akp_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for turbulent generic statistical field.'
#  ifdef TKE_DIF2
          WRITE (out,210) tkenu2(ng), 'tkenu2',                         &
     &          'Horizontal, harmonic mixing coefficient (m2/s)',       &
     &          'for turbulent energy.'
#  endif
#  ifdef TKE_DIF4
          WRITE (out,210) tkenu4(ng), 'tkenu4',                         &
     &          'Horizontal, biharmonic mixing coefficient (m4/s)',     &
     &          'for turbulent energy.'
#  endif
# endif
# ifdef GLS_MIXING
          WRITE (out,140) gls_p(ng), 'gls_p',                           &
     &          'GLS stability exponent.'
          WRITE (out,140) gls_m(ng), 'gls_m',                           &
     &          'GLS turbulent kinetic energy exponent.'
          WRITE (out,140) gls_n(ng), 'gls_n',                           &
     &          'GLS turbulent length scale exponent.'
          WRITE (out,200) gls_Kmin(ng), 'gls_Kmin',                     &
     &          'GLS minimum value of turbulent kinetic energy.'
          WRITE (out,200) gls_Pmin(ng), 'gls_Pmin',                     &
     &          'GLS minimum value of dissipation.'
          WRITE (out,200) gls_cmu0(ng), 'gls_cmu0',                     &
     &          'GLS stability coefficient.'
          WRITE (out,200) gls_c1(ng), 'gls_c1',                         &
     &          'GLS shear production coefficient.'
          WRITE (out,200) gls_c2(ng), 'gls_c2',                         &
     &          'GLS dissipation coefficient.'
          WRITE (out,200) gls_c3m(ng), 'gls_c3m',                       &
     &          'GLS stable buoyancy production coefficient.'
          WRITE (out,200) gls_c3p(ng), 'gls_c3p',                       &
     &          'GLS unstable buoyancy production coefficient.'
          WRITE (out,200) gls_sigk(ng), 'gls_sigk',                     &
     &          'GLS constant Schmidt number for TKE.'
          WRITE (out,200) gls_sigp(ng), 'gls_sigp',                     &
     &          'GLS constant Schmidt number for PSI.'
          WRITE (out,140) charnok_alpha(ng), 'charnok_alpha',           &
     &          'Charnok factor for Zos calculation.'
          WRITE (out,140) zos_hsig_alpha(ng), 'zos_hsig_alpha',         &
     &          'Factor for Zos calculation using Hsig(Awave).'
          WRITE (out,140) sz_alpha(ng), 'sz_alpha',                     &
     &          'Factor for Wave dissipation surface tke flux .'
          WRITE (out,140) crgban_cw(ng), 'crgban_cw',                   &
     &          'Factor for Craig/Banner surface tke flux.'
          WRITE (out,140) wec_alpha(ng), 'wec_alpha',                   &
     &          'WEC factor for roller/breaking energy distribution.'
# endif
# ifdef FORWARD_MIXING
          DO i=1,NAT+NPT
            itrc=i
#  ifdef T_PASSIVE
            IF (i.gt.NAT) itrc=inert(i-NAT)
#  endif
#  ifdef ADJOINT
            WRITE (out,190) ad_Akt_fac(itrc,ng), 'ad_Akt_fac', itrc,    &
     &            'ADM basic state vertical mixing scale factor',       &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
#  endif
#  if defined TANGENT || defined TL_IOMS
            WRITE (out,190) tl_Akt_fac(itrc,ng), 'tl_Akt_fac', itrc,    &
     &            'TLM basic state vertical mixing scale factor',       &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
#  endif
          END DO
#  ifdef ADJOINT
          WRITE (out,210) ad_Akv_fac(ng), 'ad_Akv_fac',                 &
     &          'ADM basic state vertical mixing scale factor',         &
     &          'for momentum.'
#  endif
#  if defined TANGENT || defined TL_IOMS
          WRITE (out,210) tl_Akv_fac(ng), 'tl_Akv_fac',                 &
     &          'TLM basic state vertical mixing scale factor',         &
     &          'for momentum.'
#  endif
# endif
#endif
          WRITE (out,200) rdrg(ng), 'rdrg',                             &
     &          'Linear bottom drag coefficient (m/s).'
          WRITE (out,200) rdrg2(ng), 'rdrg2',                           &
     &          'Quadratic bottom drag coefficient.'
          WRITE (out,200) Zob(ng), 'Zob',                               &
     &          'Bottom roughness (m).'
#ifdef BBL_MODEL
          IF (Zob(ng).le.0.0_r8) THEN
            WRITE (out,265) 'Zob = ', Zob(ng),                          &
     &            'It must be greater than zero when BBL is activated.'
            exit_flag=5
            RETURN
          END IF
#endif
#ifdef SOLVE3D
# ifdef GLS_MIXING
          WRITE (out,200) Zos(ng), 'Zos',                               &
     &          'Surface roughness (m).'
# endif
# ifdef BULK_FLUXES
          WRITE (out,200) blk_ZQ(ng), 'blk_ZQ',                         &
     &          'Height (m) of surface air humidity measurement.'
          IF (blk_ZQ(ng).le.0.0_r8) THEN
            WRITE (out,265) 'blk_ZQ = ', blk_ZQ(ng),                    &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZT(ng), 'blk_ZT',                         &
     &          'Height (m) of surface air temperature measurement.'
          IF (blk_ZT(ng).le.0.0_r8) THEN
            WRITE (out,265) 'blk_ZT = ', blk_ZT(ng),                    &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZW(ng), 'blk_ZW',                         &
     &          'Height (m) of surface winds measurement.'
          IF (blk_ZW(ng).le.0.0_r8) THEN
            WRITE (out,265) 'blk_ZW = ', blk_ZW(ng),                    &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
# endif
#endif
#if defined WET_DRY
          WRITE (out,200) Dcrit(ng), 'Dcrit',                           &
     &          'Minimum depth for wetting and drying (m).'
#endif
#ifdef SOLVE3D
# if defined LMD_SKPP || defined SOLAR_SOURCE
          WRITE (out,120) lmd_Jwt(ng), 'lmd_Jwt',                       &
     &          'Jerlov water type.'
          IF ((lmd_Jwt(ng).lt.1).or.(lmd_Jwt(ng).gt.9)) THEN
            WRITE (out,260) 'lmd_Jwt = ', lmd_Jwt(ng),                  &
     &            'It must between one and nine.'
            exit_flag=5
            RETURN
          END IF
# endif
# ifdef BODYFORCE
          WRITE (out,130) levsfrc(ng), 'levsfrc',                       &
     &          'Deepest level to apply surface stress as a',           &
     &          'bodyforce.'
          IF ((levsfrc(ng).lt.1).or.(levsfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levsfrc = ', levsfrc(ng),                  &
     &            'Out of range surface bodyforce level.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,130) levbfrc(ng), 'levbfrc',                       &
     &          'Shallowest level to apply bottom stress as a',         &
     &          'bodyforce.'
          IF ((levbfrc(ng).lt.1).or.(levbfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levbfrc = ', levbfrc(ng),                  &
     &            'Out of range bottom bodyforce level.'
            exit_flag=5
            RETURN
          END IF
# endif
#endif
#ifdef SOLVE3D
          WRITE (out,120) Vtransform(ng), 'Vtransform',                 &
     &          'S-coordinate transformation equation.'
          WRITE (out,120) Vstretching(ng), 'Vstretching',               &
     &          'S-coordinate stretching function.'
          WRITE (out,200) theta_s(ng), 'theta_s',                       &
     &          'S-coordinate surface control parameter.'
          WRITE (out,200) theta_b(ng), 'theta_b',                       &
     &          'S-coordinate bottom  control parameter.'
          IF (Tcline(ng).gt.1.0E+5_r8) THEN
            WRITE (out,210) Tcline(ng), 'Tcline',                       &
     &            'S-coordinate surface/bottom layer width (m) used',   &
     &            'in vertical coordinate stretching.'
          ELSE
            WRITE (out,160) Tcline(ng), 'Tcline',                       &
     &            'S-coordinate surface/bottom layer width (m) used',   &
     &            'in vertical coordinate stretching.'
          END IF
#endif
          WRITE (out,140) rho0, 'rho0',                                 &
     &          'Mean density (kg/m3) for Boussinesq approximation.'
#if defined SOLVE3D && (defined FOUR_DVAR || defined PROPAGATOR)
          WRITE (out,200) bvf_bak, 'bvf_bak',                           &
     &          'Background Brunt-Vaisala frequency squared (1/s2).'
#endif
          WRITE (out,140) dstart, 'dstart',                             &
     &          'Time-stamp assigned to model initialization (days).'
#if defined SSH_TIDES || defined UV_TIDES
          WRITE (out,140) tide_start, 'tide_start',                     &
     &          'Reference time origin for tidal forcing (days).'
#endif
          WRITE (out,150) time_ref, 'time_ref',                         &
     &          'Reference time for units attribute (yyyymmdd.dd)'
#ifdef SOLVE3D
          DO i=1,NAT+NPT
            itrc=i
#  ifdef T_PASSIVE
            IF (i.gt.NAT) itrc=inert(i-NAT)
#  endif
            WRITE (out,190) Tnudg(itrc,ng), 'Tnudg', itrc,              &
     &            'Nudging/relaxation time scale (days)',               &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Tnudg_SSS(ng), 'Tnudg_SSS',                   &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for sea surface salinity.'
# ifdef SCORRECTION
          IF (Tnudg_SSS(ng).le.0.0_r8) THEN
            WRITE (out,265) 'Tnudg_SSS = ', Tnudg_SSS(ng),              &
     &            'Must be greater than zero for salt flux correction.'
            exit_flag=5
            RETURN
          END IF
#  ifdef SSSC_THRESHOLD
          WRITE (out,210) sss_mismatch_threshold(ng),                   &
     &          'sss_mismatch_threshold',                               &
     &          'difference in SSS that triggers correction',           &
     &          'for sea surface salinity.'
#  endif
# endif
#endif
          WRITE (out,210) Znudg(ng), 'Znudg',                           &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for free-surface.'
          WRITE (out,210) M2nudg(ng), 'M2nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 2D momentum.'
#ifdef SOLVE3D
          WRITE (out,210) M3nudg(ng), 'M3nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 3D momentum.'
#endif
          WRITE (out,210) obcfac(ng), 'obcfac',                         &
     &          'Factor between passive and active',                    &
     &          'open boundary conditions.'
          WRITE (out,170) VolCons(1,ng), 'VolCons(1)',                  &
     &          'NLM western  edge boundary volume conservation.'
          WRITE (out,170) VolCons(2,ng), 'VolCons(2)',                  &
     &          'NLM southern edge boundary volume conservation.'
          WRITE (out,170) VolCons(3,ng), 'VolCons(3)',                  &
     &          'NLM eastern  edge boundary volume conservation.'
          WRITE (out,170) VolCons(4,ng), 'VolCons(4)',                  &
     &          'NLM northern edge boundary volume conservation.'
#ifdef ADJOINT
          WRITE (out,170) ad_VolCons(1,ng), 'ad_VolCons(1)',            &
     &          'ADM western  edge boundary volume conservation.'
          WRITE (out,170) ad_VolCons(2,ng), 'ad_VolCons(2)',            &
     &          'ADM southern edge boundary volume conservation.'
          WRITE (out,170) ad_VolCons(3,ng), 'ad_VolCons(3)',            &
     &          'ADM eastern  edge boundary volume conservation.'
          WRITE (out,170) ad_VolCons(4,ng), 'ad_VolCons(4)',            &
     &          'ADM northern edge boundary volume conservation.'
#endif
#if defined TANGENT || defined TL_IOMS
          WRITE (out,170) tl_VolCons(1,ng), 'tl_VolCons(1)',            &
     &          'TLM western  edge boundary volume conservation.'
          WRITE (out,170) tl_VolCons(2,ng), 'tl_VolCons(2)',            &
     &          'TLM southern edge boundary volume conservation.'
          WRITE (out,170) tl_VolCons(3,ng), 'tl_VolCons(3)',            &
     &          'TLM eastern  edge boundary volume conservation.'
          WRITE (out,170) tl_VolCons(4,ng), 'tl_VolCons(4)',            &
     &          'TLM northern edge boundary volume conservation.'
#endif
#ifdef SOLVE3D
          WRITE (out,140) T0(ng), 'T0',                                 &
     &          'Background potential temperature (C) constant.'
          WRITE (out,140) S0(ng), 'S0',                                 &
     &          'Background salinity (PSU) constant.'
# ifndef NONLIN_EOS
          WRITE (out,160) R0(ng), 'R0',                                 &
     &          'Background density (kg/m3) used in linear Equation',   &
     &          'of State.'
# endif
# if !defined NONLIN_EOS || defined FOUR_DVAR || defined PROPAGATOR
          WRITE (out,200) Tcoef(ng), 'Tcoef',                           &
     &          'Thermal expansion coefficient (1/Celsius).'
          WRITE (out,200) Scoef(ng), 'Scoef',                           &
     &          'Saline contraction coefficient (1/PSU).'
# endif
#endif
          WRITE (out,160) gamma2(ng), 'gamma2',                         &
     &          'Slipperiness variable: free-slip (1.0) or ',           &
     &          '                     no-slip (-1.0).'
          IF (LuvSrc(ng)) THEN
            WRITE (out,170) LuvSrc(ng), 'LuvSrc',                       &
     &          'Turning ON  momentum point Sources/Sinks.'
          ELSE
            WRITE (out,170) LuvSrc(ng), 'LuvSrc',                       &
     &          'Turning OFF momentum point Sources/Sinks.'
          END IF
          IF (LwSrc(ng)) THEN
            WRITE (out,170) LwSrc(ng), 'LwSrc',                         &
     &          'Turning ON  volume influx point Sources/Sinks.'
          ELSE
            WRITE (out,170) LwSrc(ng), 'LwSrc',                         &
     &          'Turning OFF volume influx point Sources/Sinks.'
          END IF
#ifdef SOLVE3D
          DO itrc=1,NAT
            IF (LtracerSrc(itrc,ng)) THEN
              WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,  &
     &            'Turning ON  point Sources/Sinks on tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
            ELSE
              WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,  &
     &            'Turning OFF point Sources/Sinks on tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
            END IF
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (LtracerSrc(itrc,ng)) THEN
              WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,  &
     &            'Turning ON  point Sources/Sinks on tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
            ELSE
              WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,  &
     &            'Turning OFF  point Sources/Sinks on tracer ', itrc,  &
     &            TRIM(Vname(1,idTvar(itrc)))
            END IF
          END DO
# endif
#endif
          IF (LsshCLM(ng)) THEN
            WRITE (out,170) LsshCLM(ng), 'LsshCLM',                     &
     &          'Turning ON  processing of SSH climatology.'
          ELSE
            WRITE (out,170) LsshCLM(ng), 'LsshCLM',                     &
     &          'Turning OFF processing of SSH climatology.'
          END IF
          IF (Lm2CLM(ng)) THEN
            WRITE (out,170) Lm2CLM(ng), 'Lm2CLM',                       &
     &          'Turning ON  processing of 2D momentum climatology.'
          ELSE
            WRITE (out,170) Lm2CLM(ng), 'Lm2CLM',                       &
     &          'Turning OFF processing of 2D momentum climatology.'
          END IF
#ifdef SOLVE3D
          IF (Lm3CLM(ng)) THEN
            WRITE (out,170) Lm3CLM(ng), 'Lm3CLM',                       &
     &          'Turning ON  processing of 3D momentum climatology.'
          ELSE
            WRITE (out,170) Lm3CLM(ng), 'Lm3CLM',                       &
     &          'Turning OFF processing of 3D momentum climatology.'
          END IF
          DO i=1,NAT
            IF (LtracerCLM(i,ng)) THEN
              WRITE (out,185) LtracerCLM(i,ng), 'LtracerCLM', i,        &
     &            'Turning ON  processing of climatology tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LtracerCLM(i,ng), 'LtracerCLM', i,        &
     &            'Turning OFF processing of climatology tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
# ifdef T_PASSIVE
          DO itrc=1,NPT
            i=inert(itrc)
            IF (LtracerCLM(i,ng)) THEN
              WRITE (out,185) LtracerCLM(i,ng), 'LtracerCLM', i,        &
     &            'Turning ON  processing of climatology tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LtracerCLM(i,ng), 'LtracerCLM', i,        &
     &            'Turning OFF processing of climatology tracer ', i,   &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
# endif
#endif
          IF (LnudgeM2CLM(ng)) THEN
            WRITE (out,170) LnudgeM2CLM(ng), 'LnudgeM2CLM',             &
     &          'Turning ON  nudging of 2D momentum climatology.'
          ELSE
            WRITE (out,170) LnudgeM2CLM(ng), 'LnudgeM2CLM',             &
     &          'Turning OFF nudging of 2D momentum climatology.'
          END IF
#ifdef SOLVE3D
          IF (LnudgeM3CLM(ng)) THEN
            WRITE (out,170) LnudgeM3CLM(ng), 'LnudgeM3CLM',             &
     &          'Turning ON  nudging of 3D momentum climatology.'
          ELSE
            WRITE (out,170) LnudgeM3CLM(ng), 'LnudgeM3CLM',             &
     &          'Turning OFF nudging of 3D momentum climatology.'
          END IF
          DO i=1,NAT
            IF (LnudgeTCLM(i,ng)) THEN
              WRITE (out,185) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,        &
     &            'Turning ON  nudging of climatology tracer ', i,      &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,        &
     &            'Turning OFF nudging of climatology tracer ', i,      &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
# ifdef T_PASSIVE
          DO itrc=1,NPT
            i=inert(itrc)
            IF (LnudgeTCLM(i,ng)) THEN
              WRITE (out,185) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,        &
     &            'Turning ON  nudging of climatology tracer ', i,      &
     &            TRIM(Vname(1,idTvar(i)))
            ELSE
              WRITE (out,185) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,        &
     &            'Turning OFF nudging of climatology tracer ', i,      &
     &            TRIM(Vname(1,idTvar(i)))
            END IF
          END DO
# endif
#endif
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
          WRITE (out,140) DstrS(ng), 'DstrS',                           &
     &          'Starting day for ADM sensitivity forcing.'
          WRITE (out,140) DendS(ng), 'DendS',                           &
     &          'Ending day for ADM sensitivity forcing.'
# ifdef SOLVE3D
          WRITE (out,120) KstrS(ng), 'KstrS',                           &
     &          'Deepest level whose ADM sensitivity is required.'
          IF ((KstrS(ng).lt.1).or.(KstrS(ng).gt.N(ng))) THEN
            WRITE (out,260) 'KstrS = ', KstrS(ng),                      &
     &            'Out of range ADM sensitivity starting level.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,120) KendS(ng), 'KendS',                           &
     &          'Shallowest level whose ADM sensitivity is required.'
          IF ((KendS(ng).lt.1).or.(KendS(ng).gt.N(ng))) THEN
            WRITE (out,260) 'KendS = ', KendS(ng),                      &
     &            'Out of range ADM sensitivity level.'
            exit_flag=5
            RETURN
          END IF
# endif
          IF (SCALARS(ng)%Lstate(isFsur))                               &
     &       WRITE (out,170) SCALARS(ng)%Lstate(isFsur),                &
     &       'Lstate(isFsur)',                                          &
     &       'Adjoint sensitivity on free-surface.'
          IF (SCALARS(ng)%Lstate(isUbar))                               &
     &       WRITE (out,170) SCALARS(ng)%Lstate(isUbar),                &
     &       'Lstate(isUbar)',                                          &
     &       'Adjoint sensitivity on 2D U-momentum component.'
          IF (SCALARS(ng)%Lstate(isVbar))                               &
     &       WRITE (out,170) SCALARS(ng)%Lstate(isVbar),                &
     &       'Lstate(isVbar)',                                          &
     &       'Adjoint sensitivity on 2D V-momentum component.'
# ifdef SOLVE3D
          IF (SCALARS(ng)%Lstate(isUvel))                               &
     &       WRITE (out,170) SCALARS(ng)%Lstate(isUvel),                &
     &       'Lstate(isUvel)',                                          &
     &       'Adjoint sensitivity on 3D U-momentum component.'
          IF (SCALARS(ng)%Lstate(isVvel))                               &
     &       WRITE (out,170) SCALARS(ng)%Lstate(isVvel),                &
     &       'Lstate(isVvel)',                                          &
     &       'Adjoint sensitivity on 3D V-momentum component.'
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%Lstate(isTvar(itrc)))                       &
     &         WRITE (out,180) SCALARS(ng)%Lstate(isTvar(itrc)),        &
     &         'Lstate(idTvar)',                                        &
     &         'Adjoint sensitivity on tracer ',                        &
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#if defined FORCING_SV || defined SO_SEMI || defined STOCHASTIC_OPT
          IF (SCALARS(ng)%Fstate(isFsur))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isFsur),                &
     &       'Fstate(isFsur)',                                          &
# ifdef FORCING_SV
     &       'Singular Forcing Vectors on free-surface.'
# else
     &       'Stochastic Optimals on free-surface.'
# endif
# ifdef SOLVE3D
          IF (SCALARS(ng)%Fstate(isUvel))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isUvel),                &
     &       'Fstate(isUvel)',                                          &
#  ifdef FORCING_SV
     &       'Singular Forcing Vectors on 3D U-momentum component.'
#  else
     &       'Stochastic Optimals on 3D U-momentum component.'
#  endif
          IF (SCALARS(ng)%Fstate(isVvel))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isVvel),                &
     &       'Fstate(isVvel)',                                          &
#  ifdef FORCING_SV
     &       'Singular Forcing Vectors on 3D V-momentum component.'
#  else
     &       'Stochastic Optimals on 3D V-momentum component.'
#  endif
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%Fstate(isTvar(itrc)))                       &
     &         WRITE (out,180) SCALARS(ng)%Fstate(isTvar(itrc)),        &
     &         'Fstate(idTVar)',                                        &
#  ifdef FORCING_SV
     &         'Singular Forcing Vectors on tracer',                    &
#  else
     &         'Stochastic Optimals on tracer',                         &
#  endif
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# else
          IF (SCALARS(ng)%Fstate(isUbar))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isUbar),                &
     &       'Fstate(isUbar)',                                          &
#  ifdef FORCING_SV
     &       'Singular Forcing Vectors on 2D U-momentum component.'
#  else
     &       'Stochastic Optimals on 2D U-momentum component.'
#  endif
          IF (SCALARS(ng)%Fstate(isVbar))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isVbar),                &
     &       'Fstate(isVbar)',                                          &
#  ifdef FORCING_SV
     &       'Singular Forcing Vectors on 2D V-momentum component.'
#  else
     &       'Stochastic Optimals on 2D V-momentum component.'
#  endif
# endif
          IF (SCALARS(ng)%Fstate(isUstr))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isUstr),                &
     &       'Fstate(isUstr)',                                          &
# ifdef FORCING_SV
     &       'Singular Forcing Vectors on surface U-stress.'
# else
     &       'Stochastic Optimals on surface U-stress.'
# endif
          IF (SCALARS(ng)%Fstate(isVstr))                               &
     &       WRITE (out,170) SCALARS(ng)%Fstate(isVstr),                &
     &       'Fstate(isVstr)',                                          &
# ifdef FORCING_SV
     &       'Singular Forcing Vectors on surface V-stress.'
# else
     &       'Stochastic Optimals on surface V-stress.'
# endif
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%Fstate(isTsur(itrc)))                       &
     &         WRITE (out,180) SCALARS(ng)%Fstate(isTsur(itrc)),        &
     &         'Fstate(idTsur)',                                        &
#  ifdef FORCING_SV
     &         'Singular Forcing Vectors on surface flux of tracer',    &
#  else
     &         'Stochastic Optimals on surface flux of tracer',         &
#  endif
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#ifdef SO_SEMI
# ifndef SO_SEMI_WHITE
          WRITE (out,140) SO_decay(ng), 'SO_decay',                     &
     &          'Stochastic optimals time decorrelation scale (days).'
# endif
          IF (SCALARS(ng)%Fstate(isFsur))                               &
            WRITE (out,200) SO_sdev(isFsur,ng), 'SO_sdev(isFsur)',      &
     &          'Stochastic Optimals scale, free-surface'
# ifdef SOLVE3D
          IF (SCALARS(ng)%Fstate(isUvel))                               &
            WRITE (out,200) SO_sdev(isUvel,ng), 'SO_sdev(isUvel)',      &
     &          'Stochastic Optimals scale, 3D U-momentum'
          IF (SCALARS(ng)%Fstate(isVvel))                               &
            WRITE (out,200) SO_sdev(isVvel,ng), 'SO_sdev(isVvel)',      &
     &          'Stochastic Optimals scale, 3D V-momentum'
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%Fstate(isTvar(itrc)))                       &
     &         WRITE (out,195) SO_sdev(isTvar(itrc),ng),                &
     &         'SO_sdev(idTvar)',                                       &
     &         'Stochastic Optimals scale, tracer',                     &
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# else
          IF (SCALARS(ng)%Fstate(isUbar))                               &
            WRITE (out,200) SO_sdev(isUbar,ng), 'SO_sdev(isUbar)',      &
     &          'Stochastic Optimals scale, 2D U-momentum'
          IF (SCALARS(ng)%Fstate(isVbar))                               &
            WRITE (out,200) SO_sdev(isVbar,ng), 'SO_sdev(isVbar)',      &
     &          'Stochastic Optimals scale, 2D V-momentum'
# endif
          IF (SCALARS(ng)%Fstate(isUstr))                               &
            WRITE (out,200) SO_sdev(isUstr,ng), 'SO_sdev(isUstr)',      &
     &          'Stochastic Optimals scale, surface U-stress'
          IF (SCALARS(ng)%Fstate(isVstr))                               &
            WRITE (out,200) SO_sdev(isVstr,ng), 'SO_sdev(isVstr)',      &
     &          'Stochastic Optimals scale, surface V-stress'
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%Fstate(isTsur(itrc)))                       &
     &         WRITE (out,195) SO_sdev(isTsur(itrc),ng),                &
     &         'SO_sdev(idTsur)',                                       &
     &         'Stochastic Optimals scale, surface flux of tracer',     &
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#if defined SEDIMENT && defined SED_MORPH
          IF (Hout(idBath,ng)) WRITE (out,170) Hout(idBath,ng),         &
     &       'Hout(idBath)',                                            &
     &       'Write out time-dependent bathymetry.'
#endif
          IF (Hout(idFsur,ng)) WRITE (out,170) Hout(idFsur,ng),         &
     &       'Hout(idFsur)',                                            &
     &       'Write out free-surface.'
          IF (Hout(idUbar,ng)) WRITE (out,170) Hout(idUbar,ng),         &
     &       'Hout(idUbar)',                                            &
     &       'Write out 2D U-momentum component.'
          IF (Hout(idVbar,ng)) WRITE (out,170) Hout(idVbar,ng),         &
     &       'Hout(idVbar)',                                            &
     &       'Write out 2D V-momentum component.'
          IF (Hout(idu2dE,ng)) WRITE (out,170) Hout(idu2dE,ng),         &
     &       'Hout(idu2dE)',                                            &
     &       'Write out 2D U-eastward  at RHO-points.'
          IF (Hout(idv2dN,ng)) WRITE (out,170) Hout(idv2dN,ng),         &
     &       'Hout(idv2dN)',                                            &
     &       'Write out 2D V-northward at RHO-points.'
#ifdef SOLVE3D
          IF (Hout(idUvel,ng)) WRITE (out,170) Hout(idUvel,ng),         &
     &       'Hout(idUvel)',                                            &
     &       'Write out 3D U-momentum component.'
          IF (Hout(idVvel,ng)) WRITE (out,170) Hout(idVvel,ng),         &
     &       'Hout(idVvel)',                                            &
     &       'Write out 3D V-momentum component.'
          IF (Hout(idu3dE,ng)) WRITE (out,170) Hout(idu3dE,ng),         &
     &       'Hout(idu3dE)',                                            &
     &       'Write out 3D U-wastward  component at RHO-points.'
          IF (Hout(idv3dN,ng)) WRITE (out,170) Hout(idv3dN,ng),         &
     &       'Hout(idv3dN)',                                            &
     &       'Write out 3D V-northward component at RHO-points.'
          IF (Hout(idWvel,ng)) WRITE (out,170) Hout(idWvel,ng),         &
     &       'Hout(idWvel)',                                            &
     &       'Write out W-momentum component.'
          IF (Hout(idOvel,ng)) WRITE (out,170) Hout(idOvel,ng),         &
     &       'Hout(idOvel)',                                            &
     &       'Write out omega vertical velocity.'
          DO itrc=1,NAT
            IF (Hout(idTvar(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idTvar(itrc),ng), 'Hout(idTvar)',                  &
     &          'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          IF (Hout(idpthR,ng)) WRITE (out,170) Hout(idpthR,ng),         &
     &       'Hout(idpthR)',                                            &
     &       'Write out time-varying dephts of RHO-points.'
          IF (Hout(idpthU,ng)) WRITE (out,170) Hout(idpthU,ng),         &
     &       'Hout(idpthU)',                                            &
     &       'Write out time-varying dephts of U-points.'
          IF (Hout(idpthV,ng)) WRITE (out,170) Hout(idpthV,ng),         &
     &       'Hout(idpthV)',                                            &
     &       'Write out time-varying dephts of V-points.'
          IF (Hout(idpthW,ng)) WRITE (out,170) Hout(idpthW,ng),         &
     &       'Hout(idpthW)',                                            &
     &       'Write out time-varying dephts of W-points.'
#endif
          IF (Hout(idUsms,ng)) WRITE (out,170) Hout(idUsms,ng),         &
     &       'Hout(idUsms)',                                            &
     &       'Write out surface U-momentum stress.'
          IF (Hout(idVsms,ng)) WRITE (out,170) Hout(idVsms,ng),         &
     &       'Hout(idVsms)',                                            &
     &       'Write out surface V-momentum stress.'
          IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idUbms,ng),         &
     &       'Hout(idUbms)',                                            &
     &       'Write out bottom U-momentum stress.'
          IF (Hout(idVbms,ng)) WRITE (out,170) Hout(idVbms,ng),         &
     &       'Hout(idVbms)',                                            &
     &       'Write out bottom V-momentum stress.'
#ifdef BBL_MODEL
          IF (Hout(idUbrs,ng)) WRITE (out,170) Hout(idUbrs,ng),         &
     &       'Hout(idUbrs)',                                            &
     &       'Write out bottom U-current stress.'
          IF (Hout(idVbrs,ng)) WRITE (out,170) Hout(idVbrs,ng),         &
     &       'Hout(idVbrs)',                                            &
     &       'Write out bottom V-current stress.'
          IF (Hout(idUbws,ng)) WRITE (out,170) Hout(idUbws,ng),         &
     &       'Hout(idUbws)',                                            &
     &       'Write out wind-induced, bottom U-wave stress.'
          IF (Hout(idVbws,ng)) WRITE (out,170) Hout(idVbws,ng),         &
     &       'Hout(idVbws)',                                            &
     &       'Write out wind-induced, bottom V-wave stress.'
          IF (Hout(idUbcs,ng)) WRITE (out,170) Hout(idUbcs,ng),         &
     &       'Hout(idUbcs)',                                            &
     &       'Write out max wind + current, bottom U-wave stress.'
          IF (Hout(idVbcs,ng)) WRITE (out,170) Hout(idVbcs,ng),         &
     &       'Hout(idVbcs)',                                            &
     &       'Write out max wind + current, bottom V-wave stress.'
          IF (Hout(idUVwc,ng)) WRITE (out,170) Hout(idUVwc,ng),         &
     &       'Hout(idUVwc)',                                            &
     &       'Write out max wind + current, bottom UV-wave stress.'
          IF (Hout(idUbot,ng)) WRITE (out,170) Hout(idUbot,ng),         &
     &       'Hout(idUbot)',                                            &
     &       'Write out bed wave orbital U-velocity.'
          IF (Hout(idVbot,ng)) WRITE (out,170) Hout(idVbot,ng),         &
     &       'Hout(idVbot)',                                            &
     &       'Write out bed wave orbital V-velocity.'
          IF (Hout(idUbur,ng)) WRITE (out,170) Hout(idUbur,ng),         &
     &       'Hout(idUbur)',                                            &
     &       'Write out bottom U-momentum above bed.'
          IF (Hout(idVbvr,ng)) WRITE (out,170) Hout(idVbvr,ng),         &
     &       'Hout(idVbvr)',                                            &
     &       'Write out bottom V-momentum above bed.'
#endif
#if defined WEC_MELLOR
          IF (Hout(idW2xx,ng)) WRITE (out,170) Hout(idW2xx,ng),         &
     &       'Hout(idW2xx)',                                            &
     &       'Write out 2D radiation stress, Sxx.'
          IF (Hout(idW2xy,ng)) WRITE (out,170) Hout(idW2xy,ng),         &
     &       'Hout(idW2xy)',                                            &
     &       'Write out 2D radiation stress, Sxy.'
          IF (Hout(idW2yy,ng)) WRITE (out,170) Hout(idW2yy,ng),         &
     &       'Hout(idW2yy)',                                            &
     &       'Write out 2D radiation stress, Syy.'
#endif
#if defined WEC
          IF (Hout(idU2rs,ng)) WRITE (out,170) Hout(idU2rs,ng),         &
     &       'Hout(idU2rs)',                                            &
     &       'Write out total 2D u-radiation stress.'
          IF (Hout(idV2rs,ng)) WRITE (out,170) Hout(idV2rs,ng),         &
     &       'Hout(idV2rs)',                                            &
     &       'Write out total 2D v-radiation stress.'
          IF (Hout(idU2Sd,ng)) WRITE (out,170) Hout(idU2Sd,ng),         &
     &       'Hout(idU2Sd)',                                            &
     &       'Write out 2D u-momentum stokes velocity.'
          IF (Hout(idV2Sd,ng)) WRITE (out,170) Hout(idV2Sd,ng),         &
     &       'Hout(idV2Sd)',                                            &
     &       'Write out 2D v-momentum stokes velocity.'
#endif
#ifdef SOLVE3D
# ifdef WEC_MELLOR
          IF (Hout(idW3xx,ng)) WRITE (out,170) Hout(idW3xx,ng),         &
     &       'Hout(idW3xx)',                                            &
     &       'Write out 3D radiation stress, Sxx.'
          IF (Hout(idW3xy,ng)) WRITE (out,170) Hout(idW3xy,ng),         &
     &       'Hout(idW3xy)',                                            &
     &       'Write out 3D radiation stress, Sxy.'
          IF (Hout(idW3yy,ng)) WRITE (out,170) Hout(idW3yy,ng),         &
     &       'Hout(idW3yy)',                                            &
     &       'Write out 3D radiation stress, Syy.'
          IF (Hout(idW3zx,ng)) WRITE (out,170) Hout(idW3zx,ng),         &
     &       'Hout(idW3zx)',                                            &
     &       'Write out 3D radiation stress, Szx.'
          IF (Hout(idW3zy,ng)) WRITE (out,170) Hout(idW3zy,ng),         &
     &       'Hout(idW3zy)',                                            &
     &       'Write out 3D radiation stress, Szy.'
#  endif
#  ifdef WEC
          IF (Hout(idU3rs,ng)) WRITE (out,170) Hout(idU3rs,ng),         &
     &       'Hout(idU3rs)',                                            &
     &       'Write out total 3D u-radiation stress.'
          IF (Hout(idV3rs,ng)) WRITE (out,170) Hout(idV3rs,ng),         &
     &       'Hout(idV3rs)',                                            &
     &       'Write out total 3D v-radiation stress.'
          IF (Hout(idU3Sd,ng)) WRITE (out,170) Hout(idU3Sd,ng),         &
     &       'Hout(idU3Sd)',                                            &
     &       'Write out 3D u-momentum stokes velocity.'
          IF (Hout(idV3Sd,ng)) WRITE (out,170) Hout(idV3Sd,ng),         &
     &       'Hout(idV3Sd)',                                            &
     &       'Write out 3D v-momentum stokes velocity.'
          IF (Hout(idW3Sd,ng)) WRITE (out,170) Hout(idW3Sd,ng),         &
     &       'Hout(idW3Sd)',                                            &
     &       'Write out 3D omega-momentum stokes velocity.'
          IF (Hout(idW3St,ng)) WRITE (out,170) Hout(idW3St,ng),         &
     &       'Hout(idW3St)',                                            &
     &       'Write out 3D w-momentum stokes velocity.'
# endif
# ifdef WEC_VF
          IF (Hout(idWztw,ng)) WRITE (out,170) Hout(idWztw,ng),         &
     &       'Hout(idWztw)',                                            &
     &       'Write out WEC quasi-static sea level adjustment.'
          IF (Hout(idWqsp,ng)) WRITE (out,170) Hout(idWqsp,ng),         &
     &       'Hout(idWqsp)',                                            &
     &       'Write out WEC quasi-static pressure.'
          IF (Hout(idWbeh,ng)) WRITE (out,170) Hout(idWbeh,ng),         &
     &       'Hout(idWbeh)',                                            &
     &       'Write out WEC Bernoulli head.'
# endif
#endif
#ifdef WAVES_HEIGHT
          IF (Hout(idWamp,ng)) WRITE (out,170) Hout(idWamp,ng),         &
     &       'Hout(idWamp)',                                            &
     &       'Write out wave height.'
#endif
#ifdef WAVES_LENGTH
          IF (Hout(idWlen,ng)) WRITE (out,170) Hout(idWlen,ng),         &
     &       'Hout(idWlen)',                                            &
     &       'Write out wavelength.'
#endif
#ifdef WAVES_LENGTHP
          IF (Hout(idWlep,ng)) WRITE (out,170) Hout(idWlep,ng),         &
     &       'Hout(idWlep)',                                            &
     &       'Write out peak wavelength.'
#endif
#ifdef WAVES_DIR
          IF (Hout(idWdir,ng)) WRITE (out,170) Hout(idWdir,ng),         &
     &       'Hout(idWdir)',                                            &
     &       'Write out wave direction.'
#endif
#ifdef WAVES_TOP_PERIOD
          IF (Hout(idWptp,ng)) WRITE (out,170) Hout(idWptp,ng),         &
     &       'Hout(idWptp)',                                            &
     &       'Write out wave surface period.'
#endif
#ifdef WAVES_BOT_PERIOD
          IF (Hout(idWpbt,ng)) WRITE (out,170) Hout(idWpbt,ng),         &
     &       'Hout(idWpbt)',                                            &
     &       'Write out wave bottom period.'
#endif
#ifdef WAVES_OCEAN
          IF (Hout(idWorb,ng)) WRITE (out,170) Hout(idWorb,ng),         &
     &       'Hout(idWorb)',                                            &
     &       'Write out wave bottom orbital velocity.'
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          IF (Hout(idWdif,ng)) WRITE (out,170) Hout(idWdif,ng),         &
     &       'Hout(idWdif)',                                            &
     &       'Write out wave dissipation due to bottom friction.'
#endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          IF (Hout(idWdib,ng)) WRITE (out,170) Hout(idWdib,ng),         &
     &       'Hout(idWdib)',                                            &
     &       'Write out wave dissipation due to breaking.'
          IF (Hout(idWdiw,ng)) WRITE (out,170) Hout(idWdiw,ng),         &
     &       'Hout(idWdiw)',                                            &
     &       'Write out wave dissipation due to white capping.'
#endif
#ifdef ROLLER_SVENDSEN
          IF (Hout(idWbrk,ng)) WRITE (out,170) Hout(idWbrk,ng),         &
     &       'Hout(idWbrk)',                                            &
     &       'Write out percent wave breaking.'
#endif
#ifdef WEC_ROLLER
          IF (Hout(idWdis,ng)) WRITE (out,170) Hout(idWdis,ng),         &
     &       'Hout(idWdis)',                                            &
     &       'Write out wave roller dissipation.'
#endif
#ifdef ROLLER_RENIERS
          IF (Hout(idWrol,ng)) WRITE (out,170) Hout(idWrol,ng),         &
     &       'Hout(idWrol)',                                            &
     &       'Write out wave roller action density.'
#endif
#ifdef WAVES_DSPR
          IF (Hout(idWvds,ng)) WRITE (out,170) Hout(idWvds,ng),         &
     &       'Hout(idWvds)',                                            &
     &       'Write out wave directional spread.'
          IF (Hout(idWvqp,ng)) WRITE (out,170) Hout(idWvqp,ng),         &
     &       'Hout(idWvqp)',                                            &
     &       'Write out wave spectrum peakedness.'
#endif
#ifdef UV_KIRBY
          IF (Hout(idUwav,ng)) WRITE (out,170) Hout(idUwav,ng),         &
     &       'Hout(idUwav)',                                            &
     &       'Wave-avg surface u-velocity.'
          IF (Hout(idVwav,ng)) WRITE (out,170) Hout(idVwav,ng),         &
     &       'Hout(idVwav)',                                            &
     &       'Wave-avg surface v-velocity.'
#endif
#ifdef INWAVE_MODEL
          IF (Hout(idACen,ng)) WRITE (out,170) Hout(idACen,ng),         &
     &       'Hout(idACen)',                                            &
     &       'Write out 3D wave action'
          IF (Hout(idACCt,ng)) WRITE (out,170) Hout(idACct,ng),         &
     &       'Hout(idACct)',                                            &
     &       'Wave group celerity in the theta coordinate.'
          IF (Hout(idACcx,ng)) WRITE (out,170) Hout(idACcx,ng),         &
     &       'Hout(idACcx)',                                            &
     &       'Wave group celerity in the xi coordinate.'
          IF (Hout(idACcy,ng)) WRITE (out,170) Hout(idACcy,ng),         &
     &       'Hout(idACcy)',                                            &
     &       'Wave group celerity in the eta space.'
          IF (Hout(idACtp,ng)) WRITE (out,170) Hout(idACtp,ng),         &
     &       'Hout(idACtp)',                                            &
     &       'Wave group peak period.'
#endif
#if defined SOLVE3D && defined T_PASSIVE
          DO itrc=1,NPT
            IF (Hout(idTvar(inert(itrc)),ng)) WRITE (out,180)           &
     &          Hout(idTvar(inert(itrc)),ng), 'Hout(inert)',            &
     &          'Write out inert passive tracer ', itrc,                &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
#endif
#ifdef SOLVE3D
# if defined RUNOFF
          IF (Hout(idRunoff,ng)) WRITE (out,170) Hout(idRunoff,ng),     &
     &       'Hout(idRunoff)',                                          &
     &       'Write out surface fresh water input via runoff.'
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
  || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
          IF (Hout(idPair,ng)) WRITE (out,170) Hout(idPair,ng),         &
     &       'Hout(idPair)',                                            &
     &       'Write out surface air pressure.'
# endif
# if defined BULK_FLUXES || defined ECOSIM
          IF (Hout(idUair,ng)) WRITE (out,170) Hout(idUair,ng),         &
     &       'Hout(idUair)',                                            &
     &       'Write out surface U-wind component.'
          IF (Hout(idVair,ng)) WRITE (out,170) Hout(idVair,ng),         &
     &       'Hout(idVair)',                                            &
     &       'Write out surface V-wind component.'
          IF (Hout(idUairE,ng)) WRITE (out,170) Hout(idUairE,ng),        &
     &       'Hout(idUairE)',                                            &
     &       'Write out surface Eastward U-wind component.'
          IF (Hout(idVairN,ng)) WRITE (out,170) Hout(idVairN,ng),        &
     &       'Hout(idVairN)',                                            &
     &       'Write out surface Northward V-wind component.'
# endif
          IF (Hout(idTsur(itemp),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(itemp),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net heat flux.'
          IF (Hout(idTsur(isalt),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(isalt),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net salt flux.'
# ifdef SHORTWAVE
          IF (Hout(idSrad,ng)) WRITE (out,170) Hout(idSrad,ng),         &
     &       'Hout(idSrad)',                                            &
     &       'Write out shortwave radiation flux.'
# endif
# ifdef BULK_FLUXES
          IF (Hout(idLrad,ng)) WRITE (out,170) Hout(idLrad,ng),         &
     &       'Hout(idLrad)',                                            &
     &       'Write out longwave radiation flux.'
          IF (Hout(idLhea,ng)) WRITE (out,170) Hout(idLhea,ng),         &
     &       'Hout(idLhea)',                                            &
     &       'Write out latent heat flux.'
          IF (Hout(idShea,ng)) WRITE (out,170) Hout(idShea,ng),         &
     &       'Hout(idShea)',                                            &
     &       'Write out sensible heat flux.'
#  ifdef EMINUSP
          IF (Hout(idEmPf,ng)) WRITE (out,170) Hout(idEmPf,ng),         &
     &       'Hout(idEmPf)',                                            &
     &       'Write out E-P flux.'
          IF (Hout(idevap,ng)) WRITE (out,170) Hout(idevap,ng),         &
     &       'Hout(idevap)',                                            &
     &       'Write out evaporation rate.'
          IF (Hout(idrain,ng)) WRITE (out,170) Hout(idrain,ng),         &
     &       'Hout(idrain)',                                            &
     &       'Write out rain rate.'
#  endif
# endif
          IF (Hout(idDano,ng)) WRITE (out,170) Hout(idDano,ng),         &
     &       'Hout(idDano)',                                            &
     &       'Write out density anomaly.'
          IF (Hout(idVvis,ng)) WRITE (out,170) Hout(idVvis,ng),         &
     &       'Hout(idVvis)',                                            &
     &       'Write out vertical viscosity: AKv.'
          IF (Hout(idTdif,ng)) WRITE (out,170) Hout(idTdif,ng),         &
     &       'Hout(idTdif)',                                            &
     &       'Write out vertical diffusion: AKt(itemp).'
          IF (Hout(idSdif,ng)) WRITE (out,170) Hout(idSdif,ng),         &
     &       'Hout(idSdif)',                                            &
     &       'Write out vertical diffusion: AKt(isalt).'
# ifdef LMD_SKPP
          IF (Hout(idHsbl,ng)) WRITE (out,170) Hout(idHsbl,ng),         &
     &       'Hout(idHsbl)',                                            &
     &       'Write out depth of surface boundary layer.'
# endif
# ifdef LMD_BKPP
          IF (Hout(idHbbl,ng)) WRITE (out,170) Hout(idHbbl,ng),         &
     &       'Hout(idHbbl)',                                            &
     &       'Write out depth of bottom boundary layer.'
# endif
# if defined GLS_MIXING || defined MY25_MIXING
          IF (Hout(idMtke,ng)) WRITE (out,170) Hout(idMtke,ng),         &
     &       'Hout(idMtke)',                                            &
     &       'Write out turbulent kinetic energy.'
          IF (Hout(idMtls,ng)) WRITE (out,170) Hout(idMtls,ng),         &
     &       'Hout(idMtls)',                                            &
     &       'Write out turbulent generic length-scale.'
# endif
#endif
#if defined HISTORY2
          WRITE (out,'(1x)')
          IF (Hout2(idFsur,ng)) WRITE (out,170) Hout2(idFsur,ng),       &
     &       'Hout2(idFsur)',                                           &
     &       'Write out averaged free-surface.'
          IF (Hout2(idUbar,ng)) WRITE (out,170) Hout2(idUbar,ng),       &
     &       'Hout2(idUbar)',                                           &
     &       'Write out averaged 2D U-momentum component.'
          IF (Hout2(idVbar,ng)) WRITE (out,170) Hout2(idVbar,ng),       &
     &       'Hout2(idVbar)',                                           &
     &       'Write out averaged 2D V-momentum component.'
          IF (Hout2(idu2dE,ng)) WRITE (out,170) Hout2(idu2dE,ng),       &
     &       'Hout2(idu2dE)',                                           &
     &       'Write out averaged 2D U-eastward  at RHO-points.'
          IF (Hout2(idv2dN,ng)) WRITE (out,170) Hout2(idv2dN,ng),       &
     &       'Hout2(idv2dN)',                                           &
     &       'Write out averaged 2D V-northward at RHO-points.'
# ifdef SOLVE3D
          IF (Hout2(idUvel,ng)) WRITE (out,170) Hout2(idUvel,ng),       &
     &       'Hout2(idUvel)',                                           &
     &       'Write out averaged 3D U-momentum component.'
          IF (Hout2(idVvel,ng)) WRITE (out,170) Hout2(idVvel,ng),       &
     &       'Hout2(idVvel)',                                           &
     &       'Write out averaged 3D V-momentum component.'
          IF (Hout2(idu3dE,ng)) WRITE (out,170) Hout2(idu3dE,ng),       &
     &       'Hout2(idu3dE)',                                           &
     &       'Write out averaged 3D U-eastward  at RHO-points.'
          IF (Hout2(idv3dN,ng)) WRITE (out,170) Hout2(idv3dN,ng),       &
     &       'Hout2(idv3dN)',                                           &
     &       'Write out averaged 3D V-northward at RHO-points.'
          IF (Hout2(idWvel,ng)) WRITE (out,170) Hout2(idWvel,ng),       &
     &       'Hout2(idWvel)',                                           &
     &       'Write out averaged W-momentum component.'
          IF (Hout2(idOvel,ng)) WRITE (out,170) Hout2(idOvel,ng),       &
     &       'Hout2(idOvel)',                                           &
     &       'Write out averaged omega vertical velocity.'
          DO itrc=1,NAT
            IF (Hout2(idTvar(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(idTvar(itrc),ng), 'Hout2(idTvar)',                &
     &          'Write out averaged tracer ', itrc,                     &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          IF (Hout2(idUsms,ng)) WRITE (out,170) Hout2(idUsms,ng),       &
     &       'Hout2(idUsms)',                                           &
     &       'Write out averaged surface U-momentum stress.'
          IF (Hout2(idVsms,ng)) WRITE (out,170) Hout2(idVsms,ng),       &
     &       'Hout2(idVsms)',                                           &
     &       'Write out averaged surface V-momentum stress.'
          IF (Hout2(idUbms,ng)) WRITE (out,170) Hout2(idUbms,ng),       &
     &       'Hout2(idUbms)',                                           &
     &       'Write out averaged bottom U-momentum stress.'
          IF (Hout2(idVbms,ng)) WRITE (out,170) Hout2(idVbms,ng),       &
     &       'Hout2(idVbms)',                                           &
     &       'Write out averaged bottom V-momentum stress.'
# if defined NEARSHORE_MELLOR
          IF (Hout2(idW2xx,ng)) WRITE (out,170) Hout2(idW2xx,ng),       &
     &       'Hout2(idW2xx)',                                           &
     &       'Write out averaged 2D radiation stress, Sxx.'
          IF (Hout2(idW2xy,ng)) WRITE (out,170) Hout2(idW2xy,ng),       &
     &       'Hout2(idW2xy)',                                           &
     &       'Write out averaged 2D radiation stress, Sxy.'
          IF (Hout2(idW2yy,ng)) WRITE (out,170) Hout2(idW2yy,ng),       &
     &       'Hout2(idW2yy)',                                           &
     &       'Write out averaged 2D radiation stress, Syy.'
          IF (Hout2(idU2rs,ng)) WRITE (out,170) Hout2(idU2rs,ng),       &
     &       'Hout2(idU2rs)',                                           &
     &       'Write out averaged total 2D u-radiation stress.'
          IF (Hout2(idV2rs,ng)) WRITE (out,170) Hout2(idV2rs,ng),       &
     &       'Hout2(idV2rs)',                                           &
     &       'Write out averaged total 2D v-radiation stress.'
          IF (Hout2(idU2Sd,ng)) WRITE (out,170) Hout2(idU2Sd,ng),       &
     &       'Hout2(idU2Sd)',                                           &
     &       'Write out averaged 2D u-momentum stokes velocity.'
          IF (Hout2(idV2Sd,ng)) WRITE (out,170) Hout2(idV2Sd,ng),       &
     &       'Hout2(idV2Sd)',                                           &
     &       'Write out averaged 2D v-momentum stokes velocity.'
#  ifdef SOLVE3D
          IF (Hout2(idW3xx,ng)) WRITE (out,170) Hout2(idW3xx,ng),       &
     &       'Hout2(idW3xx)',                                           &
     &       'Write out averaged 3D radiation stress, Sxx.'
          IF (Hout2(idW3xy,ng)) WRITE (out,170) Hout2(idW3xy,ng),       &
     &       'Hout2(idW3xy)',                                           &
     &       'Write out averaged 3D radiation stress, Sxy.'
          IF (Hout2(idW3yy,ng)) WRITE (out,170) Hout2(idW3yy,ng),       &
     &       'Hout2(idW3yy)',                                           &
     &       'Write out averaged 3D radiation stress, Syy.'
          IF (Hout2(idW3zx,ng)) WRITE (out,170) Hout2(idW3zx,ng),       &
     &       'Hout2(idW3zx)',                                           &
     &       'Write out averaged 3D radiation stress, Szx.'
          IF (Hout2(idW3zy,ng)) WRITE (out,170) Hout2(idW3zy,ng),       &
     &       'Hout2(idW3zy)',                                           &
     &       'Write out averaged 3D radiation stress, Szy.'
          IF (Hout2(idU3rs,ng)) WRITE (out,170) Hout2(idU3rs,ng),       &
     &       'Hout2(idU3rs)',                                           &
     &       'Write out averaged total 3D u-radiation stress.'
          IF (Hout2(idV3rs,ng)) WRITE (out,170) Hout2(idV3rs,ng),       &
     &       'Hout2(idV3rs)',                                           &
     &       'Write out averaged total 3D v-radiation stress.'
          IF (Hout2(idU3Sd,ng)) WRITE (out,170) Hout2(idU3Sd,ng),       &
     &       'Hout2(idU3Sd)',                                           &
     &       'Write out averaged 3D u-momentum stokes velocity.'
          IF (Hout2(idV3Sd,ng)) WRITE (out,170) Hout2(idV3Sd,ng),       &
     &       'Hout2(idV3Sd)',                                           &
     &       'Write out averaged 3D v-momentum stokes velocity.'
#  endif
# endif
# if defined SOLVE3D && defined T_PASSIVE
          DO itrc=1,NPT
            IF (Hout2(idTvar(inert(itrc)),ng)) WRITE (out,180)          &
     &          Hout2(idTvar(inert(itrc)),ng), 'Hout2(inert)',          &
     &          'Write out  averaged inert passive tracer ', itrc,      &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
# endif
# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          IF (Hout2(idPair,ng)) WRITE (out,170) Hout2(idPair,ng),       &
     &       'Hout2(idPair)',                                           &
     &       'Write out averaged surface air pressure.'
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          IF (Hout2(idUair,ng)) WRITE (out,170) Hout2(idUair,ng),       &
     &       'Hout2(idUair)',                                           &
     &       'Write out averaged surface U-wind component.'
          IF (Hout2(idVair,ng)) WRITE (out,170) Hout2(idVair,ng),       &
     &       'Hout2(idVair)',                                           &
     &       'Write out averaged surface V-wind component.'
          IF (Hout2(idUairE,ng)) WRITE (out,170) Hout2(idUairE,ng),     &
     &       'Hout2(idUairE)',                                          &
     &       'Write out averaged Eastward surface U-wind component.'
          IF (Hout2(idVairN,ng)) WRITE (out,170) Hout2(idVairN,ng),     &
     &       'Hout2(idVairN)',                                          &
     &       'Write out averaged Northward surface V-wind component.'
#  endif
          IF (Hout2(idTsur(itemp),ng)) WRITE (out,170)                  &
     &        Hout2(idTsur(itemp),ng), 'Hout2(idTsur)',                 &
     &       'Write out averaged surface net heat flux.'
          IF (Hout2(idTsur(isalt),ng)) WRITE (out,170)                  &
     &        Hout2(idTsur(isalt),ng), 'Hout2(idTsur)',                 &
     &       'Write out averaged surface net salt flux.'
#  ifdef SHORTWAVE
          IF (Hout2(idSrad,ng)) WRITE (out,170) Hout2(idSrad,ng),       &
     &       'Hout2(idSrad)',                                           &
     &       'Write out averaged shortwave radiation flux.'
#  endif
#  ifdef BULK_FLUXES
          IF (Hout2(idLrad,ng)) WRITE (out,170) Hout2(idLrad,ng),       &
     &       'Hout2(idLrad)',                                           &
     &       'Write out averaged longwave radiation flux.'
          IF (Hout2(idLhea,ng)) WRITE (out,170) Hout2(idLhea,ng),       &
     &       'Hout2(idLhea)',                                           &
     &       'Write out averaged latent heat flux.'
          IF (Hout2(idShea,ng)) WRITE (out,170) Hout2(idShea,ng),       &
     &       'Hout2(idShea)',                                           &
     &       'Write out averaged sensible heat flux.'
#   ifdef EMINUSP
          IF (Hout2(idevap,ng)) WRITE (out,170) Hout2(idevap,ng),       &
     &       'Hout2(idevap)',                                           &
     &       'Write out averaged evaporation rate.'
          IF (Hout2(idrain,ng)) WRITE (out,170) Hout2(idrain,ng),       &
     &       'Hout2(idrain)',                                           &
     &       'Write out averaged rain rate.'
#   endif
#  endif
          IF (Hout2(idDano,ng)) WRITE (out,170) Hout2(idDano,ng),       &
     &       'Hout2(idDano)',                                           &
     &       'Write out averaged density anomaly.'
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          IF (Hout2(idVvis,ng)) WRITE (out,170) Hout2(idVvis,ng),       &
     &       'Hout2(idVvis)',                                           &
     &       'Write out averaged vertical viscosity: AKv.'
          IF (Hout2(idTdif,ng)) WRITE (out,170) Hout2(idTdif,ng),       &
     &       'Hout2(idTdif)',                                           &
     &       'Write out averaged vertical diffusion: AKt(itemp).'
          IF (Hout2(idSdif,ng)) WRITE (out,170) Hout2(idSdif,ng),       &
     &       'Hout2(idSdif)',                                           &
     &       'Write out averaged vertical diffusion: AKt(isalt).'
#  endif
#  ifdef LMD_SKPP
          IF (Hout2(idHsbl,ng)) WRITE (out,170) Hout2(idHsbl,ng),       &
     &       'Hout2(idHsbl)',                                           &
     &       'Write out averaged depth of surface boundary layer.'
#  endif
#  ifdef LMD_BKPP
          IF (Hout2(idHbbl,ng)) WRITE (out,170) Hout2(idHbbl,ng),       &
     &       'Hout2(idHbbl)',                                           &
     &       'Write out averaged depth of bottom boundary layer.'
#  endif
# endif
          IF (Hout2(id2dRV,ng)) WRITE (out,170) Hout2(id2dRV,ng),       &
     &       'Hout2(id2dRV)',                                           &
     &       'Write out averaged 2D relative vorticity.'
          IF (Hout2(id2dPV,ng)) WRITE (out,170) Hout2(id2dPV,ng),       &
     &       'Hout2(id2dPV)',                                           &
     &       'Write out averaged 2D potential vorticity.'
# ifdef SOLVE3D
          IF (Hout2(id3dRV,ng)) WRITE (out,170) Hout2(id3dRV,ng),       &
     &       'Hout2(id3dRV)',                                           &
     &       'Write out averaged 3D relative vorticity.'
          IF (Hout2(id3dPV,ng)) WRITE (out,170) Hout2(id3dPV,ng),       &
     &       'Hout2(id3dPV)',                                           &
     &       'Write out averaged 3D potential vorticity.'
# endif
          IF (Hout2(idZZav,ng)) WRITE (out,170) Hout2(idZZav,ng),       &
     &       'Hout2(idZZav)',                                           &
     &       'Write out averaged quadratic <zeta*zeta> term.'
          IF (Hout2(idU2av,ng)) WRITE (out,170) Hout2(idU2av,ng),       &
     &       'Hout2(idU2av)',                                           &
     &       'Write out averaged quadratic <ubar*ubar> term.'
          IF (Hout2(idV2av,ng)) WRITE (out,170) Hout2(idV2av,ng),       &
     &       'Hout2(idV2av)',                                           &
     &       'Write out averaged quadratic <vbar*vbar> term.'
# ifdef SOLVE3D
          IF (Hout2(idHUav,ng)) WRITE (out,170) Hout2(idHUav,ng),       &
     &       'Hout2(idHUav)',                                           &
     &       'Write out averaged u-volume flux, Huon.'
          IF (Hout2(idHVav,ng)) WRITE (out,170) Hout2(idHVav,ng),       &
     &       'Hout2(idHVav)',                                           &
     &       'Write out averaged v-volume flux, Hvom.'
          IF (Hout2(idUUav,ng)) WRITE (out,170) Hout2(idUUav,ng),       &
     &       'Hout2(idUUav)',                                           &
     &       'Write out averaged quadratic <u*u> term.'
          IF (Hout2(idUVav,ng)) WRITE (out,170) Hout2(idUVav,ng),       &
     &       'Hout2(idUVav)',                                           &
     &       'Write out averaged quadratic <u*v> term.'
          IF (Hout2(idVVav,ng)) WRITE (out,170) Hout2(idVVav,ng),       &
     &       'Hout2(idVVav)',                                           &
     &       'Write out averaged quadratic <v*v> term.'
          DO itrc=1,NAT
            IF (Hout2(idTTav(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(idTTav(itrc),ng), 'Hout2(idTTav)',                &
     &          'Write out averaged <t*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Hout2(idUTav(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(idUTav(itrc),ng), 'Hout2(idUTav)',                &
     &          'Write out averaged <u*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Hout2(idVTav(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(idVTav(itrc),ng), 'Hout2(idVTav)',                &
     &          'Write out averaged <v*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Hout2(iHUTav(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(iHUTav(itrc),ng), 'Hout2(iHUTav)',                &
     &          'Write out averaged <Huon*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Hout2(iHVTav(itrc),ng)) WRITE (out,180)                 &
     &          Hout2(iHVTav(itrc),ng), 'Hout2(iHVTav)',                &
     &          'Write out averaged <Hvom*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          WRITE (out,'(1x)')
          IF (Aout(idFsur,ng)) WRITE (out,170) Aout(idFsur,ng),         &
     &       'Aout(idFsur)',                                            &
     &       'Write out averaged free-surface.'
          IF (Aout(idUbar,ng)) WRITE (out,170) Aout(idUbar,ng),         &
     &       'Aout(idUbar)',                                            &
     &       'Write out averaged 2D U-momentum component.'
          IF (Aout(idVbar,ng)) WRITE (out,170) Aout(idVbar,ng),         &
     &       'Aout(idVbar)',                                            &
     &       'Write out averaged 2D V-momentum component.'
          IF (Aout(idu2dE,ng)) WRITE (out,170) Aout(idu2dE,ng),         &
     &       'Aout(idu2dE)',                                            &
     &       'Write out averaged 2D U-eastward  at RHO-points.'
          IF (Aout(idv2dN,ng)) WRITE (out,170) Aout(idv2dN,ng),         &
     &       'Aout(idv2dN)',                                            &
     &       'Write out averaged 2D V-northward at RHO-points.'
# ifdef SOLVE3D
          IF (Aout(idUvel,ng)) WRITE (out,170) Aout(idUvel,ng),         &
     &       'Aout(idUvel)',                                            &
     &       'Write out averaged 3D U-momentum component.'
          IF (Aout(idVvel,ng)) WRITE (out,170) Aout(idVvel,ng),         &
     &       'Aout(idVvel)',                                            &
     &       'Write out averaged 3D V-momentum component.'
          IF (Aout(idu3dE,ng)) WRITE (out,170) Aout(idu3dE,ng),         &
     &       'Aout(idu3dE)',                                            &
     &       'Write out averaged 3D U-eastward  at RHO-points.'
          IF (Aout(idv3dN,ng)) WRITE (out,170) Aout(idv3dN,ng),         &
     &       'Aout(idv3dN)',                                            &
     &       'Write out averaged 3D V-northward at RHO-points.'
          IF (Aout(idWvel,ng)) WRITE (out,170) Aout(idWvel,ng),         &
     &       'Aout(idWvel)',                                            &
     &       'Write out averaged W-momentum component.'
          IF (Aout(idOvel,ng)) WRITE (out,170) Aout(idOvel,ng),         &
     &       'Aout(idOvel)',                                            &
     &       'Write out averaged omega vertical velocity.'
          DO itrc=1,NAT
            IF (Aout(idTvar(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idTvar(itrc),ng), 'Aout(idTvar)',                  &
     &          'Write out averaged tracer ', itrc,                     &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          IF (Aout(idUsms,ng)) WRITE (out,170) Aout(idUsms,ng),         &
     &       'Aout(idUsms)',                                            &
     &       'Write out averaged surface U-momentum stress.'
          IF (Aout(idVsms,ng)) WRITE (out,170) Aout(idVsms,ng),         &
     &       'Aout(idVsms)',                                            &
     &       'Write out averaged surface V-momentum stress.'
          IF (Aout(idUbms,ng)) WRITE (out,170) Aout(idUbms,ng),         &
     &       'Aout(idUbms)',                                            &
     &       'Write out averaged bottom U-momentum stress.'
          IF (Aout(idVbms,ng)) WRITE (out,170) Aout(idVbms,ng),         &
     &       'Aout(idVbms)',                                            &
     &       'Write out averaged bottom V-momentum stress.'
# ifdef BBL_MODEL
          IF (Aout(idUbrs,ng)) WRITE (out,170) Aout(idUbrs,ng),         &
     &       'Aout(idUbrs)',                                            &
     &       'Write out averaged bottom U-current stress.'
          IF (Aout(idVbrs,ng)) WRITE (out,170) Aout(idVbrs,ng),         &
     &       'Aout(idVbrs)',                                            &
     &       'Write out averaged bottom V-current stress.'
          IF (Aout(idUbws,ng)) WRITE (out,170) Aout(idUbws,ng),         &
     &       'Aout(idUbws)',                                            &
     &       'Write out averaged wind-induced, bottom U-wave stress.'
          IF (Aout(idVbws,ng)) WRITE (out,170) Aout(idVbws,ng),         &
     &       'Aout(idVbws)',                                            &
     &       'Write out averaged wind-induced, bottom V-wave stress.'
          IF (Aout(idUbcs,ng)) WRITE (out,170) Aout(idUbcs,ng),         &
     &       'Aout(idUbcs)',                                            &
     &       'Write out averaged max wind + current, bottom U-wave stress.'
          IF (Aout(idVbcs,ng)) WRITE (out,170) Aout(idVbcs,ng),         &
     &       'Aout(idVbcs)',                                            &
     &       'Write out averaged max wind + current, bottom V-wave stress.'
          IF (Aout(idUVwc,ng)) WRITE (out,170) Aout(idUVwc,ng),         &
     &       'Aout(idUVwc)',                                            &
     &       'Write out averaged max wind + current, bottom UV-wave stress.'
          IF (Aout(idUbot,ng)) WRITE (out,170) Aout(idUbot,ng),         &
     &       'Aout(idUbot)',                                            &
     &       'Write out averaged bed wave orbital U-velocity.'
          IF (Aout(idVbot,ng)) WRITE (out,170) Aout(idVbot,ng),         &
     &       'Aout(idVbot)',                                            &
     &       'Write out averaged bed wave orbital V-velocity.'
          IF (Aout(idUbur,ng)) WRITE (out,170) Aout(idUbur,ng),         &
     &       'Aout(idUbur)',                                            &
     &       'Write out averaged bottom U-momentum above bed.'
          IF (Aout(idVbvr,ng)) WRITE (out,170) Aout(idVbvr,ng),         &
     &       'Aout(idVbvr)',                                            &
     &       'Write out averaged bottom V-momentum above bed.'
# endif
# if defined WEC_MELLOR
          IF (Aout(idW2xx,ng)) WRITE (out,170) Aout(idW2xx,ng),         &
     &       'Aout(idW2xx)',                                            &
     &       'Write out averaged 2D radiation stress, Sxx.'
          IF (Aout(idW2xy,ng)) WRITE (out,170) Aout(idW2xy,ng),         &
     &       'Aout(idW2xy)',                                            &
     &       'Write out averaged 2D radiation stress, Sxy.'
          IF (Aout(idW2yy,ng)) WRITE (out,170) Aout(idW2yy,ng),         &
     &       'Aout(idW2yy)',                                            &
     &       'Write out averaged 2D radiation stress, Syy.'
# endif
# if defined WEC
          IF (Aout(idU2rs,ng)) WRITE (out,170) Aout(idU2rs,ng),         &
     &       'Aout(idU2rs)',                                            &
     &       'Write out averaged total 2D u-wec stress.'
          IF (Aout(idV2rs,ng)) WRITE (out,170) Aout(idV2rs,ng),         &
     &       'Aout(idV2rs)',                                            &
     &       'Write out averaged total 2D v-wec stress.'
          IF (Aout(idU2Sd,ng)) WRITE (out,170) Aout(idU2Sd,ng),         &
     &       'Aout(idU2Sd)',                                            &
     &       'Write out averaged 2D u-momentum stokes velocity.'
          IF (Aout(idV2Sd,ng)) WRITE (out,170) Aout(idV2Sd,ng),         &
     &       'Aout(idV2Sd)',                                            &
     &       'Write out averaged 2D v-momentum stokes velocity.'
# endif
# if defined WEC_MELLOR
#  ifdef SOLVE3D
          IF (Aout(idW3xx,ng)) WRITE (out,170) Aout(idW3xx,ng),         &
     &       'Aout(idW3xx)',                                            &
     &       'Write out averaged 3D radiation stress, Sxx.'
          IF (Aout(idW3xy,ng)) WRITE (out,170) Aout(idW3xy,ng),         &
     &       'Aout(idW3xy)',                                            &
     &       'Write out averaged 3D radiation stress, Sxy.'
          IF (Aout(idW3yy,ng)) WRITE (out,170) Aout(idW3yy,ng),         &
     &       'Aout(idW3yy)',                                            &
     &       'Write out averaged 3D radiation stress, Syy.'
          IF (Aout(idW3zx,ng)) WRITE (out,170) Aout(idW3zx,ng),         &
     &       'Aout(idW3zx)',                                            &
     &       'Write out averaged 3D radiation stress, Szx.'
          IF (Aout(idW3zy,ng)) WRITE (out,170) Aout(idW3zy,ng),         &
     &       'Aout(idW3zy)',                                            &
     &       'Write out averaged 3D radiation stress, Szy.'
#  endif
# endif
# if defined WEC
#  ifdef SOLVE3D
          IF (Aout(idU3rs,ng)) WRITE (out,170) Aout(idU3rs,ng),         &
     &       'Aout(idU3rs)',                                            &
     &       'Write out averaged total 3D u-wec stress.'
          IF (Aout(idV3rs,ng)) WRITE (out,170) Aout(idV3rs,ng),         &
     &       'Aout(idV3rs)',                                            &
     &       'Write out averaged total 3D v-wec stress.'
          IF (Aout(idU3Sd,ng)) WRITE (out,170) Aout(idU3Sd,ng),         &
     &       'Aout(idU3Sd)',                                            &
     &       'Write out averaged 3D u-momentum stokes velocity.'
          IF (Aout(idV3Sd,ng)) WRITE (out,170) Aout(idV3Sd,ng),         &
     &       'Aout(idV3Sd)',                                            &
     &       'Write out averaged 3D v-momentum stokes velocity.'
          IF (Aout(idW3Sd,ng)) WRITE (out,170) Aout(idW3Sd,ng),         &
     &       'Aout(idW3Sd)',                                            &
     &       'Write out averaged 3D Omega-momentum stokes velocity.'
          IF (Aout(idW3St,ng)) WRITE (out,170) Aout(idW3St,ng),         &
     &       'Aout(idW3St)',                                            &
     &       'Write out averaged 3D w-momentum stokes velocity.'
#  endif
          IF (Aout(idWztw,ng)) WRITE (out,170) Aout(idWztw,ng),         &
     &       'Aout(idWztw)',                                            &
     &       'Write out averaged total 3D u-wec stress.'


# endif
# if defined SOLVE3D && defined T_PASSIVE
          DO itrc=1,NPT
            IF (Aout(idTvar(inert(itrc)),ng)) WRITE (out,180)           &
     &          Aout(idTvar(inert(itrc)),ng), 'Aout(inert)',            &
     &          'Write out  averaged inert passive tracer ', itrc,      &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
# endif
# ifdef SOLVE3D
#  if defined RUNOFF
          IF (Aout(idRunoff,ng)) WRITE (out,170) Aout(idRunoff,ng),     &
     &       'Aout(idRunoff)',                                          &
     &       'Write out averaged surface fresh water input via runoff.'
#  endif
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS \
   || defined ATM2OCN_FLUXES || defined CCSM_FLUXES2D
          IF (Aout(idPair,ng)) WRITE (out,170) Aout(idPair,ng),         &
     &       'Aout(idPair)',                                            &
     &       'Write out averaged surface air pressure.'
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          IF (Aout(idUair,ng)) WRITE (out,170) Aout(idUair,ng),         &
     &       'Aout(idUair)',                                            &
     &       'Write out averaged surface U-wind component.'
          IF (Aout(idVair,ng)) WRITE (out,170) Aout(idVair,ng),         &
     &       'Aout(idVair)',                                            &
     &       'Write out averaged surface V-wind component.'
          IF (Aout(idUairE,ng)) WRITE (out,170) Aout(idUairE,ng),        &
     &       'Aout(idUairE)',                                            &
     &       'Write out averaged Eastward surface U-wind component.'
          IF (Aout(idVairN,ng)) WRITE (out,170) Aout(idVairN,ng),        &
     &       'Aout(idVairN)',                                            &
     &       'Write out averaged Northward surface V-wind component.'
#  endif
          IF (Aout(idTsur(itemp),ng)) WRITE (out,170)                   &
     &        Aout(idTsur(itemp),ng), 'Aout(idTsur)',                   &
     &       'Write out averaged surface net heat flux.'
          IF (Aout(idTsur(isalt),ng)) WRITE (out,170)                   &
     &        Aout(idTsur(isalt),ng), 'Aout(idTsur)',                   &
     &       'Write out averaged surface net salt flux.'
#  ifdef SHORTWAVE
          IF (Aout(idSrad,ng)) WRITE (out,170) Aout(idSrad,ng),         &
     &       'Aout(idSrad)',                                            &
     &       'Write out averaged shortwave radiation flux.'
#  endif
#  ifdef BULK_FLUXES
          IF (Aout(idLrad,ng)) WRITE (out,170) Aout(idLrad,ng),         &
     &       'Aout(idLrad)',                                            &
     &       'Write out averaged longwave radiation flux.'
          IF (Aout(idLhea,ng)) WRITE (out,170) Aout(idLhea,ng),         &
     &       'Aout(idLhea)',                                            &
     &       'Write out averaged latent heat flux.'
          IF (Aout(idShea,ng)) WRITE (out,170) Aout(idShea,ng),         &
     &       'Aout(idShea)',                                            &
     &       'Write out averaged sensible heat flux.'
#   ifdef EMINUSP
          IF (Aout(idevap,ng)) WRITE (out,170) Aout(idevap,ng),         &
     &       'Aout(idevap)',                                            &
     &       'Write out averaged evaporation rate.'
          IF (Aout(idrain,ng)) WRITE (out,170) Aout(idrain,ng),         &
     &       'Aout(idrain)',                                            &
     &       'Write out averaged rain rate.'
#   endif
#  endif
          IF (Aout(idDano,ng)) WRITE (out,170) Aout(idDano,ng),         &
     &       'Aout(idDano)',                                            &
     &       'Write out averaged density anomaly.'
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          IF (Aout(idVvis,ng)) WRITE (out,170) Aout(idVvis,ng),         &
     &       'Aout(idVvis)',                                            &
     &       'Write out averaged vertical viscosity: AKv.'
          IF (Aout(idTdif,ng)) WRITE (out,170) Aout(idTdif,ng),         &
     &       'Aout(idTdif)',                                            &
     &       'Write out averaged vertical diffusion: AKt(itemp).'
          IF (Aout(idSdif,ng)) WRITE (out,170) Aout(idSdif,ng),         &
     &       'Aout(idSdif)',                                            &
     &       'Write out averaged vertical diffusion: AKt(isalt).'
#  endif
#  ifdef LMD_SKPP
          IF (Aout(idHsbl,ng)) WRITE (out,170) Aout(idHsbl,ng),         &
     &       'Aout(idHsbl)',                                            &
     &       'Write out averaged depth of surface boundary layer.'
#  endif
#  ifdef LMD_BKPP
          IF (Aout(idHbbl,ng)) WRITE (out,170) Aout(idHbbl,ng),         &
     &       'Aout(idHbbl)',                                            &
     &       'Write out averaged depth of bottom boundary layer.'
#  endif
# endif
          IF (Aout(id2dRV,ng)) WRITE (out,170) Aout(id2dRV,ng),         &
     &       'Aout(id2dRV)',                                            &
     &       'Write out averaged 2D relative vorticity.'
          IF (Aout(id2dPV,ng)) WRITE (out,170) Aout(id2dPV,ng),         &
     &       'Aout(id2dPV)',                                            &
     &       'Write out averaged 2D potential vorticity.'
# ifdef SOLVE3D
          IF (Aout(id3dRV,ng)) WRITE (out,170) Aout(id3dRV,ng),         &
     &       'Aout(id3dRV)',                                            &
     &       'Write out averaged 3D relative vorticity.'
          IF (Aout(id3dPV,ng)) WRITE (out,170) Aout(id3dPV,ng),         &
     &       'Aout(id3dPV)',                                            &
     &       'Write out averaged 3D potential vorticity.'
# endif
          IF (Aout(idZZav,ng)) WRITE (out,170) Aout(idZZav,ng),         &
     &       'Aout(idZZav)',                                            &
     &       'Write out averaged quadratic <zeta*zeta> term.'
          IF (Aout(idU2av,ng)) WRITE (out,170) Aout(idU2av,ng),         &
     &       'Aout(idU2av)',                                            &
     &       'Write out averaged quadratic <ubar*ubar> term.'
          IF (Aout(idV2av,ng)) WRITE (out,170) Aout(idV2av,ng),         &
     &       'Aout(idV2av)',                                            &
     &       'Write out averaged quadratic <vbar*vbar> term.'
# ifdef SOLVE3D
          IF (Aout(idHUav,ng)) WRITE (out,170) Aout(idHUav,ng),         &
     &       'Aout(idHUav)',                                            &
     &       'Write out averaged u-volume flux, Huon.'
          IF (Aout(idHVav,ng)) WRITE (out,170) Aout(idHVav,ng),         &
     &       'Aout(idHVav)',                                            &
     &       'Write out averaged v-volume flux, Hvom.'
          IF (Aout(idUUav,ng)) WRITE (out,170) Aout(idUUav,ng),         &
     &       'Aout(idUUav)',                                            &
     &       'Write out averaged quadratic <u*u> term.'
          IF (Aout(idUVav,ng)) WRITE (out,170) Aout(idUVav,ng),         &
     &       'Aout(idUVav)',                                            &
     &       'Write out averaged quadratic <u*v> term.'
          IF (Aout(idVVav,ng)) WRITE (out,170) Aout(idVVav,ng),         &
     &       'Aout(idVVav)',                                            &
     &       'Write out averaged quadratic <v*v> term.'
          DO itrc=1,NAT+NPT
            IF (Aout(idTTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idTTav(itrc),ng), 'Aout(idTTav)',                  &
     &          'Write out averaged <t*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT+NPT
            IF (Aout(idUTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idUTav(itrc),ng), 'Aout(idUTav)',                  &
     &          'Write out averaged <u*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT+NPT
            IF (Aout(idVTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idVTav(itrc),ng), 'Aout(idVTav)',                  &
     &          'Write out averaged <v*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT+NPT
            IF (Aout(iHUTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(iHUTav(itrc),ng), 'Aout(iHUTav)',                  &
     &          'Write out averaged <Huon*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT+NPT
            IF (Aout(iHVTav(itrc),ng)) WRITE (out,180)                  &
     &          Aout(iHVTav(itrc),ng), 'Aout(iHVTav)',                  &
     &          'Write out averaged <Hvom*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          WRITE (out,'(1x)')
          IF (Aout(idFsuD,ng)) WRITE (out,170) Aout(idFsuD,ng),         &
     &       'Aout(idFsuD)',                                            &
     &       'Write out detided free-surface.'
          IF (Aout(idu2dD,ng)) WRITE (out,170) Aout(idu2dD,ng),         &
     &       'Aout(idu2dD)',                                            &
     &       'Write out detided 2D U-velocity.'
          IF (Aout(idv2dD,ng)) WRITE (out,170) Aout(idv2dD,ng),         &
     &       'Aout(idv2dD)',                                            &
     &       'Write out detided 2D V-velocity.'
#  ifdef SOLVE3D
          IF (Aout(idu3dD,ng)) WRITE (out,170) Aout(idu3dD,ng),         &
     &       'Aout(idu3dD)',                                            &
     &       'Write out detided 3D U-velocity.'
          IF (Aout(idv3dD,ng)) WRITE (out,170) Aout(idv3dD,ng),         &
     &       'Aout(idv3dD)',                                            &
     &       'Write out detided 3D V-velocity.'
          DO itrc=1,NAT
            IF (Aout(idTrcD(itrc),ng)) WRITE (out,180)                  &
     &          Aout(idTrcD(itrc),ng), 'Aout(idTrcD)',                  &
     &          'Write out detided tracer ', itrc,                      &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
#  endif
# endif
#endif
#if defined AVERAGES2
          WRITE (out,'(1x)')
          IF (Aout2(idFsur,ng)) WRITE (out,170) Aout2(idFsur,ng),       &
     &       'Aout2(idFsur)',                                           &
     &       'Write out averaged free-surface.'
          IF (Aout2(idUbar,ng)) WRITE (out,170) Aout2(idUbar,ng),       &
     &       'Aout2(idUbar)',                                           &
     &       'Write out averaged 2D U-momentum component.'
          IF (Aout2(idVbar,ng)) WRITE (out,170) Aout2(idVbar,ng),       &
     &       'Aout2(idVbar)',                                           &
     &       'Write out averaged 2D V-momentum component.'
          IF (Aout2(idu2dE,ng)) WRITE (out,170) Aout2(idu2dE,ng),       &
     &       'Aout2(idu2dE)',                                           &
     &       'Write out averaged 2D U-eastward  at RHO-points.'
          IF (Aout2(idv2dN,ng)) WRITE (out,170) Aout2(idv2dN,ng),       &
     &       'Aout2(idv2dN)',                                           &
     &       'Write out averaged 2D V-northward at RHO-points.'
# ifdef SOLVE3D
          IF (Aout2(idUvel,ng)) WRITE (out,170) Aout2(idUvel,ng),       &
     &       'Aout2(idUvel)',                                           &
     &       'Write out averaged 3D U-momentum component.'
          IF (Aout2(idVvel,ng)) WRITE (out,170) Aout2(idVvel,ng),       &
     &       'Aout2(idVvel)',                                           &
     &       'Write out averaged 3D V-momentum component.'
          IF (Aout2(idu3dE,ng)) WRITE (out,170) Aout2(idu3dE,ng),       &
     &       'Aout2(idu3dE)',                                           &
     &       'Write out averaged 3D U-eastward  at RHO-points.'
          IF (Aout2(idv3dN,ng)) WRITE (out,170) Aout2(idv3dN,ng),       &
     &       'Aout2(idv3dN)',                                           &
     &       'Write out averaged 3D V-northward at RHO-points.'
          IF (Aout2(idWvel,ng)) WRITE (out,170) Aout2(idWvel,ng),       &
     &       'Aout2(idWvel)',                                           &
     &       'Write out averaged W-momentum component.'
          IF (Aout2(idOvel,ng)) WRITE (out,170) Aout2(idOvel,ng),       &
     &       'Aout2(idOvel)',                                           &
     &       'Write out averaged omega vertical velocity.'
          DO itrc=1,NAT
            IF (Aout2(idTvar(itrc),ng)) WRITE (out,180)                 &
     &          Aout2(idTvar(itrc),ng), 'Aout2(idTvar)',                &
     &          'Write out averaged tracer ', itrc,                     &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          IF (Aout2(idUsms,ng)) WRITE (out,170) Aout2(idUsms,ng),       &
     &       'Aout2(idUsms)',                                           &
     &       'Write out averaged surface U-momentum stress.'
          IF (Aout2(idVsms,ng)) WRITE (out,170) Aout2(idVsms,ng),       &
     &       'Aout2(idVsms)',                                           &
     &       'Write out averaged surface V-momentum stress.'
          IF (Aout2(idUbms,ng)) WRITE (out,170) Aout2(idUbms,ng),       &
     &       'Aout2(idUbms)',                                           &
     &       'Write out averaged bottom U-momentum stress.'
          IF (Aout2(idVbms,ng)) WRITE (out,170) Aout2(idVbms,ng),       &
     &       'Aout2(idVbms)',                                           &
     &       'Write out averaged bottom V-momentum stress.'
# if defined NEARSHORE_MELLOR
          IF (Aout2(idW2xx,ng)) WRITE (out,170) Aout2(idW2xx,ng),       &
     &       'Aout2(idW2xx)',                                           &
     &       'Write out averaged 2D radiation stress, Sxx.'
          IF (Aout2(idW2xy,ng)) WRITE (out,170) Aout2(idW2xy,ng),       &
     &       'Aout2(idW2xy)',                                           &
     &       'Write out averaged 2D radiation stress, Sxy.'
          IF (Aout2(idW2yy,ng)) WRITE (out,170) Aout2(idW2yy,ng),       &
     &       'Aout2(idW2yy)',                                           &
     &       'Write out averaged 2D radiation stress, Syy.'
          IF (Aout2(idU2rs,ng)) WRITE (out,170) Aout2(idU2rs,ng),       &
     &       'Aout2(idU2rs)',                                           &
     &       'Write out averaged total 2D u-radiation stress.'
          IF (Aout2(idV2rs,ng)) WRITE (out,170) Aout2(idV2rs,ng),       &
     &       'Aout2(idV2rs)',                                           &
     &       'Write out averaged total 2D v-radiation stress.'
          IF (Aout2(idU2Sd,ng)) WRITE (out,170) Aout2(idU2Sd,ng),       &
     &       'Aout2(idU2Sd)',                                           &
     &       'Write out averaged 2D u-momentum stokes velocity.'
          IF (Aout2(idV2Sd,ng)) WRITE (out,170) Aout2(idV2Sd,ng),       &
     &       'Aout2(idV2Sd)',                                           &
     &       'Write out averaged 2D v-momentum stokes velocity.'
#  ifdef SOLVE3D
          IF (Aout2(idW3xx,ng)) WRITE (out,170) Aout2(idW3xx,ng),       &
     &       'Aout2(idW3xx)',                                           &
     &       'Write out averaged 3D radiation stress, Sxx.'
          IF (Aout2(idW3xy,ng)) WRITE (out,170) Aout2(idW3xy,ng),       &
     &       'Aout2(idW3xy)',                                           &
     &       'Write out averaged 3D radiation stress, Sxy.'
          IF (Aout2(idW3yy,ng)) WRITE (out,170) Aout2(idW3yy,ng),       &
     &       'Aout2(idW3yy)',                                           &
     &       'Write out averaged 3D radiation stress, Syy.'
          IF (Aout2(idW3zx,ng)) WRITE (out,170) Aout2(idW3zx,ng),       &
     &       'Aout2(idW3zx)',                                           &
     &       'Write out averaged 3D radiation stress, Szx.'
          IF (Aout2(idW3zy,ng)) WRITE (out,170) Aout2(idW3zy,ng),       &
     &       'Aout2(idW3zy)',                                           &
     &       'Write out averaged 3D radiation stress, Szy.'
          IF (Aout2(idU3rs,ng)) WRITE (out,170) Aout2(idU3rs,ng),       &
     &       'Aout2(idU3rs)',                                           &
     &       'Write out averaged total 3D u-radiation stress.'
          IF (Aout2(idV3rs,ng)) WRITE (out,170) Aout2(idV3rs,ng),       &
     &       'Aout2(idV3rs)',                                           &
     &       'Write out averaged total 3D v-radiation stress.'
          IF (Aout2(idU3Sd,ng)) WRITE (out,170) Aout2(idU3Sd,ng),       &
     &       'Aout2(idU3Sd)',                                           &
     &       'Write out averaged 3D u-momentum stokes velocity.'
          IF (Aout2(idV3Sd,ng)) WRITE (out,170) Aout2(idV3Sd,ng),       &
     &       'Aout2(idV3Sd)',                                           &
     &       'Write out averaged 3D v-momentum stokes velocity.'
#  endif
# endif
# if defined SOLVE3D && defined T_PASSIVE
          DO itrc=1,NPT
            IF (Aout2(idTvar(inert(itrc)),ng)) WRITE (out,180)          &
     &          Aout2(idTvar(inert(itrc)),ng), 'Aout2(inert)',          &
     &          'Write out  averaged inert passive tracer ', itrc,      &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
# endif
# ifdef SOLVE3D
#  if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
          IF (Aout2(idPair,ng)) WRITE (out,170) Aout2(idPair,ng),       &
     &       'Aout2(idPair)',                                           &
     &       'Write out averaged surface air pressure.'
#  endif
#  if defined BULK_FLUXES || defined ECOSIM
          IF (Aout2(idUair,ng)) WRITE (out,170) Aout2(idUair,ng),       &
     &       'Aout2(idUair)',                                           &
     &       'Write out averaged surface U-wind component.'
          IF (Aout2(idVair,ng)) WRITE (out,170) Aout2(idVair,ng),       &
     &       'Aout2(idVair)',                                           &
     &       'Write out averaged surface V-wind component.'
          IF (Aout2(idUairE,ng)) WRITE (out,170) Aout2(idUairE,ng),     &
     &       'Aout2(idUairE)',                                          &
     &       'Write out averaged Eastward surface U-wind component.'
          IF (Aout2(idVairN,ng)) WRITE (out,170) Aout2(idVairN,ng),     &
     &       'Aout2(idVairN)',                                          &
     &       'Write out averaged Northward surface V-wind component.'
#  endif
          IF (Aout2(idTsur(itemp),ng)) WRITE (out,170)                  &
     &        Aout2(idTsur(itemp),ng), 'Aout2(idTsur)',                 &
     &       'Write out averaged surface net heat flux.'
          IF (Aout2(idTsur(isalt),ng)) WRITE (out,170)                  &
     &        Aout2(idTsur(isalt),ng), 'Aout2(idTsur)',                 &
     &       'Write out averaged surface net salt flux.'
#  ifdef SHORTWAVE
          IF (Aout2(idSrad,ng)) WRITE (out,170) Aout2(idSrad,ng),       &
     &       'Aout2(idSrad)',                                           &
     &       'Write out averaged shortwave radiation flux.'
#  endif
#  ifdef BULK_FLUXES
          IF (Aout2(idLrad,ng)) WRITE (out,170) Aout2(idLrad,ng),       &
     &       'Aout2(idLrad)',                                           &
     &       'Write out averaged longwave radiation flux.'
          IF (Aout2(idLhea,ng)) WRITE (out,170) Aout2(idLhea,ng),       &
     &       'Aout2(idLhea)',                                           &
     &       'Write out averaged latent heat flux.'
          IF (Aout2(idShea,ng)) WRITE (out,170) Aout2(idShea,ng),       &
     &       'Aout2(idShea)',                                           &
     &       'Write out averaged sensible heat flux.'
#   ifdef EMINUSP
          IF (Aout2(idevap,ng)) WRITE (out,170) Aout2(idevap,ng),       &
     &       'Aout2(idevap)',                                           &
     &       'Write out averaged evaporation rate.'
          IF (Aout2(idrain,ng)) WRITE (out,170) Aout2(idrain,ng),       &
     &       'Aout2(idrain)',                                           &
     &       'Write out averaged rain rate.'
#   endif
#  endif
          IF (Aout2(idDano,ng)) WRITE (out,170) Aout2(idDano,ng),       &
     &       'Aout2(idDano)',                                           &
     &       'Write out averaged density anomaly.'
#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
          IF (Aout2(idVvis,ng)) WRITE (out,170) Aout2(idVvis,ng),       &
     &       'Aout2(idVvis)',                                           &
     &       'Write out averaged vertical viscosity: AKv.'
          IF (Aout2(idTdif,ng)) WRITE (out,170) Aout2(idTdif,ng),       &
     &       'Aout2(idTdif)',                                           &
     &       'Write out averaged vertical diffusion: AKt(itemp).'
          IF (Aout2(idSdif,ng)) WRITE (out,170) Aout2(idSdif,ng),       &
     &       'Aout2(idSdif)',                                           &
     &       'Write out averaged vertical diffusion: AKt(isalt).'
#  endif
#  ifdef LMD_SKPP
          IF (Aout2(idHsbl,ng)) WRITE (out,170) Aout2(idHsbl,ng),       &
     &       'Aout2(idHsbl)',                                           &
     &       'Write out averaged depth of surface boundary layer.'
#  endif
#  ifdef LMD_BKPP
          IF (Aout2(idHbbl,ng)) WRITE (out,170) Aout2(idHbbl,ng),       &
     &       'Aout2(idHbbl)',                                           &
     &       'Write out averaged depth of bottom boundary layer.'
#  endif
# endif
          IF (Aout2(id2dRV,ng)) WRITE (out,170) Aout2(id2dRV,ng),       &
     &       'Aout2(id2dRV)',                                           &
     &       'Write out averaged 2D relative vorticity.'
          IF (Aout2(id2dPV,ng)) WRITE (out,170) Aout2(id2dPV,ng),       &
     &       'Aout2(id2dPV)',                                           &
     &       'Write out averaged 2D potential vorticity.'
# ifdef SOLVE3D
          IF (Aout2(id3dRV,ng)) WRITE (out,170) Aout2(id3dRV,ng),       &
     &       'Aout2(id3dRV)',                                           &
     &       'Write out averaged 3D relative vorticity.'
          IF (Aout2(id3dPV,ng)) WRITE (out,170) Aout2(id3dPV,ng),       &
     &       'Aout2(id3dPV)',                                           &
     &       'Write out averaged 3D potential vorticity.'
# endif
          IF (Aout2(idZZav,ng)) WRITE (out,170) Aout2(idZZav,ng),       &
     &       'Aout2(idZZav)',                                           &
     &       'Write out averaged quadratic <zeta*zeta> term.'
          IF (Aout2(idU2av,ng)) WRITE (out,170) Aout2(idU2av,ng),       &
     &       'Aout2(idU2av)',                                           &
     &       'Write out averaged quadratic <ubar*ubar> term.'
          IF (Aout2(idV2av,ng)) WRITE (out,170) Aout2(idV2av,ng),       &
     &       'Aout2(idV2av)',                                           &
     &       'Write out averaged quadratic <vbar*vbar> term.'
# ifdef SOLVE3D
          IF (Aout2(idHUav,ng)) WRITE (out,170) Aout2(idHUav,ng),       &
     &       'Aout2(idHUav)',                                           &
     &       'Write out averaged u-volume flux, Huon.'
          IF (Aout2(idHVav,ng)) WRITE (out,170) Aout2(idHVav,ng),       &
     &       'Aout2(idHVav)',                                           &
     &       'Write out averaged v-volume flux, Hvom.'
          IF (Aout2(idUUav,ng)) WRITE (out,170) Aout2(idUUav,ng),       &
     &       'Aout2(idUUav)',                                           &
     &       'Write out averaged quadratic <u*u> term.'
          IF (Aout2(idUVav,ng)) WRITE (out,170) Aout2(idUVav,ng),       &
     &       'Aout2(idUVav)',                                           &
     &       'Write out averaged quadratic <u*v> term.'
          IF (Aout2(idVVav,ng)) WRITE (out,170) Aout2(idVVav,ng),       &
     &       'Aout2(idVVav)',                                           &
     &       'Write out averaged quadratic <v*v> term.'
          DO itrc=1,NAT
            IF (Aout2(idTTav(itrc),ng)) WRITE (out,180)                 &
     &          Aout2(idTTav(itrc),ng), 'Aout2(idTTav)',                &
     &          'Write out averaged <t*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout2(idUTav(itrc),ng)) WRITE (out,180)                 &
     &          Aout2(idUTav(itrc),ng), 'Aout2(idUTav)',                &
     &          'Write out averaged <u*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout2(idVTav(itrc),ng)) WRITE (out,180)                 &
     &          Aout2(idVTav(itrc),ng), 'Aout2(idVTav)',                &
     &          'Write out averaged <v*t> for tracer ', itrc,           &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout2(iHUTav(itrc),ng)) WRITE (out,180)                 &
     &          Aout2(iHUTav(itrc),ng), 'Aout2(iHUTav)',                &
     &          'Write out averaged <Huon*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
          DO itrc=1,NAT
            IF (Aout2(iHVTav(itrc),ng)) WRITE (out,180)                 &
     &          Aout2(iHVTav(itrc),ng), 'Aout2(iHVTav)',                &
     &          'Write out averaged <Hvom*t> for tracer ', itrc,        &
     &          TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#ifdef DIAGNOSTICS_UV
          WRITE (out,'(1x)')
          IF (Dout(idDu2d(M2rate),ng).or.Dout(idDv2d(M2rate),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2rate)',                    &
     &       'Write out 2D momentum acceleration.'
          IF (Dout(idDu2d(M2pgrd),ng).or.Dout(idDv2d(M2pgrd),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2pgrd)',                    &
     &       'Write out 2D momentum pressure gradient.'
# ifdef UV_COR
          IF (Dout(idDu2d(M2fcor),ng).or.Dout(idDv2d(M2fcor),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2fcor)',                    &
     &       'Write out 2D momentum Coriolis force.'
# endif
# ifdef UV_ADV
          IF (Dout(idDu2d(M2hadv),ng).or.Dout(idDv2d(M2hadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2hadv)',                    &
     &       'Write out 2D momentum horizontal advection.'
          IF (Dout(idDu2d(M2xadv),ng).or.Dout(idDv2d(M2xadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2xadv)',                    &
     &       'Write out 2D momentum horizontal X-advection.'
          IF (Dout(idDu2d(M2yadv),ng).or.Dout(idDv2d(M2yadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2yadv)',                    &
     &       'Write out 2D momentum horizontal Y-advection.'
# endif
# ifdef NEARSHORE_MELLOR
          IF (Dout(idDu2d(M2hrad),ng).or.Dout(idDv2d(M2hrad),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2hrad)',                    &
     &       'Write out 2D momentum horizontal radiation stress.'
# endif
# if defined UV_VIS2 || defined UV_VIS4
          IF (Dout(idDu2d(M2hvis),ng).or.Dout(idDv2d(M2hvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2hvis)',                    &
     &       'Write out 2D momentum horizontal viscosity.'
          IF (Dout(idDu2d(M2xvis),ng).or.Dout(idDv2d(M2xvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2xvis)',                    &
     &       'Write out 2D momentum horizontal X-viscosity.'
          IF (Dout(idDu2d(M2yvis),ng).or.Dout(idDv2d(M2yvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2yvis)',                    &
     &       'Write out 2D momentum horizontal Y-viscosity.'
# endif
          IF (Dout(idDu2d(M2sstr),ng).or.Dout(idDv2d(M2sstr),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2sstr)',                    &
     &       'Write out 2D momentum surface stress.'
          IF (Dout(idDu2d(M2bstr),ng).or.Dout(idDv2d(M2bstr),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M2bstr)',                    &
     &       'Write out 2D momentum bottom stress.'
# ifdef SOLVE3D
          WRITE (out,'(1x)')
          IF (Dout(idDu3d(M3rate),ng).or.Dout(idDv3d(M3rate),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3rate)',                    &
     &       'Write out 3D momentum acceleration.'
          IF (Dout(idDu3d(M3pgrd),ng).or.Dout(idDv3d(M3pgrd),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3pgrd)',                    &
     &       'Write out 3D momentum pressure gradient.'
#  ifdef UV_COR
          IF (Dout(idDu3d(M3fcor),ng).or.Dout(idDv3d(M3fcor),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3fcor)',                    &
     &       'Write out 3D momentum Coriolis force.'
#  endif
#  ifdef UV_ADV
          IF (Dout(idDu3d(M3hadv),ng).or.Dout(idDv3d(M3hadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3hadv)',                    &
     &       'Write out 3D momentum horizontal advection.'
          IF (Dout(idDu3d(M3xadv),ng).or.Dout(idDv3d(M3xadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3xadv)',                    &
     &       'Write out 3D momentum horizontal X-advection.'
          IF (Dout(idDu3d(M3yadv),ng).or.Dout(idDv3d(M3yadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3yadv)',                    &
     &       'Write out 3D momentum horizontal Y-advection.'
          IF (Dout(idDu3d(M3vadv),ng).or.Dout(idDv3d(M3vadv),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3vadv)',                    &
     &       'Write out 3D momentum vertical advection.'
#  endif
#  ifdef NEARSHORE_MELLOR
          IF (Dout(idDu3d(M3hrad),ng).or.Dout(idDv3d(M3hrad),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3hrad)',                    &
     &       'Write out 3D momentum horizontal radiation stress.'
          IF (Dout(idDu3d(M3vrad),ng).or.Dout(idDv3d(M3vrad),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3vrad)',                    &
     &       'Write out 3D momentum vertical radiation stress.'
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
          IF (Dout(idDu3d(M3hvis),ng).or.Dout(idDv3d(M3hvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3hvis)',                    &
     &       'Write out 3D momentum horizontal viscosity.'
          IF (Dout(idDu3d(M3xvis),ng).or.Dout(idDv3d(M3xvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3xvis)',                    &
     &       'Write out 3D momentum horizontal X-viscosity.'
          IF (Dout(idDu3d(M3yvis),ng).or.Dout(idDv3d(M3yvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3yvis)',                    &
     &       'Write out 3D momentum horizontal Y-viscosity.'
#  endif
          IF (Dout(idDu3d(M3vvis),ng).or.Dout(idDv3d(M3vvis),ng))       &
     &       WRITE (out,170) .TRUE., 'Dout(M3vvis)',                    &
     &       'Write out 3D momentum vertical viscosity.'
# endif
#endif
#if defined DIAGNOSTICS_TS && defined SOLVE3D
          WRITE (out,'(1x)')
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTrate),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTrate)',                   &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTrate),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTrate)',                   &
     &            'Write out rate of change of tracer ', itrc,          &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iThadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThadv)',                   &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iThadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThadv)',                   &
     &            'Write out horizontal advection, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTxadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTxadv)',                   &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTxadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTxadv)',                   &
     &            'Write out horizontal X-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTyadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTyadv)',                   &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTyadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTyadv)',                   &
     &            'Write out horizontal Y-advection, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTvadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvadv)',                   &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTvadv),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvadv)',                   &
     &            'Write out vertical advection, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
# if defined TS_DIF2 || defined TS_DIF4
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iThdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThdif)',                   &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#  ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iThdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iThdif)',                   &
     &            'Write out horizontal diffusion, tracer ', itrc,      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#  endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTxdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTxdif)',                   &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#  ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(i,iTxdif),ng))                              &
     &        WRITE (out,180) .TRUE., 'Dout(iTxdif)',                   &
     &            'Write out horizontal X-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#  endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTydif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTydif)',                   &
     &            'Write out horizontal Y-diffusion , tracer ', itrc,   &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#  ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTydif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTydif)',                   &
     &            'Write out horizontal Y-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#  endif
#  if defined MIX_GEO_TS || defined MIX_ISO_TS
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTsdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTsdif)',                   &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#   ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTsdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTsdif)',                   &
     &            'Write out horizontal S-diffusion, tracer ', itrc,    &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#   endif
#  endif
# endif
          DO itrc=1,NAT
            IF (Dout(idDtrc(itrc,iTvdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvdif)',                   &
     &            'Write out vertical diffusion, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO i=1,NPT
            itrc=inert(i)
            IF (Dout(idDtrc(itrc,iTvdif),ng))                           &
     &        WRITE (out,180) .TRUE., 'Dout(iTvdif)',                   &
     &            'Write out vertical diffusion, tracer ', itrc,        &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#if defined HDF5 && defined DEFLATE
          WRITE (out,120) shuffle, 'shuffle',                           &
     &          'NetCDF-4/HDF5 file format shuffle filer flag.'
          WRITE (out,120) deflate, 'deflate',                           &
     &          'NetCDF-4/HDF5 file format deflate filer flag.'
          WRITE (out,120) deflate_level, 'deflate_level',               &
     &          'NetCDF-4/HDF5 file format deflate level parameter.'
#endif
!
!-----------------------------------------------------------------------
!  Report output/input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,220)
#if defined FOUR_DVAR || defined ENKF_RESTART
          WRITE (out,230) 'Output DA Initial/Restart File:  ',          &
     &                    TRIM(DAI(ng)%name)
#endif
#ifdef PROPAGATOR
          WRITE (out,230) '       Output GST Restart File:  ',          &
     &                    TRIM(GST(ng)%name)
#endif
          WRITE (out,230) '           Output Restart File:  ',          &
     &                    TRIM(RST(ng)%name)
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).eq.0) THEN
              WRITE (out,230) '           Output History File:  ',      &
     &                        TRIM(HIS(ng)%name)
            ELSE
              WRITE (out,230) '      Prefix for History Files:  ',      &
     &                        TRIM(HIS(ng)%base)
            END IF
          END IF
#ifdef TANGENT
          IF (ndefTLM(ng).eq.0) THEN
            WRITE (out,230) '           Output Tangent File:  ',        &
     &                      TRIM(TLM(ng)%name)
          ELSE
            WRITE (out,230) '      Prefix for Tangent Files:  ',        &
     &                      TRIM(TLM(ng)%base)
          END IF
#endif
#ifdef WEAK_CONSTRAINT
          WRITE (out,230) '   Output Impulse Forcing File:  ',          &
     &                    TRIM(TLF(ng)%name)
#endif
#ifdef ADJOINT
          IF (ndefADJ(ng).eq.0) THEN
            WRITE (out,230) '           Output Adjoint File:  ',        &
     &                      TRIM(ADM(ng)%name)
          ELSE
            WRITE (out,230) '      Prefix for Adjoint Files:  ',        &
     &                      TRIM(ADM(ng)%base)
          END IF
#endif
#if defined FORWARD_WRITE && !defined FOUR_DVAR
          WRITE (out,230) '     Output Forward State File:  ',          &
     &                    TRIM(FWD(ng)%name)
#endif
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          IF (ndefAVG(ng).eq.0) THEN
            WRITE (out,230) '          Output Averages File:  ',        &
     &                      TRIM(AVG(ng)%name)
          ELSE
            WRITE (out,230) '     Prefix for Averages Files:  ',        &
     &                      TRIM(AVG(ng)%base)
          END IF
#endif
#if defined AVERAGES2
          IF (ndefAVG2(ng).eq.0) THEN
            WRITE (out,230) '          Output Averages File:  ',        &
     &                      TRIM(AVG2(ng)%name)
          ELSE
            WRITE (out,230) '     Prefix for Averages Files:  ',        &
     &                      TRIM(AVG2(ng)%base)
          END IF
#endif
#if defined HISTORY2
          IF (ndefHIS2(ng).eq.0) THEN
            WRITE (out,230) '          Output Surface File:  ',        &
     &                      TRIM(HIS2(ng)%name)
          ELSE
            WRITE (out,230) '     Prefix for Surface Files:  ',        &
     &                      TRIM(HIS2(ng)%base)
          END IF
#endif
#ifdef DIAGNOSTICS
          IF (ndefDIA(ng).eq.0) THEN
            WRITE (out,230) '       Output Diagnostics File:  ',        &
     &                      TRIM(DIA(ng)%name)
          ELSE
            WRITE (out,230) '  Prefix for Diagnostics Files:  ',        &
     &                      TRIM(DIA(ng)%base)
          END IF
#endif
#ifdef STATIONS
          WRITE (out,230) '          Output Stations File:  ',          &
     &                    TRIM(STA(ng)%name)
#endif
#ifdef FLOATS
          WRITE (out,230) '            Output Floats File:  ',          &
     &                    TRIM(FLT(ng)%name)
#endif
#if defined MODEL_COUPLING || defined COAWST_COUPLING
          WRITE (out,230) '      Physical parameters File:  ',          &
     &                    TRIM(Iname)
#endif
#ifndef ANA_GRID
          fname=GRD(ng)%name
          IF (.not.find_file(ng, fname, 'GRDNAME')) GO TO 30
          WRITE (out,230) '               Input Grid File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef NESTING
          fname=NGCname
          IF (.not.find_file(ng, fname, 'NGCNAME')) GO TO 30
          WRITE (out,230) 'Nesting grid connectivity File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef INI_FILE
# ifdef NONLINEAR
          fname=INI(ng)%name
          IF (.not.find_file(ng, fname, 'ININAME')) GO TO 30
          WRITE (out,230) '  Input Nonlinear Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined TANGENT && \
    !(defined FOUR_DVAR         || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS  || defined SANITY_CHECK        || \
      defined SENSITIVITY_4DVAR || defined TLM_CHECK)
          fname=ITL(ng)%name
          IF (.not.find_file(ng, fname, 'ITLNAME')) GO TO 30
          WRITE (out,230) '    Input Tangent Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined WEAK_CONSTRAINT && \
    !(defined W4DPSAS        || defined W4DPSAS_SENSITIVITY)
          fname=IRP(ng)%name
          IF (.not.find_file(ng, fname, 'IRPNAME')) GO TO 30
          WRITE (out,230) 'Input Representer Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined ADJOINT && \
    !(defined AD_SENSITIVITY      || defined FOUR_DVAR         || \
      defined IS4DVAR_SENSITIVITY || defined FORCING_SV        || \
      defined OPT_OBSERVATIONS    || defined OPT_PERTURBATION  || \
      defined SANITY_CHECK        || defined SENSITIVITY_4DVAR || \
      defined SO_SEMI             || defined STOCHASTIC_OPT    || \
      defined HESSIAN_SV          || defined HESSIAN_SO        || \
      defined HESSIAN_FSV)
          fname=IAD(ng)%name
          IF (.not.find_file(ng, fname, 'IADNAME')) GO TO 30
          WRITE (out,230) '    Input Adjoint Initial File:  ',          &
     &                    TRIM(fname)
# endif
#endif
#ifndef ANA_PSOURCE
          IF (LuvSrc(ng).or.LwSrc(ng).or.(ANY(LtracerSrc(:,ng)))) THEN
            fname=SSF(ng)%name
            IF (.not.find_file(ng, fname, 'SSFNAME')) GO TO 30
            WRITE (out,230) '      Input Sources/Sinks File:  ',        &
     &                      TRIM(fname)
          END IF
#endif
#ifdef FRC_FILE
          DO i=1,nFfiles(ng)
            DO ifile=1,FRC(i,ng)%Nfiles
              fname=FRC(i,ng)%files(ifile)
              IF (.not.find_file(ng, fname, 'FRCNAME')) GO TO 30
              IF (ifile.eq.1) THEN
                WRITE (out,310) '         Input Forcing File ', i,      &
     &                          ':  ', TRIM(fname)
              ELSE
                WRITE (out,'(35x,a)') TRIM(fname)
              END IF
            END DO
          END DO
#endif
          DO i=1,nCLMfiles(ng)
            IF (CLM_FILE(ng)) THEN
              DO ifile=1,CLM(i,ng)%Nfiles
                fname=CLM(i,ng)%files(ifile)
                IF (.not.find_file(ng, fname, 'CLMNAME')) GO TO 30
                IF (ifile.eq.1) THEN
                  WRITE (out,230) '        Input Climatology File:  ',  &
     &                            TRIM(fname)
                ELSE
                  WRITE (out,'(35x,a)') TRIM(fname)
                END IF
              END DO
            END IF
          END DO
#ifndef ANA_NUDGCOEF
          IF (LnudgeM2CLM(ng).or.LnudgeM3CLM(ng).or.                    &
     &        (ANY(LnudgeTCLM(:,ng)))) THEN
            fname=NUD(ng)%name
            IF (.not.find_file(ng, fname, 'NUDNAME')) GO TO 30
            WRITE (out,230) ' Input Nudge Coefficients File:  ',        &
     &                      TRIM(fname)
          END IF
#endif
#if defined FORWARD_READ && !(defined FOUR_DVAR || defined PICARD_TEST)
          fname=FWD(ng)%name
          IF (.not.find_file(ng, fname, 'FWDNAME')) GO TO 30
          WRITE (out,230) '      Input Forward State File:  ',          &
     &                    TRIM(fname)
#endif
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
          fname=ADS(ng)%name
          IF (.not.find_file(ng, fname, 'ADSNAME')) GO TO 30
          WRITE (out,230) 'Input Adjoint Sensitivity File:  ',          &
     &                    TRIM(fname)
#endif
!
          IF (ObcData(ng)) THEN
            DO i=1,nBCfiles(ng)
              DO ifile=1,BRY(i,ng)%Nfiles
                fname=BRY(i,ng)%files(ifile)
                IF (.not.find_file(ng, fname, 'BRYNAME')) GO TO 30
                IF (ifile.eq.1) THEN
                  WRITE (out,310) '       Input Boundary File ', i,   &
     &                            ':  ', TRIM(fname)
                ELSE
                  WRITE (out,'(35x,a)') TRIM(fname)
                END IF
              END DO
            END DO
          END IF
#ifdef STATIONS
          fname=sposnam
          IF (.not.find_file(ng, fname, 'SPOSNAM')) GO TO 30
          WRITE (out,230) '        Station positions File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef FOUR_DVAR
          fname=aparnam
          IF (.not.find_file(ng, fname, 'APARNAM')) GO TO 30
          WRITE (out,230) '  Assimilation Parameters File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef FLOATS
          fname=fposnam
          IF (.not.find_file(ng, fname, 'FPOSNAM')) GO TO 30
          WRITE (out,230) ' Initial Floats Positions File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef BIOLOGY
          fname=bparnam
          IF (.not.find_file(ng, fname, 'BPARNAM')) GO TO 30
          WRITE (out,230) '       Biology Parameters File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef VEGETATION
          fname=vegnam
          IF (.not.find_file(ng, fname, 'VEGNAM')) GO TO 30
          WRITE (out,230) '       Vegetation Parameters File:  ',       &
     &                    TRIM(fname)
#endif
#ifdef ICE_MODEL
          fname=iparnam
          IF (.not.find_file(ng, fname, 'IPARNAM')) GO TO 30
          WRITE (out,230) '     Ice Model Parameters File:  ',          &
     &                    TRIM(fname)
#endif
          fname=varname
          IF (.not.find_file(ng, fname, 'VARNAME')) GO TO 30
          GO TO 40
  30      IF (Master) WRITE (out,270) TRIM(fname)
          exit_flag=4
          RETURN
  40      CONTINUE
        END DO
        IF (Nuser.gt.0) THEN
          WRITE (out,230) '        Input/Output USER File:  ',          &
     &                    TRIM(USRname)
        END IF
!
!-----------------------------------------------------------------------
!  Report generic USER parameters.
!-----------------------------------------------------------------------
!
        IF (Nuser.gt.0) THEN
          WRITE (out,240)
          DO i=1,Nuser
            WRITE (out,250) user(i), i, i
          END DO
        END IF
      END IF

#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_EOFS    || defined POSTERIOR_ERROR_F || \
    defined POSTERIOR_ERROR_I)
!
!-----------------------------------------------------------------------
!  If weak constraint and estimating posterior analysis error
!  covariance matrix, stop it using Nouter > 1. Currently, the
!  analysis is only possible for Nouter = 1.
!-----------------------------------------------------------------------
!
      IF (Nouter.gt.1) THEN
        IF (Master) THEN
          WRITE (out,330) 'Nouter = ', Nouter,                          &
     &        'Posterior analysis error available for Nouter=1 only.'
        END IF
        exit_flag=5
        RETURN
      END IF
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Rescale active tracer parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NAT+NPT
          itrc=i
# ifdef T_PASSIVE
          IF (i.gt.NAT) itrc=inert(i-NAT)
# endif
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
          nl_tnu4(itrc,ng)=SQRT(ABS(nl_tnu4(itrc,ng)))
#ifdef ADJOINT
          ad_tnu4(itrc,ng)=SQRT(ABS(ad_tnu4(itrc,ng)))
#endif
#if defined TANGENT || defined TL_IOMS
          tl_tnu4(itrc,ng)=SQRT(ABS(tl_tnu4(itrc,ng)))
#endif
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(itrc,ng).gt.0.0_r8) THEN
            Tnudg(itrc,ng)=1.0_r8/(Tnudg(itrc,ng)*86400.0_r8)
          ELSE
            Tnudg(itrc,ng)=0.0_r8
          END IF
        END DO
        IF (Tnudg_SSS(ng).gt.0.0_r8) THEN
          Tnudg_SSS(ng)=1.0_r8/(Tnudg_SSS(ng)*86400.0_r8)
        END IF
      END DO
#endif

  50  FORMAT (/,' READ_PhyPar - Error while processing line: ',/,a)
#ifdef DISTRIBUTE
  60  FORMAT (/,1x,a,/,                                                 &
     &        /,1x,'Operating system : ',a,                             &
     &        /,1x,'CPU/hardware     : ',a,                             &
     &        /,1x,'Compiler system  : ',a,                             &
     &        /,1x,'Compiler command : ',a,                             &
     &        /,1x,'Compiler flags   : ',a,/,                           &
     &        /,1x,'Input Script  : ',a,/,                              &
     &        /,1x,'SVN Root URL  : ',a,                                &
     &        /,1x,'SVN Revision  : ',a,/,                              &
     &        /,1x,'Local Root    : ',a,                                &
     &        /,1x,'Header Dir    : ',a,                                &
     &        /,1x,'Header file   : ',a,                                &
     &        /,1x,'Analytical Dir: ',a)
#else
  60  FORMAT (/,1x,a,/,                                                 &
     &        /,1x,'Operating system : ',a,                             &
     &        /,1x,'CPU/hardware     : ',a,                             &
     &        /,1x,'Compiler system  : ',a,                             &
     &        /,1x,'Compiler command : ',a,                             &
     &        /,1x,'Compiler flags   : ',a,/,                           &
     &        /,1x,'SVN Root URL  : ',a,                                &
     &        /,1x,'SVN Revision  : ',a,/,                              &
     &        /,1x,'Local Root    : ',a,                                &
     &        /,1x,'Header Dir    : ',a,                                &
     &        /,1x,'Header file   : ',a,                                &
     &        /,1x,'Analytical Dir: ',a)
#endif
  70  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Nodes: ',i3,',',2x,'Tiling: ',i3.3,      &
     &        'x',i3.3)
  80  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i2.2,1x,          &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI * NtileJ  must be equal to the number of ', &
     &        'parallel nodes.',                                        &
     &        /,12x,'Change -np value to mpirun or',                    &
     &        /,12x,'change domain partition in input script.')
  90  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Threads: ',i2,',',2x,'Tiling: ',i3.3,    &
     &        'x',i3.3)
 100  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i3.3,1x,          &
     &        'partition or number of parallel threads.',               &
     &        /,12x,'NtileI*NtileJ must be a positive multiple of the', &
     &        ' number of threads.',                                    &
     &        /,12x,'Change number of threads (environment variable) ', &
     &        'or',/,12x,'change domain partition in input script.')
 110  FORMAT (/,/,' Physical Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
 120  FORMAT (1x,i10,2x,a,t32,a)
 130  FORMAT (1x,i10,2x,a,t32,a,/,t34,a)
 140  FORMAT (f11.3,2x,a,t32,a)
 150  FORMAT (f11.2,2x,a,t32,a)
 160  FORMAT (f11.3,2x,a,t32,a,/,t34,a)
 170  FORMAT (10x,l1,2x,a,t32,a)
 180  FORMAT (10x,l1,2x,a,t32,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t32,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t32,a,/,t34,a,i2.2,':',1x,a)
 195  FORMAT (1p,e11.4,2x,a,t32,a,i2.2,':',1x,a)
 200  FORMAT (1p,e11.4,2x,a,t32,a)
 210  FORMAT (1p,e11.4,2x,a,t32,a,/,t34,a)
 220  FORMAT (/,' Output/Input Files:',/)
 230  FORMAT (2x,a,a)
 240  FORMAT (/,' Generic User Parameters:',/)
 250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t32,                         &
     &        'User parameter ',i2.2,'.')
 260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
     &        i4,/,15x,a)
 265  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,           &
     &        1p,e11.4,/,15x,a)
 270  FORMAT (/,' READ_PHYPAR - could not find input file:  ',a)
 280  FORMAT (/,' READ_PHYPAR - variable info not yet loaded, ', a)
 290  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,i4,    &
     &        /,15x,a)
 300  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,'(',   &
     &        i2.2,')',/,15x,a)
 310  FORMAT (2x,a,i2.2,a,a)
 320  FORMAT (/,' READ_PHYPAR - could not find input parameter: ', a,   &
     &        /,15x,'in ROMS standard input script.',/,15x,a)
 330  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i4,/,15x,a)

      RETURN
      END SUBROUTINE read_PhyPar
