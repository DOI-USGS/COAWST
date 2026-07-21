#include "cppdefs.h"
      MODULE inp_par_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This routine reads in input model parameters from standard input.   !
!  It also writes out these parameters to standard output.             !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
#ifdef GRID_EXTRACT
      USE mod_netcdf
# if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
# endif
#endif
      USE mod_scalars
#ifdef DISTRIBUTE
      USE mod_strings
#endif
# ifdef UV_WAVEDRAG
      USE mod_tides, ONLY : nstep_hour
# endif
!
      USE dateclock_mod,    ONLY : get_date
#ifdef DISTRIBUTE
      USE distribute_mod,   ONLY : mp_bcasti, mp_bcasts
#endif
      USE lbc_mod,          ONLY : lbc_report
      USE ran_state,        ONLY : ran_seed
#ifdef NESTING
      USE set_contact_mod,  ONLY : set_contact
#endif
      USE stdinp_mod,       ONLY : stdinp_unit
      USE strings_mod,      ONLY : FoundError
#ifdef SOLVE3D
      USE tadv_mod,         ONLY : tadv_report
#endif
      USE tile_indices_mod, ONLY : tile_indices, tile_obs_bounds
!
      implicit none
!
      PUBLIC  :: inp_par
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE inp_par (model)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      logical :: GotFile, Lwrite
!
      integer :: Nghost, tile
      integer :: Imin, Imax, Jmin, Jmax
#ifdef GRID_EXTRACT
      integer :: I_padd, J_padd
#endif
#ifdef DISTRIBUTE
      integer :: MaxHaloLenI, MaxHaloLenJ
#endif
      integer :: ibry, inp, out, i, ic, ifield, itrc, j, ng, npts
      integer :: sequence, varid
!
      real(r8) :: cff
      real(r8), parameter :: epsilon = 1.0E-8_r8
      real(r8), parameter :: spv = 0.0_r8
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
!
#if defined DISTRIBUTE || defined COAWST_COUPLING
!
!  Get in ROMS standard input script filename (Iname) and and open it
!  as a regular formatted file in distributed-memory configurations.
!
      inp=stdinp_unit(Master, GotFile)
      out=stdout
      Lwrite=Master
!
      IF (.not.GotFile) THEN
        IF (Master) WRITE (out,10)
 10     FORMAT (/,' INP_PAR - Unable to ROMS standard input file, ',   &
                'Iname')
        exit_flag=2
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#else
!
!  Set standard inpur and output units.
!
      inp=stdinp
      out=stdout
      Lwrite=Master
#endif
#if defined SPLIT_4DVAR && SUPPRESS_REPORT
!
!  Supress reporting the information in the split 4D-Var algorithm when
!  appending into standard output.
!
      IF (Lappend) THEN
        Lwrite=.FALSE.
      END IF
#endif
!
!  Get current date.
!
#ifdef DISTRIBUTE
      IF (Master) CALL get_date (date_str)
      CALL mp_bcasts (1, model, date_str)
#else
      CALL get_date (date_str)
#endif
!
!-----------------------------------------------------------------------
!  Read in physical model input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) WRITE (out,20) version, TRIM(date_str)
 20   FORMAT (80('-'),/,                                                &
              ' Model Input Parameters:  ROMS version ',a,/,            &
     &        26x,a,/,80('-'))
!
!  Process ROMS standard input Iname script.
!
      CALL read_PhyPar (model, inp, out, Lwrite)
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
#endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#if defined SPLIT_4DVAR && SUPPRESS_REPORT
!
!  If appending into standard output file, supress the reporting of
!  information. Turn of "LwrtInfo" switch.
!
      IF (Lappend) THEN
        DO ng=1,Ngrids
          LwrtInfo(ng)=.FALSE.
        END DO
      END IF
#endif
#ifdef SEAICE
!
!-----------------------------------------------------------------------
!  Read in sea-ice model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (15, FILE=TRIM(iparnam), FORM='formatted', STATUS='old')

      CALL read_IcePar (model, 15, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Read in biological model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (25, FILE=TRIM(bparnam), FORM='formatted', STATUS='old')

      CALL read_BioPar (model, 25, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Read in sediment model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (35, FILE=TRIM(sparnam), FORM='formatted', STATUS='old')

      CALL read_SedPar (model, 35, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef NESTING
!
!-----------------------------------------------------------------------
!  Read in nesting contact points NetCDF file and allocate and
!  initialize several structures and variables.
!-----------------------------------------------------------------------
!
      CALL set_contact (1, model)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Set application domain parameters and switches.
!-----------------------------------------------------------------------
!
!  Set switch for three ghost-points in the halo region.
!
#ifdef SOLVE3D
      ThreeGhostPoints=ANY(Hadvection(:,:)%MPDATA).or.                  &
     &                 ANY(Hadvection(:,:)%HSIMT)
#endif
#ifdef UV_VIS4
      ThreeGhostPoints=.TRUE.
#endif
!
!  Determine the number of ghost-points in the halo region.
!
      IF (ThreeGhostPoints) THEN
        NghostPoints=3
      ELSE
        NghostPoints=2
      END IF
      IF (ANY(CompositeGrid).or.ANY(RefinedGrid)) THEN
        NghostPoints=MAX(3,NghostPoints)
      END IF
!
!  Determine the switch to process input open boundary conditions data.
!
!  In nesting applications, the lateral boundary conditions data is
!  is needed only by the main coarser grid (RefineScale(ng)=0).
!
      DO ng=1,Ngrids
        IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
          LprocessOBC(ng)=.TRUE.
        END IF
      END DO

#if defined SSH_TIDES || defined UV_TIDES
!
!  Determine the switch to process input tidal forcing data.
!
!  In nesting applications, the tides are processed only by the main
!  coarser grid (RefineScale(ng)=0) and the other grids get tidal
!  forcing from the contact areas.
!
      DO ng=1,Ngrids
        IF (.not.(RefinedGrid(ng).and.RefineScale(ng).gt.0)) THEN
          LprocessTides(ng)=.TRUE.
        END IF
      END DO
#endif
      CALL tile_indices (model, Im, Jm, Lm, Mm,                         &
     &                   BOUNDS, DOMAIN, IOBOUNDS)

#ifdef GRID_EXTRACT
!
!-----------------------------------------------------------------------
!  If extracting output solution, set application domain decomposition
!  bounds, indices, and switches per tile partition for all grids.
!-----------------------------------------------------------------------
!
!  Inquire about the extract grid dimensions.
!
     DO ng=1,Ngrids
        SELECT CASE (GRX(ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_dim (ng, iNLM, TRIM(GRX(ng)%name))
#  if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_dim (ng, iNLM, TRIM(GRX(ng)%name))
#  endif
          CASE DEFAULT
            IF (Master) WRITE (stdout,40) GRX(ng)%IOtype
 40         FORMAT (/,' INP_PAR - Illegal output type,',                &
     &              ' io_type = ',i0)
            exit_flag=3
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        xtr_Lm(ng)=0
        xtr_Mm(ng)=0
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'xi_rho') THEN
            xtr_Lm(ng)=dim_size(i)-2            ! I-computational points
            I_padd=(Lm(ng)+2)/2-(Lm(ng)+1)/2
            xtr_Im(ng)=xtr_Lm(ng)+I_padd
          ELSE IF (TRIM(dim_name(i)).eq.'eta_rho') THEN
            xtr_Mm(ng)=dim_size(i)-2            ! J-computational points
            J_padd=(Mm(ng)+2)/2-(Mm(ng)+1)/2
            xtr_Jm(ng)=xtr_Mm(ng)+J_padd
          END IF
        END DO
        IF (xtr_Lm(ng).eq.0) THEN
          WRITE (stdout,20) 'xi_rho', TRIM(GRX(ng)%name)
          exit_flag=2
          RETURN
        END IF
        IF (xtr_Jm(ng).eq.0) THEN
          WRITE (stdout,20) 'eta_rho', TRIM(GRX(ng)%name)
          exit_flag=2
          RETURN
        END IF
 45     FORMAT (/,' INP_PAR - error inquiring dimension: ',a,2x,        &
     &          'in input NetCDF file: ',a)
      END DO
!
!  Set decomposition bounds, indices, and switches per tile partition
!  for all grids.
!
      CALL tile_indices (model, xtr_Im, xtr_Jm, xtr_Lm, xtr_Mm,         &
     &                   xtr_BOUNDS, xtr_DOMAIN, xtr_IOBOUNDS)
#endif
!
!-----------------------------------------------------------------------
!  Set minimum and maximum fractional coordinates for processing
!  observations.
!-----------------------------------------------------------------------
!
      CALL tile_obs_bounds (model, Im, Jm, Lm, Mm,                      &
     &                      DOMAIN)
!
!-----------------------------------------------------------------------
!  Check tile partition starting and ending (I,J) indices for illegal
!  domain decomposition parameters NtileI and NtileJ in standard input
!  file.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        DO ng=1,Ngrids
#ifdef SOLVE3D
          WRITE (stdout,50) ng, Lm(ng), Mm(ng), N(ng),                  &
     &                      NtileI(ng), NtileJ(ng)
#else
          WRITE (stdout,50) ng, Lm(ng), Mm(ng),                         &
     &                      NtileI(ng), NtileJ(ng)
#endif
#if !defined DISTRIBUTE && defined ADJOINT
          IF ((NtileI(ng).ne.1).or.(NtileJ(ng).ne.1)) THEN
            WRITE (stdout,60)
            exit_flag=6
            RETURN
          END IF
#endif
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
#ifdef SOLVE3D
            npts=(BOUNDS(ng)%Iend(tile)-                                &
     &            BOUNDS(ng)%Istr(tile)+1)*                             &
     &           (BOUNDS(ng)%Jend(tile)-                                &
     &            BOUNDS(ng)%Jstr(tile)+1)*N(ng)
#else
            npts=(BOUNDS(ng)%Iend(tile)-                                &
     &            BOUNDS(ng)%Istr(tile)+1)*                             &
     &           (BOUNDS(ng)%Jend(tile)-                                &
     &            BOUNDS(ng)%Jstr(tile)+1)
#endif
            WRITE (stdout,70) tile,                                     &
     &                        BOUNDS(ng)%Istr(tile),                    &
     &                        BOUNDS(ng)%Iend(tile),                    &
     &                        BOUNDS(ng)%Jstr(tile),                    &
     &                        BOUNDS(ng)%Jend(tile),                    &
     &                        npts
            IF ((BOUNDS(ng)%Iend(tile)-                                 &
     &           BOUNDS(ng)%Istr(tile)+1).lt.2) THEN
              WRITE (stdout,80) ng, 'NtileI = ', NtileI(ng),            &
     &                              'Lm = ', Lm(ng),                    &
     &                              'Istr = ', BOUNDS(ng)%Istr(tile),   &
     &                              '  Iend = ', BOUNDS(ng)%Iend(tile), &
     &                              'NtileI'
              exit_flag=6
              RETURN
            END IF
            IF ((BOUNDS(ng)%Jend(tile)-                                 &
     &           BOUNDS(ng)%Jstr(tile)+1).lt.2) THEN
              WRITE (stdout,80) ng, 'NtileJ = ', NtileJ(ng),            &
     &                              'Mm = ', Mm(ng),                    &
     &                              'Jstr = ', BOUNDS(ng)%Jstr(tile),   &
     &                              '  Jend = ', BOUNDS(ng)%Jend(tile), &
     &                              'NtileJ'
              exit_flag=6
              RETURN
            END IF
          END DO
        END DO
#ifdef SOLVE3D
 50     FORMAT (/,' Tile partition information for Grid ',i2.2,':',2x,  &
     &          i0,'x',i0,'x',i0,2x,'tiling: ',i0,'x',i0,/,/,           &
     &          5x,'tile',5x,'Istr',5x,'Iend',5x,'Jstr',5x,'Jend',      &
     &          5x,'Npts',/)
#else
 50     FORMAT (/,' Tile partition information for Grid ',i2.2,':',2x,  &
     &          i0,'x',i0,2x,'tiling: ',i0,'x',i0,/,/,                  &
     &          5x,'tile',5x,'Istr',5x,'Iend',5x,'Jstr',5x,'Jend',      &
     &          5x,'Npts',/)
#endif
#if !defined DISTRIBUTE && defined ADJOINT
 60     FORMAT (/,' INP_PAR - illegal domain decomposition for the ',   &
     &                       'Adjoint model.',/,11x,'Partitions are ',  &
     &          'allowed in distributed-menory (MPI) applications.'/)
#endif
 70     FORMAT (5(4x,i5),1x,i8)
 80     FORMAT (/,' INP_PAR - domain decomposition error in input ',    &
     &                        'script file for grid: ',i2.2,/,          &
     &          /,11x,'The domain partition parameter, ',a,i0,          &
     &          /,11x,'is incompatible with grid size, ',a,i0,          &
     &          /,11x,'because it yields too small tile, ',a,i0,a,i0,   &
     &          /,11x,'Decrease partition parameter: ',a)
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
#endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Report tile minimum and maximum fractional grid coordinates.
!
      DO ng=1,Ngrids
        IF (Master.and.Lwrite) THEN
          WRITE (stdout,90) ng
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            WRITE (stdout,100) tile,                                    &
     &                         DOMAIN(ng)%Xmin_rho(tile),               &
     &                         DOMAIN(ng)%Xmax_rho(tile),               &
     &                         DOMAIN(ng)%Ymin_rho(tile),               &
     &                         DOMAIN(ng)%Ymax_rho(tile), 'RHO-points'
          END DO
          WRITE (stdout,'(1x)')
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            WRITE (stdout,100) tile,                                    &
     &                         DOMAIN(ng)%Xmin_u(tile),                 &
     &                         DOMAIN(ng)%Xmax_u(tile),                 &
     &                         DOMAIN(ng)%Ymin_u(tile),                 &
     &                         DOMAIN(ng)%Ymax_u(tile), '  U-points'
          END DO
          WRITE (stdout,'(1x)')
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            WRITE (stdout,100) tile,                                    &
     &                         DOMAIN(ng)%Xmin_v(tile),                 &
     &                         DOMAIN(ng)%Xmax_v(tile),                 &
     &                         DOMAIN(ng)%Ymin_v(tile),                 &
     &                         DOMAIN(ng)%Ymax_v(tile), '  V-points'
          END DO
 90       FORMAT (/,' Tile minimum and maximum fractional coordinates', &
     &            ' for Grid ',i2.2,':'/,                               &
#ifdef FULL_GRID
     &            '   (interior and boundary points)',/,/,              &
#else
     &            '   (interior points only)',/,/,                      &
#endif
     &            5x,'tile',5x,'Xmin',5x,'Xmax',5x,'Ymin',5x,'Ymax',    &
     &            5x,'grid',/)
 100      FORMAT (5x,i4,4f9.2,2x,a)
        END IF
      END DO

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Determine the maximum tile lengths in XI and ETA directions for
!  distributed-memory communications.  Notice that halo size are
!  increased by few points to allow exchanging of private arrays.
!-----------------------------------------------------------------------
!
      IF (ANY(EWperiodic).or.ANY(NSperiodic)) THEN
        Nghost=NghostPoints+1
      ELSE
        Nghost=NghostPoints
      END IF

      DO ng=1,Ngrids
        MaxHaloLenI=0
        MaxHaloLenJ=0
        HaloBry(ng)=Nghost
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          Imin=BOUNDS(ng)%LBi(tile)-1
          Imax=BOUNDS(ng)%UBi(tile)+1
          Jmin=BOUNDS(ng)%LBj(tile)-1
          Jmax=BOUNDS(ng)%UBj(tile)+1
          MaxHaloLenI=MAX(MaxHaloLenI,(Imax-Imin+1))
          MaxHaloLenJ=MAX(MaxHaloLenJ,(Jmax-Jmin+1))
        END DO
        HaloSizeI(ng)=Nghost*MaxHaloLenI+6*Nghost
        HaloSizeJ(ng)=Nghost*MaxHaloLenJ+6*Nghost
        TileSide(ng)=MAX(MaxHaloLenI,MaxHaloLenJ)
        TileSize(ng)=MaxHaloLenI*MaxHaloLenJ
        IF (Master.and.Lwrite) THEN
          WRITE (stdout,110) ng, HaloSizeI(ng), ng, HaloSizeJ(ng),      &
     &                       ng, TileSide(ng),  ng, TileSize(ng)
 110      FORMAT (/,' Maximum halo size in XI and ETA directions:',/,   &
     &            /,'               HaloSizeI(',i1,') = ',i7,           &
     &            /,'               HaloSizeJ(',i1,') = ',i7,           &
     &            /,'                TileSide(',i1,') = ',i7,           &
     &            /,'                TileSize(',i1,') = ',i7,/)
        END IF
      END DO
#endif

#if defined FOUR_DVAR || defined VERIFICATION
!
!-----------------------------------------------------------------------
!  Read in input assimilation parameters.
!-----------------------------------------------------------------------
!
      OPEN (35, FILE=TRIM(aparnam), FORM='formatted', STATUS='old')

      CALL read_AssPar (model, 35, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Read in floats input parameters.
!-----------------------------------------------------------------------
!
      OPEN (45, FILE=TRIM(fposnam), FORM='formatted', STATUS='old')

      CALL read_FltPar (model, 45, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#if defined FLOATS && defined FLOAT_BIOLOGY
!
!-----------------------------------------------------------------------
!  Read in biological float behavior model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (50, FILE=TRIM(fbionam), FORM='formatted', STATUS='old')

      CALL read_FltBioPar (model, 50, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef STATIONS
!
!-----------------------------------------------------------------------
!  Read in stations input parameters.
!-----------------------------------------------------------------------
!
      OPEN (55, FILE=TRIM(sposnam), FORM='formatted', STATUS='old')

      CALL read_StaPar (model, 55, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
!
#ifdef VEGETATION
!
!-----------------------------------------------------------------------
!  Read in vegetation model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (65, FILE=TRIM(vegname), FORM='formatted', STATUS='old')

      CALL read_VegPar (model, 65, out, Lwrite)
#endif 
#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Report tracer advection scheme.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        WRITE (out,120) 'NLM'
 120    FORMAT (/,1x,'Tracer Advection Scheme: ',a,/,1x,24('='),/,      &
     &          /,1x,'Variable',t25,'Grid',t31,'Horizontal',            &
     &          t50,'Vertical', /,1x,'---------',t25,'----',            &
     &          t31,2('------------',7x))
      END IF
      CALL tadv_report (out, iNLM, Hadvection, Vadvection, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined ADJOINT || defined TANGENT || defined TL_IOMS
!
      IF (Master.and.Lwrite) THEN
        WRITE (out,120) 'TLM, RPM, and ADM'
      END IF
      CALL tadv_report (out, iADM, ad_Hadvection, ad_Vadvection, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
#endif
#if defined TANGENT || defined TL_IOMS
!
!  Set tracer advection scheme switches for the tangent linear models
!  (TLM and RPM) to the same values as the adjoint model.
!
      DO ng=1,Ngrids
        DO i=1,NT(ng)
          tl_Hadvection(i,ng)%AKIMA4    = ad_Hadvection(i,ng)%AKIMA4
          tl_Hadvection(i,ng)%CENTERED2 = ad_Hadvection(i,ng)%CENTERED2
          tl_Hadvection(i,ng)%CENTERED4 = ad_Hadvection(i,ng)%CENTERED4
          tl_Hadvection(i,ng)%HSIMT     = ad_Hadvection(i,ng)%HSIMT
          tl_Hadvection(i,ng)%MPDATA    = ad_Hadvection(i,ng)%MPDATA
          tl_Hadvection(i,ng)%SPLINES   = ad_Hadvection(i,ng)%SPLINES
          tl_Hadvection(i,ng)%SPLIT_U3  = ad_Hadvection(i,ng)%SPLIT_U3
          tl_Hadvection(i,ng)%UPSTREAM3 = ad_Hadvection(i,ng)%UPSTREAM3
!
          tl_Vadvection(i,ng)%AKIMA4    = ad_Vadvection(i,ng)%AKIMA4
          tl_Vadvection(i,ng)%CENTERED2 = ad_Vadvection(i,ng)%CENTERED2
          tl_Vadvection(i,ng)%CENTERED4 = ad_Vadvection(i,ng)%CENTERED4
          tl_Vadvection(i,ng)%HSIMT     = ad_Vadvection(i,ng)%HSIMT
          tl_Vadvection(i,ng)%MPDATA    = ad_Vadvection(i,ng)%MPDATA
          tl_Vadvection(i,ng)%SPLINES   = ad_Vadvection(i,ng)%SPLINES
          tl_Vadvection(i,ng)%SPLIT_U3  = ad_Vadvection(i,ng)%SPLIT_U3
          tl_Vadvection(i,ng)%UPSTREAM3 = ad_Vadvection(i,ng)%UPSTREAM3
        END DO
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Report lateral boundary conditions.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        WRITE (out,130) 'NLM'
 130    FORMAT (/,1x,'Lateral Boundary Conditions: ',a,/,1x,28('='),/,  &
     &          /,1x,'Variable',t25,'Grid',t31,'West Edge',             &
     &          t44,'South Edge', t57,'East Edge',t70,'North Edge',     &
     &          /,1x,'---------',t25,'----',t31,4('----------',3x))
        DO ifield=1,nLBCvar
          IF (idBvar(ifield).gt.0) THEN
            CALL lbc_report (out, ifield, LBC)
          END IF
        END DO

#if defined ADJOINT || defined TANGENT || defined TL_IOMS
!
        WRITE (out,130) 'TLM, RPM, and ADM'
        DO ifield=1,nLBCvar
          IF (idBvar(ifield).gt.0) THEN
            CALL lbc_report (out, ifield, ad_LBC)
          END IF
        END DO
#endif
      END IF
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Compute various constants.
!-----------------------------------------------------------------------
!
      gorho0=g/rho0
      DO ng=1,Ngrids
        dtfast(ng)=dt(ng)/REAL(ndtfast(ng),r8)
#ifdef UV_WAVEDRAG
        IF (MOD(3600,INT(dt(ng))).ne.0) THEN
          IF (Master) print *,                                          &
     &                'Model timesteps must fall on hour interval'
          exit_flag=2
        END IF
        nstep_hour = nw*3600/INT(dt(ng))
#endif
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
        nl_visc4(ng)=SQRT(ABS(nl_visc4(ng)))
#ifdef ADJOINT
        ad_visc4(ng)=SQRT(ABS(ad_visc4(ng)))
#endif
#if defined TANGENT || defined TL_IOMS
        tl_visc4(ng)=SQRT(ABS(tl_visc4(ng)))
#endif
        tkenu4(ng)=SQRT(ABS(tkenu4(ng)))
!
!  Set internal switch for activating sponge areas.
!
#ifdef SOLVE3D
        IF (LuvSponge(ng).or.                                           &
     &      ANY(LtracerSponge(:,ng))) THEN
          Lsponge(ng)=.TRUE.
        END IF
#else
        IF (LuvSponge(ng)) THEN
          Lsponge(ng)=.TRUE.
        END IF
#endif
!
!  Set switch to processing nudging coefficients for passive/active
!  boundary conditions.
!
        NudgingCoeff(ng)=ANY(LBC(:,:,ng)%nudging)
#if defined ADJOINT || defined TANGENT || defined TL_IOMS
        NudgingCoeff(ng)=NudgingCoeff(ng).or.ANY(ad_LBC(:,:,ng)%nudging)
#endif
!
!  Set internal switch for processing climatology data.
!
#ifdef SOLVE3D
# if defined TS_MIX_CLIMA && (defined TS_DIF2 || defined TS_DIF4)
        Lclimatology(ng)=.TRUE.
# endif
        IF (LsshCLM(ng).or.                                             &
            Lm2CLM (ng).or.LnudgeM2CLM(ng).or.                          &
            Lm3CLM (ng).or.LnudgeM3CLM(ng).or.                          &
            ANY(LtracerCLM(:,ng)).or.ANY(LnudgeTCLM(:,ng))) THEN
          Lclimatology(ng)=.TRUE.
        END IF
#else
        IF (LsshCLM(ng).or.                                             &
            Lm2CLM (ng).or.LnudgeM2CLM(ng)) THEN
          Lclimatology(ng)=.TRUE.
        END IF
#endif
!
!  Set internal switch for nudging to climatology fields.
!
#ifdef SOLVE3D
        IF (LnudgeM2CLM(ng).or.                                         &
     &      LnudgeM3CLM(ng).or.                                         &
     &      ANY(LnudgeTCLM(:,ng))) THEN
          Lnudging(ng)=.TRUE.
        END IF
#else
        IF (LnudgeM2CLM(ng)) THEN
          Lnudging(ng)=.TRUE.
        END IF
#endif
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
        IF (Znudg(ng).gt.0.0_r8) THEN
          Znudg(ng)=1.0_r8/(Znudg(ng)*86400.0_r8)
        ELSE
          Znudg(ng)=0.0_r8
        END IF
!
        IF (M2nudg(ng).gt.0.0_r8) THEN
          M2nudg(ng)=1.0_r8/(M2nudg(ng)*86400.0_r8)
        ELSE
          M2nudg(ng)=0.0_r8
        END IF
#ifdef SOLVE3D
!
        IF (M3nudg(ng).gt.0.0_r8) THEN
          M3nudg(ng)=1.0_r8/(M3nudg(ng)*86400.0_r8)
        ELSE
          M3nudg(ng)=0.0_r8
        END IF
#endif
!
!  Set nudging coefficients (1/s) for passive/active (outflow/inflow)
!  open boundary conditions.  Weak nudging is expected in passive
!  outflow conditions and strong nudging is expected in active inflow
!  conditions. If nudging to climatology fields, these values are
!  replaced by spatial nudging coefficients distribution in the
!  open boundary condition routines.
!
        IF (NudgingCoeff(ng)) THEN
          DO ibry=1,4
            IF (LBC(ibry,isFsur,ng)%nudging) THEN
              FSobc_out(ng,ibry)=Znudg(ng)
              FSobc_in (ng,ibry)=obcfac(ng)*Znudg(ng)
            END IF
!
            IF (LBC(ibry,isUbar,ng)%nudging.or.                         &
     &          LBC(ibry,isVbar,ng)%nudging) THEN
              M2obc_out(ng,ibry)=M2nudg(ng)
              M2obc_in (ng,ibry)=obcfac(ng)*M2nudg(ng)
            END IF
#ifdef SOLVE3D
!
            IF (LBC(ibry,isUvel,ng)%nudging.or.                         &
     &          LBC(ibry,isVvel,ng)%nudging) THEN
              M3obc_out(ng,ibry)=M3nudg(ng)
              M3obc_in (ng,ibry)=obcfac(ng)*M3nudg(ng)
            END IF
!
            DO itrc=1,NT(ng)
              IF (LBC(ibry,isTvar(itrc),ng)%nudging) THEN
                Tobc_out(itrc,ng,ibry)=Tnudg(itrc,ng)
                Tobc_in (itrc,ng,ibry)=obcfac(ng)*Tnudg(itrc,ng)
              END IF
            END DO
#endif
          END DO
        END IF

#if defined SO_SEMI        || \
   (defined STOCHASTIC_OPT && !defined STOCH_OPT_WHITE)
       SO_decay(ng)=SO_decay(ng)*86400.0_r8
#endif
!
!  Convert momentum stresses and tracer flux scales to kinematic
!  Values. Recall, that all the model fluxes are kinematic.
!
        cff=1.0_r8/rho0
        Fscale(idUsms,ng)=cff*Fscale(idUsms,ng)
        Fscale(idVsms,ng)=cff*Fscale(idVsms,ng)
        Fscale(idUbms,ng)=cff*Fscale(idUbms,ng)
        Fscale(idVbms,ng)=cff*Fscale(idVbms,ng)
        Fscale(idUbrs,ng)=cff*Fscale(idUbrs,ng)
        Fscale(idVbrs,ng)=cff*Fscale(idVbrs,ng)
        Fscale(idUbws,ng)=cff*Fscale(idUbws,ng)
        Fscale(idVbws,ng)=cff*Fscale(idVbws,ng)
        Fscale(idUbcs,ng)=cff*Fscale(idUbcs,ng)
        Fscale(idVbcs,ng)=cff*Fscale(idVbcs,ng)
        cff=1.0_r8/(rho0*Cp)
        Fscale(idTsur(itemp),ng)=cff*Fscale(idTsur(itemp),ng)
        Fscale(idTbot(itemp),ng)=cff*Fscale(idTbot(itemp),ng)
        Fscale(idSrad,ng)=cff*Fscale(idSrad,ng)
        Fscale(idLdwn,ng)=cff*Fscale(idLdwn,ng)
        Fscale(idLrad,ng)=cff*Fscale(idLrad,ng)
        Fscale(idLhea,ng)=cff*Fscale(idLhea,ng)
        Fscale(idShea,ng)=cff*Fscale(idShea,ng)
        Fscale(iddQdT,ng)=cff*Fscale(iddQdT,ng)

#ifdef SOLVE3D
!
!  Determine the number of climatology tracers to process.
!
        IF (ANY(LtracerCLM(:,ng)).or.ANY(LnudgeTCLM(:,ng))) THEN
          ic=0
          DO itrc=1,NT(ng)
            IF (LtracerCLM(itrc,ng)) THEN
              ic=ic+1
            END IF
          END DO
          NTCLM(ng)=ic
        END IF
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Set lateral boundary condition switches for the tangent linear
!  models (TLM and RPM) to the same values as the adjoint model.
!
        DO j=1,nLBCvar
          DO i=1,4
            tl_LBC(i,j,ng)%acquire     = ad_LBC(i,j,ng)%acquire
            tl_LBC(i,j,ng)%Chapman_explicit =                           &
     &                                   ad_LBC(i,j,ng)%Chapman_explicit
            tl_LBC(i,j,ng)%Chapman_implicit =                           &
     &                                   ad_LBC(i,j,ng)%Chapman_implicit
            tl_LBC(i,j,ng)%clamped     = ad_LBC(i,j,ng)%clamped
            tl_LBC(i,j,ng)%closed      = ad_LBC(i,j,ng)%closed
            tl_LBC(i,j,ng)%Flather     = ad_LBC(i,j,ng)%Flather
            tl_LBC(i,j,ng)%gradient    = ad_LBC(i,j,ng)%gradient
            tl_LBC(i,j,ng)%nested      = ad_LBC(i,j,ng)%nested
            tl_LBC(i,j,ng)%nudging     = ad_LBC(i,j,ng)%nudging
            tl_LBC(i,j,ng)%periodic    = ad_LBC(i,j,ng)%periodic
            tl_LBC(i,j,ng)%radiation   = ad_LBC(i,j,ng)%radiation
            tl_LBC(i,j,ng)%reduced     = ad_LBC(i,j,ng)%reduced
            tl_LBC(i,j,ng)%Shchepetkin = ad_LBC(i,j,ng)%Shchepetkin
          END DO
        END DO
#endif

      END DO

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Set climatology tracers (active and passive) metadata.  It needs to
!  be done here because information is needed from all input scripts.
!  The variable name and units are the same as the basic tracers. The
!  default time-variable name is the same as the variable name but with
!  the "_time" suffix.  Recall that other time-variables names are
!  allowed provided that the input NetCDF variable has the "time"
!  attribute with the appropriate value.
!-----------------------------------------------------------------------
!
      varid=last_varid
      IF (ANY(LtracerCLM).or.ANY(LnudgeTCLM)) THEN
        DO i=1,MT
          varid=varid+1
          IF (varid.gt.MV) THEN
            WRITE (stdout,130) MV, varid
            STOP
          END IF
          idTclm(i)=varid
          DO ng=1,Ngrids
            Fscale(varid,ng)=1.0_r8
            Iinfo(1,varid,ng)=r3dvar
          END DO
          WRITE (Vname(1,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(1,idTvar(i))))
          WRITE (Vname(2,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(2,idTvar(i)))), ' climatology'
          WRITE (Vname(3,varid),'(a)')                                  &
     &          TRIM(ADJUSTL(Vname(3,idTvar(i))))
          WRITE (Vname(4,varid),'(a,a)')                                &
     &          TRIM(Vname(1,varid)), ', scalar, series'
          WRITE (Vname(5,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(1,idTvar(i)))), '_time'
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Set tracers inverse nudging coeffcients metadata.  It needs to be
!  done here because information is needed from all input scripts.
!  The variable name is the same as the basic tracer but with the
!  "_NudgeCoef" suffix.
!-----------------------------------------------------------------------
!
      DO i=1,MT
        IF (ANY(LnudgeTCLM(i,:))) THEN
          varid=varid+1
          IF (varid.gt.MV) THEN
            WRITE (stdout,140) MV, varid
 140        FORMAT (/,' INP_PAR - too small dimension ',                &
     &              'parameter, MV = ',2i5,/,15x,                       &
     &              'change file  mod_ncparam.F  and recompile.')
            STOP
          END IF
          idTnud(i)=varid
          DO ng=1,Ngrids
            Fscale(varid,ng)=1.0_r8/86400        ! default units: 1/day
            Iinfo(1,varid,ng)=r3dvar
          END DO
          WRITE (Vname(1,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(1,idTvar(i)))), '_NudgeCoef'
          WRITE (Vname(2,varid),'(a,a)')                                &
     &          TRIM(ADJUSTL(Vname(2,idTvar(i)))),                      &
     &          ', inverse nudging coefficients'
          WRITE (Vname(3,varid),'(a,1x,a)')                             &
     &          TRIM(ADJUSTL(Vname(3,idTvar(i)))), 'day-1'
          WRITE (Vname(4,varid),'(a,a)')                                &
     &        TRIM(Vname(1,varid)), ', scalar'
          WRITE (Vname(5,varid),'(a)') 'nulvar'
        ELSE
          idTnud(i)=0
        END IF
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options and definitions.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        CALL checkdefs
        FLUSH (out)
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
      CALL mp_bcasts (1, model, Coptions)
#endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Initialize random number sequence so we can get identical results
!  everytime that we run the same solution.
!-----------------------------------------------------------------------
!
      sequence=759
      CALL ran_seed (sequence)
!
      RETURN
      END SUBROUTINE inp_par
!
      END MODULE inp_par_mod
