#include "cppdefs.h"
      SUBROUTINE inp_par (model)
!
!svn $Id: inp_par.F 838 2008-11-17 04:22:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in input model parameters from standard input.   !
!  It also writes out these parameters to standard output.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
#ifdef DISTRIBUTE
      USE mod_strings
#endif
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE ran_state, ONLY: ran_seed
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      logical :: Lwrite

      integer :: Itile, Jtile, Nghost, Ntiles, tile
      integer :: Imin, Imax, Jmin, Jmax
#ifdef DISTRIBUTE
      integer :: MaxHaloLenI, MaxHaloLenJ
#endif
      integer :: inp, out, itrc, ng, npts, sequence

      real(r8) :: cff
      real(r8), parameter :: spv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
!
!  Set input units.
!
#if defined DISTRIBUTE || defined MODEL_COUPLING
      Lwrite=Master
      inp=1
      out=stdout
#else
      Lwrite=Master
      inp=stdinp
      out=stdout
#endif
!
!  Get current data.
!
#ifdef DISTRIBUTE
      IF (Master) CALL get_date (date_str)
      DO ng=1,Ngrids
        CALL mp_bcasts (ng, model, date_str)
      END DO
#else
      CALL get_date (date_str)
#endif
!
!-----------------------------------------------------------------------
!  Read in physical model input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) WRITE (out,10) version, TRIM(date_str)
 10   FORMAT (/,' Model Input Parameters:  ROMS/TOMS version ',a,/,     &
     &       26x,a,/,1x,77('-'))

#ifdef MODEL_COUPLING
!
!  In multiple model coupling, the ocean input physical parameter
!  script needs to be opened and processed as a regular file.
!  Its file name is specified in the coupling standard input script.
!
      OPEN (inp, FILE=TRIM(Iname), FORM='formatted', STATUS='old',      &
     &      ERR=20)
      GO TO 40
 20   IF (Master) WRITE (stdout,30)
      exit_flag=2
      RETURN
 30   FORMAT (/,' INP_PAR - Unable to open ROMS/TOMS input script ',    &
     &              'file.')
 40   CONTINUE

#else
# ifdef DISTRIBUTE
!
!  In distributed-memory configurations, the input physical parameters
!  script is opened as a regular file.  It is read and processed by all
!  parallel nodes.  This is to avoid a very complex broadcasting of the
!  input parameters to all nodes.
!
!!    CALL my_getarg (1,Iname)
      IF (Master) CALL my_getarg (1,Iname)
      DO ng=1,Ngrids
        CALL mp_bcasts (ng, model, Iname)
      END DO
      OPEN (inp, FILE=TRIM(Iname), FORM='formatted', STATUS='old',      &
     &      ERR=20)
      GO TO 40
 20   IF (Master) WRITE (stdout,30)
      exit_flag=2
      RETURN
 30   FORMAT (/,' INP_PAR - Unable to open ROMS/TOMS input script ',    &
     &              'file.',                                            &
     &        /,11x,'In distributed-memory applications, the input',    &
     &        /,11x,'script file is processed in parallel. The Unix',   &
     &        /,11x,'routine GETARG is used to get script file name.',  &
     &        /,11x,'For example, in MPI applications make sure that',  &
     &        /,11x,'command line is something like:',/,                &
     &        /,11x,'mpirun -np 4 ocean ocean.in',/,/,11x,'and not',/,  &
     &        /,11x,'mpirun -np 4 ocean < ocean.in',/)
 40   CONTINUE
# endif
#endif

      CALL read_PhyPar (model, inp, out, Lwrite)
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
#endif
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Set lower and upper bounds indices per domain partition for all
!  nested grids.
!-----------------------------------------------------------------------
!
!  Allocate structure.
!
      IF (.not.allocated(BOUNDS)) THEN
        allocate ( BOUNDS(Ngrids) )
        DO ng=1,Ngrids
          Ntiles=NtileI(ng)*NtileJ(ng)-1
          allocate ( BOUNDS(ng) % tile (-1:Ntiles) )
          allocate ( BOUNDS(ng) % LBi  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % UBi  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % LBj  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % UBj  (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Iend (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Istr (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jstr (-1:Ntiles) )
          allocate ( BOUNDS(ng) % Jend (-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IstrU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrV(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendR(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendT(-1:Ntiles) )
          allocate ( BOUNDS(ng) % Imin (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Imax (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Jmin (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % Jmax (4,0:1,0:Ntiles) )
          allocate ( BOUNDS(ng) % IstrTU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrTV(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendTU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendTV(-1:Ntiles) )
        END DO
      END IF
!
!  Set boundary edge I- or J-indices for each variable type.
!
      DO ng=1,Ngrids
        BOUNDS(ng) % edge(iwest ,p2dvar) = 1
        BOUNDS(ng) % edge(iwest ,r2dvar) = 0
        BOUNDS(ng) % edge(iwest ,u2dvar) = 1
        BOUNDS(ng) % edge(iwest ,v2dvar) = 0

        BOUNDS(ng) % edge(ieast ,p2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(ieast ,r2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(ieast ,u2dvar) = Lm(ng)+1
        BOUNDS(ng) % edge(ieast ,v2dvar) = Lm(ng)+1

        BOUNDS(ng) % edge(isouth,p2dvar) = 1
        BOUNDS(ng) % edge(isouth,u2dvar) = 0
        BOUNDS(ng) % edge(isouth,r2dvar) = 0
        BOUNDS(ng) % edge(isouth,v2dvar) = 1

        BOUNDS(ng) % edge(inorth,p2dvar) = Mm(ng)+1
        BOUNDS(ng) % edge(inorth,r2dvar) = Mm(ng)+1
        BOUNDS(ng) % edge(inorth,u2dvar) = Mm(ng)+1
        BOUNDS(ng) % edge(inorth,v2dvar) = Mm(ng)+1
      END DO
!
!  Set tile computational indices and arrays allocation bounds.
!
      Nghost=GHOST_POINTS
      DO ng=1,Ngrids
        BOUNDS(ng) % LBij = 0
        BOUNDS(ng) % UBij = MAX(Lm(ng)+1,Mm(ng)+1)
        DO tile=-1,NtileI(ng)*NtileJ(ng)-1
          BOUNDS(ng) % tile(tile) = tile
          CALL get_tile (ng, tile, Itile, Jtile,                        &
     &                   BOUNDS(ng) % Istr(tile),                       &
     &                   BOUNDS(ng) % Iend(tile),                       &
     &                   BOUNDS(ng) % Jstr(tile),                       &
     &                   BOUNDS(ng) % Jend(tile),                       &
     &                   BOUNDS(ng) % IstrR(tile),                      &
     &                   BOUNDS(ng) % IstrT(tile),                      &
     &                   BOUNDS(ng) % IstrU(tile),                      &
     &                   BOUNDS(ng) % IendR(tile),                      &
     &                   BOUNDS(ng) % IendT(tile),                      &
     &                   BOUNDS(ng) % JstrR(tile),                      &
     &                   BOUNDS(ng) % JstrT(tile),                      &
     &                   BOUNDS(ng) % JstrV(tile),                      &
     &                   BOUNDS(ng) % JendR(tile),                      &
     &                   BOUNDS(ng) % JendT(tile),                      &
     &                   BOUNDS(ng) % IstrTU(tile),                     &
     &                   BOUNDS(ng) % JstrTV(tile),                     &
     &                   BOUNDS(ng) % IendTU(tile),                     &
     &                   BOUNDS(ng) % JendTV(tile))
          CALL get_bounds (ng, tile, 0, Nghost, Itile, Jtile,           &
     &                     BOUNDS(ng) % LBi(tile),                      &
     &                     BOUNDS(ng) % UBi(tile),                      &
     &                     BOUNDS(ng) % LBj(tile),                      &
     &                     BOUNDS(ng) % UBj(tile))
        END DO
      END DO
!
!  Set I/O processing minimum (Imin, Jmax) and maximum (Imax, Jmax)
!  indices for non-overlapping (Nghost=0) and overlapping (Nghost>0)
!  tiles for each C-grid type variable.
!
      Nghost=GHOST_POINTS
      DO ng=1,Ngrids
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_bounds (ng, tile, p2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(1,0,tile),                 &
     &                     BOUNDS(ng) % Imax(1,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(1,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(1,0,tile))
          CALL get_bounds (ng, tile, p2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(1,1,tile),                 &
     &                     BOUNDS(ng) % Imax(1,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(1,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(1,1,tile))

          CALL get_bounds (ng, tile, r2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(2,0,tile),                 &
     &                     BOUNDS(ng) % Imax(2,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(2,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(2,0,tile))
          CALL get_bounds (ng, tile, r2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(2,1,tile),                 &
     &                     BOUNDS(ng) % Imax(2,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(2,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(2,1,tile))

          CALL get_bounds (ng, tile, u2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(3,0,tile),                 &
     &                     BOUNDS(ng) % Imax(3,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(3,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(3,0,tile))
          CALL get_bounds (ng, tile, u2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(3,1,tile),                 &
     &                     BOUNDS(ng) % Imax(3,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(3,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(3,1,tile))

          CALL get_bounds (ng, tile, v2dvar, 0     , Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(4,0,tile),                 &
     &                     BOUNDS(ng) % Imax(4,0,tile),                 &
     &                     BOUNDS(ng) % Jmin(4,0,tile),                 &
     &                     BOUNDS(ng) % Jmax(4,0,tile))
          CALL get_bounds (ng, tile, v2dvar, Nghost, Itile, Jtile,      &
     &                     BOUNDS(ng) % Imin(4,1,tile),                 &
     &                     BOUNDS(ng) % Imax(4,1,tile),                 &
     &                     BOUNDS(ng) % Jmin(4,1,tile),                 &
     &                     BOUNDS(ng) % Jmax(4,1,tile))
        END DO
      END DO
!
!  Set NetCDF IO bounds.
!
      DO ng=1,Ngrids
        CALL get_iobounds (ng)
      END DO
!
!-----------------------------------------------------------------------
!  Check tile partition starting and ending (I,J) indices for illegal
!  domain decomposition parameters NtileI and NtileJ in standard input
!  file.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
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
     &          i4.4,'x',i4.4,'x',i4.4,2x,'tiling: ',i3.3,'x',i3.3,/,/, &
     &          5x,'tile',5x,'Istr',5x,'Iend',5x,'Jstr',5x,'Jend',      &
     &          5x,'Npts',/)
#else
 50     FORMAT (/,' Tile partition information for Grid ',i2.2,':',2x,  &
     &          i4.4,'x',i4.4,2x,'tiling: ',i3.3,'x',i3.3,/,/,          &
     &          5x,'tile',5x,'Istr',5x,'Iend',5x,'Jstr',5x,'Jend',      &
     &          5x,'Npts',/)
#endif
#if !defined DISTRIBUTE && defined ADJOINT
 60     FORMAT (/,' INP_PAR - illegal domain decomposition for the ',   &
     &                       'Adjoint model.',/,11x,'Partitions are ',  &
     &          'allowed in distributed-menory (MPI) applications.'/)
#endif
 70     FORMAT (5(5x,i4),2x,i7)
 80     FORMAT (/,' INP_PAR - domain decomposition error in input ',    &
     &                        'script file for grid: ',i2,/,            &
     &          /,11x,'The domain partition parameter, ',a,i3,          &
     &          /,11x,'is incompatible with grid size, ',a,i4,          &
     &          /,11x,'because it yields too small tile, ',a,i3,a,i3,   &
     &          /,11x,'Decrease partition parameter: ',a)
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
#endif
      IF (exit_flag.ne.NoError) RETURN
!
!  Report tile minimum and maximum fractional grid coordinates.
!
      DO ng=1,Ngrids
        IF (Master) THEN
          WRITE (stdout,90)
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
 90       FORMAT (/,' Tile minimum and maximum fractional grid ',       &
     &            'coordinates:',/,                                     &
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
# if defined EW_PERIODIC || defined NS_PERIODIC
      Nghost=1+GHOST_POINTS
# else
      Nghost=GHOST_POINTS
# endif
      DO ng=1,Ngrids
        MaxHaloLenI=0
        MaxHaloLenJ=0
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
        IF (Master) THEN
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

#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Read in biological model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (15, FILE=TRIM(bparnam), FORM='formatted', STATUS='old')

      CALL read_BioPar (model, 15, out, Lwrite)
#endif
#ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Read in sediment model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (25, FILE=TRIM(sparnam), FORM='formatted', STATUS='old')

      CALL read_SedPar (model, 25, out, Lwrite)
#endif
#if defined ASSIMILATION || defined FOUR_DVAR || defined NUDGING || \
    defined VERIFICATION
!
!-----------------------------------------------------------------------
!  Read in input assimilation parameters.
!-----------------------------------------------------------------------
!
      OPEN (35, FILE=TRIM(aparnam), FORM='formatted', STATUS='old')

      CALL read_AssPar (model, 35, out, Lwrite)
#endif
#ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Read in floats input parameters.
!-----------------------------------------------------------------------
!
      OPEN (45, FILE=TRIM(fposnam), FORM='formatted', STATUS='old')

      CALL read_FloatsPar (model, 45, out, Lwrite)
#endif
#ifdef STATIONS
!
!-----------------------------------------------------------------------
!  Read in stations input parameters.
!-----------------------------------------------------------------------
!
      OPEN (55, FILE=TRIM(sposnam), FORM='formatted', STATUS='old')

      CALL read_StaPar (model, 55, out, Lwrite)
#endif
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options and definitions.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        CALL checkdefs
        CALL my_flush (out)
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
      CALL mp_bcasts (1, model, Coptions)
#endif
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Compute various constants.
!-----------------------------------------------------------------------
!
      gorho0=g/rho0
      DO ng=1,Ngrids
        dtfast(ng)=dt(ng)/REAL(ndtfast(ng),r8)
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
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
        IF (Znudg(ng).gt.0.0_r8) THEN
          Znudg(ng)=1.0_r8/(Znudg(ng)*86400.0_r8)
        ELSE
          Znudg(ng)=0.0_r8
        END IF
        IF (M2nudg(ng).gt.0.0_r8) THEN
          M2nudg(ng)=1.0_r8/(M2nudg(ng)*86400.0_r8)
        ELSE
          M2nudg(ng)=0.0_r8
        END IF
#ifdef SOLVE3D
        IF (M3nudg(ng).gt.0.0_r8) THEN
          M3nudg(ng)=1.0_r8/(M3nudg(ng)*86400.0_r8)
        ELSE
          M3nudg(ng)=0.0_r8
        END IF
#endif
#ifdef SO_SEMI
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
      END DO
!
!-----------------------------------------------------------------------
!  Initialize random number sequence so we can get identical results
!  everytime that we run the same solution.
!-----------------------------------------------------------------------
!
      sequence=759
      CALL ran_seed (sequence)

      RETURN
      END SUBROUTINE inp_par

#ifdef BIOLOGY
# if defined BIO_FENNEL
#  include <fennel_inp.h>
# elif defined ECOSIM
#  include <ecosim_inp.h>
# elif defined NEMURO
#  include <nemuro_inp.h>
# elif defined NPZD_FRANKS
#  include <npzd_Franks_inp.h>
# elif defined NPZD_IRON
#  include <npzd_iron_inp.h>
# elif defined NPZD_POWELL
#  include <npzd_Powell_inp.h>
# endif
#endif

#ifdef SEDIMENT
# include <sediment_inp.h>
#endif

#ifdef MODEL_COUPLING
      SUBROUTINE read_CouplePar (model)
!
!=======================================================================
!                                                                      !
!  This routine reads in physical model input parameters.              !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasts
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      logical :: Lwrite

      integer :: Npts, Nval, i, ic, j, inp, ng, out, status

      integer :: decode_line, load_i, load_l, load_r

      real(r8), dimension(100) :: Rval

!     real(r8), allocatable :: MyRval(:)

      character (len=1 ), parameter :: blank = ' '
      character (len=40) :: KeyWord
      character (len=80) :: Cname
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Determine coupling standard input file name.  In distributed-memory,
!  this name is assigned at the executtion command line and processed
!  with the Unix routine GETARG.  The ROMS/TOMS input parameter script
!  name is specified in this coupling script.
!-----------------------------------------------------------------------
!
# ifdef DISTRIBUTE
      Lwrite=Master
      inp=1
      out=stdout
!
      IF (MyRank.eq.0) CALL my_getarg (1,Cname)
      DO ng=1,Ngrids
        CALL mp_bcasts (ng, model, Cname)
      END DO
      IF (MyRank.eq.0) THEN
        WRITE(stdout,*) 'Coupled Input File name = ', TRIM(Cname)
      END IF
      OPEN (inp, FILE=TRIM(Cname), FORM='formatted', STATUS='old',      &
     &      ERR=10)
      GO TO 30
 10   WRITE (stdout,20)
      IF (Master) WRITE(stdout,*) 'MyRank = ', MyRank, TRIM(Cname)
      exit_flag=4
      RETURN
 20   FORMAT (/,' INP_PAR - Unable to open coupling input script.',     &
     &        /,11x,'In distributed-memory applications, the input',    &
     &        /,11x,'script file is processed in parallel. The Unix',   &
     &        /,11x,'routine GETARG is used to get script file name.',  &
     &        /,11x,'For example, in MPI applications make sure that',  &
     &        /,11x,'command line is something like:',/,                &
     &        /,11x,'mpirun -np 4 masterM coupling.in',/,               &
     &        /,11x,'and not',/,                                        &
     &        /,11x,'mpirun -np 4 masterM < coupling.in',/)
 30   CONTINUE
# else
      Lwrite=Master
      inp=stdinp
      out=stdout
# endif
!
!-----------------------------------------------------------------------
!  Read in multiple models coupling parameters. Then, load input
!  data into module. Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=40,END=50) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'VARNAME') THEN
            DO i=1,LEN(USRname)
              varname(i:i)=blank
            END DO
            varname=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NnodesATM') THEN
            Npts=load_i(Nval, Rval, 1, NnodesATM)
          ELSE IF (TRIM(KeyWord).eq.'NnodesWAV') THEN
            Npts=load_i(Nval, Rval, 1, NnodesWAV)
          ELSE IF (TRIM(KeyWord).eq.'NnodesOCN') THEN
            Npts=load_i(Nval, Rval, 1, NnodesOCN)
          ELSE IF (TRIM(KeyWord).eq.'TI_ATM_WAV') THEN
            Npts=load_r(Nval, Rval, 1, TI_ATM_WAV)
          ELSE IF (TRIM(KeyWord).eq.'TI_ATM_OCN') THEN
            Npts=load_r(Nval, Rval, 1, TI_ATM_OCN)
          ELSE IF (TRIM(KeyWord).eq.'TI_WAV_OCN') THEN
            Npts=load_r(Nval, Rval, 1, TI_WAV_OCN)
          ELSE IF (TRIM(KeyWord).eq.'WRF_CPL_GRID') THEN
            Npts=load_i(Nval, Rval, 1, WRF_CPL_GRID)
          ELSE IF (TRIM(KeyWord).eq.'OCN_name') THEN
            DO i=1,LEN(Iname)
              Iname(i:i)=blank
            END DO
            Iname=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'WAV_name') THEN
            DO i=1,LEN(Wname(Nval))
              Wname(Nval)(i:i)=blank
            END DO
            Wname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ATM_name') THEN
            DO i=1,LEN(Aname)
              Aname(i:i)=blank
            END DO
            Aname=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'O2ANAME') THEN
            DO i=1,LEN(O2Aname(Nval))
              O2Aname(Nval)(i:i)=blank
            END DO
            O2Aname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'O2WNAME') THEN
            DO i=1,LEN(O2Wname(Nval))
              O2Wname(Nval)(i:i)=blank
            END DO
            O2Wname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'A2ONAME') THEN
            DO i=1,LEN(A2Oname(Nval))
              A2Oname(Nval)(i:i)=blank
            END DO
            A2Oname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'A2WNAME') THEN
            DO i=1,LEN(A2Wname(Nval))
              A2Wname(Nval)(i:i)=blank
            END DO
            A2Wname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'W2ANAME') THEN
            DO i=1,LEN(W2Aname(Nval))
              W2Aname(Nval)(i:i)=blank
            END DO
            W2Aname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'W2ONAME') THEN
            DO i=1,LEN(W2Oname(Nval))
              W2Oname(Nval)(i:i)=blank
            END DO
            W2Oname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
  40  IF (Master) WRITE (out,60) line
      exit_flag=4
      RETURN
  50  CLOSE (inp)

  60  FORMAT (/,' READ_CouplePar - Error while processing line: ',/,a)

      RETURN
      END SUBROUTINE read_CouplePar
#endif

      SUBROUTINE read_PhyPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in physical model input parameters.              !
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
      USE mct_coupler_params
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
#ifdef PROPAGATOR
      USE mod_storage
#endif
      USE mod_strings
#ifdef INWAVE_COUPLING
      USE mod_inwave_params
#endif
#ifdef WAVES_OCEAN
      USE mct_coupler_params
#endif
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
      logical :: inhere
#if defined SOLVE3D && defined SEDIMENT
      logical :: LreadNCS = .FALSE.
      logical :: LreadNNS = .FALSE.
#endif
      integer :: Lstr, Npts, Nval, i, itrc, k, ng, status

      integer :: decode_line, load_i, load_l, load_r

#if defined SOLVE3D && defined T_PASSIVE
      logical, allocatable :: Linert(:,:)
#endif
#if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
      logical, dimension(MBOTP,Ngrids) :: Lbottom
#endif
      logical, allocatable :: Ltracer(:,:)
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
      logical, allocatable :: Ladsen(:)
#endif

      real(r8), allocatable :: Rtracer(:,:)
      real(r8), allocatable :: tracer(:,:)

      real(r8), dimension(100) :: Rval

      character (len=1 ), parameter :: blank = ' '
      character (len=19) :: ref_att
      character (len=40) :: KeyWord
      character (len=160) :: fname, line
      character (len=160), dimension(100) :: Cval
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
          IF (TRIM(KeyWord).eq.'TITLE') THEN
            IF (Nval.eq.1) THEN
              title=TRIM(ADJUSTL(Cval(Nval)))
            ELSE
              WRITE(title,'(a,1x,a)') TRIM(ADJUSTL(title)),             &
     &                                TRIM(ADJUSTL(Cval(Nval)))
            END IF
          ELSE IF (TRIM(KeyWord).eq.'MyAppCPP') THEN
            DO i=1,LEN(MyAppCPP)
              MyAppCPP(i:i)=blank
            END DO
            MyAppCPP=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'VARNAME') THEN
            DO i=1,LEN(varname)
              varname(i:i)=blank
            END DO
            varname=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'Lm') THEN
            Npts=load_i(Nval, Rval, Ngrids, Lm)
            DO ng=1,Ngrids
              IF (Lm(ng).le.0) THEN
                IF (Master) WRITE (out,300) 'Lm', ng,                   &
     &                                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Mm') THEN
            Npts=load_i(Nval, Rval, Ngrids, Mm)
            DO ng=1,Ngrids
              IF (Mm(ng).le.0) THEN
                IF (Master) WRITE (out,300) 'Mm', ng,                   &
     &                                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'N') THEN
            Npts=load_i(Nval, Rval, Ngrids, N)
            DO ng=1,Ngrids
              IF (N(ng).lt.0) THEN
                IF (Master) WRITE (out,300) 'N', ng,                    &
     &                                      'must be greater than zero.'
                exit_flag=5
                RETURN
              END IF
            END DO
#ifdef INWAVE_COUPLING
          ELSE IF (TRIM(KeyWord).eq.'ND') THEN
            Npts=load_i(Nval, Rval, Ngrids, ND)
            DO ng=1,Ngrids
              IF (ND(ng).lt.0) THEN
                WRITE (out,300) 'ND', ng,                               &
     &                          'must be greater than zero.'
                STOP
              END IF
            END DO
#endif
#if defined SEDIMENT && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Nbed') THEN
            Npts=load_i(Nval, Rval, 1, Nbed)
            IF (Nbed.le.0) THEN
              IF (Master) WRITE (out,290) 'Nbed = ', Nbed,              &
     &                                    'must be greater than zero.'
              exit_flag=5
              RETURN
            END IF
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'NAT') THEN
            Npts=load_i(Nval, Rval, 1, NAT)
            IF ((NAT.lt.1).or.(NAT.gt.2)) THEN
              IF (Master) WRITE (out,290) 'NAT = ', NAT,                &
     &                          'make sure that NAT is either 1 or 2.'
              exit_flag=5
              RETURN
            END IF
# ifdef SALINITY
            IF (NAT.ne.2) THEN
              IF (Master) WRITE (out,290) 'NAT = ', NAT,                &
     &                          'make sure that NAT is equal to 2.'
              exit_flag=5
              RETURN
            END IF
# endif
#endif
#if defined T_PASSIVE && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'NPT') THEN
            Npts=load_i(Nval, Rval, 1, NPT)
            IF (NPT.le.0) THEN
              IF (Master) WRITE (out,290) 'NPT = ', NPT,                &
     &                                    'must be greater than zero.'
              exit_flag=5
              RETURN
            END IF
#endif
#if defined SEDIMENT && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'NCS') THEN
            Npts=load_i(Nval, Rval, 1, NCS)
            IF (NCS.lt.0) THEN
              IF (Master) WRITE (out,290) 'NCS = ', NCS,                &
     &                                    'must be greater than zero.'
              exit_flag=5
              RETURN
            END IF
            LreadNCS=.TRUE.
            IF (LreadNNS.and.((NCS+NNS).le.0)) THEN
              IF (Master) WRITE (out,290) 'NST = ', NCS+NNS,            &
     &                    'either NCS or NNS must be greater than zero.'
              exit_flag=5
              RETURN
            END IF
            NST=NST+NCS
          ELSE IF (TRIM(KeyWord).eq.'NNS') THEN
            Npts=load_i(Nval, Rval, 1, NNS)
            IF (NNS.lt.0) THEN
              IF (Master) WRITE (out,290) 'NNS = ',                     &
     &                                    'must be greater than zero.'
              exit_flag=5
              RETURN
            END IF
            LreadNNS=.TRUE.
            IF (LreadNCS.and.((NCS+NNS).le.0)) THEN
              IF (Master) WRITE (out,290) 'NST = ', NCS+NNS,            &
     &                    'either NCS or NNS must be greater than zero.'
              exit_flag=5
              RETURN
            END IF
            NST=NST+NNS
#endif
          ELSE IF (TRIM(KeyWord).eq.'NtileI') THEN
            Npts=load_i(Nval, Rval, Ngrids, NtileI)
#ifdef DISTRIBUTE
            NtileX(1:Ngrids)=1
#else
            NtileX(1:Ngrids)=NtileI(1:Ngrids)
#endif
          ELSE IF (TRIM(KeyWord).eq.'NtileJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, NtileJ)
#ifdef DISTRIBUTE
            NtileE(1:Ngrids)=1
#else
            NtileE(1:Ngrids)=NtileJ(1:Ngrids)
#endif
#ifdef BIOLOGY
            CALL initialize_biology
#endif
#if defined SEDIMENT || defined BBL_MODEL
            CALL initialize_sediment
#endif
            CALL initialize_param
            CALL initialize_scalars
            CALL initialize_ncparam
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
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
#endif
            IF (.not.allocated(Rtracer)) THEN
              allocate (Rtracer(NAT+NPT,Ngrids))
            END IF
            IF (.not.allocated(tracer)) THEN
              allocate (tracer(MT,Ngrids))
            END IF
          ELSE IF (TRIM(KeyWord).eq.'NTIMES') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntimes)
          ELSE IF (TRIM(KeyWord).eq.'DT') THEN
            Npts=load_r(Nval, Rval, Ngrids, dt)
#ifdef AIR_OCEAN
            DO ng=1,Ngrids
              nOCN_ATM(ng)=MAX(1,INT(TI_ATM_OCN/dt(ng)))
            END DO
#endif
#ifdef WAVES_OCEAN
            DO ng=1,Ngrids
              nOCN_WAV(ng)=MAX(1,INT(TI_WAV_OCN/dt(ng)))
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'NDTFAST') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndtfast)
          ELSE IF (TRIM(KeyWord).eq.'ERstr') THEN
            Npts=load_i(Nval, Rval, 1, ERstr)
          ELSE IF (TRIM(KeyWord).eq.'ERend') THEN
            Npts=load_i(Nval, Rval, 1, ERend)
          ELSE IF (TRIM(KeyWord).eq.'Nouter') THEN
            Npts=load_i(Nval, Rval, 1, Nouter)
          ELSE IF (TRIM(KeyWord).eq.'Ninner') THEN
            Npts=load_i(Nval, Rval, 1, Ninner)
          ELSE IF (TRIM(KeyWord).eq.'Nintervals') THEN
            Npts=load_i(Nval, Rval, 1, Nintervals)
#ifdef PROPAGATOR
          ELSE IF (TRIM(KeyWord).eq.'NEV') THEN
            Npts=load_i(Nval, Rval, 1, NEV)
          ELSE IF (TRIM(KeyWord).eq.'NCV') THEN
            Npts=load_i(Nval, Rval, 1, NCV)
# if defined FT_EIGENMMODES || defined AFT_EIGENMODES
            IF (NCV.lt.(2*NEV+1)) THEN
              IF (Master) WRITE (out,260) 'NCV',                        &
     &                       'Must be greater than or equal to 2*NEV+1'
              exit_flag=5
              RETURN
            END IF
# elif defined OPT_PERTURBATION
            IF (NCV.lt.(2*NEV)) THEN
              IF (Master) WRITE (out,260) 'NCV',                        &
     &                       'Must be greater than or equal to 2*NEV'
              exit_flag=5
              RETURN
            END IF
# else
            IF (NCV.lt.(2*NEV)) THEN
              IF (Master) WRITE (out,260) 'NCV',                        &
     &                                    'Must be greater than NEV'
              exit_flag=5
              RETURN
            END IF
# endif
#endif
          ELSE IF (TRIM(KeyWord).eq.'NRREC') THEN
            Npts=load_i(Nval, Rval, Ngrids, nrrec)
            DO ng=1,Ngrids
              IF (nrrec(ng).lt.0) THEN
                LastRec(ng)=.TRUE.
              ELSE
                LastRec(ng)=.FALSE.
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'LcycleRST') THEN
            Npts=load_l(Nval, Cval, Ngrids, LcycleRST)
          ELSE IF (TRIM(KeyWord).eq.'NRST') THEN
            Npts=load_i(Nval, Rval, Ngrids, nRST)
          ELSE IF (TRIM(KeyWord).eq.'NSTA') THEN
            Npts=load_i(Nval, Rval, Ngrids, nSTA)
          ELSE IF (TRIM(KeyWord).eq.'NFLT') THEN
            Npts=load_i(Nval, Rval, Ngrids, nFLT)
          ELSE IF (TRIM(KeyWord).eq.'NINFO') THEN
            Npts=load_i(Nval, Rval, Ngrids, ninfo)
          ELSE IF (TRIM(KeyWord).eq.'LDEFOUT') THEN
            Npts=load_l(Nval, Cval, Ngrids, ldefout)
          ELSE IF (TRIM(KeyWord).eq.'NHIS') THEN
            Npts=load_i(Nval, Rval, Ngrids, nHIS)
          ELSE IF (TRIM(KeyWord).eq.'NDEFHIS') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefHIS)
          ELSE IF (TRIM(KeyWord).eq.'NTSAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntsAVG)
#ifdef ADJOINT
            DO ng=1,Ngrids
              IF (ntsAVG(ng).eq.1) ntsAVG(ng)=ntimes(ng)
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'NAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, nAVG)
          ELSE IF (TRIM(KeyWord).eq.'NDEFAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefAVG)
          ELSE IF (TRIM(KeyWord).eq.'NTSDIA') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntsDIA)
          ELSE IF (TRIM(KeyWord).eq.'NDIA') THEN
            Npts=load_i(Nval, Rval, Ngrids, nDIA)
          ELSE IF (TRIM(KeyWord).eq.'NDEFDIA') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefDIA)
          ELSE IF (TRIM(KeyWord).eq.'LcycleTLM') THEN
            Npts=load_l(Nval, Cval, Ngrids, LcycleTLM)
          ELSE IF (TRIM(KeyWord).eq.'NTLM') THEN
            Npts=load_i(Nval, Rval, Ngrids, nTLM)
          ELSE IF (TRIM(KeyWord).eq.'NDEFTLM') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefTLM)
          ELSE IF (TRIM(KeyWord).eq.'LcycleADJ') THEN
            Npts=load_l(Nval, Cval, Ngrids, LcycleADJ)
          ELSE IF (TRIM(KeyWord).eq.'NADJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, nADJ)
          ELSE IF (TRIM(KeyWord).eq.'NDEFADJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefADJ)
          ELSE IF (TRIM(KeyWord).eq.'NOBC') THEN
            Npts=load_i(Nval, Rval, Ngrids, nOBC)
#ifdef ADJUST_BOUNDARY
            DO ng=1,Ngrids
              Nbrec(ng)=1+ntimes(ng)/nOBC(ng)
            END DO
            allocate ( OBC_time(MAXVAL(Nbrec),Ngrids) )
#endif
          ELSE IF (TRIM(KeyWord).eq.'NSFF') THEN
            Npts=load_i(Nval, Rval, Ngrids, nSFF)
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
            DO ng=1,Ngrids
              Nfrec(ng)=1+ntimes(ng)/nSFF(ng)
            END DO
            allocate ( SF_time(MAXVAL(Nfrec),Ngrids) )
#endif
          ELSE IF (TRIM(KeyWord).eq.'LrstGST') THEN
            Npts=load_l(Nval, Cval, 1, LrstGST)
          ELSE IF (TRIM(KeyWord).eq.'MaxIterGST') THEN
            Npts=load_i(Nval, Rval, 1, MaxIterGST)
          ELSE IF (TRIM(KeyWord).eq.'NGST') THEN
            Npts=load_i(Nval, Rval, 1, nGST)
#ifdef PROPAGATOR
          ELSE IF (TRIM(KeyWord).eq.'Ritz_tol') THEN
            Npts=load_r(Nval, Rval, 1, Ritz_tol)
#endif
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                nl_tnu2(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                nl_tnu4(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_TNU2') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                ad_tnu2(itrc,ng)=Rtracer(itrc,ng)
                tl_tnu2(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_TNU4') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                ad_tnu4(itrc,ng)=Rtracer(itrc,ng)
                tl_tnu4(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VISC2') THEN
            Npts=load_r(Nval, Rval, Ngrids, nl_visc2)
          ELSE IF (TRIM(KeyWord).eq.'VISC4') THEN
            Npts=load_r(Nval, Rval, Ngrids, nl_visc4)
          ELSE IF (TRIM(KeyWord).eq.'ad_VISC2') THEN
            Npts=load_r(Nval, Rval, Ngrids, ad_visc2)
            DO ng=1,Ngrids
              tl_visc2(ng)=ad_visc2(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VISC4') THEN
            Npts=load_r(Nval, Rval, Ngrids, ad_visc4)
            DO ng=1,Ngrids
              tl_visc4(ng)=ad_visc4(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                Akt_bak(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ad_AKT_fac') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                ad_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
                tl_Akt_fac(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKV_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
          ELSE IF (TRIM(KeyWord).eq.'ad_AKV_fac') THEN
            Npts=load_r(Nval, Rval, Ngrids, ad_Akv_fac)
            DO ng=1,Ngrids
              tl_Akv_fac(ng)=ad_AKv_fac(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKK_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akk_bak)
          ELSE IF (TRIM(KeyWord).eq.'AKP_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akp_bak)
          ELSE IF (TRIM(KeyWord).eq.'TKENU2') THEN
            Npts=load_r(Nval, Rval, Ngrids, tkenu2)
          ELSE IF (TRIM(KeyWord).eq.'TKENU4') THEN
            Npts=load_r(Nval, Rval, Ngrids, tkenu4)
          ELSE IF (TRIM(KeyWord).eq.'GLS_P') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_p)
          ELSE IF (TRIM(KeyWord).eq.'GLS_M') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_m)
          ELSE IF (TRIM(KeyWord).eq.'GLS_N') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_n)
          ELSE IF (TRIM(KeyWord).eq.'GLS_Kmin') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_Kmin)
          ELSE IF (TRIM(KeyWord).eq.'GLS_Pmin') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_Pmin)
          ELSE IF (TRIM(KeyWord).eq.'GLS_CMU0') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_cmu0)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C1') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c1)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C2') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c2)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C3M') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c3m)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C3P') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c3p)
          ELSE IF (TRIM(KeyWord).eq.'GLS_SIGK') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_sigk)
          ELSE IF (TRIM(KeyWord).eq.'GLS_SIGP') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_sigp)
          ELSE IF (TRIM(KeyWord).eq.'CHARNOK_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, charnok_alpha)
          ELSE IF (TRIM(KeyWord).eq.'ZOS_HSIG_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, zos_hsig_alpha)
          ELSE IF (TRIM(KeyWord).eq.'SZ_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, sz_alpha)
          ELSE IF (TRIM(KeyWord).eq.'CRGBAN_CW') THEN
            Npts=load_r(Nval, Rval, Ngrids, crgban_cw)
          ELSE IF (TRIM(KeyWord).eq.'WEC_ALPHA') THEN
            Npts=load_r(Nval, Rval, Ngrids, wec_alpha)
          ELSE IF (TRIM(KeyWord).eq.'RDRG') THEN
            Npts=load_r(Nval, Rval, Ngrids, rdrg)
          ELSE IF (TRIM(KeyWord).eq.'RDRG2') THEN
            Npts=load_r(Nval, Rval, Ngrids, rdrg2)
          ELSE IF (TRIM(KeyWord).eq.'Zob') THEN
            Npts=load_r(Nval, Rval, Ngrids, Zob)
          ELSE IF (TRIM(KeyWord).eq.'Zos') THEN
            Npts=load_r(Nval, Rval, Ngrids, Zos)
#ifdef BULK_FLUXES
          ELSE IF (TRIM(KeyWord).eq.'BLK_ZQ') THEN
            Npts=load_r(Nval, Rval, Ngrids, blk_ZQ)
          ELSE IF (TRIM(KeyWord).eq.'BLK_ZT') THEN
            Npts=load_r(Nval, Rval, Ngrids, blk_ZT)
          ELSE IF (TRIM(KeyWord).eq.'BLK_ZW') THEN
            Npts=load_r(Nval, Rval, Ngrids, blk_ZW)
#endif
          ELSE IF (TRIM(KeyWord).eq.'DCRIT') THEN
            Npts=load_r(Nval, Rval, Ngrids, Dcrit)
          ELSE IF (TRIM(KeyWord).eq.'WTYPE') THEN
            Npts=load_i(Nval, Rval, Ngrids, lmd_Jwt)
          ELSE IF (TRIM(KeyWord).eq.'LEVSFRC') THEN
            Npts=load_i(Nval, Rval, Ngrids, levsfrc)
          ELSE IF (TRIM(KeyWord).eq.'LEVBFRC') THEN
            Npts=load_i(Nval, Rval, Ngrids, levbfrc)
          ELSE IF (TRIM(KeyWord).eq.'Vtransform') THEN
            Npts=load_i(Nval, Rval, Ngrids, Vtransform)
            DO ng=1,Ngrids
              IF ((Vtransform(ng).lt.0).or.                               &
     &            (Vtransform(ng).gt.2)) THEN
                IF (Master) WRITE (out,260) 'Vtransform = ',              &
     &                                      Vtransform(ng),               &
     &                                      'Must be either 1 or 2'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Vstretching') THEN
            Npts=load_i(Nval, Rval, Ngrids, Vstretching)
            DO ng=1,Ngrids
              IF ((Vstretching(ng).lt.0).or.                              &
     &            (Vstretching(ng).gt.3)) THEN
                IF (Master) WRITE (out,260) 'Vstretching = ',             &
     &                                      Vstretching(ng),              &
     &                                      'Must between 1 and 3'
                exit_flag=5
                RETURN
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'THETA_S') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_s)
          ELSE IF (TRIM(KeyWord).eq.'THETA_B') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_b)
          ELSE IF (TRIM(KeyWord).eq.'TCLINE') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcline)
            DO ng=1,Ngrids
              hc(ng)=Tcline(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'RHO0') THEN
            Npts=load_r(Nval, Rval, 1, rho0)
          ELSE IF (TRIM(KeyWord).eq.'BVF_BAK') THEN
            Npts=load_r(Nval, Rval, 1, bvf_bak)
          ELSE IF (TRIM(KeyWord).eq.'DSTART') THEN
            Npts=load_r(Nval, Rval, 1, dstart)
          ELSE IF (TRIM(KeyWord).eq.'TIDE_START') THEN
            Npts=load_r(Nval, Rval, 1, tide_start)
          ELSE IF (TRIM(KeyWord).eq.'TIME_REF') THEN
            Npts=load_r(Nval, Rval, 1, time_ref)
            r_text=ref_att(time_ref,r_date)
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                Tnudg(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZNUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, Znudg)
          ELSE IF (TRIM(KeyWord).eq.'M2NUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, M2nudg)
          ELSE IF (TRIM(KeyWord).eq.'M3NUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, M3nudg)
          ELSE IF (TRIM(KeyWord).eq.'OBCFAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, obcfac)
          ELSE IF (TRIM(KeyWord).eq.'R0') THEN
            Npts=load_r(Nval, Rval, Ngrids, R0)
            DO ng=1,Ngrids
              IF (R0(ng).lt.100.0_r8) R0(ng)=R0(ng)+1000.0_r8
            END DO
          ELSE IF (TRIM(KeyWord).eq.'T0') THEN
            Npts=load_r(Nval, Rval, Ngrids, T0)
          ELSE IF (TRIM(KeyWord).eq.'S0') THEN
            Npts=load_r(Nval, Rval, Ngrids, S0)
          ELSE IF (TRIM(KeyWord).eq.'TCOEF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcoef)
            DO ng=1,Ngrids
              Tcoef(ng)=ABS(Tcoef(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SCOEF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Scoef)
            DO ng=1,Ngrids
              Scoef(ng)=ABS(Scoef(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'GAMMA2') THEN
            Npts=load_r(Nval, Rval, Ngrids, gamma2)
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
          ELSE IF (TRIM(KeyWord).eq.'DstrS') THEN
            Npts=load_r(Nval, Rval, Ngrids, DstrS)
          ELSE IF (TRIM(KeyWord).eq.'DendS') THEN
            Npts=load_r(Nval, Rval, Ngrids, DendS)
          ELSE IF (TRIM(KeyWord).eq.'KstrS') THEN
            Npts=load_i(Nval, Rval, Ngrids, KstrS)
          ELSE IF (TRIM(KeyWord).eq.'KendS') THEN
            Npts=load_i(Nval, Rval, Ngrids, KendS)
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isFsur)') THEN
            IF (isFsur.eq.0) THEN
              IF (Master) WRITE (out,280) 'isFsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isFsur)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isUbar)') THEN
            IF (isUbar.eq.0) THEN
              IF (Master) WRITE (out,280) 'isUbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isUbar)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isVbar)') THEN
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
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isUvel)') THEN
            IF (isUvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'isUvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isUvel)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isVvel)') THEN
            IF (isVvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'isVvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isVvel)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isTvar)') THEN
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
# ifdef SO_SEMI
          ELSE IF (TRIM(KeyWord).eq.'SO_decay') THEN
            Npts=load_r(Nval, Rval, Ngrids, SO_decay)
          ELSE IF (TRIM(KeyWord).eq.'SOstate(isUstr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%SOstate(isUstr)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SOstate(isVstr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%SOstate(isVstr)=Ladsen(ng)
            END DO
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'SOstate(isTsur)') THEN
            Npts=load_l(Nval, Cval, MT*Ngrids, Ladsen)
            k=0
            DO ng=1,Ngrids
              DO itrc=1,MT
                k=k+1
                i=isTsur(itrc)
                SCALARS(ng)%SOstate(i)=Ladsen(k)
              END DO
            END DO
#  endif
          ELSE IF (TRIM(KeyWord).eq.'SO_sdev(isUstr)') THEN
            IF (isUstr.eq.0) THEN
              IF (Master) WRITE (out,280) 'isUstr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isUstr,1))
          ELSE IF (TRIM(KeyWord).eq.'SO_sdev(isVstr)') THEN
            IF (isUstr.eq.0) THEN
              IF (Master) WRITE (out,280) 'isVstr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isVstr,1))
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'SO_sdev(isTsur)') THEN
            IF (MAXVAL(isTsur).eq.0) THEN
              IF (Master) WRITE (out,280) 'isTsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_r(Nval, Rval, MT*Ngrids, tracer)
            k=0
            DO ng=1,Ngrids
              DO itrc=1,MT
                k=k+1
                i=isTsur(itrc)
                SO_sdev(i,ng)=tracer(k,ng)
              END DO
            END DO
#  endif
# endif
#endif
#if defined SOLVE3D && defined TS_PSOURCE
          ELSE IF (TRIM(KeyWord).eq.'LtracerSrc') THEN
            Npts=load_l(Nval, Cval, (NAT+NPT)*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                LtracerSrc(itrc,ng)=Ltracer(itrc,ng)
              END DO
            END DO
# ifdef T_PASSIVE
            IF (MAXVAL(inert).eq.0) THEN
              IF (Master) WRITE (out,280) 'inert'
              exit_flag=5
              RETURN
            END IF
            DO ng=1,Ngrids
              DO itrc=1,NPT
                i=inert(itrc)
                LtracerSrc(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
# endif
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUvel)') THEN
            IF (idUvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvel)') THEN
            IF (idVvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWvel)') THEN
            IF (idWvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idOvel)') THEN
            IF (idOvel.eq.0) THEN
              IF (Master) WRITE (out,280) 'idOvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idOvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbar)') THEN
            IF (idUbar.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbar,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbar)') THEN
            IF (idVbar.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbar,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idFsur)') THEN
            IF (idFsur.eq.0) THEN
              IF (Master) WRITE (out,280) 'idFsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idFsur,1:Ngrids))
#if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(KeyWord).eq.'Hout(idBath)') THEN
            IF (idbath.eq.0) THEN
              IF (Master) WRITE (out,280) 'idbath'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idBath,1:Ngrids))
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
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
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUair)') THEN
            IF (idUair.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUair'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUair,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVair)') THEN
            IF (idVair.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVair'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVair,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUsms)') THEN
            IF (idUsms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUsms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUsms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVsms)') THEN
            IF (idVsms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVsms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVsms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbms)') THEN
            IF (idUbms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbms)') THEN
            IF (idVbms.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbms'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbrs)') THEN
            IF (idUbrs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbrs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbrs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbrs)') THEN
            IF (idVbrs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbrs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbrs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbws)') THEN
            IF (idUbws.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbws'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbws,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbws)') THEN
            IF (idVbws.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbws'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbws,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbcs)') THEN
            IF (idUbcs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbcs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbcs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbcs)') THEN
            IF (idVbcs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbcs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbcs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUVwc)') THEN
            IF (idUVwc.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUVwc'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUVwc,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbot)') THEN
            IF (idUbot.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbot'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbot,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbot)') THEN
            IF (idVbot.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbot'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbot,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbur)') THEN
            IF (idUbur.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUbur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbur,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbvr)') THEN
            IF (idVbvr.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVbvr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbvr,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2xx)') THEN
            IF (idW2xx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW2xx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xx,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2xy)') THEN
            IF (idW2xy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW2xy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2yy)') THEN
            IF (idW2yy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW2yy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2yy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU2rs)') THEN
            IF (idU2rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU2rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU2rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV2rs)') THEN
            IF (idV2rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV2rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV2rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU2Sd)') THEN
            IF (idU2Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU2Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU2Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV2Sd)') THEN
            IF (idV2Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV2Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV2Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3xx)') THEN
            IF (idW3xx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3xx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xx,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3xy)') THEN
            IF (idW3xy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3xy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3yy)') THEN
            IF (idW3yy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3yy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3yy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3zx)') THEN
            IF (idW3zx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3zx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zx,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3zy)') THEN
            IF (idW3zy.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3zy'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU3rs)') THEN
            IF (idU3rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU3rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU3rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV3rs)') THEN
            IF (idV3rs.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV3rs'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV3rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU3Sd)') THEN
            IF (idU3Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idU3Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU3Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV3Sd)') THEN
            IF (idV3Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idV3Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV3Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3Sd)') THEN
            IF (idW3Sd.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3Sd'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3St)') THEN
            IF (idW3St.eq.0) THEN
              IF (Master) WRITE (out,280) 'idW3St'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3St,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWamp)') THEN
            IF (idWamp.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWamp'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWamp,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWlen)') THEN
            IF (idWlen.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWlen'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWlen,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWlep)') THEN
            IF (idWlep.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWlep'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWlep,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdir)') THEN
            IF (idWdir.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWdir'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdir,1:Ngrids))
#ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWptp)') THEN
            IF (idWptp.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWptp'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWptp,1:Ngrids))
#endif
#ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWpbt)') THEN
            IF (idWpbt.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWpbt'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWpbt,1:Ngrids))
#endif
#ifdef WAVES_UB
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWorb)') THEN
            IF (idWorb.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWorb'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWorb,1:Ngrids))
#endif
#if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdif)') THEN
            IF (idWdif.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWdif'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdif,1:Ngrids))
#endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdib)') THEN
            IF (idWdib.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWdib'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdib,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdiw)') THEN
            IF (idWdiw.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWdiw'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdiw,1:Ngrids))
#endif
#ifdef ROLLER_SVENDSEN
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWbrk)') THEN
            IF (idWbrk.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWbrk'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWbrk,1:Ngrids))
#endif
#ifdef WEC_ROLLER
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdis)') THEN
            IF (idWdis.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWdis'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdis,1:Ngrids))
#endif
#ifdef ROLLER_RENIERS
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWrol)') THEN
            IF (idWrol.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWrol'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWrol,1:Ngrids))
#endif
#ifdef WEC_VF
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWztw)') THEN
            IF (idWztw.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWztw'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWztw,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWqsp)') THEN
            IF (idWqsp.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWqsp'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWqsp,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWbeh)') THEN
            IF (idWbeh.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWbeh'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWbeh,1:Ngrids))
#endif
#if defined UV_KIRBY
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUwav)') THEN
            IF (idUwav.eq.0) THEN
              WRITE (out,280) 'idUwav'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUwav,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVwav)') THEN
            IF (idVwav.eq.0) THEN
              WRITE (out,280) 'idVwav'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVwav,1:Ngrids))
#endif 
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
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
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLhea)') THEN
            IF (idLhea.eq.0) THEN
              IF (Master) WRITE (out,280) 'idLhea'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLhea,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idShea)') THEN
            IF (idShea.eq.0) THEN
              IF (Master) WRITE (out,280) 'idShea'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idShea,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLrad)') THEN
            IF (idLrad.eq.0) THEN
              IF (Master) WRITE (out,280) 'idLrad'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLrad,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSrad)') THEN
            IF (idSrad.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSrad'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idEmPf)') THEN
            IF (idEmPf.eq.0) THEN
              IF (Master) WRITE (out,280) 'idEmPf'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idevap)') THEN
            IF (idevap.eq.0) THEN
              IF (Master) WRITE (out,280) 'idevap'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idevap,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idrain)') THEN
            IF (idrain.eq.0) THEN
              IF (Master) WRITE (out,280) 'idrain'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idrain,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idDano)') THEN
            IF (idDano.eq.0) THEN
              IF (Master) WRITE (out,280) 'idDano'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idDano,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvis)') THEN
            IF (idVvis.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVvis'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvis,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTdif)') THEN
            IF (idTdif.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTdif'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTdif,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSdif)') THEN
            IF (idSdif.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSdif'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSdif,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHsbl)') THEN
            IF (idHsbl.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHsbl'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHsbl,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHbbl)') THEN
            IF (idHbbl.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHbbl'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHbbl,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtke)') THEN
            IF (idMtke.eq.0) THEN
              IF (Master) WRITE (out,280) 'idMtke'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtke,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtls)') THEN
            IF (idMtls.eq.0) THEN
              IF (Master) WRITE (out,280) 'idMtls'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtls,1:Ngrids))
#if defined SOLVE3D && defined T_PASSIVE
          ELSE IF (TRIM(KeyWord).eq.'Hout(inert)') THEN
            Npts=load_l(Nval, Cval, NPT*Ngrids, Linert)
            DO ng=1,Ngrids
              DO itrc=1,NPT
                i=idTvar(inert(itrc))
                Hout(i,ng)=Linert(itrc,ng)
              END DO
            END DO
#endif
#if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
          ELSE IF (TRIM(KeyWord).eq.'Hout(idBott)') THEN
            IF (MAXVAL(idBott).eq.0) THEN
              WRITE (out,280) 'idBott'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, MBOTP*Ngrids, Lbottom)
            DO ng=1,Ngrids
              DO itrc=1,MBOTP
                i=idBott(itrc)
                Hout(i,ng)=Lbottom(itrc,ng)
              END DO
            END DO
#endif
          ELSE IF (TRIM(KeyWord).eq.'NUSER') THEN
            Npts=load_i(Nval, Rval, 1, Nuser)
          ELSE IF (TRIM(KeyWord).eq.'USER') THEN
            Npts=load_r(Nval, Rval, MAX(1,Nuser), user)
          ELSE IF (TRIM(KeyWord).eq.'NC_SHUFFLE') THEN
            Npts=load_i(Nval, Rval, 1, shuffle)
          ELSE IF (TRIM(KeyWord).eq.'NC_DEFLATE') THEN
            Npts=load_i(Nval, Rval, 1, deflate)
          ELSE IF (TRIM(KeyWord).eq.'NC_DLEVEL') THEN
            Npts=load_i(Nval, Rval, 1, deflate_level)
          ELSE IF (TRIM(KeyWord).eq.'GSTNAME') THEN
            DO i=1,LEN(GSTname(Nval))
              GSTname(Nval)(i:i)=blank
            END DO
            GSTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'RSTNAME') THEN
            DO i=1,LEN(RSTname(Nval))
              RSTname(Nval)(i:i)=blank
            END DO
            RSTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'HISNAME') THEN
            DO i=1,LEN(HISname(Nval))
              HISname(Nval)(i:i)=blank
              HISbase(Nval)(i:i)=blank
            END DO
            HISname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            HISbase(Nval)=TRIM(ADJUSTL(HISname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'TLMNAME') THEN
            DO i=1,LEN(TLMname(Nval))
              TLMname(Nval)(i:i)=blank
              TLMbase(Nval)(i:i)=blank
            END DO
            TLMname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            TLMbase(Nval)=TRIM(ADJUSTL(TLMname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'TLFNAME') THEN
            DO i=1,LEN(TLMname(Nval))
              TLFname(Nval)(i:i)=blank
            END DO
            TLFname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ADJNAME') THEN
            DO i=1,LEN(ADJname(Nval))
              ADJname(Nval)(i:i)=blank
              ADJbase(Nval)(i:i)=blank
            END DO
            ADJname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            ADJbase(Nval)=TRIM(ADJUSTL(ADJname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'AVGNAME') THEN
            DO i=1,LEN(AVGname(Nval))
              AVGname(Nval)(i:i)=blank
              AVGbase(Nval)(i:i)=blank
            END DO
            AVGname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            AVGbase(Nval)=TRIM(ADJUSTL(AVGname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'DIANAME') THEN
            DO i=1,LEN(DIAname(Nval))
              DIAname(Nval)(i:i)=blank
              DIAbase(Nval)(i:i)=blank
            END DO
            DIAname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            DIAbase(Nval)=TRIM(ADJUSTL(DIAname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'STANAME') THEN
            DO i=1,LEN(STAname(Nval))
              STAname(Nval)(i:i)=blank
            END DO
            STAname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FLTNAME') THEN
            DO i=1,LEN(FLTname(Nval))
              FLTname(Nval)(i:i)=blank
            END DO
            FLTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'GRDNAME') THEN
            DO i=1,LEN(GRDname(Nval))
              GRDname(Nval)(i:i)=blank
            END DO
            GRDname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ININAME') THEN
            DO i=1,LEN(INIname(Nval))
              INIname(Nval)(i:i)=blank
            END DO
            INIname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'IRPNAME') THEN
            DO i=1,LEN(ITLname(Nval))
              IRPname(Nval)(i:i)=blank
            END DO
            IRPname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ITLNAME') THEN
            DO i=1,LEN(ITLname(Nval))
              ITLname(Nval)(i:i)=blank
            END DO
            ITLname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'IADNAME') THEN
            DO i=1,LEN(IADname(Nval))
              IADname(Nval)(i:i)=blank
            END DO
            IADname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NFFILES') THEN
            Npts=load_i(Nval, Rval, Ngrids, nFfiles)
            DO ng=1,Ngrids
              IF (nFfiles(ng).le.0) THEN
                IF (Master) WRITE (out,260) 'NFFILES',                  &
     &                            'Must be equal or greater than one.'
                exit_flag=4
                RETURN
              END IF
            END DO
            Npts=MAXVAL(nFfiles)
            allocate ( FRCids (Npts,Ngrids) )
            allocate ( FRCname(Npts,Ngrids) )
            FRCids(1:Npts,1:Ngrids)=-1
            DO ng=1,Ngrids
              DO k=1,Npts
                DO i=1,LEN(FRCname(k,ng))
                  FRCname(k,ng)(i:i)=blank
                END DO
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'FRCNAME') THEN
            DO ng=1,Ngrids
              DO i=1,nFfiles(ng)
                IF (Nval.eq.(i+(ng-1)*nFfiles(ng))) THEN
                  FRCname(i,ng)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'CLMNAME') THEN
            DO i=1,LEN(CLMname(Nval))
              CLMname(Nval)(i:i)=blank
            END DO
            CLMname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'BRYNAME') THEN
            DO i=1,LEN(BRYname(Nval))
              BRYname(Nval)(i:i)=blank
            END DO
            BRYname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
#ifdef INWAVE_COUPLING
          ELSE IF (TRIM(KeyWord).eq.'IWININAME') THEN
            DO i=1,LEN(IWINIname(Nval))
              IWINIname(Nval)(i:i)=blank
            END DO
            IWINIname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'IWBNDNAME') THEN
            DO i=1,LEN(IWBNDname(Nval))
              IWBNDname(Nval)(i:i)=blank
            END DO
            IWBNDname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
#endif
          ELSE IF (TRIM(KeyWord).eq.'FWDNAME') THEN
            DO i=1,LEN(FWDname(Nval))
              FWDname(Nval)(i:i)=blank
              FWDbase(Nval)(i:i)=blank
            END DO
            FWDname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            FWDbase(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ADSNAME') THEN
            DO i=1,LEN(ADSname(Nval))
              ADSname(Nval)(i:i)=blank
            END DO
            ADSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'APARNAM') THEN
            DO i=1,LEN(aparnam)
              aparnam(i:i)=blank
            END DO
            aparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SPOSNAM') THEN
            DO i=1,LEN(sposnam)
              sposnam(i:i)=blank
            END DO
            sposnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FPOSNAM') THEN
            DO i=1,LEN(fposnam)
              fposnam(i:i)=blank
            END DO
            fposnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'BPARNAM') THEN
            DO i=1,LEN(bparnam)
              bparnam(i:i)=blank
            END DO
            bparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SPARNAM') THEN
            DO i=1,LEN(sparnam)
              sparnam(i:i)=blank
            END DO
            sparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'USRNAME') THEN
            DO i=1,LEN(USRname)
              USRname(i:i)=blank
            END DO
            USRname=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
  10  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)

#if defined BULK_FLUXES && defined NL_BULK_FLUXES
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
!
!  Set switch to create NetCDF file.
!
      DO ng=1,Ngrids
#ifdef SOLVE3D
        DO i=1,NV
          IF (Hout(i,ng)) LdefHIS(ng)=.TRUE.
        END DO
#else
        IF (Hout(idFsur,ng).or.Hout(idUbar,ng).or.Hout(idVbar,ng)) THEN
          LdefHIS(ng)=.TRUE.
        END IF
#endif
#if defined IS4DVAR || defined WEAK_CONSTRAINT
!
!  If weak constraint, disallow recycling of the adjoint model.
!
        LcycleADJ(ng)=.FALSE.
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
#if defined FOUR_DVAR || defined IMPULSE
!
!  Set size of additonal dimension for error covariance normalization
!  and standard deviation factors.
!
# ifdef WEAK_CONSTRAINT
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
!  If weak constraint assimilation, ensure that nADJ=nHIS.
!
        IF (nADJ(ng).lt.ntimes(ng)) THEN
          nADJ(ng)=nHIS(ng)
        END IF
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
          WRITE (out,70) ng, Lm(ng), Mm(ng), N(ng), numnodes,           &
     &                   NtileI(ng), NtileJ(ng)
          IF ((NtileI(ng)*NtileJ(ng)).ne.numnodes) THEN
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
#ifdef REFINED_GRID
          WRITE (out,120) NestedGrids, 'Numgrids',                      &
     &          'Number of refined grids.'
#endif
#ifdef AIR_OCEAN
          WRITE (out,140) TI_ATM_OCN, 'TI_ATM_OCN',                     &
     &          'Time interval (s) between coupling ATM-OCN models.'
          WRITE (out,120) nOCN_ATM(ng), 'nOCN_ATM',                     &
     &          'Number of OCN timesteps between coupling to ATM.'
#endif
#ifdef WAVES_OCEAN
          WRITE (out,140) TI_WAV_OCN, 'TI_WAV_OCN',                     &
     &          'Time interval (s) between coupling WAV-OCN models.'
          WRITE (out,120) nOCN_WAV(ng), 'nOCN_WAV',                     &
     &          'Number of OCN timesteps between coupling to WAV.'
#endif
#ifdef WRF_COUPLING
          WRITE (out,120) WRF_CPL_GRID, 'WRF_CPLG_GRID',                &
     &          'WRF grid number that roms/swan will couple to.'
#endif
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
    defined W4DPSAS           || defined W4DVAR
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
#ifdef AVERAGES
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
          WRITE (out,130) nSFF(ng), 'nOBC',                             &
     &          'Number of timesteps between 4DVAR adjustment of',      &
     &          'open boundaries.'
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          WRITE (out,130) nSFF(ng), 'nSFF',                             &
     &          'Number of timesteps between 4DVAR adjustment of',      &
     &          'surface forcing fields.'
#endif
#ifdef PROPAGATOR
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
          DO itrc=1,NAT+NPT
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
          DO itrc=1,NAT+NPT
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
#ifdef SOLVE3D
          DO itrc=1,NAT+NPT
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
# endif
# ifdef WEC
          WRITE (out,140) wec_alpha(ng), 'wec_alpha',                   &
     &          'Factor for WEC roller/breaking dissipation.'
# endif
# ifdef FORWARD_MIXING
          DO itrc=1,NAT+NPT
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
            WRITE (out,260) 'Zob.',                                     &
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
            WRITE (out,260) 'blk_ZQ.',                                  &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZT(ng), 'blk_ZT',                         &
     &          'Height (m) of surface air temperature measurement.'
          IF (blk_ZT(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZT.',                                  &
     &            'It must be greater than zero.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,200) blk_ZW(ng), 'blk_ZW',                         &
     &          'Height (m) of surface winds measurement.'
          IF (blk_ZW(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZW.',                                  &
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
          IF ((lmd_Jwt(ng).lt.1).or.(lmd_Jwt(ng).gt.5)) THEN
            WRITE (out,260) 'lmd_Jwt.',                                 &
     &            'It must between one and five.'
            exit_flag=5
            RETURN
          END IF
# endif
# ifdef BODYFORCE
          WRITE (out,130) levsfrc(ng), 'levsfrc',                       &
     &          'Deepest level to apply surface stress as a',           &
     &          'bodyforce.'
          IF ((levsfrc(ng).lt.1).or.(levsfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levsfrc.',                                 &
     &            'Out of range surface bodyforce level.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,130) levbfrc(ng), 'levbfrc',                       &
     &          'Shallowest level to apply bottom stress as a',         &
     &          'bodyforce.'
          IF ((levbfrc(ng).lt.1).or.(levbfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levbfrc.',                                 &
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
          WRITE (out,160) Tcline(ng), 'Tcline',                         &
     &          'S-coordinate surface/bottom layer width (m) used',     &
     &          'in vertical coordinate stretching.'
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
          DO itrc=1,NAT+NPT
            WRITE (out,190) Tnudg(itrc,ng), 'Tnudg', itrc,              &
     &            'Nudging/relaxation time scale (days)',               &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef SCORRECTION
          IF (Tnudg(isalt,ng).le.0.0_r8) THEN
            WRITE (out,260) 'Tnudg(isalt).',                            &
     &            'Must be greater than zero for salt flux correction.'
            exit_flag=5
            RETURN
          END IF
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
            WRITE (out,260) 'KstrS',                                    &
     &            'Out of range ADM sensitivity starting level.'
            exit_flag=5
            RETURN
          END IF
          WRITE (out,120) KendS(ng), 'KendS',                           &
     &          'Shallowest level whose ADM sensitivity is required.'
          IF ((KendS(ng).lt.1).or.(KendS(ng).gt.N(ng))) THEN
            WRITE (out,260) 'KendS',                                    &
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
# ifdef SO_SEMI
#  ifndef SO_SEMI_WHITE
          WRITE (out,140) SO_decay(ng), 'SO_decay',                     &
     &          'Stochastic optimals time decorrelation scale (days).'
#  endif
          IF (SCALARS(ng)%SOstate(isUstr))                              &
     &       WRITE (out,170) SCALARS(ng)%SOstate(isUstr),               &
     &       'SOstate(isUstr)',                                         &
     &       'Stochastic optimals on surface U-stress.'
          IF (SCALARS(ng)%SOstate(isVstr))                              &
     &       WRITE (out,170) SCALARS(ng)%SOstate(isVstr),               &
     &       'SOstate(isVstr)',                                         &
     &       'Stochastic optimals on surface V-stress.'
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%SOstate(isTsur(itrc)))                      &
     &         WRITE (out,180) SCALARS(ng)%SOstate(isTsur(itrc)),       &
     &         'SOstate(idTsur)',                                       &
     &         'Stochastic optimals on surface flux of tracer',         &
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
#  endif
          IF (SCALARS(ng)%SOstate(isUstr))                              &
            WRITE (out,200) SO_sdev(isUstr,ng), 'SO_sdev(isUstr)',      &
     &          'Stochastic optimals scale, surface U-stress'
          IF (SCALARS(ng)%SOstate(isVstr))                              &
            WRITE (out,200) SO_sdev(isVstr,ng), 'SO_sdev(isVstr)',      &
     &          'Stochastic optimals scale, surface V-stress'
#  ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (SCALARS(ng)%SOstate(isTsur(itrc)))                      &
     &         WRITE (out,195) SO_sdev(isTsur(itrc),ng),                &
     &         'SO_sdev(idTsur)',                                       &
     &         'Stochastic optimals scale, surface flux of tracer',     &
     &         itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
#  endif
# endif
#endif
#if defined SOLVE3D && defined TS_PSOURCE
          DO itrc=1,NAT
            WRITE (out,185) LtracerSrc(itrc,ng), 'LtracerSrc', itrc,    &
     &            'Processing point sources/Sink on tracer ', itrc,     &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef T_PASSIVE
          DO itrc=1,NPT
            i=inert(itrc)
            WRITE (out,185) LtracerSrc(i,ng), 'LtracerSrc', i,          &
     &            'Processing point sources/Sink on tracer ', i,        &
     &            TRIM(Vname(1,idTvar(i)))
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
#ifdef SOLVE3D
          IF (Hout(idUvel,ng)) WRITE (out,170) Hout(idUvel,ng),         &
     &       'Hout(idUvel)',                                            &
     &       'Write out 3D U-momentum component.'
          IF (Hout(idVvel,ng)) WRITE (out,170) Hout(idVvel,ng),         &
     &       'Hout(idVvel)',                                            &
     &       'Write out 3D V-momentum component.'
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
#endif
          IF (Hout(idUair,ng)) WRITE (out,170) Hout(idUair,ng),         &
     &       'Hout(idUair)',                                            &
     &       'Write out surface U-wind.'
          IF (Hout(idVair,ng)) WRITE (out,170) Hout(idVair,ng),         &
     &       'Hout(idVair)',                                            &
     &       'Write out surface V-wind.'
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
          IF (Hout(idU2rs,ng)) WRITE (out,170) Hout(idU2rs,ng),         &
     &       'Hout(idU2rs)',                                            &
     &       'Write out total 2D u-radiation stress.'
          IF (Hout(idV2rs,ng)) WRITE (out,170) Hout(idV2rs,ng),         &
     &       'Hout(idV2rs)',                                            &
     &       'Write out total 2D v-radiation stress.'
# ifdef SOLVE3D
          IF (Hout(idW3xx,ng)) WRITE (out,170) Hout(idW3xx,ng),         &
     &       'Hout(idW3xx)',                                            &
     &       'Write out 3D horizonrtal radiation stress, Sxx.'
          IF (Hout(idW3xy,ng)) WRITE (out,170) Hout(idW3xy,ng),         &
     &       'Hout(idW3xy)',                                            &
     &       'Write out 3D horizonrtal radiation stress, Sxy.'
          IF (Hout(idW3yy,ng)) WRITE (out,170) Hout(idW3yy,ng),         &
     &       'Hout(idW3yy)',                                            &
     &       'Write out 3D horizonrtal radiation stress, Syy.'
          IF (Hout(idW3zx,ng)) WRITE (out,170) Hout(idW3zx,ng),         &
     &       'Hout(idW3zx)',                                            &
     &       'Write out 3D vertical radiation stress, Szx.'
          IF (Hout(idW3zy,ng)) WRITE (out,170) Hout(idW3zy,ng),         &
     &       'Hout(idW3zy)',                                            &
     &       'Write out 3D vertical radiation stress, Szy.'
          IF (Hout(idU3rs,ng)) WRITE (out,170) Hout(idU3rs,ng),         &
     &       'Hout(idU3rs)',                                            &
     &       'Write out total 3D u-radiation stress.'
          IF (Hout(idV3rs,ng)) WRITE (out,170) Hout(idV3rs,ng),         &
     &       'Hout(idV3rs)',                                            &
     &       'Write out total 3D v-radiation stress.'
# endif
#endif
#ifdef WEC
          IF (Hout(idU2Sd,ng)) WRITE (out,170) Hout(idU2Sd,ng),         &
     &       'Hout(idU2Sd)',                                            &
     &       'Write out 2D u-momentum stokes velocity.'
          IF (Hout(idV2Sd,ng)) WRITE (out,170) Hout(idV2Sd,ng),         &
     &       'Hout(idV2Sd)',                                            &
     &       'Write out 2D v-momentum stokes velocity.'
          IF (Hout(idU3Sd,ng)) WRITE (out,170) Hout(idU3Sd,ng),         &
     &       'Hout(idU3Sd)',                                            &
     &       'Write out 3D u-momentum stokes velocity.'
          IF (Hout(idV3Sd,ng)) WRITE (out,170) Hout(idV3Sd,ng),         &
     &       'Hout(idV3Sd)',                                            &
     &       'Write out 3D v-momentum stokes velocity.'
#endif
#ifdef WAVES_HEIGHT
          IF (Hout(idWamp,ng)) WRITE (out,170) Hout(idWamp,ng),         &
     &       'Hout(idWamp)',                                            &
     &       'Write out wave height.'
#endif
#ifdef WAVES_LENGTH
          IF (Hout(idWlen,ng)) WRITE (out,170) Hout(idWlen,ng),         &
     &       'Hout(idWlen)',                                            &
     &       'Write out mean wavelength.'
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
#ifdef WAVES_UB
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
#ifdef WEC_VF
          IF (Hout(idWztw,ng)) WRITE (out,170) Hout(idWztw,ng),         &
     &       'Hout(idWztw)',                                            &
     &       'Write out WEC quasi-static sea level adjustment.'
          IF (Hout(idWqsp,ng)) WRITE (out,170) Hout(idWqsp,ng),         &
     &       'Hout(idWqsp)',                                            &
     &       'Write out WEC quasi-static pressure.'
          IF (Hout(idWbeh,ng)) WRITE (out,170) Hout(idWbeh,ng),         &
     &       'Hout(idWbeh)',                                            &
     &       'Write out WEC Bernoulli head.'
#endif
#if defined SOLVE3D && defined T_PASSIVE
          DO itrc=1,NPT
            IF (Hout(idTvar(inert(itrc)),ng)) WRITE (out,180)           &
     &          Hout(idTvar(inert(itrc)),ng), 'Hout(inert)',            &
     &          'Write out inert passive tracer ', itrc,                &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
#endif
#if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
          DO itrc=1,MBOTP
            IF (Hout(idBott(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idBott(itrc),ng), 'Hout(idBott)',                  &
     &          'Write out bottom property ', itrc,                     &
     &          TRIM(Vname(1,idBott(itrc)))
          END DO
#endif
#ifdef SOLVE3D
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
     &       'Write out vertical viscosity coefficient.'
          IF (Hout(idTdif,ng)) WRITE (out,170) Hout(idTdif,ng),         &
     &       'Hout(idTdif)',                                            &
     &       'Write out vertical T-diffusion coefficient.'
          IF (Hout(idSdif,ng)) WRITE (out,170) Hout(idSdif,ng),         &
     &       'Hout(idSdif)',                                            &
     &       'Write out vertical S-diffusion coefficient.'
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
#if defined NETCDF4 && defined DEFLATE
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
#ifdef PROPAGATOR
          WRITE (out,230) '       Output GST Restart File:  ',          &
     &                    TRIM(GSTname(ng))
#endif
          WRITE (out,230) '           Output Restart File:  ',          &
     &                    TRIM(RSTname(ng))
          IF (LdefHIS(ng)) THEN
            IF (ndefHIS(ng).eq.0) THEN
              WRITE (out,230) '           Output History File:  ',      &
     &                        TRIM(HISname(ng))
            ELSE
              Lstr=LEN_TRIM(HISname(ng))
              WRITE (out,230) '      Prefix for History Files:  ',      &
     &                        HISname(ng)(1:Lstr-3)
            END IF
          END IF
#ifdef TANGENT
          IF (ndefTLM(ng).eq.0) THEN
            WRITE (out,230) '           Output Tangent File:  ',        &
     &                      TRIM(TLMname(ng))
          ELSE
            Lstr=LEN_TRIM(TLMname(ng))
            WRITE (out,230) '      Prefix for Tangent Files:  ',        &
     &                      TLMname(ng)(1:Lstr-3)
          END IF
#endif
#ifdef WEAK_CONSTRAINT
          WRITE (out,230) '   Output Impulse Forcing File:  ',          &
     &                    TRIM(TLFname(ng))
#endif
#ifdef ADJOINT
          IF (ndefADJ(ng).eq.0) THEN
            WRITE (out,230) '           Output Adjoint File:  ',        &
     &                      TRIM(ADJname(ng))
          ELSE
            Lstr=LEN_TRIM(ADJname(ng))
            WRITE (out,230) '      Prefix for Adjoint Files:  ',        &
     &                      ADJname(ng)(1:Lstr-3)
          END IF
#endif
#if defined FORWARD_WRITE && !defined FOUR_DVAR
          WRITE (out,230) '     Output Forward State File:  ',          &
     &                    TRIM(FWDname(ng))
#endif
#ifdef AVERAGES
          IF (ndefAVG(ng).eq.0) THEN
            WRITE (out,230) '          Output Averages File:  ',        &
     &                      TRIM(AVGname(ng))
          ELSE
            Lstr=LEN_TRIM(AVGname(ng))
            WRITE (out,230) '     Prefix for Averages Files:  ',        &
     &                      AVGname(ng)(1:Lstr-3)
          END IF
#endif
#ifdef DIAGNOSTICS
          IF (ndefDIA(ng).eq.0) THEN
            WRITE (out,230) '       Output Diagnostics File:  ',        &
     &                      TRIM(DIAname(ng))
          ELSE
            Lstr=LEN_TRIM(DIAname(ng))
            WRITE (out,230) '  Prefix for Diagnostics Files:  ',        &
     &                      DIAname(ng)(1:Lstr-3)
          END IF
#endif
#ifdef STATIONS
          WRITE (out,230) '          Output Stations File:  ',          &
     &                    TRIM(STAname(ng))
#endif
#ifdef FLOATS
          WRITE (out,230) '            Output Floats File:  ',          &
     &                    TRIM(FLTname(ng))
#endif
#ifdef MODEL_COUPLING
          WRITE (out,230) '      Physical parameters File:  ',          &
     &                    TRIM(Iname)
#endif
#ifndef ANA_GRID
          fname=GRDname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '               Input Grid File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef INI_FILE
# ifdef NONLINEAR
          fname=INIname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '  Input Nonlinear Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined TANGENT && \
    !(defined FOUR_DVAR         || defined IS4DVAR_SENSITIVITY || \
      defined OPT_OBSERVATIONS  || defined SANITY_CHECK        || \
      defined SENSITIVITY_4DVAR || defined TLM_CHECK)
          fname=ITLname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '    Input Tangent Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined WEAK_CONSTRAINT && \
    !(defined W4DPSAS        || defined W4DPSAS_SENSITIVITY)
          fname=IRPname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) 'Input Representer Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined ADJOINT && \
    !(defined AD_SENSITIVITY      || defined FOUR_DVAR         || \
      defined IS4DVAR_SENSITIVITY || defined OPT_OBSERVATIONS  || \
      defined SANITY_CHECK        || defined SENSITIVITY_4DVAR || \
      defined SO_SEMI)
          fname=IADname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '    Input Adjoint Initial File:  ',          &
     &                    TRIM(fname)
# endif
#endif
#ifdef FRC_FILE
          DO i=1,nFfiles(ng)
            fname=FRCname(i,ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,310) '         Input Forcing File ', i, ':  ',   &
     &                      TRIM(fname)
          END DO
#endif
#ifdef CLM_FILE
          fname=CLMname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '        Input Climatology File:  ',          &
     &                    TRIM(fname)
#endif
#if defined FORWARD_READ && !(defined FOUR_DVAR || defined PICARD_TEST)
          fname=FWDname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '      Input Forward State File:  ',          &
     &                    TRIM(fname)
#endif
#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
          fname=ADSname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) 'Input Adjoint Sensitivity File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef OBC_DATA
          fname=BRYname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '           Input Boundary File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef STATIONS
          fname=sposnam
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '        Station positions File:  ',          &
     &                    TRIM(fname)
#endif
#if defined ASSIMILATION || defined NUDGING || defined FOUR_DVAR
          fname=aparnam
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '  Assimilation Parameters File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef FLOATS
          fname=fposnam
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) ' Initial Floats Positions File:  ',          &
     &                    TRIM(fname)
#endif
#ifdef BIOLOGY
          fname=bparnam
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '       Biology Parameters File:  ',          &
     &                    TRIM(fname)
#endif
          fname=varname
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
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
          exit_flag=5
          IF (Master) THEN
            WRITE (out,320) 'Nouter = ', Nouter,                        &
     &        'Posterior analysis error available for Nouter=1 only.'
          END IF
        END IF
#endif

#ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Rescale active tracer parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO itrc=1,NAT+NPT
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
 120  FORMAT (1x,i10,2x,a,t30,a)
 130  FORMAT (1x,i10,2x,a,t30,a,/,t32,a)
 140  FORMAT (f11.3,2x,a,t30,a)
 150  FORMAT (f11.2,2x,a,t30,a)
 160  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
 170  FORMAT (10x,l1,2x,a,t30,a)
 180  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 185  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 195  FORMAT (1p,e11.4,2x,a,t30,a,i2.2,':',1x,a)
 200  FORMAT (1p,e11.4,2x,a,t30,a)
 210  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 220  FORMAT (/,' Output/Input Files:',/)
 230  FORMAT (2x,a,a)
 240  FORMAT (/,' Generic User Parameters:',/)
 250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t30,                         &
     &        'User parameter ',i2.2,'.')
 260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,/,15x,a)
 270  FORMAT (/,' READ_PHYPAR - could not find input file:  ',a)
 280  FORMAT (/,' READ_PHYPAR - variable info not yet loaded, ', a)
 290  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,i4,    &
     &        /,15x,a)
 300  FORMAT (/,' READ_PHYPAR - Invalid dimension parameter, ',a,'(',   &
     &        i2.2,')',/,15x,a)
 310  FORMAT (2x,a,i2.2,a,a)
 320  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,i4,/,15x,a)

      RETURN
      END SUBROUTINE read_PhyPar

#if defined ASSIMILATION || defined FOUR_DVAR || defined NUDGING || \
    defined VERIFICATION
      SUBROUTINE read_AssPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This subroutine reads in input model assimilation parameters.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# if defined FOUR_DVAR || defined VERIFICATION
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
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
      logical :: inhere

      integer :: Npts, Nval
      integer :: i, ib, igrid, itrc, k, ng, status

      integer :: decode_line, load_i, load_l, load_r

# if defined ASSIMILATION || defined NUDGING
      logical, dimension(MT,Ngrids) :: Lassi

      real(r8), dimension(MT,Ngrids) :: Rassi
# endif
# if defined FOUR_DVAR || defined VERIFICATION
      logical, dimension(MT) :: Ltracer
#  if defined ADJUST_STFLUX && defined SOLVE3D
      logical, dimension(MT,Ngrids) :: Ltsur
#  endif
#  ifdef ADJUST_BOUNDARY
      logical, dimension(4,Ngrids) :: Lbry
#   ifdef SOLVE3D
      logical, dimension(4,MT,Ngrids) :: Lbry_trc
      logical, dimension(MT,4) :: Lboundary
      real(r8), dimension(MT,4,Ngrids) :: Rboundary
#   endif
#  endif
      real(r8), dimension(MT,Ngrids) :: Rtracer
# endif
      real(r8), dimension(100) :: Rval

      character (len=1 ), parameter :: blank = ' '
# ifdef ADJUST_BOUNDARY
      character (len=7) :: Text
# endif
      character (len=40) :: KeyWord
      character (len=160) :: fname, line
      character (len=160), dimension(100) :: Cval

# if defined ASSIMILATION || defined NUDGING
!
!-----------------------------------------------------------------------
!  Read in assimilation parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lassimilate') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lassimilate)
          ELSE IF (TRIM(KeyWord).eq.'Emod0') THEN
            Npts=load_r(Nval, Rval, 1, Emod0)
          ELSE IF (TRIM(KeyWord).eq.'Tgrowth') THEN
            Npts=load_r(Nval, Rval, 1, Tgrowth)
          ELSE IF (TRIM(KeyWord).eq.'cor') THEN
            Npts=load_r(Nval, Rval, 1, cor)
          ELSE IF (TRIM(KeyWord).eq.'assi_SSH') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_SSH)
          ELSE IF (TRIM(KeyWord).eq.'assi_SST') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_SST)
          ELSE IF (TRIM(KeyWord).eq.'assi_UVsur') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_UVsur)
          ELSE IF (TRIM(KeyWord).eq.'assi_UV') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_UV)
          ELSE IF (TRIM(KeyWord).eq.'assi_T') THEN
            Npts=load_l(Nval, Cval, NAT*Ngrids, Lassi)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                assi_T(itrc,ng)=Lassi(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'M2nudass') THEN
            Npts=load_r(Nval, Rval, Ngrids, M2nudass)
          ELSE IF (TRIM(KeyWord).eq.'M3nudass') THEN
            Npts=load_r(Nval, Rval, Ngrids, M3nudass)
          ELSE IF (TRIM(KeyWord).eq.'Tnudass') THEN
            Npts=load_r(Nval, Rval, NAT*Ngrids, Rassi)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                Tnudass(itrc,ng)=Rassi(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'extend_SST') THEN
            Npts=load_l(Nval, Cval, Ngrids, extend_SST)
          ELSE IF (TRIM(KeyWord).eq.'ZmSST') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZmSST)
            DO ng=1,Ngrids
              ZmSST(ng)=-ABS(ZmSST(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZoSST') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZoSST)
            DO ng=1,Ngrids
              ZoSST(ng)=-ABS(ZoSST(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'npSST') THEN
            Npts=load_i(Nval, Rval, 1, npSST)
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_SST') THEN
            DO ng=1,Ngrids
              DO i=0,npSST(ng)
                READ (inp,*) igrid, k, pcoef_SST(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Perr_SST') THEN
            DO ng=1,Ngrids
              DO i=0,npSST(ng)
                READ (inp,*) igrid, k, perr_SST(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'extend_UV') THEN
            Npts=load_l(Nval, Cval, Ngrids, extend_UV)
          ELSE IF (TRIM(KeyWord).eq.'ZmUV') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZmUV)
            DO ng=1,Ngrids
              ZmUV(ng)=-ABS(ZmUV(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZoUV') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZoUV)
            DO ng=1,Ngrids
              ZoUV(ng)=-ABS(ZoUV(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'npUV') THEN
            Npts=load_i(Nval, Rval, 1, npUV)
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_U') THEN
            DO ng=1,Ngrids
              DO i=0,npUV(ng)
                READ (inp,*) igrid, k, pcoef_U(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_V') THEN
            DO ng=1,Ngrids
              DO i=0,npUV(ng)
                READ (inp,*) igrid, k, pcoef_V(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Perr_UV') THEN
            DO ng=1,Ngrids
              DO i=0,npUV(ng)
                READ (inp,*) igrid, k, perr_V(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SSHname') THEN
            SSHname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SSTname') THEN
            SSTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'VSURname') THEN
            VSURname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'VOBSname') THEN
            VOBSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'TOBSname') THEN
            TOBSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
 10   IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
 20   CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lassimilate(ng)) THEN
            WRITE (out,60) ng
#  ifdef ASSIMILATION
            WRITE (out,100) Emod0(ng), 'Emod0',                         &
     &            'Initial model error (percentage).'
            WRITE (out,100) Tgrowth(ng), 'Tgrowth',                     &
     &            'Empirical model error growth scale (days).'
            WRITE (out,100) cor(ng), 'cor',                             &
     &            'Correlation between model and observations.'
            WRITE (out,80) assi_SSH(ng), 'assi_SSH',                    &
     &            'Activate assimilation of SSH.'
#   ifdef SOLVE3D
            WRITE (out,80) assi_SST(ng), 'assi_SST',                    &
     &            'Activate assimilation of SST.'
            DO itrc=1,NAT
              WRITE (out,90) assi_T(itrc,ng), 'assi_T', itrc,           &
     &              'Activate assimilation of tracer', itrc,            &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            WRITE (out,80) assi_UVsur(ng), 'assi_UVsur',                &
     &            'Activate assimilation of surface currents.'
            WRITE (out,80) assi_UV(ng), 'assi_UV',                      &
     &            'Activate assimilation of currents.'
#   endif
#  endif
#  ifdef SOLVE3D
#   if defined NUDGING_T || defined NUDGING_SST
            DO itrc=1,NAT
              IF (assi_T(itrc,ng)) THEN
                WRITE (out,120) Tnudass(itrc,ng), 'Tnudass', itrc,      &
     &                'Nudging assimilation time scale (days)',         &
     &                'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
              END IF
            END DO
#   endif
#  endif
#  ifdef SOLVE3D
#   if defined NUDGING_UV || defined NUDGING_UVsur
            WRITE (out,135) M3nudass(ng), 'M3nudass',                   &
     &            'Nudging assimilation time scale (days)',             &
     &            'for 3D momentum.'
#   endif
#   if defined NUDGING_SST || defined ASSIMILATION_SST
            IF (extend_SST(ng)) THEN
              WRITE (out,80) extend_SST(ng), 'extend_SST',              &
     &              'Extend SST vertically.'
              WRITE (out,70) npSST(ng), 'npSST',                        &
     &              'Order of polynomial for SST extension.'
              WRITE (out,100) ZmSST(ng), 'zmSST',                       &
     &              'SST, maximum extension depth (m).'
              WRITE (out,100) ZoSST(ng), 'zoSST',                       &
     &              'E-folding depth (m) to extend SST error variance.'
              IF (KP.lt.npSST(ng)) THEN
                WRITE (out,160) 'KP: ', KP, npSST(ng)
                exit_flag=5
                RETURN
              END IF
            END IF
            IF (npSST(ng).gt.0) THEN
              WRITE (out,*)
              DO k=0,npSST(ng)
                WRITE (out,140) k, pcoef_SST(k,ng), perr_SST(k,ng)
              END DO
            END IF
#   endif
#   if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
            IF (extend_UV(ng)) THEN
              WRITE (out,80) extend_UV(ng), 'extend_UV',                &
     &              'Extend surface currents vertically.'
              WRITE (out,70) npUV(ng), 'npUV',                          &
     &             'Order of polynomial for surface currents extension.'
              WRITE (out,100) ZmUV(ng), 'zmUV',                         &
     &              'Surface currents, maximum extension depth (m).'
              WRITE (out,110) ZoSST(ng), 'zoSST',                       &
     &              'E-folding depth (m) to extend surface currents',   &
     &              'error variance.'
              IF (KP.lt.npUV(ng)) THEN
                WRITE (out,160) 'KP: ', KP, npUV(ng)
                exit_flag=5
                RETURN
              END IF
            END IF
            IF (npUV(ng).gt.0) THEN
              WRITE (out,*)
              DO k=0,npUV(ng)
                WRITE (out,150) k, pcoef_U(k,ng), pcoef_V(k,ng),        &
     &                             perr_V(k,ng)
              END DO
            END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Report input files and check availability of input files.
!-----------------------------------------------------------------------
!
            WRITE (out,170)
            WRITE (out,180) ' Assimilation parameters File:  ',         &
     &                      TRIM(aparnam)
#  ifdef ASSIMILATION_SSH
            fname=SSHname(ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,180) '      Sea Surface Height File:  ',         &
     &                      TRIM(fname)
#  endif
#  ifdef SOLVE3D
#   if defined ASSIMILATION_SST || defined NUDGING_SST
            fname=SSTname(ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,180) ' Sea Surface Temperature File:  ',         &
     &                      TRIM(fname)
#   endif
#   if defined ASSIMILATION_T || defined NUDGING_T
            fname=TOBSname(ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,180) '                 Tracers File:  ',         &
     &                      TRIM(fname)
#   endif
#   if defined ASSIMILATION_UVsur || defined NUDGING_UVsur
            fname=VSURname(ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,180) '        Surface Currents File:  ',         &
     &                      TRIM(fname)
#   endif
#   if defined ASSIMILATION_UV || defined NUDGING_UV
            fname=VOBSname(ng)
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
            WRITE (out,180) '     Horizontal Currents File:  ',         &
     &                      TRIM(fname)
#   endif
#  endif
            GO TO 40
  30        WRITE (out,190) TRIM(fname)
            exit_flag=4
            RETURN
  40        CONTINUE
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant parameters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        Emod0(ng)=1.0_r8-0.01_r8*Emod0(ng)
# ifdef NUDGING
        IF (M2nudass(ng).gt.0.0_r8) THEN
          M2nudass(ng)=1.0_r8/(M2nudass(ng)*86400.0_r8)
        ELSE
          M2nudass(ng)=0.0_r8
        END IF
#  ifdef SOLVE3D
        IF (M3nudass(ng).gt.0.0_r8) THEN
          M3nudass(ng)=1.0_r8/(M3nudass(ng)*86400.0_r8)
        ELSE
          M3nudass(ng)=0.0_r8
        END IF
        DO itrc=1,NAT
          IF (Tnudass(itrc,ng).gt.0.0_r8) THEN
            Tnudass(itrc,ng)=1.0_r8/(Tnudass(itrc,ng)*86400.0_r8)
          ELSE
            Tnudass(itrc,ng)=0.0_r8
          END IF
        END DO
#  endif
# endif
      END DO

  50  FORMAT (/,' READ_AssPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Assimilation Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  70  FORMAT (1x,i10,2x,a,t30,a)
  80  FORMAT (10x,l1,2x,a,t30,a)
  90  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,1x,i2.2,':',1x,a)
 100  FORMAT (f11.3,2x,a,t30,a)
 110  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
 120  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 130  FORMAT (1p,e11.4,2x,a,t30,a)
 135  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 140  FORMAT (1x,'SST,Err polynomial order ',                           &
     &        i2.2,1x,1p,e15.8,1x,1p,e15.8)
 150  FORMAT (1x,'U,V,Err polynomial order: ',                          &
     &        i2.2,1x,1p,e15.8,1x,1p,e15.8,1x,1p,e15.8)
 160  FORMAT (/,' READ_ASSPAR - too small dimension parameter, ',a,     &
     &        2i4,/,15x,'change file  mod_scalars.F  and recompile.')
 170  FORMAT (/,' Input Assimilation Files:',/)
 180  FORMAT (2x,a,a)
 190  FORMAT (/,' READ_ASSPAR - could not find input file:  ',a)

# elif defined FOUR_DVAR || defined VERIFICATION
!
!-----------------------------------------------------------------------
!  Read in 4DVAR assimilation parameters. Then, load input data into
!  module.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'dTdz_min') THEN
            Npts=load_r(Nval, Rval, Ngrids, dTdz_min)
          ELSE IF (TRIM(KeyWord).eq.'ml_depth') THEN
            Npts=load_r(Nval, Rval, Ngrids, ml_depth)
            DO ng=1,Ngrids
              ml_depth(ng)=ABS(ml_depth(ng))
            END DO
#  if defined BALANCE_OPERATOR && defined ZETA_ELLIPTIC
          ELSE IF (TRIM(KeyWord).eq.'Nbico') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nbico)
#  endif
          ELSE IF (TRIM(KeyWord).eq.'LNM_depth') THEN
            Npts=load_r(Nval, Rval, Ngrids, LNM_depth)
            DO ng=1,Ngrids
              LNM_depth(ng)=ABS(LNM_depth(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'LNM_flag') THEN
            Npts=load_i(Nval, Rval, 1, LNM_flag)
          ELSE IF (TRIM(KeyWord).eq.'GradErr') THEN
            Npts=load_r(Nval, Rval, 1, GradErr)
          ELSE IF (TRIM(KeyWord).eq.'HevecErr') THEN
            Npts=load_r(Nval, Rval, 1, HevecErr)
          ELSE IF (TRIM(KeyWord).eq.'LhessianEV') THEN
            Npts=load_l(Nval, Cval, 1, LhessianEV)
          ELSE IF (TRIM(KeyWord).eq.'LhotStart') THEN
            Npts=load_l(Nval, Cval, 1, LhotStart)
          ELSE IF (TRIM(KeyWord).eq.'Lprecond') THEN
            Npts=load_l(Nval, Cval, 1, Lprecond)
#  if defined WEAK_CONSTRAINT
            IF ( LhessianEV.and.Lprecond ) THEN
              LhessianEV=.FALSE.
            END IF
#  endif
          ELSE IF (TRIM(KeyWord).eq.'Lritz') THEN
            Npts=load_l(Nval, Cval, 1, Lritz)
          ELSE IF (TRIM(KeyWord).eq.'NritzEV') THEN
            Npts=load_i(Nval, Rval, 1, NritzEV)
          ELSE IF (TRIM(KeyWord).eq.'NpostI') THEN
            Npts=load_i(Nval, Rval, 1, NpostI)
          ELSE IF (TRIM(KeyWord).eq.'tl_M2diff') THEN
            Npts=load_r(Nval, Rval, Ngrids, tl_M2diff)
          ELSE IF (TRIM(KeyWord).eq.'tl_M3diff') THEN
            Npts=load_r(Nval, Rval, Ngrids, tl_M3diff)
          ELSE IF (TRIM(KeyWord).eq.'tl_Tdiff') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                tl_Tdiff(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'LdefNRM') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, LdefNRM)
          ELSE IF (TRIM(KeyWord).eq.'LwrtNRM') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, LwrtNRM)
          ELSE IF (TRIM(KeyWord).eq.'CnormM(isFsur)') THEN
            IF (isFsur.eq.0) THEN
              IF (Master) WRITE (out,190) 'isFsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(2,isFsur))
          ELSE IF (TRIM(KeyWord).eq.'CnormM(isUbar)') THEN
            IF (isUbar.eq.0) THEN
              IF (Master) WRITE (out,190) 'isUbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(2,isUbar))
          ELSE IF (TRIM(KeyWord).eq.'CnormM(isVbar)') THEN
            IF (isVbar.eq.0) THEN
              IF (Master) WRITE (out,190) 'isVbar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(2,isVbar))
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'CnormM(isUvel)') THEN
            IF (isUvel.eq.0) THEN
              IF (Master) WRITE (out,190) 'isUvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(2,isUvel))
          ELSE IF (TRIM(KeyWord).eq.'CnormM(isVvel)') THEN
            IF (isVvel.eq.0) THEN
              IF (Master) WRITE (out,190) 'isVvel'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(2,isVvel))
          ELSE IF (TRIM(KeyWord).eq.'CnormM(isTvar)') THEN
            IF (MAXVAL(isTvar).eq.0) THEN
              IF (Master) WRITE (out,190) 'isTvar'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, MT, Ltracer)
            DO itrc=1,MT
              i=isTvar(itrc)
              Cnorm(2,i)=Ltracer(itrc)
            END DO
#  endif
          ELSE IF (TRIM(KeyWord).eq.'CnormI(isFsur)') THEN
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isFsur))
          ELSE IF (TRIM(KeyWord).eq.'CnormI(isUbar)') THEN
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isUbar))
          ELSE IF (TRIM(KeyWord).eq.'CnormI(isVbar)') THEN
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isVbar))
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'CnormI(isUvel)') THEN
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isUvel))
          ELSE IF (TRIM(KeyWord).eq.'CnormI(isVvel)') THEN
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isVvel))
          ELSE IF (TRIM(KeyWord).eq.'CnormI(isTvar)') THEN
            Npts=load_l(Nval, Cval, MT, Ltracer)
            DO itrc=1,MT
              i=isTvar(itrc)
              Cnorm(1,i)=Ltracer(itrc)
            END DO
#  endif
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(KeyWord).eq.'CnormB(isFsur)') THEN
            Npts=load_l(Nval, Cval, 4, CnormB(isFsur,:))
          ELSE IF (TRIM(KeyWord).eq.'CnormB(isUbar)') THEN
            Npts=load_l(Nval, Cval, 4, CnormB(isUbar,:))
          ELSE IF (TRIM(KeyWord).eq.'CnormB(isVbar)') THEN
            Npts=load_l(Nval, Cval, 4, CnormB(isVbar,:))
#   ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'CnormB(isUvel)') THEN
            Npts=load_l(Nval, Cval, 4, CnormB(isUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'CnormB(isVvel)') THEN
            Npts=load_l(Nval, Cval, 4, CnormB(isVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'CnormB(isTvar)') THEN
            Npts=load_l(Nval, Cval, MT*4, Lboundary)
            DO ib=1,4
              DO itrc=1,MT
                i=isTvar(itrc)
                CnormB(i,ib)=Lboundary(itrc,ib)
              END DO
            END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(KeyWord).eq.'CnormF(isUstr)') THEN
            IF (isUstr.eq.0) THEN
              IF (Master) WRITE (out,190) 'isUstr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isUstr))
          ELSE IF (TRIM(KeyWord).eq.'CnormF(isVstr)') THEN
            IF (isVstr.eq.0) THEN
              IF (Master) WRITE (out,190) 'isVstr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(1,isVstr))
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'CnormF(isTsur)') THEN
            IF (MAXVAL(isTsur).eq.0) THEN
              IF (Master) WRITE (out,190) 'isTsur'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, MT, Ltracer)
            DO itrc=1,MT
              i=isTsur(itrc)
              Cnorm(1,i)=Ltracer(itrc)
            END DO
#  endif
          ELSE IF (TRIM(KeyWord).eq.'balance(isSalt)') THEN
            Npts=load_l(Nval, Cval, 1, balance(isTvar(isalt)))
          ELSE IF (TRIM(KeyWord).eq.'balance(isFsur)') THEN
            Npts=load_l(Nval, Cval, 1, balance(isFsur))
          ELSE IF (TRIM(KeyWord).eq.'balance(isVbar)') THEN
            Npts=load_l(Nval, Cval, 1, balance(isVbar))
          ELSE IF (TRIM(KeyWord).eq.'balance(isVvel)') THEN
            Npts=load_l(Nval, Cval, 1, balance(isVvel))
          ELSE IF (TRIM(KeyWord).eq.'Nmethod') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nmethod)
          ELSE IF (TRIM(KeyWord).eq.'Rscheme') THEN
            Npts=load_i(Nval, Rval, Ngrids, Rscheme)
          ELSE IF (TRIM(KeyWord).eq.'Nrandom') THEN
            Npts=load_i(Nval, Rval, 1, Nrandom)
          ELSE IF (TRIM(KeyWord).eq.'Hgamma') THEN
            Npts=load_r(Nval, Rval, 4, Hgamma)
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Vgamma') THEN
            Npts=load_r(Nval, Rval, 4, Vgamma)
#  endif
          ELSE IF (TRIM(KeyWord).eq.'HdecayM(isFsur)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isFsur,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayM(isUbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isUbar,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayM(isVbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isVbar,:))
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'HdecayM(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayM(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(2,isVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayM(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Hdecay(2,isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VdecayM(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vdecay(2,isUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'VdecayM(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vdecay(2,isVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'VdecayM(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Vdecay(2,isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#  endif
          ELSE IF (TRIM(KeyWord).eq.'HdecayI(isFsur)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isFsur,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayI(isUbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUbar,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayI(isVbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVbar,:))
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'HdecayI(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayI(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayI(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Hdecay(1,isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VdecayI(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vdecay(1,isUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'VdecayI(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vdecay(1,isVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'VdecayI(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Vdecay(1,isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#  endif
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(KeyWord).eq.'HdecayB(isFsur)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isFsur,:,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayB(isUbar)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isUbar,:,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayB(isVbar)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isVbar,:,:))
#   ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'HdecayB(isUvel)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isUvel,:,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayB(isVvel)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, HdecayB(isVvel,:,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayB(isTvar)') THEN
            Npts=load_r(Nval, Rval, 4*MT*Ngrids, Rboundary)
            DO ng=1,Ngrids
              DO ib=1,4
                DO itrc=1,NT(ng)
                  HdecayB(isTvar(itrc),ib,ng)=Rboundary(itrc,ib,ng)
                END DO
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VdecayB(isUvel)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, VdecayB(isUvel,:,:))
          ELSE IF (TRIM(KeyWord).eq.'VdecayB(isVvel)') THEN
            Npts=load_r(Nval, Rval, 4*Ngrids, VdecayB(isVvel,:,:))
          ELSE IF (TRIM(KeyWord).eq.'VdecayB(isTvar)') THEN
            Npts=load_r(Nval, Rval, 4*MT*Ngrids, Rboundary)
            DO ng=1,Ngrids
              DO ib=1,4
                DO itrc=1,NT(ng)
                  VdecayB(isTvar(itrc),ib,ng)=Rboundary(itrc,ib,ng)
                END DO
              END DO
            END DO
#   endif
#  endif
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(KeyWord).eq.'HdecayF(isUstr)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isUstr,:))
          ELSE IF (TRIM(KeyWord).eq.'HdecayF(isVstr)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(1,isVstr,:))
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'HdecayF(isTsur)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Hdecay(1,isTsur(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#  endif
#  if defined ADJUST_STFLUX && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Lstflux') THEN
            Npts=load_l(Nval, Cval, MT*Ngrids, Ltsur)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Lstflux(itrc,ng)=Ltsur(itrc,ng)
              END DO
            END DO
#  endif
#  ifdef ADJUST_BOUNDARY
          ELSE IF (TRIM(KeyWord).eq.'Lobc(isFsur)') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
            DO ng=1,Ngrids
              DO ib=1,4
                Lobc(ib,isFsur,ng)=Lbry(ib,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lobc(isUbar)') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
            DO ng=1,Ngrids
              DO ib=1,4
                Lobc(ib,isUbar,ng)=Lbry(ib,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lobc(isVbar)') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
            DO ng=1,Ngrids
              DO ib=1,4
                Lobc(ib,isVbar,ng)=Lbry(ib,ng)
              END DO
            END DO
#   ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Lobc(isUvel)') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
            DO ng=1,Ngrids
              DO ib=1,4
                Lobc(ib,isUvel,ng)=Lbry(ib,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lobc(isVvel)') THEN
            Npts=load_l(Nval, Cval, 4*Ngrids, Lbry)
            DO ng=1,Ngrids
              DO ib=1,4
                Lobc(ib,isVvel,ng)=Lbry(ib,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lobc(isTvar)') THEN
            Npts=load_l(Nval, Cval, 4*MT*Ngrids, Lbry_trc)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                i=isTvar(itrc)
                DO ib=1,4
                  Lobc(ib,i,ng)=Lbry_trc(ib,itrc,ng)
                END DO
              END DO
            END DO
#   endif
#  endif
          ELSE IF (TRIM(KeyWord).eq.'STDnameI') THEN
            DO i=1,LEN(STDname(1,Nval))
              STDname(1,Nval)(i:i)=blank
            END DO
            STDname(1,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'STDnameM') THEN
            DO i=1,LEN(STDname(2,Nval))
              STDname(2,Nval)(i:i)=blank
            END DO
            STDname(2,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'STDnameB') THEN
            DO i=1,LEN(STDname(3,Nval))
              STDname(3,Nval)(i:i)=blank
            END DO
            STDname(3,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'STDnameF') THEN
            DO i=1,LEN(STDname(4,Nval))
              STDname(4,Nval)(i:i)=blank
            END DO
            STDname(4,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NRMnameI') THEN
            DO i=1,LEN(NRMname(1,Nval))
              NRMname(1,Nval)(i:i)=blank
            END DO
            NRMname(1,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NRMnameM') THEN
            DO i=1,LEN(NRMname(2,Nval))
              NRMname(2,Nval)(i:i)=blank
            END DO
            NRMname(2,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NRMnameB') THEN
            DO i=1,LEN(NRMname(3,Nval))
              NRMname(3,Nval)(i:i)=blank
            END DO
            NRMname(3,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NRMnameF') THEN
            DO i=1,LEN(NRMname(4,Nval))
              NRMname(4,Nval)(i:i)=blank
            END DO
            NRMname(4,Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'OBSname') THEN
            DO i=1,LEN(OBSname(Nval))
              OBSname(Nval)(i:i)=blank
            END DO
            OBSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'HSSname') THEN
            DO i=1,LEN(HSSname(Nval))
              HSSname(Nval)(i:i)=blank
              HSSbase(Nval)(i:i)=blank
            END DO
            HSSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            HSSbase(Nval)=TRIM(ADJUSTL(HSSname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'LCZname') THEN
            DO i=1,LEN(LCZname(Nval))
              LCZname(Nval)(i:i)=blank
            END DO
            LCZname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'MODname') THEN
            DO i=1,LEN(MODname(Nval))
              MODname(Nval)(i:i)=blank
            END DO
            MODname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ERRname') THEN
            DO i=1,LEN(ERRname(Nval))
              ERRname(Nval)(i:i)=blank
            END DO
            ERRname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
 10   IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
 20   CONTINUE

#  ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  Check switches to adjust boundaries for consistency.
!-----------------------------------------------------------------------
!
!  Make sure that both momentum components are activated for processing.
!  If adjusting 2D momentum in 3D applications, make sure that the
!  free-surface and 3D momentum switches are activated. This is because
!  the 2D momentum adjustments are computed from the vertical integral
!  of the 3D momentum increments.
!
      DO ng=1,Ngrids
        DO ib=1,4
          IF (.not.Lobc(ib,isUbar,ng).and.Lobc(ib,isVbar,ng)) THEN
            Lobc(ib,isUbar,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isVbar,ng).and.Lobc(ib,isUbar,ng)) THEN
            Lobc(ib,isVbar,ng)=.TRUE.
          END IF
#   ifdef SOLVE3D
          IF (.not.Lobc(ib,isUvel,ng).and.Lobc(ib,isVvel,ng)) THEN
            Lobc(ib,isUvel,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isVvel,ng).and.Lobc(ib,isUvel,ng)) THEN
            Lobc(ib,isVvel,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isFsur,ng).and.Lobc(ib,isUbar,ng)) THEN
            Lobc(ib,isFsur,ng)=.TRUE.
          END IF
          IF (.not.Lobc(ib,isUvel,ng).and.Lobc(ib,isUbar,ng)) THEN
            Lobc(ib,isUvel,ng)=.TRUE.
            Lobc(ib,isVvel,ng)=.TRUE.
          END IF
#   endif
        END DO
      END DO
#  endif

#  ifdef BALANCE_OPERATOR
!
!-----------------------------------------------------------------------
!  Check balance operator switches for consitency.
!-----------------------------------------------------------------------
!
#   ifdef SOLVE3D
      IF (.not.balance(isTvar(isalt)).and.balance(isFsur)) THEN
        balance(isTvar(isalt))=.TRUE.
      END IF
      IF (.not.balance(isTvar(isalt)).and.balance(isVvel)) THEN
        balance(isTvar(isalt))=.TRUE.
      END IF
      IF (balance(isTvar(isalt))) THEN
        balance(isTvar(itemp))=.TRUE.
      END IF
      IF (balance(isVvel)) THEN
        balance(isUvel)=.TRUE.
      END IF
      IF (balance(isVbar)) THEN
        balance(isVbar)=.FALSE.
      END IF
#   else
      IF (balance(isVbar)) THEN
        balance(isUbar)=.TRUE.
      END IF
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
#  if defined FOUR_DVAR
          WRITE (out,60) ng
#   ifdef WEAK_CONSTRAINT
#    ifdef RPM_RELAXATION
          WRITE (out,120) tl_M2diff, 'tl_M2diff',                       &
     &            'RPM 2D momentum diffusive relaxation coefficient.'
#     ifdef SOLVE3D
          WRITE (out,130) tl_M3diff, 'tl_M3diff',                       &
     &            'RPM 3D momentum diffusive relaxation coefficient.'
          DO itrc=1,NT(ng)
            WRITE (out,130) tl_Tdiff(itrc,ng), 'tl_Tdiff',              &
     &            'RPM tracer diffusive relaxation coefficient, ',      &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#   endif
#   ifndef IS4DVAR_SENSITIVITY
#    ifdef IS4DVAR
          WRITE (out,100) GradErr, 'GradErr',                           &
     &            'Upper bound on relative error of the gradient.'
          WRITE (out,100) HevecErr, 'HevecErr',                         &
     &            'Accuracy required for Hessian eigenvectors.'
          WRITE (out,70) LhessianEV, 'LhessianEV',                      &
     &            'Switch to compute Hessian eigenvectors.'
#    endif
#    ifdef WEAK_CONSTRAINT
          WRITE (out,70) LhotStart, 'LhotStart',                        &
     &            'Switch for hot start of subsequent outer loops.'
#    endif
          WRITE (out,70) Lprecond, 'Lprecond',                          &
     &            'Switch for conjugate gradient preconditioning.'
          WRITE (out,70) Lritz, 'Lritz',                                &
     &            'Switch for Ritz limited-memory preconditioning.'
#    ifdef WEAK_CONSTRAINT
          IF (Lprecond.and.(NritzEV.gt.0)) THEN
            WRITE (out,80) NritzEV, 'NritzEV',                          &
     &            'Number of preconditioning eigenpairs to use.'
          END IF
#    endif
#   endif
#   ifdef BALANCE_OPERATOR
#    ifdef ZETA_ELLIPTIC
          WRITE (out,80) Nbico(ng), 'Nbico',                            &
     &            'Number of iterations in SSH elliptic equation.'
#    endif
          WRITE (out,100) dTdz_min(ng), 'dTdz_min',                     &
     &            'Minimum dTdz (C/m) used in balanced salinity.'
          WRITE (out,100) ml_depth(ng), 'ml_depth',                     &
     &            'Mixed-layer depth (m) used in balanced salinity.'
          IF (balance(isFsur)) THEN
            WRITE (out,100) LNM_depth(ng), 'LNM_depth',                 &
     &            'Level of no motion (m) in balanced free-sruface.'
            WRITE (out,80) LNM_flag, 'LNM_flag',                        &
     &            'Level of no motion integration flag.'
          END IF
          WRITE (out,70) balance(isFsur), 'balance(isFsur)',            &
                  'Switch to include free-surface in balance operator.'
#    ifdef SOLVE3D
          WRITE (out,70) balance(isTvar(isalt)), 'balance(isSalt)',     &
                  'Switch to include salinity in balance operator.'
          WRITE (out,70) balance(isVvel), 'balance(isVvel)',            &
                  'Switch to include 3D momentum in balance operator.'
#    else
          WRITE (out,70) balance(isVbar), 'balance(isVbar)',            &
                  'Switch to include 2D momentum in balance operator.'
#    endif
#   endif
#   if defined POSTERIOR_EOFS && defined WEAK_CONSTRAINT
          WRITE (out,80) NpostI, 'NpostI',                              &
     &            'Number of Lanczos iterations in posterior analysis.'
#   endif
#   ifndef TLM_CHECK
#    ifndef IS4DVAR_SENSITIVITY
          WRITE (out,170) LdefNRM(1:4,ng), 'LdefNRM',                   &
     &            'Switch to create a normalization NetCDF file.'
          WRITE (out,170) LwrtNRM(1:4,ng), 'LwrtNRM',                   &
     &            'Switch to write out normalization factors.'
          IF (ANY(LwrtNRM(:,ng))) THEN
            IF (Nmethod(ng).eq.0) THEN
              WRITE (out,80) Nmethod(ng), 'Nmethod',                    &
     &            'Correlation normalization method: Exact.'
            ELSE IF (Nmethod(ng).eq.1) THEN
              WRITE (out,80) Nmethod(ng), 'Nmethod',                    &
     &            'Correlation normalization method: Randomization.'
              WRITE (out,80) Rscheme(ng), 'Rscheme',                    &
     &            'Random number generation scheme'
              WRITE (out,80) Nrandom, 'Nrandom',                        &
     &            'Number of iterations for randomization.'
            END IF
          END IF
#     if defined SENSITIVITY_4DVAR || \
         defined TL_W4DPSAS        || defined TL_W4DVAR || \
         defined W4DPSAS           || defined W4DVAR
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(2,isFsur), 'CnormM(isFsur)',           &
     &            'Compute model 2D RHO-normalization factors.'
            WRITE (out,70) Cnorm(2,isUbar), 'CnormM(isUbar)',           &
     &            'Compute model 2D U-normalization factors.'
            WRITE (out,70) Cnorm(2,isVbar), 'CnormM(isVbar)',           &
     &            'Compute model 2D V-normalization factors.'
#      ifdef SOLVE3D
            WRITE (out,70) Cnorm(2,isUvel), 'CnormM(isUvel)',           &
     &            'Compute model 3D U-normalization factors.'
            WRITE (out,70) Cnorm(2,isVvel), 'CnormM(isVvel)',           &
     &            'Compute model 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(2,isTvar(itrc)), 'CnormM(isTvar)',  &
     &            'Compute model normalization factors for tracer ',    &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#      endif
          END IF
#     endif
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(1,isFsur), 'CnormI(isFsur)',           &
     &            'Compute initial 2D RHO-normalization factors.'
            WRITE (out,70) Cnorm(1,isUbar), 'CnormI(isUbar)',           &
     &            'Compute initial 2D U-normalization factors.'
            WRITE (out,70) Cnorm(1,isVbar), 'CnormI(isVbar)',           &
     &            'Compute initial 2D V-normalization factors.'
#     ifdef SOLVE3D
            WRITE (out,70) Cnorm(1,isUvel), 'CnormI(isUvel)',           &
     &            'Compute initial 3D U-normalization factors.'
            WRITE (out,70) Cnorm(1,isVvel), 'CnormI(isVvel)',           &
     &            'Compute initial 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(1,isTvar(itrc)), 'CnormI(isTvar)',  &
     &            'Compute initial normalization factors for tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
          END IF
#     ifdef ADJUST_BOUNDARY
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,170) CnormB(isFsur,1:4), 'CnormB(isFsur)',       &
     &            'Compute boundary 2D RHO-normalization factors.'
            WRITE (out,170) CnormB(isUbar,1:4), 'CnormB(isUbar)',       &
     &            'Compute boundary 2D U-normalization factors.'
            WRITE (out,170) CnormB(isVbar,1:4), 'CnormB(isVbar)',       &
     &            'Compute initial 2D V-normalization factors.'
#      ifdef SOLVE3D
            WRITE (out,170) CnormB(isUvel,1:4), 'CnormB(isUvel)',       &
     &            'Compute initial 3D U-normalization factors.'
            WRITE (out,170) CnormB(isVvel,1:4), 'CnormI(isVvel)',       &
     &            'Compute initial 3D V-normalization factors.'
            DO itrc=1,NT(ng)
              WRITE (out,180) CnormB(isTvar(itrc),1:4),'CnormI(isTvar)',&
     &            'Compute initial normalization factors for tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#      endif
          END IF
#     endif
#     ifdef ADJUST_WSTRESS
          IF (ANY(LwrtNRM(:,ng))) THEN
            WRITE (out,70) Cnorm(1,isUstr), 'CnormF(isUstr)',           &
     &            'Compute normalization factors at surface U-stress.'
            WRITE (out,70) Cnorm(1,isVstr), 'CnormF(isVstr)',           &
     &            'Compute normalization factors at surface V-stress.'
          END IF
#     endif
#     if defined ADJUST_STFLUX && defined SOLVE3D
          IF (ANY(LwrtNRM(:,ng))) THEN
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(1,isTsur(itrc)), 'CnormF(isTsur)',  &
     &            'Compute normalization factors for flux of tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
          END IF
#     endif
#    endif
          WRITE (out,100) Hgamma(1), 'Hgamma',                          &
     &            'Horizontal diffusion factor, initial conditions.'
#    ifdef WEAK_CONSTRAINT
          WRITE (out,100) Hgamma(2), 'HgammaM',                         &
     &            'Horizontal diffusion factor, model error.'
#    endif
#    ifdef ADJUST_BOUNDARY
          WRITE (out,100) Hgamma(3), 'HgammaB',                         &
     &            'Horizontal diffusion factor, boundary conditions.'
#    endif
#    ifdef ADJUST_STFLUX
          WRITE (out,100) Hgamma(4), 'HgammaF',                         &
     &            'Horizontal diffusion factor, surface forcing.'
#    endif
#    ifdef SOLVE3D
          WRITE (out,100) Vgamma(1), 'Vgamma',                          &
     &            'Vertical diffusion factor, initial conditions.'
#     ifdef WEAK_CONSTRAINT
          WRITE (out,100) Vgamma(2), 'VgammaM',                         &
     &            'Vertical diffusion factor, model error.'
#     endif
#     ifdef ADJUST_BOUNDARY
          WRITE (out,100) Vgamma(3), 'VgammaB',                         &
     &            'Vertical diffusion factor, boundary conditions.'
#     endif
#    endif
#    if defined SENSITIVITY_4DVAR || \
        defined TL_W4DPSAS        || defined TL_W4DVAR || \
        defined W4DPSAS           || defined W4DVAR
         IF (nADJ(ng).lt.ntimes(ng)) THEN
            WRITE (out,120) Hdecay(2,isFsur,ng), 'HdecayM(isFsur)',     &
     &            'Model decorrelation H-scale (m), free-surface.'
            WRITE (out,120) Hdecay(2,isUbar,ng), 'HdecayM(isUbar)',     &
     &            'Model decorrelation H-scale (m), 2D U-momentum.'
            WRITE (out,120) Hdecay(2,isVbar,ng), 'HdecayM(isVbar)',     &
     &            'Model decorrelation H-scale (m), 2D V-momentum.'
#     ifdef SOLVE3D
            WRITE (out,120) Hdecay(2,isUvel,ng), 'HdecayM(isUvel)',     &
     &            'Model decorrelation H-scale (m), 3D U-momentum.'
            WRITE (out,120) Hdecay(2,isVvel,ng), 'HdecayM(isVvel)',     &
     &            'Model decorrelation H-scale (m), 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,130) Hdecay(2,isTvar(itrc),ng),                &
     &            'HdecayM(idTvar)',                                    &
     &            'Model decorrelation H-scale (m), ',                  &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
            WRITE (out,120) Vdecay(2,isUvel,ng), 'VdecayM(isUvel)',     &
     &            'Model decorrelation V-scale (m), 3D U-momentum.'
            WRITE (out,120) Vdecay(2,isVvel,ng), 'VdecayM(isVvel)',     &
     &            'Model decorrelation V-scale (m), 3D V-momentum.'
            DO itrc=1,NT(ng)
              WRITE (out,130) Vdecay(2,isTvar(itrc),ng),                &
     &            'VdecayM(idTvar)',                                    &
     &            'Model decorrelation V-scale (m), ',                  &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
          END IF
#    endif
          WRITE (out,120) Hdecay(1,isFsur,ng), 'HdecayI(isFsur)',       &
     &            'Initial decorrelation H-scale (m), free-surface.'
          WRITE (out,120) Hdecay(1,isUbar,ng), 'HdecayI(isUbar)',       &
     &            'Initial decorrelation H-scale (m), 2D U-momentum.'
          WRITE (out,120) Hdecay(1,isVbar,ng), 'HdecayI(isVbar)',       &
     &            'Initial decorrelation H-scale (m), 2D V-momentum.'
#    ifdef SOLVE3D
          WRITE (out,120) Hdecay(1,isUvel,ng), 'HdecayI(isUvel)',       &
     &            'Initial decorrelation H-scale (m), 3D U-momentum.'
          WRITE (out,120) Hdecay(1,isVvel,ng), 'HdecayI(isVvel)',       &
     &            'Initial decorrelation H-scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(1,isTvar(itrc),ng),                  &
     &            'HdecayI(idTvar)',                                    &
     &            'Initial decorrelation H-scale (m), ',                &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,120) Vdecay(1,isUvel,ng), 'VdecayI(isUvel)',       &
     &            'Initial decorrelation V-scale (m), 3D U-momentum.'
          WRITE (out,120) Vdecay(1,isVvel,ng), 'VdecayI(isVvel)',       &
     &            'Initial decorrelation V-scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Vdecay(1,isTvar(itrc),ng),                  &
     &            'VdecayI(idTvar)',                                    &
     &            'Initial decorrelation V-scale (m), ',                &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    ifdef ADJUST_BOUNDARY
          DO ib=1,4
            IF (ib.eq.iwest) THEN
              Text='W-bry  '
            ELSE IF (ib.eq.isouth) THEN
              Text='S-bry  '
            ELSE IF (ib.eq.ieast) THEN
              Text='E-bry  '
            ELSE IF (ib.eq.inorth) THEN
              Text='N-bry  '
            END IF
            IF (Lobc(ib,isFsur,ng)) THEN
              WRITE (out,120) HdecayB(isFsur,ib,ng), 'HdecayB(isFsur)', &
     &              Text//' decorrelation H-scale (m), free-surface.'
            END IF
            IF (Lobc(ib,isUbar,ng)) THEN
              WRITE (out,120) HdecayB(isUbar,ib,ng), 'HdecayB(isUbar)', &
     &              Text//' decorrelation H-scale (m), 2D U-momentum.'
            END IF
            IF (Lobc(ib,isVbar,ng)) THEN
              WRITE (out,120) HdecayB(isVbar,ib,ng), 'HdecayB(isVbar)', &
     &              Text//' decorrelation H-scale (m), 2D V-momentum.'
            END IF
#     ifdef SOLVE3D
            IF (Lobc(ib,isUvel,ng)) THEN
              WRITE (out,120) HdecayB(isUvel,ib,ng), 'HdecayB(isUvel)', &
     &              Text//' decorrelation H-scale (m), 3D U-momentum.'
            END IF
            IF (Lobc(ib,isVvel,ng)) THEN
              WRITE (out,120) HdecayB(isVvel,ib,ng), 'HdecayB(isVvel)', &
     &              Text//' decorrelation H-scale (m), 3D V-momentum.'
            END IF
            DO i=1,NT(ng)
              IF (Lobc(ib,isTvar(i),ng)) THEN
                WRITE(out,130) HdecayB(isTvar(i),ib,ng),                &
     &              'HdecayB(idTvar)',                                  &
     &              Text//' decorrelation H-scale (m), ',               &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            IF (Lobc(ib,isUvel,ng)) THEN
              WRITE (out,120) VdecayB(isUvel,ib,ng), 'VdecayB(isUvel)', &
     &              Text//' decorrelation V-scale (m), 3D U-momentum.'
            END IF
            IF (Lobc(ib,isVvel,ng)) THEN
              WRITE (out,120) VdecayB(isVvel,ib,ng), 'VdecayB(isVvel)', &
     &              Text//' decorrelation V-scale (m), 3D V-momentum.'
            END IF
            DO i=1,NT(ng)
              IF (Lobc(ib,isTvar(i),ng)) THEN
                WRITE(out,130) VdecayB(isTvar(i),ib,ng),                &
     &              'VdecayB(idTvar)',                                  &
     &              Text//' decorrelation V-scale (m), ',               &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
#     endif
          END DO
#    endif
#    ifdef ADJUST_WSTRESS
          WRITE (out,120) Hdecay(1,isUstr,ng), 'HdecayF(isUstr)',       &
     &            'Forcing decorrelation H-scale (m), U-stress.'
          WRITE (out,120) Hdecay(1,isVstr,ng), 'HdecayF(isVstr)',       &
     &            'Forcing decorrelation H-scale (m), V-stress.'
#    endif
#    if defined ADJUST_STFLUX && defined SOLVE3D
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(1,isTsur(itrc),ng),                  &
     &            'HdecayF(idTsur)',                                    &
     &            'Forcing decorrelation H-scale (m), ',                &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    if defined ADJUST_STFLUX && defined SOLVE3D
          DO itrc=1,NT(ng)
            WRITE (out,110) Lstflux(itrc,ng), 'Lstflux(itrc)',          &
     &            'Adjusting surface flux of tracer ', itrc,            &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#    ifdef ADJUST_BOUNDARY
          WRITE (out,170) Lobc(1:4,isFsur,ng), 'Lobc(isFsur)',          &
     &            'Adjusting free-surface boundaries.'
          WRITE (out,170) Lobc(1:4,isUbar,ng), 'Lobc(isUbar)',          &
     &            'Adjusting 2D U-momentum boundaries.'
          WRITE (out,170) Lobc(1:4,isVbar,ng), 'Lobc(isVbar)',          &
     &            'Adjusting 2D V-momentum boundaries.'
#     ifdef SOLVE3D
          WRITE (out,170) Lobc(1:4,isUvel,ng), 'Lobc(isUvel)',          &
     &            'Adjusting 3D U-momentum boundaries.'
          WRITE (out,170) Lobc(1:4,isVvel,ng), 'Lobc(isVvel)',          &
     &            'Adjusting 3D V-momentum boundaries.'
          DO itrc=1,NT(ng)
            WRITE (out,180) Lobc(1:4,isTvar(itrc),ng),'Lobc(isTvar)',   &
     &            'Adjusting boundaries for tracer ', itrc,             &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#   endif
#  endif
!
!-----------------------------------------------------------------------
!  Report input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,150)
#  ifdef VERIFICATION
          WRITE (out,160) ' Verification Parameters File:  ',           &
     &                    TRIM(aparnam)
#  else
          WRITE (out,160) ' Assimilation Parameters File:  ',           &
     &                    TRIM(aparnam)
#  endif
#  ifdef FOUR_DVAR
#   if defined IS4DVAR          || defined OBS_SENSITIBITY || \
       defined OPT_OBSERVATIONS || defined WEAK_CONSTRAINT
#    if defined SENSITIVITY_4DVAR || \
        defined TL_W4DPSAS        || defined TL_W4DVAR || \
        defined W4DPSAS           || defined W4DVAR
          fname=STDname(2,ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,160) '               Model STD File:  ',           &
     &                    TRIM(fname)
#    endif
          fname=STDname(1,ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,160) '  Initial Conditions STD File:  ',           &
     &                    TRIM(fname)
#    ifdef ADJUST_BOUNDARY
          fname=STDname(3,ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,160) ' Boundary Conditions STD File:  ',           &
     &                    TRIM(fname)
#    endif
#    if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          fname=STDname(4,ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,160) '     Surface Forcing STD File:  ',           &
     &                    TRIM(fname)
#    endif
#   endif
#   if defined SENSITIVITY_4DVAR || \
       defined TL_W4DPSAS        || defined TL_W4DVAR || \
       defined W4DPSAS           || defined W4DVAR
          fname=NRMname(2,ng)
          WRITE (out,160) '              Model Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(2,ng)) THEN
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
          END IF
#   endif
          fname=NRMname(1,ng)
          WRITE (out,160) ' Initial Conditions Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(1,ng)) THEN
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
          END IF
#   ifdef ADJUST_BOUNDARY
          fname=NRMname(3,ng)
          WRITE (out,160) 'Boundary Conditions Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(3,ng)) THEN
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
          END IF
#   endif
#   if defined ADJUST_WSTRESS || defined ADJUST_STFLUX
          fname=NRMname(4,ng)
          WRITE (out,160) '    Surface Forcing Norm File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(4,ng)) THEN
            INQUIRE (FILE=TRIM(fname), EXIST=inhere)
            IF (.not.inhere) GO TO 30
          END IF
#   endif
#   if !(defined CORRELATION || defined CONVOLUTION || \
         defined OPT_OBSERVATIONS)
          fname=OBSname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,160) '            Observations File:  ',           &
     &                    TRIM(fname)
#   endif
#  endif
          WRITE (out,160) '    Input/Output Lanczos File:  ',           &
     &                    TRIM(LCZname(ng))
#  ifndef IS4DVAR_SENSITIVITY
          WRITE (out,160) '    Input/Output Hessian File:  ',           &
     &                    TRIM(HSSname(ng))
#  endif
#  ifdef VERIFICATION
          WRITE (out,160) '     Output verification File:  ',           &
     &                    TRIM(MODname(ng))
#  else
          WRITE (out,160) '            Output 4DVAR File:  ',           &
     &                    TRIM(MODname(ng))
#  endif
#  if defined WEAK_CONSTRAINT   && \
     (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
          WRITE (out,160) '  Output Posterior Error File:  ',           &
     &                    TRIM(ERRname(ng))
#  endif
          GO TO 40
  30      WRITE (out,200) TRIM(fname)
          exit_flag=4
          RETURN
  40      CONTINUE
        END DO
      END IF

  50  FORMAT (/,' READ_AssPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Assimilation Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)
  70  FORMAT (10x,l1,2x,a,t30,a)
  80  FORMAT (1x,i10,2x,a,t30,a)
  90  FORMAT (1x,i10,2x,a,t30,a,/,t32,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 120  FORMAT (f11.3,2x,a,t30,a)
 130  FORMAT (f11.3,2x,a,t30,a,a,'.')
 150  FORMAT (/,' Input Assimilation Files:',/)
 160  FORMAT (2x,a,a)
 170  FORMAT (3x,4(1x,l1),2x,a,t30,a)
 180  FORMAT (3x,4(1x,l1),2x,a,t30,a,i2.2,':',1x,a)
 190  FORMAT (/,' READ_ASSPAR - variable info not yet loaded, ', a)
 200  FORMAT (/,' READ_ASSPAR - could not find input file:  ',a)
# endif

      RETURN
      END SUBROUTINE read_AssPar
#endif

#ifdef FLOATS
      SUBROUTINE read_FloatsPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in input station parameters.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
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
      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status

      integer, dimension(Ngrids) :: ncount, nentry

      integer, allocatable :: Fcoor(:,:), Fcount(:,:), Ftype(:,:)

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: xfloat, yfloat, zfloat

      real(r8), dimension(100) :: Rval

      real(r8), allocatable :: Ft0(:,:), Fx0(:,:), Fy0(:,:), Fz0(:,:)
      real(r8), allocatable :: Fdt(:,:), Fdx(:,:), Fdy(:,:), Fdz(:,:)

      character (len=35) :: frmt
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in initial float locations.
!-----------------------------------------------------------------------
!
!  Notice I added one when allocating local scratch arrays to avoid
!  out of bounds in some compilers when reading the last blank line
!  which signal termination of input data.
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lfloats') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lfloats)
          ELSE IF (TRIM(KeyWord).eq.'FRREC') THEN
            Npts=load_i(Nval, Rval, Ngrids, frrec)
          ELSE IF (TRIM(KeyWord).eq.'NFLOATS') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nfloats)
          ELSE IF (TRIM(KeyWord).eq.'POS') THEN
            Npts=Nfloats(1)+1
            IF (Ngrids.gt.1) Npts=MAXVAL(Nfloats)+1
            IF (.not.allocated(Fcoor))  allocate ( Fcoor (Npts,Ngrids) )
            IF (.not.allocated(Fcount)) allocate ( Fcount(Npts,Ngrids) )
            IF (.not.allocated(Ftype))  allocate ( Ftype (Npts,Ngrids) )
            IF (.not.allocated(Ft0))    allocate ( Ft0(Npts,Ngrids) )
            IF (.not.allocated(Fx0))    allocate ( Fx0(Npts,Ngrids) )
            IF (.not.allocated(Fy0))    allocate ( Fy0(Npts,Ngrids) )
            IF (.not.allocated(Fz0))    allocate ( Fz0(Npts,Ngrids) )
            IF (.not.allocated(Fdt))    allocate ( Fdt(Npts,Ngrids) )
            IF (.not.allocated(Fdx))    allocate ( Fdx(Npts,Ngrids) )
            IF (.not.allocated(Fdy))    allocate ( Fdy(Npts,Ngrids) )
            IF (.not.allocated(Fdz))    allocate ( Fdz(Npts,Ngrids) )
            DO ng=1,Ngrids
              CALL allocate_floats (ng)
            END DO
            ncount(1:Ngrids)=0
            nentry(1:Ngrids)=0
            DO WHILE (.TRUE.)
              READ (inp,*,ERR=30,END=30) igrid,                         &
     &                                   Fcoor (nentry(igrid)+1,igrid), &
     &                                   Ftype (nentry(igrid)+1,igrid), &
     &                                   Fcount(nentry(igrid)+1,igrid), &
     &                                   Ft0(nentry(igrid)+1,igrid),    &
     &                                   Fx0(nentry(igrid)+1,igrid),    &
     &                                   Fy0(nentry(igrid)+1,igrid),    &
     &                                   Fz0(nentry(igrid)+1,igrid),    &
     &                                   Fdt(nentry(igrid)+1,igrid),    &
     &                                   Fdx(nentry(igrid)+1,igrid),    &
     &                                   Fdy(nentry(igrid)+1,igrid),    &
     &                                   Fdz(nentry(igrid)+1,igrid)
              IF (igrid.gt.Ngrids) THEN
                IF (Lwrite) WRITE (out,40) fposnam
                exit_flag=4
                RETURN
              END IF
              ncount(igrid)=ncount(igrid)+Fcount(nentry(igrid)+1,igrid)
              nentry(igrid)=nentry(igrid)+1
            END DO
          END IF
        END IF
      END DO
  20  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!  Turn off the processing of floats if not running long enough to
!  create a floats file (LdefFLT=.FALSE. because nFLT < ntimes or
!  nFLT = 0 when nrrec = 0).
!
      DO ng=1,Ngrids
        IF (.not.LdefFLT(ng).and.Lfloats(ng)) THEN
          Lfloats(ng)=.FALSE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lfloats(ng)) THEN
            IF (ncount(ng).ne.Nfloats(ng)) THEN
              WRITE (stdout,60) ncount(ng), Nfloats(ng)
              exit_flag=4
              RETURN
            END IF
            WRITE (out,70) ng
            DO i=1,nentry(ng)
              IF (.not.spherical.and.(Fcoor(i,ng).eq.0)) THEN
                frmt='(i1,i2,i5,f10.4,2f8.2,f8.2,4f9.3)'
              ELSE
                frmt='(i1,i2,i5,f10.4,3f8.2,4f9.3)'
              END IF
              WRITE (out,frmt) Fcoor(i,ng), Ftype(i,ng), Fcount(i,ng),  &
     &                         Ft0(i,ng), Fx0(i,ng), Fy0(i,ng),         &
     &                         Fz0(i,ng), Fdt(i,ng), Fdx(i,ng),         &
     &                         Fdy(i,ng), Fdz(i,ng)
            END DO
            WRITE (out,80) Nfloats(ng),                                 &
     &            'Nfloats',                                            &
     &            'Number of float trajectories to compute.'
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Process initial float locations.
!-----------------------------------------------------------------------
!
!  Set time of float release (seconds after model initialization) and
!  initial float horizontal positions (grid units).  Fill the initial
!  vertical level or depth position.
!
      DO ng=1,Ngrids
        mc=0
        nc=0
        IF (Lfloats(ng)) THEN
          DO i=1,nentry(ng)
            IF (Fcount(i,ng).eq.1) THEN
              nc=nc+1
              FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
              FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
              FLT(ng)%Ftype(nc)=Ftype(i,ng)
              IF (Fcoor(i,ng).eq.0) THEN
                FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),      &
     &                                  REAL(Lm(ng),r8)+0.5_r8)
                FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),      &
     &                                  REAL(Mm(ng),r8)+0.5_r8)
              ELSE
                mc=mc+1
                FLT(ng)%Flon(mc)=Fx0(i,ng)
                FLT(ng)%Flat(mc)=Fy0(i,ng)
                FLT(ng)%Findex(mc)=nc
              END IF
            ELSE IF (Fcount(i,ng).gt.1) THEN
              DO j=1,Fcount(i,ng)
                nc=nc+1
                IF (Fdt(i,ng).gt.0.0_r8) THEN
                  FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+            &
     &                                     REAL(j-1,r8)*Fdt(i,ng))*     &
     &                                    day2sec
                  FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                  FLT(ng)%Ftype(nc)=Ftype(i,ng)
                  IF (Fcoor(i,ng).eq.0) THEN
                    FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),  &
     &                                      REAL(Lm(ng),r8)+0.5_r8)
                    FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),  &
     &                                      REAL(Mm(ng),r8)+0.5_r8)
                  ELSE
                    mc=mc+1
                    FLT(ng)%Flon(mc)=Fx0(i,ng)
                    FLT(ng)%Flat(mc)=Fy0(i,ng)
                    FLT(ng)%Findex(mc)=nc
                  END IF
                ELSE
                  FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
                  IF (Fdz(i,ng).eq.0.0_r8) THEN
                    FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                  ELSE
                    IF (Fz0(i,ng).gt.0.0_r8) THEN
                      zfloat=Fz0(i,ng)+REAL(j-1,r8)*Fdz(i,ng)
                      FLT(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,zfloat),   &
     &                                        REAL(N(ng),r8))
                    ELSE
                      FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+                &
     &                                        REAL(j-1,r8)*Fdz(i,ng)
                    END IF
                  END IF
                  FLT(ng)%Ftype(nc)=Ftype(i,ng)
                  IF (Fcoor(i,ng).eq.0) THEN
                    xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                    FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,xfloat),     &
     &                                      REAL(Lm(ng),r8)+0.5_r8)
                    yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                    FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,yfloat),     &
     &                                      REAL(Mm(ng),r8)+0.5_r8)
                  ELSE
                    mc=mc+1
                    FLT(ng)%Flon(mc)=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                    FLT(ng)%Flat(mc)=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                    FLT(ng)%Findex(mc)=nc
                  END IF
                END IF
              END DO
            END IF
          END DO
          FLT(ng)%Findex(0)=mc
        END IF
      END DO
!
!  Deallocate local arrays.
!
      IF (allocated(Fcoor))  deallocate ( Fcoor )
      IF (allocated(Fcount)) deallocate ( Fcount )
      IF (allocated(Ftype))  deallocate ( Ftype )
      IF (allocated(Ft0))    deallocate ( Ft0 )
      IF (allocated(Fx0))    deallocate ( Fx0 )
      IF (allocated(Fy0))    deallocate ( Fy0 )
      IF (allocated(Fz0))    deallocate ( Fz0 )
      IF (allocated(Fdt))    deallocate ( Fdt )
      IF (allocated(Fdx))    deallocate ( Fdx )
      IF (allocated(Fdy))    deallocate ( Fdy )
      IF (allocated(Fdz))    deallocate ( Fdz )
!
  40  FORMAT (/,' READ_FloatsPar - Error while reading floats',         &
     &          ' locations in input script: ',a)
  50  FORMAT (/,' READ_FloatsPar - Error while processing line: ',/,a)
  60  FORMAT (/,' READ_FloatsPar - Inconsistent number of floats to',   &
     &          ' process: ', 2i6,/,18x,'change input script.')
  70  FORMAT (/,/,' Floats Initial Locations, Grid: ',i2.2,             &
     &        /,  ' ==================================',/,/,            &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 80   FORMAT (/,1x,i10,2x,a,t30,a)

      RETURN
      END SUBROUTINE read_FloatsPar
#endif

#ifdef STATIONS
      SUBROUTINE read_StaPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in input station parameters.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
      USE mod_scalars
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
      integer :: Mstation, Npts, Nval
      integer :: flag, i, igrid, ista, itrc, ng, status

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: Xpos, Ypos

# if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
      logical, dimension(MBOTP,Ngrids) :: Lbottom
# endif
      logical, dimension(MT,Ngrids) :: Lsta

      integer, dimension(Ngrids) :: is

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in stations parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lstations') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lstations)
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idOvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idOvel,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbar,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbar,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idFsur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idFsur,:))
# if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(KeyWord).eq.'Sout(idBath)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idBath,:))
# endif
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTvar)') THEN
            Npts=load_l(Nval, Cval, MT*Ngrids, Lsta)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Sout(idTvar(itrc),ng)=Lsta(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUsms,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVsms,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbms,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbms,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbrs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbrs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbrs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbrs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbws)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbws,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbws)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbws,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbcs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbcs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbcs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbcs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbot,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbot,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbur,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbvr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbvr,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW2xx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW2xx,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW2xy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW2xy,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW2yy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW2yy,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU2rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU2rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV2rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV2rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU2Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU2Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV2Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV2Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3xx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3xx,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3xy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3xy,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3yy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3yy,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3zx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3zx,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3zy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3zy,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU3rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU3rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV3rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV3rs,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU3Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU3Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV3Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV3Sd,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWamp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWamp,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWlen)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWlen,:))
          ELSE IF (TRIM(Keyword).eq.'Sout(idWlep)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWlep,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWdir)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWdir,:))
# ifdef WAVES_TOP_PERIOD
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWptp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWptp,:))
# endif
# ifdef WAVES_BOT_PERIOD
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWpbt)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWpbt,:))
# endif
# ifdef WAVES_UB
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWorb)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWorb,:))
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWdif,:))
# endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWdib)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWdib,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWdiw)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWdiw,:))
# endif
#ifdef ROLLER_SVENDSEN
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWbrk)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWbrk,:))
# endif
#ifdef WEC_ROLLER
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWdis)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWdis,:))
# endif
#ifdef ROLLER_RENIERS
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWrol)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWrol,:))
# endif
#ifdef WEC_VF
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWztw)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWztw,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWqsp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWqsp,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWbeh)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWbeh,:))
# endif
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NAT*Ngrids, Lsta)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                Sout(idTsur(itrc),ng)=Lsta(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLhea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLhea,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idShea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idShea,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLrad,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSrad,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idEmPf)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idEmPf,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idevap)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idevap,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idrain)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idrain,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idDano)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idDano,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVvis)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVvis,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTdif,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSdif,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHsbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHsbl,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHbbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHbbl,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idMtke)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idMtke,:))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idMtls)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idMtls,:))
# if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
          ELSE IF (TRIM(KeyWord).eq.'Sout(idBott)') THEN
            Npts=load_l(Nval, Cval, MBOTP*Ngrids, Lbottom)
            DO ng=1,Ngrids
              DO itrc=1,MBOTP
                i=idBott(itrc)
                Sout(i,ng)=Lbottom(itrc,ng)
              END DO
            END DO
# endif
          ELSE IF (TRIM(KeyWord).eq.'NSTATION') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nstation)
          ELSE IF (TRIM(KeyWord).eq.'POS') THEN
            DO ng=1,Ngrids
              allocate ( SCALARS(ng) % Sflag(Nstation(ng)) )
              allocate ( SCALARS(ng) % SposX(Nstation(ng)) )
              allocate ( SCALARS(ng) % SposY(Nstation(ng)) )
            END DO
            is(1:Ngrids)=0
            DO WHILE (.TRUE.)
              READ (inp,*,ERR=10,END=10) igrid, flag, Xpos, Ypos
              ng=MAX(1,ABS(igrid))
              is(ng)=is(ng)+1
              SCALARS(ng)%Sflag(is(ng))=flag
              SCALARS(ng)%SposX(is(ng))=Xpos
              SCALARS(ng)%SposY(is(ng))=Ypos
            END DO
 10         DO ng=1,Ngrids
              IF (Nstation(ng).ne.is(ng)) THEN
                IF (Master) WRITE (out,40) Nstation(ng), is(ng)
                exit_flag=4
                RETURN
              END IF
            END DO
          END IF
        END IF
      END DO
 20   IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
 30   CONTINUE
!
!  Turn off the processing of stations if not running long enough to
!  create a stations file (LdefSTA=.FALSE. because nSTA < ntimes or
!  nSTA = 0 when nrrec = 0).
!
      DO ng=1,Ngrids
        IF (.not.LdefSTA(ng).and.Lstations(ng)) THEN
          Lstations(ng)=.FALSE.
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lstations(ng)) THEN
            WRITE (out,60) ng
            WRITE (out,70) Nstation(ng), 'Nstation',                    &
     &            'Number of stations to write out into stations file.'
#if defined SEDIMENT && defined SED_MORPH
            IF (Sout(idbath,ng)) WRITE (out,80) Sout(idbath,ng),        &
     &          'Sout(idbath)',                                         &
     &          'Write out free-surface.'
#endif
            IF (Sout(idFsur,ng)) WRITE (out,80) Sout(idFsur,ng),        &
     &          'Sout(idFsur)',                                         &
     &          'Write out free-surface.'
            IF (Sout(idUbar,ng)) WRITE (out,80) Sout(idUbar,ng),        &
     &          'Sout(idUbar)',                                         &
     &          'Write out 2D U-momentum component.'
            IF (Sout(idVbar,ng)) WRITE (out,80) Sout(idVbar,ng),        &
     &          'Sout(idVbar)',                                         &
     &          'Write out 2D V-momentum component.'
# ifdef SOLVE3D
            IF (Sout(idUvel,ng)) WRITE (out,80) Sout(idUvel,ng),        &
     &          'Sout(idUvel)',                                         &
     &          'Write out 3D U-momentum component.'
            IF (Sout(idVvel,ng)) WRITE (out,80) Sout(idVvel,ng),        &
     &          'Sout(idVvel)',                                         &
     &          'Write out 3D V-momentum component.'
            IF (Sout(idWvel,ng)) WRITE (out,80) Sout(idWvel,ng),        &
     &          'Sout(idWvel)',                                         &
     &          'Write out W-momentum component.'
            IF (Sout(idOvel,ng)) WRITE (out,80) Sout(idOvel,ng),        &
     &          'Sout(idOvel)',                                         &
     &          'Write out omega vertical velocity.'
            DO itrc=1,NT(ng)
              IF (Sout(idTvar(itrc),ng)) WRITE (out,90)                 &
     &            Sout(idTvar(itrc),ng), 'Sout(idTvar)',                &
     &            'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
# endif
            IF (Sout(idUsms,ng)) WRITE (out,80) Sout(idUsms,ng),        &
     &          'Sout(idUsms)',                                         &
     &          'Write out surface U-momentum stress.'
            IF (Sout(idVsms,ng)) WRITE (out,80) Sout(idVsms,ng),        &
     &          'Sout(idVsms)',                                         &
     &          'Write out surface V-momentum stress.'
            IF (Sout(idUbms,ng)) WRITE (out,80) Sout(idUbms,ng),        &
     &          'Sout(idUbms)',                                         &
     &          'Write out bottom U-momentum stress.'
            IF (Sout(idVbms,ng)) WRITE (out,80) Sout(idVbms,ng),        &
     &          'Sout(idVbms)',                                         &
     &          'Write out bottom V-momentum stress.'
# ifdef BBL_MODEL
            IF (Sout(idUbrs,ng)) WRITE (out,80) Sout(idUbrs,ng),        &
     &          'Sout(idUbrs)',                                         &
     &          'Write out bottom U-current stress.'
            IF (Sout(idVbrs,ng)) WRITE (out,80) Sout(idVbrs,ng),        &
     &          'Sout(idVbrs)',                                         &
     &          'Write out bottom V-current stress.'
            IF (Sout(idUbws,ng)) WRITE (out,80) Sout(idUbws,ng),        &
     &          'Sout(idUbws)',                                         &
     &          'Write out wind-induced, bottom U-wave stress.'
            IF (Sout(idVbws,ng)) WRITE (out,80) Sout(idVbws,ng),        &
     &          'Sout(idVbws)',                                         &
     &          'Write out wind-induced, bottom V-wave stress.'
            IF (Sout(idUbcs,ng)) WRITE (out,80) Sout(idUbcs,ng),        &
     &          'Sout(idUbcs)',                                         &
     &          'Write out max wind + current, bottom U-wave stress.'
            IF (Sout(idVbcs,ng)) WRITE (out,80) Sout(idVbcs,ng),        &
     &          'Sout(idVbcs)',                                         &
     &          'Write out max wind + current, bottom V-wave stress.'
            IF (Sout(idUbot,ng)) WRITE (out,80) Sout(idUbot,ng),        &
     &          'Sout(idUbot)',                                         &
     &          'Write out bed wave orbital U-velocity.'
            IF (Sout(idVbot,ng)) WRITE (out,80) Sout(idVbot,ng),        &
     &          'Sout(idVbot)',                                         &
     &          'Write out bed wave orbital V-velocity.'
            IF (Sout(idUbur,ng)) WRITE (out,80) Sout(idUbur,ng),        &
     &          'Sout(idUbur)',                                         &
     &          'Write out bottom U-velocity above bed.'
            IF (Sout(idVbvr,ng)) WRITE (out,80) Sout(idVbvr,ng),        &
     &          'Sout(idVbvr)',                                         &
     &          'Write out bottom V-velocity above bed.'
# endif
# if defined WEC_MELLOR
            IF (Sout(idW2xx,ng)) WRITE (out,80) Sout(idW2xx,ng),        &
     &         'Sout(idW2xx)',                                          &
     &         'Write out 2D radiation stress, Sxx.'
            IF (Sout(idW2xy,ng)) WRITE (out,80) Sout(idW2xy,ng),        &
     &         'Sout(idW2xy)',                                          &
     &         'Write out 2D radiation stress, Sxy.'
            IF (Sout(idW2yy,ng)) WRITE (out,80) Sout(idW2yy,ng),        &
     &         'Sout(idW2yy)',                                          &
     &         'Write out 2D radiation stress, Syy.'
            IF (Sout(idU2rs,ng)) WRITE (out,80) Sout(idU2rs,ng),        &
     &         'Sout(idU2rs)',                                          &
     &         'Write out total 2D u-radiation stress.'
            IF (Sout(idV2rs,ng)) WRITE (out,80) Sout(idV2rs,ng),        &
     &         'Sout(idV2rs)',                                          &
     &         'Write out total 2D v-radiation stress.'
#  ifdef SOLVE3D
            IF (Sout(idW3xx,ng)) WRITE (out,80) Sout(idW3xx,ng),        &
     &         'Sout(idW3xx)',                                          &
     &         'Write out 3D horizonrtal radiation stress, Sxx.'
            IF (Sout(idW3xy,ng)) WRITE (out,80) Sout(idW3xy,ng),        &
     &         'Sout(idW3xy)',                                          &
     &         'Write out 3D horizonrtal radiation stress, Sxy.'
            IF (Sout(idW3yy,ng)) WRITE (out,80) Sout(idW3yy,ng),        &
     &         'Sout(idW3yy)',                                          &
     &         'Write out 3D horizonrtal radiation stress, Syy.'
            IF (Sout(idW3zx,ng)) WRITE (out,80) Sout(idW3zx,ng),        &
     &         'Sout(idW3zx)',                                          &
     &         'Write out 3D vertical radiation stress, Spx.'
            IF (Sout(idW3zy,ng)) WRITE (out,80) Sout(idW3zy,ng),        &
     &         'Sout(idW3zy)',                                          &
     &         'Write out 3D vertical radiation stress, Spy.'
            IF (Sout(idU3rs,ng)) WRITE (out,80) Sout(idU3rs,ng),        &
     &         'Sout(idU3rs)',                                          &
     &         'Write out total 3D u-radiation stress.'
            IF (Sout(idV3rs,ng)) WRITE (out,80) Sout(idV3rs,ng),        &
     &         'Sout(idV3rs)',                                          &
     &         'Write out total 3D v-radiation stress.'
#  endif
# endif
# ifdef WEC
            IF (Sout(idU2Sd,ng)) WRITE (out,80) Sout(idU2Sd,ng),        &
     &         'Sout(idU2Sd)',                                          &
     &         'Write out 2D u-momentum stokes velocity.'
            IF (Sout(idV2Sd,ng)) WRITE (out,80) Sout(idV2Sd,ng),        &
     &         'Sout(idV2Sd)',                                          &
     &         'Write out 2D v-momentum stokes velocity.'
            IF (Sout(idU3Sd,ng)) WRITE (out,80) Sout(idU3Sd,ng),        &
     &         'Sout(idU3Sd)',                                          &
     &         'Write out 3D u-momentum stokes velocity.'
            IF (Sout(idV3Sd,ng)) WRITE (out,80) Sout(idV3Sd,ng),        &
     &         'Sout(idV3Sd)',                                          &
     &         'Write out 3D v-momentum stokes velocity.'
# endif
# ifdef WAVES_HEIGHT
            IF (Sout(idWamp,ng)) WRITE (out,80) Sout(idWamp,ng),        &
     &         'Sout(idWamp)',                                          &
     &         'Write out wave height.'
# endif
# ifdef WAVES_LENGTH
            IF (Sout(idWlen,ng)) WRITE (out,80) Sout(idWlen,ng),        &
     &         'Sout(idWlen)',                                          &
     &         'Write out wave length.'
# endif
# ifdef WAVES_LENGTHP
            IF (Sout(idWlep,ng)) WRITE (out,80) Sout(idWlep,ng),        &
     &         'Sout(idWlep)',                                          &
     &         'Write out peak wave length.'
# endif
# ifdef WAVES_DIR
            IF (Sout(idWdir,ng)) WRITE (out,80) Sout(idWdir,ng),        &
     &         'Sout(idWdir)',                                          &
     &         'Write out wave direction.'
# endif
# ifdef WAVES_TOP_PERIOD
            IF (Sout(idWptp,ng)) WRITE (out,80) Sout(idWptp,ng),        &
     &         'Sout(idWptp)',                                          &
     &         'Write out wave surface period.'
# endif
# ifdef WAVES_BOT_PERIOD
            IF (Sout(idWpbt,ng)) WRITE (out,80) Sout(idWpbt,ng),        &
     &         'Sout(idWpbt)',                                          &
     &         'Write out wave bottom period.'
# endif
# ifdef WAVES_UB
            IF (Sout(idWorb,ng)) WRITE (out,80) Sout(idWorb,ng),        &
     &         'Sout(idWorb)',                                          &
     &         'Write out wave bottom orbital velocity.'
# endif
# if defined WAVES_OCEAN || (defined WEC_VF && defined BOTTOM_STREAMING)
            IF (Sout(idWdif,ng)) WRITE (out,80) Sout(idWdif,ng),        &
     &         'Sout(idWdif)',                                          &
     &         'Write out wave dissipation due to bottom friction.'
# endif
#if defined WAVES_OCEAN || defined TKE_WAVEDISS || \
    defined WDISS_THORGUZA || defined WDISS_CHURTHOR
            IF (Sout(idWdib,ng)) WRITE (out,80) Sout(idWdib,ng),        &
     &         'Sout(idWdib)',                                          &
     &         'Write out wave dissipation due to breaking.'
            IF (Sout(idWdiw,ng)) WRITE (out,80) Sout(idWdiw,ng),        &
     &         'Sout(idWdiw)',                                          &
     &         'Write out wave dissipation due to white capping.'
# endif
#ifdef ROLLER_SVENDSEN
            IF (Sout(idWbrk,ng)) WRITE (out,80) Sout(idWbrk,ng),        &
     &         'Sout(idWbrk)',                                          &
     &         'Write out percent wave breaking.'
# endif
#ifdef WEC_ROLLER
            IF (Sout(idWdis,ng)) WRITE (out,80) Sout(idWdis,ng),        &
     &         'Sout(idWdis)',                                          &
     &         'Write out wave roller dissipation.'
# endif
#ifdef ROLLER_RENIERS
            IF (Sout(idWrol,ng)) WRITE (out,80) Sout(idWrol,ng),        &
     &         'Sout(idWrol)',                                          &
     &         'Write out wave roller action density.'
# endif
#ifdef WEC_VF
            IF (Sout(idWztw,ng)) WRITE (out,80) Sout(idWztw,ng),        &
     &         'Sout(idWztw)',                                          &
     &         'Write out WEC quasi-static sea level adjustment.'
            IF (Sout(idWqsp,ng)) WRITE (out,80) Sout(idWqsp,ng),        &
     &         'Sout(idWqsp)',                                          &
     &         'Write out WEC quasi-static pressure.'
            IF (Sout(idWbeh,ng)) WRITE (out,80) Sout(idWbeh,ng),        &
     &         'Sout(idWbeh)',                                          &
     &         'Write out WEC Bernoulli head.'
# endif
# if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
            DO itrc=1,MBOTP
              IF (Sout(idBott(itrc),ng)) WRITE (out,80)                 &
     &            Sout(idBott(itrc),ng), 'Sout(idBott)',                &
     &            'Write out bottom property ', itrc,                   &
     &            TRIM(Vname(1,idBott(itrc)))
            END DO
# endif
# ifdef SOLVE3D
            IF (Sout(idTsur(itemp),ng)) WRITE (out,80)                  &
     &          Sout(idTsur(itemp),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net heat flux.'
            IF (Sout(idTsur(isalt),ng)) WRITE (out,80)                  &
     &          Sout(idTsur(isalt),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net salt flux.'
#  ifdef SHORTWAVE
            IF (Sout(idSrad,ng)) WRITE (out,80) Sout(idSrad,ng),        &
     &          'Sout(idSrad)',                                         &
     &          'Write out shortwave radiation flux.'
#  endif
#  ifdef BULK_FLUXES
            IF (Sout(idLrad,ng)) WRITE (out,80) Sout(idLrad,ng),        &
     &          'Sout(idLrad)',                                         &
     &          'Write out longwave radiation flux.'
            IF (Sout(idLhea,ng)) WRITE (out,80) Sout(idLhea,ng),        &
     &          'Sout(idLhea)',                                         &
     &          'Write out latent heat flux.'
            IF (Sout(idShea,ng)) WRITE (out,80) Sout(idShea,ng),        &
     &          'Sout(idShea)',                                         &
     &          'Write out sensible heat flux.'
#   ifdef EMINUSP
            IF (Sout(idEmPf,ng)) WRITE (out,80) Sout(idEmPf,ng),        &
     &         'Sout(idEmPf)',                                          &
     &         'Write out E-P flux.'
            IF (Sout(idevap,ng)) WRITE (out,80) Sout(idevap,ng),        &
     &         'Sout(idevap)',                                          &
     &         'Write out evaporation rate.'
            IF (Sout(idrain,ng)) WRITE (out,80) Sout(idrain,ng),        &
     &         'Sout(idrain)',                                          &
     &         'Write out rain rate.'
#   endif
#  endif
            IF (Sout(idDano,ng)) WRITE (out,80) Sout(idDano,ng),        &
     &          'Sout(idDano)',                                         &
     &          'Write out density anomaly.'
            IF (Sout(idVvis,ng)) WRITE (out,80) Sout(idVvis,ng),        &
     &          'Sout(idVvis)',                                         &
     &          'Write out vertical viscosity coefficient.'
            IF (Sout(idTdif,ng)) WRITE (out,80) Sout(idTdif,ng),        &
     &          'Sout(idTdif)',                                         &
     &          'Write out vertical T-diffusion coefficient.'
            IF (Sout(idSdif,ng)) WRITE (out,80) Sout(idSdif,ng),        &
     &          'Sout(idSdif)',                                         &
     &          'Write out vertical S-diffusion coefficient.'
#  ifdef LMD_SKPP
            IF (Sout(idHsbl,ng)) WRITE (out,80) Sout(idHsbl,ng),        &
     &          'Sout(idHsbl)',                                         &
     &          'Write out depth of surface boundary layer.'
#  endif
#  ifdef LMD_BKPP
            IF (Sout(idHbbl,ng)) WRITE (out,80) Sout(idHbbl,ng),        &
     &          'Sout(idHbbl)',                                         &
     &          'Write out depth of bottom boundary layer.'
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
            IF (Sout(idMtke,ng)) WRITE (out,80) Sout(idMtke,ng),        &
     &          'Sout(idMtke)',                                         &
     &          'Write out turbulent kinetic energy.'
            IF (Sout(idMtls,ng)) WRITE (out,80) Sout(idMtls,ng),        &
     &          'Sout(idMtls)',                                         &
     &          'Write out turbulent generic length-scale.'
#  endif
# endif
            WRITE (out,*)
            DO i=1,Nstation(ng)
              WRITE (out,100) i, SCALARS(ng)%Sflag(i),                  &
     &                           SCALARS(ng)%SposX(i),                  &
     &                           SCALARS(ng)%SposY(i)
            END DO
          END IF
        END DO
      END IF

  40  FORMAT (/,' READ_StaPar - Inconsistent number of stations, ',     &
     &        'Nstation = ',2i8,/,15x,'change input script values.')
  50  FORMAT (/,' READ_StaPar - Error while processing line: ',/,a)

  60  FORMAT (/,/,' Stations Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  70  FORMAT (1x,i10,2x,a,t30,a)
  80  FORMAT (10x,l1,2x,a,t30,a)
  90  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 100  FORMAT (13x,'Flag and positions for station ',i4.4,':',           &
     &        i3,1x,2f10.4)
 110  FORMAT (/,' READ_StaPAR - variable info not yet loaded, ', a)

      RETURN
      END SUBROUTINE read_StaPar
#endif

      FUNCTION decode_line (line_text, KeyWord, Nval, Cval, Rval)
!
!=======================================================================
!                                                                      !
!  This function decodes lines of text from input script files.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
! Imported variable declarations.
!
      character (len=*), intent(in) :: line_text

      character (len=40), intent(inout) :: KeyWord

      integer, intent(inout) :: Nval

      character (len=160), dimension(100), intent(inout) :: Cval

      real(r8), dimension(100), intent(inout) :: Rval
!
! Local variable declarations
!
      logical :: IsString, Kextract, decode, nested
      integer :: Iblank, Icmm, Kstr, Kend, Linp
      integer :: Lend, LenS, Lstr, Lval, Nmul, Schar
      integer :: copies, i, ic, ie, is, j, status

      integer, dimension(20) :: Imul

      integer :: decode_line

      character (len=1 ), parameter :: blank = ' '
      character (len=160) :: Vstring, line, string
!
!------------------------------------------------------------------------
!  Decode input line.
!------------------------------------------------------------------------
!
!  Initialize.
!
      DO i=1,LEN(line)
        line(i:i)=blank
        Vstring(i:i)=blank
        string(i:i)=blank
      END DO
!
!  Get length of "line". Remove comment after the KEYWORD, if any.
!  Then, remove leading and trailing blanks.
!
      Linp=LEN(line_text)
      IF ((Linp.gt.0).and.(line_text(1:1).ne.CHAR(33))) THEN
        Icmm=INDEX(line_text,CHAR(33),BACK=.FALSE.)
        IF (Icmm.gt.0) Linp=Icmm-1
        line=TRIM(ADJUSTL(line_text(1:Linp)))
        Linp=LEN_TRIM(line)
      ELSE
        line=TRIM(ADJUSTL(line_text))
        Linp=LEN_TRIM(line)
      END IF
!
!  If not a blank or comment line [char(33)=!], decode and extract input
!  values.  Find equal sign [char(61)].
!
      status=-1
      nested=.FALSE.
      IF ((Linp.gt.0).and.(line(1:1).ne.CHAR(33))) THEN
        status=1
        Kstr=1
        Kend=INDEX(line,CHAR(61),BACK=.FALSE.)-1
        Lstr=INDEX(line,CHAR(61),BACK=.TRUE.)+1
!
! Determine if KEYWORD is followed by double equal sign (==) indicating
! nested parameter.
!
        IF ((Lstr-Kend).eq.3) nested=.TRUE.
!
! Extract KEYWORD, trim leading and trailing blanks.
!
        Kextract=.FALSE.
        IF (Kend.gt.0) THEN
          Lend=Linp
          KeyWord=line(Kstr:Kend)
          Nval=0
          Kextract=.TRUE.
        ELSE
          Lstr=1
          Lend=Linp
          Kextract=.TRUE.
        END IF
!
! Extract parameter values string.  Remove continuation symbol
! [char(92)=\], if any.  Trim leading trailing blanks.
!
        IF (Kextract) THEN
          Icmm=INDEX(line,CHAR(92),BACK=.FALSE.)
          IF (Icmm.gt.0) Lend=Icmm-1
          Vstring=ADJUSTL(line(Lstr:Lend))
          Lval=LEN_TRIM(Vstring)
!
! The TITLE KEYWORD is a special one since it can include strings,
! numbers, spaces, and continuation symbol.
!
          IsString=.FALSE.
          IF (TRIM(KeyWord).eq.'TITLE') THEN
            Nval=Nval+1
            Cval(Nval)=Vstring(1:Lval)
            IsString=.TRUE.
          ELSE
!
! Check if there is a multiplication symbol [char(42)=*] in the variable
! string indicating repetition of input values.
!
            Nmul=0
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(42)) THEN
                Nmul=Nmul+1
                Imul(Nmul)=i
              END IF
            END DO
            ic=1
!
! Check for blank spaces [char(32)=' '] between entries and decode.
!
            is=1
            ie=Lval
            Iblank=0
            decode=.FALSE.
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(32)) THEN
                IF (Vstring(i+1:i+1).ne.CHAR(32)) decode=.TRUE.
                Iblank=i
              ELSE
                ie=i
              ENDIF
              IF (decode.or.(i.eq.Lval)) THEN
                Nval=Nval+1
!
! Processing numeric values.  Check starting character to determine
! if numeric or character values. It is possible to have both when
! processing repetitions via the multiplication symbol.
!
                Schar=ICHAR(Vstring(is:is))
                IF (((48.le.Schar).and.(Schar.le.57)).or.               &
     &              (Schar.eq.43).or.(Schar.eq.45)) THEN
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Schar=ICHAR(Vstring(Imul(ic)+1:Imul(ic)+1))
                    IF ((43.le.Schar).and.(Schar.le.57)) THEN
                      READ (Vstring(Imul(ic)+1:ie),*) Rval(Nval)
                      DO j=1,copies-1
                        Rval(Nval+j)=Rval(Nval)
                      END DO
                    ELSE
                      string=Vstring(Imul(ic)+1:ie)
                      LenS=LEN_TRIM(string)
                      Cval(Nval)=string(1:LenS)
                      DO j=1,copies-1
                        Cval(Nval+j)=Cval(Nval)
                      END DO
                    END IF
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    LenS=LEN_TRIM(string)
                    READ (string(1:LenS),*) Rval(Nval)
                  END IF
                ELSE
!
! Processing character values (logicals and strings).
!
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Cval(Nval)=Vstring(Imul(ic)+1:ie)
                    DO j=1,copies-1
                      Cval(Nval+j)=Cval(Nval)
                    END DO
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    Cval(Nval)=TRIM(ADJUSTL(string))
                  END IF
                  IsString=.TRUE.
                END IF
                is=Iblank+1
                ie=Lval
                decode=.FALSE.
              END IF
            END DO
          END IF
        END IF
        status=Nval
      END IF
      decode_line=status
      RETURN
      END FUNCTION decode_line

      FUNCTION load_i (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model integer     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     load_i     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(r8), intent(in) :: Vinp(Ninp)
      integer, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_i
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(Ninp))
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
      END IF
      load_i=ic

      RETURN
      END FUNCTION load_i

      FUNCTION load_l (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model logical     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     load_l     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      character (len=*), intent(in) :: Vinp(Ninp)
      logical, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_l
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.(Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          IF ((Vinp(Ninp)(1:1).eq.'T').or.(Vinp(Ninp)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.(Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
      END IF
      load_l=ic

      RETURN
      END FUNCTION load_l

      FUNCTION load_r (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!                                                                      !
!  This function loads input values into a requested model real        !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output real variable.                                 !
!     load_r     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(r8), intent(in) :: Vinp(Ninp)
      real(r8), intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_r
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=Vinp(Ninp)
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
      END IF
      load_r=ic

      RETURN
      END FUNCTION load_r
