#include "cppdefs.h"
      MODULE inp_par_mod
!
!svn $Id: inp_par.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in input model parameters from standard input.   !
!  It also writes out these parameters to standard output.             !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE inp_par (model)
!***********************************************************************
!

#ifdef FILTERED
      USE mod_filter, ONLY: nfil,fcount,nw,dtfil,FILN
#endif
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
#ifdef DISTRIBUTE
      USE mod_strings
#endif
# ifdef UV_WAVEDRAG
      USE mod_tides, ONLY : nstep_hour
# endif
!
      USE dateclock_mod,  ONLY : get_date
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE ran_state,      ONLY : ran_seed
      USE strings_mod,    ONLY : FoundError
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
!
      integer :: Itile, Jtile, Nghost, Ntiles, tile
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Uoff, Voff
#ifdef DISTRIBUTE
      integer :: MaxHaloLenI, MaxHaloLenJ
#endif
      integer :: ibry, inp, out, i, ic, ifield, itrc, j, ng, npts
      integer :: ifile, sequence, varid

      real(r8) :: cff
      real(r8), parameter :: epsilon = 1.0E-8_r8
      real(r8), parameter :: spv = 0.0_r8
      character*2 :: ngnumber
!
      character (len=10) :: stdout_file

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------

#ifdef ROMS_STDOUT
!
!  Change default Fortran standard out unit, so ROMS run information is
!  directed to a file. This is advantageous in coupling applications to
!  ROMS information separated from other models. It overwite Fortran
!  default unit 6.
!
# if defined DISTRIBUTE && defined DISJOINTED
       stdout=60+ForkColor
       WRITE (stdout_file,'(a,i2.2,a)') 'log', ForkColor+1, '.roms'
# else
       stdout=60
       stdout_file='log.roms'
# endif
!
!  Open ROMS standard output file.
!
       IF (Master) THEN
         IF (Lappend) THEN
           OPEN (stdout, FILE=TRIM(stdout_file), FORM='formatted',      &
     &           STATUS='old', POSITION='append', ACTION='write')
         ELSE
           OPEN (stdout, FILE=TRIM(stdout_file), FORM='formatted',      &
     &           STATUS='replace')
         END IF
       END IF
#endif
!
!  Set input units.
!
#if defined DISTRIBUTE || defined COAWST_COUPLING
      Lwrite=Master
      inp=1
      out=stdout
#else
      Lwrite=Master
      inp=stdinp
      out=stdout
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
#ifndef NO_DATE_STAMP
# ifdef DISTRIBUTE
      IF (Master) CALL get_date (date_str)
      CALL mp_bcasts (1, model, date_str)
# else
      CALL get_date (date_str)
# endif
#endif
!
!-----------------------------------------------------------------------
!  Read in physical model input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) WRITE (out,10) version, TRIM(date_str)
 10   FORMAT (80('-'),/,                                                &
              ' Model Input Parameters:  ROMS/TOMS version ',a,/,       &
     &        26x,a,/,80('-'))

#if defined COAWST_COUPLING || defined JEDI
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
!!    CALL my_getarg (1, Iname)
      IF (Master) CALL my_getarg (1, Iname)
      CALL mp_bcasts (1, model, Iname)
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
     &        /,11x,'mpirun -np 4 romsM roms.in',/,/,11x,'and not',/,   &
     &        /,11x,'mpirun -np 4 romsM < roms.in',/)
 40   CONTINUE
# endif
#endif
!
      CALL read_PhyPar (model, inp, out, Lwrite)
#ifdef DISTRIBUTE
      CALL mp_bcasti (1, model, exit_flag)
#endif
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

#if defined SPLIT_4DVAR && SUPPRESS_REPORT
!
!  If appending into standard output file, supress the reporting of
!  information.  Turn of "LwrtInfo" switch.
!  appending into standard output.
!
      IF (Lappend) THEN
        DO ng=1,Ngrids
          LwrtInfo(ng)=.FALSE.
        END DO
      END IF
#endif
#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Read in biological model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (15, FILE=TRIM(bparnam), FORM='formatted', STATUS='old')

      CALL read_BioPar (model, 15, out, Lwrite)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
#endif
#ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Read in sediment model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (25, FILE=TRIM(sparnam), FORM='formatted', STATUS='old')

      CALL read_SedPar (model, 25, out, Lwrite)
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
!  Set lower and upper bounds indices per domain partition for all
!  nested grids.
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
!  Set logical switches needed when processing variables in tiles
!  adjacent to the domain boundary edges or corners.  This needs to
!  be computed first since these switches are used in "get_tile".
!
      DO ng=1,Ngrids
        DO tile=-1,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain_edges (ng, tile,                              &
     &                           DOMAIN(ng) % Eastern_Edge    (tile),   &
     &                           DOMAIN(ng) % Western_Edge    (tile),   &
     &                           DOMAIN(ng) % Northern_Edge   (tile),   &
     &                           DOMAIN(ng) % Southern_Edge   (tile),   &
     &                           DOMAIN(ng) % NorthEast_Corner(tile),   &
     &                           DOMAIN(ng) % NorthWest_Corner(tile),   &
     &                           DOMAIN(ng) % SouthEast_Corner(tile),   &
     &                           DOMAIN(ng) % SouthWest_Corner(tile),   &
     &                           DOMAIN(ng) % NorthEast_Test  (tile),   &
     &                           DOMAIN(ng) % NorthWest_Test  (tile),   &
     &                           DOMAIN(ng) % SouthEast_Test  (tile),   &
     &                           DOMAIN(ng) % SouthWest_Test  (tile))
        END DO
      END DO
!
!  Set tile computational indices and arrays allocation bounds
!
      Nghost=NghostPoints
      DO ng=1,Ngrids
        BOUNDS(ng) % LBij = 0
        BOUNDS(ng) % UBij = MAX(Lm(ng)+1,Mm(ng)+1)
        DO tile=-1,NtileI(ng)*NtileJ(ng)-1
          BOUNDS(ng) % tile(tile) = tile
          CALL get_tile (ng, tile, Itile, Jtile,                        &
     &                   BOUNDS(ng) % Istr   (tile),                    &
     &                   BOUNDS(ng) % Iend   (tile),                    &
     &                   BOUNDS(ng) % Jstr   (tile),                    &
     &                   BOUNDS(ng) % Jend   (tile),                    &
     &                   BOUNDS(ng) % IstrM  (tile),                    &
     &                   BOUNDS(ng) % IstrR  (tile),                    &
     &                   BOUNDS(ng) % IstrU  (tile),                    &
     &                   BOUNDS(ng) % IendR  (tile),                    &
     &                   BOUNDS(ng) % JstrM  (tile),                    &
     &                   BOUNDS(ng) % JstrR  (tile),                    &
     &                   BOUNDS(ng) % JstrV  (tile),                    &
     &                   BOUNDS(ng) % JendR  (tile),                    &
     &                   BOUNDS(ng) % IstrB  (tile),                    &
     &                   BOUNDS(ng) % IendB  (tile),                    &
     &                   BOUNDS(ng) % IstrP  (tile),                    &
     &                   BOUNDS(ng) % IendP  (tile),                    &
     &                   BOUNDS(ng) % IstrT  (tile),                    &
     &                   BOUNDS(ng) % IendT  (tile),                    &
     &                   BOUNDS(ng) % JstrB  (tile),                    &
     &                   BOUNDS(ng) % JendB  (tile),                    &
     &                   BOUNDS(ng) % JstrP  (tile),                    &
     &                   BOUNDS(ng) % JendP  (tile),                    &
     &                   BOUNDS(ng) % JstrT  (tile),                    &
     &                   BOUNDS(ng) % JendT  (tile),                    &
     &                   BOUNDS(ng) % Istrm3 (tile),                    &
     &                   BOUNDS(ng) % Istrm2 (tile),                    &
     &                   BOUNDS(ng) % Istrm1 (tile),                    &
     &                   BOUNDS(ng) % IstrUm2(tile),                    &
     &                   BOUNDS(ng) % IstrUm1(tile),                    &
     &                   BOUNDS(ng) % Iendp1 (tile),                    &
     &                   BOUNDS(ng) % Iendp2 (tile),                    &
     &                   BOUNDS(ng) % Iendp2i(tile),                    &
     &                   BOUNDS(ng) % Iendp3 (tile),                    &
     &                   BOUNDS(ng) % Jstrm3 (tile),                    &
     &                   BOUNDS(ng) % Jstrm2 (tile),                    &
     &                   BOUNDS(ng) % Jstrm1 (tile),                    &
     &                   BOUNDS(ng) % JstrVm2(tile),                    &
     &                   BOUNDS(ng) % JstrVm1(tile),                    &
     &                   BOUNDS(ng) % Jendp1 (tile),                    &
     &                   BOUNDS(ng) % Jendp2 (tile),                    &
     &                   BOUNDS(ng) % Jendp2i(tile),                    &
     &                   BOUNDS(ng) % Jendp3 (tile))

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
      Nghost=NghostPoints
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
!  Set minimum and maximum fractional coordinates for processing
!  observations. Either the full grid or only interior points will
!  be considered.  The strategy here is to add a small value (epsilon)
!  to the eastern and northern boundary values of Xmax and Ymax so
!  observations at such boundaries locations are processed. This
!  is needed because the .lt. operator in the following conditional:
!
!     IF (...
!    &    ((Xmin.le.Xobs(iobs)).and.(Xobs(iobs).lt.Xmax)).and.          &
!    &    ((Ymin.le.Yobs(iobs)).and.(Yobs(iobs).lt.Ymax))) THEN
!-----------------------------------------------------------------------
!
!  Set RHO-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
#ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, r2dvar, 0, Itile, Jtile,           &
     &                   rILB(ng), rIUB(ng), rJLB(ng), rJUB(ng))
# ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          rILB(ng)=rILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          rIUB(ng)=rIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          rJLB(ng)=rJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          rJUB(ng)=rJUB(ng)-1
        END IF
# endif
#else
# ifdef FULL_GRID
        rILB(ng)=0
        rIUB(ng)=Lm(ng)+1
        rJLB(ng)=0
        rJUB(ng)=Mm(ng)+1
# else
        rILB(ng)=1
        rIUB(ng)=Lm(ng)
        rJLB(ng)=1
        rJUB(ng)=Mm(ng)
# endif
#endif
!
!  Minimum and maximum fractional coordinates for RHO-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, r2dvar, 0, epsilon,                &
#ifdef FULL_GRID
     &                     .TRUE.,                                      &
#else
     &                     .FALSE.,                                     &
#endif
     &                     DOMAIN(ng) % Xmin_rho(tile),                 &
     &                     DOMAIN(ng) % Xmax_rho(tile),                 &
     &                     DOMAIN(ng) % Ymin_rho(tile),                 &
     &                     DOMAIN(ng) % Ymax_rho(tile))
        END DO
#ifdef DISTRIBUTE
        rXmin(ng)=DOMAIN(ng)%Xmin_rho(MyRank)
        rXmax(ng)=DOMAIN(ng)%Xmax_rho(MyRank)
        rYmin(ng)=DOMAIN(ng)%Ymin_rho(MyRank)
        rYmax(ng)=DOMAIN(ng)%Ymax_rho(MyRank)
#else
        rXmin(ng)=DOMAIN(ng)%Xmin_rho(0)
        rXmax(ng)=DOMAIN(ng)%Xmax_rho(0)
        rYmin(ng)=DOMAIN(ng)%Ymin_rho(0)
        rYmax(ng)=DOMAIN(ng)%Ymax_rho(0)
#endif
      END DO
!
!  Set U-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
        IF (EWperiodic(ng)) THEN
          Uoff=0
        ELSE
          Uoff=1
        END IF
#ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, u2dvar, 0, Itile, Jtile,           &
     &                   uILB(ng), uIUB(ng), uJLB(ng), uJUB(ng))
# ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          uILB(ng)=uILB(ng)+Uoff
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          uIUB(ng)=uIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          uJLB(ng)=uJLB(ng)+1
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          uJUB(ng)=uJUB(ng)-1
        END IF
# endif
#else
# ifdef FULL_GRID
        uILB(ng)=1
        uIUB(ng)=Lm(ng)+1
        uJLB(ng)=0
        uJUB(ng)=Mm(ng)+1
# else
        uILB(ng)=1+Uoff
        uIUB(ng)=Lm(ng)
        uJLB(ng)=1
        uJUB(ng)=Mm(ng)
# endif
#endif
!
!  Minimum and maximum fractional coordinates for U-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, u2dvar, 0, epsilon,                &
#ifdef FULL_GRID
     &                     .TRUE.,                                      &
#else
     &                     .FALSE.,                                     &
#endif
     &                     DOMAIN(ng) % Xmin_u(tile),                   &
     &                     DOMAIN(ng) % Xmax_u(tile),                   &
     &                     DOMAIN(ng) % Ymin_u(tile),                   &
     &                     DOMAIN(ng) % Ymax_u(tile))
        END DO
#ifdef DISTRIBUTE
        uXmin(ng)=DOMAIN(ng)%Xmin_u(MyRank)
        uXmax(ng)=DOMAIN(ng)%Xmax_u(MyRank)
        uYmin(ng)=DOMAIN(ng)%Ymin_u(MyRank)
        uYmax(ng)=DOMAIN(ng)%Ymax_u(MyRank)
#else
        uXmin(ng)=DOMAIN(ng)%Xmin_u(0)
        uXmax(ng)=DOMAIN(ng)%Xmax_u(0)
        uYmin(ng)=DOMAIN(ng)%Ymin_u(0)
        uYmax(ng)=DOMAIN(ng)%Ymax_u(0)
#endif
      END DO
!
!  Set V-points domain lower and upper bounds (integer).
!
      DO ng=1,Ngrids
        IF (NSperiodic(ng)) THEN
          Voff=0
        ELSE
          Voff=1
        END IF
#ifdef DISTRIBUTE
        CALL get_bounds (ng, MyRank, v2dvar, 0, Itile, Jtile,           &
     &                   vILB(ng), vIUB(ng), vJLB(ng), vJUB(ng))
# ifndef FULL_GRID
        IF (Itile.eq.0) THEN
          vILB(ng)=vILB(ng)+1
        END IF
        IF (Itile.eq.(NtileI(ng)-1)) THEN
          vIUB(ng)=vIUB(ng)-1
        END IF
        IF (Jtile.eq.0) THEN
          vJLB(ng)=vJLB(ng)+Voff
        END IF
        IF (Jtile.eq.(NtileJ(ng)-1)) THEN
          vJUB(ng)=vJUB(ng)-1
        END IF
# endif
#else
#  ifdef FULL_GRID
        vILB(ng)=0
        vIUB(ng)=Lm(ng)+1
        vJLB(ng)=1
        vJUB(ng)=Mm(ng)+1
# else
        vILB(ng)=1
        vIUB(ng)=Lm(ng)
        vJLB(ng)=1+Voff
        vJUB(ng)=Mm(ng)
# endif
#endif
!
!  Minimum and maximum fractional coordinates for V-points.
!
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
          CALL get_domain (ng, tile, v2dvar, 0, epsilon,                &
#ifdef FULL_GRID
     &                     .TRUE.,                                      &
#else
     &                     .FALSE.,                                     &
#endif
     &                     DOMAIN(ng) % Xmin_v(tile),                   &
     &                     DOMAIN(ng) % Xmax_v(tile),                   &
     &                     DOMAIN(ng) % Ymin_v(tile),                   &
     &                     DOMAIN(ng) % Ymax_v(tile))
        END DO
#ifdef DISTRIBUTE
        vXmin(ng)=DOMAIN(ng)%Xmin_v(MyRank)
        vXmax(ng)=DOMAIN(ng)%Xmax_v(MyRank)
        vYmin(ng)=DOMAIN(ng)%Ymin_v(MyRank)
        vYmax(ng)=DOMAIN(ng)%Ymax_v(MyRank)
#else
        vXmin(ng)=DOMAIN(ng)%Xmin_v(0)
        vXmax(ng)=DOMAIN(ng)%Xmax_v(0)
        vYmin(ng)=DOMAIN(ng)%Ymin_v(0)
        vYmax(ng)=DOMAIN(ng)%Ymax_v(0)
#endif
      END DO
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
      OPEN (65, FILE=TRIM(vegnam), FORM='formatted', STATUS='old')

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
          CALL lbc_report (out, ifield, LBC)
        END DO

#if defined ADJOINT || defined TANGENT || defined TL_IOMS
!
        WRITE (out,130) 'TLM, RPM, and ADM'
        DO ifield=1,nLBCvar
          CALL lbc_report (out, ifield, ad_LBC)
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
#ifdef FILTERED
!
!  Check various parameters for writing time-filtered output
!
        IF (ntsAVG(ng).ne.1) THEN
          ntsAVG(ng) = 1
          IF (Master) print *, 'Resetting ntsAVG to 1'
        END IF
        IF (MOD(dt(ng),1.0_r8).ne.0) THEN
          IF (Master) print *, 'Model timestep must be a whole number'
          exit_flag=2
        END IF
        IF (MOD(dtfil,INT(dt(ng))).ne.0) THEN
          IF (Master) print *,                                          &
     &                'Model timesteps must fall on filter interval'
          exit_flag=2
        END IF
        nfil = nw*dtfil/INT(dt(ng))
        IF (MOD(nfil,navg(ng)).ne.0) THEN
          IF (Master) print *,                                          &
     &                'Must write out in sync with filter length'
          exit_flag=2
        END IF
        IF (nfil/navg(ng).ne.FILN) THEN
          IF (Master) print *, 'Number of filter sums is wrong'
          exit_flag=2
        END IF
        DO itrc=1,FILN
          fcount(itrc) = (1-itrc)*nw/FILN
        END DO
# ifdef FILTERED_RST
!
!  Hardcode the names of the filter restart files.
!
        WRITE(ngnumber,'(i2.2)') ng
        FIL(1,ng)%name=TRIM(ADJUSTL('ocean_fil2d_'//ngnumber//'.nc'))
        FIL(2,ng)%name=TRIM(ADJUSTL('ocean_filu_'//ngnumber//'.nc'))
        FIL(3,ng)%name=TRIM(ADJUSTL('ocean_filv_'//ngnumber//'.nc'))
        FIL(4,ng)%name=TRIM(ADJUSTL('ocean_filakt_'//ngnumber//'.nc'))
        FIL(5,ng)%name=TRIM(ADJUSTL('ocean_filaks_'//ngnumber//'.nc'))
        FIL(6,ng)%name=TRIM(ADJUSTL('ocean_filakv_'//ngnumber//'.nc'))
        FIL(7,ng)%name=TRIM(ADJUSTL('ocean_filr_'//ngnumber//'.nc'))
        FIL(8,ng)%name=TRIM(ADJUSTL('ocean_filt_'//ngnumber//'.nc'))
        FIL(9,ng)%name=TRIM(ADJUSTL('ocean_filw_'//ngnumber//'.nc'))

        DO ifile=1,nfile
          allocate ( FIL(ifile,ng)%Nrec(1) )
          allocate ( FIL(ifile,ng)%Vid(NV) )
          allocate ( FIL(ifile,ng)%Tid(MT) )
          allocate ( FIL(ifile,ng)%time_min(1) )
          allocate ( FIL(ifile,ng)%time_max(1) )
          allocate ( FIL(ifile,ng)%files(1) )
          FIL(ifile,ng)%Nfiles=1
          FIL(ifile,ng)%Rindex=0
          FIL(ifile,ng)%Fcount=1
          FIL(ifile,ng)%Nrec=0
          FIL(ifile,ng)%ncid=-1
          FIL(ifile,ng)%Vid=-1
          FIL(ifile,ng)%Tid=-1
          FIL(ifile,ng)%time_min=0.0_r8
          FIL(ifile,ng)%time_max=0.0_r8
          FIL(ifile,ng)%label='FILTER - tidal filtering'
        END DO
# endif
#endif
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
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options and definitions.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        CALL checkdefs
        CALL my_flush (out)
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
