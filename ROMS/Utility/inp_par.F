#include "cppdefs.h"
      SUBROUTINE inp_par (model)
!
!svn $Id: inp_par.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in input model parameters from standard input.   !
!  It also writes out these parameters to standard output.             !
!                                                                      !
!=======================================================================
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
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
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
      IF (Lwrite) WRITE (out,10) version, TRIM(date_str)
 10   FORMAT (/,' Model Input Parameters:  ROMS/TOMS version ',a,/,     &
     &       26x,a,/,1x,77('-'))

#ifdef COAWST_COUPLING
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
#ifdef NESTING
!
!  Read in nesting contact points NetCDF file and allocate and
!  initialize several structures and variables.
!
      CALL set_contact (1, model)
      IF (exit_flag.ne.NoError) RETURN
#endif
!
!-----------------------------------------------------------------------
!  Set lower and upper bounds indices per domain partition for all
!  nested grids.
!-----------------------------------------------------------------------
!
!  Determine the number of ghost-points in the halo region.
!
#if defined TS_MPDATA || defined UV_VIS4
      NghostPoints=3
#else
      NghostPoints=2
#endif
      IF (ANY(CompositeGrid).or.ANY(RefinedGrid)) THEN
        NghostPoints=MAX(3,NghostPoints)
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
      IF (Master) THEN
        DO ng=1,Ngrids
#ifdef SOLVE3D
          WRITE (stdout,50) ng, Lm(ng), Mm(ng), N(ng),                  &
     &                      NtileI(ng), NtileJ(ng)
#else
          WRITE (stdout,50) ng, Lm(ng), Mm(ng),                         &
     &                      NtileI(ng), NtileJ(ng)
#endif
          Print *,'Number of tracers: ', NT(ng)
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
#ifdef ICE_MODEL
!
!-----------------------------------------------------------------------
!  Read in ice model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (inp, FILE=TRIM(iparnam), FORM='formatted', STATUS='old')

      CALL read_IcePar (inp, out, Lwrite)
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
#if defined FOUR_DVAR || defined VERIFICATION
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

      CALL read_FltPar (model, 45, out, Lwrite)
#endif
#if defined FLOATS && defined FLOAT_BIOLOGY
!
!-----------------------------------------------------------------------
!  Read in biological float behavior model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (50, FILE=TRIM(fbionam), FORM='formatted', STATUS='old')

      CALL read_FltBioPar (model, 50, out, Lwrite)
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
#ifdef VEGETATION
!
!-----------------------------------------------------------------------
!  Read in vegetation model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (65, FILE=TRIM(vegnam), FORM='formatted', STATUS='old')

      CALL read_VegPar (model, 65, out, Lwrite)
#endif 
!
!-----------------------------------------------------------------------
!  Report lateral boundary conditions.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        WRITE (out,120) 'NLM'
 120    FORMAT (/,1x,'Lateral Boundary Conditions: ',a,/,1x,28('='),/,  &
     &          /,1x,'Variable',t25,'Grid',t31,'West Edge',             &
     &          t44,'South Edge', t57,'East Edge',t70,'North Edge',     &
     &          /,1x,'---------',t25,'----',t31,4('----------',3x))
        DO ifield=1,nLBCvar
          CALL lbc_report (out, ifield, LBC)
        END DO

#if defined ADJOINT || defined TANGENT || defined TL_IOMS
        WRITE (out,120) 'TLM, RPM, and ADM'
        DO ifield=1,nLBCvar
          CALL lbc_report (out, ifield, ad_LBC)
        END DO
#endif
      END IF
      IF (exit_flag.ne.NoError) RETURN
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
#ifdef ICE_MODEL
        FIL(10,ng)%name=TRIM(ADJUSTL('ocean_filice_'//ngnumber//'.nc'))
#endif
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
# if defined CLIMA_TS_MIX && (defined TS_DIF2 || defined TS_DIF4)
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
#ifdef NCEP_FLUXES
        cff=1.0_r8
#else
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
#endif
#ifdef NCEP_FLUXES
        cff=1.0_r8
#else
        cff=1.0_r8/(rho0*Cp)
        Fscale(idTsur(itemp),ng)=cff*Fscale(idTsur(itemp),ng)
        Fscale(idTbot(itemp),ng)=cff*Fscale(idTbot(itemp),ng)
        Fscale(idSrad,ng)=cff*Fscale(idSrad,ng)
        Fscale(idLdwn,ng)=cff*Fscale(idLdwn,ng)
        Fscale(idLrad,ng)=cff*Fscale(idLrad,ng)
        Fscale(idLhea,ng)=cff*Fscale(idLhea,ng)
        Fscale(idShea,ng)=cff*Fscale(idShea,ng)
        Fscale(iddQdT,ng)=cff*Fscale(iddQdT,ng)
#endif

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
            WRITE (stdout,130) MV, varid
 130        FORMAT (/,' INP_PAR - too small dimension ',                &
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
!  Initialize random number sequence so we can get identical results
!  everytime that we run the same solution.
!-----------------------------------------------------------------------
!
      sequence=759
      CALL ran_seed (sequence)

      RETURN
      END SUBROUTINE inp_par

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

      character (len=256), dimension(200), intent(inout) :: Cval

      real(r8), dimension(200), intent(inout) :: Rval
!
! Local variable declarations
!
      logical :: IsString, Kextract, decode, nested
      integer :: Iblank, Icomm, Icont, Ipipe, Kstr, Kend, Linp
      integer :: Lend, LenS, Lstr, Lval, Nmul, Schar
      integer :: copies, i, ic, ie, ierr, is, j, status

      integer, dimension(20) :: Imul

      integer :: decode_line

      character (len=1 ), parameter :: blank = ' '
      character (len=256) :: Vstring, line, string
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
        Icomm=INDEX(line_text,CHAR(33),BACK=.FALSE.)
        IF (Icomm.gt.0) Linp=Icomm-1
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
! [char(92)=\] or multi-line value [char(124)=|], if any.  Trim
! leading trailing blanks.
!
        IF (Kextract) THEN
          Icont=INDEX(line,CHAR(92 ),BACK=.FALSE.)
          Ipipe=INDEX(line,CHAR(124),BACK=.FALSE.)
          IF (Icont.gt.0) Lend=Icont-1
          IF (Ipipe.gt.0) Lend=Ipipe-1
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
!                   READ (string(1:LenS),*) Rval(Nval)
                    READ (string(1:LenS),*,IOSTAT=ierr) Rval(Nval)
                    if (ierr.ne.0) then
                      write(*,*)                                        &
     &     '#### ERROR :: Cannot interpret string ', string(1:LenS),    &
     &' as a REAL number.'
                    end if
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

      FUNCTION find_file (ng, fname, KeyWord) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This function checks if provided input file exits.                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     fname      File name (path and name).                            !
!     KeyWord    Keyword associated with file name (string,OPTIONAL).  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     foundit    The value of the result is TRUE/FALSE if the file     !
!                  was found or not.                                   !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      character (len=*), intent(in) :: fname
      character (len=*), intent(in) :: KeyWord
!
!  Local variable declarations.
!
      logical :: foundit, isURL

      integer :: ncid
!
      SourceFile='inp_par.F, find_file'
!
!-----------------------------------------------------------------------
!  Check if the file exit.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
!
!  Check for empty file name string.
!
      IF (LEN_TRIM(fname).eq.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(KeyWord)
 10       FORMAT (/,' INP_PAR - empty file name string for ',           &
     &            'input script KeyWord: ',a)
        END IF
        exit_flag=5
      END IF
!
!  Check if provided file is a URL.  This implies the file is a NetCDF
!  file on Data Access Protocol (DAP) server (like OPeNDAP).
!
      isURL=.FALSE.
      IF (INDEX(TRIM(fname),'http:').ne.0) THEN
        isURL=.TRUE.
      END IF
!
!  Use F90 intrinsic function for non URL files.
!
      IF (.not.isURL) THEN
        INQUIRE (FILE=TRIM(fname), EXIST=foundit)
!
!  Use NetCDF library (version 4.1.1 or higher) to check URL NetCDF
!  files.
!
      ELSE
        CALL netcdf_open (ng, iNLM, fname, 0, ncid)
        IF (exit_flag.eq.NoError) THEN
          foundit=.TRUE.
        END IF
      END IF

      RETURN
      END FUNCTION find_file

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

      FUNCTION load_lbc (Ninp, Vinp, line, nline, ifield, igrid,        &
     &                   iTrcStr, iTrcEnd, svname, S)
!
!=======================================================================
!                                                                      !
!  This function sets lateral boundary conditions logical switches     !
!  according to input string keywords.                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable (integer)                      !
!     Vinp       Input values (string)                                 !
!     line       Current input line (string)                           !
!     nline      Multi-line counter (integer)                          !
!     ifield     Lateral boundary variable index (integer)             !
!     igrid      Nested grid counter (integer)                         !
!     iTrcStr    Starting tracer index to process (integer)            !
!     iTrcEnd    Ending   tracer index to process (integer)            !
!     svname     State variable name (string)                          !
!     S          Derived type structure, TYPE(T_LBC)                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nline      Updated multi-line counter (integer)                  !
!     igrid      Updated nested grid counter (integer)                 !
!     S          Updated derived type structure, TYPE(T_LBC)           !
!     load_lbc   Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, ifield, iTrcStr, iTrcEnd
      integer, intent(inout) :: igrid, nline

      character (len=256), intent(in) :: line
      character (len=256), intent(in) :: Vinp(Ninp)
      character (len=*  ), intent(in) :: svname

      TYPE(T_LBC), intent(inout) :: S(4,nLBCvar,Ngrids)
!
!  Local variable declarations.
!
      integer :: Icont, i, ibry, ic
      integer :: load_lbc

      character (len=10) :: Bstring(4), string
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions switches in structure.
!-----------------------------------------------------------------------
!
!  Check current line for the continuation symbol [char(92)=\].
!
      Icont=INDEX(TRIM(line),CHAR(92) ,BACK=.FALSE.)
!
!  Extract lateral boundary condition keywords from Vinp. Notice that
!  additional array elements are added to Vinp during continuation
!  lines.
!
      i=nline*4
      Bstring(1)=TRIM(Vinp(i+1))
      Bstring(2)=TRIM(Vinp(i+2))
      Bstring(3)=TRIM(Vinp(i+3))
      Bstring(4)=TRIM(Vinp(i+4))
!
!  Advance or reset entry lines counter.
!
      IF (Icont.gt.0) THEN
        nline=nline+1
      ELSE
        nline=0
      END IF
!
!  Set switches for each boundary segment.
!
      ic=1
      IF ((0.lt.ifield).and.(ifield.le.nLBCvar)) THEN
        DO ibry=1,4
          string=uppercase(Bstring(ibry))
          SELECT CASE (TRIM(string))
            CASE ('CHA')
              S(ibry,ifield,igrid)%Chapman_implicit = .TRUE.
            CASE ('CHE')
              S(ibry,ifield,igrid)%Chapman_explicit = .TRUE.
            CASE ('CLA')
              S(ibry,ifield,igrid)%clamped = .TRUE.
              S(ibry,ifield,igrid)%acquire = .TRUE.
            CASE ('CLO')
              S(ibry,ifield,igrid)%closed = .TRUE.
            CASE ('FLA')
              S(ibry,ifield,igrid)%Flather = .TRUE.
              S(ibry,ifield,igrid)%acquire = .TRUE.
              S(ibry,isFsur,igrid)%acquire = .TRUE.
            CASE ('GRA')
              S(ibry,ifield,igrid)%gradient = .TRUE.
            CASE ('MIX')
              S(ibry,ifield,igrid)%mixed   = .TRUE.
              S(ibry,ifield,igrid)%acquire = .TRUE.
            CASE ('NES')
              S(ibry,ifield,igrid)%nested = .TRUE.
            CASE ('PER')
              S(ibry,ifield,igrid)%periodic = .TRUE.
              IF ((ibry.eq.ieast).or.(ibry.eq.iwest)) THEN
                EWperiodic(igrid)=.TRUE.
              ELSE IF ((ibry.eq.inorth).or.(ibry.eq.isouth)) THEN
                NSperiodic(igrid)=.TRUE.
              END IF
            CASE ('RAD')
              S(ibry,ifield,igrid)%radiation = .TRUE.
            CASE ('RADNUD')
              S(ibry,ifield,igrid)%radiation = .TRUE.
              S(ibry,ifield,igrid)%nudging = .TRUE.
              S(ibry,ifield,igrid)%acquire = .TRUE.
            CASE ('RED')
              S(ibry,ifield,igrid)%reduced = .TRUE.
#if defined FSOBC_REDUCED
              S(ibry,isFsur,igrid)%acquire = .TRUE.
#endif
            CASE ('SHC')
              S(ibry,ifield,igrid)%Shchepetkin = .TRUE.
              S(ibry,ifield,igrid)%acquire = .TRUE.
              S(ibry,isFsur,igrid)%acquire = .TRUE.
            CASE DEFAULT
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vinp(ibry)), TRIM(line)
              END IF
              exit_flag=2
              RETURN
          END SELECT
        END DO

#ifdef SOLVE3D
!
!  If processing tracers and last standard input entry (Icont=0), set
!  unspecified tracer values to the last tracer entry.
!
        IF ((iTrcStr.gt.0).and.(iTrcEnd.gt.0)) THEN
          IF ((Icont.eq.0).and.(ifield.lt.isTvar(iTrcEnd))) THEN
            DO i=ifield+1,isTvar(iTrcEnd)
              DO ibry=1,4
                S(ibry,i,igrid)%clamped   = S(ibry,ifield,igrid)%clamped
                S(ibry,i,igrid)%closed    = S(ibry,ifield,igrid)%closed
                S(ibry,i,igrid)%gradient  = S(ibry,ifield,igrid)%gradient
                S(ibry,i,igrid)%nested    = S(ibry,ifield,igrid)%nested
                S(ibry,i,igrid)%periodic  = S(ibry,ifield,igrid)%periodic
                S(ibry,i,igrid)%radiation = S(ibry,ifield,igrid)%radiation
                S(ibry,i,igrid)%nudging   = S(ibry,ifield,igrid)%nudging
                S(ibry,i,igrid)%acquire   = S(ibry,ifield,igrid)%acquire
              END DO
              ic=ic+1
            END DO
          END IF
        END IF
#endif
      END IF
!
!  If appropriate, increase or reset nested grid counter.
!
      IF ((Icont.gt.0).and.(Ngrids.gt.1)) THEN
        IF ((iTrcStr.gt.0).and.(iTrcEnd.gt.0)) THEN
          IF ((ifield.eq.isTvar(iTrcEnd)).or.(ic.gt.1)) THEN
            igrid=igrid+MIN(1,Icont)
          END IF
        ELSE
          igrid=igrid+MIN(1,Icont)
        END IF
        IF (igrid.gt.Ngrids) THEN
          IF (Master) THEN
            WRITE (stdout,20) TRIM(line)
          END IF
          exit_flag=2
          RETURN
        END IF
      ELSE
        igrid=1
      END IF
      load_lbc=ic

 10   FORMAT (/,' INP_PAR - illegal lateral boundary condition ',       &
     &        'keyword: ',a,/,11x,a)
 20   FORMAT (/,' INP_PAR - incorrect continuation symbol in line:',/,  &
     &        11x,a,/,11x,'number of nested grid values exceeded.')

      RETURN
      END FUNCTION load_lbc

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

      FUNCTION load_s1d (Nval, Fname, Fdim, line, label, igrid, Nfiles, &
     &                   S)
!
!=======================================================================
!                                                                      !
!  This function loads input values into requested 1D structure        !
!  containing information about I/O files.                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Nval       Number of values processed (integer)                  !
!     Fname      File name(s) processed (string array)                 !
!     Fdim       File name(s) dimension in calling program (integer)   !
!     line       Current input line (string)                           !
!     label      I/O structure label (string)                          !
!     igrid      Nested grid counter (integer)                         !
!     Nfiles     Number of files per grid (integer array)              !
!     S          Derived type structure, TYPE(T_IO)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     igrid      Updated nested grid counter.                          !
!     S          Updated derived type structure, TYPE(T_IO).           !
!     load_s1d   Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in)    :: Nval, Fdim
      integer, intent(inout) :: igrid
      integer, intent(inout) :: Nfiles(Ngrids)

      character (len=*),   intent(in) :: line
      character (len=256), intent(in) :: Fname(Fdim)
      character (len=*),   intent(inout) :: label

      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: load, persist

      integer :: Icont, Ipipe, i, j, lstr, my_Ngrids, ng
      integer :: load_s1d

      character (len=1 ), parameter :: blank = ' '
!
!-----------------------------------------------------------------------
!  Count files for all grids and activate load switch.
!-----------------------------------------------------------------------
!
!  Check current line for the continuation symbol [char(92)=\] or pipe
!  symbol [char(124)=|]. The continuation symbol is used to separate
!  string values for different grid, whereas the pipe symbol is used
!  to separate multi-string values for split input files. User may
!  split the records for a particular input field into several files.
!
      Icont=INDEX(TRIM(line),CHAR(92) ,BACK=.FALSE.)
      Ipipe=INDEX(TRIM(line),CHAR(124),BACK=.FALSE.)
      IF ((Icont.eq.0).and.(Ipipe.eq.0)) THEN
        load=.TRUE.                           ! last input string
      ELSE
        load=.FALSE.                          ! process next string
      END IF
!
!  Accumulate number of multi-files per each grid.
!
      Nfiles(igrid)=Nfiles(igrid)+1
!
!  Set grid counter.
!
      IF (.not.load) THEN
        igrid=igrid+MIN(1,Icont)
      END IF
      IF (igrid.gt.Ngrids) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(line)
        END IF
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Load I/O information into structure.
!-----------------------------------------------------------------------
!
      IF (load) THEN
!
!  If nesting and the number of file name entries is less than Ngrids,
!  persist the last values provided.  This is the case when not enough
!  entries are provided by "==" plural symbol after the KEYWORD.
!
        IF (igrid.lt.Ngrids) THEN
          DO i=igrid+1,Ngrids
            Nfiles(i)=Nfiles(igrid)
          END DO
          my_Ngrids=igrid
          persist=.TRUE.
        ELSE
          my_Ngrids=Ngrids
          persist=.FALSE.
        END IF
!
!  Allocate various fields in structure, if not continuation or pipe
!  symbol is found which indicates end of input data.
!
        DO ng=1,Ngrids
          allocate ( S(ng)%Nrec(Nfiles(ng)) )
          allocate ( S(ng)%time_min(Nfiles(ng)) )
          allocate ( S(ng)%time_max(Nfiles(ng)) )
          IF (label(1:3).eq.'FLT') THEN
#ifdef FLOAT_BIOLOGY
            allocate ( S(ng)%Vid(-10:NV) )
#else
            allocate ( S(ng)%Vid(-6:NV) )
#endif
          ELSE
            allocate ( S(ng)%Vid(NV) )
          END IF
          allocate ( S(ng)%Tid(MT) )
          allocate ( S(ng)%files(Nfiles(ng)) )
        END DO
!
!  Intialize strings to blank to facilitate processing.
!
        DO ng=1,Ngrids
          lstr=LEN(S(ng)%name)
          DO i=1,lstr
            S(ng)%base(i:i)=blank
            S(ng)%name(i:i)=blank
          END DO
          DO j=1,Nfiles(ng)
            DO i=1,lstr
              S(ng)%files(j)(i:i)=blank
            END DO
          END DO
        END DO
!
!  Initialize and load fields into structure.
!
        i=0
        DO ng=1,my_Ngrids
          S(ng)%Nfiles=Nfiles(ng)              ! number of multi-files
          S(ng)%Fcount=1                       ! multi-file counter
          S(ng)%Rindex=0                       ! time index
          S(ng)%ncid=-1                        ! closed NetCDF state
          S(ng)%Vid=-1                         ! NetCDF variables IDs
          S(ng)%Tid=-1                         ! NetCDF tracers IDs
          DO j=1,Nfiles(ng)
            i=i+1
            S(ng)%files(j)=TRIM(Fname(i))      ! load multi-files
            S(ng)%Nrec(j)=0                    ! record counter
            S(ng)%time_min(j)=0.0_r8           ! starting time
            S(ng)%time_max(j)=0.0_r8           ! ending time
          END DO
          S(ng)%label=TRIM(label)              ! structure label
          S(ng)%name=TRIM(S(ng)%files(1))      ! load first file
          lstr=LEN_TRIM(S(ng)%name)
          S(ng)%base=S(ng)%name(1:lstr-3)      ! do not include ".nc"
          Nfiles(ng)=0                         ! clean file counter
        END DO
!
!  If appropriate, persist last value(s).
!
        IF (persist) THEN
          DO ng=igrid+1,Ngrids
            S(ng)%Nfiles=S(igrid)%Nfiles
            S(ng)%Fcount=1
            S(ng)%Rindex=0
            S(ng)%ncid=-1
            S(ng)%Vid=-1
            S(ng)%Tid=-1
            DO j=1,S(igrid)%Nfiles
              S(ng)%files(j)=S(igrid)%files(j)
              S(ng)%Nrec(j)=0
              S(ng)%time_min(j)=0.0_r8
              S(ng)%time_max(j)=0.0_r8
            END DO
            S(ng)%label=TRIM(label)
            S(ng)%name=S(igrid)%name
            S(ng)%base=S(igrid)%base
            Nfiles(ng)=0
          END DO
        END IF
!
!  Reset counters and clean label.
!
        igrid=1
        DO ng=1,Ngrids
          Nfiles(ng)=0
        END DO
        DO i=1,LEN(label)
          label(i:i)=blank
        END DO
      END IF
      load_s1d=Nval

 10   FORMAT (/,' INP_PAR - incorrect continuation symbol in line:',/,  &
     &        11x,a,/,11x,'number of nested grid values exceeded.')

      RETURN
      END FUNCTION load_s1d

      FUNCTION load_s2d (Nval, Fname, Fdim, line, label, ifile, igrid,  &
     &                   Nfiles, Ncount, idim, S)
!
!=======================================================================
!                                                                      !
!  This function loads input values into requested 2D structure        !
!  containing information about input forcing files.                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Nval       Number of values processed (integer)                  !
!     Fname      File name(s) processed (string array)                 !
!     Fdim       File name(s) dimension in calling program (integer)   !
!     line       Current input line (string)                           !
!     label      I/O structure label (string)                          !
!     ifile      File structure counter (integer)                      !
!     igrid      Nested grid counter (integer)                         !
!     Nfiles     Number of input files per grid (integer vector)       !
!     Ncount     Number of files per grid counter (integer array)      !
!     idim       Size of structure inner dimension (integer)           !
!     S          Derived type structure, TYPE(T_IO)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ifile      Updated file counter.                                 !
!     igrid      Updated nested grid counter.                          !
!     S          Updated derived type structure, TYPE(T_IO).           !
!     load_s2d   Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in)    :: Nval, Fdim, idim
      integer, intent(in)    :: Nfiles(Ngrids)
      integer, intent(inout) :: ifile, igrid
      integer, intent(inout) :: Ncount(idim,Ngrids)

      character (len=*),   intent(in) :: line
      character (len=256), intent(in) :: Fname(Fdim)
      character (len=*),   intent(inout) :: label

      TYPE(T_IO), intent(inout) :: S(idim,Ngrids)
!
!  Local variable declarations.
!
      logical :: load, persist

      integer :: Icont, Ipipe, i, j, k, lstr, my_Ngrids, ng
      integer :: load_s2d

      character (len=1 ), parameter :: blank = ' '
!
!-----------------------------------------------------------------------
!  Count files for all grids and activate load switch.
!-----------------------------------------------------------------------
!
!  Check current line for the continuation symbol [char(92)=\] or pipe
!  symbol [char(124)=|]. The continuation symbol is used to separate
!  string values for different grid, whereas the pipe symbol is used
!  to separate multi-string values for split input files. User may
!  split the records for a particular input field into several files.
!
      Icont=INDEX(TRIM(line),CHAR(92) ,BACK=.FALSE.)
      Ipipe=INDEX(TRIM(line),CHAR(124),BACK=.FALSE.)
      IF ((Icont.eq.0).and.(Ipipe.eq.0)) THEN
        load=.TRUE.                           ! last input string
      ELSE
        load=.FALSE.                          ! process next string
      END IF
!
!  Accumulate number of multi-files per each grid.
!
      Ncount(ifile,igrid)=Ncount(ifile,igrid)+1
!
!  Set counters for next processing file, if any.  The continuation
!  symbol in the input "line" is used to advance the counters.
!
      IF (.not.load) THEN
        IF ((ifile.lt.Nfiles(igrid)).or.(Ipipe.ne.0)) THEN
          ifile=ifile+MIN(1,Icont)
        ELSE
          ifile=1
          igrid=igrid+MIN(1,Icont)
        END IF
      END IF
      IF (ifile.gt.idim) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(line)
        END IF
        exit_flag=2
        RETURN
      END IF
      IF (igrid.gt.Ngrids) THEN
        IF (Master) THEN
          WRITE (stdout,20) TRIM(line)
        END IF
        exit_flag=2
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Load I/O information into structure.
!-----------------------------------------------------------------------
!
      IF (load) THEN
!
!  If nesting and the number of file name entries is less than Ngrids,
!  persist the last values provided.  This is the case when not enough
!  entries are provided by "==" plural symbol after the KEYWORD.
!
        IF (igrid.lt.Ngrids) THEN
          DO j=igrid+1,Ngrids
            DO i=1,idim
              Ncount(i,j)=Ncount(i,igrid)
            END DO
          END DO
          my_Ngrids=igrid
          persist=.TRUE.
        ELSE
          my_Ngrids=Ngrids
          persist=.FALSE.
        END IF

!
!  Allocate various fields in structure, if not continuation or pipe
!  symbol is found which indicates end of input data.
!
        DO ng=1,Ngrids
          DO i=1,idim
            allocate ( S(i,ng)%Nrec(Ncount(i,ng)) )
            allocate ( S(i,ng)%time_min(Ncount(i,ng)) )
            allocate ( S(i,ng)%time_max(Ncount(i,ng)) )
            allocate ( S(i,ng)%Vid(NV) )
            allocate ( S(i,ng)%Tid(MT) )
            allocate ( S(i,ng)%files(Ncount(i,ng)) )
          END DO
        END DO
!
!  Intialize strings to blank to facilitate processing.
!
        DO ng=1,Ngrids
          DO i=1,idim
            lstr=LEN(S(i,ng)%name)
            DO j=1,lstr
              S(i,ng)%base(j:j)=blank
              S(i,ng)%name(j:j)=blank
            END DO
            DO k=1,Ncount(i,ng)
              DO j=1,lstr
                S(i,ng)%files(k)(j:j)=blank
              END DO
            END DO
          END DO
        END DO
!
!  Initialize and load fields into structure.
!
        k=0
        DO ng=1,my_Ngrids
          DO i=1,Nfiles(ng)
            S(i,ng)%Nfiles=Ncount(i,ng)         ! number of multi-files
            S(i,ng)%Fcount=1                    ! multi-file counter
            S(i,ng)%Rindex=0                    ! time index
            S(i,ng)%ncid=-1                     ! closed NetCDF state
            S(i,ng)%Vid=-1                      ! NetCDF variables IDs
            S(i,ng)%Tid=-1                      ! NetCDF tracers IDs
            DO j=1,Ncount(i,ng)
              k=k+1
              S(i,ng)%files(j)=TRIM(Fname(k))   ! load multi-files
              S(i,ng)%Nrec(j)=0                 ! record counter
              S(i,ng)%time_min(j)=0.0_r8        ! starting time
              S(i,ng)%time_max(j)=0.0_r8        ! ending time
            END DO
            S(i,ng)%label=TRIM(label)           ! structure label
            S(i,ng)%name=TRIM(S(i,ng)%files(1)) ! load first file
            lstr=LEN_TRIM(S(i,ng)%name)
            S(i,ng)%base=S(i,ng)%name(1:lstr-3) ! do not include ".nc"
          END DO
        END DO
!
!  If appropriate, persist last value(s).
!
        IF (persist) THEN
          DO ng=igrid+1,Ngrids
            DO i=1,Nfiles(ng)
              S(i,ng)%Nfiles=S(i,igrid)%Nfiles
              S(i,ng)%Fcount=1
              S(i,ng)%Rindex=0
              S(i,ng)%ncid=-1
              S(i,ng)%Vid=-1
              S(i,ng)%Tid=-1
              DO j=1,S(i,igrid)%Nfiles
                S(i,ng)%files(j)=S(i,igrid)%files(j)
                S(i,ng)%Nrec(j)=0
                S(i,ng)%time_min(j)=0.0_r8
                S(i,ng)%time_max(j)=0.0_r8
              END DO
              S(i,ng)%label=TRIM(label)
              S(i,ng)%name=S(i,igrid)%name
              S(i,ng)%base=S(i,igrid)%base
              Ncount(i,ng)=0
            END DO
          END DO
        END IF
!
!  Reset counters and clean label.
!
        igrid=1
        ifile=1
        DO ng=1,Ngrids
          DO i=1,idim
            Ncount(i,ng)=0
          END DO
        END DO
        DO i=1,LEN(label)
          label(i:i)=blank
        END DO
      END IF
      load_s2d=Nval

 10   FORMAT (/,' INP_PAR - incorrect continuation symbol in line:',/,  &
     &        11x,a,/,11x,'inner dimension of structure exceeded.')
 20   FORMAT (/,' INP_PAR - incorrect continuation symbol in line:',/,  &
     &        11x,a,/,11x,'number of nested grid values exceeded.')

      RETURN
      END FUNCTION load_s2d
