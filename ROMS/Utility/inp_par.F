#include "cppdefs.h"
      SUBROUTINE inp_par (model)
!
!svn $Id: inp_par.F 838 2008-11-17 04:22:18Z jcwarner $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasts
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
        CALL mp_bcasts (ng, model, date_str, LEN(date_str))
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
 20   WRITE (stdout,30)
      STOP
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
        CALL mp_bcasts (ng, model, Iname, LEN(Iname))
      END DO
      OPEN (inp, FILE=TRIM(Iname), FORM='formatted', STATUS='old',      &
     &      ERR=20)
      GO TO 40
 20   WRITE (stdout,30)
      STOP
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
          allocate ( BOUNDS(ng) % IstrTU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JstrTV(-1:Ntiles) )
          allocate ( BOUNDS(ng) % IendTU(-1:Ntiles) )
          allocate ( BOUNDS(ng) % JendTV(-1:Ntiles) )
        END DO
      END IF
!
!  Set bounds.
!
      Nghost=GHOST_POINTS
      DO ng=1,Ngrids
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
          WRITE (stdout,90) ng, HaloSizeI(ng), ng, HaloSizeJ(ng),        &
     &                      ng, TileSide(ng),  ng, TileSize(ng)
 90       FORMAT (/,' Maximum halo size in XI and ETA directions:',/,    &
     &            /,'               HaloSizeI(',i1,') = ',i7,            &
     &            /,'               HaloSizeJ(',i1,') = ',i7,            &
     &            /,'                TileSide(',i1,') = ',i7,            &
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
        visc4(ng)=SQRT(ABS(visc4(ng)))
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

#ifdef MODEL_COUPLING
      SUBROUTINE read_CouplePar (model)
!
!=======================================================================
!                                                                      !
!  This routine reads in physical model input parameters.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
!     USE mod_coupler
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
        CALL mp_bcasts (ng, model, Cname, LEN(Cname))
      END DO
      IF (MyRank.eq.0) THEN
        WRITE(stdout,*) 'Coupled Input File name = ', TRIM(Cname)
      END IF
      OPEN (inp, FILE=TRIM(Cname), FORM='formatted', STATUS='old',      &
     &      ERR=10)
      GO TO 30
 10   WRITE (stdout,20)
      WRITE(stdout,*) 'MyRank = ', MyRank, TRIM(Cname)
      STOP
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
  40  WRITE (out,60) line
      STOP
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
#ifdef MODEL_COUPLING
      USE mod_coupler
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
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
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
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
                WRITE (out,300) 'Lm', ng,                               &
     &                          'must be greater than zero.'
                STOP
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Mm') THEN
            Npts=load_i(Nval, Rval, Ngrids, Mm)
            DO ng=1,Ngrids
              IF (Mm(ng).le.0) THEN
                WRITE (out,300) 'Mm', ng,                               &
     &                          'must be greater than zero.'
                STOP
              END IF
            END DO
          ELSE IF (TRIM(KeyWord).eq.'N') THEN
            Npts=load_i(Nval, Rval, Ngrids, N)
            DO ng=1,Ngrids
              IF (N(ng).lt.0) THEN
                WRITE (out,300) 'N', ng,                                &
     &                          'must be greater than zero.'
                STOP
              END IF
            END DO
#if defined SEDIMENT && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Nbed') THEN
            Npts=load_i(Nval, Rval, 1, Nbed)
            IF (Nbed.le.0) THEN
              WRITE (out,290) 'Nbed = ', Nbed,                          &
     &                        'must be greater than zero.'
              STOP
            END IF
#endif
#ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'NAT') THEN
            Npts=load_i(Nval, Rval, 1, NAT)
            IF ((NAT.lt.0).or.(NAT.gt.2)) THEN
              WRITE (out,290) 'NAT = ', NAT,                            &
     &                        'make sure that NAT is either 1 or 2.'
              STOP
            END IF
# ifdef SALINITY
            IF (NAT.ne.2) THEN
              WRITE (out,290) 'NAT = ', NAT,                            &
     &                        'make sure that NAT is equal to 2.'
              STOP
            END IF
# endif
#endif
#if defined T_PASSIVE && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'NPT') THEN
            Npts=load_i(Nval, Rval, 1, NPT)
            IF (NPT.le.0) THEN
              WRITE (out,290) 'NPT = ', NPT,                            &
     &                        'must be greater than zero.'
              STOP
            END IF
#endif
#if defined SEDIMENT && defined SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'NCS') THEN
            Npts=load_i(Nval, Rval, 1, NCS)
            IF (NCS.lt.0) THEN
              WRITE (out,290) 'NCS = ', NCS,                            &
     &                        'must be greater than zero.'
              STOP
            END IF
            LreadNCS=.TRUE.
            IF (LreadNNS.and.((NCS+NNS).le.0)) THEN
              WRITE (out,290) 'NST = ', NCS+NNS,                        &
     &              'either NCS or NNS must be greater than zero.'
              STOP
            END IF
            NST=NST+NCS
          ELSE IF (TRIM(KeyWord).eq.'NNS') THEN
            Npts=load_i(Nval, Rval, 1, NNS)
            IF (NNS.lt.0) THEN
              WRITE (out,290) 'NNS = ',                                 &
     &                        'must be greater than zero.'
              STOP
            END IF
            LreadNNS=.TRUE.
            IF (LreadNCS.and.((NCS+NNS).le.0)) THEN
              WRITE (out,290) 'NST = ', NCS+NNS,                        &
     &              'either NCS or NNS must be greater than zero.'
              STOP
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
            CALL initialize_param
            CALL initialize_scalars
            CALL initialize_ncparam
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
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
              WRITE (out,260) 'NCV',                                    &
     &                        'Must be greater than or equal to 2*NEV+1'
              STOP
            END IF
# elif defined OPT_PERTURBATION
            IF (NCV.lt.(2*NEV)) THEN
              WRITE (out,260) 'NCV',                                    &
     &                        'Must be greater than or equal to 2*NEV'
              STOP
            END IF
# else
            IF (NCV.lt.(2*NEV)) THEN
              WRITE (out,260) 'NCV', 'Must be greater than NEV'
              STOP
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
                tnu2(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                tnu4(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VISC2') THEN
            Npts=load_r(Nval, Rval, Ngrids, visc2)
          ELSE IF (TRIM(KeyWord).eq.'VISC4') THEN
            Npts=load_r(Nval, Rval, Ngrids, visc4)
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, (NAT+NPT)*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT+NPT
                Akt_bak(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKV_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
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
          ELSE IF (TRIM(KeyWord).eq.'THETA_S') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_s)
          ELSE IF (TRIM(KeyWord).eq.'THETA_B') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_b)
          ELSE IF (TRIM(KeyWord).eq.'TCLINE') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcline)
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
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
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
              WRITE (out,280) 'isFsur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isFsur)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isUbar)') THEN
            IF (isUbar.eq.0) THEN
              WRITE (out,280) 'isUbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isUbar)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isVbar)') THEN
            IF (isVbar.eq.0) THEN
              WRITE (out,280) 'isVbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isVbar)=Ladsen(ng)
            END DO
# ifdef SOLVE3D  
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isUvel)') THEN
            IF (isUvel.eq.0) THEN
              WRITE (out,280) 'isUvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isUvel)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isVvel)') THEN
            IF (isVvel.eq.0) THEN
              WRITE (out,280) 'isVvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Ladsen)
            DO ng=1,Ngrids
              SCALARS(ng)%Lstate(isVvel)=Ladsen(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Lstate(isTvar)') THEN
            IF (MAXVAL(isTvar).eq.0) THEN
              WRITE (out,280) 'isTvar'
              STOP
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
              WRITE (out,280) 'isUstr'
              STOP
            END IF
            Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isUstr,1))
          ELSE IF (TRIM(KeyWord).eq.'SO_sdev(isVstr)') THEN
            IF (isUstr.eq.0) THEN
              WRITE (out,280) 'isVstr'
              STOP
            END IF
            Npts=load_r(Nval, Rval, Ngrids, SO_sdev(isVstr,1))
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'SO_sdev(isTsur)') THEN
            IF (MAXVAL(isTsur).eq.0) THEN
              WRITE (out,280) 'isTsur'
              STOP
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
#if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(KeyWord).eq.'Hout(idBath)') THEN
            IF (idbath.eq.0) THEN
              WRITE (out,280) 'idbath'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idBath,1:Ngrids))
#endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idFsur)') THEN
            IF (idFsur.eq.0) THEN
              WRITE (out,280) 'idFsur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idFsur,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbar)') THEN
            IF (idUbar.eq.0) THEN
              WRITE (out,280) 'idUbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbar,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbar)') THEN
            IF (idVbar.eq.0) THEN
              WRITE (out,280) 'idVbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbar,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUvel)') THEN
            IF (idUvel.eq.0) THEN
              WRITE (out,280) 'idUvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvel)') THEN
            IF (idVvel.eq.0) THEN
              WRITE (out,280) 'idVvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWvel)') THEN
            IF (idWvel.eq.0) THEN
              WRITE (out,280) 'idWvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idOvel)') THEN
            IF (idOvel.eq.0) THEN
              WRITE (out,280) 'idOvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idOvel,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            IF (MAXVAL(idTvar).eq.0) THEN
              WRITE (out,280) 'idTvar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                i=idTvar(itrc)
                Hout(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUsms)') THEN
            IF (idUsms.eq.0) THEN
              WRITE (out,280) 'idUsms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUsms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVsms)') THEN
            IF (idVsms.eq.0) THEN
              WRITE (out,280) 'idVsms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVsms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbms)') THEN
            IF (idUbms.eq.0) THEN
              WRITE (out,280) 'idUbms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbms)') THEN
            IF (idVbms.eq.0) THEN
              WRITE (out,280) 'idVbms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbms,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbrs)') THEN
            IF (idUbrs.eq.0) THEN
              WRITE (out,280) 'idUbrs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbrs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbrs)') THEN
            IF (idVbrs.eq.0) THEN
              WRITE (out,280) 'idVbrs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbrs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbws)') THEN
            IF (idUbws.eq.0) THEN
              WRITE (out,280) 'idUbws'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbws,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbws)') THEN
            IF (idVbws.eq.0) THEN
              WRITE (out,280) 'idVbws'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbws,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbcs)') THEN
            IF (idUbcs.eq.0) THEN
              WRITE (out,280) 'idUbcs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbcs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbcs)') THEN
            IF (idVbcs.eq.0) THEN
              WRITE (out,280) 'idVbcs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbcs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbot)') THEN
            IF (idUbot.eq.0) THEN
              WRITE (out,280) 'idUbot'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbot,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbot)') THEN
            IF (idVbot.eq.0) THEN
              WRITE (out,280) 'idVbot'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbot,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbur)') THEN
            IF (idUbur.eq.0) THEN
              WRITE (out,280) 'idUbur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbur,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbvr)') THEN
            IF (idVbvr.eq.0) THEN
              WRITE (out,280) 'idVbvr'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbvr,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2xx)') THEN
            IF (idW2xx.eq.0) THEN
              WRITE (out,280) 'idW2xx'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xx,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2xy)') THEN
            IF (idW2xy.eq.0) THEN
              WRITE (out,280) 'idW2xy'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2xy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW2yy)') THEN
            IF (idW2yy.eq.0) THEN
              WRITE (out,280) 'idW2yy'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW2yy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3xx)') THEN
            IF (idW3xx.eq.0) THEN
              WRITE (out,280) 'idW3xx'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xx,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3xy)') THEN
            IF (idW3xy.eq.0) THEN
              WRITE (out,280) 'idW3xy'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3xy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3yy)') THEN
            IF (idW3yy.eq.0) THEN
              WRITE (out,280) 'idW3yy'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3yy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3zx)') THEN
            IF (idW3zx.eq.0) THEN
              WRITE (out,280) 'idW3zx'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zx,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idW3zy)') THEN
            IF (idW3zy.eq.0) THEN
              WRITE (out,280) 'idW3zy'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idW3zy,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU2rs)') THEN
            IF (idU2rs.eq.0) THEN
              WRITE (out,280) 'idU2rs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU2rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV2rs)') THEN
            IF (idV2rs.eq.0) THEN
              WRITE (out,280) 'idV2rs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV2rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU2Sd)') THEN
            IF (idU2Sd.eq.0) THEN
              WRITE (out,280) 'idU2Sd'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU2Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV2Sd)') THEN
            IF (idV2Sd.eq.0) THEN
              WRITE (out,280) 'idV2Sd'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV2Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU3rs)') THEN
            IF (idU3rs.eq.0) THEN
              WRITE (out,280) 'idU3rs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU3rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV3rs)') THEN
            IF (idV3rs.eq.0) THEN
              WRITE (out,280) 'idV3rs'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV3rs,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idU3Sd)') THEN
            IF (idU3Sd.eq.0) THEN
              WRITE (out,280) 'idU3Sd'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idU3Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idV3Sd)') THEN
            IF (idV3Sd.eq.0) THEN
              WRITE (out,280) 'idV3Sd'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idV3Sd,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWamp)') THEN
            IF (idWamp.eq.0) THEN
              WRITE (out,280) 'idWamp'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWamp,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWlen)') THEN
            IF (idWlen.eq.0) THEN
              WRITE (out,280) 'idWlen'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWlen,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWdir)') THEN
            IF (idWdir.eq.0) THEN
              WRITE (out,280) 'idWdir'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWdir,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            IF (idTsur(itemp).eq.0) THEN
              WRITE (out,280) 'idTsur'
              STOP
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
              WRITE (out,280) 'idLhea'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLhea,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idShea)') THEN
            IF (idShea.eq.0) THEN
              WRITE (out,280) 'idShea'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idShea,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLrad)') THEN
            IF (idLrad.eq.0) THEN
              WRITE (out,280) 'idLrad'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLrad,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSrad)') THEN
            IF (idSrad.eq.0) THEN
              WRITE (out,280) 'idSrad'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idevap)') THEN
            IF (idevap.eq.0) THEN
              WRITE (out,280) 'idevap'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idevap,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idrain)') THEN
            IF (idrain.eq.0) THEN
              WRITE (out,280) 'idrain'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idrain,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idDano)') THEN
            IF (idDano.eq.0) THEN
              WRITE (out,280) 'idDano'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idDano,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvis)') THEN
            IF (idVvis.eq.0) THEN
              WRITE (out,280) 'idVvis'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvis,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTdif)') THEN
            IF (idTdif.eq.0) THEN
              WRITE (out,280) 'idTdif'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTdif,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSdif)') THEN
            IF (idSdif.eq.0) THEN
              WRITE (out,280) 'idSdif'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSdif,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHsbl)') THEN
            IF (idHsbl.eq.0) THEN
              WRITE (out,280) 'idHsbl'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHsbl,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHbbl)') THEN
            IF (idHbbl.eq.0) THEN
              WRITE (out,280) 'idHbbl'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHbbl,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtke)') THEN
            IF (idMtke.eq.0) THEN
              WRITE (out,280) 'idMtke'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtke,1:Ngrids))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtls)') THEN
            IF (idMtls.eq.0) THEN
              WRITE (out,280) 'idMtls'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtls,1:Ngrids))
#if defined SOLVE3D && defined T_PASSIVE
          ELSE IF (TRIM(KeyWord).eq.'Hout(inert)') THEN
            IF (MAXVAL(inert).eq.0) THEN
              WRITE (out,280) 'inert'
              STOP
            END IF
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
              STOP
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
                WRITE (out,260) 'NFFILES',                              &
     &                          'Must be equal or greater than one.'
                STOP
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
                IF (Nval.eq.(i+(Ngrids-1)*nFfiles(ng))) THEN
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
  10  WRITE (out,50) line
      STOP
  20  CLOSE (inp)
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
#if defined IS4DVAR || defined IOM || defined WEAK_CONSTRAINT
!
!  If weak constraint, disallow recycling of the adjoint model.
!
        LcycleADJ(ng)=.FALSE.
#endif
#if defined IS4DVAR || defined IS4DVAR_OLD || defined S4DVAR
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
          WRITE (out,120) ERstr, 'ERstr',                               &
     &          'Starting ensemble/perturbation run number.'
          WRITE (out,120) ERend, 'ERend',                               &
     &          'Ending ensemble/perturbation run number.'
#ifdef FOUR_DVAR
          WRITE (out,120) Nouter, 'Nouter',                             &
     &          'Maximun number of 4DVAR outer loop iterations.'
#endif
#if defined IS4DVAR         || defined IS4DVAR_OLD || \
    defined OBS_SENSITIVITY || defined W4DPSAS     || \
    defined W4DVAR
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
# if defined WEAK_CONSTRAINT || defined IOM
          LcycleTLM(ng)=.FALSE.       ! disallow recycling
          nTLM(ng)=nHIS(ng)           ! same as basic state
# endif
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
            WRITE (out,190) tnu2(itrc,ng), 'tnu2', itrc,                &
     &            'Horizontal, harmonic mixing coefficient (m2/s)',     &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
# ifdef TS_DIF4
          DO itrc=1,NAT+NPT
            WRITE (out,190) tnu4(itrc,ng), 'tnu4', itrc,                &
     &            'Horizontal, biharmonic mixing coefficient (m4/s)',   &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#ifdef UV_VIS2
          WRITE (out,210) visc2(ng), 'visc2',                           &
     &          'Horizontal, harmonic mixing coefficient (m2/s)',       &
     &          'for momentum.'
#endif
#ifdef UV_VIS4
          WRITE (out,210) visc4(ng), 'visc4',                           &
     &          'Horizontal, biharmonic mixing coefficient (m4/s)',     &
     &          'for momentum.'
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
            STOP
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
            STOP
          END IF
          WRITE (out,200) blk_ZT(ng), 'blk_ZT',                         &
     &          'Height (m) of surface air temperature measurement.'
          IF (blk_ZT(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZT.',                                  &
     &            'It must be greater than zero.'
            STOP
          END IF
          WRITE (out,200) blk_ZW(ng), 'blk_ZW',                         &
     &          'Height (m) of surface winds measurement.'
          IF (blk_ZW(ng).le.0.0_r8) THEN
            WRITE (out,260) 'blk_ZW.',                                  &
     &            'It must be greater than zero.'
            STOP
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
            STOP
          END IF
# endif
# ifdef BODYFORCE
          WRITE (out,130) levsfrc(ng), 'levsfrc',                       &
     &          'Deepest level to apply surface stress as a',           &
     &          'bodyforce.'
          IF ((levsfrc(ng).lt.1).or.(levsfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levsfrc.',                                 &
     &            'Out of range surface bodyforce level.'
            STOP
          END IF
          WRITE (out,130) levbfrc(ng), 'levbfrc',                       &
     &          'Shallowest level to apply bottom stress as a',         &
     &          'bodyforce.'
          IF ((levbfrc(ng).lt.1).or.(levbfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levbfrc.',                                 &
     &            'Out of range bottom bodyforce level.'
            STOP
          END IF
# endif
#endif
#ifdef SOLVE3D
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
            STOP
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
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
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
            STOP
          END IF
          WRITE (out,120) KendS(ng), 'KendS',                           &
     &          'Shallowest level whose ADM sensitivity is required.'
          IF ((KendS(ng).lt.1).or.(KendS(ng).gt.N(ng))) THEN
            WRITE (out,260) 'KendS',                                    &
     &            'Out of range ADM sensitivity level.'
            STOP
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
#if defined NEARSHORE_MELLOR
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
          IF (Hout(idU2Sd,ng)) WRITE (out,170) Hout(idU2Sd,ng),         &
     &       'Hout(idU2Sd)',                                            &
     &       'Write out 2D u-momentum stokes velocity.'
          IF (Hout(idV2Sd,ng)) WRITE (out,170) Hout(idV2Sd,ng),         &
     &       'Hout(idV2Sd)',                                            &
     &       'Write out 2D v-momentum stokes velocity.'
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
          IF (Hout(idU3Sd,ng)) WRITE (out,170) Hout(idU3Sd,ng),         &
     &       'Hout(idU3Sd)',                                            &
     &       'Write out 3D u-momentum stokes velocity.'
          IF (Hout(idV3Sd,ng)) WRITE (out,170) Hout(idV3Sd,ng),         &
     &       'Hout(idV3Sd)',                                            &
     &       'Write out 3D v-momentum stokes velocity.'
# endif
#endif
          IF (Hout(idWamp,ng)) WRITE (out,170) Hout(idWamp,ng),         &
     &       'Hout(idWamp)',                                            &
     &       'Write out wave height.'
          IF (Hout(idWlen,ng)) WRITE (out,170) Hout(idWlen,ng),         &
     &       'Hout(idWlen)',                                            &
     &       'Write out wave length.'
          IF (Hout(idWdir,ng)) WRITE (out,170) Hout(idWdir,ng),         &
     &       'Hout(idWdir)',                                            &
     &       'Write out wave direction.'
#if defined SOLVE3D && defined T_PASSIVE
          DO itrc=1,NPT
            IF (Hout(idTvar(inert(itrc)),ng)) WRITE (out,180)           &
     &          Hout(idTvar(inert(itrc)),ng), 'Hout(inert)',            &
     &          'Write out inert passive tracer ', itrc,                &
     &          TRIM(Vname(1,idTvar(inert(itrc))))
          END DO
#endif
#if defined SEDIMENT && defined BEDLOAD
          DO itrc=1,NST
            IF (Hout(idUbld(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idUbld(itrc),ng), 'Hout(idUbld)',                  &
     &          'Write out bed load u-direction ', itrc,                &
     &          TRIM(Vname(1,idUbld(itrc)))
            IF (Hout(idVbld(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idVbld(itrc),ng), 'Hout(idVbld)',                  &
     &          'Write out bed load v-direction ', itrc,                &
     &          TRIM(Vname(1,idVbld(itrc)))
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
#if defined WEAK_CONSTRAINT || defined IOM
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
    !(defined FOUR_DVAR       || defined GRADIENT_CHECK   || \
      defined OBS_SENSITIVITY || defined OPT_OBSERVATIONS || \
      defined SANITY_CHECK    || defined TLM_CHECK)
          fname=ITLname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) '    Input Tangent Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if (defined WEAK_CONSTRAINT || defined IOM) && !defined W4DPSAS
          fname=IRPname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,230) 'Input Representer Initial File:  ',          &
     &                    TRIM(fname)
# endif
# if defined ADJOINT && \
    !(defined AD_SENSITIVITY  || defined FOUR_DVAR        || \
      defined OBS_SENSITIVITY || defined OPT_OBSERVATIONS || \
      defined SANITY_CHECK    || defined SO_SEMI)
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
#if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SO_SEMI
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
  30      WRITE (out,270) TRIM(fname)
          STOP
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
          tnu4(itrc,ng)=SQRT(ABS(tnu4(itrc,ng)))
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

      RETURN
      END SUBROUTINE read_PhyPar

#ifdef BIOLOGY
      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
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
      integer :: Npts, Nval, i, is, itrc, ng, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(Ngrids) :: Lbio
      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio
 
      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval

# if defined BIO_FASHAM
!
!-----------------------------------------------------------------------
!  Read in Fasham biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN 
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'AttSW') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSW)
          ELSE IF (TRIM(KeyWord).eq.'AttChl') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttChl)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'Vp0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vp0)
          ELSE IF (TRIM(KeyWord).eq.'I_thNH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, I_thNH4)
          ELSE IF (TRIM(KeyWord).eq.'D_p5NH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, D_p5NH4)
          ELSE IF (TRIM(KeyWord).eq.'NitriR') THEN
            Npts=load_r(Nval, Rval, Ngrids, NitriR)
          ELSE IF (TRIM(KeyWord).eq.'K_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NO3)
          ELSE IF (TRIM(KeyWord).eq.'K_NH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NH4)
          ELSE IF (TRIM(KeyWord).eq.'K_Phy') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_Phy)
          ELSE IF (TRIM(KeyWord).eq.'Chl2C_m') THEN
            Npts=load_r(Nval, Rval, Ngrids, Chl2C_m)
          ELSE IF (TRIM(KeyWord).eq.'ChlMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, ChlMin)
          ELSE IF (TRIM(KeyWord).eq.'PhyCN') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyCN)
          ELSE IF (TRIM(KeyWord).eq.'PhyIP') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIP)
          ELSE IF (TRIM(KeyWord).eq.'PhyIS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIS)
          ELSE IF (TRIM(KeyWord).eq.'PhyMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMin)
          ELSE IF (TRIM(KeyWord).eq.'PhyMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMR)
          ELSE IF (TRIM(KeyWord).eq.'ZooAE_N') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooAE_N)
          ELSE IF (TRIM(KeyWord).eq.'ZooBM') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooBM)
          ELSE IF (TRIM(KeyWord).eq.'ZooCN') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooCN)
          ELSE IF (TRIM(KeyWord).eq.'ZooER') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooER)
          ELSE IF (TRIM(KeyWord).eq.'ZooGR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGR)
          ELSE IF (TRIM(KeyWord).eq.'ZooMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMin)
          ELSE IF (TRIM(KeyWord).eq.'ZooMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMR)
          ELSE IF (TRIM(KeyWord).eq.'LDeRRN') THEN
            Npts=load_r(Nval, Rval, Ngrids, LDeRRN)
          ELSE IF (TRIM(KeyWord).eq.'LDeRRC') THEN
            Npts=load_r(Nval, Rval, Ngrids, LDeRRC)
          ELSE IF (TRIM(KeyWord).eq.'CoagR') THEN
            Npts=load_r(Nval, Rval, Ngrids, CoagR)
          ELSE IF (TRIM(KeyWord).eq.'SDeRRN') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeRRN)
          ELSE IF (TRIM(KeyWord).eq.'SDeRRC') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeRRC)
          ELSE IF (TRIM(KeyWord).eq.'wPhy') THEN
            Npts=load_r(Nval, Rval, Ngrids, wPhy)
          ELSE IF (TRIM(KeyWord).eq.'wLDet') THEN
            Npts=load_r(Nval, Rval, Ngrids, wLDet)
          ELSE IF (TRIM(KeyWord).eq.'wSDet') THEN
            Npts=load_r(Nval, Rval, Ngrids, wSDet)
          ELSE IF (TRIM(KeyWord).eq.'pCO2air') THEN
            Npts=load_r(Nval, Rval, Ngrids, pCO2air)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTvar(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTsur(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#  ifdef DIAGNOSTICS_BIO
#   ifdef CARBON
          ELSE IF (TRIM(KeyWord).eq.'Hout(iCOfx)') THEN
            IF (iDbio2(iCOfx).eq.0) THEN
              WRITE (out,130) 'iDbio2(iCOfx)'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(iCOfx)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
#   endif
#   ifdef DENITRIFICATION
          ELSE IF (TRIM(KeyWord).eq.'Hout(iDNIT)') THEN
            IF (iDbio2(iDNIT).eq.0) THEN
              WRITE (out,130) 'iDbio2(iDNIT)'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(iDNIT)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
#   endif
#   ifdef CARBON
          ELSE IF (TRIM(KeyWord).eq.'Hout(ipCO2)') THEN
            IF (iDbio2(ipCO2).eq.0) THEN
              WRITE (out,130) 'iDbio2(ipCO2)'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(ipCO2)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
#   endif
#   ifdef OXYGEN
          ELSE IF (TRIM(KeyWord).eq.'Hout(iO2fx)') THEN
            IF (iDbio2(iO2fx).eq.0) THEN
              WRITE (out,130) 'iDbio2(iO2fx)'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio2(iO2fx)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
#   endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(iPPro)') THEN
            IF (iDbio3(iPPro).eq.0) THEN
              WRITE (out,130) 'iDbio3(iPPro)'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio3(iPPro)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iNO3u)') THEN
            IF (iDbio3(iNO3u).eq.0) THEN
              WRITE (out,130) 'iDbio3(iNO3u)'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Lbio)
            i=iDbio3(iNO3u)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbio(ng)
            END DO
#  endif
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,100) AttSW(ng), 'AttSW',                         &
     &            'Light attenuation of seawater (m-1).'
            WRITE (out,100) AttChl(ng), 'AttChl',                       &
     &            'Light attenuation by chlorophyll (1/(mg_Chl m-2)).'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) Vp0(ng), 'Vp0',                             &
     &            'Eppley temperature-limited growth parameter',        &
     &            '(nondimensional).'
            WRITE (out,100) I_thNH4(ng), 'I_thNH4',                     &
     &            'Radiation threshold for nitrification (W/m2).'
            WRITE (out,100) D_p5NH4(ng), 'D_p5NH4',                     &
     &            'Half-saturation radiation for nitrification (W/m2).'
            WRITE (out,100) NitriR(ng), 'NitriR',                       &
     &            'Nitrification rate (day-1).'
            WRITE (out,110) K_NO3(ng), 'K_NO3',                         &
     &            'Inverse half-saturation for phytoplankton NO3',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,110) K_NH4(ng), 'K_NH4',                         &
     &            'Inverse half-saturation for phytoplankton NH4',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,110) K_Phy(ng), 'K_Phy',                         &
     &            'Zooplankton half-saturation constant for ingestion', &
     &            '(mmol_N m-3)^2.'
            WRITE (out,100) Chl2C_m(ng), 'Chl2C_m',                     &
     &            'Maximum chlorophyll to carbon ratio (mg_Chl/mg_C).'
            WRITE (out,100) ChlMin(ng), 'ChlMin',                       &
     &            'Chlorophyll minimum threshold (mg_Chl/m3).'
            WRITE (out,100) PhyCN(ng), 'PhyCN',                         &
     &            'Phytoplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,100) PhyIP(ng), 'PhyIP',                         &
     &            'Phytoplankton NH4 inhibition parameter (1/mmol_N).'
            WRITE (out,110) PhyIS(ng), 'PhyIS',                         &
     &            'Phytoplankton growth, initial slope of P-I curve',   &
     &            '(mg_C/(mg_Chl Watts m-2 day)).'
            WRITE (out,100) PhyMin(ng), 'PhyMin',                       &
     &            'Phytoplankton minimum threshold (mmol_N/m3).'
            WRITE (out,100) PhyMR(ng), 'PhyMR',                         &
     &            'Phytoplankton mortality rate (day-1).'
            WRITE (out,110) ZooAE_N(ng), 'ZooAE_N',                     &
     &            'Zooplankton nitrogen assimilation efficiency',       &
     &            '(nondimensional).'
            WRITE (out,100) ZooBM(ng), 'ZooBM',                         &
     &            'Rate for zooplankton basal metabolism (1/day).'
            WRITE (out,100) ZooCN(ng), 'ZooCN',                         &
     &            'Zooplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,100) ZooER(ng), 'ZooER',                         &
     &            'Zooplankton specific excretion rate (day-1).'
            WRITE (out,100) ZooGR(ng), 'ZooGR',                         &
     &            'Zooplankton maximum growth rate (day-1).'
            WRITE (out,100) ZooMin(ng), 'ZooMin',                       &
     &            'Zooplankton minimum threshold (mmol_N/m3).'
            WRITE (out,100) ZooER(ng), 'ZooMR',                         &
     &            'Zooplankton mortality rate (day-1).'
            WRITE (out,100) LDeRRN(ng), 'LDeRRN',                       &
     &            'Large detritus N re-mineralization rate (day-1).'
            WRITE (out,100) LDeRRC(ng), 'LDeRRC',                       &
     &            'Large detritus C re-mineralization rate (day-1).'
            WRITE (out,100) CoagR(ng), 'CoagR',                         &
     &            'Coagulation rate (day-1).'
            WRITE (out,100) SDeRRN(ng), 'SDeRRN',                       &
     &            'Remineralization rate for small detritus N (day-1).'
            WRITE (out,100) SDeRRC(ng), 'SDeRRC',                       &
     &            'Remineralization rate for small detritus C (day-1).'
            WRITE (out,100) wPhy(ng), 'wPhy',                           &
     &            'Phytoplankton sinking velocity (m/day).'
            WRITE (out,100) wLDet(ng), 'wLDet',                         &
     &            'Large detritus sinking velocity (m/day).'
            WRITE (out,100) wSDet(ng), 'wSDet',                         &
     &            'Small detritus sinking velocity (m/day).'
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  ifdef DIAGNOSTICS_BIO
#   if !(defined CARBON || defined OXYGEN || defined DENITRIFICATION) 
            DO itrc=1,NDbio2d
              i=iDbio2(itrc)
              IF (Hout(i,ng)) WRITE (out,140)                           &
     &            Hout(i,ng), 'Hout(iDbio2)', 'Diagnostics for',        &
     &            TRIM(Vname(1,i))
            END DO
#   endif
            DO itrc=1,NDbio3d
              i=iDbio3(itrc)
              IF (Hout(i,ng)) WRITE (out,140)                           &
     &            Hout(i,ng), 'Hout(iDbio3)', 'Diagnostics for',        &
     &            TRIM(Vname(1,i))
            END DO
#  endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Fasham Model Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t30,a)
  80  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',       &
     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t30,a,1x,a)

# elif defined BIO_LIMADONEY
!
!-----------------------------------------------------------------------
!  Read in biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'AttSW') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSW)
          ELSE IF (TRIM(KeyWord).eq.'AttChl') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttChl)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'alphaPI') THEN
            Npts=load_r(Nval, Rval, Ngrids, alphaPI)
          ELSE IF (TRIM(KeyWord).eq.'thetaN0') THEN
            Npts=load_r(Nval, Rval, Ngrids, thetaN0)
          ELSE IF (TRIM(KeyWord).eq.'tref') THEN
            Npts=load_r(Nval, Rval, Ngrids, tref)
          ELSE IF (TRIM(KeyWord).eq.'e_a') THEN
            Npts=load_r(Nval, Rval, Ngrids, e_a)
          ELSE IF (TRIM(KeyWord).eq.'p_cref') THEN
            Npts=load_r(Nval, Rval, Ngrids, p_cref)
          ELSE IF (TRIM(KeyWord).eq.'z_umax1') THEN
            Npts=load_r(Nval, Rval, Ngrids, z_umax1)
          ELSE IF (TRIM(KeyWord).eq.'z_umax2') THEN
            Npts=load_r(Nval, Rval, Ngrids, z_umax2)
          ELSE IF (TRIM(KeyWord).eq.'p1_excr') THEN
            Npts=load_r(Nval, Rval, Ngrids, p1_excr)
          ELSE IF (TRIM(KeyWord).eq.'p2_excr') THEN
            Npts=load_r(Nval, Rval, Ngrids, p2_excr)
          ELSE IF (TRIM(KeyWord).eq.'p1_mort2') THEN
            Npts=load_r(Nval, Rval, Ngrids, p1_mort2)
          ELSE IF (TRIM(KeyWord).eq.'p2_mort2') THEN
            Npts=load_r(Nval, Rval, Ngrids, p2_mort2)
          ELSE IF (TRIM(KeyWord).eq.'z_mort') THEN
            Npts=load_r(Nval, Rval, Ngrids, z_mort)
          ELSE IF (TRIM(KeyWord).eq.'z_mort2') THEN
            Npts=load_r(Nval, Rval, Ngrids, z_mort2)
          ELSE IF (TRIM(KeyWord).eq.'d1_rem') THEN
            Npts=load_r(Nval, Rval, Ngrids, d1_rem)
          ELSE IF (TRIM(KeyWord).eq.'d2_rem') THEN
            Npts=load_r(Nval, Rval, Ngrids, d2_rem)
          ELSE IF (TRIM(KeyWord).eq.'wDet2') THEN
            Npts=load_r(Nval, Rval, Ngrids, wDet2)
          ELSE IF (TRIM(KeyWord).eq.'p1_kno3') THEN
            Npts=load_r(Nval, Rval, Ngrids, p1_kno3)
          ELSE IF (TRIM(KeyWord).eq.'p2_kno3') THEN
            Npts=load_r(Nval, Rval, Ngrids, p2_kno3)
          ELSE IF (TRIM(KeyWord).eq.'p1_knh4') THEN
            Npts=load_r(Nval, Rval, Ngrids, p1_knh4)
          ELSE IF (TRIM(KeyWord).eq.'p2_knh4') THEN
            Npts=load_r(Nval, Rval, Ngrids, p2_knh4)
          ELSE IF (TRIM(KeyWord).eq.'p2_si') THEN
            Npts=load_r(Nval, Rval, Ngrids, p2_si)
          ELSE IF (TRIM(KeyWord).eq.'p_nitr') THEN
            Npts=load_r(Nval, Rval, Ngrids, p_nitr)
          ELSE IF (TRIM(KeyWord).eq.'falloc0') THEN
            Npts=load_r(Nval, Rval, Ngrids, falloc0)
          ELSE IF (TRIM(KeyWord).eq.'si_alloc') THEN
            Npts=load_r(Nval, Rval, Ngrids, si_alloc)
          ELSE IF (TRIM(KeyWord).eq.'z_grz') THEN
            Npts=load_r(Nval, Rval, Ngrids,  z_grz)
          ELSE IF (TRIM(KeyWord).eq.'zegest0') THEN
            Npts=load_r(Nval, Rval, Ngrids, zegest0)
          ELSE IF (TRIM(KeyWord).eq.'r_ref') THEN
            Npts=load_r(Nval, Rval, Ngrids, r_ref)
          ELSE IF (TRIM(KeyWord).eq.'q_min') THEN
            Npts=load_r(Nval, Rval, Ngrids, q_min)
          ELSE IF (TRIM(KeyWord).eq.'q_max') THEN
            Npts=load_r(Nval, Rval, Ngrids, q_max)
          ELSE IF (TRIM(KeyWord).eq.'q_simin') THEN
            Npts=load_r(Nval, Rval, Ngrids, q_simin)
          ELSE IF (TRIM(KeyWord).eq.'q_simax') THEN
            Npts=load_r(Nval, Rval, Ngrids, q_simax)
          ELSE IF (TRIM(KeyWord).eq.'vncref') THEN
            Npts=load_r(Nval, Rval, Ngrids, vncref)
          ELSE IF (TRIM(KeyWord).eq.'biolambda') THEN
            Npts=load_r(Nval, Rval, Ngrids, biolambda)          
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTvar(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTsur(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,100) AttSW(ng), 'AttSW',                         &
     &            'Light attenuation of seawater (m-1).'
            WRITE (out,100) AttChl(ng), 'AttChl',                       &
     &            'Light attenuation by chlorophyll (1/(mg_Chl m-2)).'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) alphaPI(ng), 'alphaPI',                     &
     &            '',            &
     &            ''
            WRITE (out,110) thetaN0(ng), 'thetaN0',                     &
     &            '',            &
     &            ''
            WRITE (out,110) tref(ng), 'tref',                           &
     &            '',            &
     &            ''
            WRITE (out,110) e_a(ng), 'e_a',                             &
     &            '',            &
     &            ''
            WRITE (out,110) p_cref(ng), 'p_cref',                       &
     &            '',            &
     &            ''
            WRITE (out,110) z_umax1(ng), 'z_umax1',                     &
     &            '',            &
     &            ''
            WRITE (out,110) z_umax2(ng), 'z_umax2',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p1_excr(ng), 'p1_excr',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p2_excr(ng), 'p2_excr',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p1_mort2(ng), 'p1_mort2',                   &
     &            '',            &
     &            ''
            WRITE (out,110) p2_mort2(ng), 'p2_mort2',                   &
     &            '',            &
     &            ''
            WRITE (out,110) z_mort(ng), 'z_mort',                       &
     &            '',            &
     &            ''
            WRITE (out,110) z_mort2(ng), 'z_mort2',                     &
     &            '',            &
     &            ''
            WRITE (out,110) d1_rem(ng), 'd1_rem',                       &
     &            '',            &
     &            ''
            WRITE (out,110) d2_rem(ng), 'd2_rem',                       &
     &            '',            &
     &            ''
            WRITE (out,110) wDet2(ng), 'wDet2',                         &
     &            '',            &
     &            ''
            WRITE (out,110) p1_kno3(ng), 'p1_kno3',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p2_kno3(ng), 'p2_kno3',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p1_knh4(ng), 'p1_knh4',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p2_knh4(ng), 'p2_knh4',                     &
     &            '',            &
     &            ''
            WRITE (out,110) p2_si(ng), 'p2_si',                         &
     &            '',            &
     &            ''
            WRITE (out,110) p_nitr(ng), 'p_nitr',                       &
     &            '',            &
     &            ''
            WRITE (out,110) falloc0(ng), 'falloc0',                     &
     &            '',            &
     &            ''
            WRITE (out,110) si_alloc(ng), 'si_alloc',                   &
     &            '',            &
     &            ''
            WRITE (out,110) z_grz(ng), 'z_grz',                         &
     &            '',            &
     &            ''
            WRITE (out,110) zegest0(ng), 'zegest0',                     &
     &            '',            &
     &            ''
            WRITE (out,110) r_ref(ng), 'r_ref',                         &
     &            '',            &
     &            ''
            WRITE (out,110) q_min(ng), 'q_min',                         &
     &            '',            &
     &            ''
            WRITE (out,110) q_max(ng), 'q_max',                         &
     &            '',            &
     &            ''
            WRITE (out,110) q_simin(ng), 'q_simin',                     &
     &            '',            &
     &            ''
            WRITE (out,110) q_simax(ng), 'q_simax',                     &
     &            '',            &
     &            ''
            WRITE (out,110) vncref(ng), 'vncref',                       &
     &            '',            &
     &            ''
            WRITE (out,110) biolambda(ng), 'biolambda',                 &
     &            '',            &
     &            ''
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' LimaDoney Model Parameters, Grid: ',i2.2,              &
     &        /,  ' =================================',/)  
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t28,a)
  80  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t28,a)
 110  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',       &
     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t28,a,1x,a)

# elif defined NEMURO
!
!-----------------------------------------------------------------------
!  Read in Fasham biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN 
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'AttSW') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSW)
          ELSE IF (TRIM(KeyWord).eq.'AttPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttPS)
          ELSE IF (TRIM(KeyWord).eq.'AttPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttPL)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'AlphaPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaPS)
          ELSE IF (TRIM(KeyWord).eq.'AlphaPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaPL)
          ELSE IF (TRIM(KeyWord).eq.'BetaPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaPS)
          ELSE IF (TRIM(KeyWord).eq.'BetaPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaPL)
          ELSE IF (TRIM(KeyWord).eq.'VmaxS') THEN
            Npts=load_r(Nval, Rval, Ngrids, VmaxS)
          ELSE IF (TRIM(KeyWord).eq.'VmaxL') THEN
            Npts=load_r(Nval, Rval, Ngrids, VmaxL)
          ELSE IF (TRIM(KeyWord).eq.'KNO3S') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNO3S)
          ELSE IF (TRIM(KeyWord).eq.'KNO3L') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNO3L)
          ELSE IF (TRIM(KeyWord).eq.'KNH4S') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNH4S)
          ELSE IF (TRIM(KeyWord).eq.'KNH4L') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNH4L)
          ELSE IF (TRIM(KeyWord).eq.'KSiL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KSiL)
          ELSE IF (TRIM(KeyWord).eq.'PusaiS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiS)
          ELSE IF (TRIM(KeyWord).eq.'PusaiL') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiL)
          ELSE IF (TRIM(KeyWord).eq.'KGppS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGppS)
          ELSE IF (TRIM(KeyWord).eq.'KGppL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGppL)
          ELSE IF (TRIM(KeyWord).eq.'ResPS0') THEN
            Npts=load_r(Nval, Rval, Ngrids, ResPS0)
          ELSE IF (TRIM(KeyWord).eq.'ResPL0') THEN
            Npts=load_r(Nval, Rval, Ngrids, ResPL0)
          ELSE IF (TRIM(KeyWord).eq.'KResPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KResPS)
          ELSE IF (TRIM(KeyWord).eq.'KResPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KResPL)
          ELSE IF (TRIM(KeyWord).eq.'GammaS') THEN
            Npts=load_r(Nval, Rval, Ngrids, GammaS)
          ELSE IF (TRIM(KeyWord).eq.'GammaL') THEN
            Npts=load_r(Nval, Rval, Ngrids, GammaL)
          ELSE IF (TRIM(KeyWord).eq.'MorPS0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorPS0)
          ELSE IF (TRIM(KeyWord).eq.'MorPL0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorPL0)
          ELSE IF (TRIM(KeyWord).eq.'KMorPS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorPS)
          ELSE IF (TRIM(KeyWord).eq.'KMorPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorPL)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxSps') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxSps)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxLps') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxLps)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxLpl') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxLpl)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxLzs') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxLzs)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxPpl') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxPpl)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxPzs') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxPzs)
          ELSE IF (TRIM(KeyWord).eq.'GRmaxPzl') THEN
            Npts=load_r(Nval, Rval, Ngrids, GRmaxPzl)
          ELSE IF (TRIM(KeyWord).eq.'KGraS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGraS)
          ELSE IF (TRIM(KeyWord).eq.'KGral') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGral)
          ELSE IF (TRIM(KeyWord).eq.'KGraP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KGraP)
          ELSE IF (TRIM(KeyWord).eq.'LamS') THEN
            Npts=load_r(Nval, Rval, Ngrids, LamS)
          ELSE IF (TRIM(KeyWord).eq.'LamL') THEN
            Npts=load_r(Nval, Rval, Ngrids, LamL)
          ELSE IF (TRIM(KeyWord).eq.'LamP') THEN
            Npts=load_r(Nval, Rval, Ngrids, LamP)
          ELSE IF (TRIM(KeyWord).eq.'KPS2ZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPS2ZS)
          ELSE IF (TRIM(KeyWord).eq.'KPS2ZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPS2ZL)
          ELSE IF (TRIM(KeyWord).eq.'KPL2ZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPL2ZL)
          ELSE IF (TRIM(KeyWord).eq.'KZS2ZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KZS2ZL)
          ELSE IF (TRIM(KeyWord).eq.'KPL2ZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KPL2ZP)
          ELSE IF (TRIM(KeyWord).eq.'KZS2ZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KZS2ZP)
          ELSE IF (TRIM(KeyWord).eq.'KZL2ZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KZL2ZP)
          ELSE IF (TRIM(KeyWord).eq.'PS2ZSstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PS2ZSstar)
          ELSE IF (TRIM(KeyWord).eq.'PS2ZLstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PS2ZLstar)
          ELSE IF (TRIM(KeyWord).eq.'PL2ZLstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PL2ZLstar)
          ELSE IF (TRIM(KeyWord).eq.'ZS2ZLstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZS2ZLstar)
          ELSE IF (TRIM(KeyWord).eq.'PL2ZPstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, PL2ZPstar)
          ELSE IF (TRIM(KeyWord).eq.'ZS2ZPstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZS2ZPstar)
          ELSE IF (TRIM(KeyWord).eq.'ZL2ZPstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZL2ZPstar)
          ELSE IF (TRIM(KeyWord).eq.'PusaiPL') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiPL)
          ELSE IF (TRIM(KeyWord).eq.'PusaiZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PusaiZS)
          ELSE IF (TRIM(KeyWord).eq.'MorZS0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorZS0)
          ELSE IF (TRIM(KeyWord).eq.'MorZL0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorZL0)
          ELSE IF (TRIM(KeyWord).eq.'MorZP0') THEN
            Npts=load_r(Nval, Rval, Ngrids, MorZP0)
          ELSE IF (TRIM(KeyWord).eq.'KMorZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorZS)
          ELSE IF (TRIM(KeyWord).eq.'KMorZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorZL)
          ELSE IF (TRIM(KeyWord).eq.'KMorZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, KMorZP)
          ELSE IF (TRIM(KeyWord).eq.'AlphaZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaZS)
          ELSE IF (TRIM(KeyWord).eq.'AlphaZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, AlphaZL)
          ELSE IF (TRIM(KeyWord).eq.'BetaZS') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaZS)
          ELSE IF (TRIM(KeyWord).eq.'BetaZL') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaZL)
          ELSE IF (TRIM(KeyWord).eq.'BetaZP') THEN
            Npts=load_r(Nval, Rval, Ngrids, BetaZP)
          ELSE IF (TRIM(KeyWord).eq.'Nit0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Nit0)
          ELSE IF (TRIM(KeyWord).eq.'VP2N0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VP2N0)
          ELSE IF (TRIM(KeyWord).eq.'VP2D0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VP2D0)
          ELSE IF (TRIM(KeyWord).eq.'VD2N0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VD2N0)
          ELSE IF (TRIM(KeyWord).eq.'VO2S0') THEN
            Npts=load_r(Nval, Rval, Ngrids, VO2S0)
          ELSE IF (TRIM(KeyWord).eq.'KNit') THEN
            Npts=load_r(Nval, Rval, Ngrids, KNit)
          ELSE IF (TRIM(KeyWord).eq.'KP2D') THEN
            Npts=load_r(Nval, Rval, Ngrids, KP2D)
          ELSE IF (TRIM(KeyWord).eq.'KP2N') THEN
            Npts=load_r(Nval, Rval, Ngrids, KP2N)
          ELSE IF (TRIM(KeyWord).eq.'KD2N') THEN
            Npts=load_r(Nval, Rval, Ngrids, KD2N)
          ELSE IF (TRIM(KeyWord).eq.'KO2S') THEN
            Npts=load_r(Nval, Rval, Ngrids, KO2S)
          ELSE IF (TRIM(KeyWord).eq.'RSiN') THEN
            Npts=load_r(Nval, Rval, Ngrids, RSiN)
          ELSE IF (TRIM(KeyWord).eq.'setVPON') THEN
            Npts=load_r(Nval, Rval, Ngrids, setVPON)
          ELSE IF (TRIM(KeyWord).eq.'setVOpal') THEN
            Npts=load_r(Nval, Rval, Ngrids, setVOpal)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTvar(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTsur(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,100) AttSW(ng), 'AttSW',                         &
     &            'Light attenuation due to seawater (m-1)'
            WRITE (out,110) AttPS(ng), 'AttPS',                         &
     &            'Light attenuation due to small phytoplankton',       &
     &            '(m2/millimole_N).'
            WRITE (out,110) AttPL(ng), 'AttPL',                         &
     &            'Light attenuation due to large phytoplankton',       &
     &            '(m2/millimole_N).'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) AlphaPS(ng), 'AlphaPS',                     &
     &            'Small phytoplankton initial slope of the P-I curve', &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,110) AlphaPL(ng), 'AlphaPL',                     &
     &            'Small phytoplankton initial slope of the P-I curve', &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,110) BetaPS(ng), 'BetaPS',                       &
     &            'Small phytoplankton photoinhibition coefficient',    &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,110) BetaPL(ng), 'BetaPL',                       &
     &            'Large phytoplankton photoinhibition coefficient',    &
     &            '(1/(W/m2) 1/day).'
            WRITE (out,110) VmaxS(ng), 'VmaxS',                         &
     &            'Small phytoplankton maximum photosynthetic rate',    &
     &            '(1/day).'
            WRITE (out,110) VmaxL(ng), 'VmaxL',                         &
     &            'Large phytoplankton maximum photosynthetic rate',    &
     &            '(1/day).'
            WRITE (out,110) KNO3S(ng), 'KNO3S',                         &
     &            'Small phytoplankton NO3 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,110) KNO3L(ng), 'KNO3L',                         &
     &            'Large phytoplankton NO3 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,110) KNH4S(ng), 'KNH4S',                         &
     &            'Small phytoplankton NH4 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,110) KNH4L(ng), 'KNH4L',                         &
     &            'Large phytoplankton NH4 half saturation constant',   &
     &            '(millimole_N/m3).'
            WRITE (out,110) KSiL(ng), 'KSiL',                           &
     &            'Small phytoplankton SiOH4 half saturation constant', &
     &            '(millimole_Si/m3).'
            WRITE (out,110) PusaiS(ng), 'PusaiS',                       &
     &            'Small phytoplankton NH4 inhibition coefficient',     &
     &            '(m3/millimole_N).'
            WRITE (out,110) PusaiL(ng), 'PusaiL',                       &
     &            'Large phytoplankton NH4 inhibition coefficient',     &
     &            '(m3/millimole_N).'
            WRITE (out,110) KGppS(ng), 'KGppS',                         &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'photosynthetic rate (1/Celsius).'
            WRITE (out,110) KGppL(ng), 'KGppL',                         &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'photosynthetic rate (1/Celsius).'
            WRITE (out,100) ResPS0(ng), 'ResPS0',                       &
     &            'Small phytoplankton respiration rate (1/day).'
            WRITE (out,100) ResPL0(ng), 'ResPL0',                       &
     &            'Large phytoplankton respiration rate (1/day).'
            WRITE (out,110) KResPS(ng), 'KResPS',                       &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'respiration (1/Celsius).'
            WRITE (out,110) KResPL(ng), 'KResPL',                       &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'respiration (1/Celsius).'
            WRITE (out,110) GammaS(ng), 'GammaS',                       &
     &            'Small phytoplankton ratio of extracellular',         &
     &            'excretion to photosynthesis (nondimensional).'
            WRITE (out,110) GammaL(ng), 'GammaL',                       &
     &            'Large phytoplankton ratio of extracellular',         &
     &            'excretion to photosynthesis (nondimensional).'
            WRITE (out,110) MorPS0(ng), 'MorPS0',                       &
     &            'Small phytoplankton mortality rate',                 &
     &            '(m3/millimole_N/day).'
            WRITE (out,110) MorPS0(ng), 'MorPL0',                       &
     &            'Large phytoplankton mortality rate',                 &
     &            '(m3/millimole_N/day).'
            WRITE (out,110) KMorPS(ng), 'KMorPS',                       &
     &            'Small phytoplankton temperature coefficient for',    &
     &            'mortality (1/Celsius).'
            WRITE (out,110) KMorPL(ng), 'KMorPL',                       &
     &            'Large phytoplankton temperature coefficient for',    &
     &            'mortality (1/Celsius).'
            WRITE (out,110) GRmaxSps(ng), 'GRmaxSps',                   &
     &            'Small zooplankton grazing rate on small',            &
     &            'phytoplankton (1/day).'
            WRITE (out,110) GRmaxLps(ng), 'GRmaxLps',                   &
     &            'Large zooplankton grazing rate on small',            &
     &            'phytoplankton (1/day).'
            WRITE (out,110) GRmaxLpl(ng), 'GRmaxLpl',                   &
     &            'Large zooplankton grazing rate on large',            &
     &            'phytoplankton (1/day).'
            WRITE (out,110) GRmaxLzs(ng), 'GRmaxLzs',                   &
     &            'Large zooplankton grazing rate on small',            &
     &            'zooplankton (1/day).'
            WRITE (out,110) GRmaxPpl(ng), 'GRmaxPpl',                   &
     &            'Predator zooplankton grazing rate on large',         &
     &            'phytoplankton (1/day).'
            WRITE (out,110) GRmaxPzs(ng), 'GRmaxPzs',                   &
     &            'Predator zooplankton grazing rate on small',         &
     &            'zooplankton (1/day).'
            WRITE (out,110) GRmaxPzl(ng), 'GRmaxPzl',                   &
     &            'Predator zooplankton grazing rate on large',         &
     &            'zooplankton (1/day).'
            WRITE (out,110) KGraS(ng), 'KGraS',                         &
     &            'Small zooplankton temperature coefficient for',      &
     &            'grazing (1/Celsius).'
            WRITE (out,110) KGraL(ng), 'KGraL',                         &
     &            'Large zooplankton temperature coefficient for',      &
     &            'grazing (1/Celsius).'
            WRITE (out,110) KGraP(ng), 'KGraP',                         &
     &            'Predator zooplankton temperature coefficient for',   &
     &            'grazing (1/Celsius).'
            WRITE (out,110) LamS(ng), 'LamS',                           &
     &            'Small zooplankton grazing Ivlev constant',           &
     &            '(m3/millimole_N).'
            WRITE (out,110) LamL(ng), 'LamL',                           &
     &            'Large zooplankton grazing Ivlev constant',           &
     &            '(m3/millimole_N).'
            WRITE (out,110) LamP(ng), 'LamP',                           &
     &            'Preditor zooplankton grazing Ivlev constant',        &
     &            '(m3/millimole_N).'
#  ifdef HOLLING_GRAZING
            WRITE (out,110) KPS2ZS(ng), 'KPS2ZS',                       &
     &            'Half-saturation constant for small zooplankton',     &
     &            'grazing on small phytoplankton (millimole_N/m3)^2.'
            WRITE (out,110) KPS2ZL(ng), 'KPS2ZL',                       &
     &            'Half-saturation constant for large zooplankton',     &
     &            'grazing on small phytoplankton (millimole_N/m3)^2.'
            WRITE (out,110) KPL2ZL(ng), 'KPL2ZL',                       &
     &            'Half-saturation constant for large zooplankton',     &
     &            'grazing on large phytoplankton (millimole_N/m3)^2.'
            WRITE (out,110) KPL2ZP(ng), 'KPL2ZP',                       &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on large phytoplankton (millimole_N/m3)^2.'
            WRITE (out,110) KZS2ZP(ng), 'KZS2ZP',                       &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on small zooplankton (millimole_N/m3)^2.'
            WRITE (out,110) KZL2ZP(ng), 'KZL2ZP',                       &
     &            'Half-saturation constant for predator zooplankton',  &
     &            'grazing on large zooplankton (millimole_N/m3)^2.'
#  else
            WRITE (out,110) PS2ZSstar(ng), 'PS2ZSstar',                 &
     &            'Small zooplankton threshold for grazing on small',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,110) PS2ZLstar(ng), 'PS2ZLstar',                 &
     &            'Large zooplankton threshold for grazing on small',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,110) PL2ZLstar(ng), 'PL2ZLstar',                 &
     &            'Large zooplankton threshold for grazing on large',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,110) PL2ZLstar(ng), 'PL2ZLstar',                 &
     &            'Large zooplankton threshold for grazing on small',   &
     &            'zooplankton (millimole_N/m3).'
            WRITE (out,110) PL2ZPstar(ng), 'PL2ZPstar',                 &
     &            'Predator zooplankton threshold for grazing on large',&
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,110) ZS2ZPstar(ng), 'ZS2ZPstar',                 &
     &            'Predator zooplankton threshold for grazing on small',&
     &            'zooplankton (millimole_N/m3).'
            WRITE (out,110) ZL2ZPstar(ng), 'ZL2ZPstar',                 &
     &            'Predator zooplankton threshold for grazing on large',&
     &            'zooplankton (millimole_N/m3).'
#  endif
            WRITE (out,110) PusaiPL(ng), 'PusauPL',                     &
     &            'Predator zooplankton grazing inhibition on large',   &
     &            'phytoplankton (millimole_N/m3).'
            WRITE (out,110) PusaiZS(ng), 'PusauZS',                     &
     &            'Predator zooplankton grazing inhibition on small',   &
     &            'zootoplankton (millimole_N/m3).'
            WRITE (out,110) MorZS0(ng), 'MorZS0',                       &
     &            'Small zooplankton mortality rate at 0 Celsius',      &
     &            '(m3/millimole_N/day).'
            WRITE (out,110) MorZL0(ng), 'MorZL0',                       &
     &            'Large zooplankton mortality rate at 0 Celsius',      &
     &            '(m3/millimole_N/day).'
            WRITE (out,110) MorZP0(ng), 'MorZP0',                       &
     &            'Predator zooplankton mortality rate at 0 Celsius',   &
     &            '(m3/millimole_N/day).'
            WRITE (out,110) KMorZS(ng), 'KMorZS',                       &
     &            'Small zooplankton temperature coefficient for',      &
     &            'mortality (1/Celsius).'
            WRITE (out,110) KMorZL(ng), 'KMorZL',                       &
     &            'Large zooplankton temperature coefficient for',      &
     &            'mortality (1/Celsius).'
            WRITE (out,110) KMorZP(ng), 'KMorZP',                       &
     &            'Predator zooplankton temperature coefficient for',   &
     &            'mortality (1/Celsius).'
            WRITE (out,110) AlphaZS(ng), 'AlphaZS',                     &
     &            'Small zooplankton assimilation efficiency',          &
     &            '(nondimensional).'
            WRITE (out,110) AlphaZS(ng), 'AlphaZL',                     &
     &            'Large zooplankton assimilation efficiency',          &
     &            '(nondimensional).'
            WRITE (out,110) AlphaZP(ng), 'AlphaZP',                     &
     &            'Predator zooplankton assimilation efficiency',       &
     &            '(nondimensional).'
            WRITE (out,110) BetaZS(ng), 'BetaZS',                       &
     &            'Small zooplankton growth efficiency',                &
     &            '(nondimensional).'
            WRITE (out,110) BetaZL(ng), 'BetaZL',                       &
     &            'Large zooplankton growth efficiency',                &
     &            '(nondimensional).'
            WRITE (out,110) BetaZP(ng), 'BetaZP',                       &
     &            'Predator zooplankton growth efficiency',             &
     &            '(nondimensional).'
            WRITE (out,100) Nit0(ng), 'Nit0',                           &
     &            'NH4 to NO3 decomposition rate (1/day).'
            WRITE (out,100) VP2N0(ng), 'VP2N0',                         &
     &            'PON to NH4 decomposition rate (1/day).'
            WRITE (out,100) VP2D0(ng), 'VP2D0',                         &
     &            'PON to DON decomposition rate (1/day).'
            WRITE (out,100) VD2N0(ng), 'VD2N0',                         &
     &            'DON to NH4 decomposition rate (1/day).'
            WRITE (out,100) VP2N0(ng), 'VO2S0',                         &
     &            'Opal to SiOH4 decomposition rate (1/day).'
            WRITE (out,110) KNit(ng), 'KNit',                           &
     &            'Temperature coefficient for NH4 to NO3',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,110) KP2D(ng), 'KP2D',                           &
     &            'Temperature coefficient for PON to DON',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,110) KP2N(ng), 'KP2N',                           &
     &            'Temperature coefficient for PON to NH4',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,110) KP2N(ng), 'KD2N',                           &
     &            'Temperature coefficient for DON to NH4',             &
     &            'decomposition (1/Celsius).'
            WRITE (out,110) KO2S(ng), 'KO2S',                           &
     &            'Temperature coefficient for Opal to SiOH4',          &
     &            'decomposition (1/Celsius).'
            WRITE (out,100) RSiN(ng), 'RSiN',                           &
     &            'Si:N ratio (millimole_Si/millimole_N)'
            WRITE (out,100) setVPON(ng), 'setVPON',                     &
     &            'PON sinking velocity (m/day).'
            WRITE (out,100) setVOpal(ng), 'setVOpal',                   &
     &            'Opal sinking velocity (m/day).'
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Nemuro Model Parameters, Grid: ',i2.2,                &
     &        /,  ' =============================',/)
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)


# elif defined NPZD_FRANKS
!
!-----------------------------------------------------------------------
!  Read in NPZD biological model parameters.
!-----------------------------------------------------------------------
!
#  ifdef ANA_BIOLOGY
      IF (.not.allocated(BioIni)) allocate ( BioIni(MT,Ngrids) )
#  endif
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN 
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
#  ifdef ANA_BIOLOGY
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iNO3_)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iNO3_,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iPhyt)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iPhyt,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iZoop)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iZoop,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iSDet)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iSDet,1))
#  endif
          ELSE IF (TRIM(KeyWord).eq.'K_ext') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_ext)
          ELSE IF (TRIM(KeyWord).eq.'K_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NO3)
          ELSE IF (TRIM(KeyWord).eq.'K_Phy') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_Phy)
          ELSE IF (TRIM(KeyWord).eq.'Vm_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vm_NO3)
          ELSE IF (TRIM(KeyWord).eq.'PhyMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMR)
          ELSE IF (TRIM(KeyWord).eq.'ZooGR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGR)
          ELSE IF (TRIM(KeyWord).eq.'ZooMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMR)
          ELSE IF (TRIM(KeyWord).eq.'ZooMD') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMD)
          ELSE IF (TRIM(KeyWord).eq.'ZooGA') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGA)
          ELSE IF (TRIM(KeyWord).eq.'ZooEC') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooEC)
          ELSE IF (TRIM(KeyWord).eq.'DetRR') THEN
            Npts=load_r(Nval, Rval, Ngrids, DetRR)
          ELSE IF (TRIM(KeyWord).eq.'wDet') THEN
            Npts=load_r(Nval, Rval, Ngrids, wDet)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTvar(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTsur(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
#  ifdef ANA_BIOLOGY
            WRITE (out,100) BioIni(iNO3_,ng), 'BioIni(iNO3_)',          &
     &            'Nitrate initial concentration (mmol/m3).'
            WRITE (out,100) BioIni(iPhyt,ng), 'BioIni(iPhyt)',          &
     &            'Phytoplankton initial concentration (mmol/m3).'
            WRITE (out,100) BioIni(iZoop,ng), 'BioIni(iZoop)',          &
     &            'Zooplankton initial concentration (mmol/m3).'
            WRITE (out,100) BioIni(iSDet,ng), 'BioIni(iSDet)',          &
     &            'Small detritus initial concentration (mmol/m3).'
#  endif
            WRITE (out,100) K_ext(ng), 'K_ext',                         &
     &            'Light extinction coefficient (m-1).'
            WRITE (out,110) K_NO3(ng), 'K_NO3',                         &
     &            'Inverse half-saturation for phytoplankton NO3',      &
     &            'uptake (1/(mmol m-3)).'
            WRITE (out,110) K_Phy(ng), 'K_Phy',                         &
     &            'Phytoplankton saturation coefficient',               &
     &            '(mmol/m3)^2.'
            WRITE (out,100) Vm_NO3(ng), 'Vm_NO3',                       &
     &            'Nitrate upatake rate (day-1).'
            WRITE (out,100) PhyMR(ng), 'PhyMR',                         &
     &            'Phytoplankton senescence/mortality rate (day-1)'
            WRITE (out,100) ZooGR(ng), 'ZooGR',                         &
     &            'Zooplankton maximum growth rate (day-1).'
            WRITE (out,100) ZooMR(ng), 'ZooMR',                         &
     &            'Zooplankton mortality rate (day-1).'
            WRITE (out,100) ZooMD(ng), 'ZooMD',                         &
     &            'Zooplankton death bits rate (day-1).'
            WRITE (out,100) ZooGA(ng), 'ZooGA',                         &
     &            'Zooplankton grazing inefficiency (nondimensional).'
            WRITE (out,100) ZooEC(ng), 'ZooEC',                         &
     &            'Zooplankton excreted fraction (nondimensional).'
            WRITE (out,100) DetRR(ng), 'DetRR',                         &
     &            'Detritus remineralization rate (day-1).'
            WRITE (out,100) wDet(ng), 'wDet',                           &
     &            'Detrital sinking rate (m/day).'
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' NPZD Model Parameters, Grid: ',i2.2,                &
     &        /,  ' ===============================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)

# elif defined NPZD_POWELL
!
!-----------------------------------------------------------------------
!  Read in NPZD biological model parameters.
!-----------------------------------------------------------------------
!
#  ifdef ANA_BIOLOGY
      IF (.not.allocated(BioIni)) allocate ( BioIni(MT,Ngrids) )
#  endif
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN 
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
#  ifdef ANA_BIOLOGY
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iNO3_)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iNO3_,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iPhyt)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iPhyt,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iZoop)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iZoop,1))
          ELSE IF (TRIM(KeyWord).eq.'BioIni(iSDet)') THEN
            Npts=load_r(Nval, Rval, Ngrids, BioIni(iSDet,1))
#  endif
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'AttSW') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSW)
          ELSE IF (TRIM(KeyWord).eq.'AttPhy') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttPhy)
          ELSE IF (TRIM(KeyWord).eq.'PhyIS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIS)
          ELSE IF (TRIM(KeyWord).eq.'Vm_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vm_NO3)
          ELSE IF (TRIM(KeyWord).eq.'PhyMRD') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMRD)
          ELSE IF (TRIM(KeyWord).eq.'PhyMRN') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMRN)
          ELSE IF (TRIM(KeyWord).eq.'K_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NO3)
          ELSE IF (TRIM(KeyWord).eq.'Ivlev') THEN
            Npts=load_r(Nval, Rval, Ngrids, Ivlev)
          ELSE IF (TRIM(KeyWord).eq.'ZooGR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGR)
          ELSE IF (TRIM(KeyWord).eq.'ZooEED') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooEED)
          ELSE IF (TRIM(KeyWord).eq.'ZooEEN') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooEEN)
          ELSE IF (TRIM(KeyWord).eq.'ZooMRD') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMRD)
          ELSE IF (TRIM(KeyWord).eq.'ZooMRN') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMRN)
          ELSE IF (TRIM(KeyWord).eq.'DetRR') THEN
            Npts=load_r(Nval, Rval, Ngrids, DetRR)
          ELSE IF (TRIM(KeyWord).eq.'wPhy') THEN
            Npts=load_r(Nval, Rval, Ngrids, wPhy)
          ELSE IF (TRIM(KeyWord).eq.'wDet') THEN
            Npts=load_r(Nval, Rval, Ngrids, wDet)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTvar(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTsur(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
#  ifdef ANA_BIOLOGY
            WRITE (out,100) BioIni(iNO3_,ng), 'BioIni(iNO3_)',          &
     &            'Nitrate initial concentration (mmol/m3).'
            WRITE (out,100) BioIni(iPhyt,ng), 'BioIni(iPhyt)',          &
     &            'Phytoplankton initial concentration (mmol/m3).'
            WRITE (out,100) BioIni(iZoop,ng), 'BioIni(iZoop)',          &
     &            'Zooplankton initial concentration (mmol/m3).'
            WRITE (out,100) BioIni(iSDet,ng), 'BioIni(iSDet)',          &
     &            'Small detritus initial concentration (mmol/m3).'
#  endif
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,100) AttSW(ng), 'AttSW',                         &
     &            'Light attenuation of seawater (m-1).'
            WRITE (out,100) AttPhy(ng), 'AttPhy',                       &
     &            'Light attenuation by phytoplankton (m2/mmole_N).'
            WRITE (out,110) PhyIS(ng), 'PhyIS',                         &
     &            'Phytoplankton growth, initial slope of P-I curve',   &
     &            '(m2/W).'
            WRITE (out,100) Vm_NO3(ng), 'Vm_NO3',                       &
     &            'Nitrate upatake rate (day-1).'
            WRITE (out,100) PhyMRD(ng), 'PhyMRD',                       &
     &            'Phytoplankton mortality rate to Detritus (day-1)'
            WRITE (out,100) PhyMRN(ng), 'PhyMRN',                       &
     &            'Phytoplankton mortality rate to Nitrogen (day-1)'
            WRITE (out,110) K_NO3(ng), 'K_NO3',                         &
     &            'Inverse half-saturation for phytoplankton NO3',      &
     &            'uptake (1/(mmol m-3)).'
            WRITE (out,110) Ivlev(ng), 'Ivlev',                         &
     &            'Ivlev constant for zooplankton grazing',             &
     &            '(nondimensional).'
            WRITE (out,100) ZooGR(ng), 'ZooGR',                         &
     &            'Zooplankton maximum growth rate (day-1).'
            WRITE (out,110) ZooEED(ng), 'ZooEED',                       &
     &            'Zooplankton excretion efficiency to Detritus',       &
     &            'pool (nondimensional).'
            WRITE (out,110) ZooEEN(ng), 'ZooEEN',                       &
     &            'Zooplankton excretion efficiency to Nitrogen',       &
     &            'pool (nondimensional).'
            WRITE (out,100) ZooMRD(ng), 'ZooMRD',                       &
     &            'Zooplankton mortality rate to Detritus (day-1).'
            WRITE (out,100) ZooMRN(ng), 'ZooMRN',                       &
     &            'Zooplankton mortality rate to Nitrogen (day-1).'
            WRITE (out,100) DetRR(ng), 'DetRR',                         &
     &            'Detritus remineralization rate (day-1).'
            WRITE (out,100) wPhy(ng), 'wPhy',                           &
     &            'Phytoplankton sinking rate (m/day).'
            WRITE (out,100) wDet(ng), 'wDet',                           &
     &            'Detrital sinking rate (m/day).'
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' NPZD Model Parameters, Grid: ',i2.2,                &
     &        /,  ' ===============================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)

# elif defined ECOSIM
!
!-----------------------------------------------------------------------
!  Read in EcoSim bio-optical model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'RtUVR_flag') THEN
            Npts=load_l(Nval, Cval, Ngrids, RtUVR_flag)
          ELSE IF (TRIM(KeyWord).eq.'NFIX_flag') THEN
            Npts=load_l(Nval, Cval, Ngrids, NFIX_flag)
          ELSE IF (TRIM(KeyWord).eq.'Regen_flag') THEN
            Npts=load_l(Nval, Cval, Ngrids, Regen_flag)
          ELSE IF (TRIM(KeyWord).eq.'HsNO3') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsNO3)
          ELSE IF (TRIM(KeyWord).eq.'HsNH4') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsNH4)
          ELSE IF (TRIM(KeyWord).eq.'HsSiO') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsSiO)
          ELSE IF (TRIM(KeyWord).eq.'HsPO4') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsPO4)
          ELSE IF (TRIM(KeyWord).eq.'HsFe') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsFe)
          ELSE IF (TRIM(KeyWord).eq.'GtALG_max') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, GtALG_max)
          ELSE IF (TRIM(KeyWord).eq.'PhyTbase') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, PhyTbase)
          ELSE IF (TRIM(KeyWord).eq.'PhyTfac') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, PhyTfac)
          ELSE IF (TRIM(KeyWord).eq.'BET_') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, BET_)
          ELSE IF (TRIM(KeyWord).eq.'maxC2nALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2nALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2nALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2nALG)
          ELSE IF (TRIM(KeyWord).eq.'C2nALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2nALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'maxC2SiALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2SiALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2SiALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2SiALG)
          ELSE IF (TRIM(KeyWord).eq.'C2SiALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2SiALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'maxC2pALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2pALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2pALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2pALG)
          ELSE IF (TRIM(KeyWord).eq.'C2pALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2pALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'maxC2FeALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, maxC2FeALG)
          ELSE IF (TRIM(KeyWord).eq.'minC2FeALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, minC2FeALG)
          ELSE IF (TRIM(KeyWord).eq.'C2FeALGminABS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2FeALGminABS)
          ELSE IF (TRIM(KeyWord).eq.'qu_yld') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, qu_yld)
          ELSE IF (TRIM(KeyWord).eq.'E0_comp') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, E0_comp)
          ELSE IF (TRIM(KeyWord).eq.'E0_inhib') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, E0_inhib)
          ELSE IF (TRIM(KeyWord).eq.'inhib_fac') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, inhib_fac)
          ELSE IF (TRIM(KeyWord).eq.'C2CHL_max') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2CHL_max)
          ELSE IF (TRIM(KeyWord).eq.'mxC2Cl') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxC2Cl)
          ELSE IF (TRIM(KeyWord).eq.'b_C2Cl') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_C2Cl)
          ELSE IF (TRIM(KeyWord).eq.'mxC2Cn') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxC2Cn)
          ELSE IF (TRIM(KeyWord).eq.'b_C2Cn') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_C2Cn)
          ELSE IF (TRIM(KeyWord).eq.'mxPacEff') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPacEff)
          ELSE IF (TRIM(KeyWord).eq.'b_PacEff') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PacEff)
          ELSE IF (TRIM(KeyWord).eq.'mxChlB') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxChlB)
          ELSE IF (TRIM(KeyWord).eq.'b_ChlB') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_ChlB)
          ELSE IF (TRIM(KeyWord).eq.'mxChlC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxChlC)
          ELSE IF (TRIM(KeyWord).eq.'b_ChlC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_ChlC)
          ELSE IF (TRIM(KeyWord).eq.'mxPSC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPSC)
          ELSE IF (TRIM(KeyWord).eq.'b_PSC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PSC)
          ELSE IF (TRIM(KeyWord).eq.'mxPPC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxPPC)
          ELSE IF (TRIM(KeyWord).eq.'b_PPC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_PPC)
          ELSE IF (TRIM(KeyWord).eq.'mxLPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxLPUb)
          ELSE IF (TRIM(KeyWord).eq.'b_LPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_LPUb)
          ELSE IF (TRIM(KeyWord).eq.'mxHPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, mxHPUb)
          ELSE IF (TRIM(KeyWord).eq.'b_HPUb') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, b_HPUb)
          ELSE IF (TRIM(KeyWord).eq.'FecDOC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, FecDOC)
          ELSE IF (TRIM(KeyWord).eq.'FecPEL') THEN
            Npts=load_r(Nval, Rval, Nphy*Nfec*Ngrids, FecPEL)
          ELSE IF (TRIM(KeyWord).eq.'FecCYC') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, FecCYC)
          ELSE IF (TRIM(KeyWord).eq.'ExALG') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, ExALG)
          ELSE IF (TRIM(KeyWord).eq.'WS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, WS)
          ELSE IF (TRIM(KeyWord).eq.'HsGRZ') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsGRZ)
          ELSE IF (TRIM(KeyWord).eq.'MinRefuge') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, MinRefuge)
          ELSE IF (TRIM(KeyWord).eq.'RefugeDep') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, RefugeDep)
          ELSE IF (TRIM(KeyWord).eq.'Norm_Vol') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, Norm_Vol)
          ELSE IF (TRIM(KeyWord).eq.'Norm_Surf') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, Norm_Surf)
          ELSE IF (TRIM(KeyWord).eq.'HsDOP') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsDOP)
          ELSE IF (TRIM(KeyWord).eq.'C2pALKPHOS') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, C2pALKPHOS)
          ELSE IF (TRIM(KeyWord).eq.'HsDON') THEN
            Npts=load_r(Nval, Rval, Nphy*Ngrids, HsDON)
          ELSE IF (TRIM(KeyWord).eq.'C2nNupDON') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2nNupDON)
          ELSE IF (TRIM(KeyWord).eq.'C2nBAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2nBAC)
          ELSE IF (TRIM(KeyWord).eq.'C2pBAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2pBAC)
          ELSE IF (TRIM(KeyWord).eq.'C2FeBAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, C2FeBAC)
          ELSE IF (TRIM(KeyWord).eq.'HsDOC_ba') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, HsDOC_ba)
          ELSE IF (TRIM(KeyWord).eq.'GtBAC_max') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, GtBAC_max)
          ELSE IF (TRIM(KeyWord).eq.'BacTbase') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, BacTbase)
          ELSE IF (TRIM(KeyWord).eq.'BacTfac') THEN
            Npts=load_r(Nval, Rval, Nbac*Ngrids, BacTfac)
          ELSE IF (TRIM(KeyWord).eq.'BacDOC') THEN
            Npts=load_r(Nval, Rval, Ngrids, BacDOC)
          ELSE IF (TRIM(KeyWord).eq.'BacPEL') THEN
            Npts=load_r(Nval, Rval, Ngrids, BacPEL)
          ELSE IF (TRIM(KeyWord).eq.'BacCYC') THEN
            Npts=load_r(Nval, Rval, Ngrids, BacCYC)
          ELSE IF (TRIM(KeyWord).eq.'ExBAC_c') THEN
            Npts=load_r(Nval, Rval, Ngrids, ExBAC_c)
          ELSE IF (TRIM(KeyWord).eq.'ExBacC2N') THEN
            Npts=load_r(Nval, Rval, Ngrids, ExBacC2N)
          ELSE IF (TRIM(KeyWord).eq.'Bac_Ceff') THEN
            Npts=load_r(Nval, Rval, Ngrids, Bac_Ceff)
          ELSE IF (TRIM(KeyWord).eq.'RtNIT') THEN
            Npts=load_r(Nval, Rval, Ngrids, RtNIT)
          ELSE IF (TRIM(KeyWord).eq.'HsNIT') THEN
            Npts=load_r(Nval, Rval, Ngrids, HsNIT)
          ELSE IF (TRIM(KeyWord).eq.'cDOCfrac_c') THEN
            Npts=load_r(Nval, Rval, Ndom*Ngrids, cDOCfrac_c)
          ELSE IF (TRIM(KeyWord).eq.'RtUVR_DIC') THEN
            Npts=load_r(Nval, Rval, Ngrids, RtUVR_DIC)
          ELSE IF (TRIM(KeyWord).eq.'RtUVR_DOC') THEN
            Npts=load_r(Nval, Rval, Ngrids, RtUVR_DOC)
          ELSE IF (TRIM(KeyWord).eq.'WF') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, WF)
          ELSE IF (TRIM(KeyWord).eq.'RegTbase') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegTbase)
          ELSE IF (TRIM(KeyWord).eq.'RegTfac') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegTfac)
          ELSE IF (TRIM(KeyWord).eq.'RegCR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegCR)
          ELSE IF (TRIM(KeyWord).eq.'RegNR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegNR)
          ELSE IF (TRIM(KeyWord).eq.'RegSR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegSR)
          ELSE IF (TRIM(KeyWord).eq.'RegPR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegPR)
          ELSE IF (TRIM(KeyWord).eq.'RegFR') THEN
            Npts=load_r(Nval, Rval, Nfec*Ngrids, RegFR)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,150) 'idTvar(idbio(', itrc, '))'
                  STOP
                END IF
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  WRITE (out,120) 'idTsur(idbio(', itrc, '))'
                  STOP
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Initialize secondary EcoSim parameters.
!-----------------------------------------------------------------------
!
      CALL initialize_biology
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,60) RtUVR_flag(ng), 'RtUVR_flag',                &
     &            'Switch to calculate CDOC UV photolysis.'
            WRITE (out,60) NFIX_flag(ng), 'NFIX_flag',                  &
     &            'Switch to calculate temperature based N fixation.'
            WRITE (out,60) Regen_flag(ng), 'Regen_flag',                &
     &            'Switch to calculate fecal matter regeneration.'
            WRITE (out,70) 'HsNO3',                                     &
     &            'Half-saturation for phytoplankton NO3 uptake',       &
     &            '(micromole_NO3/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsNO3(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsNH4',                                     &
     &            'Half-saturation for phytoplankton NH4 uptake',       &
     &            '(micromole_NH4/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsNH4(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsSiO',                                     &
     &            'Half-saturation for phytoplankton SiO uptake',       &
     &            '(micromole_SiO/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsSiO(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsPO4',                                     &
     &            'Half-saturation for phytoplankton PO4 uptake',       &
     &            '(micromole_PO4/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsPO4(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsFe',                                      &
     &            'Half-saturation for phytoplankton Fe uptake',        &
     &            '(micromole_Fe/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsFe(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'GtALG_max',                                 &
     &            'Maximum 24 hour growth rate (1/day):'
              DO is=1,Nphy
                WRITE (out,90) GtALG_max(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'PhyTbase',                                  &
     &            'Temperature base for exponential response to',       &
     &            'temperature (Celsius):'
              DO is=1,Nphy
                WRITE (out,90) PhyTbase(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'PhyTfac',                                   &
     &            'Phytoplankton exponential temperature factor',       &
     &            '(1/Celsius):'
              DO is=1,Nphy
                WRITE (out,90) PhyTfac(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'BET_',                                      &
     &            'Nitrate uptake inhibition for NH4 (l/micromole):'
              DO is=1,Nphy
                WRITE (out,90) BET_(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'maxC2nALG',                                 &
     &            'Maximum phytoplankton C:N ratio',                    &
     &            '(micromole_C/micromole_N):'   
              DO is=1,Nphy
                WRITE (out,90) maxC2nALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'minC2nALG',                                 &
     &            'Balanced phytoplankton C:N ratio',                   &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) minC2nALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2nALGminABS',                              &
     &            'Absolute minimum phytoplankton C:N ratio',           &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) C2nALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'maxC2SiALG',                                &
     &            'Maximum phytoplankton C:Si ratio',                   &
     &            '(micromole_C/micromole_Si)'
              DO is=1,Nphy
                WRITE (out,90) maxC2SiALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'minC2SiALG',                                &
     &            'Balanced phytoplankton C:Si ratio',                  &
     &            '(micromole_C/micromole_Si):'
              DO is=1,Nphy
                WRITE (out,90) minC2SiALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2SiALGminABS',                             &
     &            'Absolute minimum phytoplankton C:Si ratio',          &
     &            '(micromole_C/micromole_Si):'
              DO is=1,Nphy
                WRITE (out,90) C2SiALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'maxC2pALG',                                 &
     &            'Maximum phytoplankton C:P ratio',                    &
     &            '(micromole_C/micromole_P):'
              DO is=1,Nphy
                WRITE (out,90) maxC2pALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'minC2pALG',                                 &
     &            'Balanced phytoplankton C:P ratio',                   &
     &            '(micromole_C/micromole_P):'
              DO is=1,Nphy
                WRITE (out,90) minC2pALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2pALGminABS',                              &
     &            'Absolute minimum phytoplankton C:P ratio',           &
     &            '(micromole_C/micromole_P)'
              DO is=1,Nphy
                WRITE (out,90) C2pALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'maxC2FeALG',                                &
     &            'Maximum phytoplankton C:Fe ratio',                   &
     &            '(micromole_C/micromole_Fe):'
              DO is=1,Nphy
                WRITE (out,90) maxC2FeALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'minC2FeALG',                                &
     &            'Balanced phytoplankton C:Fe ratio',                  &
     &            '(micromole_C/micromole_Fe):'
              DO is=1,Nphy
                WRITE (out,90) minC2FeALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2FeALGminABS',                             &
     &            'Absolute minimum phytoplankton C:Fe ratio',          &
     &            '(micromole_C/micromole_Fe):'
              DO is=1,Nphy
                WRITE (out,90) C2FeALGminABS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'qu_yld',                                    &
     &            'Maximum quantum yield',                              &
     &            '(micromole_C/micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) qu_yld(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'E0_comp',                                   &
     &            'Compensation light level (micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) E0_comp(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'E0_inhib',                                  &
     &            'Light level for onset of photoinhibition',           &
     &            '(micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) E0_inhib(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'inhib_fac',                                 &
     &            'Exponential decay factor for light limited growth',  &
     &            '(1/micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) inhib_fac(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2CHL_max',                                 &
     &            'Maximum lighted limited C:Chl ratio',                &
     &            '(microgram_C/microgram_Chl):'
              DO is=1,Nphy
                WRITE (out,90) C2CHL_max(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxC2Cl',                                    &
     &            'Rate of change in light limited C:Chl ratio',        &
     &            '(microgram_C/microgram_Chl/micromole_quanta):'
              DO is=1,Nphy
                WRITE (out,90) mxC2Cl(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_C2Cl',                                    &
     &            'Minimum lighted limited C:Chl ratio',                &
     &            '(microgram_C/microgram_Chl):'
              DO is=1,Nphy
                WRITE (out,90) b_C2Cl(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxC2Cn',                                    &
     &            'Rate of change in nutient limited C:Chl ratio',      &
     &            '[(ug_C/ug_Chl)/(umole_C/umole_N)]:'
              DO is=1,Nphy
                WRITE (out,90) mxC2Cn(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_C2Cn',                                    &
     &            'Minimum nutrient limited C:Chl ratio',               &
     &            '(microgram_C/microgram_Chl):'
              DO is=1,Nphy
                WRITE (out,90) b_C2Cn(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxPacEff',                                  &
     &            'Rate of change in package effect',                   &
     &            '[1/(microgram_C/microgram_Chl)]:'
              DO is=1,Nphy
                WRITE (out,90) mxPacEff(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_PacEff',                                  &
     &            'Maximum package effect',                             &
     &            '[1/(microgram_C/microgram_Chl)]:'
              DO is=1,Nphy
                WRITE (out,90) b_PacEff(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxChlB',                                    &
     &            'Rate of change in the Chl_b:Chl_a ratio',            &
     &            '[(ug_Chl_b/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxChlB(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_ChlB',                                    &
     &            'Maximum Chl_b:Chl_a ratio',                          &
     &            '(microgram_Chl_b/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_ChlB(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxChlC',                                    &
     &            'Rate of change in the Chl_c:Chl_a ratio',            &
     &            '[(ug_Chl_c/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxChlC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_ChlC',                                    &
     &            'Maximum Chl_c:Chl_a ratio',                          &
     &            '(microgram_Chl_c/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_ChlC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxPSC',                                     &
     &            'Rate of change in the PSC:Chl_a ratio',              &
     &            '[(ug_PSC/ug_Chl_a)/ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxPSC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_PSC',                                     &
     &            'Maximum PSC:Chl_a ratio',                            &
     &            '(microgram_PSC/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_PSC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxPPC',                                     &
     &            'Rate of change in the PPC:Chl_a ratio',              &
     &            '[(ug_PPC/ug_Chl_a)/(ug_C/ug_Chl_ a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxPPC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_PPC',                                     &
     &            'Maximum PPC:Chl_a ratio',                            &
     &            '(microgram_Chl_c/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_PPC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxLPUb',                                    &
     &            'Rate of change in the LPUb:Chl_a ratio',             &
     &            '[(ug_LPUb/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxLPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'b_LPUb',                                    &
     &            'Maximum LPUb:Chl_a ratio',                           &
     &            '(migrogram_HPUb/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_LPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'mxHPUb',                                    &
     &            'Rate of change in the HPUb:Chl_a ratio',             &
     &            '[(ug_HPUb/ug_Chl_a)/(ug_C/ug_Chl_a)]:'
              DO is=1,Nphy
                WRITE (out,90) mxHPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70)'b_HPUb',                                     &
     &            'Maximum HPUb:Chl_a ratio',                           &
     &            '(microgram_HPUb/microgram_Chl_a):'
              DO is=1,Nphy
                WRITE (out,90) b_HPUb(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'FecDOC',                                    &
     &            'Proportion of grazing stress apportioned to DOM',    &
     &            '(nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) FecDOC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'FecPEL',                                    &
     &            'Proportion of grazing stress apportioned to fecal',  &
     &            '(nondimensional):'
              DO i=1,Nfec
                DO is=1,Nphy
                  WRITE (out,100) FecPEL(is,i,ng), i, TRIM(PhyName(is))
                END DO
              END DO
            WRITE (out,70) 'FecCYC',                                    &
     &            'Proportion of grazing stress that is recycled',      &
     &            '(nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) FecCYC(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'ExALG',                                     &
     &            'Proportion of daily production lost to excretion',   &
     &            '(nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) ExALG(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'WS',                                        &
     &            'Phytoplankton sinking speed (meters/day):'
              DO is=1,Nphy
                WRITE (out,90) WS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'HsGRZ',                                     &
     &            'Phytoplankton grazing parameter (nondimensional):'     
              DO is=1,Nphy
                WRITE (out,90) HsGRZ(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'MinRefuge',                                 &
     &            'Refuge Phytoplankton population (micromole_C/liter):'
              DO is=1,Nphy
                WRITE (out,90) MinRefuge(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'RefugeDep',                                 &
     &            'Maximum Refuge Phytoplankton depth (meters):'
              DO is=1,Nphy
                WRITE (out,90) RefugeDep(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'Norm_Vol',                                  &
     &            'Normalized Volume factor (nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) Norm_Vol(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,80) 'Norm_Surf',                                 &
     &            'Normalized Surface Area factor (nondimensional):'
              DO is=1,Nphy
                WRITE (out,90) Norm_Surf(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsDOP',                                     &
     &            'Half Saturation Constant for DOP uptake',            &
     &            '(micromole_DOP/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsDOP(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2pALKPHOS',                                &
     &            'C:P ratio where DOP uptake begins',                  &
     &            '(micromole_C/micromole_P):'
              DO is=1,Nphy
                WRITE (out,90) C2pALKPHOS(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsDON',                                     &
     &            'Half Saturation Constant for DON uptake',            &
     &            '(micromole_DON/liter):'
              DO is=1,Nphy
                WRITE (out,90) HsDON(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'C2nNupDON',                                 &
     &            'C:N ratio where DON uptake begins',                  &
     &            '(micromole_C/micromole_N):'
              DO is=1,Nphy
                WRITE (out,90) C2nNupDON(is,ng), TRIM(PhyName(is))
              END DO
            WRITE (out,70) 'HsDOC_ba',                                  &
     &            'Half saturation constant for bacteria DOC uptake',   &
     &            '(micromole_DOC/liter):'
              DO is=1,Nbac
                WRITE (out,90) HsDOC_ba(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,80) 'GtBAC_max',                                 &
     &            'Maximum 24 hour bacterial growth rate (1/day):'
              DO is=1,Nbac
                WRITE (out,90) GtBAC_max(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,70) 'BacTbase',                                  &
     &            'Temperature base for exponetial response to',        &
     &            'temperature (Celsius):'
              DO is=1,Nbac
                WRITE (out,90) BacTbase(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,70) 'BacTfac',                                   &
     &            'Bacteria exponential temperature factor',            &
     &            '(1/Celsius):'
              DO is=1,Nbac
                WRITE (out,90) BacTfac(is,ng), TRIM(BacName(is))
              END DO
            WRITE (out,110) C2nBAC(ng), 'C2nBAC',                       &
     &            'Carbon to Nitrogen ratio of Bacteria',               &
     &            '(micromole_C/micromole_N).'
            WRITE (out,110) C2pBAC(ng), 'C2pBAC',                       &
     &            'Carbon to Phosphorus ratio of Bacteria',             &
     &            '(micromole_C/micromole_P).'
            WRITE (out,110) C2FeBAC(ng), 'C2FeBAC',                     &
     &            'Carbon to Iron ratio of Bacteria',                   &
     &            '(micromole_C/micromole_Fe).'
            WRITE (out,110) BacDOC(ng), 'BacDOC',                       &
     &            'Proportion of bacteria grazing stress apportioned',  &
     &            'to DOM (nondimensional).'
            WRITE (out,110) BacPEL(ng), 'BacPEL',                       &
     &            'Proportion of bacteria grazing stress apportioned',  &
     &            'to fecal (nondimensional).'
            WRITE (out,110) BacCYC(ng), 'BacCYC',                       &
     &            'Proportion of bacteria grazing stress recycled',     &
     &            '(nondimensional).'
            WRITE (out,110) ExBAC_c(ng), 'ExBAC_c',                     &
     &            'Bacterial recalcitrant C excretion as proportion',   &
     &            'of uptake (nondimensional).'
            WRITE (out,110) ExBacC2N(ng), 'ExBacC2N',                   &
     &            'Bacterial recalcitrant excretion carbon:nitrogen',   &
     &            'ratio (micromole_C/micromole_N).'
            WRITE (out,110) Bac_Ceff(ng), 'Bac_Ceff',                   &
     &            'Bacterial gross growth carbon efficiency',           &
     &            '(nondimensional).'
            WRITE (out,120) RtNIT(ng), 'RtNIT',                         &
     &            'Maximum nitrification rate (1/day).'
            WRITE (out,110) HsNIT(ng), 'HsNIT',                         &
     &            'Half saturation constant for bacteria nitrification',&
     &            '(micromole_NH4/liter).'
            WRITE (out,70) 'cDOCfrac_c',                                &
     &            'Colored fraction of DOC from phytoplakton and',      &
     &            'bacterial losses (nondimensional):'
              DO is=1,Ndom
                WRITE (out,90) cDOCfrac_c(is,ng), TRIM(DomName(is))
              END DO
            WRITE (out,110) RtUVR_DIC(ng), 'RtUVR_DIC',                 &
     &            'UV degradation of DOC into DIC at 410 nm',           &
     &            '(micromole/meter/liter/hour).'
            WRITE (out,110) RtUVR_DOC(ng), 'RtUVR_DOC',                 &
     &            'UV degradation of DOC to colorless labile DOC',      &
     &            'at 410 nm (micromole/meter/liter/hour).'
            WRITE (out,80) 'WF',                                        &
     &            'Fecal sinking flux (meters/day):'
              DO is=1,Nfec
                WRITE (out,90) WF(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,70) 'RegTbase',                                  &
     &            'Fecal regeneration temperature base for exponential',&
     &            'response to temperature (Celsius):'
              DO is=1,Nfec
                WRITE (out,90) RegTbase(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,70) 'RegTfac',                                   &
     &            'Fecal regeneration exponential temperature factor',  &
     &            '(1/Celsius):'
              DO is=1,Nfec
                WRITE (out,90) RegTfac(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegCR',                                     &
     &            'Fecal carbon regeneration rate (1/day):'
              DO is=1,Nfec
                WRITE (out,90) RegCR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegNR',                                     &
     &            'Fecal nitrogen regeneration rate (1/day)'
              DO is=1,Nfec
                WRITE (out,90) RegNR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegSR',                                     &
     &            'Fecal silica regeneration rate (1/day):'
              DO is=1,Nfec
                WRITE (out,90) RegSR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegPR',                                     &
     &            'Fecal phosphorus regeneration rate (1/day):'
              DO is=1,Nfec
                WRITE (out,90) RegPR(is,ng), TRIM(FecName(is))
              END DO
            WRITE (out,80) 'RegFR',                                     &
     &            'Fecal iron regeneration rate (1/day)'
              DO is=1,Nfec
                WRITE (out,90) RegFR(is,ng), TRIM(FecName(is))
              END DO
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,130) tnu2(i,ng), 'tnu2', i,                    &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,130) tnu4(i,ng), 'tnu4', i,                    &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,130) Akt_bak(i,ng), 'Akt_bak', i,               &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,130) Tnudg(i,ng), 'Tnudg', i,                  &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,140)                   &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' EcoSim Parameters, Grid: ',i2.2,                    &
     &        /,  ' ===========================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a)
  70  FORMAT ('...........',2x,a,t30,a,/,t32,a)
  80  FORMAT ('...........',2x,a,t30,a)
  90  FORMAT (1p,e11.4,t33,a) 
 100  FORMAT (1p,e11.4,t33,'Fecal Group ',i1,', ',a) 
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a) 
 120  FORMAT (1p,e11.4,2x,a,t30,a)
 130  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 140  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
 150  FORMAT (/,' read_BioPar - variable info not yet loaded, ',       &
     &        a,i2.2,a)

# endif

      RETURN
      END SUBROUTINE read_BioPar
#endif

#ifdef SEDIMENT
      SUBROUTINE read_SedPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in cohesive and non-cohesive sediment model      !
!  parameters.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
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
      integer :: Npts, Nval, i, ng, itrc, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(Ngrids) :: Lbed
      logical, dimension(MBOTP,Ngrids) :: Lbottom
      logical, dimension(NCS,Ngrids) :: Lmud
      logical, dimension(NNS,Ngrids) :: Lsand

      real(r8), dimension(Ngrids) :: Rbed
      real(r8), dimension(NCS,Ngrids) :: Rmud
      real(r8), dimension(NNS,Ngrids) :: Rsand
 
      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in cohesive and non-cohesive model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lsediment') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lsediment)
          ELSE IF (TRIM(KeyWord).eq.'NEWLAYER_THICK') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              newlayer_thick(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MINLAYER_THICK') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              minlayer_thick(ng)=Rbed(ng)
            END DO
# ifdef MIXED_BED
          ELSE IF (TRIM(KeyWord).eq.'TRANSC') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              transC(ng)=Rbed(ng)
          END DO
          ELSE IF (TRIM(KeyWord).eq.'TRANSN') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              transN(ng)=Rbed(ng)
          END DO
# endif
          ELSE IF (TRIM(KeyWord).eq.'BEDLOAD_COEFF') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              bedload_coeff(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_SD50') THEN 
            IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Sd50(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_CSED') THEN
            IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud )
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Csed(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_SRHO') THEN
            IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Srho(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_WSED') THEN
            IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Wsed(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_ERATE') THEN
            IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Erate(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAU_CE') THEN
            IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                tau_ce(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAU_CD') THEN
            IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                tau_cd(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_POROS') THEN
            IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                poros(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNU2') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                tnu2(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNU4') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                tnu4(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                Akt_bak(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNUDG') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                Tnudg(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_MORPH_FAC') THEN
            IF (.not.allocated(morph_fac)) THEN
              allocate (morph_fac(NST,Ngrids))
            END IF
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                morph_fac(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
#if defined COHESIVE_BED || defined MIXED_BED
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_MIN') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_min(ng)=Rbed(ng)
          END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_MAX') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_max(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_SLOPE') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_slp(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_OFF') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_off(ng)=Rbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAUCR_TIME') THEN 
            Npts=load_r(Nval, Rval, Ngrids, Rbed)
            DO ng=1,Ngrids
              tcr_tim(ng)=Rbed(ng)
          END DO 
# endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idmud)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idTvar(idsed(itrc))
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMfrac)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idfrac(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMmass)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idBmas(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
# ifdef BEDLOAD
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMUbld)') THEN
            DO ng=1,Ngrids
              DO itrc=1,NCS
                IF (idUbld(itrc).eq.0) THEN
                  WRITE (out,100) 'idUbld'
                 STOP
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idUbld(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMVbld)') THEN
            DO ng=1,Ngrids
              DO itrc=1,NCS
                IF (idVbld(itrc).eq.0) THEN
                  WRITE (out,100) 'idVbld'
                  STOP
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idVbld(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
# endif
          ELSE IF (TRIM(KeyWord).eq.'SAND_SD50') THEN 
            IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Sd50(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_CSED') THEN
            IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand )
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Csed(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_SRHO') THEN
            IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Srho(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_WSED') THEN
            IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Wsed(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_ERATE') THEN
            IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Erate(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TAU_CE') THEN
            IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                tau_ce(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TAU_CD') THEN
            IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                tau_cd(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_POROS') THEN
            IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                poros(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNU2') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                tnu2(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNU4') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                tnu4(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                Akt_bak(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNUDG') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                Tnudg(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_MORPH_FAC') THEN
            IF (.not.allocated(morph_fac)) THEN
              allocate (morph_fac(NST,Ngrids))
            END IF
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                morph_fac(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idsand)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idTvar(idsed(NCS+itrc))
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSfrac)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idfrac(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSmass)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idBmas(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
# ifdef BEDLOAD
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSUbld)') THEN
            DO ng=1,Ngrids
              DO itrc=NCS+1,NST
                IF (idUbld(itrc).eq.0) THEN
                  WRITE (out,100) 'idUbld'
                 STOP
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idUbld(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSVbld)') THEN
            DO ng=1,Ngrids
              DO itrc=NCS+1,NST
                IF (idVbld(itrc).eq.0) THEN
                  WRITE (out,100) 'idVbld'
                  STOP
                END IF
              END DO
            END DO
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idVbld(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
# endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(ithck)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(ithck)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iaged)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(iaged)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iporo)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(iporo)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
# if defined COHESIVE_BED || defined SED_BIODIFF || defined MIXED_BED
          ELSE IF (TRIM(KeyWord).eq.'Hout(ibtcr)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbed) 
            i=idSbed(ibtcr) 
            DO ng=1,Ngrids 
              Hout(i,ng)=Lbed(ng) 
            END DO
# endif
          ELSE IF (TRIM(KeyWord).eq.'Hout(idiff)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(idiff)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,40) line
      STOP
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lsediment(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60)
            DO itrc=1,NST
              WRITE (out,70) itrc, Sd50(itrc,ng), Csed(itrc,ng),        &
     &                       Srho(itrc,ng), Wsed(itrc,ng),              &
     &                       Erate(itrc,ng), poros(itrc,ng)
            END DO
            WRITE (out,80)
            DO itrc=1,NST
              i=idsed(itrc)
              WRITE (out,70) itrc, tau_ce(itrc,ng), tau_cd(itrc,ng),    &
     &                       tnu2(i,ng), tnu4(i,ng), Akt_bak(i,ng),     &
     &                       Tnudg(i,ng)
            END DO
            WRITE (out,85)
            DO itrc=1,NST
              WRITE (out,70) itrc,  morph_fac(itrc,ng)
            END DO
            WRITE (out,110) newlayer_thick(ng)
            WRITE (out,111) minlayer_thick(ng)
            WRITE (out,120) bedload_coeff(ng)
# ifdef MIXED_BED
            WRITE (out,121) transC(ng)
            WRITE (out,122) transN(ng)
# endif
            DO itrc=1,NST
              i=idTvar(idsed(itrc))
              IF (Hout(i,ng)) WRITE (out,90) Hout(i,ng),                &
     &            'Hout(idTvar)',                                       &
     &            'Write out sediment', itrc, TRIM(Vname(1,i))
            END DO
            DO itrc=1,NST
              i=idfrac(itrc)
              IF (Hout(i,ng)) WRITE (out,90) Hout(i,ng),                &
     &            'Hout(idfrac)',                                       &
     &            'Write out bed fraction, sediment ', itrc,            &
     &            TRIM(Vname(1,i))
            END DO
            DO itrc=1,NST
              i=idBmas(itrc)
              IF (Hout(i,ng)) WRITE (out,90) Hout(i,ng),                &
     &            'Hout(idfrac)',                                       &
     &            'Write out mass, sediment ', itrc,                    &
     &            TRIM(Vname(1,i))
            END DO
# ifdef BEDLOAD
            DO itrc=1,NST
              i=idUbld(itrc)
              IF (Hout(i,ng)) WRITE (out,90) Hout(i,ng),                &
     &            'Hout(idUbld)',                                       &
     &            'Write out bed load at U-points, sediment ', itrc,    &
     &            TRIM(Vname(1,i))
              i=idVbld(itrc)
              IF (Hout(i,ng)) WRITE (out,90) Hout(i,ng),                &
     &            'Hout(idVbld)',                                       &
     &            'Write out bed load at V-points, sediment ', itrc,    &
     &            TRIM(Vname(1,i))
            END DO
# endif
            DO itrc=1,MBEDP
              i=idSbed(itrc)
              IF (Hout(i,ng)) WRITE (out,90) Hout(i,ng),                &
     &            'Hout(idSbed)',                                       &
     &            'Write out BED property ', itrc, TRIM(Vname(1,i))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant input parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NST
          Sd50(i,ng)=Sd50(i,ng)*0.001_r8
          Wsed(i,ng)=Wsed(i,ng)*0.001_r8
          tau_ce(i,ng)=tau_ce(i,ng)/rho0
          tau_cd(i,ng)=tau_cd(i,ng)/rho0
          tnu4(idsed(i),ng)=SQRT(ABS(tnu4(idsed(i),ng)))
          IF (Tnudg(idsed(i),ng).gt.0.0_r8) THEN
            Tnudg(idsed(i),ng)=1.0_r8/(Tnudg(idsed(i),ng)*86400.0_r8)
          ELSE
            Tnudg(idsed(i),ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' READ_SedPar - variable info not yet loaded, ', a)
  40  FORMAT (/,' READ_SedPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' Sediment Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  60  FORMAT (/,1x,'Size',5x,'Sd50',8x,'Csed',8x,'Srho',8x,'Wsed',      &
     &        8x,'Erate',7x,'poros',/,1x,'Class',4x,'(mm)',7x,          &
     &        '(kg/m3)',5x,'(kg/m3)',5x,'(mm/s)',5x,'(kg/m2/s)',4x,     &
     &        '(nondim)',/)
  70  FORMAT (2x,i2,2x,6(1x,1p,e11.4))
  80  FORMAT (/,9x,'tau_ce',6x,'tau_cd',7x,'tnu2',8x,'tnu4',7x,         &
     &        'Akt_bak',6x,'Tnudg',/,9x,'(N/m2)',6x,'(N/m2)',6x,        &
     &        '(m2/s)',6x,'(m4/s)',7x,'(m2/s)',6x,'(day)',/)
  85  FORMAT (/,9x,'morph_fac',/,9x,'(nondim)',/)
  90  FORMAT (10x,l1,2x,a,t29,a,i2.2,':',1x,a)
 100  FORMAT (/,' READ_SEDPAR - variable info not yet loaded, ', a)
 110  FORMAT (/,' New bed layer formed when deposition exceeds ',e11.5, &
     &        ' (m).')
 111  FORMAT (' Two first layers are combined when 2nd layer smaller ', &
     &         'than ',e11.5,' (m).')
 120  FORMAT (' Rate coefficient for bed load transport = ',e11.5,/)
 121  FORMAT (' Transition for mixed sediment =',e11.5,/)
 122  FORMAT (' Transition for cohesive sediment =',e11.5,/)
      RETURN
      END SUBROUTINE read_SedPar
#endif

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
      integer :: i, igrid, itrc, k, ng, status

      integer :: decode_line, load_i, load_l, load_r

# if defined ASSIMILATION || defined NUDGING
      logical, dimension(MT,Ngrids) :: Lassi

      real(r8), dimension(MT,Ngrids) :: Rassi
# endif
# if defined FOUR_DVAR || defined VERIFICATION
      logical, dimension(MT) :: Ltracer

      real(r8), dimension(MT,Ngrids) :: Rtracer
# endif
      real(r8), dimension(100) :: Rval

      character (len=1 ), parameter :: blank = ' '
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
          ELSE IF (TRIM(KeyWord).eq.'Znudass') THEN
            Npts=load_r(Nval, Rval, Ngrids, Znudass)
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
 10   WRITE (out,50) line
      STOP
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
#  ifdef NUDGING_SSH
            WRITE (out,135) Znudass(ng), 'Znudass',                     &
     &            'Nudging assimilation time scale (days)',             &
     &            'for free-surface.'
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
                STOP
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
                STOP
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
#  if defined ASSIMILATION_SSH || defined NUDGING_SSH
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
            STOP
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
        IF (Znudass(ng).gt.0.0_r8) THEN
          Znudass(ng)=1.0_r8/(Znudass(ng)*86400.0_r8)
        ELSE
          Znudass(ng)=0.0_r8
        END IF
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
          IF (TRIM(KeyWord).eq.'LcycleINI') THEN 
            Npts=load_l(Nval, Cval, Ngrids, LcycleINI)
          ELSE IF (TRIM(KeyWord).eq.'ICG') THEN
            Npts=load_i(Nval, Rval, 1, ICG)
          ELSE IF (TRIM(KeyWord).eq.'NiterSD') THEN
            Npts=load_i(Nval, Rval, 1, NiterSD)
          ELSE IF (TRIM(KeyWord).eq.'CGtol') THEN
            Npts=load_r(Nval, Rval, 1, CGtol)
          ELSE IF (TRIM(KeyWord).eq.'CGstepI') THEN
            Npts=load_r(Nval, Rval, 1, CGstepI)
          ELSE IF (TRIM(KeyWord).eq.'CGeps') THEN
            Npts=load_r(Nval, Rval, 1, CGeps)
          ELSE IF (TRIM(KeyWord).eq.'CostFunFac') THEN
            Npts=load_r(Nval, Rval, 1, CostFunFac)
          ELSE IF (TRIM(KeyWord).eq.'GradErr') THEN
            Npts=load_r(Nval, Rval, 1, GradErr)
          ELSE IF (TRIM(KeyWord).eq.'HevecErr') THEN
            Npts=load_r(Nval, Rval, 1, HevecErr)
          ELSE IF (TRIM(KeyWord).eq.'LhessianEV') THEN
            Npts=load_l(Nval, Cval, 1, LhessianEV)
          ELSE IF (TRIM(KeyWord).eq.'Lprecond') THEN
            Npts=load_l(Nval, Cval, 1, Lprecond)
#  if defined WEAK_CONSTRAINT
            IF ( LhessianEV.and.Lprecond ) THEN
              LhessianEV=.FALSE.
            END IF
#  endif
          ELSE IF (TRIM(KeyWord).eq.'NritzEV') THEN
            Npts=load_i(Nval, Rval, 1, NritzEV)
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
            Npts=load_l(Nval, Cval, Ngrids, LdefNRM)
          ELSE IF (TRIM(KeyWord).eq.'LwrtNRM') THEN 
            Npts=load_l(Nval, Cval, Ngrids, LwrtNRM)
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isFsur)') THEN
            IF (isFsur.eq.0) THEN
              WRITE (out,170) 'isFsur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isFsur))
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isUbar)') THEN
            IF (isUbar.eq.0) THEN
              WRITE (out,170) 'isUbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isUbar))
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isVbar)') THEN
            IF (isVbar.eq.0) THEN
              WRITE (out,170) 'isVbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isVbar))
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isUstr)') THEN
            IF (isUstr.eq.0) THEN
              WRITE (out,170) 'isUstr'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isUstr))
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isVstr)') THEN
            IF (isVstr.eq.0) THEN
              WRITE (out,170) 'isVstr'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isVstr))
#  endif
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isUvel)') THEN
            IF (isUvel.eq.0) THEN
              WRITE (out,170) 'isUvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isUvel))
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isVvel)') THEN
            IF (isVvel.eq.0) THEN
              WRITE (out,170) 'isVvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, 1, Cnorm(isVvel))
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isTvar)') THEN
            IF (MAXVAL(isTvar).eq.0) THEN
              WRITE (out,170) 'isTvar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, MT, Ltracer)
            DO itrc=1,MT
              i=isTvar(itrc)
              Cnorm(i)=Ltracer(itrc)
            END DO
#  ifdef ADJUST_STFLUX
          ELSE IF (TRIM(KeyWord).eq.'Cnorm(isTsur)') THEN
            IF (MAXVAL(isTsur).eq.0) THEN
              WRITE (out,170) 'isTsur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, MT, Ltracer)
            DO itrc=1,MT
              i=isTsur(itrc)
              Cnorm(i)=Ltracer(itrc)
            END DO
#  endif
          ELSE IF (TRIM(KeyWord).eq.'Nmethod') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nmethod)
          ELSE IF (TRIM(KeyWord).eq.'Rscheme') THEN
            Npts=load_i(Nval, Rval, Ngrids, Rscheme)
          ELSE IF (TRIM(KeyWord).eq.'Nrandom') THEN
            Npts=load_i(Nval, Rval, 1, Nrandom)
          ELSE IF (TRIM(KeyWord).eq.'Hgamma') THEN
            Npts=load_r(Nval, Rval, 1, Hgamma)
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Vgamma') THEN
            Npts=load_r(Nval, Rval, 1, Vgamma)
#  endif
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isFsur)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isFsur,1))
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isUbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isUbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isVbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isVbar,1))
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isUstr)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isUstr,1))
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isVstr)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isVstr,1))
#  endif
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isUvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Hdecay(isVvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Hdecay(isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#   ifdef ADJUST_STFLUX
          ELSE IF (TRIM(KeyWord).eq.'Hdecay(isTsur)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Hdecay(isTsur(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#   endif
          ELSE IF (TRIM(KeyWord).eq.'Vdecay(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vdecay(isUvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Vdecay(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vdecay(isVvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Vdecay(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Vdecay(isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#  endif
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isFsur)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isFsur,1))
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isUbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isUbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isVbar)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isVbar,1))
#  ifdef ADJUST_WSTRESS
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isUstr)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isUstr,1))
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isVstr)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isVstr,1))
#  endif
#  ifdef SOLVE3D
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isUvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isUvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isVvel)') THEN
            Npts=load_r(Nval, Rval, Ngrids, Cfscale(isVvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isTvar)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Cfscale(isTvar(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#   ifdef ADJUST_STFLUX
          ELSE IF (TRIM(KeyWord).eq.'Cfscale(isTsur)') THEN
            Npts=load_r(Nval, Rval, MT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Cfscale(isTsur(itrc),ng)=Rtracer(itrc,ng)
              END DO
            END DO
#   endif
#  endif
          ELSE IF (TRIM(KeyWord).eq.'npN2') THEN
            Npts=load_i(Nval, Rval, 1, npN2)
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_N2') THEN
            allocate ( pcoef_N2(0:npN2) )
            DO i=0,npN2
              READ (inp,*) k, pcoef_N2(k)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'STDname') THEN
            DO i=1,LEN(STDname(Nval))
              STDname(Nval)(i:i)=blank
            END DO
            STDname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NRMname') THEN
            DO i=1,LEN(NRMname(Nval))
              NRMname(Nval)(i:i)=blank
            END DO
            NRMname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'OBSname') THEN
            DO i=1,LEN(OBSname(Nval))
              OBSname(Nval)(i:i)=blank
            END DO
            OBSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'HSSname') THEN
            DO i=1,LEN(HSSname(Nval))
              HSSname(Nval)(i:i)=blank
            END DO
            HSSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'MODname') THEN
            DO i=1,LEN(MODname(Nval))
              MODname(Nval)(i:i)=blank
            END DO
            MODname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
 10   WRITE (out,50) line
      STOP
 20   CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
#  if defined FOUR_DVAR && !defined OBS_SENSITIVITY
          WRITE (out,60) ng
          WRITE (out,70) LcycleINI(ng), 'LcycleINI',                    &
     &            'Switch to recycle time-record in initial file.'
#   if !(defined WEAK_CONSTRAINT || defined IOM)
#    ifndef IS4DVAR
          WRITE (out,90) ICG, 'ICG',                                    &
     &            'Conjugate gradient descent scheme: ',                &
     &            '[0] Fletcher-reeves, [1] Polak-Ribiere'
          WRITE (out,80) NiterSD, 'NiterSD',                            &
     &            'Number of iterations between a steepest descent.'
          WRITE (out,100) CGtol, 'CGtol',                               &
     &            'Conjugate gradient orthogonality tolerance.'
#    endif
          WRITE (out,100) CGstepI, 'CGstepI',                           &
     &            'Initial descent algorithm step size.'
          WRITE (out,100) CostFunFac, 'CostFunFac',                     &
     &            'Cost function change convergence factor.'
#   else
          WRITE (out,100) CGeps, 'CGeps',                               &
     &            'Conjugate gradient convergence criteria.'
#    ifdef RPM_RELAXATION
          WRITE (out,120) tl_M2diff, 'tl_M2diff',                       &
     &            'RPM 2D momentum diffusive relaxation coefficient.'
#     ifdef SOLVE3D
          WRITE (out,130) tl_M3diff, 'tl_M3diff',                       &
     &            'RPM 3D momentum diffusive relaxation coefficient.'
          DO itrc=1,NT(ng)
            WRITE (out,130) tl_Tdiff(itrc,ng), 'tl_Tdiff',              &
     &            'RPM tracers diffusive relaxation coefficient, ',     &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     endif
#    endif
#   endif
#   if defined IS4DVAR && defined LANCZOS
          WRITE (out,100) GradErr, 'GradErr',                           &
     &            'Upper bound on relative error of the gradient.'
          WRITE (out,100) HevecErr, 'HevecErr',                         &
     &            'Accuracy required for Hessian eigenvectors.'
          WRITE (out,70) LhessianEV, 'LhessianEV',                      &
     &            'Switch to compute Hessian eigenvectors.'
#   endif
          WRITE (out,70) Lprecond, 'Lprecond',                          &
     &           'Switch to preconditioning with Hessian eigenpairs.'
#   ifndef LANCZOS
          IF (Lprecond.and.(NritzEV.gt.0)) THEN
            WRITE (out,80) NritzEV, 'NritzEV',                          &
     &           'Number of preconditioning eigenpairs to use.'
          END IF
#   endif
#   if !(defined GRADIENT_CHECK || defined TLM_CHECK)
          WRITE (out,70) LdefNRM(ng), 'LdefNRM',                        &
     &            'Switch to create a normalization NetCDF file.'
          WRITE (out,70) LwrtNRM(ng), 'LwrtNRM',                        &
     &            'Switch to write out normalization factors.'
          IF (LwrtNRM(ng)) THEN
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
          IF (LwrtNRM(ng)) THEN
            DO i=1,NstateVar(ng)
              Cnorm(i)=.TRUE.         ! force computation of all factors
            END DO
          END IF
          IF (LwrtNRM(ng)) THEN
            WRITE (out,70) Cnorm(isFsur), 'Cnorm(isFsur)',              &
     &            'Compute 2D normalization factors at RHO-points.'
            WRITE (out,70) Cnorm(isUbar), 'Cnorm(isUbar)',              &
     &            'Compute 2D normalization factors at U-points.'
            WRITE (out,70) Cnorm(isVbar), 'Cnorm(isVbar)',              &
     &            'Compute 2D normalization factors at V-points.'
#    ifdef ADJUST_WSTRESS
            WRITE (out,70) Cnorm(isUstr), 'Cnorm(isUstr)',              &
     &            'Compute normalization factors at surface U-stress.'
            WRITE (out,70) Cnorm(isVstr), 'Cnorm(isVstr)',              &
     &            'Compute normalization factors at surface V-stress.'
#    endif
#    ifdef SOLVE3D
            WRITE (out,70) Cnorm(isUvel), 'Cnorm(isUvel)',              &
     &            'Compute 3D normalization factors at U-points.'
            WRITE (out,70) Cnorm(isVvel), 'Cnorm(isVvel)',              &
     &            'Compute 3D normalization factors at V-points.'
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(isTvar(itrc)), 'Cnorm(isTvar)',     &
     &            'Compute normalization factors for tracer ',itrc,     &
     &            TRIM(Vname(1,idTvar(itrc)))
            END DO
#     ifdef ADJUST_STFLUX
            DO itrc=1,NT(ng)
              WRITE (out,110) Cnorm(isTsur(itrc)), 'Cnorm(isTsur)',     &
     &            'Compute normalization factors for flux of tracer ',  &
     &            itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
#     endif
#    endif
          END IF
          WRITE (out,100) Hgamma, 'Hgamma',                             &
     &            'Horizontal diffusion stability/accuracy factor.'
#    ifdef SOLVE3D
          WRITE (out,100) Vgamma, 'Vgamma',                             &
     &            'Vertical diffusion stability/accuracy factor.'
#    endif
          WRITE (out,120) Hdecay(isFsur,ng), 'Hdecay(isFsur)',          &
     &            'Horizontal decorrelation scale (m), free-surface.'
          WRITE (out,120) Hdecay(isUbar,ng), 'Hdecay(isUbar)',          &
     &            'Horizontal decorrelation scale (m), 2D U-momentum.'
          WRITE (out,120) Hdecay(isVbar,ng), 'Hdecay(isVbar)',          &
     &            'Horizontal decorrelation scale (m), 2D V-momentum.'
#    ifdef ADJUST_WSTRESS
          WRITE (out,120) Hdecay(isUstr,ng), 'Hdecay(isUstr)',          &
     &            'Horizontal decorrelation scale (m), U-stress.'
          WRITE (out,120) Hdecay(isVstr,ng), 'Hdecay(isVstr)',          &
     &            'Horizontal decorrelation scale (m), V-stress.'
#    endif
#    ifdef SOLVE3D
          WRITE (out,120) Hdecay(isUvel,ng), 'Hdecay(isUvel)',          &
     &            'Horizontal decorrelation scale (m), 3D U-momentum.'
          WRITE (out,120) Hdecay(isVvel,ng), 'Hdecay(isVvel)',          &
     &            'Horizontal decorrelation scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(isTvar(itrc),ng), 'Hdecay(idTvar)',  &
     &            'Horizontal decorrelation scale (m), ',               &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#     ifdef ADJUST_STFLUX
          DO itrc=1,NT(ng)
            WRITE (out,130) Hdecay(isTsur(itrc),ng), 'Hdecay(idTsur)',  &
     &            'Horizontal decorrelation scale (m), ',               &
     &            TRIM(Vname(1,idTsur(itrc)))
          END DO
#     endif
          WRITE (out,120) Vdecay(isUvel,ng), 'Vdecay(isUvel)',          &
     &            'Vertical decorrelation scale (m), 3D U-momentum.'
          WRITE (out,120) Vdecay(isVvel,ng), 'Vdecay(isVvel)',          &
     &            'Vertical decorrelation scale (m), 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Vdecay(isTvar(itrc),ng), 'Vdecay(idTvar)',  &
     &            'Vertical decorrelation scale (m), ',                 &
     &            TRIM(Vname(1,idTvar(itrc)))
          END DO
#    endif
#   endif
          WRITE (out,120) Cfscale(isFsur,ng), 'Cfscale(isFsur)',        &
     &            'Forcing error scale factor, free-surface.'
          WRITE (out,120) Cfscale(isUbar,ng), 'Cfscale(isUbar)',        &
     &            'Forcing error scale factor, 2D U-momentum.'
          WRITE (out,120) Cfscale(isVbar,ng), 'Cfscale(isVbar)',        &
     &            'Forcing error scale factor, 2D V-momentum.'
#    ifdef ADJUST_WSTRESS
          WRITE (out,120) Cfscale(isUstr,ng), 'Cfscale(isUstr)',        &
     &            'Forcing error scale factor, surface U-stress.'
          WRITE (out,120) Cfscale(isVstr,ng), 'Cfscale(isVstr)',        &
     &            'Forcing error scale factor, surface V-stress.'
#    endif
#    ifdef SOLVE3D
          WRITE (out,120) Cfscale(isUvel,ng), 'Cfscale(isUvel)',        &
     &            'Forcing error scale factor, 3D U-momentum.'
          WRITE (out,120) Cfscale(isVvel,ng), 'Cfscale(isVvel)',        &
     &            'Forcing error scale factor, 3D V-momentum.'
          DO itrc=1,NT(ng)
            WRITE (out,130) Cfscale(isTvar(itrc),ng),                   &
     &              'Cfscale(idTvar)',                                  &
     &              'Forcing error scale factor, ',                     &
     &              TRIM(Vname(1,idTvar(itrc)))
          END DO
#     ifdef ADJUST_STFLUX
          DO itrc=1,NT(ng)
            WRITE (out,130) Cfscale(isTsur(itrc),ng),                   &
     &              'Cfscale(idTsur)',                                  &
     &              'Forcing error scale factor flux of tracer, ',      &
     &              TRIM(Vname(1,idTsur(itrc)))
          END DO
#     endif
#    endif
#   if defined SOLVE3D && defined N2NORM_PROFILE
          WRITE (out,80) npN2, 'npN2',                                  &
     &            'Order of polynomial for N2(z) profile.'
          IF (npN2.gt.0) THEN
            WRITE (out,*)
            DO k=0,npN2
              WRITE (out,140) k, pcoef_N2(k)
            END DO
          END IF
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
#   if defined IS4DVAR || defined IS4DVAR_OLD || defined IOM || \
       defined OPT_OBSERVATIONS || defined WEAK_CONSTRAINT
          fname=STDname(ng)
          INQUIRE (FILE=TRIM(fname), EXIST=inhere)
          IF (.not.inhere) GO TO 30
          WRITE (out,160) '      Standard Deviation File:  ',           &
     &                    TRIM(fname)
#   endif
#   if !(defined OBS_SENSITIVITY || S4DVAR)
          fname=NRMname(ng)
          WRITE (out,160) 'Correlation Norm Factors File:  ',           &
     &                    TRIM(fname)
          IF (.not.LdefNRM(ng)) THEN
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
#  ifdef OBS_SENSITIVITY
          WRITE (out,160) '    Input/Output Lanczos File:  ',           &
     &                    TRIM(HSSname(ng))
#  else
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
          GO TO 40
  30      WRITE (out,170) TRIM(fname)
          STOP
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
 140  FORMAT (1x,'N2(z) polynomial coefficient: ',i2.2,1x,1p,e15.8)
 150  FORMAT (/,' Input Assimilation Files:',/)
 160  FORMAT (2x,a,a)
 170  FORMAT (/,' READ_ASSPAR - could not find input file:  ',a)
 180  FORMAT (/,' READ_ASSPAR - variable info not yet loaded, ', a)
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
            Npts=Nfloats(1)
            IF (Ngrids.gt.1) Npts=MAXVAL(Nfloats)
            allocate ( Fcoor (Npts,Ngrids) )
            allocate ( Fcount(Npts,Ngrids) )
            allocate ( Ftype (Npts,Ngrids) )
            allocate ( Ft0(Npts,Ngrids) )
            allocate ( Fx0(Npts,Ngrids) )
            allocate ( Fy0(Npts,Ngrids) )
            allocate ( Fz0(Npts,Ngrids) )
            allocate ( Fdt(Npts,Ngrids) )
            allocate ( Fdx(Npts,Ngrids) )
            allocate ( Fdy(Npts,Ngrids) )
            allocate ( Fdz(Npts,Ngrids) )
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
                STOP
              END IF
              ncount(igrid)=ncount(igrid)+Fcount(nentry(igrid)+1,igrid)
              nentry(igrid)=nentry(igrid)+1        
            END DO
          END IF
        END IF
      END DO
  20  WRITE (out,50) line
      STOP
  30  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (ncount(ng).ne.Nfloats(ng)) THEN
            WRITE (stdout,60) ncount(ng), Nfloats(ng)
            STOP
          END IF
          IF (Lfloats(ng)) THEN
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
        DO i=1,nentry(ng)
          IF (Fcount(i,ng).eq.1) THEN
            nc=nc+1
            FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
            FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
            FLT(ng)%Ftype(nc)=Ftype(i,ng)
            IF (Fcoor(i,ng).eq.0) THEN
              FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),        &
     &                                REAL(Lm(ng),r8)+0.5_r8)
              FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),        &
     &                                REAL(Mm(ng),r8)+0.5_r8)
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
                FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+              &
     &                                   REAL(j-1,r8)*Fdt(i,ng))*       &
     &                                  day2sec
                FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                FLT(ng)%Ftype(nc)=Ftype(i,ng)
                IF (Fcoor(i,ng).eq.0) THEN
                  FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),    &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),    &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
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
                    FLT(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,zfloat),     &
     &                                      REAL(N(ng),r8))
                  ELSE
                    FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+                  &
     &                                      REAL(j-1,r8)*Fdz(i,ng)
                  END IF
                END IF
                FLT(ng)%Ftype(nc)=Ftype(i,ng)
                IF (Fcoor(i,ng).eq.0) THEN
                  xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,xfloat),       &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,yfloat),       &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
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
      END DO
!
!  Deallocate local arrays.
!
      deallocate ( Fcoor )
      deallocate ( Fcount )
      deallocate ( Ftype )
      deallocate ( Ft0 )
      deallocate ( Fx0 )
      deallocate ( Fy0 )
      deallocate ( Fz0 )
      deallocate ( Fdt )
      deallocate ( Fdx )
      deallocate ( Fdy )
      deallocate ( Fdz )
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
#if defined SEDIMENT && defined SED_MORPH
          ELSE IF (TRIM(KeyWord).eq.'Sout(idbath)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idbath,1))
#endif
          ELSE IF (TRIM(KeyWord).eq.'Sout(idFsur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idFsur,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idOvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idOvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTvar)') THEN
            Npts=load_l(Nval, Cval, MT*Ngrids, Lsta)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Sout(idTvar(itrc),ng)=Lsta(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUsms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVsms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbrs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbrs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbrs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbrs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbws)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbws,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbws)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbws,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbcs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbcs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbcs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbcs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbot,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbot,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbur,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbvr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbvr,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW2xx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW2xx,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW2xy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW2xy,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW2yy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW2yy,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU2rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU2rs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV2rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV2rs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU2Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU2Sd,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV2Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV2Sd,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3xx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3xx,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3xy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3xy,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3yy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3yy,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3zx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3zx,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idW3zy)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idW3zy,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU3rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU3rs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV3rs)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV3rs,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idU3Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idU3Sd,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idV3Sd)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idV3Sd,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWamp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWamp,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWlen)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWlen,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWdir)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWdir,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NAT*Ngrids, Lsta)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                Sout(idTsur(itrc),ng)=Lsta(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLhea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLhea,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idShea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idShea,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLrad,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSrad,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idevap)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idevap,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idrain)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idrain,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idDano)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idDano,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVvis)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVvis,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTdif,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSdif,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHsbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHsbl,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHbbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHbbl,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idMtke)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idMtke,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idMtls)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idMtls,1))
# if defined SOLVE3D && (defined BBL_MODEL || defined SEDIMENT)
          ELSE IF (TRIM(KeyWord).eq.'Sout(idBott)') THEN
            IF (MAXVAL(idBott).eq.0) THEN
              WRITE (out,110) 'idBott'
              STOP
            END IF
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
                WRITE (out,40) Nstation(ng), is(ng)
                STOP
              END IF
            END DO
          END IF
        END IF
      END DO
 20   WRITE (out,50) line
      STOP
 30   CONTINUE
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
# if defined NEARSHORE_MELLOR
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
            IF (Sout(idU2Sd,ng)) WRITE (out,80) Sout(idU2Sd,ng),        &
     &         'Sout(idU2Sd)',                                          &
     &         'Write out 2D u-momentum stokes velocity.'
            IF (Sout(idV2Sd,ng)) WRITE (out,80) Sout(idV2Sd,ng),        &
     &         'Sout(idV2Sd)',                                          &
     &         'Write out 2D v-momentum stokes velocity.'
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
            IF (Sout(idU3Sd,ng)) WRITE (out,80) Sout(idU3Sd,ng),        &
     &         'Sout(idU3Sd)',                                          &
     &         'Write out 3D u-momentum stokes velocity.'
            IF (Sout(idV3Sd,ng)) WRITE (out,80) Sout(idV3Sd,ng),        &
     &         'Sout(idV3Sd)',                                          &
     &         'Write out 3D v-momentum stokes velocity.'
#  endif
# endif
            IF (Sout(idWamp,ng)) WRITE (out,80) Sout(idWamp,ng),        &
     &         'Sout(idWamp)',                                          &
     &         'Write out wave height.'
            IF (Sout(idWlen,ng)) WRITE (out,80) Sout(idWlen,ng),        &
     &         'Sout(idWlen)',                                          &
     &         'Write out wave length.'
            IF (Sout(idWdir,ng)) WRITE (out,80) Sout(idWdir,ng),        &
     &         'Sout(idWdir)',                                          &
     &         'Write out wave direction.'
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
!  Get length of "line", remove leading and trailing blanks.
!
      line=TRIM(ADJUSTL(line_text))
      Linp=LEN_TRIM(line)
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
! Extract parameter values string.  Remove comments [char(33)=!] or
! continuation symbol [char(92)=\], if any.  Trim leading trailing
! blanks.
!
        IF (Kextract) THEN
          Icmm=INDEX(line,CHAR(33),BACK=.FALSE.)
          IF (Icmm.gt.0) Lend=Icmm-1
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
