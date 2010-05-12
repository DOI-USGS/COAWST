1) Computation with unstructured grids is parallelized using OpenMP directives.
2) For writing block output to a Matlab binary file, the old format (Level 4)
   is replaced by the new one (Level 5 MAT-file format).

--- Makefile	2009-08-24 13:08:38.000000000 +0200
+++ Makefile	2009-08-24 13:06:37.000000000 +0200
@@ -76,6 +76,7 @@
 SwanInterpolateOutput.$(EXTO) \
 SwanConvAccur.$(EXTO) \
 SwanConvStopc.$(EXTO) \
+SwanThreadBounds.$(EXTO) \
 SwanFindObstacles.$(EXTO) \
 SwanCrossObstacle.$(EXTO) \
 SwanComputeForce.$(EXTO) \
@@ -87,7 +88,7 @@
 HCAT_EXE = hcat.exe
 HCAT_OBJS = swanhcat.$(EXTO)
 
-.SUFFIXES: .f .F .for .f90 .F90
+.SUFFIXES: .f .for .f90
 
 .PHONEY: help
 
@@ -115,7 +116,7 @@
 	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" $(HCAT_EXE)
 
 omp:
-	@perl switch.pl $(swch) -omp *.ftn *.ftn90
+	@perl switch.pl $(swch) *.ftn *.ftn90
 	$(MAKE) FOR=$(F90_OMP) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_OMP)" \
                 INCS="$(INCS_OMP)" LIBS="$(LIBS_OMP)" $(SWAN_EXE)
 	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" $(HCAT_EXE)
@@ -145,12 +146,6 @@
 .f90.o:
 	$(FOR) $< -c $(FFLAGS) $(INCS)
 
-.F.o:
-	$(FOR) $< -c $(FFLAGS) $(INCS)
-
-.F90.o:
-	$(FOR) $< -c $(FFLAGS) $(INCS)
-
 .for.o:
 	$(FOR) $< -c $(FFLAGS) $(INCS)
 
@@ -164,10 +159,10 @@
 	$(RM) *.$(EXTO) *.mod
 
 clobber:
-	$(RM) *.$(EXTO) *.mod *.f *.F *.F90 *.for *.f90 $(SWAN_EXE) $(HCAT_EXE)
+	$(RM) *.$(EXTO) *.mod *.f *.for *.f90 $(SWAN_EXE) $(HCAT_EXE)
 
 allclean:
-	$(RM) *.$(EXTO) *.mod *.f *.F *.F90 *.for *.f90 $(SWAN_EXE) $(HCAT_EXE)
+	$(RM) *.$(EXTO) *.mod *.f *.for *.f90 $(SWAN_EXE) $(HCAT_EXE)
 
 cleandoc:
 	$(MAKE) -f Makefile.latex TARGET=swanuse cleandoc
--- ocpmix.ftn	2009-08-24 13:08:38.000000000 +0200
+++ ocpmix.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -1461,6 +1461,8 @@
 !
 !  8. SUBROUTINE USED
 !
+!$    LOGICAL,EXTERNAL :: OMP_IN_PARALLEL
+!
 !  9. SUBROUTINES CALLING
 !
 ! 10. ERROR MESSAGES
@@ -1481,6 +1483,7 @@
         STPNOW = .FALSE.
       END IF
       IF (MAXERR.EQ.-1) STPNOW = .FALSE.
+!$    IF (OMP_IN_PARALLEL()) STPNOW = .FALSE.
 !
       RETURN
       END
--- platform.pl	2009-08-24 13:08:38.000000000 +0200
+++ platform.pl	2009-08-24 13:06:37.000000000 +0200
@@ -154,18 +154,18 @@
     print OUTFILE "F90_SER = ifort\n";
     print OUTFILE "F90_OMP = ifort\n";
     print OUTFILE "F90_MPI = mpif90\n";
-    if ($cpu =~ /i686/i) {
-      print OUTFILE "FLAGS_OPT = -O2 -xN -mp1\n";
-    }
-    elsif ($cpu =~ /x86_64/i) {
-      print OUTFILE "FLAGS_OPT = -O2 -ipo -xW -mp1\n";
-    }
-    else {
+#    if ($cpu =~ /i686/i) {
+#      print OUTFILE "FLAGS_OPT = -O2 -xN -mp1\n";
+#    }
+#    elsif ($cpu =~ /x86_64/i) {
+#      print OUTFILE "FLAGS_OPT = -O2 -ipo -xW -mp1\n";
+#    }
+#    else {
       print OUTFILE "FLAGS_OPT = -O2\n";
-    }
+#    }
     print OUTFILE "FLAGS_MSC = -W0 -assume byterecl -traceback\n";
     print OUTFILE "FLAGS_SER =\n";
-    print OUTFILE "FLAGS_OMP = -openmp -assume cc_omp -fpp2\n";
+    print OUTFILE "FLAGS_OMP = -openmp\n";
     print OUTFILE "FLAGS_MPI =\n";
     print OUTFILE "INCS_SER =\n";
     print OUTFILE "INCS_OMP =\n";
@@ -192,7 +192,7 @@
     print OUTFILE "FLAGS_OPT = -O2 -tpp7\n";
     print OUTFILE "FLAGS_MSC = -W0 -auto\n";
     print OUTFILE "FLAGS_SER =\n";
-    print OUTFILE "FLAGS_OMP = -openmp -fpp2\n";
+    print OUTFILE "FLAGS_OMP = -openmp\n";
     print OUTFILE "FLAGS_MPI =\n";
     print OUTFILE "INCS_SER =\n";
     print OUTFILE "INCS_OMP =\n";
@@ -219,7 +219,7 @@
     print OUTFILE "FLAGS_OPT = -O2 -tpp1\n";
     print OUTFILE "FLAGS_MSC = -W0 -auto\n";
     print OUTFILE "FLAGS_SER =\n";
-    print OUTFILE "FLAGS_OMP = -openmp -fpp2\n";
+    print OUTFILE "FLAGS_OMP = -openmp\n";
     print OUTFILE "FLAGS_MPI =\n";
     print OUTFILE "INCS_SER =\n";
     print OUTFILE "INCS_OMP =\n";
@@ -238,27 +238,26 @@
     system 'rm pgf90';
     print OUTFILE "##############################################################################\n";
     print OUTFILE "# IA32_PGF:		Intel Pentium with Linux using Portland Group compiler\n";
-    print OUTFILE "# Note: -lpgc is required for OpenMP.\n";
     print OUTFILE "##############################################################################\n";
     print OUTFILE "F90_SER = pgf90\n";
     print OUTFILE "F90_OMP = pgf90\n";
     print OUTFILE "F90_MPI = mpif90\n";
-    print OUTFILE "FLAGS_OPT =\n";
-    print OUTFILE "FLAGS_MSC = -mcmodel=medium\n";
-    print OUTFILE "FLAGS_SER = -fast\n";
+    print OUTFILE "FLAGS_OPT = -fast\n";
+    print OUTFILE "FLAGS_MSC =\n";
+    print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP = -mp\n";
-    print OUTFILE "FLAGS_MPI = -fast\n";
+    print OUTFILE "FLAGS_MPI =\n";
     print OUTFILE "INCS_SER =\n";
     print OUTFILE "INCS_OMP =\n";
     print OUTFILE "INCS_MPI =\n";
     print OUTFILE "LIBS_SER =\n";
-    print OUTFILE "LIBS_OMP = -lpgc\n";
-    print OUTFILE "LIBS_MPI = -lmpi\n";
+    print OUTFILE "LIBS_OMP =\n";
+    print OUTFILE "LIBS_MPI =\n";
     print OUTFILE "OUT = -o \n";
     print OUTFILE "EXTO = o\n";
     print OUTFILE "MAKE = make\n";
     print OUTFILE "RM = rm -f\n";
-    print OUTFILE "swch = -unix -impi -timg\n";
+    print OUTFILE "swch = -unix -f95 -impi -timg\n";
   }
   elsif ( -f "lf95" )
   {
@@ -298,7 +297,7 @@
     print OUTFILE "FLAGS_OPT = -O\n";
     print OUTFILE "FLAGS_MSC = -w -ffree-line-length-none\n";
     print OUTFILE "FLAGS_SER =\n";
-    print OUTFILE "FLAGS_OMP =\n";
+    print OUTFILE "FLAGS_OMP = -fopenmp\n";
     print OUTFILE "FLAGS_MPI =\n";
     print OUTFILE "INCS_SER =\n";
     print OUTFILE "INCS_OMP =\n";
@@ -353,7 +352,7 @@
   print OUTFILE "FLAGS_OPT = /optimize:2\n";
   print OUTFILE "FLAGS_MSC = /assume:byterecl /traceback /nowarn /nologo\n";
   print OUTFILE "FLAGS_SER =\n";
-  print OUTFILE "FLAGS_OMP = /Qopenmp /assume:cc_omp /fpp\n";
+  print OUTFILE "FLAGS_OMP = /Qopenmp\n";
   print OUTFILE "FLAGS_MPI =\n";
   print OUTFILE "INCS_SER =\n";
   print OUTFILE "INCS_OMP =\n";
@@ -365,7 +364,7 @@
   print OUTFILE "EXTO = obj\n";
   print OUTFILE "MAKE = nmake\n";
   print OUTFILE "RM = del\n";
-  print OUTFILE "swch = -dos -f95 -impi -cvis -timg\n";
+  print OUTFILE "swch = -dos -impi -cvis\n";
 }
 elsif ($os =~ /Darwin/i) {
   print OUTFILE "##############################################################################\n";
--- swancom1.ftn	2009-08-24 13:08:40.000000000 +0200
+++ swancom1.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -1008,7 +1008,6 @@
 !$OMP+PRIVATE(IP,IDC,ISC)                                                 40.31
 !$OMP+PRIVATE(I1GRD,I2GRD,I1MYC,I2MYC)                                    40.31
 !$OMP+PRIVATE(JDUM,JJ,III)                                                40.31
-!$OMP+PRIVATE(IXGRID,IXQUAD,IQERR,J1,J2)                                  40.41
 !$OMP+PRIVATE(IS,IE,INCI,JS,JE,INCJ,JSD,JED,JNODE,JWFRS,JWFRE)            40.31
 !$OMP+COPYIN(ICMAX,CSETUP)
 !$OMP+COPYIN(COSLAT,PROPSL)
@@ -1016,7 +1015,10 @@
 !
 !$OMP MASTER                                                              40.22
 !  Print number of threads set by environment                             40.22
-!$    WRITE(SCREEN,*)"Number of OpenMP threads=",OMP_GET_NUM_THREADS()    40.22
+!$    IF ( IT.EQ.1 )                                                      41.10
+!$   +WRITE(SCREEN,'(a,i2/)')                                             41.10 40.22
+!$   +      ' Number of threads during execution of parallel region = ',  41.10 40.22
+!$   +      OMP_GET_NUM_THREADS()                                         41.10 40.22
 !$OMP END MASTER                                                          40.22
 !TIMG
 !TIMG      CALL SWTSTA(101)                                                    40.23
@@ -1173,6 +1175,7 @@
 !     *** interaction files                                    ***        40.41
 !                                                                         40.41
 !TIMG      CALL SWTSTA(135)                                                    40.41
+!$OMP MASTER
       IF (IQUAD.EQ.51.OR.IQUAD.EQ.52.OR.IQUAD.EQ.53) THEN                 40.41
          CALL init_constants                                              40.41
          IXQUAD = IQUAD - 50                                              40.41
@@ -1186,13 +1189,11 @@
          END IF                                                           40.41
 !                                                                         40.41
          IXGRID = 3                                                       40.41
-         J1 = I1GRD                                                       40.41
-         IF (J1.EQ.1) J1 = 2                                              40.41
-         J2 = I2GRD                                                       40.41
          CALL xnl_init(SPCSIG  , SPCDIR(:,1)*180./PI , MSC , MDC ,        40.41
-     &                 -PWTAIL(1), GRAV   , COMPDA(J1:J2,JDP2) ,          40.41
-     &                 J2-J1+1   , IXQUAD , IXGRID   ,INODE,IQERR )       40.41
+     &                 -PWTAIL(1), GRAV   , COMPDA(2,JDP2) ,              40.41
+     &                 MCGRD-1   , IXQUAD , IXGRID   ,INODE,IQERR )       40.41
       END IF                                                              40.41
+!$OMP END MASTER
 !TIMG      CALL SWTSTO(135)                                                    40.41
 !
 !TIMG      CALL SWTSTA(103)                                                    40.23
@@ -1655,11 +1656,7 @@
      &            CAX1             ,CAY1                                  40.22
      &                                                                 )
 !TIMG                CALL SWTSTO(104)                                          40.23
-!               This step is disallowed within a parallel region in       40.22
-!               OpenMP.                                                   40.22
-!OMP#ifndef _OPENMP                                                           40.22
                 IF (STPNOW()) RETURN                                      34.01
-!OMP#endif                                                                    40.22
 
 !----------------------------------------------------------------------   40.22
 !               Once the computation is done for grid point (IX,IY) the   40.22
@@ -1830,15 +1827,10 @@
 !       *** compute wave-induced setup ***                                32.02
 !
 !TIMG        CALL SWTSTA(106)                                                  40.23
-        IF (LSETUP.GT.0) THEN                                             31.03
-          CALL SETUPP ( KGRPNT, MSTPDA, SETPDA, AC2, COMPDA(1,JDP2),      40.41
+        IF (LSETUP.GT.0)                                                  31.03
+     &    CALL SETUPP ( KGRPNT, MSTPDA, SETPDA, AC2, COMPDA(1,JDP2),      40.41
      &                  COMPDA(1,JDPSAV), COMPDA(1,JSETUP),               40.41
      &                  XCGRID, YCGRID, SPCSIG, SPCDIR )                  40.41
-!         This step is disallowed within a parallel region in OpenMP.     40.22
-!OMP#ifndef _OPENMP                                                           40.22
-          IF (STPNOW()) RETURN                                            34.01
-!OMP#endif                                                                    40.22
-        ENDIF                                                             31.03
 !TIMG        CALL SWTSTO(106)                                                  40.23
 !
 !----------------------------------------------------------------------   40.31
@@ -3085,25 +3077,26 @@
 !       leak and radiation stress in present grid point
 !
 !TIMG        CALL SWTSTA(124)                                                  40.23
-        CALL ADDDIS (COMPDA(1,JDISS)    ,COMPDA(1,JLEAK)    ,
-     &               AC2                ,LSWMAT(1,1,JABIN)  ,             40.02
-     &               SWMATR(1,1,JDIS0)  ,SWMATR(1,1,JDIS1)  ,
-     &               SWMATR(1,1,JGEN0)  ,SWMATR(1,1,JGEN1)  ,             40.85
-     &               SWMATR(1,1,JRED0)  ,SWMATR(1,1,JRED1)  ,             40.85
-     &               SWMATR(1,1,JTRA0)  ,SWMATR(1,1,JTRA1)  ,             40.85
-     &               SWMATR(1,1,JMATL)  ,SWMATR(1,1,JMATU)  ,             40.85
-     &               SWMATR(1,1,JMAT5)  ,SWMATR(1,1,JMAT6)  ,             40.85
-     &               COMPDA(1,JDSXB)    ,                                 40.67 40.61
-     &               COMPDA(1,JDSXS)    ,                                 40.67 40.61
-     &               COMPDA(1,JDSXW)    ,                                 40.67 40.61
-     &               COMPDA(1,JGSXW)    ,COMPDA(1,JGENR)    ,             40.85
-     &               COMPDA(1,JRSXQ)    ,COMPDA(1,JRSXT)    ,             40.85
-     &               COMPDA(1,JREDS)    ,                                 40.85
-     &               COMPDA(1,JTSXG)    ,COMPDA(1,JTSXT)    ,             40.85
-     &               COMPDA(1,JTSXS)    ,COMPDA(1,JTRAN)    ,             40.85
-     &               SWMATR(1,1,JLEK1)  ,COMPDA(1,JRADS)    ,             40.85
-     &               SPCSIG                                               30.72
-     &                                                      )
+        IF ( LADDS )                                                      40.85
+     &     CALL ADDDIS (COMPDA(1,JDISS)    ,COMPDA(1,JLEAK)    ,
+     &                  AC2                ,LSWMAT(1,1,JABIN)  ,          40.02
+     &                  SWMATR(1,1,JDIS0)  ,SWMATR(1,1,JDIS1)  ,
+     &                  SWMATR(1,1,JGEN0)  ,SWMATR(1,1,JGEN1)  ,          40.85
+     &                  SWMATR(1,1,JRED0)  ,SWMATR(1,1,JRED1)  ,          40.85
+     &                  SWMATR(1,1,JTRA0)  ,SWMATR(1,1,JTRA1)  ,          40.85
+     &                  SWMATR(1,1,JMATL)  ,SWMATR(1,1,JMATU)  ,          40.85
+     &                  SWMATR(1,1,JMAT5)  ,SWMATR(1,1,JMAT6)  ,          40.85
+     &                  COMPDA(1,JDSXB)    ,                              40.67 40.61
+     &                  COMPDA(1,JDSXS)    ,                              40.67 40.61
+     &                  COMPDA(1,JDSXW)    ,                              40.67 40.61
+     &                  COMPDA(1,JGSXW)    ,COMPDA(1,JGENR)    ,          40.85
+     &                  COMPDA(1,JRSXQ)    ,COMPDA(1,JRSXT)    ,          40.85
+     &                  COMPDA(1,JREDS)    ,                              40.85
+     &                  COMPDA(1,JTSXG)    ,COMPDA(1,JTSXT)    ,          40.85
+     &                  COMPDA(1,JTSXS)    ,COMPDA(1,JTRAN)    ,          40.85
+     &                  SWMATR(1,1,JLEK1)  ,COMPDA(1,JRADS)    ,          40.85
+     &                  SPCSIG                                            30.72
+     &                                                         )
 !TIMG        CALL SWTSTO(124)                                                  40.23
 !
 !       *** if negative action density occur rescale with a factor ***
@@ -3861,9 +3854,7 @@
       ARR(2) = SMN2                                                       40.30
       CALL SWREDUCE( ARR, 2, SWREAL, SWSUM )                              40.30
       CALL SWREDUCE( NINDX, 1, SWINT, SWSUM )                             40.41
-!OMP#ifndef _OPENMP                                                           40.31
-      IF (STPNOW()) RETURN                                                40.30
-!OMP#endif                                                                    40.31
+!MPI      IF (STPNOW()) RETURN                                                40.30
       HSMN2 = ARR(1) / REAL(NINDX)                                        40.30 30.82
        SMN2 = ARR(2) / REAL(NINDX)                                        40.30 30.82
 !$OMP END MASTER                                                          40.31
@@ -3962,9 +3953,7 @@
 !$OMP MASTER                                                              40.31
 !
       CALL SWREDUCE ( IACCUR, 1, SWINT, SWSUM )                           40.30
-!OMP#ifndef _OPENMP                                                           40.31
-      IF (STPNOW()) RETURN                                                40.30
-!OMP#endif                                                                    40.31
+!MPI      IF (STPNOW()) RETURN                                                40.30
       ACCUR  = REAL(IACCUR) * 100. / REAL(NINDX)
 !$OMP END MASTER                                                          40.31
 !$OMP BARRIER                                                             40.31
@@ -8705,9 +8694,7 @@
       IARR(1) = IACCUR
       IARR(2) = WETGRD
       CALL SWREDUCE ( IARR, 2, SWINT, SWSUM )
-!OMP#ifndef _OPENMP
-      IF (STPNOW()) RETURN
-!OMP#endif
+!MPI      IF (STPNOW()) RETURN
       IACCUR = IARR(1)
       WETGRD = IARR(2)
       ACCUR  = REAL(IACCUR) * 100. / REAL(WETGRD)
--- SwanCompdata.ftn90	2009-08-24 13:08:38.000000000 +0200
+++ SwanCompdata.ftn90	2009-08-24 13:06:40.000000000 +0200
@@ -64,6 +64,7 @@
     integer, dimension(MICMAX)                 :: vs     ! computational stencil, i.e. set of vertices
                                                          ! needed for the computation of a new value
                                                          ! in the present vertex
+!$omp threadprivate(vs)
 !
     integer, dimension(:,:), save, allocatable :: blist  ! list of boundary vertices in ascending order for each boundary polygon
     integer, dimension(:)  , save, allocatable :: vlist  ! vertex list
--- SwanCompUnstruc.ftn90	2009-08-24 13:08:40.000000000 +0200
+++ SwanCompUnstruc.ftn90	2009-08-24 13:06:49.000000000 +0200
@@ -35,6 +35,7 @@
 !   40.85: Marcel Zijlema
 !   41.02: Marcel Zijlema
 !   41.07: Marcel Zijlema
+!   41.10: Marcel Zijlema
 !
 !   Updates
 !
@@ -42,6 +43,7 @@
 !   40.85,   August 2008: add propagation, generation and redistribution terms for output purposes
 !   41.02, February 2009: implementation of diffraction
 !   41.07,   August 2009: bug fix: never-ending sweep is prevented
+!   41.10,   August 2009: parallelization using OpenMP directives
 !
 !   Purpose
 !
@@ -109,6 +111,7 @@
     integer                               :: idwmin    ! minimum counter for spectral wind direction
     integer, save                         :: ient = 0  ! number of entries in this subroutine
     integer                               :: iface     ! face index
+    integer                               :: inocnt    ! inocnv counter for calling thread
     integer                               :: inocnv    ! integer indicating number of vertices in which solver does not converged
     integer                               :: isslow    ! minimum frequency that is propagated within a sweep
     integer                               :: isstop    ! maximum frequency that is propagated within a sweep
@@ -116,6 +119,8 @@
     integer                               :: istot     ! maximum number of bins in frequency space for considered sweep
     integer                               :: iter      ! iteration counter
     integer                               :: ivert     ! vertex index
+    integer                               :: ivlow     ! lower index in range of vertices in calling thread
+    integer                               :: ivup      ! upper index in range of vertices in calling thread
     integer                               :: j         ! loop counter
     integer                               :: jc        ! loop counter
     integer                               :: k         ! loop counter
@@ -132,8 +137,10 @@
     integer                               :: npfr      ! number of vertices in which rescaling is used
     integer                               :: nupdv     ! number of updates of present vertex (based on surrounding cells)
     integer                               :: nwetp     ! total number of active vertices
+    integer                               :: swpcol    ! own sweep color (0/1=black/white) in calling thread
     integer                               :: swpdir    ! sweep counter
     integer                               :: swpnr     ! sweep number
+    integer                               :: tid       ! thread number
     integer, dimension(3)                 :: v         ! vertices in present cell
     integer                               :: vb        ! vertex of begin of present face
     integer                               :: ve        ! vertex of end of present face
@@ -147,6 +154,7 @@
     integer, dimension(:), allocatable    :: islmin    ! lowest sigma-index occured in applying limiter
     integer, dimension(:), allocatable    :: nflim     ! number of frequency use of limiter in each vertex
     integer, dimension(:), allocatable    :: nrscal    ! number of frequency use of rescaling in each vertex
+!$  integer, dimension(:), allocatable    :: tlist     ! vertex list for calling thread
     !
     real                                  :: abrbot    ! near bottom excursion
     real                                  :: accur     ! percentage of active vertices in which required accuracy has been reached
@@ -232,6 +240,9 @@
     !
     type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
     type(verttype), dimension(:), pointer :: vert      ! datastructure for vertices with their attributes
+    !
+!$  integer, external                     :: omp_get_num_threads ! number of OpenMP threads being used
+!$  integer, external                     :: omp_get_thread_num  ! get thread number
 !
 !   Structure
 !
@@ -254,6 +265,10 @@
     IXCGRD(1) = -9999  ! to be used in routines SINTGRL and SOURCE so that Ursell number and
     IYCGRD(1) = -9999  ! quadruplets are calculated once in each vertex during an iteration
     !
+    tid   = 0
+    ivlow = 1
+    ivup  = nverts
+    !
     ! print all the settings used in SWAN run
     !
     if ( it == 1 .and. ITEST > 0 ) call SWPRSET (spcsig)
@@ -292,8 +307,51 @@
        allocate(memnl4(0,0,0))
     endif
     !
+    !$ allocate(tlist(nverts))
+!TIMG    call SWTSTO(101)
+    !
+    ! initialization of shared arrays
+    !
+    hscurr = 0.
+    hsprev = 0.
+    hsdifc = 0.
+    tmcurr = 0.
+    tmprev = 0.
+    tmdifc = 0.
+    !
+    swtsda = 0.
+    !
+    ! spawn a parallel region
+    !
+    !$omp parallel default(shared) &
+    !$omp private(cad, cas, cax, cay, cgo, kwave) &
+    !$omp private(idcmin, idcmax, iscmin, iscmax, anybin) &
+    !$omp private(amat, rhs, ac2old) &
+    !$omp private(anywnd, obredf, reflso, alimw, groww, anyblk) &
+    !$omp private(disc0, disc1, genc0, genc1, redc0, redc1, trac0, trac1, leakcf) &
+    !$omp private(ue, sa1, sa2, sfnl) &
+    !$omp private(da1c, da1p, da1m, da2c, da2p, da2m, dsnl) &
+    !$omp private(tid, ivlow, ivup, iter, swpcol, kvert, kstart, kend, kstep) &
+    !$omp private(ivert, nupdv, jc, k, j, icell, v, vu, l, swpnr, rdx, rdy, lpredt, vb, ve, iface, link, inocnt) &
+    !$omp private(iddlow, iddtop, idtot, isslow, isstop, istot) &
+    !$omp private(abrbot, kmespc, idwmin, idwmax, hs, etot, qbloc, ufric, fpm, thetaw, hm, wind10, smebrk) &
+    !$omp copyin(ICMAX, COSLAT, IPTST, TESTFL)
+    !
+    ! print number of threads set by environment
+    !
+    !$omp master
+    !$    if ( it == 1 ) &
+    !$    write (SCREEN,'(a,i2/)') ' Number of threads during execution of parallel region = ', omp_get_num_threads()
+    !$omp end master
+    !
+    ! get thread number
+    !
+    !$ tid = omp_get_thread_num()
+    tid = tid + 1
+    !
     ! allocation of private arrays
     !
+!TIMG    call SWTSTA(101)
     allocate(   cad(MDC,MSC      ))
     allocate(   cas(MDC,MSC      ))
     allocate(   cax(MDC,MSC,ICMAX))
@@ -330,9 +388,11 @@
     !
     ! calculate ranges of spectral space for arrays related to 4 wave-wave interactions
     !
+    !$omp single
 !TIMG    call SWTSTA(135)
     if ( IQUAD > 0 ) call FAC4WW ( xis, snlc1, dal1, dal2, dal3, spcsig, wwint, wwawg, wwswg )
 !TIMG    call SWTSTO(135)
+    !$omp end single
     !
     if ( IQUAD > 0 ) then
        allocate(  ue(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
@@ -369,21 +429,11 @@
        allocate(da2m(0,0))
        allocate(dsnl(0,0))
     endif
-    !
-    ! initialization of shared arrays
-    !
-    hscurr = 0.
-    hsprev = 0.
-    hsdifc = 0.
-    tmcurr = 0.
-    tmprev = 0.
-    tmdifc = 0.
-    !
-    swtsda = 0.
 !TIMG    call SWTSTO(101)
     !
     ! marks vertices active and non-active
     !
+    !$omp single
     nwetp = 0
     do kvert = 1, nverts
        if ( compda(kvert,JDP2) > DEPMIN ) then
@@ -403,12 +453,40 @@
     else
        fguess = .false.
     endif
+    !$omp end single
+    !
+    ! compute load-balanced spatial loop bounds for each thread
+    !
+    !$ call SwanThreadBounds( nwetp, ivlow, ivup, tlist )
     !
 !TIMG    call SWTSTA(103)
     iterloop: do iter = 1, ITERMX
        !
        ! some initializations
        !
+       if ( IQUAD  > 2 ) memnl4(1:MDC,1:MSC,ivlow:ivup) = 0.
+       if ( ITRIAD > 0 .or. ISURF == 5 ) compda(ivlow:ivup,JURSEL) = 0.
+       !
+       compda(ivlow:ivup,JDISS) = 0.
+       compda(ivlow:ivup,JLEAK) = 0.
+       compda(ivlow:ivup,JDSXB) = 0.
+       compda(ivlow:ivup,JDSXS) = 0.
+       compda(ivlow:ivup,JDSXW) = 0.
+       compda(ivlow:ivup,JGENR) = 0.
+       compda(ivlow:ivup,JGSXW) = 0.
+       compda(ivlow:ivup,JREDS) = 0.
+       compda(ivlow:ivup,JRSXQ) = 0.
+       compda(ivlow:ivup,JRSXT) = 0.
+       compda(ivlow:ivup,JTRAN) = 0.
+       compda(ivlow:ivup,JTSXG) = 0.
+       compda(ivlow:ivup,JTSXT) = 0.
+       compda(ivlow:ivup,JTSXS) = 0.
+       compda(ivlow:ivup,JRADS) = 0.
+       compda(ivlow:ivup,JQB  ) = 0.
+       !
+       inocnt = 0
+       !
+       !$omp master
        inocnv = 0
        !
        islmin = 9999
@@ -417,26 +495,6 @@
        !
        acnrms = -9999.
        !
-       if ( IQUAD  > 2 ) memnl4 = 0.
-       if ( ITRIAD > 0 .or. ISURF == 5 ) compda(:,JURSEL) = 0.
-       !
-       compda(:,JDISS) = 0.
-       compda(:,JLEAK) = 0.
-       compda(:,JDSXB) = 0.
-       compda(:,JDSXS) = 0.
-       compda(:,JDSXW) = 0.
-       compda(:,JGENR) = 0.
-       compda(:,JGSXW) = 0.
-       compda(:,JREDS) = 0.
-       compda(:,JRSXQ) = 0.
-       compda(:,JRSXT) = 0.
-       compda(:,JTRAN) = 0.
-       compda(:,JTSXG) = 0.
-       compda(:,JTSXT) = 0.
-       compda(:,JTSXS) = 0.
-       compda(:,JRADS) = 0.
-       compda(:,JQB  ) = 0.
-       !
        ! During first iteration, first guess of action density is based on 2nd generation mode
        ! After first iteration, user-defined settings are re-used
        !
@@ -521,6 +579,10 @@
        enddo
        !
        swpdir = 0
+       !$omp end master
+       !
+       ! synchronize threads before loop over sweep directions
+       !$omp barrier
        !
        ! loop over a number of sweeps through grid
        !
@@ -528,6 +590,7 @@
           !
           ! check whether all vertices in grid are completely updated (in both geographic and spectral spaces)
           !
+          !$omp single
           swpfull = .true.
           do kvert = 1, nverts
              if ( .not.vert(kvert)%fullupdated ) then
@@ -535,8 +598,11 @@
                 exit
              endif
           enddo
+          !$omp end single
+          !
           if ( swpfull ) exit sweeploop
           !
+          !$omp single
           swpdir = swpdir + 1
           !
           if ( SCREEN /= PRINTF ) then
@@ -551,22 +617,26 @@
                 endif
              endif
           endif
+          !$omp end single
           !
           ! loop over vertices in the grid (depending on sweep direction)
           !
-          if ( mod(swpdir,2) == 1 ) then
-             kstart = 1
-             kend   = nverts
+          swpcol = mod(swpdir,2)
+          !
+          if ( mod(tid,2) == swpcol ) then
+             kstart = ivlow
+             kend   = ivup
              kstep  = 1
           else
-             kstart = nverts
-             kend   = 1
+             kstart = ivup
+             kend   = ivlow
              kstep  = -1
           endif
           !
           vertloop: do kvert = kstart, kend, kstep
             !
             ivert = vlist(kvert)
+            !$ ivert = tlist(kvert)
             !
             if ( vert(ivert)%active ) then   ! this active vertex needs to be updated
                !
@@ -858,7 +928,7 @@
 !TIMG                              call SWTSTA(120)
                               call SWSIP ( ac2        , amat(1,1,1)    , rhs            , amat(1,1,4), &
                                            amat(1,1,5), amat(1,1,2)    , amat(1,1,3)    , ac2old     , &
-                                           PNUMS(12)  , nint(PNUMS(14)), nint(PNUMS(13)), inocnv     , &
+                                           PNUMS(12)  , nint(PNUMS(14)), nint(PNUMS(13)), inocnt     , &
                                            iddlow     , iddtop         , isstop         , idcmin     , &
                                            idcmax     )
 !TIMG                              call SWTSTO(120)
@@ -881,20 +951,21 @@
                         ! store propagation, generation, dissipation, redistribution, leak and radiation stress in present vertex
                         !
 !TIMG                        call SWTSTA(124)
-                        call ADDDIS ( compda(1,JDISS), compda(1,JLEAK), ac2   , anybin , &
-                                      disc0          , disc1          ,                  &
-                                      genc0          , genc1          ,                  &
-                                      redc0          , redc1          ,                  &
-                                      trac0          , trac1          ,                  &
-                                      amat(1,1,4)    , amat(1,1,5)    ,                  &
-                                      amat(1,1,2)    , amat(1,1,3)    ,                  &
-                                      compda(1,JDSXB), compda(1,JDSXS), compda(1,JDSXW), &
-                                      compda(1,JGSXW), compda(1,JGENR),                  &
-                                      compda(1,JRSXQ), compda(1,JRSXT),                  &
-                                      compda(1,JREDS),                                   &
-                                      compda(1,JTSXG), compda(1,JTSXT),                  &
-                                      compda(1,JTSXS), compda(1,JTRAN),                  &
-                                      leakcf         , compda(1,JRADS), spcsig           )
+                        if ( LADDS )                                                        &
+                           call ADDDIS ( compda(1,JDISS), compda(1,JLEAK), ac2   , anybin , &
+                                         disc0          , disc1          ,                  &
+                                         genc0          , genc1          ,                  &
+                                         redc0          , redc1          ,                  &
+                                         trac0          , trac1          ,                  &
+                                         amat(1,1,4)    , amat(1,1,5)    ,                  &
+                                         amat(1,1,2)    , amat(1,1,3)    ,                  &
+                                         compda(1,JDSXB), compda(1,JDSXS), compda(1,JDSXW), &
+                                         compda(1,JGSXW), compda(1,JGENR),                  &
+                                         compda(1,JRSXQ), compda(1,JRSXT),                  &
+                                         compda(1,JREDS),                                   &
+                                         compda(1,JTSXG), compda(1,JTSXT),                  &
+                                         compda(1,JTSXS), compda(1,JTRAN),                  &
+                                         leakcf         , compda(1,JRADS), spcsig           )
 !TIMG                        call SWTSTO(124)
                         !
                         ! if negative action density occur rescale with a factor
@@ -937,10 +1008,22 @@
             !
           enddo vertloop
           !
+          ! synchronize threads before starting next sweep
+          !$omp barrier
+          !
        enddo sweeploop
        !
+       ! global sum to inocnv counter
+       !
+       !$omp atomic
+       inocnv = inocnv + inocnt
+       !
+       ! synchronize threads before checking stop criterion
+       !$omp barrier
+       !
        ! store the source terms assembled in test points per iteration in the files IFPAR, IFS1D and IFS2D
        !
+       !$omp master
 !TIMG       call SWTSTA(105)
        if ( NPTST > 0 .and. NSTATM == 0 ) then
           if ( IFPAR > 0 ) write (IFPAR,151) iter
@@ -1004,6 +1087,7 @@
              write (PRINTF,136) inocnv
           endif
        endif
+       !$omp end master
        !
        ! info regarding the iteration process and the accuracy
        !
@@ -1012,14 +1096,15 @@
 !TIMG          call SWTSTA(102)
           if ( PNUMS(21) == 0. ) then
              !
-             call SwanConvAccur ( accur, hscurr, tmcurr, compda(1,JDHS), compda(1,JDTM), xytst, spcsig, ac2 )
+             call SwanConvAccur ( accur, hscurr, tmcurr, compda(1,JDHS), compda(1,JDTM), xytst, spcsig, ac2, ivlow, ivup )
              !
           else
              !
-             call SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, compda(1,JDHS), compda(1,JDTM), xytst, spcsig, ac2 )
+             call SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, compda(1,JDHS), compda(1,JDTM), xytst, spcsig, ac2, ivlow, ivup )
              !
           endif
           !
+          !$omp master
           if ( iter == 1 ) then
              accur = -9999.
              write (PRINTF,142)
@@ -1028,14 +1113,19 @@
              write (PRINTF,140) accur, PNUMS(4)
              if ( NSTATC == 0 ) write (SCREEN,140) accur, PNUMS(4)
           endif
+          !$omp end master
 !TIMG          call SWTSTO(102)
           !
+          ! synchronize threads
+          !$omp barrier
+          !
           ! if accuracy has been reached then terminates iteration process
           !
           if ( accur >= PNUMS(4) ) exit iterloop
           !
        elseif ( PNUMS(21) == 2. ) then
           !
+         !$omp master
 !TIMG          call SWTSTA(102)
           write (PRINTF,141)
           if ( NSTATC == 0 ) write (SCREEN,141)
@@ -1056,6 +1146,10 @@
           endif
           acnrmo = acnrms(1)
 !TIMG          call SWTSTO(102)
+          !$omp end master
+          !
+          ! synchronize threads
+          !$omp barrier
           !
           ! if accuracy has been reached then terminates iteration process
           !
@@ -1066,20 +1160,6 @@
     enddo iterloop
 !TIMG    call SWTSTO(103)
     !
-    ! store the source terms assembled in test points per time step in the files IFPAR, IFS1D and IFS2D
-    !
-!TIMG    call SWTSTA(105)
-    if ( NPTST > 0 .and. NSTATM == 1 ) then
-       if ( IFPAR > 0 ) write (IFPAR,152) CHTIME
-       if ( IFS1D > 0 ) write (IFS1D,152) CHTIME
-       if ( IFS2D > 0 ) write (IFS2D,152) CHTIME
-       call PLTSRC ( swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), swtsda(1,1,1,JPBTFR), &
-                     swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPTRI) , &
-                     ac2                 , spcsig              , compda(1,JDP2)      , xytst               , &
-                     dummy               )
-    endif
-!TIMG    call SWTSTO(105)
-    !
     ! deallocation of private arrays
     !
 !TIMG    call SWTSTA(101)
@@ -1128,9 +1208,29 @@
     deallocate(da2p)
     deallocate(da2m)
     deallocate(dsnl)
+!TIMG    call SWTSTO(101)
+    !
+    ! end of parallel region
+    !
+    !$omp end parallel
+    !
+    ! store the source terms assembled in test points per time step in the files IFPAR, IFS1D and IFS2D
+    !
+!TIMG    call SWTSTA(105)
+    if ( NPTST > 0 .and. NSTATM == 1 ) then
+       if ( IFPAR > 0 ) write (IFPAR,152) CHTIME
+       if ( IFS1D > 0 ) write (IFS1D,152) CHTIME
+       if ( IFS2D > 0 ) write (IFS2D,152) CHTIME
+       call PLTSRC ( swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), swtsda(1,1,1,JPBTFR), &
+                     swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPTRI) , &
+                     ac2                 , spcsig              , compda(1,JDP2)      , xytst               , &
+                     dummy               )
+    endif
+!TIMG    call SWTSTO(105)
     !
     ! deallocation of shared arrays
     !
+!TIMG    call SWTSTA(101)
     deallocate(islmin)
     deallocate( nflim)
     deallocate(nrscal)
@@ -1145,6 +1245,8 @@
     deallocate(swtsda)
     !
     deallocate(memnl4)
+    !
+!$  deallocate(tlist)
 !TIMG    call SWTSTO(101)
     !
  101 format (// ' Settings of 2nd generation mode as first guess are used:')
--- SwanConvAccur.ftn90	2009-08-24 13:08:40.000000000 +0200
+++ SwanConvAccur.ftn90	2009-08-24 13:06:41.000000000 +0200
@@ -1,4 +1,4 @@
-subroutine SwanConvAccur ( accur, hscurr, tmcurr, delhs, deltm, xytst, spcsig, ac2 )
+subroutine SwanConvAccur ( accur, hscurr, tmcurr, delhs, deltm, xytst, spcsig, ac2, ivlow, ivup )
 !
 !   --|-----------------------------------------------------------|--
 !     | Delft University of Technology                            |
@@ -32,10 +32,12 @@
 !   Authors
 !
 !   40.80: Marcel Zijlema
+!   41.10: Marcel Zijlema
 !
 !   Updates
 !
 !   40.80, October 2007: New subroutine
+!   41.10,  August 2009: parallelization using OpenMP directives
 !
 !   Purpose
 !
@@ -53,6 +55,8 @@
 !
 !   Argument variables
 !
+    integer, intent(in)                         :: ivlow  ! lower index in range of vertices in calling thread
+    integer, intent(in)                         :: ivup   ! upper index in range of vertices in calling thread
     integer, dimension(NPTST), intent(in)       :: xytst  ! test points for output purposes
     !
     real, intent(out)                           :: accur  ! percentage of active vertices in which required accuracy has been reached
@@ -71,7 +75,9 @@
     integer                               :: ivert    ! loop counter over vertices
     integer                               :: j        ! loop counter
     integer                               :: npacc    ! number of vertices in which required accuracy has been reached
+    integer                               :: npacct   ! npacc counter for calling thread
     integer                               :: nwetp    ! total number of active vertices
+    integer                               :: nwetpt   ! nwetp counter for calling thread
     !
     real                                  :: fact     ! auxiliary factor
     real                                  :: hsabs    ! absolute difference in wave height between last 2 iterations
@@ -79,9 +85,11 @@
     real                                  :: hsprev   ! wave height at previous iteration level
     real                                  :: hsrel    ! required accuracy with respect to relative error in wave height
     real                                  :: hsmean   ! space-averaged wave height at previous iteration level
+    real                                  :: hsmeant  ! hsmean variable for calling thread
     real                                  :: m0       ! moment of zeroth order
     real                                  :: m1       ! moment of first order
     real                                  :: tmmean   ! space-averaged mean period at previous iteration level
+    real                                  :: tmmeant  ! tmmean variable for calling thread
     real                                  :: tmabs    ! absolute difference in mean period between last 2 iterations
     real                                  :: tmoval   ! required accuracy with respect to overall error in mean period
     real                                  :: tmprev   ! mean period at previous iteration level
@@ -91,6 +99,9 @@
     logical                               :: tstfl    ! indicates whether vertex is a test point
     !
     type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
+!   Common variables
+!
+    common/convaccur/npacc,nwetp,hsmean,tmmean
 !
 !   Structure
 !
@@ -104,9 +115,17 @@
     !
     vert => gridobject%vert_grid
     !
+    !$omp single
     npacc = 0
     nwetp = 0
     !
+    hsmean = 0.
+    tmmean = 0.
+    !$omp end single
+    !
+    npacct = 0
+    nwetpt = 0
+    !
     deltm = 0.
     delhs = 0.
     !
@@ -114,28 +133,40 @@
     !
     ! calculate space-averaged wave height and mean period
     !
-    hsmean = 0.
-    tmmean = 0.
+    hsmeant = 0.
+    tmmeant = 0.
     !
-    do ivert = 1, nverts
+    do ivert = ivlow, ivup
        !
        if ( vert(ivert)%active ) then
           !
-          nwetp = nwetp + 1
+          nwetpt = nwetpt + 1
           !
-          hsmean = hsmean + hscurr(ivert)
-          tmmean = tmmean + tmcurr(ivert)
+          hsmeant = hsmeant + hscurr(ivert)
+          tmmeant = tmmeant + tmcurr(ivert)
           !
        endif
        !
     enddo
     !
+    ! global sum to nwetp, hsmean and tmmean
+    !
+    !$omp atomic
+    nwetp  = nwetp + nwetpt
+    !$omp atomic
+    hsmean = hsmean + hsmeant
+    !$omp atomic
+    tmmean = tmmean + tmmeant
+    !
+    !$omp barrier
+    !$omp single
     hsmean = hsmean/real(nwetp)
     tmmean = tmmean/real(nwetp)
+    !$omp end single
     !
     ! calculate a set of accuracy parameters based on relative and overall accuracy measures for Hs and Tm
     !
-    do ivert = 1, nverts
+    do ivert = ivlow, ivup
        !
        if ( vert(ivert)%active ) then
           !
@@ -197,7 +228,7 @@
           !
           ! count vertices where wave height and mean period have reached required accuracies
           !
-          if ( hsabs <= max(hsrel,hsoval) .and. tmabs <= max(tmrel,tmoval) ) npacc = npacc + 1
+          if ( hsabs <= max(hsrel,hsoval) .and. tmabs <= max(tmrel,tmoval) ) npacct = npacct + 1
           !
           if (tstfl) then
              if (lhead) write(PRINTF,11)
@@ -214,9 +245,17 @@
        !
     enddo
     !
+    ! global sum to npacc
+    !
+    !$omp atomic
+    npacc = npacc + npacct
+    !
     ! compute percentage of active vertices where required accuracy has been reached
     !
+    !$omp barrier
+    !$omp single
     accur = real(npacc)*100./real(nwetp)
+    !$omp end single
     !
  11 format(11x,'dHrel          ','dHoval         ','dTm01rel       ','dTm01oval ')
  12 format(1x,ss,'k=',i5,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
--- SwanConvStopc.ftn90	2009-08-24 13:08:40.000000000 +0200
+++ SwanConvStopc.ftn90	2009-08-24 13:06:41.000000000 +0200
@@ -1,4 +1,4 @@
-subroutine SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, delhs, deltm, xytst, spcsig, ac2 )
+subroutine SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, delhs, deltm, xytst, spcsig, ac2, ivlow, ivup )
 !
 !   --|-----------------------------------------------------------|--
 !     | Delft University of Technology                            |
@@ -33,11 +33,13 @@
 !
 !   40.80: Marcel Zijlema
 !   40.93: Marcel Zijlema
+!   41.10: Marcel Zijlema
 !
 !   Updates
 !
 !   40.80,   October 2007: New subroutine
 !   40.93, September 2008: extended with curvature of mean period
+!   41.10,    August 2009: parallelization using OpenMP directives
 !
 !   Purpose
 !
@@ -55,6 +57,8 @@
 !
 !   Argument variables
 !
+    integer, intent(in)                         :: ivlow  ! lower index in range of vertices in calling thread
+    integer, intent(in)                         :: ivup   ! upper index in range of vertices in calling thread
     integer, dimension(NPTST), intent(in)       :: xytst  ! test points for output purposes
     !
     real, intent(out)                           :: accur  ! percentage of active vertices in which required accuracy has been reached
@@ -77,7 +81,9 @@
     integer                               :: ivert    ! loop counter over vertices
     integer                               :: j        ! loop counter
     integer                               :: npacc    ! number of vertices in which required accuracy has been reached
+    integer                               :: npacct   ! npacc counter for calling thread
     integer                               :: nwetp    ! total number of active vertices
+    integer                               :: nwetpt   ! nwetp counter for calling thread
     !
     real                                  :: curvah   ! required accuracy with respect to curvature in wave height
     real                                  :: curvat   ! required accuracy with respect to curvature in mean period
@@ -100,6 +106,10 @@
     !
     type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
 !
+!   Common variables
+!
+    common/convstopc/npacc,nwetp
+!
 !   Structure
 !
 !   Description of the pseudo code
@@ -112,8 +122,13 @@
     !
     vert => gridobject%vert_grid
     !
+    !$omp single
     npacc = 0
     nwetp = 0
+    !$omp end single
+    !
+    npacct = 0
+    nwetpt = 0
     !
     deltm = 0.
     delhs = 0.
@@ -122,7 +137,7 @@
     !
     ! calculate a set of accuracy parameters based on relative error and curvature for Hs and Tm
     !
-    do ivert = 1, nverts
+    do ivert = ivlow, ivup
        !
        if ( vert(ivert)%active ) then
           !
@@ -138,7 +153,7 @@
           !
           ! count active points
           !
-          nwetp = nwetp + 1
+          nwetpt = nwetpt + 1
           !
           ! store wave height and mean period of previous iteration levels
           !
@@ -205,7 +220,7 @@
           ! count vertices where wave height and period have reached required accuracies
           !
           if ( (hscurv <= curvah .and. hsabs <= max(hsrel,PNUMS(2)) ) .and. &
-               (tmcurv <= curvat .and. tmabs <= max(tmrel,PNUMS(3)) ) ) npacc = npacc + 1
+               (tmcurv <= curvat .and. tmabs <= max(tmrel,PNUMS(3)) ) ) npacct = npacct + 1
           !
           if (tstfl) then
              if (lhead) write(PRINTF,11)
@@ -222,9 +237,19 @@
        !
     enddo
     !
+    ! global sum to npacc and nwetp
+    !
+    !$omp atomic
+    npacc = npacc + npacct
+    !$omp atomic
+    nwetp = nwetp + nwetpt
+    !
     ! compute percentage of active vertices where required accuracy has been reached
     !
+    !$omp barrier
+    !$omp single
     accur = real(npacc)*100./real(nwetp)
+    !$omp end single
     !
  11 format(11x,'dHabs          ','dHrel          ','Curvature H    ','dTabs          ','dTrel          ','Curvature T    ')
  12 format(1x,ss,'k=',i5,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
--- swan.edt	2009-08-24 13:08:40.000000000 +0200
+++ swan.edt	2009-08-24 13:06:37.000000000 +0200
@@ -205,18 +205,18 @@
 !                              (SPEC  [ndec])
 !
 !   BLOCK   'sname'  HEADER | NOHEADER  'fname' (LAY-OUT [idla])            &
-!          <  DSPR|HSIGN|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DEPTH|VEL|  &
-!             FRCOEFF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS| &
-!             DHSIGN|DRTM01|LEAK|TSEC|XP|YP|DIST|SETUP|TMM10|RTMM10|        &
+!          <  HSIGN|HSWELL|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DSPR|VEL| &
+!             FRCOEF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS|  &
+!             DHSIGN|DRTM01|LEAK|TSEC|XP|YP|DIST|SETUP|TMM10|RTMM10|DEPTH|  &
 !             TMBOT|QP|BFI|WATLEV|BOTLEV|TPS|DISBOT|DISSURF|DISWCAP|        &
 !             GENE|GENW|REDI|REDQ|REDT|PROPA|PROPX|PROPT|PROPS|RADS >       &
 !             ([unit]) (OUTPUT [tbegblk] [deltblk] SEC|MIN|HR|DAY)
 !
 !   TABLE   'sname'  HEADER | NOHEADER | INDEXED 'fname'                    &
-!          <  DSPR|HSIGN|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DEPTH|VEL|  &
-!             FRCOEFF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS| &
+!          <  HSIGN|HSWELL|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DSPR|VEL| &
+!             FRCOEF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS|  &
 !             DHSIGN|DRTM01|LEAK|TIME|TSEC|XP|YP|DIST|SETUP|TMM10|RTMM10|   &
-!             TMBOT|QP|BFI|WATLEV|BOTLEV|TPS|DISBOT|DISSURF|DISWCAP|        &
+!             DEPTH|TMBOT|QP|BFI|WATLEV|BOTLEV|TPS|DISBOT|DISSURF|DISWCAP|  &
 !             GENE|GENW|REDI|REDQ|REDT|PROPA|PROPX|PROPT|PROPS|RADS >       &
 !             ([unit]) (OUTPUT [tbegtbl] [delttbl] SEC|MIN|HR|DAY)
 !
--- swanimp.tex	2009-08-24 13:08:40.000000000 +0200
+++ swanimp.tex	2009-08-24 13:06:37.000000000 +0200
@@ -15,7 +15,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72ABC}
+{\Large\bf SWAN Cycle III version 40.72ABCD}
 \end{center}
 
 \cleardoublepage
@@ -71,7 +71,7 @@
 SWAN wave model on your computer system. The total material consists of:
 \begin{itemize}
     \item the SWAN source code,
-    \item the SWAN executables for MS Windows,
+    \item the SWAN executable for Microsoft Windows,
     \item the User Manual,
     \item this Implementation Manual,
     \item the Technical documentation,
@@ -103,6 +103,8 @@
 additional files, e.g. for testing SWAN. Please read the copyright in this
 manual and in the source code with respect to the terms of usage and distribution
 of SWAN. You are permitted to implement SWAN on your computer system.
+However, for any use of the SWAN source code in your environment, proper reference
+must be made to the origin of the software!
 \\[2ex]
 \noindent
 Implementation involves the following steps:
@@ -137,12 +139,12 @@
                          \> $\,\,$ swanout2.ftn \\
 service routines         \>:       swanser.ftn \\
                          \> $\,\,$ SwanIntgratSpc.ftn90 \\
-routines for support     \> \\
-parallel MPI runs        \>:       swanparll.ftn \\
 \end{tabbing}
 \newpage
 \begin{tabbing}
 xxxxxxxxxxxxxxxxxxxxxxxx\= \kill
+routines for support     \> \\
+parallel MPI runs        \>:       swanparll.ftn \\
 routines for unstructured \> $\,\,$ \\
 grids                     \>:       SwanReadGrid.ftn90 \\
                           \> $\,\,$ SwanReadADCGrid.ftn90 \\
@@ -175,6 +177,7 @@
                           \> $\,\,$ SwanInterpolateOutput.ftn90 \\
                           \> $\,\,$ SwanConvAccur.ftn90 \\
                           \> $\,\,$ SwanConvStopc.ftn90 \\
+                          \> $\,\,$ SwanThreadBounds.ftn90 \\
                           \> $\,\,$ SwanFindObstacles.ftn90 \\
                           \> $\,\,$ SwanCrossObstacle.ftn90 \\
                           \> $\,\,$ SwanComputeForce.ftn90 \\
@@ -199,14 +202,11 @@
 
 \noindent
 The source code is written in Fortran 90. Some routines are written in fixed form and depending on your system,
-the extension may be {\tt for}, {\tt f} or {\tt F}. Other routines are written in free form and are indicated by
-extension {\tt f90} or {\tt F90}. The conversion from {\tt ftn} or {\tt ftn90} to one of these extensions
+the extension may be {\tt for} or {\tt f}. Other routines are written in free form and are indicated by
+extension {\tt f90}. The conversion from {\tt ftn} or {\tt ftn90} to one of these extensions
 can be done automatically or manually; see Chapter \ref{ch:instal}.
 \\[2ex]
 \noindent
-For any use of the SWAN source code in your environment, proper reference must be made to the origin of the software!
-\\[2ex]
-\noindent
 You are allow to make changes in the source code of SWAN, but Delft University of Technology
 will not support modified versions of SWAN. If you ever want your modifications to be
 implemented in the authorized version of SWAN (the version on the SWAN web page), you need
@@ -302,18 +302,18 @@
 
 \chap{Installation of SWAN on your computer} \label{ch:instal}
 
-The portability of the SWAN code between single processor machines is guaranteed by the use of standard ANSI
+The portability of the SWAN code between different platforms is guaranteed by the use of standard ANSI
 FORTRAN 90. Hence, virtually all Fortran compilers can be used for installing SWAN. See also the manual Programming rules.
 \\[2ex]
 \noindent
 The SWAN code is parallelized, which enables a considerable reduction in the turn-around time for relatively
 large CPU-demanding calculations. Two parallelization strategies are available:
 \begin{itemize}
-  \item A message passing modelling is employed based on the Message Passing Interface (MPI) standard that
-        enables communication between independent processors. Only simple point-to-point and collective
-        communications have been employed. Hence, users can optionally run SWAN on a cluster of PC nodes.
   \item The computational kernel of SWAN contains a number of OpenMP compiler directives, so that users can
-        optionally run SWAN on shared-memory supercomputers or a dual-core PC.
+        optionally run SWAN on a dual core PC or multiprocessor systems.
+  \item A message passing modelling is employed based on the Message Passing Interface (MPI) standard that
+        enables communication between independent processors. Hence, users can optionally run SWAN on a
+        Linux cluster.
 \end{itemize}
 \noindent
 The material on the SWAN web site provides a {\tt Makefile} and two Perl scripts ({\tt platform.pl}
@@ -370,19 +370,20 @@
 (e.g., instationary or nesting ones), two ways of parallelism for reducing the turn-around time are
 available:
 \begin{itemize}
-   \item The SWAN code contains a number of so-called OpenMP directives that tells the
-         compiler how to generate multi-threaded code on a shared memory computer. For this,
-         you need a Fortran 90 compiler having the OpenMP option. The performance is good
-         for a restricted number of threads ($< 8$).This type of parallelism can be used e.g.,
-         on symmetric multiprocessors and PC's with dual-core or quad-core processors.
+   \item The SWAN code contains a number of so-called OpenMP directives that enables the
+         compiler to generate multithreaded code on a shared memory computer. For this,
+         you need a Fortran 90 compiler supporting OpenMP 2.0. The performance is good
+         for a restricted number of threads ($< 8$). This type of parallelism can be used e.g.,
+         on (symmetric) multiprocessors and PC's with dual core or quad core processors.
    \item If the user want to run SWAN on a relative large number of processors, a message passing
          model is a good alternative. It is based on independent processors which do not share
          any memory but are connected via an interconnection network, such as Beowulf systems
-         (cluster of Linux PC's connected via fast Ethernet switches). Since, the Message Passing
-         Interface (MPI) standard (e.g., MPICH distribution, freely available for several platforms,
-         such as Linux and Windows NT/2000/XP, at \hl{http://www-unix.mcs.anl.gov/mpi/mpich}) is
-         very popular nowadays, the SWAN code contains a set of generic subroutines that call a
-         number of MPI-routines, meant for local data exchange, gathering data, global reductions,
+         (cluster of Linux PC's connected via fast Ethernet switches). The implementation is based on
+         the Message Passing Interface (MPI) standard (e.g., MPICH2 distribution, freely available
+         for several platforms, such as Linux and Windows NT/2000/XP, at
+         \hl{http://www-unix.mcs.anl.gov/mpi/mpich}).
+         The SWAN code contains a set of generic subroutines that call a
+         number of MPI routines, meant for local data exchange, gathering data, global reductions,
          etc. This technique is beneficial for larger simulations only, such that the communication
          times are relatively small compared to the computing times.
 \end{itemize}
@@ -426,7 +427,7 @@
                  (see http://support.microsoft.com/default.aspx?scid=kb;en-us;Q132084). Run the
                  downloaded file ({\tt Nmake15.exe}) to extract it. Copy both the {\tt NMAKE.EXE} and
                  the {\tt NMAKE.ERR} file to SWAN directory.
-           \item This setup does support OpenMP for Windows dual-core systems, if Intel Fortran
+           \item This setup does support OpenMP for Windows dual core systems, if Intel Fortran
                  compiler is provided.
            \item This installation currently supports MPICH for
                  Windows NT/2000/XP (Professional); Win9x/ME are not supported.
@@ -438,7 +439,7 @@
         \end{itemize}
   \item One of the commands {\tt make ser}, {\tt make omp} and {\tt make mpi}
         must be preceded by {\tt make config}.
-  \item If desirable, you may clean-up the generated object files and modules by
+  \item If desirable, you may clean up the generated object files and modules by
         typing {\tt make~clean}. If you want to go back to the original state
         with respect to the source code, i.e. removing everything that has been
         generated by the Makefile, just type {\tt make~clobber}.
@@ -452,14 +453,14 @@
 \subsection{Modifications in the source code}
 
 To compile SWAN on your computer system properly, some subroutines should be adapted first
-depending on the operating system, use of compilers and the wish to use MPI or OpenMP for
+depending on the operating system, use of compilers and the wish to use MPI for
 parallel runs. This can be done by removing the switches started with '!' followed by an
 indentifiable prefix in the first 3 or 4 columns of the subroutine. A Perl script called
 {\tt switch.pl} is provided in the material that enables the user to quickly select the
 switches to be removed. This script can be used as follows:
 \begin{verbatim}
-perl switch.pl [-dos] [-unix] [-f95] [-mpi] [-omp] [-cray]
-               [-sgi] [-cvis] [-timg] [-impi] *.ftn
+perl switch.pl [-dos] [-unix] [-f95] [-mpi] [-cray] [-sgi]
+               [-cvis] [-timg] [-impi] *.ftn
 \end{verbatim}
 where the options are all optionally. The meaning of these options are as follows.
 \begin{itemize}
@@ -480,10 +481,6 @@
   \item[{\tt -mpi}]
   For the proper use of MPI, you must remove the switch !MPI at several places in the file
   {\tt swanparll.ftn}, {\tt swancom1.ftn} and {\tt swmod1.ftn}.
-  \item[{\tt -omp}]
-  The subroutine {\tt SWCOMP} (in {\tt swancom1.ftn}) contains a number of {\tt \_OPENMP}
-  macro that needs to be set by first preprocessing it by the C (or Fortran) preprocessor.
-  In order to use this macro just remove the !OMP switch.
   \item[{\tt -cray}, {\tt -sgi}]
   If you use a Cray or SGI Fortran 90 compiler, the subroutines {\tt OCPINI} (in {\tt ocpids.ftn})
   and {\tt FOR} (in {\tt ocpmix.ftn}) should be adapted by removing the switch !/Cray or !/SGI since,
@@ -506,7 +503,7 @@
 \end{itemize}
 
 \noindent
-For example, you work on a Beowulf cluster where MPI has been installed and use the Intel Fortran
+For example, you work on a Linux cluster where MPI has been installed and use the Intel Fortran
 compiler (that can handle Fortran 95 statements), then type the following:
 \begin{verbatim}
 perl switch.pl -unix -f95 -mpi *.ftn *.ftn90
@@ -517,13 +514,11 @@
 \subsection{Compiling and linking SWAN source code}
 
 After the necessary modifications are made as described in the previous section, the source code is
-ready for compilation. All source code is written in fixed form Fortran 90 so you must have a
+ready for compilation. All source code is written in Fortran 90 so you must have a
 Fortran 90 compiler in order to compile SWAN. The source code cannot be compiled with a Fortran 77
 compiler. If you intended to use MPI for parallel runs, you must use the command {\tt mpif90} instead
 of the original compiler command or using the Integrated Development Environment e.g., for Visual
-Fortran (see Installation and User's Guide for MPICH). For parallel runs using OpenMP, the compiler
-must have an option to interpret OpenMP directives and the extension {\tt ftn} must be changed into
-{\tt F} (will be done automatically when using the script {\tt switch.pl}).
+Fortran (see Installation and User's Guide for MPICH).
 \\[2ex]
 \noindent
 The SWAN source code complies with the ANSI Fortran 90 standard, except for a few cases, where
@@ -539,7 +534,7 @@
 \\[2ex]
 \noindent
 Furthermore, SWAN can generate binary MATLAB files on request, which are unformatted. Some compilers,
-e.g. Compaq Visual Fortran and Intel Fortran version 9.x, measured record length in 4-byte units and
+e.g. Compaq Visual Fortran and Intel Fortran version 9.x, measured record length in 4-byte or longword units and
 as a consequence, these unformatted files cannot be loaded in MATLAB. Hence, in such as case a
 compiler option is needed to request 1-byte units, e.g. for Compaq Visual Fortran this is
 {\tt /assume:byterecl} and for Intel Fortran version 9.x this is {\tt -assume byterecl}.
@@ -563,7 +558,7 @@
 \newpage
 \section{Make SWAN documentation} \label{sec:makedoc}
 
-SWAN comes with 4 detailed documents which are provided as downloadable PDF files as well as browsable web-pages:
+SWAN comes with 4 detailed documents which are provided as downloadable PDF files as well as browsable web pages:
 \begin{itemize}
   \item The \underline{User Manual} describes the complete input and usage of the SWAN package.
   \item The \underline{Implementation Manual} explains the installation procedure of SWAN on a single- or multi-processor machine
@@ -675,7 +670,7 @@
 \noindent
 Two run procedures are provided among the source code, one for the Windows platform, called
 {\tt swanrun.bat}, and one for the UNIX/Linux platform, called {\tt swanrun}. Basically, the
-following actions need to be done by the run procedure:
+run procedure carries out the following actions:
 \begin{itemize}
   \item Copy the command file with extension {\tt swn} to {\tt INPUT} (assuming {\tt INPUT} is
         the standard file name for command input, see Chapter \ref{ch:usechan}).
@@ -688,22 +683,42 @@
 program {\tt mpirun} is needed and is provided in the MPICH distribution.
 \\[2ex]
 \noindent
-Before calling the run procedure, the environment variable PATH need to be adapted by
+Before calling the run procedure, the environment variable {\small \tt PATH} need to be adapted by
 including the pathname of the directory where {\tt swan.exe} can be found. In case of Windows,
 this pathname can be specified through the category {\it System} of {\it Control Panel} (on the
 {\it Advanced} tab, click {\it Environment Variables}) or by adding it in the {\tt AUTOEXEC.BAT}
-file. In case of UNIX or Linux running the Bourne shell, the environment variable PATH may be
-changed as follows:
+file. In case of UNIX or Linux running the bash shell (sh or ksh), the environment variable
+{\small \tt PATH} may be changed as follows:
 \begin{verbatim}
 export PATH=${PATH}:/usr/local/swan
 \end{verbatim}
-if {\tt /usr/local/swan} is the directory where the executable {\tt swan.exe} can be found.
-In case of the C shell, use the following command:
+if {\tt /usr/local/swan} is the directory where the executable {\tt swan.exe} is resided.
+In case of the C shell (csh), use the following command:
 \begin{verbatim}
 setenv PATH ${PATH}:/usr/local/swan
 \end{verbatim}
+
+\noindent
 If appropriate, you also need to add the directory path where the {\tt bin} directory of MPICH
-is resided to PATH to have access to the command {\tt mpirun}.
+is resided to {\small \tt PATH} to have access to the command {\tt mpirun}.
+\\[2ex]
+\noindent
+You may also specify the number of threads to be used during execution of the multithreaded
+implementation of SWAN on multiprocessor systems. The environment variable for this is
+{\small \tt OMP\_NUM\_THREADS} and can be set like
+\begin{verbatim}
+export OMP_NUM_THREADS=4
+\end{verbatim}
+or
+\begin{verbatim}
+setenv OMP_NUM_THREADS 4
+\end{verbatim}
+or, in case of Windows,
+\begin{verbatim}
+OMP_NUM_THREADS = 4
+\end{verbatim}
+When dynamic adjustment of the number of threads is enabled, the value given in\\
+{\small \tt OMP\_NUM\_THREADS} represents the maximum number of threads allowed.
 \\[2ex]
 \noindent
 The provided run procedures enable the user to properly and easily run SWAN both serial as
@@ -720,27 +735,34 @@
 where {\tt filename} is the name of your command file without extension (assuming it is {\tt swn})
 and {\tt nprocs} indicates how many processes need to be launched for a parallel MPI run (do not
 type the brackets; they just indicate that {\tt nprocs} is optional). By default, $nprocs = 1$.
-You may also run on a dual-core PC with the appropriate SWAN executable; do not set {\tt nprocs}.
+You may also run on a dual/quad core computer; do not set {\tt nprocs}.
 \\[2ex]
 \noindent
 The command line for the UNIX script {\tt swanrun} is as follows:
 \begin{verbatim}
 ./swanrun -input filename [-omp n | -mpi n]
 \end{verbatim}
-where {\tt filename} is the name of your command file without extension. The parameter {\tt -omp~n}
-specifies a parallel run on $n$ processors using OpenMP. The parameter {\tt -mpi~n} specifies a
-parallel run on $n$ processors using MPI. The parameter {\tt -input} is obliged, whereas the
-parameters {\tt -omp~n} and {\tt -mpi~n} can be omitted (default: $n=1$). Note that the script
+where {\tt filename} is the name of your command file without extension. Note that the script
 {\tt swanrun} need to be made executable first, as follows:
 \begin{verbatim}
 chmod +rx ./swanrun
 \end{verbatim}
+The parameter {\tt -omp~n}
+specifies a parallel run on $n$ cores using OpenMP. Note that the UNIX script will set
+{\small \tt OMP\_NUM\_THREADS} to $n$. The parameter {\tt -mpi~n} specifies a
+parallel run on $n$ processors using MPI. The parameter {\tt -input} is obliged, whereas the
+parameters {\tt -omp~n} and {\tt -mpi~n} can be omitted (default: $n=1$). To redirect screen
+output to a file, use the sign $>$. Use an ampersand to run SWAN in the background. An
+example:
+\begin{verbatim}
+./swanrun -input f31har01 -omp 4 > swanout &
+\end{verbatim}
 
 \noindent
 For a parallel MPI run, you may also need a {\tt machinefile} that contains the names of the
 nodes in your parallel environment. Put one node per line in the file. Lines starting with
 the \# character are comment lines. You can specify a number after the node name to indicate
-how many processes to launch on the node. This is useful e.g., for dual-processors. The run
+how many cores to launch on the node. This is useful e.g., for multi-core processors. The run
 procedure will cycle through this list until all the requested processes are launched. Example
 of such a file may look like:
 \begin{verbatim}
@@ -785,8 +807,8 @@
 
 The SWAN system consists of one executable file ({\tt swan.exe}), a command file ({\tt swan.edt}) and
 a run procedure ({\tt swanrun.bat} or {\tt swanrun}).
-Depending on your system, you may use 32-bit or 64-bit single-core or 32-/64-bit dual-core executable.
-These executables can be obtained from the SWAN web site.
+Depending on your system, you may use 32-bit or 64-bit executable.
+These executables for Windows NT/2000/XP can be obtained from the SWAN web site.
 The input and output to a number of test problems
 is provided on the SWAN web page. The files with extension {\tt swn} are the command files for these
 tests; the files with extension {\tt bot} are the bottom files for these tests, etc. This input can be
--- SwanInterpolateOutput.ftn90	2009-08-24 13:08:40.000000000 +0200
+++ SwanInterpolateOutput.ftn90	2009-08-24 13:06:43.000000000 +0200
@@ -65,8 +65,8 @@
 !
 !   Argument variables
 !
-    integer, dimension(mip), intent(in) :: kvert  ! vertex indices of output points
     integer, intent(in)                 :: mip    ! number of given points
+    integer, dimension(mip), intent(in) :: kvert  ! vertex indices of output points
     !
     real, intent(in)                    :: excval ! exception value for output quantity
     real, dimension(nverts), intent(in) :: finp   ! output quantity defined on the computational grid
--- swanmain.ftn	2009-08-24 13:08:40.000000000 +0200
+++ swanmain.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -858,6 +858,7 @@
       CALL BUGFIX ('A')
       CALL BUGFIX ('B')
       CALL BUGFIX ('C')
+      CALL BUGFIX ('D')
 !
       CALL OCPINI ('swaninit', .TRUE.,INERR)                              34.01
       IF (INERR.GT.0) RETURN                                              34.01
@@ -1031,6 +1032,7 @@
       ONED   = .FALSE.                                                    32.02
       ACUPDA = .TRUE.                                                     40.07
       OFFSRC = .FALSE.                                                    40.80
+      LADDS  = .FALSE.                                                    40.85
       HSRERR = 0.1                                                        32.01
 !
 !     higher order propagation and spherical coordinates                  33.08
@@ -4577,10 +4579,6 @@
          IF (EQREAL(PWCAP( 9),    0.0)) PWCAP( 9) = 0.                    40.53
          IF (EQREAL(PWCAP(11),    0.0)) PWCAP(11) = 0.                    40.53
       END IF                                                              40.53
-!OMP!
-!OMP      IF ( IQUAD.EQ.51 .OR. IQUAD.EQ.52 .OR. IQUAD.EQ.53 ) THEN           40.41
-!OMP         CALL MSGERR(4,'XNL not supported within OpenMP environment')     40.41
-!OMP      END IF                                                              40.41
 !
       IF ( ITEST .GE. 120 ) THEN
         WRITE(PRINTF,3000) IWIND ,IQUAD, ICUR, IWCAP, MSC
--- swanout2.ftn	2009-08-24 13:08:40.000000000 +0200
+++ swanout2.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -174,7 +174,7 @@
         IF (MATLAB) THEN                                                  40.30
            CLOSE(NREF)                                                    40.30
            OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',               40.30
-     &          ACCESS='DIRECT', RECL=1)                                  40.30
+     &          ACCESS='DIRECT', RECL=4)                                  41.08 40.30
            IREC(IRQ) = 1                                                  40.51
         END IF
       ENDIF
@@ -1849,6 +1849,7 @@
 !
 !****************************************************************
 !
+      USE OCPCOMM2                                                        41.08
       USE OCPCOMM4                                                        40.41
 !
       IMPLICIT NONE
@@ -1887,11 +1888,13 @@
 !
 !     40.30: Marcel Zijlema
 !     40.41: Marcel Zijlema
+!     41.08: Pieter Smit
 !
 !  1. Updates
 !
-!     40.30, May 03: New subroutine
+!     40.30, May  03: New subroutine
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
+!     41.08, Aug. 09: adapted to write binary file in the Level 5 MAT-file format
 !
 !  2. Purpose
 !
@@ -1904,7 +1907,7 @@
 !        statement:
 !
 !        OPEN(UNIT=IOUTMA, FILE=BINFIL, FORM='UNFORMATTED',
-!             ACCESS='DIRECT', RECL=1)
+!             ACCESS='DIRECT', RECL=4)
 !
 !        Furthermore, initialize record counter to IREC = 1
 !
@@ -1912,7 +1915,7 @@
 !        matrices to be saved
 !
 !     3) The matrix may contain signed infinity and/or Not a Numbers.
-!        According to the IEEE standard, on a 32-bit machine, the real
+!        According to the IEEE 754 standard, on a 32-bit machine, the real
 !        format has an 8-bit biased exponent (=actual exponent increased
 !        by bias=127) and a 23-bit fraction or mantissa. The leftmost
 !        bit is the sign bit. Let a fraction, biased exponent and sign
@@ -1927,14 +1930,11 @@
 !        representation of a NaN equals the representation of +Inf
 !        plus 1.
 !
-!        The Cray machine C916 at SARA Information Centre does not
-!        support the IEEE standard.
-!
 !     4) The NaN's or Inf's are indicated by a dummy value as given
 !        by dumval
 !
-!     For more information consult "Appendix - MAT-File Structure"
-!     of the MATLAB External Data Reference guide (Version 4.2)
+!     For more information on the Level 5 MAT-file format consult
+!     document "MAT-File Format" of MathWorks
 !
 !  4. Argument variables
 !
@@ -1956,67 +1956,55 @@
 !
 !  5. Parameter variables
 !
-!     ---
+!     BlockSize   size of matlab data segment
+!     DataSize    number of bytes written per write statement
+!     HeaderSize  size of the header in bytes
+!     mChar       character data
+!     mInt32      signed   integer*4
+!     mUInt32     unsigned integer*8
+!     mSingle     real
+!
+!     --- standard sizes
+!
+      INTEGER, PARAMETER :: DataSize   = 4
+      INTEGER, PARAMETER :: HeaderSize = 128
+      INTEGER, PARAMETER :: BlockSize  = 8
+!
+!     --- Matlab data types
+!
+      INTEGER, PARAMETER :: mChar      = 4
+      INTEGER, PARAMETER :: mInt32     = 5
+      INTEGER, PARAMETER :: mUInt32    = 6
+      INTEGER, PARAMETER :: mSingle    = 7
 !
 !  6. Local variables
 !
-!     BVAL  :     a byte value
-!     CHARS :     array to pass character info to MSGERR
+!     CTMP  :     a temporary character array
+!     HEADER:     header of binary MAT-file
 !     I     :     loop variable
 !     IENT  :     number of entries
-!     IF    :     first non-character in string
-!     IL    :     last non-character in string
-!     IMAGF :     a 4-byte imaginary flag. Possible values are:
-!                 0: there is only real data
-!                 1: the data has also an imaginary part
 !     IOS   :     auxiliary integer with iostat-number
-!     ITYPE :     the type flag containing a 4-byte integer whose
-!                 decimal digits encode storage information.
-!                 If the integer is represented as ABCD then:
-!                 "A" indicates the format to write the binary
-!                 data to a file on the machine. Possible values are:
-!                   0: Intel based machines (PC 386/486, Pentium)
-!                   1: Motorola 68000 based machines (Macintosh,
-!                      HP 9000, SPARC, Apollo, SGI)
-!                   2: VAX-D format
-!                   3: VAX-G format
-!                   4: Cray
-!                 "B" is always zero
-!                 "C" indicates which format the data is stored.
-!                  Possible values are:
-!                   0: double precision (64 bit) floating point numbers
-!                   1: single precision (32 bit) floating point numbers
-!                   2: 32-bit signed integers
-!                   3: 16-bit signed integers
-!                   4: 16-bit unsigned integers
-!                   5: 8-bit unsigned integers
-!                 "D" indicates the type of data (matrix).
-!                  Possible values:
-!                   0: numeric matrix
-!                   1: textual matrix
-!                   2: sparse  matrix
+!     IRECS :     size of array including tags and flags
 !     J     :     index
 !     M     :     loop variable
 !     MSGSTR:     string to pass message to call MSGERR
 !     N     :     loop variable
 !     NAMLEN:     a 4-byte integer representing the number of
-!                 characters in matrix name plus 1
+!                 characters in matrix name
 !     NANVAL:     an integer representing Not a Number
+!     NTOT  :     size of data array
 !
-      INTEGER I, J, IENT, IF, IL, IOS, M, N
-      INTEGER BVAL(4), IMAGF, ITYPE, NAMLEN, NANVAL
-      CHARACTER*20 INTSTR, CHARS
+      INTEGER I, J, IENT, IOS, M, N, NTOT
+      INTEGER NAMLEN, NANVAL
+      INTEGER, SAVE :: IRECS
       CHARACTER*80 MSGSTR
+      CHARACTER(LEN=HeaderSize) HEADER
+      CHARACTER(LEN=BlockSize) CTMP
+
 !
 !  8. Subroutines used
 !
-!     INTSTR           Converts integer to string
 !     MSGERR           Writes error message
-!     TXPBLA           Removes leading and trailing blanks in string
-!     SWI2B            Calculates 32-bit representation of an
-!                      integer number
-!     SWR2B            Calculates 32-bit representation of a
-!                      floating-point number
 !
 !  9. Subroutines calling
 !
@@ -2033,14 +2021,12 @@
 ! 12. Structure
 !
 !     set Not a Number
-!
-!     set some flags
-!
-!     write header consisting of ITYPE, MROWS, NCOLS, IMAGF, NAMLEN and
-!     name of matrix MATNAM
-!
+!     length of name matrix
+!     size of data array
+!     write header once
+!     array name
 !     write matrix
-!
+!     write the size of the array
 !     if necessary, give message that error occurred while writing file
 !
 ! 13. Source text
@@ -2049,61 +2035,121 @@
       DATA IENT/0/
       IF (LTRACE) CALL STRACE (IENT,'SWRMAT')
 
+      IOS = 0
+
 !     --- set Not a Number
 
       NANVAL = 255 * 2**23 + 1
 
-!     --- set some flags
+!     --- length of name matrix
 
-      ITYPE = 1010
-      IMAGF = 0
-      IOS   = 0
-
-!     --- write header consisting of ITYPE, MROWS, NCOLS, IMAGF,
-!         NAMLEN and name of matrix MATNAM
-!         the name should be ended by zero-byte terminator
-
-      CALL SWI2B ( ITYPE, BVAL )
-      DO I = 1, 4
-         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
-         IREC = IREC + 1
-      END DO
+      NAMLEN = LEN_TRIM(MATNAM)
 
-      CALL SWI2B ( MROWS, BVAL )
-      DO I = 1, 4
-         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
-         IREC = IREC + 1
-      END DO
+!     --- size of data array
 
-      CALL SWI2B ( NCOLS, BVAL )
-      DO I = 1, 4
-         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
-         IREC = IREC + 1
-      END DO
+      NTOT = MROWS * NCOLS
 
-      CALL SWI2B ( IMAGF, BVAL )
-      DO I = 1, 4
-         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
-         IREC = IREC + 1
-      END DO
+!     --- descriptive header
+
+      WRITE (HEADER, '(6A)') 'Data produced by SWAN version ',
+     &                       TRIM(VERTXT),'; project: ',TRIM(PROJID),
+     &                       '; run number: ',PROJNR
+
+!     --- data offset
+
+      HEADER(117:124) = CHAR(ICHAR(' '))
+
+!     --- version
+
+      HEADER(125:126) = CHAR(0) // CHAR(1)
+
+!     --- endian indicator
+
+      WRITE(HEADER(127:128),'(A)') INT(19785,KIND=2)
+
+!     --- write header once
+
+      IF ( IREC.EQ.1 ) THEN
+         DO I = 1, HeaderSize/DataSize
+            J = DataSize*(I-1) + 1
+            IF ( IOS.EQ.0 )
+     &         WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) HEADER(J:J+DataSize-1)
+            IREC = IREC + 1
+         END DO
+      END IF
+
+!     --- array tag
+
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 14                ! Matlab array
+      IREC = IREC + 1
+      IRECS = IREC
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
+      IREC = IREC + 1
+
+!     --- array flags
 
-      CALL TXPBLA(MATNAM,IF,IL)
-      NAMLEN = IL - IF + 2
-      CALL SWI2B ( NAMLEN, BVAL )
-      DO I = 1, 4
-         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mInt32
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 2*DataSize
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 7                 ! single precision array
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
+      IREC = IREC + 1
+
+      IF ( MOD(2,BlockSize/DataSize).NE.0 ) THEN
+         IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
          IREC = IREC + 1
-      END DO
+      END IF
+
+!     --- dimensions array
 
-      DO I = IF, IL
-         IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) MATNAM(I:I)
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mInt32
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 2*DataSize
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) MROWS
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) NCOLS
+      IREC = IREC + 1
+
+      IF ( MOD(2,BlockSize/DataSize).NE.0 ) THEN
+         IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) 0
          IREC = IREC + 1
-      END DO
-      IF (IOS.EQ.0) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(0)
+      END IF
+
+!     --- array name
+
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) mChar
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) NAMLEN
       IREC = IREC + 1
 
+      I = 1
+      DO
+        CTMP(1:8) = CHAR(ICHAR(' '))
+        IF ( I.GT.NAMLEN ) THEN
+           EXIT
+        ELSE IF ( I+BlockSize.LE.NAMLEN ) THEN
+           CTMP(1:8) = MATNAM(I:I+BlockSize-1)
+        ELSE
+           CTMP(1:NAMLEN-I+1) = MATNAM(I:NAMLEN)
+        END IF
+
+        IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CTMP(1:4)
+        IREC = IREC + 1
+        IF ( IOS.EQ.0 ) WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CTMP(5:8)
+        I = I + BlockSize
+        IREC = IREC + 1
+      END DO
+
 !     --- write matrix
 
+      IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) mSingle
+      IREC = IREC + 1
+      IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) NTOT*DataSize
+      IREC = IREC + 1
+
       DO M = 1, NCOLS
          DO N = 1, MROWS
             IF ( IDLA.EQ.1 ) THEN
@@ -2111,27 +2157,32 @@
             ELSE
                J = (N-1)*NCOLS + M
             END IF
-            IF (RDATA(J).EQ.DUMVAL ) THEN
-               CALL SWI2B ( NANVAL, BVAL )
+            IF (RDATA(J).NE.DUMVAL ) THEN
+               WRITE (IOUTMA,REC=IREC) RDATA(J)
             ELSE
-               CALL SWR2B ( RDATA(J), BVAL )
+               WRITE (IOUTMA,REC=IREC) NANVAL
             END IF
-            DO I = 1, 4
-              IF (IOS.EQ.0)
-     &                   WRITE(IOUTMA,REC=IREC,IOSTAT=IOS) CHAR(BVAL(I))
-              IREC = IREC + 1
-            END DO
+            IREC = IREC + 1
          END DO
       END DO
 
+      IF ( MOD(NTOT,BlockSize/DataSize).NE.0 ) THEN
+         IF ( IOS.EQ.0 ) WRITE (IOUTMA,REC=IREC,IOSTAT=IOS) 0.
+         IREC = IREC + 1
+      END IF
+
+!     --- write the size of the array
+
+      IF ( IOS.EQ.0 )
+     &   WRITE (IOUTMA,REC=IRECS,IOSTAT=IOS) (IREC-IRECS-1)*DataSize
+
 !     --- if necessary, give message that error occurred while writing file
 
       IF ( IOS.NE.0 ) THEN
-         CHARS = INTSTR(IOS)
-         CALL TXPBLA(CHARS,IF,IL)
-         MSGSTR = 'Error while writing binary MAT-file - '//
-     &            'IOSTAT number is '//CHARS(IF:IL)
-         CALL MSGERR ( 4, MSGSTR )
+         WRITE (MSGSTR, '(A,I5)')
+     &            'Error while writing binary MAT-file - '//
+     &            'IOSTAT number is ', IOS
+         CALL MSGERR( 4, TRIM(MSGSTR) )
          RETURN
       END IF
 
--- swanparll.ftn	2009-08-24 13:08:40.000000000 +0200
+++ swanparll.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -4347,7 +4347,7 @@
             IF (MATLAB) THEN
                CLOSE(NREF)
                OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',
-     &              ACCESS='DIRECT', RECL=1)
+     &              ACCESS='DIRECT', RECL=4)                              41.08
                IREC(IRQ) = 1                                              40.51
             END IF
 
--- swanpre2.ftn	2009-08-24 13:08:38.000000000 +0200
+++ swanpre2.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -1328,6 +1328,9 @@
                    JRADS  = MCMVAR                                        40.85
                    ALOCMP = .TRUE.                                        40.97
                 ENDIF                                                     40.85
+                IF (IVTYPE.EQ.7  .OR. IVTYPE.EQ.9  .OR.                   40.85
+     &              IVTYPE.EQ.54 .OR. IVTYPE.EQ.55 .OR.                   40.85
+     &              IVTYPE.EQ.56 .OR. IVTYPE.GE.60 ) LADDS = .TRUE.       40.85
              ENDIF
              GOTO 70
           ENDIF
@@ -1517,6 +1520,9 @@
                JRADS  = MCMVAR                                            40.85
                ALOCMP = .TRUE.                                            40.97
             ENDIF                                                         40.85
+            IF (IVTYPE.EQ.7  .OR. IVTYPE.EQ.9  .OR.                       40.85
+     &          IVTYPE.EQ.54 .OR. IVTYPE.EQ.55 .OR.                       40.85
+     &          IVTYPE.EQ.56 .OR. IVTYPE.GE.60 ) LADDS = .TRUE.           40.85
             CALL INKEYW ('STA', ' ')                                      40.00
             IF (KEYWIS('UNIT')) THEN
               CALL MSGERR (1, 'UNIT is ignored in this version')          40.00
--- SwanReadADCGrid.ftn90	2009-08-24 13:08:40.000000000 +0200
+++ SwanReadADCGrid.ftn90	2009-08-24 13:06:46.000000000 +0200
@@ -60,6 +60,7 @@
 !
 !   Local variables
 !
+    character(80)           :: filenm   ! file name
     integer, save           :: ient = 0 ! number of entries in this subroutine
     integer                 :: idum     ! dummy integer
     integer                 :: ii       ! auxiliary integer
@@ -88,9 +89,10 @@
     !
     ! open file fort.14
     !
+    filenm = 'fort.14'
     ndsd   = 0
     iostat = 0
-    call for (ndsd, 'fort.14', 'OF', iostat)
+    call for (ndsd, trim(filenm), 'OF', iostat)
     if (stpnow()) goto 900
     !
     ! skip first line
--- swanser.ftn	2009-08-24 13:08:38.000000000 +0200
+++ swanser.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -31,8 +31,6 @@
 !     NUMSTR                                                              40.23
 !     SWCOPI                                                              40.23
 !     SWCOPR                                                              40.23
-!     SWI2B                                                               40.30
-!     SWR2B                                                               40.30
 !
 !  functions:
 !  ----------
@@ -6887,400 +6885,3 @@
 
       RETURN
       END
-!****************************************************************
-!
-      SUBROUTINE SWI2B ( IVAL, BVAL )
-!
-!****************************************************************
-!
-      USE OCPCOMM4                                                        40.41
-!
-      IMPLICIT NONE
-!
-!
-!   --|-----------------------------------------------------------|--
-!     | Delft University of Technology                            |
-!     | Faculty of Civil Engineering and Geosciences              |
-!     | Environmental Fluid Mechanics Section                     |
-!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
-!     |                                                           |
-!     | Programmer: Marcel Zijlema                                |
-!   --|-----------------------------------------------------------|--
-!
-!
-!     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2009  Delft University of Technology
-!
-!     This program is free software; you can redistribute it and/or
-!     modify it under the terms of the GNU General Public License as
-!     published by the Free Software Foundation; either version 2 of
-!     the License, or (at your option) any later version.
-!
-!     This program is distributed in the hope that it will be useful,
-!     but WITHOUT ANY WARRANTY; without even the implied warranty of
-!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-!     GNU General Public License for more details.
-!
-!     A copy of the GNU General Public License is available at
-!     http://www.gnu.org/copyleft/gpl.html#SEC3
-!     or by writing to the Free Software Foundation, Inc.,
-!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
-!
-!
-!  0. Authors
-!
-!     40.30: Marcel Zijlema
-!     40.41: Marcel Zijlema
-!
-!  1. Updates
-!
-!     40.30, May 03: New subroutine
-!     40.41, Oct. 04: common blocks replaced by modules, include files removed
-!
-!  2. Purpose
-!
-!     Calculates 32-bit representation of an integer number
-!
-!  3. Method
-!
-!     The representation of an integer number is divided into 4 parts
-!     of 8 bits each, resulting in 32-bit word in memory. Generally,
-!     storage words are represented with bits counted from the right,
-!     making bit 0 the lower-order bit and bit 31 the high-order bit,
-!     which is also the sign bit.
-!
-!     The integer number is always an exact representation of an
-!     integer of value positive, negative, or zero. Each bit, except
-!     the leftmost bit, corresponds to the actual exponent as power
-!     of two.
-!
-!     For representing negative numbers, the method called
-!     "excess 2**(m - 1)" is used, which represents an m-bit number by
-!     storing it as the sum of itself and 2**(m - 1). For a 32-bit
-!     machine, m = 32. This results in a positive number, so the
-!     leftmost bit need to be reversed. This method is identical to the
-!     two's complement method.
-!
-!     An example:
-!
-!        the 32-bit representation of 5693 is
-!
-!        decimal    :     0        0       22       61
-!        hexidecimal:     0        0       16       3D
-!        binary     : 00000000 00000000 00010110 00111101
-!
-!        since,
-!
-!        5693 = 2^12 + 2^10 + 2^9 + 2^5 + 2^4 + 2^3 + 2^2 + 2^0
-!
-!  4. Argument variables
-!
-!     BVAL        a byte value as a part of the representation of
-!                 integer number
-!     IVAL        integer number
-!
-      INTEGER BVAL(4), IVAL
-!
-!  6. Local variables
-!
-!     I     :     loop counter
-!     IENT  :     number of entries
-!     IQUOT :     auxiliary integer with quotient
-!     M     :     maximal exponent number possible (for 32-bit machine, m=32)
-!
-      INTEGER I, IENT, IQUOT, M
-!
-! 12. Structure
-!
-!     initialise 4 parts of the representation
-!     if integer < 0, increased it by 2**(m-1)
-!     compute the actual part of the representation
-!     determine the sign bit
-!
-! 13. Source text
-!
-      SAVE IENT
-      DATA IENT/0/, M/32/
-      IF (LTRACE) CALL STRACE (IENT,'SWI2B')
-
-!     --- initialise 4 parts of the representation
-
-      DO 10 I = 1, 4
-         BVAL(I) = 0
- 10   CONTINUE
-
-      IQUOT = IVAL
-
-!     --- if integer < 0, increased it by 2*(m-1)
-
-      IF ( IVAL.LT.0 ) IQUOT = IQUOT + 2**(M-1)
-
-!     --- compute the actual part of the representation
-
-      DO 20 I = 4, 1, -1
-         BVAL(I) = MOD(IQUOT,256)
-         IQUOT = INT(IQUOT/256)
- 20   CONTINUE
-
-!     --- determine the sign bit
-
-      IF ( IVAL.LT.0 ) BVAL(1) = BVAL(1) + 128
-
-      RETURN
-      END
-!****************************************************************
-!
-      SUBROUTINE SWR2B ( RVAL, BVAL )
-!
-!****************************************************************
-!
-      USE OCPCOMM4                                                        40.41
-!
-      IMPLICIT NONE
-!
-!
-!   --|-----------------------------------------------------------|--
-!     | Delft University of Technology                            |
-!     | Faculty of Civil Engineering and Geosciences              |
-!     | Environmental Fluid Mechanics Section                     |
-!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
-!     |                                                           |
-!     | Programmer: Marcel Zijlema                                |
-!   --|-----------------------------------------------------------|--
-!
-!
-!     SWAN (Simulating WAves Nearshore); a third generation wave model
-!     Copyright (C) 2009  Delft University of Technology
-!
-!     This program is free software; you can redistribute it and/or
-!     modify it under the terms of the GNU General Public License as
-!     published by the Free Software Foundation; either version 2 of
-!     the License, or (at your option) any later version.
-!
-!     This program is distributed in the hope that it will be useful,
-!     but WITHOUT ANY WARRANTY; without even the implied warranty of
-!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-!     GNU General Public License for more details.
-!
-!     A copy of the GNU General Public License is available at
-!     http://www.gnu.org/copyleft/gpl.html#SEC3
-!     or by writing to the Free Software Foundation, Inc.,
-!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
-!
-!
-!  0. Authors
-!
-!     40.30: Marcel Zijlema
-!     40.41: Marcel Zijlema
-!
-!  1. Updates
-!
-!     40.30, May 03: New subroutine
-!     40.41, Oct. 04: common blocks replaced by modules, include files removed
-!
-!  2. Purpose
-!
-!     Calculates 32-bit representation of a floating-point number
-!
-!  3. Method
-!
-!     The representation of a floating-point number is divided into 4
-!     parts of 8 bits each, resulting in 32-bit word in memory.
-!     Generally, storage words are represented with bits counted from
-!     the right, making bit 0 the lower-order bit and bit 31 the
-!     high-order bit, which is also the sign bit.
-!
-!     The floating-point number is a processor approximation. Its format
-!     has an 8-bit biased exponent and a 23-bit fraction or mantissa. The
-!     leftmost bit is the sign bit which is zero for plus and 1 for
-!     minus. The biased exponent equals the bias and the actual exponent
-!     (power of two) of the number. For a 32-bit machine, bias=127.
-!
-!     Furthermore, the floating-point number is usually stored in the
-!     normalized form, i.e. it has a binary point to the left of the
-!     mantissa and an implied leading 1 to the left of the binary point.
-!     Thus, if X is a floating-point number, then it is calculated as
-!     follows:
-!
-!         X = (-1)**sign bit 1.fraction * 2**(biased exponent-bias)
-!
-!     There are several exceptions. Let a fraction, biased exponent
-!     and sign bit be denoted as F, E and S, respectively. The following
-!     formats adhere to IEEE standard:
-!
-!     S = 0, E = 00000000 and F  = 00 ... 0 : X = 0
-!     S = 0, E = 00000000 and F <> 00 ... 0 : X = +0.fraction * 2**(1-bias)
-!     S = 1, E = 00000000 and F <> 00 ... 0 : X = -0.fraction * 2**(1-bias)
-!     S = 0, E = 11111111 and F  = 00 ... 0 : X = +Inf
-!     S = 1, E = 11111111 and F  = 00 ... 0 : X = -Inf
-!     S = 0, E = 11111111 and F <> 00 ... 0 : X = NaN
-!
-!     A NaN (Not a Number) is a value reserved for signalling an
-!     attempted invalid operation, like 0/0. Its representation
-!     equals the representation of +Inf plus 1, i.e. 2**31 - 2**23 + 1
-!
-!     An example:
-!
-!        the 32-bit representation of 23.1 is
-!
-!        decimal    :    65      184      204      205
-!        hexidecimal:    41       B8       CC       CD
-!        binary     : 01000001 10111000 11001100 11001101
-!
-!        since,
-!
-!        23.1 = 2^4 + 2^2 + 2^1 + 2^0 + 2^-4 + 2^-5 + 2^-8 + 2^-9 +
-!               2^-12 + 2^-13 + 2^-16 + 2^-17 + 2^-19
-!
-!        so that the biased exponent = 4 + 127 = 131 = 10000011 = E
-!        and the sign bit = 0 = S. The remaining of the 32-bit word is
-!        the fraction, which is
-!
-!    3 2 1 0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19
-!
-! F= 0 1 1 1  0  0  0  1  1  0  0  1  1   0   0   1   1   0   0   1   1   0   1
-!
-!  4. Argument variables
-!
-!     BVAL        a byte value as a part of the representation of
-!                 floating-point number
-!     RVAL        floating-point number
-!
-      INTEGER BVAL(4)
-      REAL    RVAL
-!
-!  6. Local variables
-!
-!     ACTEXP:     actual exponent in the representation
-!     BEXPO :     biased exponent
-!     BIAS  :     bias (for 32-bit machine, bias=127)
-!     EXPO  :     calculated exponent of floating-point number
-!     FRAC  :     fraction of floating-point number
-!     I     :     loop counter
-!     IENT  :     number of entries
-!     IPART :     i-the part of the representation
-!     IQUOT :     auxiliary integer with quotient
-!     LEADNR:     leading number of floating-point number
-!     LFRAC :     length of fraction in representation
-!           :     (for 32-bit machine, lfrac=23)
-!     RFRAC :     auxiliary real with fraction
-!
-      INTEGER ACTEXP, EXPO, I, IENT, IPART, LFRAC, BIAS, BEXPO,
-     &        LEADNR, IQUOT
-      REAL FRAC, RFRAC
-!
-! 12. Structure
-!
-!     initialise 4 parts of the representation and biased exponent
-!     determine leading number and fraction
-!     do while leading number >= 1
-!        calculate positive exponent as power of two
-!     or while fraction > 0
-!        calculate negative exponent as power of two
-!     end do
-!     compute the actual part of the representation
-!
-! 13. Source text
-!
-      SAVE IENT
-      DATA IENT/0/, LFRAC/23/, BIAS/127/
-      IF (LTRACE) CALL STRACE (IENT,'SWR2B')
-
-!     --- initialise 4 parts of the representation and biased exponent
-
-      DO 10 I = 1, 4
-         BVAL(I) = 0
- 10   CONTINUE
-      BEXPO  = -1
-
-!     --- determine leading number and fraction
-
-      IF ( ABS(RVAL).LT.1.E-7 ) THEN
-         LEADNR = 0
-         FRAC   = 0.
-      ELSE
-         LEADNR = INT(ABS(RVAL))
-         FRAC   = ABS(RVAL) - REAL(LEADNR)
-      END IF
-
- 20   IF ( LEADNR.GE.1 ) THEN
-
-!        --- calculate positive exponent as power of two
-
-         EXPO  = 0
-         IQUOT = LEADNR
- 30      IF ( IQUOT.GE.2 ) THEN
-
-              IQUOT = INT(IQUOT/2)
-              EXPO  = EXPO + 1
-
-         GOTO 30
-         END IF
-
-      ELSE IF ( FRAC.GT.0. ) THEN
-
-!        --- calculate negative exponent as power of two
-
-         EXPO = 0
-         RFRAC = FRAC
- 40      IF ( RFRAC.LT.1. ) THEN
-
-              RFRAC = RFRAC * 2.
-              EXPO  = EXPO - 1
-
-         GOTO 40
-         END IF
-
-      ELSE
-
-         GOTO 50
-
-      END IF
-
-!     --- compute the actual part of the representation
-
-      IF ( BEXPO.EQ.-1 ) THEN
-
-!        --- determine biased exponent
-
-         BEXPO = EXPO + BIAS
-
-!        --- the first seven bits of biased exponent belong
-!            to first part of the representation
-
-         BVAL(1) = INT(BEXPO/2)
-
-!        --- determine the sign bit
-
-         IF ( RVAL.LT.0. ) BVAL(1) = BVAL(1) + 128
-
-!        --- the eighth bit of biased component is the leftmost
-!            bit of second part of the representation
-
-         BVAL(2) = MOD(BEXPO,2)*2**7
-         IPART = 2
-
-      ELSE
-
-!        --- compute the actual exponent of bit 1 in i-th part of
-!            the representation
-
-         ACTEXP = (IPART-2)*8 + 7 - BEXPO + BIAS + EXPO
-         IF ( ACTEXP.LT.0 ) THEN
-            ACTEXP = ACTEXP + 8
-            IPART = IPART + 1
-            IF ( IPART.GT.4 ) GOTO 50
-         END IF
-         BVAL(IPART) = BVAL(IPART) + 2**ACTEXP
-
-      END IF
-
-      IF ( EXPO.LT.(BEXPO-BIAS-LFRAC) ) GOTO 50
-      LEADNR = LEADNR - 2.**EXPO
-      IF ( EXPO.LT.0 ) FRAC = FRAC - 2.**EXPO
-
-      GOTO 20
- 50   CONTINUE
-
-      RETURN
-      END
--- SwanThreadBounds.ftn90	1970-01-01 01:00:00.000000000 +0100
+++ SwanThreadBounds.ftn90	2009-08-24 13:06:47.000000000 +0200
@@ -0,0 +1,160 @@
+subroutine SwanThreadBounds ( nwetp, ivlow, ivup, tlist )
+!
+!   --|-----------------------------------------------------------|--
+!     | Delft University of Technology                            |
+!     | Faculty of Civil Engineering and Geosciences              |
+!     | Environmental Fluid Mechanics Section                     |
+!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
+!     |                                                           |
+!     | Programmer: Marcel Zijlema                                |
+!   --|-----------------------------------------------------------|--
+!
+!
+!     SWAN (Simulating WAves Nearshore); a third generation wave model
+!     Copyright (C) 2009  Delft University of Technology
+!
+!     This program is free software; you can redistribute it and/or
+!     modify it under the terms of the GNU General Public License as
+!     published by the Free Software Foundation; either version 2 of
+!     the License, or (at your option) any later version.
+!
+!     This program is distributed in the hope that it will be useful,
+!     but WITHOUT ANY WARRANTY; without even the implied warranty of
+!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
+!     GNU General Public License for more details.
+!
+!     A copy of the GNU General Public License is available at
+!     http://www.gnu.org/copyleft/gpl.html#SEC3
+!     or by writing to the Free Software Foundation, Inc.,
+!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
+!
+!
+!   Authors
+!
+!   41.10: Marcel Zijlema
+!
+!   Updates
+!
+!   41.10, August 2009: New subroutine
+!
+!   Purpose
+!
+!   Determines load-balanced loop bounds for calling thread
+!
+!   Modules used
+!
+    use ocpcomm4
+    use SwanGriddata
+    use SwanGridobjects
+    use SwanCompdata
+!
+    implicit none
+!
+!   Argument variables
+!
+    integer, intent(out)                    :: ivlow     ! lower index in range of vertices in calling thread
+    integer, intent(out)                    :: ivup      ! upper index in range of vertices in calling thread
+    integer, intent(in)                     :: nwetp     ! total number of active vertices
+    integer, dimension(nverts), intent(out) :: tlist     ! vertex list for calling thread
+!
+!   Local variables
+!
+    integer                                 :: i         ! counter
+    integer, save                           :: ient = 0  ! number of entries in this subroutine
+    integer                                 :: ith       ! thread counter
+    integer                                 :: ivert     ! vertex index
+    integer                                 :: j         ! counter
+    integer                                 :: kvert     ! loop counter over vertices
+    integer, dimension(10)                  :: nacvt     ! number of active vertices for i-th thread
+    integer                                 :: ncurvt    ! number of currently assigned vertices to a thread
+    integer                                 :: nvcum     ! cumulative number of vertices
+    integer                                 :: nth       ! number of threads
+    integer                                 :: tid       ! thread number
+    !
+    type(verttype), dimension(:), pointer   :: vert      ! datastructure for vertices with their attributes
+    !
+!$  integer, external                       :: omp_get_num_threads ! number of OpenMP threads being used
+!$  integer, external                       :: omp_get_thread_num  ! get thread number
+!
+!   Structure
+!
+!   Description of the pseudo code
+!
+!   Source text
+!
+    if (ltrace) call strace (ient,'SwanThreadBounds')
+    !
+    ! point to vertex object
+    !
+    vert => gridobject%vert_grid
+    !
+    ! determine number of threads and thread number
+    !
+    nth = 1
+    tid = 0
+    !$ nth = omp_get_num_threads()
+    !$ tid = omp_get_thread_num()
+    tid = tid + 1
+    !
+    ! determine load-balanced sizes for all threads
+    !
+    nvcum = 0
+    do i = 1, nth
+       nacvt(i) = (nwetp*i)/nth - nvcum
+       nvcum    = (nwetp*i)/nth
+    enddo
+    !
+    ! determine loop bounds for calling thread
+    !
+    ivlow = nverts+1
+    ivup  = 0
+    !
+    ith    = 1
+    ncurvt = 0
+    !
+    do kvert = 1, nverts
+       !
+       ivert = vlist(kvert)
+       !
+       if ( vert(ivert)%active ) then
+          !
+          if ( ith == tid ) then
+             ivlow = min(kvert,ivlow)
+             ivup  = max(kvert,ivup )
+          endif
+          ncurvt = ncurvt + 1
+          !
+          if ( ncurvt >= nacvt(ith) ) then
+             ith    = ith + 1
+             ncurvt = 0
+          endif
+          !
+       endif
+       !
+    enddo
+    !
+    ! determine vertex list for calling thread based on vlist
+    ! first active vertices followed by inactive ones
+    !
+    i = ivlow
+    j = ivup
+    !
+    do kvert = ivlow, ivup
+       !
+       ivert = vlist(kvert)
+       !
+       if ( vert(ivert)%active ) then
+          !
+          tlist(i) = ivert
+          i = i + 1
+          !
+       else
+          !
+          tlist(j) = ivert
+          j = j - 1
+          !
+       endif
+       !
+    enddo
+    !
+end subroutine SwanThreadBounds
--- swanuse.tex	2009-08-24 13:08:40.000000000 +0200
+++ swanuse.tex	2009-08-24 13:11:58.000000000 +0200
@@ -32,7 +32,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72ABC}
+{\Large\bf SWAN Cycle III version 40.72ABCD}
 \end{center}
 
 \cleardoublepage
@@ -2246,13 +2246,13 @@
 {\tt OPEN}      \> the curve on which the boundary condition is given, is not closed.\\
 {\tt [xgc]}     \> if SWAN is used with \underline{Cartesian coordinates}:\+\\
                    longitude of south-west corner of SWAN computational grid (in degrees); if the\\
-                   south-west corner of the nest in the WAM computation is on land this value is\\
+                   south-west corner of the nest in the WW3 computation is on land this value is\\
                    required.\\
                    If SWAN is used with \underline{spherical coordinates} then {\tt [xgc]} is ignored by SWAN.\\
                    Default: the location of the first spectrum encountered in the nest file.\-\\
 {\tt [ygc]}     \> if SWAN is used with \underline{Cartesian coordinates}:\+\\
                    longitude of south-west corner of SWAN computational grid (in degrees); if the\\
-                   south-west corner of the nest in the WAM computation is on land this value is\\
+                   south-west corner of the nest in the WW3 computation is on land this value is\\
                    required.\\
                    If SWAN is used with \underline{spherical coordinates} then {\tt [ygc]} is ignored by SWAN.\\
                    Default: the location of the first spectrum encountered in the nest file.\-\\
@@ -2739,7 +2739,7 @@
 the domain, e.g. nearshore coasts and estuaries) in contrast to closed basin, e.g. lakes, where this option
 should not be used.
 Note that set-up is not computed correctly with spherical coordinates.
-Note that set-up is not supported in case of parallel runs using MPI and also not tested with OpenMP!
+Note that set-up is not supported in case of parallel runs using either MPI or OpenMP!
 \begin{tabbing}
 xxxxxxxxxxxx\= \kill
 {\tt [supcor]} \> by default the wave-induced set-up is computed with a constant added such that the\+\\
@@ -3335,8 +3335,7 @@
 \begin{itemize}
    \item SWAN accepts at the boundaries of the computational grid only the user-imposed incoming
          wave components and it replaces the user-imposed outgoing wave components with computed
-         components (propagating to the boundary from the interior region). Note that this will
-         not happen in case of triangular meshes.
+         components (propagating to the boundary from the interior region).
    \item during the computation of the parameters, SWAN adds an analytical (diagnostic) high-frequency
          tail to the discrete spectrum.
    \item SWAN has an option to only compute within a pre-set directional sector (pre-set by the user).
@@ -3613,7 +3612,9 @@
                    {\tt HEADER} is the {\tt PRINT} file. In case of {\tt NOHEADER} the filename is required. Note\\
                    that when the extension is `{\bf .mat}', a binary MATLAB file will be generated\\
                    automatically. This file requires less space on your computer and can be loaded\\
-                   in MATLAB much faster than the ASCII-file.\\
+                   in MATLAB much faster than an ASCII file. Also note that the output\\
+                   parameters are stored as single precision. (Hence, use the Matlab command\\
+                   {\small \tt double} for conversion to double precision, if necessary.)\\
                    Binary MATLAB files are particularly useful for the computation with\\
                    unstructured grids. Some MATLAB scripts are provided with the SWAN\\
                    source code that can be used to plot wave parameters as maps in a\\
@@ -3727,7 +3728,7 @@
                 Output only if command {\tt FRICTION} is used. If one wants to output {\tt UBOT} but\\
                 friction is ignored in the computation, then one should use the command\\
                 {\tt FRICTION} with the value of the friction set to zero ({\tt FRICTION COLLINS 0}).\-\\
-{\tt URMS}   \> the rms-value of the of the orbital velocity near the bottom (in m/s).\+\\
+{\tt URMS}   \> the rms-value of the orbital velocity near the bottom (in m/s).\+\\
                 If one wants to output {\tt URMS} but friction is ignored in the computation,\\
                 then one should use the command {\tt FRICTION} with the value of the friction\\
                 set to zero ({\tt FRICTION COLLINS 0}).\-\\
@@ -4795,18 +4796,18 @@
 !                              (SPEC  [ndec])
 !
 !   BLOCK   'sname'  HEADER | NOHEADER  'fname' (LAY-OUT [idla])            &
-!          <  DSPR|HSIGN|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DEPTH|VEL|  &
-!             FRCOEFF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS| &
-!             DHSIGN|DRTM01|LEAK|TSEC|XP|YP|DIST|SETUP|TMM10|RTMM10|        &
+!          <  HSIGN|HSWELL|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DSPR|VEL| &
+!             FRCOEF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS|  &
+!             DHSIGN|DRTM01|LEAK|TSEC|XP|YP|DIST|SETUP|TMM10|RTMM10|DEPTH|  &
 !             TMBOT|QP|BFI|WATLEV|BOTLEV|TPS|DISBOT|DISSURF|DISWCAP|        &
 !             GENE|GENW|REDI|REDQ|REDT|PROPA|PROPX|PROPT|PROPS|RADS >       &
 !             ([unit]) (OUTPUT [tbegblk] [deltblk] SEC|MIN|HR|DAY)
 !
 !   TABLE   'sname'  HEADER | NOHEADER | INDEXED 'fname'                    &
-!          <  DSPR|HSIGN|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DEPTH|VEL|  &
-!             FRCOEFF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS| &
+!          <  HSIGN|HSWELL|DIR|PDIR|TDIR|TM01|RTM01|RTP|TM02|FSPR|DSPR|VEL| &
+!             FRCOEF|WIND|DISSIP|QB|TRANSP|FORCE|UBOT|URMS|WLEN|STEEPNESS|  &
 !             DHSIGN|DRTM01|LEAK|TIME|TSEC|XP|YP|DIST|SETUP|TMM10|RTMM10|   &
-!             TMBOT|QP|BFI|WATLEV|BOTLEV|TPS|DISBOT|DISSURF|DISWCAP|        &
+!             DEPTH|TMBOT|QP|BFI|WATLEV|BOTLEV|TPS|DISBOT|DISSURF|DISWCAP|  &
 !             GENE|GENW|REDI|REDQ|REDT|PROPA|PROPX|PROPT|PROPS|RADS >       &
 !             ([unit]) (OUTPUT [tbegtbl] [delttbl] SEC|MIN|HR|DAY)
 !
@@ -4856,7 +4857,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72ABC
+$ Data produced by SWAN version 40.72ABCD
 $ Project:'projname'     ;   run number:'runnum'
 TIME                                    time-dependent data
      1                                  time coding option
@@ -4962,7 +4963,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72ABC
+$ Data produced by SWAN version 40.72ABCD
 $ Project:'projname'     ;   run number:'runnum'
 LOCATIONS                               locations in x-y-space
      1                                  number of locations
@@ -5142,7 +5143,7 @@
 
 \bibitem{Implman}
 {SWAN -- Implementation manual}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.swan.tudelft.nl} (Version 40.72ABC, August 2009).
+\hl{http://www.swan.tudelft.nl} (Version 40.72ABCD, August 2009).
 
 \bibitem{Progrul}
 {SWAN -- Programming rules}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
@@ -5206,7 +5207,7 @@
 		\hyperpage{26, 27}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38}, \hyperpage{45}, \hyperpage{52}, 
 		\hyperpage{65}, \hyperpage{67}, \hyperpage{79, 80}, 
-		\hyperpage{83}, \hyperpage{90}
+		\hyperpage{84}, \hyperpage{90}
   \item CURVE, \hyperpage{68}
   \item curvi-linear, \hyperpage{3}, \hyperpage{6}, \hyperpage{8, 9}, 
 		\hyperpage{13, 14}, \hyperpage{16}, \hyperpage{21, 22}, 
@@ -5236,7 +5237,7 @@
 		\hyperpage{14, 15}, \hyperpage{26}, \hyperpage{30, 31}, 
 		\hyperpage{42}, \hyperpage{45}, \hyperpage{51, 52}, 
 		\hyperpage{54--56}, \hyperpage{62}, \hyperpage{64--66}, 
-		\hyperpage{72--74}, \hyperpage{79}, \hyperpage{83}, 
+		\hyperpage{72--74}, \hyperpage{79}, \hyperpage{83, 84}, 
 		\hyperpage{86}, \hyperpage{89}, \hyperpage{91}, 
 		\hyperpage{108, 109}
   \item FRICTION, \hyperpage{55}
@@ -5371,7 +5372,7 @@
   \item TRIAD, \hyperpage{56}
   \item triads, \hyperpage{16}, \hyperpage{56}, \hyperpage{80}, 
 		\hyperpage{92}
-  \item triangular, \hyperpage{3}, \hyperpage{16}, \hyperpage{72}
+  \item triangular, \hyperpage{3}, \hyperpage{16}
 
   \indexspace
 
@@ -5385,8 +5386,8 @@
 
   \item WAM, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{18}, \hyperpage{22}, \hyperpage{37}, 
-		\hyperpage{46}, \hyperpage{48--50}, \hyperpage{54}, 
-		\hyperpage{81, 82}, \hyperpage{84, 85}, \hyperpage{87}
+		\hyperpage{46}, \hyperpage{48, 49}, \hyperpage{54}, 
+		\hyperpage{81}, \hyperpage{83--85}, \hyperpage{87}
   \item WAVEWATCH, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{22}, \hyperpage{46}, \hyperpage{49, 50}
   \item whitecapping, \hyperpage{7}, \hyperpage{16}, \hyperpage{53}, 
--- switch.pl	2008-05-25 14:59:38.000000000 +0200
+++ switch.pl	2009-08-24 13:06:37.000000000 +0200
@@ -3,7 +3,6 @@
 $tim = "FALSE";
 $mpi = "FALSE";
 $f95 = "FALSE";
-$omp = "FALSE";
 $dos = "FALSE";
 $unx = "FALSE";
 $cry = "FALSE";
@@ -16,7 +15,6 @@
    if ($ARGV[0]=~/-timg/) {$tim="TRUE";shift;}
    if ($ARGV[0]=~/-mpi/) {$mpi="TRUE";shift;}
    if ($ARGV[0]=~/-f95/) {$f95="TRUE";shift;}
-   if ($ARGV[0]=~/-omp/) {$omp="TRUE";shift;}
    if ($ARGV[0]=~/-dos/) {$dos="TRUE";shift;}
    if ($ARGV[0]=~/-unix/) {$unx="TRUE";shift;}
    if ($ARGV[0]=~/-cray/) {$cry="TRUE";shift;}
@@ -35,13 +33,7 @@
 foreach $file (@files)
 {
 # --- set output file name
-  if ($unx=~/TRUE/ && $omp=~/TRUE/)
-  {
-    ($tempf)=split(/.ftn/, $file);
-    $ext = ($file =~ m/ftn90/) ? "F90" : "F";
-    $outfile = join(".",$tempf,$ext);
-  }
-  elsif ($unx=~/TRUE/)
+  if ($unx=~/TRUE/)
   {
     ($tempf)=split(/.ftn/, $file);
     $ext = ($file =~ m/ftn90/) ? "f90" : "f";
@@ -68,7 +60,6 @@
       if ($tim=~/TRUE/) {$newline=~s/^!TIMG//;}
       if ($mpi=~/TRUE/) {$newline=~s/^!MPI//;}
       if ($f95=~/TRUE/) {$newline=~s/^!F95//;}
-      if ($omp=~/TRUE/) {$newline=~s/^!OMP//;}
       if ($dos=~/TRUE/) {$newline=~s/^!DOS//;}
       if ($unx=~/TRUE/) {$newline=~s/^!UNIX//;}
       if ($cry=~/TRUE/) {$newline=~s/^!\/Cray//;}
--- swmod1.ftn	2009-08-24 13:08:40.000000000 +0200
+++ swmod1.ftn	2009-08-24 13:06:49.000000000 +0200
@@ -2030,6 +2030,7 @@
 !                 =5; for command GEN3 YAN ..., if wind is made active,
 !                 =5; for command GROWTH G3 YAN ..., not documented in manual,
 !                 not documented in manual
+! LADDS  [.fals.] indicates whether extra output is requested or not
 ! LSETUP [     0] =0; setup is not calculated
 !                 =1; setup is calculated
 !                 =2; setup is calculated with the boundary conditions from
@@ -2338,6 +2339,7 @@
       LOGICAL             ACUPDA
       LOGICAL             BNDCHK,      BNAUT,       ONED,        BRESCL
       LOGICAL             OFFSRC                                          40.80
+      LOGICAL             LADDS                                           40.85
       LOGICAL             CSETUP
 !$OMP THREADPRIVATE(ICMAX,CSETUP)
 !
