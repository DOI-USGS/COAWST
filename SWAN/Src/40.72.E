1) Some adaption to Makefile and perl script plaform.pl.
2) The so-called CSM formulation for whitecapping is removed.
3) Reading input fields and block outputting for unstructured mesh
   cases is made more efficient.
4) Technical documentation is extended with useful information.
5) Bug fix in calculation of the orbital velocities.

--- Makefile	2009-12-17 10:30:36.000000000 +0100
+++ Makefile	2009-12-17 10:31:22.000000000 +0100
@@ -88,9 +88,9 @@
 HCAT_EXE = hcat.exe
 HCAT_OBJS = swanhcat.$(EXTO)
 
-.SUFFIXES: .f .for .f90
+.SUFFIXES: .o .f .for .f90
 
-.PHONEY: help
+.PHONEY: help clean clobber
 
 help:
 	@echo "This Makefile supports the following:"
@@ -112,20 +112,28 @@
 ser:
 	@perl switch.pl $(swch) *.ftn *.ftn90
 	$(MAKE) FOR=$(F90_SER) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
+	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_SER)" \
                 INCS="$(INCS_SER)" LIBS="$(LIBS_SER)" $(SWAN_EXE)
-	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" $(HCAT_EXE)
+	$(MAKE) hcat
 
 omp:
 	@perl switch.pl $(swch) *.ftn *.ftn90
 	$(MAKE) FOR=$(F90_OMP) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_OMP)" \
+	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_OMP)" \
                 INCS="$(INCS_OMP)" LIBS="$(LIBS_OMP)" $(SWAN_EXE)
-	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" $(HCAT_EXE)
+	$(MAKE) hcat
 
 mpi:
 	@perl switch.pl $(swch) -mpi *.ftn *.ftn90
 	$(MAKE) FOR=$(F90_MPI) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI)" \
+	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_MPI)" \
                 INCS="$(INCS_MPI)" LIBS="$(LIBS_MPI)" $(SWAN_EXE)
-	$(MAKE) FOR="$(F90_SER)" FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" $(HCAT_EXE)
+	$(MAKE) hcat
+
+hcat:
+	@perl switch.pl $(swch) swanhcat.ftn
+	$(MAKE) FOR=$(F90_SER) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
+	        FFLAGS90="$(FLAGS_OPT) $(FLAGS90_MSC) $(FLAGS_SER)" $(HCAT_EXE)
 
 doc:
 	$(MAKE) -f Makefile.latex TARGET=swanuse doc
@@ -144,7 +152,7 @@
 	$(FOR) $< -c $(FFLAGS) $(INCS)
 
 .f90.o:
-	$(FOR) $< -c $(FFLAGS) $(INCS)
+	$(FOR) $< -c $(FFLAGS90) $(INCS)
 
 .for.o:
 	$(FOR) $< -c $(FFLAGS) $(INCS)
@@ -153,7 +161,7 @@
 	$(FOR) $< -c $(FFLAGS) $(INCS)
 
 .f90.obj:
-	$(FOR) $< -c $(FFLAGS) $(INCS)
+	$(FOR) $< -c $(FFLAGS90) $(INCS)
 
 clean:
 	$(RM) *.$(EXTO) *.mod
--- ocpmix.ftn	2009-12-17 10:30:36.000000000 +0100
+++ ocpmix.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -1837,6 +1837,7 @@
                CALL MSGERR (4,
      &         'File cannot be opened/does not exist: '//DDNAME_L)        40.03
                IOSTAT = IEEXST
+               RETURN
             END IF
             IF (OPENED) THEN
                IF (IOSTAT.GT.-1)
--- platform.pl	2009-12-17 10:30:36.000000000 +0100
+++ platform.pl	2009-12-17 10:31:22.000000000 +0100
@@ -31,6 +31,7 @@
   print OUTFILE "F90_MPI = f90\n";
   print OUTFILE "FLAGS_OPT = -Ofast=IP35 -mips4 -r12000\n";
   print OUTFILE "FLAGS_MSC = -64\n";
+  print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
   print OUTFILE "FLAGS_SER =\n";
   print OUTFILE "FLAGS_OMP = -mp\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -56,6 +57,7 @@
   print OUTFILE "FLAGS_OPT = -O3 -qstrict -qarch=auto -qtune=auto -qnohot -qcache=auto \\\n";
   print OUTFILE "            -qunroll -qalign=4k -qfloat=hsflt\n";
   print OUTFILE "FLAGS_MSC = -w -qfixed -qnosave -q64\n";
+  print OUTFILE "FLAGS90_MSC = -w -qfree=f90 -qnosave -q64\n";
   print OUTFILE "FLAGS_SER =\n";
   print OUTFILE "FLAGS_OMP = -qsmp=omp\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -80,6 +82,7 @@
   print OUTFILE "F90_MPI = f90\n";
   print OUTFILE "FLAGS_OPT = -fast\n";
   print OUTFILE "FLAGS_MSC = -w -fixed\n";
+  print OUTFILE "FLAGS90_MSC = -w -free\n";
   print OUTFILE "FLAGS_SER =\n";
   print OUTFILE "FLAGS_OMP = -Wp,-C -omp\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -104,6 +107,7 @@
   print OUTFILE "F90_MPI = mpf90\n";
   print OUTFILE "FLAGS_OPT = -xO3 -xtarget=native -fsimple=1 -depend -libmil -xlibmopt -xlic_lib=sunperf\n";
   print OUTFILE "FLAGS_MSC = -w -silent\n";
+  print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
   print OUTFILE "FLAGS_SER =\n";
   print OUTFILE "FLAGS_OMP = -openmp\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -128,6 +132,7 @@
   print OUTFILE "F90_MPI =\n";
   print OUTFILE "FLAGS_OPT = +O2 +Onolimit\n";
   print OUTFILE "FLAGS_MSC =\n";
+  print OUTFILE "FLAGS90_MSC =\n";
   print OUTFILE "FLAGS_SER = +Onoopenmp\n";
   print OUTFILE "FLAGS_OMP = +Oopenmp\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -164,6 +169,7 @@
       print OUTFILE "FLAGS_OPT = -O2\n";
 #    }
     print OUTFILE "FLAGS_MSC = -W0 -assume byterecl -traceback\n";
+    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP = -openmp\n";
     print OUTFILE "FLAGS_MPI =\n";
@@ -177,7 +183,7 @@
     print OUTFILE "EXTO = o\n";
     print OUTFILE "MAKE = make\n";
     print OUTFILE "RM = rm -f\n";
-    print OUTFILE "swch = -unix -f95 -timg -impi\n";
+    print OUTFILE "swch = -unix -impi\n";
   }
   elsif ( -f "ifc" )
   {
@@ -191,6 +197,7 @@
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -O2 -tpp7\n";
     print OUTFILE "FLAGS_MSC = -W0 -auto\n";
+    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP = -openmp\n";
     print OUTFILE "FLAGS_MPI =\n";
@@ -218,6 +225,7 @@
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -O2 -tpp1\n";
     print OUTFILE "FLAGS_MSC = -W0 -auto\n";
+    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP = -openmp\n";
     print OUTFILE "FLAGS_MPI =\n";
@@ -243,10 +251,11 @@
     print OUTFILE "F90_OMP = pgf90\n";
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -fast\n";
-    print OUTFILE "FLAGS_MSC =\n";
+    print OUTFILE "FLAGS_MSC = -Mfixed\n";
+    print OUTFILE "FLAGS90_MSC = -Mfree\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP = -mp\n";
-    print OUTFILE "FLAGS_MPI =\n";
+    print OUTFILE "FLAGS_MPI = -tp barcelona-64\n";
     print OUTFILE "INCS_SER =\n";
     print OUTFILE "INCS_OMP =\n";
     print OUTFILE "INCS_MPI =\n";
@@ -257,7 +266,7 @@
     print OUTFILE "EXTO = o\n";
     print OUTFILE "MAKE = make\n";
     print OUTFILE "RM = rm -f\n";
-    print OUTFILE "swch = -unix -f95 -impi -timg\n";
+    print OUTFILE "swch = -unix -impi\n";
   }
   elsif ( -f "lf95" )
   {
@@ -270,6 +279,7 @@
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -O --tpp\n";
     print OUTFILE "FLAGS_MSC = --staticlink --nwo\n";
+    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP =\n";
     print OUTFILE "FLAGS_MPI =\n";
@@ -283,7 +293,7 @@
     print OUTFILE "EXTO = o\n";
     print OUTFILE "MAKE = make\n";
     print OUTFILE "RM = rm -f\n";
-    print OUTFILE "swch = -unix -f95 -timg\n";
+    print OUTFILE "swch = -unix\n";
   }
   elsif ( -f "gfortran" )
   {
@@ -292,10 +302,11 @@
     print OUTFILE "# IA32_GNU:		Intel Pentium with Linux using GNU compiler gfortran.\n";
     print OUTFILE "##############################################################################\n";
     print OUTFILE "F90_SER = gfortran\n";
-    print OUTFILE "F90_OMP = \n";
+    print OUTFILE "F90_OMP = gfortran\n";
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -O\n";
-    print OUTFILE "FLAGS_MSC = -w -ffree-line-length-none\n";
+    print OUTFILE "FLAGS_MSC = -w -fno-second-underscore\n";
+    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC) -ffree-line-length-none\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP = -fopenmp\n";
     print OUTFILE "FLAGS_MPI =\n";
@@ -309,7 +320,7 @@
     print OUTFILE "EXTO = o\n";
     print OUTFILE "MAKE = make\n";
     print OUTFILE "RM = rm -f\n";
-    print OUTFILE "swch = -unix -f95 -timg\n";
+    print OUTFILE "swch = -unix\n";
   }
   elsif ( -f "g95" )
   {
@@ -318,10 +329,11 @@
     print OUTFILE "# IA32_GNU:		Intel Pentium with Linux using GNU compiler g95.\n";
     print OUTFILE "##############################################################################\n";
     print OUTFILE "F90_SER = g95\n";
-    print OUTFILE "F90_OMP = \n";
+    print OUTFILE "F90_OMP = NO_OPENMP_WITH_G95\n";
     print OUTFILE "F90_MPI = mpif90\n";
     print OUTFILE "FLAGS_OPT = -O\n";
-    print OUTFILE "FLAGS_MSC = -w\n";
+    print OUTFILE "FLAGS_MSC = -fno-second-underscore\n";
+    print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC) -ffree-line-length-huge\n";
     print OUTFILE "FLAGS_SER =\n";
     print OUTFILE "FLAGS_OMP =\n";
     print OUTFILE "FLAGS_MPI =\n";
@@ -335,7 +347,34 @@
     print OUTFILE "EXTO = o\n";
     print OUTFILE "MAKE = make\n";
     print OUTFILE "RM = rm -f\n";
-    print OUTFILE "swch = -unix -f95 -timg\n";
+    print OUTFILE "swch = -unix\n";
+  }
+  elsif ( -f "xlf90" )
+  {
+    system 'rm xlf90';
+    print OUTFILE "##############################################################################\n";
+    print OUTFILE "# PPC64_IBM:		IBM Power6 with Linux using IBM Compiler.\n";
+    print OUTFILE "##############################################################################\n";
+    print OUTFILE "F90_SER = xlf90\n";
+    print OUTFILE "F90_OMP =\n";
+    print OUTFILE "F90_MPI = mpfort\n";
+    print OUTFILE "FLAGS_OPT = -O3 -qstrict -qarch=auto -qtune=auto\n";
+    print OUTFILE "FLAGS_MSC = -qfixed -qzerosize -qwarn64\n";
+    print OUTFILE "FLAGS90_MSC = -qfree=f90 -qzerosize -qwarn64\n";
+    print OUTFILE "FLAGS_SER =\n";
+    print OUTFILE "FLAGS_OMP =\n";
+    print OUTFILE "FLAGS_MPI =\n";
+    print OUTFILE "INCS_SER = \n";
+    print OUTFILE "INCS_OMP =\n";
+    print OUTFILE "INCS_MPI =\n";
+    print OUTFILE "LIBS_SER =\n";
+    print OUTFILE "LIBS_OMP =\n";
+    print OUTFILE "LIBS_MPI =\n";
+    print OUTFILE "OUT = -o \n";
+    print OUTFILE "EXTO = o\n";
+    print OUTFILE "MAKE = make\n";
+    print OUTFILE "RM = rm -f\n";
+    print OUTFILE "swch = -unix\n";
   }
   else
   {
@@ -351,6 +390,7 @@
   print OUTFILE "F90_MPI = ifort\n";
   print OUTFILE "FLAGS_OPT = /optimize:2\n";
   print OUTFILE "FLAGS_MSC = /assume:byterecl /traceback /nowarn /nologo\n";
+  print OUTFILE "FLAGS90_MSC = \$(FLAGS_MSC)\n";
   print OUTFILE "FLAGS_SER =\n";
   print OUTFILE "FLAGS_OMP = /Qopenmp\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -375,6 +415,7 @@
   print OUTFILE "F90_MPI =\n";
   print OUTFILE "FLAGS_OPT = -O3 -qstrict -qtune=auto -qcache=auto -qalign=4k\n";
   print OUTFILE "FLAGS_MSC = -w -qfixed\n";
+  print OUTFILE "FLAGS90_MSC = -w -qfree=f90\n";
   print OUTFILE "FLAGS_SER =\n";
   print OUTFILE "FLAGS_OMP =\n";
   print OUTFILE "FLAGS_MPI =\n";
@@ -388,7 +429,7 @@
   print OUTFILE "EXTO = o\n";
   print OUTFILE "MAKE = make\n";
   print OUTFILE "RM = rm -f\n";
-  print OUTFILE "swch = -unix -f95 -timg\n";
+  print OUTFILE "swch = -unix\n";
 }
 else
 {
--- swancom1.ftn	2009-12-17 10:30:36.000000000 +0100
+++ swancom1.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -1228,7 +1228,9 @@
 !       initialise Ursell number to 0 for each iteration                  40.03
 !       this is done in parallel within OpenMP environment                40.31
 !
-        IF (ITRIAD.GT.0 .OR. ISURF.EQ.5) THEN                             41.03
+        IF (ITRIAD.GT.0
+     &      .OR. ISURF.EQ.5                                               41.03
+     &                 ) THEN
            DO IP = I1GRD,I2GRD                                            40.31
             COMPDA(IP,JURSEL) = 0.
           ENDDO
@@ -1656,7 +1658,7 @@
      &            CAX1             ,CAY1                                  40.22
      &                                                                 )
 !TIMG                CALL SWTSTO(104)                                          40.23
-                IF (STPNOW()) RETURN                                      34.01
+!MPI                IF (STPNOW()) RETURN                                      34.01
 
 !----------------------------------------------------------------------   40.22
 !               Once the computation is done for grid point (IX,IY) the   40.22
@@ -2949,6 +2951,7 @@
      &  MEMNL4              ,WWINT               ,WWAWG               ,
      &  WWSWG               ,CGO                 ,COMPDA(1,JUSTAR)    ,
      &  COMPDA(1,JZEL)      ,SPCDIR              ,ANYWND              ,
+     &  CAS                 ,                                             41.11
      &  SWMATR(1,1,JDIS0)   ,SWMATR(1,1,JDIS1)   ,
      &  SWMATR(1,1,JGEN0)   ,SWMATR(1,1,JGEN1)   ,                        40.85
      &  SWMATR(1,1,JRED0)   ,SWMATR(1,1,JRED1)   ,                        40.85
@@ -3022,7 +3025,6 @@
      &                   INOCNV, IDDLOW, IDDTOP, ISSTOP, IDCMIN,          40.41
      &                   IDCMAX )
 !TIMG            CALL SWTSTO(120)                                              40.23
-            IF (STPNOW()) RETURN                                          34.01
 !
           ELSE IF (INT(PNUMS(8)).EQ.2 .OR. INT(PNUMS(8)).EQ.3) THEN       40.00
 !
@@ -3040,7 +3042,6 @@
      &                    IDDTOP                                 )        5/mar
 !TIMG            CALL SWTSTO(120)                                              40.23
 !
-!
           END IF
 !
         ELSE                                                              40.00
@@ -3406,10 +3407,6 @@
          WRITE(PRINTF,6136) PWCAP(6), PWCAP(7)
          WRITE(PRINTF,6137) PWCAP(8)
  6137    FORMAT('                      : KCONV  ',E12.4)
-      ELSEIF (IWCAP.EQ.6) THEN                                            40.30
-         WRITE(PRINTF,6138) PWCAP(12), PWCAP(13)                          40.30
- 6138    FORMAT(' W-cap cum. steepness : CST    ',                        40.30
-     &           E12.4,' POW   ',E12.4)                                   40.30
       ELSEIF (IWCAP.EQ.7) THEN                                            40.53
          WRITE(PRINTF,6139) PWCAP(1), PWCAP(12)                           40.53
          WRITE(PRINTF,6140) PWCAP(9), PWCAP(11)                           40.53
@@ -4720,7 +4717,7 @@
 !     QB      : computed in the subroutine FRABRE
 !     SIGM01  = ETOT1 / ETOT
 !     SIGM_10 = ETOT  / ACTOT
-!     UBOT    = Sqrt [ UB2 ]
+!     UBOT    = Sqrt [ 2 UB2 ] Note: this is Urms according to Madsen (1994)
 !     TMBOT   = 2*PI AB2 / UB2
 !
 !  4. Argument variables
@@ -4785,7 +4782,8 @@
 !     WH                 : fraction of breaking waves
 !
       REAL              :: AB2, EMAX, UB2
-      REAL              :: BIPH, WH
+      REAL              :: WH
+      REAL              :: BIPH
       REAL              :: FRINTF_X_DDIR
       REAL              :: BRCOEF   ! variable breaking coefficient (calc. in BRKPAR)  40.22
 !
@@ -4852,7 +4850,9 @@
 !
       ETOT = ETOT + ETOT_DSIG(MSC) * PWTAIL(6) / FRINTF
 !
-      IF (ISURF .EQ. 1 .OR. ISURF.EQ.4 .OR. ISURF.EQ.5) THEN              41.03
+      IF (ISURF .EQ. 1 .OR. ISURF.EQ.4                                    41.03
+     &    .OR. ISURF.EQ.5                                                 41.03
+     &                                ) THEN
         HM   = PSURF(2) * DEP2(KCGRD(1))                                  40.22
       ELSEIF (ISURF .EQ. 2 .OR. ISURF.EQ.3) THEN                          41.03
 !        Calulate the correct breaking coefficient BRCOEF
@@ -4949,7 +4949,7 @@
 !     Calculate the orbital velocity UBOT, orbital excursion ABRBOT and
 !     bottom wave period TMBOT                                            40.51
 !
-        IF ( UB2 .GT. 0.) UBOT(KCGRD(1)) = SQRT ( UB2 )
+        IF ( UB2 .GT. 0.) UBOT(KCGRD(1)) = SQRT ( 2. * UB2 )
         IF ( AB2 .GT. 0.) ABRBOT = SQRT (2. *  AB2)
         IF ( UB2.GT.0. .AND. AB2.GT.0. )                                  40.51
      &                               TMBOT(KCGRD(1)) = PI2*SQRT(AB2/UB2)  40.51
@@ -4959,7 +4959,9 @@
 !     *** calculate Ursell number ***
 !     *** update only for first encounter in a sweep                      40.16
 !
-      IF ( ITRIAD.GT.0 .OR. ISURF.EQ.5 ) THEN                             41.03 40.41
+      IF ( ITRIAD.GT.0                                                    40.41
+     &     .OR. ISURF.EQ.5                                                41.03
+     &                 ) THEN                                             40.41
          IF (( SWPDIR .EQ. 1) .OR.                                        40.16
      &       ( SWPDIR .EQ. 2 .AND. IXCGRD(1) .EQ. 1) .OR.                 40.16
      &       ( SWPDIR .EQ. 3 .AND. IYCGRD(1) .EQ. 1) .OR.                 40.16
@@ -5825,6 +5827,7 @@
      &                     MEMNL4     ,WWINT      ,WWAWG      ,
      &                     WWSWG      ,CGO        ,USTAR      ,
      &                     ZELEN      ,SPCDIR     ,ANYWND     ,
+     &                     CAS        ,                                   41.11
      &                     DISSC0     ,DISSC1     ,GENC0      ,           40.85
      &                     GENC1      ,REDC0      ,REDC1      ,           40.85
      &                     XIS        ,FRCOEF     ,IT         ,           30.00
@@ -6168,6 +6171,7 @@
       REAL  :: WWAWG(*)
       REAL  :: WWSWG(*)
       REAL  :: CGO(MSC,MICMAX)                                            40.22
+      REAL  :: CAS(MDC,MSC,MICMAX)                                        41.11
       REAL  :: USTAR(MCGRD)
       REAL  :: ZELEN(MCGRD)
       REAL  :: DISSC0(1:MDC,1:MSC,1:MDISP)                                40.67
@@ -6233,7 +6237,7 @@
 !         *** calculate surf breaking source term (5 formulations) ***    41.03
 !
           CALL SSURF (ETOT    ,HM      ,QBLOC   ,SMEBRK  ,
-     &                SPCSIG  ,AC2     ,IMATRA  ,                         30.81
+     &                KMESPC  ,SPCSIG  ,AC2     ,IMATRA  ,                30.81
      &                IMATDA  ,IDCMIN  ,IDCMAX  ,PLWBRK  ,
      &                ISSTOP  ,DISSC0  ,DISSC1  )                         40.67 40.61 30.21
 !
@@ -6302,6 +6306,7 @@
      &                            IDCMIN  ,IDCMAX  ,ISSTOP  ,             40.02
      &                            ETOT    ,IMATDA  ,IMATRA  ,PLWCAP  ,    40.02
      &                            CGO     ,UFRIC   ,                      40.53
+     &                            CAS     ,                               41.11
      &                            DEP2    ,DISSC1  ,DISSC0  )             40.67 40.61 40.12
 !TIMG      CALL SWTSTO(133)                                                    40.23
 !
@@ -7928,9 +7933,8 @@
             IF ( IMATDA(ID,IS).NE.0. ) THEN
                INVMDA(ID,IS) = 1./IMATDA(ID,IS)
             ELSE
-               CALL MSGERR ( 4,
+               CALL MSGERR ( 3,
      &                     'Main diagonal of spectral matrix is zero!' )
-               RETURN
             END IF
          END DO
       END DO
--- swancom2.ftn	2009-12-17 10:30:34.000000000 +0100
+++ swancom2.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -612,7 +612,7 @@
 !****************************************************************
 !
       SUBROUTINE SSURF (ETOT    ,HM      ,QB      ,SMEBRK  ,              30.81
-     &                  SPCSIG  ,AC2     ,IMATRA  ,                       30.81
+     &                  KMESPC  ,SPCSIG  ,AC2     ,IMATRA  ,              30.81
      &                  IMATDA  ,IDCMIN  ,IDCMAX  ,PLWBRK  ,              30.81
      &                  ISSTOP  ,DISSC0  ,DISSC1  )                       40.67 40.61 30.81 30.21
 !
@@ -665,6 +665,7 @@
 !     40.67: Nico Booij
 !     41.03: Andre van der Westhuysen
 !     41.06: Gerbrant van Vledder
+!     41.09: Andre van der Westhuysen
 !
 !  1. Updates
 !
@@ -678,11 +679,12 @@
 !     40.67, Jun. 07: more accurate computation of dissipation terms
 !     41.03, Feb. 09: extension to alternative surf breaking formula's
 !     41.06, Mar. 09: extension to frequency dependent surf breaking
+!     41.09, Oct. 09: extension to biphase breaker model (variable n)
 !
 !  2. Purpose
 !
-!     Computation of the source term due to wave breaking with one of three
-!     formulation:
+!     Computation of the source term due to wave breaking with one of the
+!     following formulation:
 !                 1) Battjes and Janssen (1978)
 !                 2) Thornton and Guza (1983)
 !                 3) biphase scaling of Van der Westhuysen (2009)
@@ -801,7 +803,9 @@
 !
 !          after Thornton and Guza (1983), Eqs. (24),(25)
 !
-!     For the weighting function W (or fraction of breaking waves), two choices can be made.
+!     and
+!
+!     for the weighting function W (or fraction of breaking waves), two choices can be made.
 !     The first one is after Thornton and Guza (1983):
 !
 !           Hrms  n
@@ -839,6 +843,7 @@
 !     IMATDA  output:   Coefficient of diagonal matrix (2D)
 !     IMATRA  output:   Coefficient of righthandside of matrix
 !     ISSTOP  input :   Maximum for counter IS
+!     KMESPC  input :   Mean average wavenumber according to the WAM-formulation
 !     PLWBRK  output:   array containing the surf breaking source term
 !                       for test-output
 !     QB      input :   Fraction of breaking waves
@@ -855,7 +860,7 @@
      &         PLWBRK(MDC,MSC,NPTST)                                      40.00
       REAL     SPCSIG(MSC)
 !
-      REAL     ETOT,  HM,  QB, SMEBRK                                     30.81
+      REAL     ETOT,  HM,  QB, SMEBRK, KMESPC                             30.81
 !
 !  5. Parameter variables
 !
@@ -876,8 +881,8 @@
 !     SbrD    Derivative of source term for surf breaking (Sbr) to action density
 !
       INTEGER          ID,       IDDUM,   IENT,   IS
-      REAL             BIPH,     DEPLOC,  URSLOC, WH,
-     &                 PP,       FAC,     EPTOT,  ETOT0,
+      REAL             PP,       FAC,     EPTOT,  ETOT0,
+     &                 BIPH,     DEPLOC,  URSLOC, WH, SMEAN,
      &                 ECS(MDC), FMIN,    FMAX,   FRFAC(MSC)
       DOUBLE PRECISION BB,       DIS0,    SbrD,
      &                 SURFA0,   SURFA1,  WS  ,
@@ -984,6 +989,12 @@
          URSLOC = (GRAV*2.*SQRT(ETOT))/(SQRT(2.)*SMEBRK**2*DEPLOC**2)
          BIPH = 0.5*PI*(TANH(PTRIAD(4)/URSLOC)-1.)
          IF ( BIPH.LT.0. ) THEN
+!
+            IF ( PSURF(6).GT.0. ) THEN                                    41.09
+               SMEAN    = SQRT(8.*ETOT)*KMESPC/(2.*PI)                    41.09
+               PSURF(4) = 4.-(4./PI)*ATAN(PSURF(6)*(SMEAN-0.0375))        41.09
+            ENDIF                                                         41.09
+!
             WH = (BIPH/PSURF(5))**PSURF(4)
             WH = MIN(1.,WH)
 !
@@ -1067,6 +1078,7 @@
      &                   IDCMIN  ,IDCMAX  ,ISSTOP  ,                      40.02
      &                   ETOT    ,IMATDA  ,IMATRA  ,PLWCAP  ,             40.02
      &                   CGO     ,UFRIC   ,                               40.53
+     &                   CAS     ,                                        41.11
      &                   DEP2    ,DISSC1  ,DISSC0  )                      40.67 40.61 40.12
 !
 !****************************************************************
@@ -1121,20 +1133,19 @@
 !     40.61: Marcel Zijlema
 !     40.63: Andre van der Westhuysen
 !     40.67: Nico Booij
+!     41.11: Andre van der Westhuysen
 !
 !  1. Updates
 !
 !     40.02, Jan. 00: New, based on the old SWCAP1-5 subroutines
 !     40.12, Nov. 00: Added WCAP to dissipation output (bug fix)
-!     40.30, Mar. 03: cumulative steepness method (taken from Fysica+)
-!     40.31, Feb. 04: some changes/improvements to CSM
-!     40.41, May. 04: Normalisation added to CSM
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
 !     40.53: Aug. 04: white-capping based Alves and Banner (2003) method
 !     40.61, Sep. 06: introduce DISWCP variable for output purposes
 !     40.63, Apr. 07: a correction to Alves and Banner method
 !     40.67, Jun. 07: more accurate computation of source terms
 !                     DISSIP and DISIMP renamed DISSC1 and DISSC0 (as elsewhere)
+!     41.11, Oct. 09: Enhanced dissipation in counter current
 !
 !
 !  2. Purpose
@@ -1152,8 +1163,6 @@
 !     C_wc1 = C_K  sig~ (k/k~): According to Komen (generalised)
 !     C_wc2 = C_BJ sig~ (k/k~): According to Battjes-Janssen (modified)
 !     C_wc3 = C_LH            : According to Longuet Higgins (1969)
-!     C_wc4 = C_CSM           : According to cumulative steepness
-!                               method (Ris et al. (1999) and Donelan (1999))
 !
 !     In these formulations C_K is defined as (Komen; 1994 p. 145):
 !
@@ -1190,42 +1199,6 @@
 !     and C3 can be varied
 !
 !
-!     Cumulative steepness method. Basic idea is that dissipation depends
-!     on the steepness of the wave spectrum at and below a particular
-!     frequency. Details can be found in "SWAN fysica plus" (2002),
-!     report H3937/A832, implemented by Delf Hydraulics and Alkyon by
-!     order of RIKZ/RWS as a part of the project HR-Ontwikkeling.
-!                                                                         40.41
-!     Recent modifications are described in:                              40.41
-!     Hurdle, D.P., and G.Ph. van Vledder, 2004: Improved spectral wave   40.41
-!     modelling of white-capping dissipation in swell sea systems. Proc.  40.41
-!     Offshore Mechanics, Arctic Engineering Conference, June 2004,       40.41
-!     Vancouver, Canada,  paper OMAE2004-51562.                           40.41
-!                                                                         40.41
-!     Default values: C_st = 4                                            40.41
-!                     m    = 2                                            40.41
-!                                                                         40.41
-!     C_CSM = A*C_st S(sig,th)                                            40.41
-!
-!     where
-!
-!     C_st is a tuneable coefficient and
-!
-!     S(sig,th) = integrals from 0 to sig and from 0 to 2pi of
-!
-!                 k^2 |cos(th-th')|^m E(sig,th) dsig dth'
-!
-!     where coefficient m controls the directional dependence in case
-!     of straining mechanism (if dominant then m = 2 else m > 10) and
-!     th is the actual direction where straining mechanism takes
-!     place.
-!                                                                         40.41
-!     A is a normalisation coefficient computed                           40.41
-!                                                                         40.41
-!              1     GAMMA(0.5*m+1.)                                      40.41
-!       A = -------- ----------------                                     40.41
-!           SQRT(PI) GAMMA(0.5*m+0.5)                                     40.41
-!
 !     In these equations the variables have the following meaning:
 !
 !     Hm   : Maximum wave height
@@ -1281,6 +1254,7 @@
       REAL, INTENT(INOUT) :: DISSC0(MDC,MSC,MDISP),DISSC1(MDC,MSC,MDISP)  40.67 40.12
       REAL, INTENT(IN)    :: UFRIC                                        40.53
       REAL, INTENT(IN)    :: CGO(MSC,MICMAX)                              40.53
+      REAL, INTENT(IN)    :: CAS(MDC,MSC,MICMAX)                          41.11
 !
 !  6. Local variables
 !
@@ -1289,13 +1263,6 @@
 !     C_BJ  : Whitecapping coefficient according to Battjes-Janssen
 !     C_K   : Whitecapping coefficient according to Komen
 !     C_LH  : Whitecapping coefficient according to Longuet Higgins
-!     CHARS : array to pass character info to MSGERR
-!     CUMSTP: Cumulative steepness spectrum
-!     DDIF  : Directional difference
-!     DELTA : Directional step (radians)
-!     DSIGMA: Step sizes, ratio with sigma
-!     DSTEEP: Contribution to cumulative steepness
-!     EBIN  : Total energy per bin
 !     EF    : Energy density spectrum in frequency domain
 !             (azimuth-integrated frequency spectrum)
 !     HM    : Maximum waveheight as used in the Battjes-Janssen expression
@@ -1306,8 +1273,6 @@
 !     IDDUM : Counter in directional space within sweep limits
 !     IENT  : Number of entries in this subroutine
 !     IS    : Counter in frequency space
-!     MSGSTR: string to pass message to call MSGERR
-!     MXWCP : maximum number of white-capping methods
 !     N1    : Exponent for the wavenumber term in the Komen expression
 !     N2    : Exponent for the steepness term in the Komen expression
 !     P     : Exponent of the relative saturation (B/Br)
@@ -1316,29 +1281,21 @@
 !     STP_OV: Overall steepness
 !     STP_PM: Overall steepness for a Pierson-Moskowitz spectrum
 !     WCAP  : Whitecapping source-term
+!     WCCUR : Whitecapping source term due to opposing current
 !     WCIMPL: Implicit part of the whitecapping source-term
-!     XFAC  : Factor between sigma & dsigma
 !
       INTEGER, SAVE     :: IENT = 0
-      INTEGER           :: ID, IDDUM, IS, ID1, ID2, IF, IL, MXWCP
+      INTEGER           :: ID, IDDUM, IS, ID1, ID2
 !
       REAL              :: A, C_BJ, HM, HRMS, N1, N2
       REAL              :: QB_WC, SIG0, STP_OV, STP_PM
-      REAL              :: DDIF, DELTA,                                   40.31 40.30
-     &                     DSTEEP, EBIN, XFAC                             40.30
-      REAL              :: CPOW, CTOT, GAMMA                              40.41
-      REAL              :: BINSIZE                                        40.41
 !
       REAL, ALLOCATABLE :: C_K(:), C_LH(:), WCAP(:), WCIMPL(:)
-      REAL, ALLOCATABLE :: CUMSTP(:,:), DSIGMA(:)                         40.30
-      REAL, ALLOCATABLE :: FCOS(:)                                        40.31
 
       REAL              :: B, P, BRKD, FBR                                40.63 40.53
       REAL              :: PWCAP1, PWCAP2, PWCAP9, PWCAP10, PWCAP11       40.63
       REAL              :: EF(MSC)                                        40.53
-!
-      CHARACTER*20 NUMSTR, CHARS                                          40.51
-      CHARACTER*80 MSGSTR                                                 40.51
+      REAL              :: WCCUR(MDC,MSC)                                 41.11
 !
 !     8. REMARKS
 !
@@ -1357,7 +1314,6 @@
 !     If IWCAP = 1, 2, or 5; Calculate C_K
 !     If IWCAP = 4 or 5; Calcualte C_BJ
 !     If IWCAP = 3; Calculate C_LH
-!     If IWCAP = 6; Calculate cumulative steepness spectrum
 !     IF IWCAP = 7; Calculate with Alves and Banner (2003)
 !     ------------------------------------------------------------
 !     For frequency dependent part of the spectrum
@@ -1373,19 +1329,6 @@
 !
       IF (LTRACE) CALL STRACE (IENT,'WCAP')
 !
-      MXWCP = 7
-      IF (IWCAP.GT.MXWCP) THEN                                            40.30
-!
-! Error message
-!
-         CHARS = NUMSTR(MXWCP+1,RNAN,'(I1)')
-         CALL TXPBLA(CHARS,IF,IL)
-         MSGSTR = 'Value for IWCAP should be less than '//
-     &            CHARS(IF:IL)
-         CALL MSGERR ( 4, MSGSTR )
-         RETURN
-      END IF
-!
 ! Initialisation
 !
       IF (ETOT.LE.0.) RETURN
@@ -1395,9 +1338,8 @@
       IF (EDRKTOT.LE.0.) RETURN
 !
       ALLOCATE (C_K(MSC), C_LH(MSC), WCAP(MSC), WCIMPL(MSC))
-      ALLOCATE (CUMSTP(0:MSC,MDC), DSIGMA(1:MSC))                         40.30
-      ALLOCATE (FCOS(1:MDC))                                              40.31
       WCIMPL(1:MSC) = 0.
+      WCCUR(1:MDC,1:MSC) = 0.                                             41.11
 !
 ! Calculate coefficients
 !
@@ -1458,45 +1400,6 @@
         END DO
       END IF
 !
-! In case of cumulative steepness method, calculate its spectrum          40.30
-!
-      IF ( IWCAP.EQ.6 ) THEN                                              40.30
-         XFAC      = (SPCSIG(3)-SPCSIG(1))/(2.*SPCSIG(2))                 40.30
-         DSIGMA(:) = XFAC*SPCSIG(:)                                       40.30
-         DELTA     = SPCDIR(2,1)-SPCDIR(1,1)                              40.30
-!
-!        --- precompute cos dependence                                    40.31
-!                                                                         40.31
-         DO ID1 = 1, MDC                                                  40.31
-           DDIF = REAL(ID1-1)*DELTA                                       40.31
-           FCOS(ID1) = (ABS(COS(DDIF)))**PWCAP(13)                        40.31
-         END DO                                                           40.31
-!
-!        --- compute normalisation coefficient                            40.41
-!
-         CPOW = PWCAP(13)                                                 40.41
-         IF (CPOW > 10) THEN                                              40.41
-           CTOT = SQRT(CPOW/(2.*PI))*(1.+0.25/CPOW)                       40.41
-         ELSE                                                             40.41
-           CTOT = 1./SQRT(PI)*GAMMA(0.5*CPOW+1.)/GAMMA(0.5*CPOW+0.5)      40.41
-         END IF                                                           40.41
-!
-         CUMSTP = 0.                                                      40.30
-         DO IS = 1,ISSTOP                                                 40.30
-            BINSIZE = SPCSIG(IS)*DELTA*DSIGMA(IS)                         40.41
-            DO ID1 = 1,MDC                                                40.30
-               CUMSTP(IS,ID1) = CUMSTP(IS-1,ID1)                          40.31
-               DO ID2 = 1,MDC                                             40.30
-                  EBIN= AC2(ID2,IS,KCGRD(1))*BINSIZE                      40.41 40.30
-                  DSTEEP = KWAVE(IS,1)**2*EBIN*FCOS(ABS(ID1-ID2)+1)       40.31 40.30
-                  CUMSTP(IS,ID1) = CUMSTP(IS,ID1) + DSTEEP                40.31 40.30
-               END DO                                                     40.30
-            END DO                                                        40.30
-         END DO                                                           40.30
-         CUMSTP = PWCAP(12)*CUMSTP                                        40.30
-         CUMSTP = CUMSTP*CTOT                                             40.41
-      END IF                                                              40.30
-!
 ! Calculate dissipation according to Alves & Banner (2003)                40.53
 !
       IF ( IWCAP.EQ.7 ) THEN                                              40.53
@@ -1551,11 +1454,20 @@
      &    (GRAV**(0.5)*KWAVE(IS,1)**(0.5)/SPCSIG(IS))**(PWCAP(10)/2-1) *  40.53
      &    GRAV**(0.5)*KWAVE(IS,1)**(0.5)                                  40.53
      &    + (1.-FBR) * C_K(IS) * SIGM_10 * (KWAVE(IS,1) / KM_WAM)         40.63
+!                                                                         41.11
+!  Calculate current-induced enhanced dissipation WCCUR(ID,IS)            41.11
+!                                                                         41.11
+          IF ( PWCAP(14).NE.0. ) THEN                                     41.11
+             DO ID = 1, MDC                                               41.11
+                WCCUR(ID,IS)=PWCAP(14)*MAX((CAS(ID,IS,1)/SPCSIG(IS)),0.)  41.11
+     &                       *(B/BRKD)**(P/2.)                            41.11
+             ENDDO                                                        41.11
+          ENDIF                                                           41.11
                                                                           40.53
         END DO                                                            40.53
       END IF                                                              40.53
 !
-      IF ( IWCAP.LT.6 ) THEN                                              40.51 40.30
+      IF ( IWCAP.LT.7 ) THEN                                              40.51 40.30
 !
 ! Calculate the whitecapping source term WCAP(IS)
 !
@@ -1588,37 +1500,19 @@
 !
 ! Fill the diagonal of the matrix and the PLWCAP-array
 !
-      IF ( IWCAP.NE.6 ) THEN                                              40.30
-
-        DO IS=1, ISSTOP
+      DO IS=1, ISSTOP
 !
-!         Only fill the values for the current sweep
+!        Only fill the values for the current sweep
 !
-          DO IDDUM = IDCMIN(IS), IDCMAX(IS)
+         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
             ID = MOD(IDDUM - 1 + MDC, MDC) + 1
             IMATDA(ID,IS)   = IMATDA(ID,IS)   + WCAP(IS)
             DISSC1(ID,IS,1) = DISSC1(ID,IS,1) + WCAP(IS)                  40.67 40.12
             IF (TESTFL) PLWCAP(ID,IS,IPTST) = -1.*(WCAP(IS)-WCIMPL(IS))
-          END DO
-        END DO
-
-      ELSE
-
-        DO IS=1, ISSTOP                                                   40.30
-!                                                                         40.30
-!         Only fill the values for the current sweep                      40.30
-!                                                                         40.30
-          DO IDDUM = IDCMIN(IS), IDCMAX(IS)                               40.30
-            ID = MOD(IDDUM - 1 + MDC, MDC) + 1                            40.30
-            IMATDA(ID,IS)   = IMATDA(ID,IS)   + CUMSTP(IS,ID)             40.30
-            DISSC1(ID,IS,1) = DISSC1(ID,IS,1) + CUMSTP(IS,ID)             40.67 40.30
-            IF (TESTFL) PLWCAP(ID,IS,IPTST) = -1.*CUMSTP(IS,ID)           40.30
-          END DO                                                          40.30
-        END DO                                                            40.30
-
-      END IF
+         END DO
+      END DO
 !
-! Add the implicit part to the right-hand side
+! Add the implicit part to the right-hand side, if appropriate
 !
       IF ((IWCAP.EQ.4).OR.
      &    (IWCAP.EQ.5)) THEN
@@ -1636,8 +1530,26 @@
         END DO
       END IF
 !
+! Add extra dissipation on opposing current, if appropriate               41.11
+!
+      IF ( IWCAP.EQ.7 .AND. PWCAP(14).NE.0. ) THEN
+!
+         DO IS=1, ISSTOP
+!
+!           Only fill the values for the current sweep
+!
+            DO IDDUM = IDCMIN(IS), IDCMAX(IS)
+               ID = MOD(IDDUM - 1 + MDC, MDC) + 1
+               IMATDA(ID,IS)   = IMATDA(ID,IS)   + WCCUR(ID,IS)
+               DISSC1(ID,IS,1) = DISSC1(ID,IS,1) + WCCUR(ID,IS)
+               IF (TESTFL) PLWCAP(ID,IS,IPTST) = PLWCAP(ID,IS,IPTST) -
+     &                                           WCCUR(ID,IS)
+            END DO
+         END DO
+!
+      ENDIF
+!
       DEALLOCATE (C_K, C_LH, WCAP, WCIMPL)
-      DEALLOCATE (CUMSTP, DSIGMA, FCOS)                                   40.31 40.30
 !
       RETURN
       END SUBROUTINE SWCAP
--- swancom5.ftn	2009-12-17 10:30:34.000000000 +0100
+++ swancom5.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -2563,7 +2563,7 @@
       IF (LTRACE) CALL STRACE (IENT,'SORDUP')
 !
       IF(NSTATC.EQ.1)THEN                                                 33.10
-         CALL MSGERR (4, 'SORDUP scheme is for stationary mode only.')    33.10
+         CALL MSGERR (3, 'SORDUP scheme is for stationary mode only.')    33.10
       END IF                                                              33.10
       IF (TESTFL .AND. ITEST .GE. 120) THEN
         WRITE(PRINTF,*) ' Initial matrix coefficients at SORDUP : '
--- SwanCompUnstruc.ftn90	2009-12-17 10:30:36.000000000 +0100
+++ SwanCompUnstruc.ftn90	2009-12-17 10:31:32.000000000 +0100
@@ -465,7 +465,9 @@
        ! some initializations
        !
        if ( IQUAD  > 2 ) memnl4(1:MDC,1:MSC,ivlow:ivup) = 0.
-       if ( ITRIAD > 0 .or. ISURF == 5 ) compda(ivlow:ivup,JURSEL) = 0.
+       if ( ITRIAD > 0 &
+            .or. ISURF == 5 &
+                       ) compda(ivlow:ivup,JURSEL) = 0.
        !
        compda(ivlow:ivup,JDISS) = 0.
        compda(ivlow:ivup,JLEAK) = 0.
@@ -884,7 +886,9 @@
                                       da1m                , da2c                , da2p                , da2m                , &
                                       sfnl                , dsnl                , memnl4              , wwint               , &
                                       wwawg               , wwswg               , cgo                 , compda(1,JUSTAR)    , &
-                                      compda(1,JZEL)      , spcdir              , anywnd              , disc0               , &
+                                      compda(1,JZEL)      , spcdir              , anywnd              ,                       &
+                                      cas                 ,                                                                   &
+                                      disc0               ,                                                                   &
                                       disc1               , genc0               , genc1               , redc0               , &
                                       redc1               , xis                 , compda(1,JFRC2)     , it                  , &
                                       compda(1,JURSEL)    , anybin              , reflso              )
--- SwanConvAccur.ftn90	2009-12-17 10:30:36.000000000 +0100
+++ SwanConvAccur.ftn90	2009-12-17 10:31:26.000000000 +0100
@@ -257,7 +257,7 @@
     accur = real(npacc)*100./real(nwetp)
     !$omp end single
     !
- 11 format(11x,'dHrel          ','dHoval         ','dTm01rel       ','dTm01oval ')
- 12 format(1x,ss,'k=',i5,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
+ 11 format(13x,'dHrel          ','dHoval         ','dTm01rel       ','dTm01oval ')
+ 12 format(1x,ss,'k=',i7,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
     !
 end subroutine SwanConvAccur
--- SwanConvStopc.ftn90	2009-12-17 10:30:36.000000000 +0100
+++ SwanConvStopc.ftn90	2009-12-17 10:31:26.000000000 +0100
@@ -251,7 +251,7 @@
     accur = real(npacc)*100./real(nwetp)
     !$omp end single
     !
- 11 format(11x,'dHabs          ','dHrel          ','Curvature H    ','dTabs          ','dTrel          ','Curvature T    ')
- 12 format(1x,ss,'k=',i5,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
+ 11 format(13x,'dHabs          ','dHrel          ','Curvature H    ','dTabs          ','dTrel          ','Curvature T    ')
+ 12 format(1x,ss,'k=',i7,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
     !
 end subroutine SwanConvStopc
--- swan.edt	2009-12-17 10:30:36.000000000 +0100
+++ swan.edt	2009-12-17 10:31:22.000000000 +0100
@@ -90,19 +90,17 @@
 !         | WESTHuysen                  |
 !
 !
-!         | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]  |
-!         |                                                  |
-!         |   JANSsen  [cds1]  [delta] [pwtail]              |
-!         |                                                  |
-!         |   LHIG     [cflhig]                              |
-!         |                                                  |
-!   WCAP <    BJ       [bjstp] [bjalf]                        >
-!         |                                                  |
-!         |   KBJ      [bjstp] [bjalf] [kconv]               |
-!         |                                                  |
-!         |   CSM      [cst] [pow]                           |
-!         |                                                  |
-!         |   AB       [cds2] [br] [p0] [powst] [powk]       |
+!         | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]
+!         |
+!         |   JANSsen  [cds1]  [delta] [pwtail]
+!         |
+!         |   LHIG     [cflhig]
+!         |
+!   WCAP <    BJ       [bjstp] [bjalf]
+!         |
+!         |   KBJ      [bjstp] [bjalf] [kconv]
+!         |
+!         |   AB       [cds2] [br] [p0] [powst] [powk]
 !
 !   QUADrupl [iquad] [limiter] [lambda] [cnl4] [csh1] [csh2] [csh3]
 !
@@ -113,12 +111,10 @@
 !        | -> CON [alpha] [gamma]                                      |
 !        |                                                             |
 !        |    VAR [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
-!        |                                                             |
-!   BRE <     RUE [alpha] [a] [b]                                       >   &
+!   BRE <                                                               > &
+!        |    RUE [alpha] [a] [b]                                      |
 !        |                                                             |
 !        |    TG  [alpha] [gamma] [pown]                               |
-!        |                                                             |
-!        !    BIP [alpha] [pown] [bref]                                |
 !
 !     ( FREQDep [power] [fmin] [fmax] )
 !
--- swanimp.tex	2009-12-17 10:30:36.000000000 +0100
+++ swanimp.tex	2009-12-17 10:31:23.000000000 +0100
@@ -15,7 +15,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72ABCD}
+{\Large\bf SWAN Cycle III version 40.72ABCDE}
 \end{center}
 
 \cleardoublepage
@@ -330,12 +330,13 @@
 Compaq True 64 Alpha (DEC ALFA)     & OSF1        & Compaq \\
 Sun SPARC                           & Solaris     & Sun \\
 PA-RISC (HP 9000 series 700/800)    & HP-UX v11   & HP \\
+IBM Power6 (pSeries 575)            & Linux       & IBM \\
 Intel Pentium (32-bit) PC           & Linux       & GNU (g95) \\
 Intel Pentium (32-bit) PC           & Linux       & GNU (gfortran) \\
 Intel Pentium (32-bit) PC           & Linux       & Intel \\
 Intel Pentium (64-bit) PC           & Linux       & Intel \\
 Intel Itanium (64-bit) PC           & Linux       & Intel \\
-Intel Pentium (32-bit) PC           & Linux       & Portland Group \\
+Intel Pentium (64-bit) PC           & Linux       & Portland Group \\
 Intel Pentium (32-bit) PC           & Linux       & Lahey \\
 Intel Pentium (32-bit) PC           & MS Windows  & Intel \\
 Intel Pentium (64-bit) PC           & MS Windows  & Intel \\
--- SwanInterpolateOutput.ftn90	2009-12-17 10:30:36.000000000 +0100
+++ SwanInterpolateOutput.ftn90	2009-12-17 10:31:28.000000000 +0100
@@ -58,6 +58,7 @@
     use swcomm2
     use swcomm3
     use m_obsta
+    use outp_data
     use SwanGriddata
     use SwanGridobjects
 !
@@ -127,6 +128,11 @@
 !
     if (ltrace) call strace (ient,'SwanInterpolateOutput')
     !
+    if ( LCOMPGRD .and. mip == nverts ) then
+       foutp = finp
+       return
+    endif
+    !
     ! point to vertex and cell objects
     !
     vert => gridobject%vert_grid
--- swanmain.ftn	2009-12-17 10:30:36.000000000 +0100
+++ swanmain.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -793,7 +793,6 @@
 !     40.23, Nov. 02: parameter PROPFL added
 !     40.23, Dec. 02: reset of some default variables
 !     40.30, Mar. 03: introduction distributed-memory approach using MPI
-!     40.30, Mar. 03: coefficients for cumulative steepness method added
 !     40.08, Mar. 03: Unneccessary variable deleted
 !     40.31, Nov. 03: removing POOL construction and HPGL functionality
 !     40.41, Sep. 04: output variables TMM10 and RTMM10 added
@@ -859,6 +858,7 @@
       CALL BUGFIX ('B')
       CALL BUGFIX ('C')
       CALL BUGFIX ('D')
+      CALL BUGFIX ('E')
 !
       CALL OCPINI ('swaninit', .TRUE.,INERR)                              34.01
       IF (INERR.GT.0) RETURN                                              34.01
@@ -1197,12 +1197,13 @@
       PWCAP(7) = 1.
       PWCAP(8) = 0.75
 !
-!     *** (12): Proportionality coefficient in cumulative steepness       40.30
-!               method                                                    40.30
-!         (13): Power of cosine term in cumulative steepness method       40.30
-!
-      PWCAP(12) = 0.5                                                     40.41 40.31 40.30
-      PWCAP(13) = 2.0                                                     40.31 40.30
+!     *** (13): Flag for activation of nonlinear dispersion (not used)    41.11
+!         (14): Coefficient for enhanced current-induced dissipation      41.11
+!         (15): Exponent for doppler shifting factor (must be =1.)        41.11
+!
+      PWCAP(13) = 0.                                                      41.11
+      PWCAP(14) = 0.                                                      41.11
+      PWCAP(15) = 1.0                                                     41.11
 !
 !     PBOT(1)   = 0.005          modified
       PBOT(1)   = 0.0                                                     20.68
@@ -3223,11 +3224,19 @@
 !
 !        ***** compute depth and water level *****
 !
-         DEP = SVALQI (XP, YP, 1, DEPTH, 1, 0, 0)                         40.80
+         IF ( IGTYPE(1).EQ.3 ) THEN
+            DEP = DEPTH(INDX)
+         ELSE
+            DEP = SVALQI (XP, YP, 1, DEPTH, 1, 0, 0)                      40.80
+         ENDIF
          COMPDA(INDX,JBOTLV) = DEP                                        40.80
 !
          IF (VARWLV) THEN                                                 40.80
-            WLVL = SVALQI (XP, YP, 7, WLEVL, 1, 0, 0)                     40.80
+            IF ( IGTYPE(7).EQ.3 ) THEN
+               WLVL = WLEVL(INDX)
+            ELSE
+               WLVL = SVALQI (XP, YP, 7, WLEVL, 1, 0, 0)                  40.80
+            ENDIF
             COMPDA(INDX,JWLV2) = WLVL
             IF (JWLV1.GT.1) COMPDA(INDX,JWLV1) = WLVL                     40.80
             IF (JWLV3.GT.1) COMPDA(INDX,JWLV3) = WLVL                     40.80
@@ -3245,8 +3254,16 @@
 !
          IF (ICUR.EQ.1 .AND. IGTYPE(2) .GE. 1) THEN                       40.80
             IF (DEPW.GT.0.) THEN
-              UU  = SVALQI (XP, YP, 2, UXB, 0, 0, 0)                      40.80
-              VV  = SVALQI (XP, YP, 3, UYB, 0, 0, 0)                      40.80
+              IF ( IGTYPE(2).EQ.3 ) THEN
+                 UU = UXB(INDX)
+              ELSE
+                 UU = SVALQI (XP, YP, 2, UXB, 0, 0, 0)                    40.80
+              ENDIF
+              IF ( IGTYPE(3).EQ.3 ) THEN
+                 VV = UYB(INDX)
+              ELSE
+                 VV = SVALQI (XP, YP, 3, UYB, 0, 0, 0)                    40.80
+              ENDIF
               VTOT = SQRT (UU*UU + VV*VV)
               CGMAX = PNUMS(18)*SQRT(GRAV*DEPW)
               IF (VTOT .GT. CGMAX) THEN
@@ -3273,7 +3290,11 @@
 !         ***** compute variable friction coefficient *****
 !
          IF (VARFR) THEN
-            FRI = SVALQI (XP, YP, 4, FRIC, 1, 0, 0)                       40.80
+            IF ( IGTYPE(4).EQ.3 ) THEN
+               FRI = FRIC(INDX)
+            ELSE
+               FRI = SVALQI (XP, YP, 4, FRIC, 1, 0, 0)                    40.80
+            ENDIF
             COMPDA(INDX,JFRC2) = FRI                                      40.80
             COMPDA(INDX,JFRC3) = FRI                                      40.80
          ENDIF
@@ -3281,8 +3302,16 @@
 !        ***** compute variable wind velocity *****
 !
          IF (VARWI) THEN
-            UU  = SVALQI (XP, YP, 5, WXI, 0, 0, 0)                        40.80
-            VV  = SVALQI (XP, YP, 6, WYI, 0, 0, 0)                        40.80
+            IF ( IGTYPE(5).EQ.3 ) THEN
+               UU = WXI(INDX)
+            ELSE
+               UU = SVALQI (XP, YP, 5, WXI, 0, 0, 0)                      40.80
+            ENDIF
+            IF ( IGTYPE(6).EQ.3 ) THEN
+               VV = WYI(INDX)
+            ELSE
+               VV = SVALQI (XP, YP, 6, WYI, 0, 0, 0)                      40.80
+            ENDIF
             COMPDA(INDX,JWX2) =  UU*COSWC + VV*SINWC
             COMPDA(INDX,JWY2) = -UU*SINWC + VV*COSWC
             IF (JWX3.GT.1) COMPDA(INDX,JWX3) = COMPDA(INDX,JWX2)          40.80
@@ -3292,7 +3321,11 @@
 !     ***** compute variable air-sea temperature difference *****         40.80
 !
          IF (VARAST) THEN
-            ASTD = SVALQI (XP, YP, 10, ASTDF, 1, 0, 0)                    40.80
+            IF ( IGTYPE(10).EQ.3 ) THEN
+               ASTD = ASTDF(INDX)
+            ELSE
+               ASTD = SVALQI (XP, YP, 10, ASTDF, 1, 0, 0)                 40.80
+            ENDIF
             COMPDA(INDX,JASTD2) = ASTD                                    40.80
             COMPDA(INDX,JASTD3) = ASTD                                    40.80
          ENDIF
@@ -4318,8 +4351,7 @@
 !     40.08, Mar. 03: "GE" changed to "EQ" and warning message related to
 !                     curvilinear coordinates added
 !     40.31, Dec. 03: removing POOL mechanism
-!     40.41, Jul. 04: added error concerning CSM and DIA,
-!                     added warning concerning frequency-resolution and DIA,
+!     40.41, Jul. 04: added warning concerning frequency-resolution and DIA,
 !                     added check of resonance condition for triads
 !     40.41, Oct. 04: common blocks replaced by modules, include files removed
 !     40.53, Mar. 05: change Alves and Banner parameters in case of XNL
@@ -4344,17 +4376,13 @@
 !
 !  6. Local variables
 !
-!     CHARS :     array to pass character info to MSGERR
 !     MSGSTR:     string to pass message to call MSGERR
 !
-      CHARACTER*20 NUMSTR, CHARS(3)
       CHARACTER*80 MSGSTR
 !
 !  8. Subroutines used
 !
 !     MSGERR : Handles error messages according to severity
-!     NUMSTR : Converts integer/real to string
-!     TXPBLA : Removes leading and trailing blanks in string
 !
 !  9. Subroutines calling
 !
@@ -4376,7 +4404,6 @@
 !
 ! 13. Source text
 !
-      LOGICAL  STPNOW                                                     40.00
       LOGICAL  EQREAL                                                     40.53
       SAVE IENT
       DATA IENT/0/
@@ -4415,6 +4442,22 @@
 !
 !     -----------------------------------------------------------------
 !
+!     *** Check formulation for whitecapping ***
+!
+      IF ( IWCAP.GT.7 ) THEN
+         WRITE (MSGSTR, '(A,I2,A)')
+     &               'Unknown method for whitecapping (IWCAP=',IWCAP,')'
+         CALL MSGERR( 1, TRIM(MSGSTR) )
+         CALL MSGERR( 1,
+     &     'Whitecapping according to Komen et al. (1984) will be used')
+         IWCAP     = 1
+         PWCAP(1)  = 2.36E-5
+         PWCAP(2)  = 3.02E-3
+         PWCAP(9)  = 2.
+         PWCAP(10) = 0.
+         PWCAP(11) = 1.
+      ENDIF
+!
 !     *** WAM cycle 3 physics ***
 !
       IF  ( (IWIND .EQ. 3 .OR. IWIND .EQ. 5) .AND. IWCAP .NE. 1
@@ -4526,7 +4569,7 @@
        IF (MSC .EQ. 3 ) THEN
          CALL MSGERR(4,'Do not activate quadruplets for boundary ')
          CALL MSGERR(4,'option BIN -> use other option           ')
-         WRITE(PRINTF,*)
+         RETURN
        END IF
 !
       END IF                                                              30.72
@@ -4559,16 +4602,9 @@
          CALL MSGERR(4,'When triads are active the number of     ')       30.72
          CALL MSGERR(4,'directions chosen in the CGRID command   ')       30.72
          CALL MSGERR(4,'must be less than 200                    ')       30.72
+         RETURN
       END IF                                                              30.72
 !
-!     When CSM is applied only DIA should be used                         40.41
-!
-      IF (IWCAP.EQ.6 .AND. IQUAD.GT.3 .AND. IQUAD.NE.8) THEN              40.41 40.30
-         CALL MSGERR(3,                                                   40.30
-     &       'In case of cumulative steepness method, only DIA '//        40.30
-     &       'technique is supported')                                    40.30
-      END IF                                                              40.30
-!
 !     When Alves and Banner and XNL are applied change the parameters     40.53
 !
       IF ( IWCAP.EQ.7 .AND. (IQUAD.EQ.51 .OR. IQUAD.EQ.52 .OR.            40.53
@@ -5005,7 +5041,11 @@
         DO INDX = 1, nverts                                               40.80
            XP = xcugrd(INDX)
            YP = ycugrd(INDX)
-           DEP = SVALQI (XP, YP, 1, DEPTH, 1, 0, 0)
+           IF ( IGTYPE(1).EQ.3 ) THEN
+              DEP = DEPTH(INDX)
+           ELSE
+              DEP = SVALQI (XP, YP, 1, DEPTH, 1, 0, 0)
+           ENDIF
            COMPDA(INDX,JDP1) = COMPDA(INDX,JDP2)
            WLVL = COMPDA(INDX,JWLV2)
            DEPW = DEP + WLVL + WLEV
@@ -6418,11 +6458,19 @@
       DO INDX = 1, nverts                                                 40.80
          XP = xcugrd(INDX)
          YP = ycugrd(INDX)
-         UU = SVALQI (XP, YP, IGR1, ARR, 0, 0, 0)
+         IF ( IGTYPE(IGR1).EQ.3 ) THEN
+            UU = ARR(INDX)
+         ELSE
+            UU = SVALQI (XP, YP, IGR1, ARR, 0, 0, 0)
+         ENDIF
          IF (IGR2.EQ.0) THEN
             COMPDA(INDX,JX3) = UU
          ELSE
-            VV = SVALQI (XP, YP, IGR2, ARR2, 0, 0, 0)
+            IF ( IGTYPE(IGR2).EQ.3 ) THEN
+               VV = ARR2(INDX)
+            ELSE
+               VV = SVALQI (XP, YP, IGR2, ARR2, 0, 0, 0)
+            ENDIF
             COMPDA(INDX,JX3) =  UU*COSFC + VV*SINFC
             COMPDA(INDX,JY3) = -UU*SINFC + VV*COSFC
          ENDIF
--- swanout1.ftn	2009-12-17 10:30:34.000000000 +0100
+++ swanout1.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -769,6 +769,7 @@
       USE SWCOMM1                                                         40.41
       USE SWCOMM3                                                         40.41
       USE SWCOMM4                                                         40.41
+      USE OUTP_DATA
 !
 !
 !   --|-----------------------------------------------------------|--
@@ -1000,6 +1001,12 @@
         WRITE (PRTEST,'(A)') ' error SWODDC: no PSTYPE defined'           40.31
       ENDIF
 !
+      IF (PSNAME.EQ.'COMPGRID') THEN
+         LCOMPGRD=.TRUE.
+      ELSE
+         LCOMPGRD=.FALSE.
+      ENDIF
+!
       IF (ITEST.GE.100 .OR. IOUTES .GE. 30) THEN
         WRITE (PRTEST, 91) PSNAME, PSTYPE, MIP
   91    FORMAT (' Exit SWODDC  ', A16, 2X, A1, 3I6)
@@ -1351,6 +1358,7 @@
       USE TIMECOMM                                                        40.41
       USE M_PARALL                                                        40.31
       USE M_DIFFR                                                         40.21
+      USE OUTP_DATA
       USE SwanGriddata                                                    40.80
       USE SwanGridobjects                                                 40.91
 !
@@ -1501,9 +1509,13 @@
 !        ---find closest vertex for given point in case of unstructured grid
 !
          ALLOCATE(KVERT(MIP))
-         DO IP = 1, MIP
-            CALL SwanFindPoint ( VOQ(IP,1), VOQ(IP,2), KVERT(IP) )        41.07
-         ENDDO
+         IF (.NOT.LCOMPGRD) THEN
+            DO IP = 1, MIP
+               CALL SwanFindPoint ( VOQ(IP,1), VOQ(IP,2), KVERT(IP) )     41.07
+            ENDDO
+         ELSE
+            KVERT = 0
+         ENDIF
 !
       ENDIF
 !
@@ -2680,6 +2692,7 @@
       USE SWCOMM2                                                         40.80
       USE SWCOMM3                                                         40.41
       USE SWCOMM4                                                         40.41
+      USE OUTP_DATA
 !
 !
 !   --|-----------------------------------------------------------|--
@@ -2867,9 +2880,13 @@
            CALL SWOINA (XC(IP), YC(IP), AC2, ACLOC, KGRPNT, DEPXY,        40.86 30.50
      &                  CROSS(1,IP), EXCPT)                               40.86
         ELSE                                                              40.80
-           IF (.NOT.EQREAL(VOQ(IP,1),OVEXCV(1))) XP = VOQ(IP,1) - XOFFS   40.80
-           IF (.NOT.EQREAL(VOQ(IP,2),OVEXCV(2))) YP = VOQ(IP,2) - YOFFS   40.80
-           CALL SwanInterpolateAc ( ACLOC, XP, YP, AC2, EXCPT )           40.80
+           IF (.NOT.LCOMPGRD) THEN
+              IF (.NOT.EQREAL(VOQ(IP,1),OVEXCV(1))) XP=VOQ(IP,1)-XOFFS    40.80
+              IF (.NOT.EQREAL(VOQ(IP,2),OVEXCV(2))) YP=VOQ(IP,2)-YOFFS    40.80
+              CALL SwanInterpolateAc ( ACLOC, XP, YP, AC2, EXCPT )        40.80
+           ELSE
+              ACLOC(:,:) = AC2(:,:,IP)
+           ENDIF
         ENDIF                                                             40.80
         IF (EXCPT) GOTO 700                                               40.86
 !
--- swanparll.ftn	2009-12-17 10:30:36.000000000 +0100
+++ swanparll.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -3651,7 +3651,7 @@
 
 !        --- if it is not open or already open thru another request
 
-         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN
+         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN+NREOQ ) THEN
 
 !           --- open generic output file or reset reference number
 
@@ -3972,7 +3972,7 @@
 
 !        --- if it is not open or already open thru another request
 
-         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN
+         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN+NREOQ ) THEN
 
 !           --- open generic output file or reset reference number
 
@@ -4330,7 +4330,7 @@
 
 !        --- if it is not open or already open thru another request
 
-         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN ) THEN
+         IF ( .NOT.OPENED .OR. NREF.LE.HIOPEN+NREOQ ) THEN
 
 !           --- open generic output file or reset reference number
 
@@ -4344,7 +4344,7 @@
 
             MATLAB = INDEX( FILENM, '.MAT' ).NE.0 .OR.                    40.41
      &               INDEX (FILENM, '.mat' ).NE.0                         40.41
-            IF (MATLAB) THEN
+            IF (MATLAB .AND. .NOT.OPENED) THEN
                CLOSE(NREF)
                OPEN(UNIT=NREF, FILE=FILENM, FORM='UNFORMATTED',
      &              ACCESS='DIRECT', RECL=4)                              41.08
--- swanpre1.ftn	2009-12-17 10:30:34.000000000 +0100
+++ swanpre1.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -213,7 +213,6 @@
 !     40.23, Sep. 02: coefficient urslim must be stored in PTRIAD(5)
 !     40.23, Nov. 02: keyword FLUXLIM added
 !     40.30, Feb. 03: introduction distributed-memory approach using MPI
-!     40.30, Mar. 03: extension to cumulative steepness method
 !     40.08, Mar. 03: warning message added for use of curvilinear
 !                     coordinate system
 !     40.31, Dec. 03: removing POOL-mechanism
@@ -1963,12 +1962,12 @@
 !           | -> CONstant [alpha] [gamma]                                      |
 !           |                                                                  |
 !           |    VARiable [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
-!           |                                                                  |
-! BREaking <     RUEssink [alpha] [a] [b]                                       > &
+! BREaking <                                                                    > &
+!           |    RUEssink [alpha] [a] [b]                                      |
 !           |                                                                  |
 !           |    TG       [alpha] [gamma] [pown]                               |
 !           |                                                                  |
-!           !    BIPhase  [alpha] [pown] [bref]                                |
+!           !    WESTHuys [alpha] [pown] [bref] [shfac]                        |
 !
 !
 !       ( FREQDep [power] [fmin] [fmax] ) (fmin and fmax not documented)
@@ -1999,11 +1998,12 @@
           CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.0)
           CALL INREAL ('GAMMA',  PSURF(2), 'STA', 0.42)
           CALL INREAL ('POWN' ,  PSURF(4), 'STA', 4.0)
-        ELSE IF (KEYWIS('BIP')) THEN                                      41.03
+        ELSE IF (KEYWIS('WESTH')) THEN                                    41.09 41.03
           ISURF = 5
-          CALL INREAL ('ALPHA',  PSURF(1), 'STA', 1.0)
-          CALL INREAL ('POWN' ,  PSURF(4), 'STA', 2.4)
-          CALL INREAL ('BREF' ,  PSURF(5), 'STA', -1.5708)
+          CALL INREAL ('ALPHA',  PSURF(1), 'STA', 0.98)                   41.09
+          CALL INREAL ('POWN' ,  PSURF(4), 'STA', 2.5)                    41.09
+          CALL INREAL ('BREF' ,  PSURF(5), 'STA', -1.3963)                41.09
+          CALL INREAL ('SHFAC',  PSURF(6), 'STA', 500.)                   41.09
         ENDIF
 !
         CALL INKEYW ('STA', '  ')                                         41.06
@@ -2018,23 +2018,25 @@
 !
 !     ------------------------------------------------------------------
 !
-!     WCAP        parameters whitecapping
+!     WCAP        parameters whitecapping (NOT documented)
 !
 ! =============================================================
 !
-!        | ->KOMen    [cds2] [stpm] [powst] [delta] [powk]  |             34.00
-!        |                                                  |
-!        |   JANSsen  [cds1]  [delta] [pwtail]              |
-! WCAP  <                                                    >    (NOT documented)
-!        |   LHIG     [cflhig]                              |
-!        |                                                  |
-!        |   BJ       [bjstp] [bjalf]                       |
-!        |                                                  |
-!        |   KBJ      [bjstp] [bjalf] [kconv]               |
-!        |                                                  |             40.30
-!        |   CSM      [cst] [pow]                           |             40.30
-!        |                                                  |
-!        |   AB       [cds2] [br] [p0] [powst] [powk]       |             40.53
+!        | ->KOMen    [cds2] [stpm] [powst] [delta] [powk]                34.00
+!        |
+!        |   JANSsen  [cds1]  [delta] [pwtail]
+!        |
+!        |   LHIG     [cflhig]
+!        |
+! WCAP  <    BJ       [bjstp] [bjalf]
+!        |
+!        |   KBJ      [bjstp] [bjalf] [kconv]
+!        |                                                                40.30
+!        |   AB       [cds2] [br] [p0] [powst] [powk]                     40.53
+!        |
+!        |   WESTHuys [cds2] [br] [p0] [powst] [powk] &                   41.11
+!        |
+!        |            [nldisp] [cds3] [powfsh]                            41.11
 !
 ! =============================================================
 !
@@ -2083,11 +2085,6 @@
           CALL INREAL ('BJSTP' , PWCAP(6), 'UNC', 0.)
           CALL INREAL ('BJALF' , PWCAP(7), 'UNC', 0.)
           CALL INREAL ('KCONV' , PWCAP(8), 'UNC', 0.)
-        ELSE IF ( KEYWIS ('CSM')) THEN                                    40.30
-!         *** whitecapping according to cumulative steepness method  ***  40.30
-          IWCAP = 6                                                       40.30
-          CALL INREAL ('CST', PWCAP(12), 'UNC', 0.)                       40.30
-          CALL INREAL ('POW', PWCAP(13), 'UNC', 0.)                       40.30
         ELSE IF ( KEYWIS ('AB')) THEN                                     40.53
 !         *** whitecapping according to Alves and Banner (2003) ***       40.53
           IWCAP = 7                                                       40.53
@@ -2096,6 +2093,18 @@
           CALL INREAL ('P0',    PWCAP(10), 'STA', 4.)                     40.53
           CALL INREAL ('POWST', PWCAP(9),  'STA', 0.)                     40.53
           CALL INREAL ('POWK',  PWCAP(11), 'STA', 0.)                     40.53
+        ELSE IF ( KEYWIS ('WESTH')) THEN                                  41.11
+!         *** whitecapping according to Vd Westhuysen et al. (2007) ***   41.11
+!         *** incl. enhanced dissipation on current (yet unpublished)     41.11
+          IWCAP = 7                                                       41.11
+          CALL INREAL ('CDS2',  PWCAP(1),  'STA', 5.0E-5)                 41.11
+          CALL INREAL ('BR',    PWCAP(12), 'STA', 1.75E-3)                41.11
+          CALL INREAL ('P0',    PWCAP(10), 'STA', 4.)                     41.11
+          CALL INREAL ('POWST', PWCAP(9),  'STA', 0.)                     41.11
+          CALL INREAL ('POWK',  PWCAP(11), 'STA', 0.)                     41.11
+          CALL INREAL ('NLDISP',PWCAP(13), 'STA', 0.)                     41.11
+          CALL INREAL ('CDS3'  ,PWCAP(14), 'STA', 1.5)                    41.11
+          CALL INREAL ('POWFSH',PWCAP(15), 'STA', 1.)                     41.11
         ELSE
           CALL WRNKEY
         ENDIF
--- SwanPrintGridInfo.ftn90	2009-12-17 10:30:32.000000000 +0100
+++ SwanPrintGridInfo.ftn90	2009-12-17 10:31:29.000000000 +0100
@@ -166,13 +166,13 @@
     ! format statements
     !
  100 format(// ' The unstructured grid contains ',a,' generated by ',a//      &
-               ' Number of vertices          = ',i6//                         &
-               ' Number of cells             = ', i6/                         &
-               '    Number of internal cells = ', i6/                         &
-               '    Number of boundary cells = ',i6//                         &
-               ' Number of faces             = ',i6/                          &
-               '    Number of internal faces = ',i6/                          &
-               '    Number of boundary faces = ',i6//)
+               ' Number of vertices          = ',i7//                         &
+               ' Number of cells             = ', i7/                         &
+               '    Number of internal cells = ', i7/                         &
+               '    Number of boundary cells = ',i7//                         &
+               ' Number of faces             = ',i7/                          &
+               '    Number of internal faces = ',i7/                          &
+               '    Number of boundary faces = ',i7//)
  200 format(' The minimum gridsize =',f12.5)
  300 format(' The maximum gridsize =',f12.5)
 
--- SwanReadADCGrid.ftn90	2009-12-17 10:30:36.000000000 +0100
+++ SwanReadADCGrid.ftn90	2009-12-17 10:31:30.000000000 +0100
@@ -60,7 +60,7 @@
 !
 !   Local variables
 !
-    character(80)           :: filenm   ! file name
+    character(7)            :: grdfil   ! file name (fort.14)
     integer, save           :: ient = 0 ! number of entries in this subroutine
     integer                 :: idum     ! dummy integer
     integer                 :: ii       ! auxiliary integer
@@ -89,10 +89,10 @@
     !
     ! open file fort.14
     !
-    filenm = 'fort.14'
     ndsd   = 0
     iostat = 0
-    call for (ndsd, trim(filenm), 'OF', iostat)
+    grdfil = 'fort.14'
+    call for (ndsd, grdfil, 'OF', iostat)
     if (stpnow()) goto 900
     !
     ! skip first line
--- swanser.ftn	2009-12-17 10:30:36.000000000 +0100
+++ swanser.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -3621,8 +3621,6 @@
 !                                                                         40.18
 !  8. Subroutines used                                                    40.09
 !                                                                         40.09
-      LOGICAL :: STPNOW                                                   40.09
-!
 !  9. Subroutines calling                                                 40.09
 !                                                                         40.09
 !     SWTRCF                                                              40.09
@@ -3808,7 +3806,6 @@
       END DO                                                              40.09
 !                                                                         40.18
       IF (LREFDIFF .GT. 0) DEALLOCATE (PRDIF)
-      IF (STPNOW()) RETURN                                                40.09
       RETURN                                                              40.09
 !     End of subroutine REFLECT                                           40.09
       END                                                                 40.09
--- swantech.tex	2009-12-17 10:30:34.000000000 +0100
+++ swantech.tex	2009-12-17 10:37:33.000000000 +0100
@@ -17,7 +17,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72ABC}
+{\Large\bf SWAN Cycle III version 40.72ABCDE}
 \end{center}
 
 \cleardoublepage
@@ -260,7 +260,7 @@
 This is called the random-phase model.
 \\[2ex]
 \noindent
-In the presene of the ambient current, it is assumed that it is uniform with respect to the vertical
+In the presence of the ambient current, it is assumed that it is uniform with respect to the vertical
 co-ordinate and the changes in the mean flow within a wave length are so small that they affect only
 negligibly the dispersion relation. The absolute radian frequency $\omega$ then equals the sum of
 the relative radian frequency $\sigma$ and the multiplication of the wave number and ambient current
@@ -505,7 +505,8 @@
 
 First, in Section~\ref{sec:gencon} general concepts of the physical processes of generation, dissipation and
 nonlinear wave-wave interactions that are implemented in SWAN are outlined. Next, complete expressions for
-these physical processes are given in subsequent sections.
+these physical processes are given in subsequent sections. Finally, for completeness, the first- and second-generation
+formulations as employed in SWAN are outlined in Section~\ref{sec:relaxmod}.
 
 \subsection{General concepts} \label{sec:gencon}
 
@@ -597,27 +598,20 @@
 from alternative calibrations of the Komen~{\it et~al} (1984) expression, e.g. Rogers~{\it et~al} (2003), to
 alternative ways of calculating mean spectral steepness, e.g. Van Vledder and Hurdle (2002).
 In SWAN, another alternative is presented.
-%In SWAN, two alternatives are presented.
 \nocite{Vle02H}
 \\[2ex]
 \noindent
-%An alternative formulation for whitecapping is based on the Cumulative Steepness Method as described
-%in Hurdle and Van Vledder (2004). With this method dissipation due to whitecapping depends on the steepness of
-%the wave spectrum at and below a particular frequency.
-%\\[2ex]
-%\noindent
-%A second
-An
-alternative of the whitecapping expression is based on Alves and Banner (2003). This expression is based
+This alternative is proposed by Van der Westhuysen~{\it et~al} (2007) and Van der Westhuysen (2007),
+based on the whitecapping expression of Alves and Banner (2003). This expression is based
 on experimental findings that whitecapping dissipation appears to be related to the nonlinear hydrodynamics within
 wave groups. This yields a dissipation term that primarily depends on quantities that are local in the frequency
 spectrum, as opposed to ones that are distributed over the spectrum, as in the expression of Komen~{\it et~al} (1984).
 However, the final whitecapping expression proposed by Alves and Banner (2003) features additional dependencies on
 the spectral mean wavenumber and steepness, which is problematic in situations of mixed sea and swell often
-encountered in the nearshore. Therefore, their whitecapping expression is applied here without these mean spectral
-dependencies. This adapted whitecapping expression is used together with a wind input term that is based on that
-of Yan (1987). Further information and details can be found in Van der Westhuysen~{\it et~al} (2007).
-\nocite{Alv03B,Yan87,Wes07ZB}
+encountered in the nearshore. Therefore, their whitecapping expression is applied in Van der Westhuysen (2007)
+without these mean spectral dependencies. This adapted whitecapping expression is used together with a wind input
+term that is based on that of Yan (1987).
+\nocite{Alv03B,Yan87,Wes07ZB,Wes07}
 \\[2ex]
 \noindent
 In shallow water the orbital motions of the water particles, induced by surface waves, extend
@@ -670,6 +664,14 @@
 \end{equation}
 in which $E_{\rm tot}$ is the total wave energy and $D_{\rm tot}~<~0$ is the rate of dissipation of
 the total energy due to wave breaking according to Battjes and Janssen (1978).
+%Adding a quadratic dependency
+%on frequency as suggested by Mase and Kirby (1992) and supported by Elgar~{\it et~al}. (1997)
+%and Chen~{\it et~al}. (1997)
+%may have effect on the SWAN results. This option is available in SWAN.
+%%Chen~{\it et~al}. (1997) inferred from observations and simulations with a
+%%Boussinesq model that the high-frequency levels are insensitive to such frequency dependency because
+%%an increased dissipation at high frequencies is compensated approximately by increased nonlinear
+%%energy transfer (but they did find the frequency dependency to be relevant in time domain).
 The value of $D_{\rm tot}$ depends critically on the breaker parameter $\gamma = H_{\max}/d$
 (in which $H_{\max}$ is the maximum possible individual wave height in the local water depth $d$).
 In SWAN, both a constant value and a variable value are available. Examples of a variable breaker parameter
@@ -890,40 +892,15 @@
 Cycle~4; Komen~{\it et~al}., 1994).
 \\[2ex]
 \noindent
-%\underline{Whitecapping: CSM formulation}\\[2ex]
-%An alternative formulation for whitecapping is based on the Cumulative Steepness Method as described
-%in Hurdle and Van Vledder (2004). With this method dissipation due to whitecapping depends on the steepness of
-%the wave spectrum at and below a particular frequency. It is defined as (directionally dependent):
-%\begin{equation}
-%  S_{\rm st} (\sigma,\theta) = A_m \int_{0}^{\sigma} \int_{0}^{2\pi} k^2 | \cos (\theta - \theta')|^m E
-%  (\sigma,\theta)d\sigma d\theta
-%  \label{eq2-5}
-%\end{equation}
-%with $A_m$ the normalisation coefficient as determined by
-%\begin{equation}
-%  \int_{0}^{2\pi} A_m \cos^m (\theta) d\theta = 1
-%\end{equation}
-%In expression (\ref{eq2-5}) the coefficient $m$ controls the directional dependence. It is expected that this
-%coefficient will be order 1 if the straining mechanism is dominant, $m$ is more than 10 if other mechanism play
-%a role (e.g. instability that occurs when vertical acceleration in the waves becomes greater than gravity).
-%Default in SWAN is $m=2$. The alternative whitecapping source term is given by
-%\begin{equation}
-%  S_{\rm wc}^{\rm st} = -C_{\rm wc}^{\rm st} \left( S_{\rm st} (\sigma,\theta) \right)^p E (\sigma,\theta)
-%  \label{eq2-6}
-%\end{equation}
-%with $C_{wc}^{st}$ a tuneable coefficient and $p$ a parameter controlling the proportionality of the dissipation
-%rate on the steepness. In SWAN, $p=1$ is assumed.
-%\nocite{Hur04V}
-%\\[2ex]
-%\noindent
 \underline{Whitecapping: saturation-based model}\\[2ex]
-The whitecapping formulation used in SWAN is an adapted form of the expression of Alves and Banner (2003), which is based
+An alternative description for whitecapping in SWAN is given by Van der Westhuysen~{\it et~al} (2007) and
+Van der Westhuysen (2007), which is an adapted form of the expression of Alves and Banner (2003). The latter is based
 on the apparent relationship between wave groups and whitecapping dissipation. This adaption is due to the fact that it
 can also be applied to mixed sea-swell conditions and in shallow water. This was done by removing the dependencies on mean
 spectral steepness and wavenumber in the original expression, and by applying source term scaling arguments for its
 calibration (see below). This led to the following expression for whitecapping dissipation
 \begin{equation}
-  S_{\rm ds,w} (\sigma,\theta) = -C_{\rm ds} \left( \frac{B(k)}{B_{\rm r}} \right)^{p/2} (\tanh(kh))^{(2-p_0)/4}
+  S_{\rm ds,break} (\sigma,\theta) = -C'_{\rm ds} \left( \frac{B(k)}{B_{\rm r}} \right)^{p/2} (\tanh(kh))^{(2-p_0)/4}
                                   \sqrt{gk} E(\sigma,\theta)
   \label{eq:satu}
 \end{equation}
@@ -932,7 +909,8 @@
 \begin{equation}
   B(k) = \int_{0}^{2\pi} c_g k^3 E(\sigma,\theta) d\theta
 \end{equation}
-and $B_{\rm r}$ is a threshold saturation level. When $B(k) > B_{\rm r}$, waves break and the exponent $p$ is set equal
+and $B_{\rm r} = 1.75 \times 10^{-3}$ is a threshold saturation level. The proportionality coefficient is set to
+$C'_{\rm ds} = 5.0 \times 10^{-5}$. When $B(k) > B_{\rm r}$, waves break and the exponent $p$ is set equal
 to a calibration parameter $p_0$. For $B(k) \leq B_{\rm r}$ there is no breaking, but some residual dissipation proved
 necessary. This is obtained by setting $p = 0$. A smooth transition between these two situations is achieved by
 (Alves and Banner, 2003);
@@ -941,6 +919,28 @@
   \label{eq:expo}
 \end{equation}
 \noindent
+In Van der Westhuysen (2007) the dissipation modes of breaking and non-breaking waves are separated, so that they
+are active over different parts of the spectrum:
+\begin{equation}
+  S_{\rm ds,w}(\sigma,\theta) = f_{\rm br}(\sigma)S_{\rm ds,break} + \left[ 1 - f_{\rm br}(\sigma) \right]
+  S_{\rm ds,non-break} \ \ ,
+  \label{eq:wctot}
+\end{equation}
+\noindent
+where $S_{\rm ds,break}$ is the contribution by breaking waves (\ref{eq:satu}), and $S_{\rm ds,non-break}$ dissipation
+by means other than breaking (e.g. turbulence). The changeover between the two modes is made with a smooth transition
+function $f_{\rm br}$ similar to (\ref{eq:expo}):
+\begin{equation}
+  f_{\rm br}(\sigma) = \frac{1}{2} + \frac{1}{2} \tanh \left[ 10 \left( \sqrt{\frac{B(k)}{B_{\rm r}}} - 1 \right) \right]
+  \label{eq:fbr}
+\end{equation}
+\noindent
+Since relatively little is known about the dissipation mechanisms of the non-breaking low-frequency waves, their dissipation
+is not modelled in detail. Instead, the expression (\ref{eq3-11}) is used for $S_{\rm ds,non-break}$, to provide general
+background dissipation of non-breaking waves. For this component, the parameter settings of Komen~{\it et~al}. (1984)
+are applied.
+\\[2ex]
+\noindent
 The wind input expression used in saturation-based model is based on that by Yan (1987). This expression embodies experimental
 findings that for strong wind forcing, $u_*/c > 0.1$ say, the wind-induced growth rate of waves depends quadratically on
 $u_*/c$ (e.g. Plant 1982), whereas for weaker forcing, $u_*/c < 0.1$ say, the growth rate depends linearly on $u_*/c$
@@ -989,7 +989,7 @@
 in which $C_{\rm b}$ is a bottom friction coefficient that generally depends on the bottom orbital motion
 represented by $U_{\rm rms}$:
 \begin{equation}
-  U^2_{\rm rms} = \int_{0}^{2\pi} \int_{0}^{\infty} \frac{\sigma^2}{g^2 \sinh^2 kd} E(\sigma,\theta) d\sigma d\theta
+  U^2_{\rm rms} = 2\,\int_{0}^{2\pi} \int_{0}^{\infty} \frac{\sigma^2}{g^2 \sinh^2 kd} E(\sigma,\theta) d\sigma d\theta
   \label{eq3-18}
 \end{equation}
 Hasselmann~{\it et~al}. (1973) found $C_{\rm b}=C_{\rm JON}=0.038$m$^{2}$s$^{-3}$ which is in
@@ -1410,6 +1410,130 @@
   \label{eq3-38}
 \end{equation}
 
+\subsection{First- and second-generation model formulations in SWAN} \label{sec:relaxmod}
+
+The source term $S_{\rm tot}$ for the first- and second-generation formulation (relaxation model) of
+SWAN is (Holthuijsen and de Boer, 1988):
+\begin{equation}
+   S_{\rm tot} =
+    \left\{
+      \begin{array}{lll}
+         S_{\rm in} = A + B E
+         \, , & \mbox{if } E < E_{\rm lim} & \mbox{and } |\theta - \theta_w| < \frac{\pi}{2} \\
+         \\
+         S_{\rm ds,w} = \frac{E_{\rm lim} - E}{\tau}
+         \, , & \mbox{if } E > E_{\rm lim} & \mbox{and } |\theta - \theta_w| < \frac{\pi}{2} \\
+         \\
+         0 \, , & \mbox{if } E > E_{\rm lim} & \mbox{and } |\theta - \theta_w| > \frac{\pi}{2} \\
+      \end{array}
+    \right.
+  \label{eq:relaxmod}
+\end{equation}
+where $S_{\rm in}$ and $S_{\rm ds,w}$ represent input by wind and decay for over-developed sea states,
+respectively, $A$ and $B E$ are a linear and exponential growth term, respectively, $E$ is the
+spectral density, $E_{\rm lim}$ is the saturated spectrum, $\tau$ is a time scale and $\theta$ and
+$\theta_w$ are the discrete spectral wave direction and the wind direction, respectively.
+The expressions for the source terms $S_{\rm in}$ and $S_{\rm ds,w}$ have been modified for shallow water applications
+(N. Booij and L.H. Holthuijsen, personal communication, 1996) and are given below.
+\nocite{Hol88B}
+\\[2ex]
+\noindent
+The distinction between first- and second-generation is only in the formulation of the
+saturated spectrum $E_{\rm lim}$ as outlined below.
+\\[2ex]
+\noindent
+\underline{Linear and exponential growth}\\[2ex]
+The linear growth term $A$ is given by  an expression due to Cavaleri and Malanotte-Rizolli
+(1981) as adapted by Holthuijsen and de Boer (1988) and Holthuijsen et al. (1996):
+\begin{equation}
+   A =
+    \left\{
+      \begin{array}{ll}
+         \frac{\beta_1}{2\pi} \frac{\pi}{g^2} C^2_{\rm drag} \left(  \frac{\rho_a}{\rho_w} \right)^2 (U_{10} \max [0,\cos(\theta-\theta_w)])^4
+         \, , &  \sigma \geq 0.7 \sigma_{\rm PM,d} \\
+         \\
+         0 \, , & \sigma < 0.7 \sigma_{\rm PM,d} \\
+      \end{array}
+    \right.
+\end{equation}
+where $\beta_1$ is a coefficient that has been tuned to be $\beta_1 = 188$,
+$C_{\rm drag}$ is a drag coefficient equal to $C_{\rm drag} = 0.0012$ and
+$\sigma_{\rm PM,d}$ is the fully developed peak frequency including the effect of shallow
+water and is estimated from the depth dependent
+relation of the Shore Protection Manual (1973):
+\begin{equation}
+  \sigma_{\rm PM,d} = \frac{\sigma_{\rm PM}}{\tanh(0.833 {\tilde d}^{0.375})}
+\end{equation}
+with the dimensionless depth
+\begin{equation}
+  {\tilde d} = \frac{gd}{U^2_{10}}
+\end{equation}
+The Pierson-Moskowitz (1964) frequency is
+\begin{equation}
+  \sigma_{\rm PM} = \frac{0.13g}{U_{10}} 2\pi
+\end{equation}
+\nocite{Hol96BB,CER73}
+\\[2ex]
+\noindent
+The exponential growth term $B E$ is due to Snyder et al. (1981) rescaled in terms of $U_{10}$ as
+adapted by Holthuijsen and de Boer (1988) and Holthuijsen et al. (1996):
+\begin{equation}
+  B = \max [0, \beta_2 \frac{5}{2\pi} \frac{\rho_a}{\rho_w} (\frac{U_{10}}{\sigma/k}\cos(\theta-\theta_w) -\beta_3)]\sigma
+\end{equation}
+in which the coefficients $\beta_2$ and $\beta_3$ have been tuned to be $\beta_2 = 0.59$ and $\beta_3 = 0.12$.
+\\[2ex]
+\noindent
+\underline{Decay}\\[2ex]
+If the spectral densities are larger than the wind-dependent saturation spectrum $E_{\rm lim}$, (e.g.,
+when the wind decreases), energy is dissipated with a relaxation model:
+\begin{equation}
+  S_{\rm ds,w} (\sigma,\theta)= \frac{E_{\rm lim}(\sigma,\theta) - E(\sigma,\theta)}{\tau(\sigma)}
+\end{equation}
+where $\tau(\sigma)$ is a time scale given by
+\begin{equation}
+  \tau(\sigma) = \beta_4 \left(\frac{2\pi}{\sigma}\right)^2 \frac{g}{U_{10}\cos(\theta - \theta_w)}
+\end{equation}
+in which the coefficient $\beta_4$ has been tuned to be $\beta_4 = 250$.
+\\[2ex]
+\noindent
+\underline{Saturated spectrum}\\[2ex]
+The saturated spectrum has been formulated in term of wave number with a $\cos^2-$directional
+distribution centred at the local wind direction $\theta_w$. It is essentially an adapted
+Pierson-Moskowitz (1964) spectrum:
+\begin{equation}
+   S_{\rm tot} =
+    \left\{
+      \begin{array}{ll}
+         \frac{\alpha k^{-3}}{2c_g} \exp{ \left \{ -\frac{5}{4}(\frac{\sigma}{\sigma_{\rm PM,d}})^{-4} \right\} }
+         \, \frac{2}{\pi} \, \cos^2(\theta-\theta_w)
+         \, , & \mbox{for } |\theta - \theta_w| < \frac{\pi}{2} \\
+         \\
+         0 \, , & \mbox{for } |\theta - \theta_w| \geq \frac{\pi}{2} \\
+      \end{array}
+    \right.
+\end{equation}
+For the \underline{first-generation} formulation, the scale factor $\alpha$ is a constant and equals
+\begin{equation}
+  \alpha = 0.0081
+\end{equation}
+For the \underline{second-generation} formulation, the scale factor $\alpha$ depends on the total dimensionless
+wave energy ${\tilde E}_{\rm tot,sea}$ of the wind sea part of the spectrum and the dimensionless depth ${\tilde d}$:
+\begin{equation}
+  \alpha = \max [ (0.0081 + (0.013 - 0.0081)e^{-{\tilde d}}), 0.0023 {\tilde E}^{-0.223}_{\rm tot,sea} ]
+\end{equation}
+where the total dimensionless wind sea wave energy ${\tilde E}_{\rm tot,sea}$ is given by:
+\begin{equation}
+  {\tilde E}_{\rm tot,sea} = \frac{g^2 E_{\rm tot,sea}}{U^4_{10}}
+\end{equation}
+with
+\begin{equation}
+  E_{\rm tot,sea} = \int_{\theta_w - \pi/2}^{\theta_w + \pi/2} \int_{0.7 \sigma_{\rm PM,d}}^{\infty}
+  E(\sigma,\theta) d\sigma d\theta
+\end{equation}
+The maximum value of $\alpha$ is taken to be 0.155. This dependency of $\alpha$ on the local
+dimensionless energy of the wind sea permits an overshoot in the wave spectrum under wave
+generation conditions. For deep water $\alpha = 0.0081$ as proposed by Pierson and Moskowitz (1964).
+
 \section{The influence of ambient current on waves} \label{sec:curre}
 
 Waves are subject to the influence of ambient current, when they propagate on it. The ambient current
@@ -4496,6 +4620,10 @@
 2nd International Symposium on Ocean Wave Measurement and Analysis}, New Orleans,
 Louisiana, July 25$-$28, 1993: New York, pp. 630$-$641
 
+\bibitem{Hol96BB}
+Holthuijsen, L.H., N. Booij and L. Bertotti, 1996: The propagation of wind errors through
+ocean wave hindcasts, {\it J. Offshore Mech. And Arctic Eng.}, {\bf 118}, 184-189
+
 \bibitem{Hol97BP}
 Holthuijsen, L.H., N. Booij and R. Padilla-Hernandez, 1997: A curvi-linear, third-generation coastal wave
 model, {\it Conf. Coastal Dynamics '97}, Plymouth, 128-136
@@ -4539,9 +4667,6 @@
 Holthuijsen, L.H., 2005: {\it Waves in oceanic and coastal waters},
 Cambridge University Press, in press.
 
-\bibitem{Hur04V}
-Hurdle, D.P. and G. Ph. van Vledder, 2004. {\it Proc. 23rd Int. Conf. on Offshore Mech. and Arctic. Eng.}
-
 \bibitem{Hsu05OL}
 Hsu, T.-W., S.-H. Ou and J.-M. Liau, 2005: Hindcasting nearshore wind waves using a FEM code for SWAN,
 {\it Coastal Engineering}, {\bf 52}, 177-195
@@ -4734,6 +4859,10 @@
 Ris, R.C., L.H. Holthuijsen and N. Booij, 1994: A spectral model for waves in the near shore zone, {\it Proc.
 24th Int. Conf. Coastal Engng}, Kobe, Japan, pp. 68-78
 
+\bibitem{Ris96H}
+Ris, R.C. and L.H. Holthuijsen, 1996: Spectral modelling of current wave-blocking.
+{\it Proc. 25$^{{\rm th}}$ Int. Conf. Coastal Engng.}, ASCE, 1247-1254.
+
 \bibitem{Ris97H}
 Ris, R.C. and L.H. Holthuijsen, 1997: Modelling of current induced wave-blocking in a spectral wave
 model, {\it 8$^{{\rm th}}$ International Biennal Conference on  Physics of Estuaries and Coastal Seas}, J.
@@ -4887,6 +5016,10 @@
 P.~Wesseling, An introduction to multigrid methods, John Wiley and Sons,
   Chichester, 1992.
 
+\bibitem{Wes07}
+Van der Westhuysen, A.J., 2007: Advances in the spectral modelling of wind waves in the nearshore, {\it Ph.D. thesis},
+Delft University of Technology, Department of Civil Engineering, The Netherlands
+
 \bibitem{Wes07ZB}
 Van der Westhuysen, A.J., M. Zijlema and J.A. Battjes, 2007: Nonlinear saturation-based whitecapping dissipation
 in SWAN for deep and shallow water, {\it Coast. Engng}., {\bf 54}, 151-170
@@ -4971,196 +5104,196 @@
 \begin{theindex}
 \addcontentsline{toc}{chapter}{Index}
   \item ambient, \hyperpage{2}, \hyperpage{7}, \hyperpage{11, 12}, 
-		\hyperpage{29, 30}, \hyperpage{101}
+		\hyperpage{32}, \hyperpage{103}
 
   \indexspace
 
-  \item bathymetry, \hyperpage{23}, \hyperpage{38}, \hyperpage{81}, 
-		\hyperpage{83}, \hyperpage{103}
-  \item bottom, \hyperpage{1--3}, \hyperpage{12}, \hyperpage{14, 15}, 
-		\hyperpage{22, 23}, \hyperpage{30, 31}, 
-		\hyperpage{46, 47}, \hyperpage{83}, \hyperpage{87, 88}, 
-		\hyperpage{98}, \hyperpage{100}, \hyperpage{103--108}
+  \item bathymetry, \hyperpage{24}, \hyperpage{40}, \hyperpage{83}, 
+		\hyperpage{85}, \hyperpage{105}
+  \item bottom, \hyperpage{1--3}, \hyperpage{13--15}, 
+		\hyperpage{22--24}, \hyperpage{32, 33}, 
+		\hyperpage{48, 49}, \hyperpage{85}, \hyperpage{89, 90}, 
+		\hyperpage{100}, \hyperpage{102}, \hyperpage{105--110}
   \item boundary, \hyperpage{3, 4}, \hyperpage{13}, \hyperpage{15}, 
-		\hyperpage{44}, \hyperpage{63, 64}, \hyperpage{66, 67}, 
-		\hyperpage{69--73}, \hyperpage{75}, \hyperpage{82}, 
-		\hyperpage{84}, \hyperpage{88}, \hyperpage{102, 103}
-  \item breaking, \hyperpage{1}, \hyperpage{3}, \hyperpage{12}, 
-		\hyperpage{14--17}, \hyperpage{20, 21}, 
-		\hyperpage{23--25}, \hyperpage{43}, \hyperpage{60}, 
-		\hyperpage{62}, \hyperpage{95}, \hyperpage{98, 99}, 
-		\hyperpage{102}, \hyperpage{107, 108}
+		\hyperpage{46}, \hyperpage{65, 66}, \hyperpage{68, 69}, 
+		\hyperpage{71--75}, \hyperpage{77}, \hyperpage{84}, 
+		\hyperpage{86}, \hyperpage{90}, \hyperpage{104, 105}
+  \item breaking, \hyperpage{1}, \hyperpage{3}, \hyperpage{13--17}, 
+		\hyperpage{20, 21}, \hyperpage{23--25}, \hyperpage{45}, 
+		\hyperpage{62}, \hyperpage{97}, \hyperpage{100, 101}, 
+		\hyperpage{104}, \hyperpage{109}
 
   \indexspace
 
-  \item Cartesian, \hyperpage{3}, \hyperpage{10--12}, \hyperpage{68}, 
-		\hyperpage{85}, \hyperpage{87, 88}
+  \item Cartesian, \hyperpage{3}, \hyperpage{10--12}, \hyperpage{70}, 
+		\hyperpage{87}, \hyperpage{89, 90}
   \item co-ordinate, \hyperpage{3}, \hyperpage{7}, \hyperpage{10--12}, 
-		\hyperpage{27, 28}, \hyperpage{52--54}, \hyperpage{70}, 
-		\hyperpage{72}
-  \item coastal, \hyperpage{1, 2}, \hyperpage{15}, \hyperpage{40}, 
-		\hyperpage{49}, \hyperpage{77}, \hyperpage{81}, 
-		\hyperpage{96, 97}, \hyperpage{99}, 
-		\hyperpage{101--103}, \hyperpage{105, 106}
-  \item convergence, \hyperpage{36, 37}, \hyperpage{42}, \hyperpage{44}, 
-		\hyperpage{49--52}, \hyperpage{75, 76}, 
-		\hyperpage{79, 80}, \hyperpage{88}, \hyperpage{105}, 
-		\hyperpage{109}
-  \item Courant, \hyperpage{40}
+		\hyperpage{27, 28}, \hyperpage{54--56}, \hyperpage{72}, 
+		\hyperpage{74}
+  \item coastal, \hyperpage{1, 2}, \hyperpage{15}, \hyperpage{42}, 
+		\hyperpage{51}, \hyperpage{79}, \hyperpage{83}, 
+		\hyperpage{98, 99}, \hyperpage{101}, 
+		\hyperpage{103--105}, \hyperpage{107, 108}
+  \item convergence, \hyperpage{38, 39}, \hyperpage{44}, \hyperpage{46}, 
+		\hyperpage{51--54}, \hyperpage{77, 78}, 
+		\hyperpage{81, 82}, \hyperpage{90}, \hyperpage{107}, 
+		\hyperpage{111}
+  \item Courant, \hyperpage{42}
   \item current, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{10--13}, 
-		\hyperpage{15}, \hyperpage{29, 30}, \hyperpage{32--34}, 
-		\hyperpage{37, 38}, \hyperpage{41}, \hyperpage{44--47}, 
-		\hyperpage{54}, \hyperpage{59}, \hyperpage{65}, 
-		\hyperpage{78}, \hyperpage{80}, \hyperpage{90}, 
-		\hyperpage{98}, \hyperpage{100, 101}, 
-		\hyperpage{105--108}
-  \item curvi-linear, \hyperpage{2, 3}, \hyperpage{52--54}, 
-		\hyperpage{56}, \hyperpage{66}, \hyperpage{72, 73}, 
-		\hyperpage{81}, \hyperpage{96}, \hyperpage{101}
+		\hyperpage{15}, \hyperpage{32}, \hyperpage{34}, 
+		\hyperpage{36}, \hyperpage{39, 40}, \hyperpage{43}, 
+		\hyperpage{46--49}, \hyperpage{56}, \hyperpage{61}, 
+		\hyperpage{67}, \hyperpage{80}, \hyperpage{82}, 
+		\hyperpage{92}, \hyperpage{100}, \hyperpage{102, 103}, 
+		\hyperpage{107--111}
+  \item curvi-linear, \hyperpage{2, 3}, \hyperpage{54--56}, 
+		\hyperpage{58}, \hyperpage{68}, \hyperpage{74, 75}, 
+		\hyperpage{83}, \hyperpage{98}, \hyperpage{103}
 
   \indexspace
 
-  \item dam, \hyperpage{16}, \hyperpage{30, 31}, \hyperpage{101}, 
-		\hyperpage{109}
-  \item diffraction, \hyperpage{3}, \hyperpage{10}, \hyperpage{30}, 
-		\hyperpage{32, 33}, \hyperpage{92}, \hyperpage{96}, 
-		\hyperpage{102, 103}
-  \item diffusion, \hyperpage{38--41}, \hyperpage{92}
+  \item dam, \hyperpage{16}, \hyperpage{32, 33}, \hyperpage{103}, 
+		\hyperpage{111}
+  \item diffraction, \hyperpage{3}, \hyperpage{10}, \hyperpage{32}, 
+		\hyperpage{34, 35}, \hyperpage{94}, \hyperpage{98}, 
+		\hyperpage{104, 105}
+  \item diffusion, \hyperpage{40--43}, \hyperpage{94}
   \item dissipation, \hyperpage{1}, \hyperpage{3}, \hyperpage{12}, 
 		\hyperpage{14--16}, \hyperpage{20--24}, 
-		\hyperpage{60, 61}, \hyperpage{82}, \hyperpage{95}, 
-		\hyperpage{97}, \hyperpage{100}, \hyperpage{105}, 
-		\hyperpage{107, 108}
+		\hyperpage{62, 63}, \hyperpage{84}, \hyperpage{97}, 
+		\hyperpage{99}, \hyperpage{102}, \hyperpage{107}, 
+		\hyperpage{109, 110}
 
   \indexspace
 
-  \item filter, \hyperpage{13}, \hyperpage{18}, \hyperpage{33}
-  \item flow, \hyperpage{7}, \hyperpage{65}, \hyperpage{97}, 
-		\hyperpage{102}, \hyperpage{104}, \hyperpage{109}
-  \item force, \hyperpage{33, 34}, \hyperpage{54}, \hyperpage{65--67}, 
-		\hyperpage{72, 73}, \hyperpage{83}, \hyperpage{87}, 
-		\hyperpage{90}, \hyperpage{98}
+  \item filter, \hyperpage{13}, \hyperpage{18}, \hyperpage{35}
+  \item flow, \hyperpage{7}, \hyperpage{67}, \hyperpage{99}, 
+		\hyperpage{104}, \hyperpage{106}, \hyperpage{111}
+  \item force, \hyperpage{36}, \hyperpage{56}, \hyperpage{67--69}, 
+		\hyperpage{74, 75}, \hyperpage{85}, \hyperpage{89}, 
+		\hyperpage{92}, \hyperpage{100}
   \item frequency, \hyperpage{7--9}, \hyperpage{11}, \hyperpage{13, 14}, 
-		\hyperpage{16}, \hyperpage{18--21}, \hyperpage{23}, 
-		\hyperpage{25}, \hyperpage{28}, \hyperpage{30}, 
-		\hyperpage{35--37}, \hyperpage{41, 42}, \hyperpage{44}, 
-		\hyperpage{46--50}, \hyperpage{57--60}, 
-		\hyperpage{63, 64}, \hyperpage{95}, \hyperpage{103}, 
-		\hyperpage{109}
-  \item friction, \hyperpage{1}, \hyperpage{3}, \hyperpage{12--15}, 
-		\hyperpage{18}, \hyperpage{21, 22}, \hyperpage{100}, 
-		\hyperpage{102, 103}, \hyperpage{105}, 
-		\hyperpage{107, 108}
+		\hyperpage{16}, \hyperpage{18--23}, \hyperpage{25}, 
+		\hyperpage{28}, \hyperpage{30}, \hyperpage{32}, 
+		\hyperpage{37--39}, \hyperpage{43, 44}, \hyperpage{46}, 
+		\hyperpage{48--52}, \hyperpage{59--62}, 
+		\hyperpage{65, 66}, \hyperpage{97}, \hyperpage{105}, 
+		\hyperpage{111}
+  \item friction, \hyperpage{1}, \hyperpage{3}, \hyperpage{13--15}, 
+		\hyperpage{18}, \hyperpage{22, 23}, \hyperpage{102}, 
+		\hyperpage{104, 105}, \hyperpage{107}, 
+		\hyperpage{109, 110}
 
   \indexspace
 
-  \item garden-sprinkler, \hyperpage{39, 40}, \hyperpage{92}
+  \item garden-sprinkler, \hyperpage{41, 42}, \hyperpage{94}
 
   \indexspace
 
-  \item harbour, \hyperpage{3}, \hyperpage{33}
+  \item harbour, \hyperpage{3}, \hyperpage{35}
 
   \indexspace
 
   \item initial, \hyperpage{13}, \hyperpage{15}, \hyperpage{28}, 
-		\hyperpage{45}, \hyperpage{63, 64}, \hyperpage{76}
-  \item island, \hyperpage{81, 82}
+		\hyperpage{47}, \hyperpage{65, 66}, \hyperpage{78}
+  \item island, \hyperpage{83, 84}
 
   \indexspace
 
-  \item Jonswap, \hyperpage{63, 64}
+  \item Jonswap, \hyperpage{65, 66}
 
   \indexspace
 
   \item latitude, \hyperpage{11, 12}
-  \item limiter, \hyperpage{36}, \hyperpage{49, 50}, \hyperpage{52}
+  \item limiter, \hyperpage{38}, \hyperpage{51, 52}, \hyperpage{54}
   \item longitude, \hyperpage{11, 12}
 
   \indexspace
 
-  \item obstacle, \hyperpage{3}, \hyperpage{30--33}, \hyperpage{38, 39}, 
-		\hyperpage{55--57}
+  \item obstacle, \hyperpage{3}, \hyperpage{32--35}, \hyperpage{40, 41}, 
+		\hyperpage{57--59}
   \item ocean, \hyperpage{1}, \hyperpage{7, 8}, \hyperpage{11}, 
-		\hyperpage{29}, \hyperpage{33}, \hyperpage{40}, 
-		\hyperpage{96}, \hyperpage{99--105}, 
-		\hyperpage{108, 109}
+		\hyperpage{32}, \hyperpage{35}, \hyperpage{42}, 
+		\hyperpage{98}, \hyperpage{101--106}, \hyperpage{108}, 
+		\hyperpage{110, 111}
 
   \indexspace
 
   \item propagation, \hyperpage{1--4}, \hyperpage{10--12}, 
-		\hyperpage{32}, \hyperpage{37}, \hyperpage{39--43}, 
-		\hyperpage{45--47}, \hyperpage{53}, \hyperpage{56}, 
-		\hyperpage{77, 78}, \hyperpage{81}, \hyperpage{83}, 
-		\hyperpage{85--88}, \hyperpage{96}, \hyperpage{98, 99}, 
-		\hyperpage{105}, \hyperpage{107}
+		\hyperpage{34, 35}, \hyperpage{39}, \hyperpage{41--45}, 
+		\hyperpage{47--49}, \hyperpage{55}, \hyperpage{58}, 
+		\hyperpage{79, 80}, \hyperpage{83}, \hyperpage{85}, 
+		\hyperpage{87--90}, \hyperpage{98}, 
+		\hyperpage{100, 101}, \hyperpage{103}, \hyperpage{107}, 
+		\hyperpage{109}
 
   \indexspace
 
-  \item quadruplets, \hyperpage{16}, \hyperpage{25}, \hyperpage{48}
+  \item quadruplets, \hyperpage{16}, \hyperpage{25}, \hyperpage{50}
 
   \indexspace
 
-  \item recti-linear, \hyperpage{73}, \hyperpage{81}
-  \item reflection, \hyperpage{3}, \hyperpage{30}, \hyperpage{32, 33}, 
-		\hyperpage{35}, \hyperpage{56}
+  \item recti-linear, \hyperpage{75}, \hyperpage{83}
+  \item reflection, \hyperpage{3}, \hyperpage{32}, \hyperpage{34, 35}, 
+		\hyperpage{37}, \hyperpage{58}
   \item refraction, \hyperpage{1}, \hyperpage{3}, \hyperpage{11}, 
-		\hyperpage{32}, \hyperpage{38}, \hyperpage{44}, 
-		\hyperpage{46}, \hyperpage{63}, \hyperpage{87}, 
-		\hyperpage{96, 97}, \hyperpage{102, 103}, 
-		\hyperpage{106}, \hyperpage{108}
+		\hyperpage{34}, \hyperpage{40}, \hyperpage{46}, 
+		\hyperpage{48}, \hyperpage{65}, \hyperpage{89}, 
+		\hyperpage{98, 99}, \hyperpage{104, 105}, 
+		\hyperpage{108}, \hyperpage{111}
   \item regular, \hyperpage{2}, \hyperpage{4}, \hyperpage{7}, 
-		\hyperpage{15}, \hyperpage{30}, \hyperpage{87}, 
-		\hyperpage{106}, \hyperpage{108}
+		\hyperpage{15}, \hyperpage{32}, \hyperpage{89}, 
+		\hyperpage{108}, \hyperpage{110, 111}
 
   \indexspace
 
-  \item set-up, \hyperpage{3, 4}, \hyperpage{23}, \hyperpage{33, 34}, 
-		\hyperpage{65}, \hyperpage{73}, \hyperpage{95}
+  \item set-up, \hyperpage{3, 4}, \hyperpage{24}, \hyperpage{36}, 
+		\hyperpage{67}, \hyperpage{75}, \hyperpage{97}
   \item shoaling, \hyperpage{1}, \hyperpage{3}, \hyperpage{15}, 
-		\hyperpage{25}, \hyperpage{98, 99}
-  \item SORDUP, \hyperpage{38}, \hyperpage{40}, \hyperpage{77}
-  \item spherical, \hyperpage{3}, \hyperpage{11, 12}, \hyperpage{88}
-  \item stability, \hyperpage{35, 36}, \hyperpage{40}, \hyperpage{44}, 
-		\hyperpage{50}, \hyperpage{61, 62}, \hyperpage{84}, 
-		\hyperpage{87}
-  \item stationary, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{36--38}, 
-		\hyperpage{40}, \hyperpage{48}, \hyperpage{53}, 
-		\hyperpage{64, 65}, \hyperpage{96}, \hyperpage{98}, 
-		\hyperpage{101}, \hyperpage{109}
-  \item steepness, \hyperpage{14}, \hyperpage{19, 20}, \hyperpage{23}, 
-		\hyperpage{30}
+		\hyperpage{100, 101}
+  \item SORDUP, \hyperpage{40}, \hyperpage{42}, \hyperpage{79}
+  \item spherical, \hyperpage{3}, \hyperpage{11, 12}, \hyperpage{90}
+  \item stability, \hyperpage{37, 38}, \hyperpage{42}, \hyperpage{46}, 
+		\hyperpage{52}, \hyperpage{63}, \hyperpage{86}, 
+		\hyperpage{89}
+  \item stationary, \hyperpage{2, 3}, \hyperpage{7}, \hyperpage{38--40}, 
+		\hyperpage{42}, \hyperpage{50}, \hyperpage{55}, 
+		\hyperpage{66, 67}, \hyperpage{98}, \hyperpage{100}, 
+		\hyperpage{103}, \hyperpage{111}
+  \item steepness, \hyperpage{14}, \hyperpage{19, 20}, \hyperpage{24}, 
+		\hyperpage{32}
   \item swell, \hyperpage{14, 15}, \hyperpage{20}, \hyperpage{22}, 
-		\hyperpage{40, 41}, \hyperpage{97}, \hyperpage{100}, 
-		\hyperpage{108}
+		\hyperpage{42, 43}, \hyperpage{99}, \hyperpage{102}, 
+		\hyperpage{110}
 
   \indexspace
 
   \item triads, \hyperpage{16}
-  \item triangular, \hyperpage{76}, \hyperpage{82}, \hyperpage{84--86}
+  \item triangular, \hyperpage{78}, \hyperpage{84}, \hyperpage{86--88}
 
   \indexspace
 
-  \item unstructured, \hyperpage{4}, \hyperpage{81, 82}, \hyperpage{87}, 
-		\hyperpage{90}, \hyperpage{92}, \hyperpage{96}, 
-		\hyperpage{109}
+  \item unstructured, \hyperpage{4}, \hyperpage{83, 84}, \hyperpage{89}, 
+		\hyperpage{92}, \hyperpage{94}, \hyperpage{98}, 
+		\hyperpage{111}
 
   \indexspace
 
   \item WAM, \hyperpage{1}, \hyperpage{13, 14}, \hyperpage{18--20}, 
-		\hyperpage{26}, \hyperpage{35, 36}, \hyperpage{49}, 
-		\hyperpage{59}, \hyperpage{63}, \hyperpage{99}, 
-		\hyperpage{101}, \hyperpage{106}, \hyperpage{108}
-  \item WAVEWATCH, \hyperpage{1}, \hyperpage{63}, \hyperpage{107}
+		\hyperpage{26}, \hyperpage{37, 38}, \hyperpage{51}, 
+		\hyperpage{61}, \hyperpage{65}, \hyperpage{101}, 
+		\hyperpage{103}, \hyperpage{108}, \hyperpage{110}
+  \item WAVEWATCH, \hyperpage{1}, \hyperpage{65}, \hyperpage{109}
   \item whitecapping, \hyperpage{3}, \hyperpage{12}, \hyperpage{14}, 
-		\hyperpage{16}, \hyperpage{19--21}, \hyperpage{30}, 
-		\hyperpage{100}, \hyperpage{102}, \hyperpage{107, 108}
+		\hyperpage{16}, \hyperpage{19, 20}, \hyperpage{22}, 
+		\hyperpage{32}, \hyperpage{102}, \hyperpage{104}, 
+		\hyperpage{109, 110}
   \item wind, \hyperpage{1--4}, \hyperpage{7}, \hyperpage{12--15}, 
-		\hyperpage{17--22}, \hyperpage{29}, \hyperpage{35}, 
-		\hyperpage{37--41}, \hyperpage{43}, \hyperpage{45}, 
-		\hyperpage{48, 49}, \hyperpage{55}, \hyperpage{64}, 
-		\hyperpage{72}, \hyperpage{81}, \hyperpage{87}, 
-		\hyperpage{95--97}, \hyperpage{100}, 
-		\hyperpage{102--109}
+		\hyperpage{17--22}, \hyperpage{30--32}, \hyperpage{37}, 
+		\hyperpage{39--43}, \hyperpage{45}, \hyperpage{47}, 
+		\hyperpage{50, 51}, \hyperpage{57}, \hyperpage{66}, 
+		\hyperpage{74}, \hyperpage{83}, \hyperpage{89}, 
+		\hyperpage{97--99}, \hyperpage{102--111}
 
 \end{theindex}
 
--- swanuse.tex	2009-12-17 10:30:36.000000000 +0100
+++ swanuse.tex	2009-12-17 10:37:33.000000000 +0100
@@ -32,7 +32,7 @@
 \end{center}
 \vfill
 \begin{center}
-{\Large\bf SWAN Cycle III version 40.72ABCD}
+{\Large\bf SWAN Cycle III version 40.72ABCDE}
 \end{center}
 
 \cleardoublepage
@@ -897,7 +897,6 @@
 {\tt GEN1}     \> SWAN runs in first generation mode\\
 {\tt GEN2}     \> SWAN runs in second generation mode\\
 {\tt GEN3}     \> SWAN runs in third generation mode\\
-%{\tt WCAPPING} \> activates cumulative steepness method for whitecapping\\
 {\tt QUAD}     \> controls the computation of quadruplets\\
 {\tt BREAKING} \> activates dissipation by depth-induced wave breaking\\
 {\tt FRICTION} \> activates dissipation by bottom friction\\
@@ -2415,26 +2414,6 @@
                  Default: {\tt [a]} = 0.0015.\-\\
 \end{tabbing}
 
-%\idxcmd{WCAPPING}
-%\linecmd
-%\begin{verbatim}
-%WCAPping CSM [cst] [pow]
-%\end{verbatim}
-%\linecmd
-%
-%\noindent
-%With this command the user wants to choose the Cumulative Steepness Method (CSM) for approximating whitecapping (see Scientific/Technical documentation)
-%and {\bf not} the formulation of Komen {\it et al.} (1984) and {\bf not} Janssen (1991a).
-%\begin{tabbing}
-%xxxxxxxxxxxx\= \kill
-%{\tt [cst]} \> the tuneable coefficient $C_{wc}^{st}$\+\\
-%               Default: {\tt [cst]} = 4.0.\-\\
-%{\tt [pow]} \> power $m$.\+\\
-%               Default: {\tt [pow]} = 2.0.\-\\
-%\end{tabbing}
-%Note that the CSM method in SWAN is still in its experimental phase. Its results are promising, but the method still suffers some numerical
-%problems.
-%
 \idxcmd{QUADRUPL}
 \linecmd
 \begin{verbatim}
@@ -2488,6 +2467,8 @@
 
 \noindent
 With this command the user can influence depth-induced wave breaking in shallow water in the SWAN model.
+%By default, the bulk dissipation rate is proportional to the energy density per wave component.
+%Optionally, this dissipation rate can also be made frequency dependent as suggested by Chen et al. (1997).
 \\[2ex]
 \noindent
 If this command is not used, \underline{SWAN will account for wave breaking anyhow} (with default options and
@@ -4681,19 +4662,17 @@
 !         | WESTHuysen                  |
 !
 !
-!         | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]  |
-!         |                                                  |
-!         |   JANSsen  [cds1]  [delta] [pwtail]              |
-!         |                                                  |
-!         |   LHIG     [cflhig]                              |
-!         |                                                  |
-!   WCAP <    BJ       [bjstp] [bjalf]                        >
-!         |                                                  |
-!         |   KBJ      [bjstp] [bjalf] [kconv]               |
-!         |                                                  |
-!         |   CSM      [cst] [pow]                           |
-!         |                                                  |
-!         |   AB       [cds2] [br] [p0] [powst] [powk]       |
+!         | -> KOMen   [cds2] [stpm] [powst] [delta] [powk]
+!         |
+!         |   JANSsen  [cds1]  [delta] [pwtail]
+!         |
+!         |   LHIG     [cflhig]
+!         |
+!   WCAP <    BJ       [bjstp] [bjalf]
+!         |
+!         |   KBJ      [bjstp] [bjalf] [kconv]
+!         |
+!         |   AB       [cds2] [br] [p0] [powst] [powk]
 !
 !   QUADrupl [iquad] [limiter] [lambda] [cnl4] [csh1] [csh2] [csh3]
 !
@@ -4704,12 +4683,10 @@
 !        | -> CON [alpha] [gamma]                                      |
 !        |                                                             |
 !        |    VAR [alpha] [gammin] [gammax] [gamneg] [coeff1] [coeff2] |
-!        |                                                             |
-!   BRE <     RUE [alpha] [a] [b]                                       >   &
+!   BRE <                                                               > &
+!        |    RUE [alpha] [a] [b]                                      |
 !        |                                                             |
 !        |    TG  [alpha] [gamma] [pown]                               |
-!        |                                                             |
-!        !    BIP [alpha] [pown] [bref]                                |
 !
 !     ( FREQDep [power] [fmin] [fmax] )
 !
@@ -4857,7 +4834,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72ABCD
+$ Data produced by SWAN version 40.72ABCDE
 $ Project:'projname'     ;   run number:'runnum'
 TIME                                    time-dependent data
      1                                  time coding option
@@ -4963,7 +4940,7 @@
 
 \begin{verbatim}
 SWAN   1                                Swan standard spectral file, version
-$ Data produced by SWAN version 40.72ABCD
+$ Data produced by SWAN version 40.72ABCDE
 $ Project:'projname'     ;   run number:'runnum'
 LOCATIONS                               locations in x-y-space
      1                                  number of locations
@@ -5143,7 +5120,7 @@
 
 \bibitem{Implman}
 {SWAN -- Implementation manual}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.swan.tudelft.nl} (Version 40.72ABCD, August 2009).
+\hl{http://www.swan.tudelft.nl} (Version 40.72ABCDE, December 2009).
 
 \bibitem{Progrul}
 {SWAN -- Programming rules}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
@@ -5151,7 +5128,7 @@
 
 \bibitem{Techdoc}
 {SWAN -- Scientific and Technical documentation}. Delft University of Technology, Environmental Fluid Mechanics Section, available from
-\hl{http://www.swan.tudelft.nl} (Version 40.72ABC, August 2009).
+\hl{http://www.swan.tudelft.nl} (Version 40.72ABCDE, December 2009).
 
 \end{thebibliography}
 
@@ -5161,29 +5138,28 @@
 
   \indexspace
 
-  \item bathymetry, \hyperpage{5}, \hyperpage{10, 11}, \hyperpage{65}, 
-		\hyperpage{86}
-  \item BLOCK, \hyperpage{75}
+  \item bathymetry, \hyperpage{5}, \hyperpage{10, 11}, \hyperpage{64}, 
+		\hyperpage{85}
+  \item BLOCK, \hyperpage{74}
   \item bottom, \hyperpage{3}, \hyperpage{8--13}, \hyperpage{16}, 
 		\hyperpage{18, 19}, \hyperpage{21, 22}, \hyperpage{24}, 
 		\hyperpage{34, 35}, \hyperpage{38--40}, \hyperpage{53}, 
-		\hyperpage{55, 56}, \hyperpage{59}, \hyperpage{66, 67}, 
+		\hyperpage{55, 56}, \hyperpage{58}, \hyperpage{66, 67}, 
 		\hyperpage{69}, \hyperpage{79, 80}, \hyperpage{86}, 
 		\hyperpage{92}, \hyperpage{94}
   \item BOUND SHAPE, \hyperpage{41}
   \item boundary, \hyperpage{1}, \hyperpage{3}, \hyperpage{6--14}, 
 		\hyperpage{22}, \hyperpage{26}, \hyperpage{28}, 
 		\hyperpage{30}, \hyperpage{32}, \hyperpage{42--51}, 
-		\hyperpage{53}, \hyperpage{62}, \hyperpage{66}, 
-		\hyperpage{70}, \hyperpage{72}, \hyperpage{83, 84}, 
-		\hyperpage{87}
+		\hyperpage{53}, \hyperpage{61}, \hyperpage{66}, 
+		\hyperpage{70--72}, \hyperpage{83, 84}, \hyperpage{87}
   \item BOUNDNEST1, \hyperpage{46}
   \item BOUNDNEST2, \hyperpage{47}
   \item BOUNDNEST3, \hyperpage{49}
   \item BOUNDSPEC, \hyperpage{42}
   \item BREAKING, \hyperpage{54}
   \item breaking, \hyperpage{11}, \hyperpage{16--18}, \hyperpage{22}, 
-		\hyperpage{53--55}, \hyperpage{57}, \hyperpage{62}, 
+		\hyperpage{53, 54}, \hyperpage{56, 57}, \hyperpage{61}, 
 		\hyperpage{80}, \hyperpage{86}, \hyperpage{92}
 
   \indexspace
@@ -5192,29 +5168,29 @@
 		\hyperpage{26, 27}, \hyperpage{29}, \hyperpage{36}, 
 		\hyperpage{41}, \hyperpage{44}, \hyperpage{46}, 
 		\hyperpage{48--51}, \hyperpage{59}, \hyperpage{67--71}, 
-		\hyperpage{74}, \hyperpage{78}, \hyperpage{86}, 
+		\hyperpage{73, 74}, \hyperpage{78}, \hyperpage{86}, 
 		\hyperpage{90}, \hyperpage{94}, \hyperpage{106}, 
 		\hyperpage{108, 109}
   \item CGRID, \hyperpage{28}
   \item co-ordinate, \hyperpage{32}
   \item coastal, \hyperpage{3--5}, \hyperpage{15}
-  \item COMPUTE, \hyperpage{87}
-  \item convergence, \hyperpage{6, 7}, \hyperpage{64}
+  \item COMPUTE, \hyperpage{86}
+  \item convergence, \hyperpage{6, 7}, \hyperpage{63}
   \item COORDINATES, \hyperpage{27}
   \item Courant, \hyperpage{14}
   \item current, \hyperpage{3}, \hyperpage{5, 6}, \hyperpage{9--12}, 
 		\hyperpage{14}, \hyperpage{16}, \hyperpage{19--21}, 
 		\hyperpage{26, 27}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38}, \hyperpage{45}, \hyperpage{52}, 
-		\hyperpage{65}, \hyperpage{67}, \hyperpage{79, 80}, 
-		\hyperpage{84}, \hyperpage{90}
+		\hyperpage{64, 65}, \hyperpage{67}, \hyperpage{79, 80}, 
+		\hyperpage{83}, \hyperpage{90}
   \item CURVE, \hyperpage{68}
   \item curvi-linear, \hyperpage{3}, \hyperpage{6}, \hyperpage{8, 9}, 
 		\hyperpage{13, 14}, \hyperpage{16}, \hyperpage{21, 22}, 
 		\hyperpage{29}, \hyperpage{32}, \hyperpage{34, 35}, 
 		\hyperpage{38, 39}, \hyperpage{43}, \hyperpage{47, 48}, 
-		\hyperpage{50}, \hyperpage{62}, \hyperpage{66}, 
-		\hyperpage{68}
+		\hyperpage{50}, \hyperpage{62}, \hyperpage{65}, 
+		\hyperpage{67}
 
   \indexspace
 
@@ -5223,21 +5199,21 @@
   \item diffraction, \hyperpage{5}, \hyperpage{22}, \hyperpage{60, 61}, 
 		\hyperpage{64}
   \item diffusion, \hyperpage{62}, \hyperpage{64, 65}
-  \item dissipation, \hyperpage{22}, \hyperpage{52--56}, \hyperpage{62}, 
-		\hyperpage{80}, \hyperpage{92}
+  \item dissipation, \hyperpage{22}, \hyperpage{52, 53}, \hyperpage{55}, 
+		\hyperpage{61}, \hyperpage{79, 80}, \hyperpage{92}
 
   \indexspace
 
-  \item filter, \hyperpage{13, 14}, \hyperpage{60, 61}
+  \item filter, \hyperpage{13, 14}, \hyperpage{60}
   \item flow, \hyperpage{6}, \hyperpage{30, 31}, \hyperpage{34}, 
-		\hyperpage{36, 37}, \hyperpage{65}
+		\hyperpage{36, 37}, \hyperpage{64}
   \item force, \hyperpage{80}, \hyperpage{93}
-  \item FRAME, \hyperpage{67}
+  \item FRAME, \hyperpage{66}
   \item frequency, \hyperpage{4}, \hyperpage{6}, \hyperpage{8}, 
 		\hyperpage{14, 15}, \hyperpage{26}, \hyperpage{30, 31}, 
 		\hyperpage{42}, \hyperpage{45}, \hyperpage{51, 52}, 
-		\hyperpage{54--56}, \hyperpage{62}, \hyperpage{64--66}, 
-		\hyperpage{72--74}, \hyperpage{79}, \hyperpage{83, 84}, 
+		\hyperpage{54--56}, \hyperpage{61}, \hyperpage{64--66}, 
+		\hyperpage{72--74}, \hyperpage{78, 79}, \hyperpage{83}, 
 		\hyperpage{86}, \hyperpage{89}, \hyperpage{91}, 
 		\hyperpage{108, 109}
   \item FRICTION, \hyperpage{55}
@@ -5249,7 +5225,7 @@
 
   \indexspace
 
-  \item garden-sprinkler, \hyperpage{62, 63}
+  \item garden-sprinkler, \hyperpage{62}
   \item GEN1, \hyperpage{52}
   \item GEN2, \hyperpage{52}
   \item GEN3, \hyperpage{53}
@@ -5278,8 +5254,8 @@
 
   \item latitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{86}, 
 		\hyperpage{108}
-  \item LIMITER, \hyperpage{57}
-  \item limiter, \hyperpage{5}, \hyperpage{7}, \hyperpage{57}, 
+  \item LIMITER, \hyperpage{56}
+  \item limiter, \hyperpage{5}, \hyperpage{7}, \hyperpage{56}, 
 		\hyperpage{64, 65}
   \item longitude, \hyperpage{8}, \hyperpage{27}, \hyperpage{49, 50}, 
 		\hyperpage{86}, \hyperpage{108}
@@ -5293,7 +5269,7 @@
   \item nautical, \hyperpage{7, 8}, \hyperpage{108, 109}
   \item NESTOUT, \hyperpage{84}
   \item NGRID, \hyperpage{70}
-  \item NUMERIC, \hyperpage{63}
+  \item NUMERIC, \hyperpage{62}
 
   \indexspace
 
@@ -5306,9 +5282,9 @@
 
   \indexspace
 
-  \item POINTS, \hyperpage{70}
+  \item POINTS, \hyperpage{69}
   \item PROJECT, \hyperpage{24}
-  \item PROP, \hyperpage{62}
+  \item PROP, \hyperpage{61}
   \item propagation, \hyperpage{6}, \hyperpage{12}, \hyperpage{22}, 
 		\hyperpage{57}, \hyperpage{61, 62}, \hyperpage{79}, 
 		\hyperpage{92}
@@ -5317,29 +5293,29 @@
 
   \item QUADRUPL, \hyperpage{53}
   \item quadruplets, \hyperpage{15, 16}, \hyperpage{22}, \hyperpage{54}, 
-		\hyperpage{57}, \hyperpage{80}, \hyperpage{92}
+		\hyperpage{56}, \hyperpage{79}, \hyperpage{92}
   \item QUANTITY, \hyperpage{72}
 
   \indexspace
 
-  \item RAY, \hyperpage{69}
+  \item RAY, \hyperpage{68}
   \item READGRID COORDINATES, \hyperpage{32}
   \item READGRID UNSTRUCTURED, \hyperpage{32}
   \item READINP, \hyperpage{37}
   \item recti-linear, \hyperpage{3}, \hyperpage{16}, \hyperpage{29}, 
-		\hyperpage{38}, \hyperpage{68}
-  \item reflection, \hyperpage{58, 59}
-  \item refraction, \hyperpage{6}, \hyperpage{11}, \hyperpage{62}, 
-		\hyperpage{64, 65}, \hyperpage{69}
+		\hyperpage{38}, \hyperpage{67}
+  \item reflection, \hyperpage{57--59}
+  \item refraction, \hyperpage{6}, \hyperpage{11}, \hyperpage{61}, 
+		\hyperpage{64}, \hyperpage{69}
   \item regular, \hyperpage{4}, \hyperpage{8, 9}, \hyperpage{16}, 
 		\hyperpage{22}, \hyperpage{28, 29}, \hyperpage{31, 32}, 
-		\hyperpage{34, 35}, \hyperpage{62}, \hyperpage{66, 67}
+		\hyperpage{34, 35}, \hyperpage{62}, \hyperpage{65, 66}
 
   \indexspace
 
   \item SET, \hyperpage{25}
   \item set-up, \hyperpage{4}, \hyperpage{6}, \hyperpage{16, 17}, 
-		\hyperpage{22}, \hyperpage{28}, \hyperpage{60}, 
+		\hyperpage{22}, \hyperpage{28}, \hyperpage{59, 60}, 
 		\hyperpage{65}
   \item SETUP, \hyperpage{59}
   \item shoaling, \hyperpage{18}
@@ -5349,20 +5325,20 @@
   \item spherical, \hyperpage{3--5}, \hyperpage{8, 9}, \hyperpage{21}, 
 		\hyperpage{25}, \hyperpage{27--29}, \hyperpage{31}, 
 		\hyperpage{36}, \hyperpage{45}, \hyperpage{48--50}, 
-		\hyperpage{59, 60}, \hyperpage{67--71}, \hyperpage{81}, 
+		\hyperpage{59}, \hyperpage{67--71}, \hyperpage{81}, 
 		\hyperpage{86}, \hyperpage{103}, \hyperpage{108}
   \item stability, \hyperpage{14}, \hyperpage{62}
   \item stationary, \hyperpage{4}, \hyperpage{6}, \hyperpage{8}, 
 		\hyperpage{14}, \hyperpage{16}, \hyperpage{18}, 
 		\hyperpage{20, 21}, \hyperpage{26, 27}, \hyperpage{34}, 
 		\hyperpage{37--39}, \hyperpage{45, 46}, \hyperpage{51}, 
-		\hyperpage{53}, \hyperpage{62--64}, \hyperpage{81}, 
-		\hyperpage{87, 88}, \hyperpage{103}, \hyperpage{106}, 
+		\hyperpage{53}, \hyperpage{62--64}, \hyperpage{80}, 
+		\hyperpage{87}, \hyperpage{103}, \hyperpage{106}, 
 		\hyperpage{108, 109}
   \item steepness, \hyperpage{53}, \hyperpage{80}, \hyperpage{93}
   \item STOP, \hyperpage{88}
-  \item swell, \hyperpage{12}, \hyperpage{15}, \hyperpage{56}, 
-		\hyperpage{63}, \hyperpage{73, 74}, \hyperpage{78}, 
+  \item swell, \hyperpage{12}, \hyperpage{15}, \hyperpage{55}, 
+		\hyperpage{62}, \hyperpage{73, 74}, \hyperpage{78}, 
 		\hyperpage{89}
 
   \indexspace
@@ -5370,7 +5346,7 @@
   \item TABLE, \hyperpage{81}
   \item TEST, \hyperpage{85}
   \item TRIAD, \hyperpage{56}
-  \item triads, \hyperpage{16}, \hyperpage{56}, \hyperpage{80}, 
+  \item triads, \hyperpage{16}, \hyperpage{56}, \hyperpage{79}, 
 		\hyperpage{92}
   \item triangular, \hyperpage{3}, \hyperpage{16}
 
@@ -5379,15 +5355,15 @@
   \item unstructured, \hyperpage{3, 4}, \hyperpage{8, 9}, 
 		\hyperpage{12}, \hyperpage{21}, \hyperpage{29}, 
 		\hyperpage{32}, \hyperpage{34, 35}, \hyperpage{44}, 
-		\hyperpage{58}, \hyperpage{61--64}, \hyperpage{66}, 
-		\hyperpage{71}, \hyperpage{78}
+		\hyperpage{58}, \hyperpage{60}, \hyperpage{62--64}, 
+		\hyperpage{66}, \hyperpage{71}, \hyperpage{78}
 
   \indexspace
 
   \item WAM, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{18}, \hyperpage{22}, \hyperpage{37}, 
 		\hyperpage{46}, \hyperpage{48, 49}, \hyperpage{54}, 
-		\hyperpage{81}, \hyperpage{83--85}, \hyperpage{87}
+		\hyperpage{81--84}, \hyperpage{87}
   \item WAVEWATCH, \hyperpage{4--7}, \hyperpage{9, 10}, \hyperpage{15}, 
 		\hyperpage{22}, \hyperpage{46}, \hyperpage{49, 50}
   \item whitecapping, \hyperpage{7}, \hyperpage{16}, \hyperpage{53}, 
@@ -5397,9 +5373,9 @@
   \item wind, \hyperpage{3}, \hyperpage{5}, \hyperpage{7--12}, 
 		\hyperpage{14--19}, \hyperpage{21}, \hyperpage{26}, 
 		\hyperpage{34, 35}, \hyperpage{39}, \hyperpage{41}, 
-		\hyperpage{51--53}, \hyperpage{56}, \hyperpage{61}, 
-		\hyperpage{63}, \hyperpage{65}, \hyperpage{79}, 
-		\hyperpage{86, 87}, \hyperpage{92}, \hyperpage{94}
+		\hyperpage{51--53}, \hyperpage{55}, \hyperpage{61, 62}, 
+		\hyperpage{64, 65}, \hyperpage{79}, \hyperpage{86, 87}, 
+		\hyperpage{92}, \hyperpage{94}
 
 \end{theindex}
 
--- swmod1.ftn	2009-12-17 10:30:36.000000000 +0100
+++ swmod1.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -1982,7 +1982,7 @@
 !                 =2; for command BRE VAR ..., surf breaking according to Nelson
 !                 =3; for command BRE RUE ..., surf breaking according to Ruessink
 !                 =4; for command BRE TG ... , surf breaking according to Thornton and Guza
-!                 =5; for command BRE BIP ..., surf breaking according to biphase scaling
+!                 =5; for command BRE WESTH ..., surf breaking according to biphase scaling
 ! ITERMX [      ] maximum number of iterations:
 !                 is set equal to MXITST in case of stationary computations
 !                 is set equal to MXITNS in case of nonstationary computations
--- swmod2.ftn	2009-12-17 10:30:32.000000000 +0100
+++ swmod2.ftn	2009-12-17 10:31:32.000000000 +0100
@@ -203,6 +203,8 @@
 
       CHARACTER (LEN=1)  :: OUT_COMMENT = '%' ! comment sign for heading lines
 
+      LOGICAL            :: LCOMPGRD
+
       ! formats for output:
       CHARACTER (LEN=40) :: FLT_BLOCK = '(6E12.4)'       ! floating point block
       CHARACTER (LEN=40) :: FLT_TABLE = '(E11.4)'        ! floating point table
