      SUBROUTINE RDMEMS(ISUB,IRET)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    RDMEMS
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS SUBROUTINE READS A PARTICULAR SUBSET FROM A BUFR
C   MESSAGE IN INTERNAL MEMORY (ARRAY MBAY IN COMMON BLOCK /BITBUF/)
C   INTO INTERNAL SUBSET ARRAYS BASED ON THE SUBSET NUMBER IN THE
C   MESSAGE.
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR
C 1998-07-08  J. WOOLLEN -- REPLACED CALL TO CRAY LIBRARY ROUTINE
C                           "ABORT" WITH CALL TO NEW INTERNAL BUFRLIB
C                           ROUTINE "BORT"
C 1998-10-27  J. WOOLLEN -- MODIFIED TO CORRECT PROBLEMS CAUSED BY IN-
C                           LINING CODE WITH FPP DIRECTIVES
C 1999-11-18  J. WOOLLEN -- THE NUMBER OF BUFR FILES WHICH CAN BE
C                           OPENED AT ONE TIME INCREASED FROM 10 TO 32
C                           (NECESSARY IN ORDER TO PROCESS MULTIPLE
C                           BUFR FILES UNDER THE MPI)
C 2000-09-19  J. WOOLLEN -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           10,000 TO 20,000 BYTES
C 2001-08-15  D. KEYSER  -- PARAMETER MAXMEM (THE MAXIMUM NUMBER OF
C                           BYTES REQUIRED TO STORE ALL MESSAGES
C                           INTERNALLY) WAS INCREASED FROM 8 MBYTES TO
C                           16 MBYTES
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- PARAMETER MAXMSG (THE MAXIMUM NUMBER OF
C                           BUFR MESSAGES WHICH CAN BE STORED
C                           INTERNALLY) INCREASED FROM 50000 TO 200000;
C                           UNIFIED/PORTABLE FOR WRF; ADDED
C                           DOCUMENTATION (INCLUDING HISTORY); OUTPUTS
C                           MORE COMPLETE DIAGNOSTIC INFO WHEN ROUTINE
C                           TERMINATES ABNORMALLY OR UNUSUAL THINGS
C                           HAPPEN
C 2004-08-09  J. ATOR    -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           20,000 TO 50,000 BYTES
C 2004-11-15  D. KEYSER  -- PARAMETER MAXMEM (THE MAXIMUM NUMBER OF
C                           BYTES REQUIRED TO STORE ALL MESSAGES
C                           INTERNALLY) WAS INCREASED FROM 16 MBYTES TO
C                           50 MBYTES
C 2009-04-21  J. ATOR    -- USE ERRWRT
C
C USAGE:    CALL RDMEMS (ISUB, IRET)
C   INPUT ARGUMENT LIST:
C     ISUB     - INTEGER: POINTER TO SUBSET NUMBER TO READ IN BUFR
C                MESSAGE
C
C   OUTPUT ARGUMENT LIST:
C     IRET     - INTEGER: RETURN CODE:
C                       0 = normal return
C                      -1 = ISUB is greater than the number of subsets
C                           in memory
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     ERRWRT   IUPB     READSB
C                               STATUS
C    THIS ROUTINE IS CALLED BY: UFBMMS   UFBMNS   UFBRMS
C                               Normally not called by any application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      CHARACTER*128 BORT_STR,ERRSTR

      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM),
     .                MDX(MXDXW),IPDXM(MXDXM),LDXM,NDXM,LDXTS,NDXTS,
     .                IFDXTS(MXDXTS),ICDXTS(MXDXTS),IPMSGS(MXDXTS)
      COMMON /MSGCWD/ NMSG(NFILES),NSUB(NFILES),MSUB(NFILES),
     .                INODE(NFILES),IDATE(NFILES)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(MXMSGLD4),MBYT(NFILES),
     .                MBAY(MXMSGLD4,NFILES)
      COMMON /UNPTYP/ MSGUNP(NFILES)
      COMMON /QUIET / IPRT

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE MESSAGE REQUEST AND FILE STATUS
C  -----------------------------------------

      CALL STATUS(MUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IF(IM.EQ.0) GOTO 902
      IF(NSUB(LUN).NE.0) GOTO 903

      IF(ISUB.GT.MSUB(LUN)) THEN
         IF(IPRT.GE.0) THEN
      CALL ERRWRT('+++++++++++++++++++++WARNING+++++++++++++++++++++++')
           WRITE ( UNIT=ERRSTR, FMT='(A,I5,A,A,I5,A)' )
     .      'BUFRLIB: RDMEMS - REQ. SUBSET #', ISUB, ' (= 1st INPUT ',
     .      'ARG.) > # OF SUBSETS IN MEMORY MESSAGE (', MSUB(LUN), ')'
           CALL ERRWRT(ERRSTR)
           CALL ERRWRT('RETURN WITH IRET = -1')
      CALL ERRWRT('+++++++++++++++++++++WARNING+++++++++++++++++++++++')
      CALL ERRWRT(' ')
         ENDIF
         IRET = -1
         GOTO 100
      ENDIF

      MBYM = MBYT(LUN)
      NBYT = 0

C  POSITION TO SUBSET NUMBER ISUB IN MEMORY MESSAGE
C  ------------------------------------------------

      IF(MSGUNP(LUN).EQ.0) THEN
         NSUB(LUN) = ISUB-1
         DO I=1,ISUB-1
         MBYT(LUN) = MBYT(LUN) + IUPB(MBAY(1,LUN),MBYT(LUN)+1,16)
         ENDDO
      ELSEIF(MSGUNP(LUN).EQ.1) THEN
c  .... message with "standard" Section 3
         DO I=1,ISUB-1
         CALL READSB(MUNIT,IRET)
         ENDDO
      ELSEIF(MSGUNP(LUN).EQ.2) THEN
c  .... compressed message
         NSUB(LUN) = ISUB-1
      ENDIF

C  NOW READ SUBSET NUMBER ISUB FROM MEMORY MESSAGE
C  -----------------------------------------------

      CALL READSB(MUNIT,IRET)
c  .... This should have already been accounted for with stmt. 902 or
c       IRET = -1 above
      IF(IRET.NE.0) GOTO 904

C  RESET SUBSET POINTER BACK TO ZERO (BEGINNING OF MESSAGE) AND RETURN
C  -------------------------------------------------------------------

      MBYT(LUN) = MBYM
      NSUB(LUN) = 0

C  EXITS
C  -----

100   RETURN
900   CALL BORT('BUFRLIB: RDMEMS - INPUT BUFR FILE IS CLOSED, IT '//
     . 'MUST BE OPEN FOR INPUT')
901   CALL BORT('BUFRLIB: RDMEMS - INPUT BUFR FILE IS OPEN FOR '//
     . 'OUTPUT, IT MUST BE OPEN FOR INPUT')
902   CALL BORT('BUFRLIB: RDMEMS - A MEMORY MESSAGE MUST BE OPEN IN '//
     . 'INPUT BUFR FILE, NONE ARE')
903   WRITE(BORT_STR,'("BUFRLIB: RDMEMS - UPON ENTRY, SUBSET POINTER '//
     . 'IN MEMORY MESSAGE IS NOT AT BEGINNING (",I3," SUBSETS HAVE '//
     . 'BEEN READ, SHOULD BE 0)")') NSUB(LUN)
      CALL BORT(BORT_STR)
904   CALL BORT('BUFRLIB: RDMEMS - CALL TO ROUTINE READSB RETURNED '//
     . 'WITH IRET = -1 (EITHER MEMORY MESSAGE NOT OPEN OR ALL '//
     . 'SUBSETS IN MESSAGE READ')
      END
