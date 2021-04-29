#include "cppdefs.h"
      MODULE stdinp_mod
!
!svn $Id: stdinp_mod.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to read requested KeyWord     !
!  parameter from ROMS standard input file. It is used to process      !
!  specific parameters before the normal call to read all of them.     !
!                                                                      !
!  getpar_i      Reads requested integer parameter.                    !
!                                                                      !
!  getpar_l      Reads requested logical parameter.                    !
!                                                                      !
!  getpar_r      Reads requested floating-point (real) parameter.      !
!                                                                      !
!  getpar_r      Reads requested string parameter.                     !
!                                                                      !
!  stdinp_unit   Determines input unit. If distributed-memory, it      !
!                  opens standard input file for reading.              !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE inp_decode_mod
!
      INTERFACE getpar_i
        MODULE PROCEDURE getpar_0d_i
        MODULE PROCEDURE getpar_1d_i
      END INTERFACE getpar_i

      INTERFACE getpar_l
        MODULE PROCEDURE getpar_0d_l
        MODULE PROCEDURE getpar_1d_l
      END INTERFACE getpar_l

      INTERFACE getpar_r
        MODULE PROCEDURE getpar_0d_r
        MODULE PROCEDURE getpar_1d_r
      END INTERFACE getpar_r

      INTERFACE getpar_s
        MODULE PROCEDURE getpar_0d_s
      END INTERFACE getpar_s
!
      CONTAINS
!
      FUNCTION stdinp_unit (localPET, GotFile) RESULT (InpUnit)
!
!***********************************************************************
!                                                                      !
!  This function determines ROMS standard input unit to process its    !
!  parameters. If distributed-memory, it gets standard input filename  !
!  and open it for processing.                                         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     GotFile    Standard input filename is known (logical)            !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits,  ONLY : Iname, SourceFile, stdinp, stdout
      USE mod_scalars,  ONLY : exit_flag

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasts
#endif
!
!  Imported variable declarations.
!
      logical, intent(out) :: GotFile
      integer, intent(in)  :: localPET
!
!  Local variable declararions
!
      integer :: InpUnit
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", stdout_unit"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Determine ROMS standard input unit.
!-----------------------------------------------------------------------
#ifdef DISTRIBUTE
!
!  In distributed-memory configurations, the input physical parameters
!  script is opened as a regular file.  It is read and processed by all
!  parallel nodes.  This is to avoid very complex broadcasting of the
!  input parameters to all nodes.
# ifdef MODEL_COUPLING
!  However, in ESM coupling, the submission argument is 'coupling.in'
!  and the Iname was read and processed in the coupling configuration
!  elsewhere.
# endif
!
      InpUnit=1
# ifndef MODEL_COUPLING
      IF (localPET.eq.0) CALL my_getarg (1, Iname)
      CALL mp_bcasts (1, 1, Iname)
# endif
      OPEN (InpUnit, FILE=TRIM(Iname), FORM='formatted', STATUS='old',  &
     &      ERR=10)
      GotFile=.TRUE.
      GO TO 30
 10   IF (localPET.eq.0) WRITE (stdout,20)
      exit_flag=2
      RETURN
 20   FORMAT (/,' STDINP_UNIT - Unable to open ROMS/TOMS input script', &
     &                        ' file.',                                 &
     &        /,11x,'In distributed-memory applications, the input',    &
     &        /,11x,'script file is processed in parallel. The Unix',   &
     &        /,11x,'routine GETARG is used to get script file name.',  &
     &        /,11x,'For example, in MPI applications make sure that',  &
     &        /,11x,'command line is something like:',/,                &
     &        /,11x,'mpirun -np 4 romsM roms.in',/,/,11x,'and not',/,   &
     &        /,11x,'mpirun -np 4 romsM < roms.in',/)
 30   CONTINUE
#else
      InpUnit=stdinp
      GotFile=.FALSE.
#endif
      END FUNCTION stdinp_unit
!
      SUBROUTINE getpar_0d_i (localPET, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a scalar integer parameter from ROMS standard input file.     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (integer)              !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      integer, intent(out) :: Value
!
      character (len=*), intent(in) :: KeyWord
      character (len=*), intent(in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile
!
      integer :: InpUnit, Npts, Nval, status
      integer :: Ivalue(1)
!
      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input integer parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_0D_i - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            Npts=load_i(Nval, Rval, 1, Ivalue)
            Value=Ivalue(1)
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_0D_I - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_0D_I - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_0d_i
!
      SUBROUTINE getpar_1d_i (localPET, Ndim, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a 1D integer parameter from ROMS standard input file.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     Ndim       Size integer variable dimension                       !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (integer 1D array)     !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      integer, intent(in)  :: Ndim
      integer, intent(out) :: Value(:)
!
      character (len=*), intent(in) :: KeyWord
      character (len=*), intent(in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile
!
      integer :: InpUnit, Npts, Nval, status
!
      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input 1D integer parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_1D_i - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            Npts=load_i(Nval, Rval, Ndim, Value)
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_1D_I - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_1D_I - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_1d_i
!
      SUBROUTINE getpar_0d_l (localPET, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a scalar logical parameter from ROMS standard input file.     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (logical)              !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      logical, intent(out) :: Value
!
      character (len=*), intent(in) :: KeyWord
      character (len=*), intent(in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile
      logical :: Lvalue(1)
!
      integer :: InpUnit, Npts, Nval, status
!
      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input logical parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_0D_L - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            Npts=load_l(Nval, Cval, 1, Lvalue)
            Value=Lvalue(1)
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_0D_L - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_0D_L - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_0d_l
!
      SUBROUTINE getpar_1d_l (localPET, Ndim, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a 1D logical parameter from ROMS standard input file.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     Ndim       Size logical variable dimension                       !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (logical 1D array)     !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      logical, intent(out) :: Value(:)
      integer, intent(in)  :: localPET
      integer, intent(in)  :: Ndim
!
      character (len=*), intent(in) :: KeyWord
      character (len=*), intent(in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile
!
      integer :: InpUnit, Npts, Nval, status
!
      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input 1D logical parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_1D_L - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            Npts=load_l(Nval, Cval, Ndim, Value)
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_1D_L - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_1D_L - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_1d_l
!
      SUBROUTINE getpar_0d_r (localPET, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a scalar floating-point parameter from ROMS standard input    !
!  file.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (real)                 !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
      real(r8), intent(out) :: Value
!
      character (len=*), intent(in) :: KeyWord
      character (len=*), intent(in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile

      integer :: InpUnit, Npts, Nval, status

      real(r8) :: Rvalue(1)

      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input floating-point parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_0D_R - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            Npts=load_r(Nval, Rval, 1, Rvalue)
            Value=Rvalue(1)
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_0D_R - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_0D_R - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_0d_r
!
      SUBROUTINE getpar_1d_r (localPET, Ndim, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a 1D floating-point parameter from ROMS standard input file.  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     Ndim       Size integer variable dimension                       !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (real 1D array)        !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: Ndim

      real(r8), intent(out) :: Value(:)
!
      character (len=*), intent(in) :: KeyWord
      character (len=*), intent(in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile
!
      integer :: InpUnit, Npts, Nval, status
!
      real(dp), dimension(nRval) :: Rval
!
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input 1D floating-point parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_1D_R - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            Npts=load_r(Nval, Rval, Ndim, Value)
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_1D_R - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_1D_R - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_1d_r
!
      SUBROUTINE getpar_0d_s (localPET, Value, KeyWord, InpName)
!
!***********************************************************************
!                                                                      !
!  Reads a scalar string from ROMS standard input file.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     localPET   Local Persistent Execution Thread (integer)           !
!     KeyWord    Keyword associated with input parameter (string)      !
!     InpName    Standard input filename (string; OPTIONAL)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Value      Standard input parameter value (string)               !
!                                                                      !
!***********************************************************************
!
      USE mod_iounits, ONLY : stdout
      USE mod_scalars, ONLY : exit_flag
!
      USE strings_mod, ONLY : FoundError
!
!  Imported variable declarations.
!
      integer, intent(in)  :: localPET
!
      character (len=*), intent( in) :: KeyWord
      character (len=*), intent(out) :: Value
      character (len=*), intent( in), optional :: InpName
!
!  Local variable declarations.
!
      logical :: foundit, GotFile
!
      integer :: InpUnit, Npts, Nval, i, status
!
      real(dp), dimension(nRval) :: Rval
!
      character (len=1  ), parameter :: blank = ' '
      character (len= 40) :: string
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Read requested ROMS standard input integer parameter.
!-----------------------------------------------------------------------
!
!  Get standard input unit.
!
      IF (PRESENT(InpName)) THEN
        InpUnit=1
        OPEN (InpUnit, FILE=TRIM(InpName), FORM='formatted',            &
     &        STATUS='old', ERR=10)
        GotFile=.TRUE.
        GO TO 30
  10    IF (localPET.eq.0) WRITE (stdout,20) TRIM(InpName)
        exit_flag=2
        RETURN
  20    FORMAT (/,' GETPAR_0D_S - Unable to input script: ',a)
      ELSE
        InpUnit=stdinp_unit(localPET, GotFile)
      END IF
!
!  Process requested parameter.
!
  30  foundit=.FALSE.
      DO WHILE (.TRUE.)
        READ (InpUnit,'(a)',ERR=40,END=60) line
        status=decode_line(line, string, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(string).eq.TRIM(KeyWord)) THEN
            DO i=1,LEN(Value)
              Value(i:i)=blank
            END DO
            Value=TRIM(ADJUSTL(Cval(Nval)))
            foundit=.TRUE.
          END IF
        END IF
      END DO
  40  IF (localPET.eq.0) THEN
        WRITE (stdout,50) line
  50    FORMAT (/,' GETPAR_0D_S - Error while processing line: ',/,a)
      END IF
      exit_flag=4
  60  CONTINUE
      IF (.not.foundit) THEN
        IF (localPET.eq.0) THEN
          WRITE (stdout,70) TRIM(KeyWord)
  70    FORMAT (/,' GETPAR_0D_S - unable to find KeyWord: ',a,          &
     &          /,15x,'in ROMS standard input file.')
        END IF
        exit_flag=5
      END IF
      IF (GotFile) THEN
        CLOSE (InpUnit)
      END IF
!
      RETURN
      END SUBROUTINE getpar_0d_s
!
      END MODULE stdinp_mod
