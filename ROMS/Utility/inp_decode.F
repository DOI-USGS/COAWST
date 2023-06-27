#include "cppdefs.h"
      MODULE inp_decode_mod
!
!git $Id$
!svn $Id: inp_decode.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several routines to process and decode ROMS    !
!  unique namelist KeyWord parameters from input script files:         !
!                                                                      !
!  decode_line    Decodes line of text from input script files for a   !
!                   particular KeyWord.                                !
!                                                                      !
!  find_file      Checks if provide input filename exits.              !
!                                                                      !
!  load_i         Processes and loads an integer parameter variable.   !
!                                                                      !
!  load_i         Processes and loads a logical parameter variable.    !
!                                                                      !
!  load_r         Processes and loads a single or double precision     !
!                   floating-point (real) parameter variable.          !
!                                                                      !
!  load_lbc       Processes and loads lateral boundary conditions      !
!                   logical switches into derived type structure,      !
!                   TYPE(T_LBC).                                       !
!                                                                      !
!  load_s1d       Processes and loads I/O parameters into 1D derived   !
!                   type structure, TYPE(T_IO).                        !
!                                                                      !
!  load_s2d       Processes and loads I/O parameters into 2D derived   !
!                   type structure, TYPE(T_IO).                        !
!                                                                      !
#ifdef SOLVE3D
!  load_tadv      Processes and loads tracers horizontal and vertical  !
!                   advection switches into derived type structure,    !
!                   TYPE(T_ADV).                                       !
#endif
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
      INTERFACE load_i
        MODULE PROCEDURE load_0d_i       ! scalar integer
        MODULE PROCEDURE load_1d_i       ! 1D integer array
        MODULE PROCEDURE load_2d_i       ! 2D integer array
        MODULE PROCEDURE load_3d_i       ! 3D integer array
      END INTERFACE load_i

      INTERFACE load_l
        MODULE PROCEDURE load_0d_l       ! scalar logical
        MODULE PROCEDURE load_1d_l       ! 1D logical array
        MODULE PROCEDURE load_2d_l       ! 2D logical array
        MODULE PROCEDURE load_3d_l       ! 3D logical array
      END INTERFACE load_l

      INTERFACE load_r
#ifdef SINGLE_PRECISION
        MODULE PROCEDURE load_0d_dp      ! scalar real(dp)
        MODULE PROCEDURE load_1d_dp      ! 1D real(dp) array
        MODULE PROCEDURE load_2d_dp      ! 2D real(dp) array
        MODULE PROCEDURE load_3d_dp      ! 3D real(dp) array
#endif
        MODULE PROCEDURE load_0d_r8      ! scalar real(r8)
        MODULE PROCEDURE load_1d_r8      ! 1D real(r8) array
        MODULE PROCEDURE load_2d_r8      ! 2D real(r8) array
        MODULE PROCEDURE load_3d_r8      ! 3D real(r8) array
      END INTERFACE load_r

      INTERFACE load_s1d
        MODULE PROCEDURE load_s1d1       ! 1D structrure, S(:)
        MODULE PROCEDURE load_s1d2       ! 2D structrure, S(Ie,:)
      END INTERFACE load_s1d
!
!  Module dimension parameters.
!
      integer, parameter :: nCval = 200  ! size of Cval character vector
      integer, parameter :: nRval = 100  ! size of Rval real vector
!
      CONTAINS
!
      FUNCTION decode_line (line_text, KeyWord, Nval, Cval, Rval)
!
!***********************************************************************
!                                                                      !
!  This function decodes lines of text from input script files. It is  !
!  to evaluate ROMS unique namelist parameters.                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!    line_text      Input scripts lines as text (string)               !
!    KeyWord        Input parameter keyword (string)                   !
!    Nval           Number of values processed (integer)               !
!    Cval           Input values as characters (string array)          !
!    Rval           Input values as mumbers (real array)               !
!                                                                      !
!***********************************************************************
!
! Imported variable declarations.
!
      character (len=*), intent(in) :: line_text

      character (len=*), intent(inout) :: KeyWord

      integer, intent(inout) :: Nval

      character (len=*), intent(inout) :: Cval(:)

      real(dp), intent(inout) :: Rval(:)
!
! Local variable declarations
!
      logical :: IsString, Kextract, decode, nested
      integer :: Iblank, Icomm, Icont, Ipipe, Kstr, Kend, Linp
      integer :: Lend, LenS, Lstr, Lval, Nmul, Schar
      integer :: copies, i, ic, ie, is, j, status

      integer, dimension(20) :: Imul

      integer :: decode_line

      character (len=256) :: Vstring, inpline, line, string
!
!------------------------------------------------------------------------
!  Decode input line.
!------------------------------------------------------------------------
!
!  Initialize. Use CHAR(32) for blank space.
!
      line=CHAR(32)
      inpline=CHAR(32)
      Vstring=CHAR(32)
      string=CHAR(32)
!!    Lstr=LEN(line)
!!    DO i=1,Lstr
!!      line(i:i)=CHAR(32)
!!      inpline(i:i)=CHAR(32)
!!      Vstring(i:i)=CHAR(32)
!!      string(i:i)=CHAR(32)
!!    END DO
!
!  Check input line and remove illegal characters.  Replace control
!  ASCII characters CHAR(0) to CHAR(31) with a blank space, CHAR(32).
!
!  Char  Dec  Key  Control Action
!  ----------------------------------------------------------------------
!  NUL   0    ^@   Null character
!  SOH   1    ^A   Start of heading, = console interrupt
!  STX   2    ^B   Start of text, maintenance mode on HP console
!  ETX   3    ^C   End of text
!  EOT   4    ^D   End of transmission, not the same as ETB
!  ENQ   5    ^E   Enquiry, goes with ACK; old HP flow control
!  ACK   6    ^F   Acknowledge, clears ENQ logon hand
!  BEL   7    ^G   Bell, rings the bell...
!  BS    8    ^H   Backspace, works on HP terminals/computers
!  HT    9    ^I   Horizontal tab, move to next tab stop
!  LF    10   ^J   Line Feed
!  VT    11   ^K   Vertical tab
!  FF    12   ^L   Form Feed, page eject
!  CR    13   ^M   Carriage Return
!  SO    14   ^N   Shift Out, alternate character set
!  SI    15   ^O   Shift In, resume defaultn character set
!  DLE   16   ^P   Data link escape
!  DC1   17   ^Q   XON, with XOFF to pause listings; ":okay to send".
!  DC2   18   ^R   Device control 2, block-mode flow control
!  DC3   19   ^S   XOFF, with XON is TERM=18 flow control
!  DC4   20   ^T   Device control 4
!  NAK   21   ^U   Negative acknowledge
!  SYN   22   ^V   Synchronous idle
!  ETB   23   ^W   End transmission block, not the same as EOT
!  CAN   24   ^X   Cancel line, MPE echoes !!!
!  EM    25   ^Y   End of medium, Control-Y interrupt
!  SUB   26   ^Z   Substitute
!  ESC   27   ^[   Escape, next character is not echoed
!  FS    28   ^\   File separator
!  GS    29   ^]   Group separator
!  RS    30   ^^   Record separator, block-mode terminator
!  US    31   ^_   Unit separator
!
!  SP    32        Space
!  !     33        Exclamation mark
!  #     35        Number sign, hash, or pound sign
!  *     42        Asterisk (star, multiply)
!  +     43        Plus
!  -     45        Hyphen, dash, minus
!  .     46        Period
!  0     48        Zero
!  1     49        One
!  2     50        Two
!  3     51        Three
!  4     52        Four
!  5     53        Five
!  6     54        Six
!  7     55        Seven
!  8     56        Eight
!  9     57        Nine
!  :     58        colon sign
!  =     61        Equals sign
!  \     92        Reverse slant (Backslash)
!  |    124        Vertical line
!
      inpline=TRIM(ADJUSTL(line_text))
      Linp=LEN_TRIM(inpline)
      DO i=1,LEN_TRIM(inpline)
        j=ICHAR(inpline(i:i))
        IF (j.lt.32) THEN
          inpline(i:i)=char(32)                           ! blank space
        END IF
      END DO
      inpline=TRIM(inpline)
!
!  Get length of "line". Remove comment after the KEYWORD, if any.
!  Then, remove leading and trailing blanks.
!
      IF ((Linp.gt.0).and.(inpline(1:1).ne.CHAR(33))) THEN
        Icomm=INDEX(inpline,CHAR(33),BACK=.FALSE.)
        IF (Icomm.gt.0) Linp=Icomm-1
        line=TRIM(ADJUSTL(inpline(1:Linp)))
        Linp=LEN_TRIM(line)
      ELSE
        line=TRIM(ADJUSTL(inpline))
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
!
      RETURN
      END FUNCTION decode_line
!
      FUNCTION find_file (ng, out, fname, KeyWord) RESULT (foundit)
!
!***********************************************************************
!                                                                      !
!  This function checks if provided input file exits.                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number                                    !
!     out        Standard output unit                                  !
!     fname      Filename (path and name)                              !
!     KeyWord    Keyword associated with file name (string,OPTIONAL)   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     foundit    The value of the result is TRUE/FALSE if the file     !
!                  was found or not                                    !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, out
!
      character (len=*), intent(in) :: fname
      character (len=*), intent(in) :: KeyWord
!
!  Local variable declarations.
!
      logical :: foundit, isURL
!
      integer :: lstr, ncid
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", find_file"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Check if the file exit.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
!
!  Check for empty file name string.
!
      lstr=LEN_TRIM(fname)
      IF (lstr.eq.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(KeyWord)
 10       FORMAT (/,' FIND_FILE - empty file name string ',             &
     &            'for standard input script KeyWord: ',a)
        END IF
        exit_flag=5
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
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
          CALL netcdf_close (ng, iNLM, ncid, fname, .FALSE.)
        END IF
      END IF
!
!  Report if not found.
!
      IF (.not.foundit) THEN
        IF (Master) WRITE (out,20) ng, TRIM(fname)
 20     FORMAT (/,' FIND_FILE - Grid ',i2.2,                           &
     &          ', cannot find input file: ',a)
        exit_flag=5
      END IF
!
      RETURN
      END FUNCTION find_file
!
      FUNCTION load_0d_i (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model scalar integer         !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Nout       Size of output integer variable dimension (not used)  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output scalar integer variable                        !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Nout
      real(dp), intent(in) :: Vinp(:)
!
      integer, intent(out) :: Vout
!
!  Local variable declarations.
!
      integer :: ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load scalar integer variable with input value.
!-----------------------------------------------------------------------
!
      ic=1
      Vout=INT(Vinp(ic))
      Nval=ic

      RETURN
      END FUNCTION load_0d_i
!
      FUNCTION load_1d_i (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 1D integer array.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Nout       Size of output integer variable dimension             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 1D integer variable                            !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Nout
      real(dp), intent(in) :: Vinp(:)
!
      integer, intent(out) :: Vout(:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load 1D integer variable with input values.
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
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vout(i)=INT(Vinp(Ninp))
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
      END IF
      Nval=ic

      RETURN
      END FUNCTION load_1d_i
!
      FUNCTION load_2d_i (Ninp, Vinp, Iout, Jout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 2D integer array.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Iout       Size of output integer variable first  I-dimension    !
!     Jout       Size of output integer variable second J-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 2D integer variable                            !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Iout, Jout
      real(dp), intent(in) :: Vinp(:)
!
      integer, intent(out) :: Vout(:,:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
      integer, dimension(Iout*Jout) :: Vwrk
!
!-----------------------------------------------------------------------
!  Load 2D integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the 2D array.
!
      ic=0
      Nout=Iout*Jout
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vwrk(i)=INT(Vinp(i))
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vwrk(i)=INT(Vinp(Ninp))
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vwrk(i)=INT(Vinp(i))
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_2d_i
!
      FUNCTION load_3d_i (Ninp, Vinp, Iout, Jout, Kout, Vout)           &
     &                    RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 3D integer array.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Iout       Size of output integer variable first  I-dimension    !
!     Jout       Size of output integer variable second J-dimension    !
!     Kout       Size of output integer variable third  K-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 3D integer variable                            !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Iout, Jout, Kout
      real(dp), intent(in) :: Vinp(:)
!
      integer, intent(out) :: Vout(:,:,:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
      integer, dimension(Iout*Jout*Kout) :: Vwrk
!
!-----------------------------------------------------------------------
!  Load 3D integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the 3D array.
!
      ic=0
      Nout=Iout*Jout*Kout
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vwrk(i)=INT(Vinp(i))
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vwrk(i)=INT(Vinp(Ninp))
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vwrk(i)=INT(Vinp(i))
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout,Kout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_3d_i
!
      FUNCTION load_0d_l (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model scalar logical         !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (character 1D array)                     !
!     Nout       Size of output logical variable dimension (not used)  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output scalar logical variable                        !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      character (len=*), intent(in) :: Vinp(:)
!
      logical, intent(out) :: Vout
!
!  Local variable declarations.
!
      integer :: ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load scalar logical variable with input value.
!-----------------------------------------------------------------------
!
      ic=1
      IF ((Vinp(ic)(1:1).eq.'T').or.                                    &
     &    (Vinp(ic)(1:1).eq.'t')) THEN
        Vout=.TRUE.
      ELSE
        Vout=.FALSE.
      END IF
      Nval=ic
!
      RETURN
      END FUNCTION load_0d_l
!
      FUNCTION load_1d_l (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 1D logical array.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (character 1D array)                     !
!     Nout       Size of output logical variable dimension             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 1D logical variable                            !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      character (len=*), intent(in) :: Vinp(:)
!
      logical, intent(out) :: Vout(:)
!
!  Local variable declarations.
!
      logical :: LastValue

      integer :: Nstr, i, ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load logical variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      LastValue=.FALSE.
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.                                 &
     &        (Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
          LastValue=Vout(i)
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vout(i)=LastValue
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.                                 &
     &        (Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.TRUE.
          ELSE
            Vout(i)=.FALSE.
          END IF
        END DO
      END IF
      Nval=ic
!
      RETURN
      END FUNCTION load_1d_l
!
      FUNCTION load_2d_l (Ninp, Vinp, Iout, Jout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 2D logical array.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (character 1D array)                     !
!     Iout       Size of output logical variable first  I-dimension    !
!     Jout       Size of output logical variable second J-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 2D logical variable                            !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Iout, Jout
      character (len=*), intent(in) :: Vinp(:)
!
      logical, intent(out) :: Vout(:,:)
!
!  Local variable declarations.
!
      logical :: LastValue
!
      logical, dimension(Iout*Jout) :: Vwrk
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
!-----------------------------------------------------------------------
!  Load 2D logical variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      Nout=Iout*Jout
      LastValue=.FALSE.
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.                                 &
     &        (Vinp(i)(1:1).eq.'t')) THEN
            Vwrk(i)=.TRUE.
          ELSE
            Vwrk(i)=.FALSE.
          END IF
          LastValue=Vwrk(i)
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vwrk(i)=LastValue
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.                                 &
     &        (Vinp(i)(1:1).eq.'t')) THEN
            Vwrk(i)=.TRUE.
          ELSE
            Vwrk(i)=.FALSE.
          END IF
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_2d_l
!
      FUNCTION load_3d_l (Ninp, Vinp, Iout, Jout, Kout, Vout)           &
     &                   RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 3D logical array.      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (character 1D array)                     !
!     Iout       Size of output logical variable first  I-dimension    !
!     Jout       Size of output logical variable second J-dimension    !
!     Kout       Size of output logical variable third  K-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 3D logical variable                            !
!     Nval       Number of output values processed                     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Iout, Jout, Kout
      character (len=*), intent(in) :: Vinp(:)
!
      logical, intent(out) :: Vout(:,:,:)
!
!  Local variable declarations.
!
      logical :: LastValue
!
      logical, dimension(Iout*Jout*Kout) :: Vwrk
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
!-----------------------------------------------------------------------
!  Load 3D logical variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      Nout=Iout*Jout*Kout
      LastValue=.FALSE.
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.                                 &
     &        (Vinp(i)(1:1).eq.'t')) THEN
            Vwrk(i)=.TRUE.
          ELSE
            Vwrk(i)=.FALSE.
          END IF
          LastValue=Vwrk(i)
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vwrk(i)=LastValue
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.                                 &
     &        (Vinp(i)(1:1).eq.'t')) THEN
            Vwrk(i)=.TRUE.
          ELSE
            Vwrk(i)=.FALSE.
          END IF
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout,Kout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_3d_l

#ifdef SINGLE_PRECISION
!
      FUNCTION load_0d_dp (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model scalar double          !
!  precision variable when numerical kernel is in single precision.    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Nout       Size of output integer variable dimension (not used)  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output scalar variable (real, KIND=dp)                !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Nout
      real(dp), intent(in) :: Vinp(:)
!
      real(dp), intent(out) :: Vout
!
!  Local variable declarations.
!
      integer :: ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load scalar floating-point variable with input value.
!-----------------------------------------------------------------------
!
      ic=1
      Vout=Vinp(ic)
      Nval=ic
!
      RETURN
      END FUNCTION load_0d_dp
!
      FUNCTION load_1d_dp (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 1D double precision    !
!  array when numerical kernel is in single precision.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Nout       Size of output integer variable dimension             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 1D variable (real, KIND=dp)                    !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Nout
      real(dp), intent(in) :: Vinp(:)
!
      real(dp), intent(out) :: Vout(:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load 1D floating-point variable with input values.
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
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vout(i)=Vinp(Ninp)
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
      END IF
      Nval=ic
!
      RETURN
      END FUNCTION load_1d_dp
!
      FUNCTION load_2d_dp (Ninp, Vinp, Iout, Jout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 2D double precision    !
!  array when numerical kernel is in single precision.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Iout       Size of output integer variable first  I-dimension    !
!     Jout       Size of output integer variable second J-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 2D variable (real, KIND=dp)                    !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Iout, Jout
      real(dp), intent(in) :: Vinp(:)
!
      real(dp), intent(out) :: Vout(:,:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
      real(dp), dimension(Iout*Jout) :: Vwrk
!
!-----------------------------------------------------------------------
!  Load 2D floating-point variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      Nout=Iout*Jout
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vwrk(i)=Vinp(i)
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vwrk(i)=Vinp(Ninp)
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vwrk(i)=Vinp(i)
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_2d_dp
!
      FUNCTION load_3d_dp (Ninp, Vinp, Iout, Jout, Kout, Vout)          &
     &                     RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 3D double precision    !
!  array when numerical kernel is in single precision.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Iout       Size of output integer variable first  I-dimension    !
!     Jout       Size of output integer variable second J-dimension    !
!     Kout       Size of output integer variable third  K-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 3D variable (real, KIND=dp)                    !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Iout, Jout, Kout
      real(dp), intent(in) :: Vinp(:)
!
      real(dp), intent(out) :: Vout(:,:,:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
      real(dp), dimension(Iout*Jout*Kout) :: Vwrk
!
!-----------------------------------------------------------------------
!  Load 3D floating-point variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      Nout=Iout*Jout*Kout
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vwrk(i)=Vinp(i)
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
            Vwrk(i)=Vinp(Ninp)
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vwrk(i)=Vinp(i)
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout,Kout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_3d_dp
#endif
!
      FUNCTION load_0d_r8 (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!=======================================================================
!                                                                      !
!  It loads input values into a requested model scalar floating-point  !
!  variable (KIND=r8).                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Nout       Size of output integer variable dimension (not used)  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output scalar variable (real, KIND=r8)                !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Nout
      real(dp), intent(in) :: Vinp(:)
!
      real(r8), intent(out) :: Vout
!
!  Local variable declarations.
!
      integer :: ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load scalar floating-point variable with input value.
!-----------------------------------------------------------------------
!
      ic=1
#ifdef SINGLE_PRECISION
      Vout=REAL(Vinp(ic),r8)
#else
      Vout=Vinp(ic)
#endif
      Nval=ic
!
      RETURN
      END FUNCTION load_0d_r8
!
      FUNCTION load_1d_r8 (Ninp, Vinp, Nout, Vout) RESULT (Nval)
!
!=======================================================================
!                                                                      !
!  It loads input values into a requested model 1D floating-point      !
!  array (KIND=r8).                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Nout       Size of output integer variable dimension             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 1D variable (real, KIND=r8)                    !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Nout
      real(dp), intent(in) :: Vinp(:)
!
      real(r8), intent(out) :: Vout(:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nval
!
!-----------------------------------------------------------------------
!  Load 1D floating-point variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
#ifdef SINGLE_PRECISION
          Vout(i)=REAL(Vinp(i),r8)
#else
          Vout(i)=Vinp(i)
#endif
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
#ifdef SINGLE_PRECISION
            Vout(i)=REAL(Vinp(Ninp),r8)
#else
            Vout(i)=Vinp(Ninp)
#endif
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
#ifdef SINGLE_PRECISION
          Vout(i)=REAL(Vinp(i),r8)
#else
          Vout(i)=Vinp(i)
#endif
        END DO
      END IF
      Nval=ic
!
      RETURN
      END FUNCTION load_1d_r8
!
      FUNCTION load_2d_r8 (Ninp, Vinp, Iout, Jout, Vout) RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 2D floating-point      !
!  array (KIND=r8).                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Iout       Size of output integer variable first  I-dimension    !
!     Jout       Size of output integer variable second J-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 2D variable (real, KIND=r8)                    !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Iout, Jout
      real(dp), intent(in) :: Vinp(:)
!
      real(r8), intent(out) :: Vout(:,:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
      real(r8), dimension(Iout*Jout) :: Vwrk
!
!-----------------------------------------------------------------------
!  Load 2D floating-point variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      Nout=Iout*Jout
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
#ifdef SINGLE_PRECISION
          Vwrk(i)=REAL(Vinp(i),r8)
#else
          Vwrk(i)=Vinp(i)
#endif
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
#ifdef SINGLE_PRECISION
            Vwrk(i)=REAL(Vinp(Ninp),r8)
#else
            Vwrk(i)=Vinp(Ninp)
#endif
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
#ifdef SINGLE_PRECISION
          Vwrk(i)=REAL(Vinp(i),r8)
#else
          Vwrk(i)=Vinp(i)
#endif
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_2d_r8
!
      FUNCTION load_3d_r8 (Ninp, Vinp, Iout, Jout, Kout, Vout)          &
     &                     RESULT (Nval)
!
!***********************************************************************
!                                                                      !
!  It loads input values into a requested model 3D floating-point      !
!  array (KIND=r8).                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Number of input elements to process in Vinp (integer) !
!     Vinp       Input values (1D real(dp) array)                      !
!     Iout       Size of output integer variable first  I-dimension    !
!     Jout       Size of output integer variable second J-dimension    !
!     Kout       Size of output integer variable third  K-dimension    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output 3D variable (real, KIND=r8)                    !
!     Nval       Number of output values processed                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,  intent(in) :: Ninp, Iout, Jout, Kout
      real(dp), intent(in) :: Vinp(:)
!
      real(r8), intent(out) :: Vout(:,:,:)
!
!  Local variable declarations.
!
      integer :: Nstr, i, ic
      integer :: Nout, Nval
!
      real(r8), dimension(Iout*Jout*Kout) :: Vwrk
!
!-----------------------------------------------------------------------
!  Load 3D floating-point variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      Nout=Iout*Jout*Kout
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
#ifdef SINGLE_PRECISION
          Vwrk(i)=REAL(Vinp(i),r8)
#else
          Vwrk(i)=Vinp(i)
#endif
        END DO
        IF (Nout.gt.Ninp) THEN
          Nstr=Ninp+1
          DO i=Nstr,Nout
            ic=ic+1
#ifdef SINGLE_PRECISION
            Vwrk(i)=REAL(Vinp(Ninp),r8)
#else
            Vwrk(i)=Vinp(Ninp)
#endif
          END DO
        END IF
      ELSE
        DO i=1,Nout
          ic=ic+1
#ifdef SINGLE_PRECISION
          Vwrk(i)=REAL(Vinp(i),r8)
#else
          Vwrk(i)=Vinp(i)
#endif
        END DO
      END IF
      Vout=RESHAPE(Vwrk,(/Iout,Jout,Kout/))
      Nval=ic
!
      RETURN
      END FUNCTION load_3d_r8
!
      FUNCTION load_lbc (Ninp, Vinp, line, nline, ifield, igrid,        &
     &                   iTrcStr, iTrcEnd, svname, S)
!
!***********************************************************************
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
!***********************************************************************
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
      integer :: Icont, i, ib, ic
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
        DO ib=1,4
          string=uppercase(Bstring(ib))
          SELECT CASE (TRIM(string))
            CASE ('CHA')
              S(ib,ifield,igrid)%Chapman_implicit = .TRUE.
            CASE ('CHE')
              S(ib,ifield,igrid)%Chapman_explicit = .TRUE.
            CASE ('CLA')
              S(ib,ifield,igrid)%clamped = .TRUE.
              S(ib,ifield,igrid)%acquire = .TRUE.
            CASE ('CLO')
              S(ib,ifield,igrid)%closed = .TRUE.
            CASE ('FLA')
              S(ib,ifield,igrid)%Flather = .TRUE.
              S(ib,ifield,igrid)%acquire = .TRUE.
              S(ib,isFsur,igrid)%acquire = .TRUE.
            CASE ('GRA')
              S(ib,ifield,igrid)%gradient = .TRUE.
            CASE ('NES')
              S(ib,ifield,igrid)%nested = .TRUE.
            CASE ('PER')
              S(ib,ifield,igrid)%periodic = .TRUE.
              IF ((ib.eq.ieast).or.(ib.eq.iwest)) THEN
                EWperiodic(igrid)=.TRUE.
              ELSE IF ((ib.eq.inorth).or.(ib.eq.isouth)) THEN
                NSperiodic(igrid)=.TRUE.
              END IF
            CASE ('RAD')
              S(ib,ifield,igrid)%radiation = .TRUE.
            CASE ('RADNUD')
              S(ib,ifield,igrid)%radiation = .TRUE.
              S(ib,ifield,igrid)%nudging = .TRUE.
              S(ib,ifield,igrid)%acquire = .TRUE.
            CASE ('RED')
              S(ib,ifield,igrid)%reduced = .TRUE.
#if defined FSOBC_REDUCED
              S(ib,isFsur,igrid)%acquire = .TRUE.
#endif
            CASE ('SHC')
              S(ib,ifield,igrid)%Shchepetkin = .TRUE.
              S(ib,ifield,igrid)%acquire = .TRUE.
              S(ib,isFsur,igrid)%acquire = .TRUE.
            CASE DEFAULT
              IF (Master) THEN
                WRITE (stdout,10) TRIM(Vinp(ib)), TRIM(line)
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
              DO ib=1,4
                S(ib,i,igrid)%clamped   = S(ib,ifield,igrid)%clamped
                S(ib,i,igrid)%closed    = S(ib,ifield,igrid)%closed
                S(ib,i,igrid)%gradient  = S(ib,ifield,igrid)%gradient
                S(ib,i,igrid)%nested    = S(ib,ifield,igrid)%nested
                S(ib,i,igrid)%periodic  = S(ib,ifield,igrid)%periodic
                S(ib,i,igrid)%radiation = S(ib,ifield,igrid)%radiation
                S(ib,i,igrid)%nudging   = S(ib,ifield,igrid)%nudging
                S(ib,i,igrid)%acquire   = S(ib,ifield,igrid)%acquire
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

 10   FORMAT (/,' LOAD_LBC - illegal lateral boundary condition ',      &
     &        'keyword: ',a,/,12x,a)
 20   FORMAT (/,' LOAD_LBC - incorrect continuation symbol in line:',/, &
     &        12x,a,/,12x,'number of nested grid values exceeded.')
!
      RETURN
      END FUNCTION load_lbc
!
      FUNCTION load_s1d1 (Nval, Fname, Fdim, line, label, igrid,        &
     &                    Mgrids, Nfiles, io_type, S)
!
!***********************************************************************
!                                                                      !
!  This function loads input values into requested 1D structure        !
!  S(Mgrids) containing information about I/O files.                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Nval        Number of values processed (integer)                 !
!     Fname       File name(s) processed (string array)                !
!     Fdim        File name(s) dimension in calling program (integer)  !
!     line        Current input line (string)                          !
!     label       I/O structure label (string)                         !
!     igrid       Nested grid counter (integer)                        !
!     Mgrids      Number of nested grids (integer)                     !
!     Nfiles      Number of files per grid (integer array)             !
!     io_type     File I/O type (integer)                              !
!     S(Mgrids)   Derived type structure array, TYPE(T_IO)             !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     igrid       Updated nested grid counter.                         !
!     S(Mgrids)   Updated derived type structure array, TYPE(T_IO).    !
!     load_s1d_1  Number of output values processed.                   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in)    :: Mgrids, Nval, Fdim, io_type
      integer, intent(inout) :: igrid
      integer, intent(inout) :: Nfiles(Mgrids)

      character (len=*),   intent(in) :: line
      character (len=256), intent(in) :: Fname(Fdim)
      character (len=*),   intent(inout) :: label

      TYPE(T_IO), intent(inout) :: S(Mgrids)
!
!  Local variable declarations.
!
      logical :: load, persist

      integer :: Icont, Ipipe, i, is, j, lstr, my_Mgrids, ng
      integer :: load_s1d1

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
      IF (igrid.gt.Mgrids) THEN
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
!  If nesting and the number of file name entries is less than Mgrids,
!  persist the last values provided.  This is the case when not enough
!  entries are provided by "==" plural symbol after the KEYWORD.
!
        IF (igrid.lt.Mgrids) THEN
          DO i=igrid+1,Mgrids
            Nfiles(i)=Nfiles(igrid)
          END DO
          my_Mgrids=igrid
          persist=.TRUE.
        ELSE
          my_Mgrids=Mgrids
          persist=.FALSE.
        END IF
!
!  Allocate various fields in structure, if not continuation or pipe
!  symbol is found which indicates end of input data.
!
        IF (label(1:3).eq.'FLT') THEN
#ifdef FLOAT_BIOLOGY
          is=-10
#else
          is=-6
#endif
        ELSE
          is=1
        END IF
!
        DO ng=1,Mgrids
          allocate ( S(ng)%Nrec(Nfiles(ng)) )
          allocate ( S(ng)%time_min(Nfiles(ng)) )
          allocate ( S(ng)%time_max(Nfiles(ng)) )
          allocate ( S(ng)%Vid(is:NV) )
          allocate ( S(ng)%Tid(MT) )
#if defined PIO_LIB && defined DISTRIBUTE
          allocate ( S(ng)%pioVar(is:NV) )
          allocate ( S(ng)%pioTrc(MT) )
#endif
          allocate ( S(ng)%files(Nfiles(ng)) )
        END DO
!
!  Intialize strings to blank to facilitate processing.
!
        DO ng=1,Mgrids
          lstr=LEN(S(ng)%name)
          DO i=1,lstr
            S(ng)%head(i:i)=blank
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
        DO ng=1,my_Mgrids
          S(ng)%IOtype=io_type                 ! file IO type
          S(ng)%Nfiles=Nfiles(ng)              ! number of multi-files
          S(ng)%Fcount=1                       ! multi-file counter
          S(ng)%load=1                         ! filename load counter
          S(ng)%Rindex=0                       ! time index
          S(ng)%ncid=-1                        ! closed NetCDF state
          S(ng)%Vid=-1                         ! NetCDF variables IDs
          S(ng)%Tid=-1                         ! NetCDF tracers IDs
#if defined PIO_LIB && defined DISTRIBUTE
          S(ng)%pioFile%fh=-1                  ! closed file handler
          DO j=1,NV
            S(ng)%pioVar(j)%vd%varID=-1        ! variables IDs
            S(ng)%pioVar(j)%dkind=-1           ! variables data kind
            S(ng)%pioVar(j)%gtype=0            ! variables C-grid type
          END DO
          DO j=1,MT
            S(ng)%pioTrc(j)%vd%varID=-1        ! tracers IDs
            S(ng)%pioTrc(j)%dkind=-1           ! tracers data kind
            S(ng)%pioTrc(j)%gtype=0            ! tracers C-grid type
          END DO
#endif
          DO j=1,Nfiles(ng)
            i=i+1
            S(ng)%files(j)=TRIM(Fname(i))      ! load multi-files
            S(ng)%Nrec(j)=0                    ! record counter
            S(ng)%time_min(j)=0.0_dp           ! starting time
            S(ng)%time_max(j)=0.0_dp           ! ending time
          END DO
          S(ng)%label=TRIM(label)              ! structure label
          S(ng)%name=TRIM(S(ng)%files(1))      ! load first file
          lstr=LEN_TRIM(S(ng)%name)
          S(ng)%head=S(ng)%name(1:lstr-3)      ! do not include ".nc"
          S(ng)%base=S(ng)%name(1:lstr-3)      ! do not include ".nc"
          Nfiles(ng)=0                         ! clean file counter
        END DO
!
!  If appropriate, persist last value(s).
!
        IF (persist) THEN
          DO ng=igrid+1,Mgrids
            S(ng)%IOtype=io_type
            S(ng)%Nfiles=S(igrid)%Nfiles
            S(ng)%Fcount=1
            S(ng)%load=1
            S(ng)%Rindex=0
            S(ng)%ncid=-1
            S(ng)%Vid=-1
            S(ng)%Tid=-1
#if defined PIO_LIB && defined DISTRIBUTE
            S(ng)%pioFile%fh=-1
            DO j=is,NV
              S(ng)%pioVar(j)%vd%varID=-1
              S(ng)%pioVar(j)%dkind=-1
              S(ng)%pioVar(j)%gtype=0
            END DO
            DO j=1,MT
              S(ng)%pioTrc(j)%vd%varID=-1
              S(ng)%pioTrc(j)%dkind=-1
              S(ng)%pioTrc(j)%gtype=0
            END DO
#endif
            DO j=1,S(igrid)%Nfiles
              S(ng)%files(j)=S(igrid)%files(j)
              S(ng)%Nrec(j)=0
              S(ng)%time_min(j)=0.0_dp
              S(ng)%time_max(j)=0.0_dp
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
        DO ng=1,Mgrids
          Nfiles(ng)=0
        END DO
        DO i=1,LEN(label)
          label(i:i)=blank
        END DO
      END IF
      load_s1d1=Nval

 10   FORMAT (/,' LOAD_S1D1 - incorrect continuation symbol in line:',  &
     &        /,14x,a,/,11x,'number of nested grid values exceeded.')
!
      RETURN
      END FUNCTION load_s1d1
!
      FUNCTION load_s1d2 (Nval, Fname, Fdim, line, label, igrid,        &
     &                    Mgrids, Nfiles, idim, Ie, io_type, S)
!
!***********************************************************************
!                                                                      !
!  This function loads input values into requested 2D structure        !
!  S(Ie,:) elemement containing information about I/O files.           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Nval         Number of values processed (integer)                !
!     Fname        File name(s) processed (string array)               !
!     Fdim         File name(s) dimension in calling program (integer) !
!     line         Current input line (string)                         !
!     label        I/O structure label (string)                        !
!     igrid        Nested grid counter (integer)                       !
!     Mgrids       Number of nested grids (integer)                    !
!     Nfiles       Number of files per grid (integer array)            !
!     idim         Size of structure inner dimension (integer)         !
!     Ie           Inner dimension element to process (integer)        !
!     io_type      File I/O type (integer)                             !
!     S(Ie,Mgrids) Derived type structure array, TYPE(T_IO)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     igrid        Updated nested grid counter.                        !
!     S(Ie,Mgrids) Updated derived type structure array, TYPE(T_IO).   !
!     load_s1d_2   Number of output values processed.                  !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in)    :: Mgrids, Nval, Fdim, Ie, idim, io_type
      integer, intent(inout) :: igrid
      integer, intent(inout) :: Nfiles(Mgrids)

      character (len=*),   intent(in) :: line
      character (len=256), intent(in) :: Fname(Fdim)
      character (len=*),   intent(inout) :: label

      TYPE(T_IO), intent(inout) :: S(idim,Mgrids)
!
!  Local variable declarations.
!
      logical :: load, persist

      integer :: Icont, Ipipe, i, is, j, lstr, my_Mgrids, ng
      integer :: load_s1d2

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
      IF (igrid.gt.Mgrids) THEN
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
!  If nesting and the number of file name entries is less than Mgrids,
!  persist the last values provided.  This is the case when not enough
!  entries are provided by "==" plural symbol after the KEYWORD.
!
        IF (igrid.lt.Mgrids) THEN
          DO i=igrid+1,Mgrids
            Nfiles(i)=Nfiles(igrid)
          END DO
          my_Mgrids=igrid
          persist=.TRUE.
        ELSE
          my_Mgrids=Mgrids
          persist=.FALSE.
        END IF
!
!  Allocate various fields in structure, if not continuation or pipe
!  symbol is found which indicates end of input data.
!
        DO ng=1,Mgrids
          allocate ( S(Ie,ng)%Nrec(Nfiles(ng)) )
          allocate ( S(Ie,ng)%time_min(Nfiles(ng)) )
          allocate ( S(Ie,ng)%time_max(Nfiles(ng)) )
          allocate ( S(Ie,ng)%Vid(NV) )
          allocate ( S(Ie,ng)%Tid(MT) )
#if defined PIO_LIB && defined DISTRIBUTE
          allocate ( S(Ie,ng)%pioVar(NV) )
          allocate ( S(Ie,ng)%pioTrc(MT) )
#endif
          allocate ( S(Ie,ng)%files(Nfiles(ng)) )
        END DO
!
!  Intialize strings to blank to facilitate processing.
!
        DO ng=1,Mgrids
          lstr=LEN(S(Ie,ng)%name)
          DO i=1,lstr
            S(Ie,ng)%head(i:i)=blank
            S(Ie,ng)%base(i:i)=blank
            S(Ie,ng)%name(i:i)=blank
          END DO
          DO j=1,Nfiles(ng)
            DO i=1,lstr
              S(Ie,ng)%files(j)(i:i)=blank
            END DO
          END DO
        END DO
!
!  Initialize and load fields into structure.
!
        i=0
        DO ng=1,my_Mgrids
          S(Ie,ng)%IOtype=io_type               ! file IO type
          S(Ie,ng)%Nfiles=Nfiles(ng)            ! number of multi-files
          S(Ie,ng)%Fcount=1                     ! multi-file counter
          S(Ie,ng)%load=1                       ! filename load counter
          S(Ie,ng)%Rindex=0                     ! time index
          S(Ie,ng)%ncid=-1                      ! closed NetCDF state
          S(Ie,ng)%Vid=-1                       ! NetCDF variables IDs
          S(Ie,ng)%Tid=-1                       ! NetCDF tracers IDs
#if defined PIO_LIB && defined DISTRIBUTE
          S(Ie,ng)%pioFile%fh=-1                ! closed file handler
          DO j=1,NV
            S(Ie,ng)%pioVar(j)%vd%varID=-1      ! variables IDs
            S(Ie,ng)%pioVar(j)%dkind=-1         ! variables data kind
            S(Ie,ng)%pioVar(j)%gtype=0          ! variables C-grid type
          END DO
          DO j=1,MT
            S(Ie,ng)%pioTrc(j)%vd%varID=-1      ! tracers IDs
            S(Ie,ng)%pioTrc(j)%dkind=-1         ! tracers data kind
            S(Ie,ng)%pioTrc(j)%gtype=0          ! tracers C-grid type
          END DO
#endif
          DO j=1,Nfiles(ng)
            i=i+1
            S(Ie,ng)%files(j)=TRIM(Fname(i))    ! load multi-files
            S(Ie,ng)%Nrec(j)=0                  ! record counter
            S(Ie,ng)%time_min(j)=0.0_dp         ! starting time
            S(Ie,ng)%time_max(j)=0.0_dp         ! ending time
          END DO
          S(Ie,ng)%label=TRIM(label)            ! structure label
          S(Ie,ng)%name=TRIM(S(Ie,ng)%files(1)) ! load first file
          lstr=LEN_TRIM(S(Ie,ng)%name)
          S(Ie,ng)%head=S(Ie,ng)%name(1:lstr-3) ! do not include ".nc"
          S(Ie,ng)%base=S(Ie,ng)%name(1:lstr-3) ! do not include ".nc"
          Nfiles(ng)=0                          ! clean file counter
        END DO
!
!  If appropriate, persist last value(s).
!
        IF (persist) THEN
          DO ng=igrid+1,Mgrids
            S(Ie,ng)%IOtype=io_type
            S(Ie,ng)%Nfiles=S(Ie,igrid)%Nfiles
            S(Ie,ng)%Fcount=1
            S(Ie,ng)%load=1
            S(Ie,ng)%Rindex=0
            S(Ie,ng)%ncid=-1
            S(Ie,ng)%Vid=-1
            S(Ie,ng)%Tid=-1
#if defined PIO_LIB && defined DISTRIBUTE
            S(Ie,ng)%pioFile%fh=-1
            DO j=1,NV
              S(Ie,ng)%pioVar(j)%vd%varID=-1
              S(Ie,ng)%pioVar(j)%dkind=-1
              S(Ie,ng)%pioVar(j)%gtype=0
            END DO
            DO j=1,MT
              S(Ie,ng)%pioTrc(j)%vd%varID=-1
              S(Ie,ng)%pioTrc(j)%dkind=-1
              S(Ie,ng)%pioTrc(j)%gtype=0
            END DO
#endif
            DO j=1,S(Ie,igrid)%Nfiles
              S(Ie,ng)%files(j)=S(Ie,igrid)%files(j)
              S(Ie,ng)%Nrec(j)=0
              S(Ie,ng)%time_min(j)=0.0_dp
              S(Ie,ng)%time_max(j)=0.0_dp
            END DO
            S(Ie,ng)%label=TRIM(label)
            S(Ie,ng)%name=S(Ie,igrid)%name
            S(Ie,ng)%base=S(Ie,igrid)%base
            Nfiles(ng)=0
          END DO
        END IF
!
!  Reset counters and clean label.
!
        igrid=1
        DO ng=1,Mgrids
          Nfiles(ng)=0
        END DO
        DO i=1,LEN(label)
          label(i:i)=blank
        END DO
      END IF
      load_s1d2=Nval

 10   FORMAT (/,' LOAD_S1D2 - incorrect continuation symbol in line:',  &
     &        /,14x,a,/,11x,'number of nested grid values exceeded.')
!
      RETURN
      END FUNCTION load_s1d2
!
      FUNCTION load_s2d (Nval, Fname, Fdim, line, label, ifile, igrid,  &
     &                   Mgrids, Nfiles, Ncount, idim, io_type, S)
!
!***********************************************************************
!                                                                      !
!  This function loads input values into requested 2D structure        !
!  containing information about input forcing files. Notice that       !
!  Mgrids is passed for flexibility in coupling algorithms.            !
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
!     Mgrids     Number of nested grids (integer)                      !
!     Nfiles     Number of input files per grid (integer vector)       !
!     Ncount     Number of files per grid counter (integer array)      !
!     idim       Size of structure inner dimension (integer)           !
!     io_type    File I/O type (integer)                               !
!     S          Derived type structure, TYPE(T_IO)                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ifile      Updated file counter.                                 !
!     igrid      Updated nested grid counter.                          !
!     S          Updated derived type structure, TYPE(T_IO).           !
!     load_s2d   Number of output values processed.                    !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in)    :: Mgrids, Nval, Fdim, idim, io_type
      integer, intent(in)    :: Nfiles(Mgrids)
      integer, intent(inout) :: ifile, igrid
      integer, intent(inout) :: Ncount(idim,Mgrids)
!
      character (len=*),   intent(in) :: line
      character (len=256), intent(in) :: Fname(Fdim)
      character (len=*),   intent(inout) :: label
!
      TYPE(T_IO), intent(inout) :: S(idim,Mgrids)
!
!  Local variable declarations.
!
      logical :: load, persist
!
      integer :: Icont, Ipipe, i, is, j, k, lstr, my_Mgrids, ng
      integer :: load_s2d
!
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
      IF (igrid.gt.Mgrids) THEN
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
!  If nesting and the number of file name entries is less than Mgrids,
!  persist the last values provided.  This is the case when not enough
!  entries are provided by "==" plural symbol after the KEYWORD.
!
        IF (igrid.lt.Mgrids) THEN
          DO j=igrid+1,Mgrids
            DO i=1,idim
              Ncount(i,j)=Ncount(i,igrid)
            END DO
          END DO
          my_Mgrids=igrid
          persist=.TRUE.
        ELSE
          my_Mgrids=Mgrids
          persist=.FALSE.
        END IF
!
!  Allocate various fields in structure, if not continuation or pipe
!  symbol is found which indicates end of input data.
!
        DO ng=1,Mgrids
          DO i=1,idim
            allocate ( S(i,ng)%Nrec(Ncount(i,ng)) )
            allocate ( S(i,ng)%time_min(Ncount(i,ng)) )
            allocate ( S(i,ng)%time_max(Ncount(i,ng)) )
            allocate ( S(i,ng)%Vid(NV) )
            allocate ( S(i,ng)%Tid(MT) )
#if defined PIO_LIB && defined DISTRIBUTE
            allocate ( S(i,ng)%pioVar(NV) )
            allocate ( S(i,ng)%pioTrc(MT) )
#endif
            allocate ( S(i,ng)%files(Ncount(i,ng)) )
          END DO
        END DO
!
!  Intialize strings to blank to facilitate processing.
!
        DO ng=1,Mgrids
          DO i=1,idim
            lstr=LEN(S(i,ng)%name)
            DO j=1,lstr
              S(i,ng)%head(j:j)=blank
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
        DO ng=1,my_Mgrids
          DO i=1,Nfiles(ng)
            S(i,ng)%IOtype=io_type              ! file IO type
            S(i,ng)%Nfiles=Ncount(i,ng)         ! number of multi-files
            S(i,ng)%Fcount=1                    ! multi-file counter
            S(i,ng)%load=1                      ! filename load counter
            S(i,ng)%Rindex=0                    ! time index
            S(i,ng)%ncid=-1                     ! closed NetCDF state
            S(i,ng)%Vid=-1                      ! NetCDF variables IDs
            S(i,ng)%Tid=-1                      ! NetCDF tracers IDs
#if defined PIO_LIB && defined DISTRIBUTE
            S(i,ng)%pioFile%fh=-1               ! closed file handler
            DO j=1,NV
              S(i,ng)%pioVar(j)%vd%varID=-1     ! variables IDs
              S(i,ng)%pioVar(j)%dkind=-1        ! variables data kind
              S(i,ng)%pioVar(j)%gtype=0         ! variables C-grid type
            END DO
            DO j=1,MT
              S(i,ng)%pioTrc(j)%vd%varID=-1     ! tracers IDs
              S(i,ng)%pioTrc(j)%dkind=-1        ! tracers data kind
              S(i,ng)%pioTrc(j)%gtype=0         ! tracers C-grid type
            END DO
#endif
            DO j=1,Ncount(i,ng)
              k=k+1
              S(i,ng)%files(j)=TRIM(Fname(k))   ! load multi-files
              S(i,ng)%Nrec(j)=0                 ! record counter
              S(i,ng)%time_min(j)=0.0_dp        ! starting time
              S(i,ng)%time_max(j)=0.0_dp        ! ending time
            END DO
            S(i,ng)%label=TRIM(label)           ! structure label
            S(i,ng)%name=TRIM(S(i,ng)%files(1)) ! load first file
            lstr=LEN_TRIM(S(i,ng)%name)
            S(i,ng)%head=S(i,ng)%name(1:lstr-3) ! do not include ".nc"
            S(i,ng)%base=S(i,ng)%name(1:lstr-3) ! do not include ".nc"
          END DO
        END DO
!
!  If appropriate, persist last value(s).
!
        IF (persist) THEN
          DO ng=igrid+1,Mgrids
            DO i=1,Nfiles(ng)
              S(i,ng)%IOtype=io_type
              S(i,ng)%Nfiles=S(i,igrid)%Nfiles
              S(i,ng)%Fcount=1
              S(i,ng)%load=1
              S(i,ng)%Rindex=0
              S(i,ng)%ncid=-1
              S(i,ng)%Vid=-1
              S(i,ng)%Tid=-1
#if defined PIO_LIB && defined DISTRIBUTE
              S(i,ng)%pioFile%fh=-1
              DO j=1,NV
                S(i,ng)%pioVar(j)%vd%varID=-1
                S(i,ng)%pioVar(j)%dkind=-1
                S(i,ng)%pioVar(j)%gtype=0
              END DO
              DO j=1,MT
                S(i,ng)%pioTrc(j)%vd%varID=-1
                S(i,ng)%pioTrc(j)%dkind=-1
                S(i,ng)%pioTrc(j)%gtype=0
              END DO
#endif
              DO j=1,S(i,igrid)%Nfiles
                S(i,ng)%files(j)=S(i,igrid)%files(j)
                S(i,ng)%Nrec(j)=0
                S(i,ng)%time_min(j)=0.0_dp
                S(i,ng)%time_max(j)=0.0_dp
              END DO
              S(i,ng)%label=TRIM(label)
              S(i,ng)%head=S(i,igrid)%head
              S(i,ng)%base=S(i,igrid)%base
              S(i,ng)%name=S(i,igrid)%name
              Ncount(i,ng)=0
            END DO
          END DO
        END IF
!
!  Reset counters and clean label.
!
        igrid=1
        ifile=1
        DO ng=1,Mgrids
          DO i=1,idim
            Ncount(i,ng)=0
          END DO
        END DO
        DO i=1,LEN(label)
          label(i:i)=blank
        END DO
      END IF
      load_s2d=Nval

 10   FORMAT (/,' LOAD_S2D - incorrect continuation symbol in line:',/, &
     &        12x,a,/,12x,'inner dimension of structure exceeded.')
 20   FORMAT (/,' LOAD_S2D - incorrect continuation symbol in line:',/, &
     &        12x,a,/,12x,'number of nested grid values exceeded.')
!
      RETURN
      END FUNCTION load_s2d
!
#ifdef SOLVE3D
      FUNCTION load_tadv (Ninp, Vinp, line, nline, itrc, igrid,         &
     &                    itracer, iTrcStr, iTrcEnd, svname, S)
!
!***********************************************************************
!                                                                      !
!  This function sets tracers advection logical switches according to  !
!  input string keywords.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable (integer)                      !
!     Vinp       Input values (string)                                 !
!     line       Current input line (string)                           !
!     nline      Multi-line counter (integer)                          !
!     itrc       Tracer array index (integer)                          !
!     itracer    Calling routine tracer counter (integer)              !
!     igrid      Nested grid counter (integer)                         !
!     iTrcStr    Starting tracer index to process (integer)            !
!     iTrcEnd    Ending   tracer index to process (integer)            !
!     svname     State variable name (string)                          !
!     S          Derived type structure, TYPE(T_ADV)                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nline      Updated multi-line counter (integer)                  !
!     itracer    Updated calling routine tracer counter (integer)      !
!     igrid      Updated nested grid counter (integer)                 !
!     S          Updated derived type structure, TYPE(T_LBC)           !
!     load_tadv  Number of output values processed.                    !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, itrc, iTrcStr, iTrcEnd
      integer, intent(inout) :: igrid, itracer, nline

      character (len=256), intent(in) :: line
      character (len=256), intent(in) :: Vinp(Ninp)
      character (len=*  ), intent(in) :: svname

      TYPE(T_ADV), intent(inout) :: S(MAXVAL(NT),Ngrids)
!
!  Local variable declarations.
!
      integer :: Icont, i, ic
      integer :: load_tadv

      character (len=10) :: Astring, string
!
!-----------------------------------------------------------------------
!  Set tracers advection switches in structure.
!-----------------------------------------------------------------------
!
!  Check current line for the continuation symbol [char(92)=\].
!
      Icont=INDEX(TRIM(line),CHAR(92) ,BACK=.FALSE.)
!
!  Extract tracer advection scheme keywords from Vinp. Notice that
!  additional array elements are added to Vinp during continuation
!  lines.
!
      i=nline
      Astring=TRIM(Vinp(i+1))
!
!  Advance or reset entry lines counter.
!
      IF (Icont.gt.0) THEN
        nline=nline+1
      ELSE
        nline=0
      END IF
!
!  Set advection switches for each tracer.
!
      ic=1
      IF ((0.lt.itrc).and.(itrc.le.iTrcEnd)) THEN
        string=uppercase(Astring)
        SELECT CASE (TRIM(string))
          CASE ('A4', 'AKIMA4')
            S(itrc,igrid) % AKIMA4 = .TRUE.
          CASE ('C2', 'CENTERED2')
            S(itrc,igrid) % CENTERED2 = .TRUE.
          CASE ('C4', 'CENTERED4')
            S(itrc,igrid) % CENTERED4 = .TRUE.
          CASE ('HS', 'HSIMT')
            S(itrc,igrid) % HSIMT = .TRUE.
          CASE ('MP', 'MPDATA')
            S(itrc,igrid) % MPDATA = .TRUE.
          CASE ('SP', 'SPLINES')
            S(itrc,igrid) % SPLINES = .TRUE.
          CASE ('SU', 'SU3', 'SPLIT_U3')
            S(itrc,igrid) % SPLIT_U3 = .TRUE.
          CASE ('U3', 'UPSTREAM3')
            S(itrc,igrid) % UPSTREAM3 = .TRUE.
          CASE DEFAULT
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Astring)
            END IF
            exit_flag=2
            RETURN
        END SELECT
!
!  If processing tracers and last standard input entry (Icont=0), set
!  unspecified tracer values to the last tracer entry.
!
        IF ((iTrcStr.gt.0).and.(iTrcEnd.gt.0)) THEN
          IF ((Icont.eq.0).and.(itracer.lt.iTrcEnd)) THEN
            DO i=itrc+1,iTrcEnd
              S(i,igrid) % AKIMA4    = S(itrc,igrid) % AKIMA4
              S(i,igrid) % CENTERED2 = S(itrc,igrid) % CENTERED2
              S(i,igrid) % CENTERED4 = S(itrc,igrid) % CENTERED4
              S(i,igrid) % HSIMT     = S(itrc,igrid) % HSIMT
              S(i,igrid) % MPDATA    = S(itrc,igrid) % MPDATA
              S(i,igrid) % SPLINES   = S(itrc,igrid) % SPLINES
              S(i,igrid) % SPLIT_U3  = S(itrc,igrid) % SPLIT_U3
              S(i,igrid) % UPSTREAM3 = S(itrc,igrid) % UPSTREAM3
            END DO
            ic=ic+1
          END IF
        END IF
      END IF
!
!  If appropriate, reset tracer grid counter.  It is done to process
!  other keywords using this function.
!
      IF ((itrc.eq.iTrcEnd).or.(ic.gt.1)) THEN
        itracer=0
      END IF
!
!  If appropriate, increase or reset nested grid counter.
!
      IF ((Icont.gt.0).and.(Ngrids.gt.1)) THEN
        IF ((iTrcStr.gt.0).and.(iTrcEnd.gt.0)) THEN
          IF ((itrc.eq.iTrcEnd).or.(ic.gt.1)) THEN
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
      load_tadv=ic

 10   FORMAT (/,' LOAD_TADV - illegal tracer advection scheme ',        &
     &        'keyword: ',a,/,13x,'Correct standard input file.',/)
 20   FORMAT (/,' LOAD_TADV - incorrect continuation symbol in line:',  &
     &        /,13x,a,/,13x,'number of nested grid values exceeded.')
!
      RETURN
      END FUNCTION load_tadv
#endif
!
      END MODULE inp_decode_mod
