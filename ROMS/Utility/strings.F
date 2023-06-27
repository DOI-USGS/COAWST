#include "cppdefs.h"
      MODULE strings_mod
!
!git $Id$
!svn $Id: strings.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several string manipulation functions:         !
!                                                                      !
!    FoundError       checks action flag against no-error code         !
!    GlobalError      global check action flag against no-error code   !
!    StandardName     sets standard_name attribute                     !
!    TaskError        task check action flag against no-error code     !
!    assign_string    allocates and assign string                      !
!    find_string      scans a character array for a specified string   !
!    join_string      concatenate character array into a single string !
!    lowercase        converts input string characters to lowercase    !
!    uppercase        converts input string characters to uppercase    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    IF (.not.find_string(var_name,n_var,'spherical',varid)) THEN      !
!      ...                                                             !
!    END IF                                                            !
!                                                                      !
!    string=lowercase('MY UPPERCASE STRING')                           !
!                                                                      !
!    string=uppercase('my lowercase string')                           !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
!
      PUBLIC :: FoundError
      PUBLIC :: GlobalError
      PUBLIC :: StandardName
      PUBLIC :: TaskError
      PUBLIC :: assign_string
      PUBLIC :: find_string
      PUBLIC :: join_string
      PUBLIC :: lowercase
      PUBLIC :: uppercase
!
      CONTAINS
!
      FUNCTION FoundError (flag, NoErr, line, routine) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This logical function checks ROMS execution flag against no-error   !
!  code and issue a message if they are not equal.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     flag         ROMS execution flag (integer)                       !
!     NoErr        No Error code (integer)                             !
!     line         Calling model routine line (integer)                !
!     routine      Calling model routine (string)                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    execution flag is in error.                       !
!                                                                      !
!=======================================================================
!
      USE mod_iounits,  ONLY : stdout
      USE mod_parallel, ONLY : Master
!
!  Imported variable declarations.
!
      integer, intent(in) :: flag, NoErr, line

      character (len=*), intent(in) :: routine
!
!  Local variable declarations.
!
      logical :: foundit
!
!-----------------------------------------------------------------------
!  Scan array for requested string.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      IF (flag.ne.NoErr) THEN
        foundit=.TRUE.
        IF (Master) THEN
          WRITE (stdout,10) flag, line, TRIM(routine)
  10      FORMAT (' Found Error: ', i2.2, t20, 'Line: ', i0,            &
     &            t35, 'Source: ', a)
        END IF
        CALL my_flush (stdout)
      END IF

      RETURN
      END FUNCTION FoundError
!
      FUNCTION GlobalError (ng, model, flag, NoErr, line, routine)      &
     &                     RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This logical function checks ROMS execution flag against no-error   !
!  code and issue a message if they are not equal. It can be used in   !
!  split, disjointed, distributed-memory communicators.  All process   !
!  in application needs to call this function and the master process   !
!  knows to flag value to broadcast.                                   !
!                                                                      !
!  If not distributed-memory, this function has the same capability    !
!  as FoundError.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     flag         ROMS execution flag (integer)                       !
!     NoErr        No Error code (integer)                             !
!     line         Calling model routine line (integer)                !
!     routine      Calling model routine (string)                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     flag         Updated flag value to all processes (integer)       !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    execution flag is in error.                       !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_iounits,    ONLY : stdout

#if defined DISTRIBUTE && defined DISJOINTED
      USE distribute_mod, ONLY : mp_barrier, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, NoErr, line
      integer, intent(inout) :: flag

      character (len=*), intent(in) :: routine
!
!  Local variable declarations.
!
      logical :: foundit
      logical :: MasterProcess
#ifdef DISTRIBUTE
      integer :: MyCOMM
#endif
!
!-----------------------------------------------------------------------
!  Set master process.
#ifdef DISTRIBUTE
!  Set distribute-memory communicator and broadcast flag value to all
!  processes.
#endif
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
# ifdef DISJOINTED
      MyCOMM=FULL_COMM_WORLD
      MasterProcess=FullRank.eq.0
# else
      MyCOMM=OCN_COMM_WORLD
      MasterProcess=MyRank.eq.0
# endif
#else
      MasterProcess=Master
#endif
#if defined DISTRIBUTE && defined DISJOINTED
!
      CALL mp_bcasti (ng, model, flag, MyCOMM)
#endif
!
!-----------------------------------------------------------------------
!  Scan array for requested string.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      IF (flag.ne.NoErr) THEN
        foundit=.TRUE.
        IF (MasterProcess) THEN
          WRITE (stdout,10) flag, line, TRIM(routine)
  10      FORMAT (' Found Error: ', i2.2, t20, 'Line: ', i0,            &
     &            t35, 'Source: ', a)
        END IF
        CALL my_flush (stdout)
      END IF

      RETURN
      END FUNCTION GlobalError
!
      FUNCTION TaskError (ng, model, flag, NoErr, line, routine)        &
     &                    RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This logical function checks ROMS execution flag against no-error   !
!  code and issue a message if they are not equal. It can be used in   !
!  split, disjointed-tasks, distributed-memory communicators.  All     !
!  processes in the task section needs to call this function and its   !
!  master process knows to flag value to broadcast.                    !
!                                                                      !
!  If not distributed-memory, this function has the same capability    !
!  as FoundError and GlobalError.                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     flag         ROMS execution flag (integer)                       !
!     NoErr        No Error code (integer)                             !
!     line         Calling model routine line (integer)                !
!     routine      Calling model routine (string)                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     flag         Updated flag value to all processes (integer)       !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    execution flag is in error.                       !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_iounits,    ONLY : stdout

#if defined DISTRIBUTE && defined DISJOINTED
      USE distribute_mod, ONLY : mp_barrier, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, NoErr, line
      integer, intent(inout) :: flag

      character (len=*), intent(in) :: routine
!
!  Local variable declarations.
!
      logical :: foundit
      logical :: MasterProcess
#ifdef DISTRIBUTE
      integer :: MyCOMM
#endif
!
!-----------------------------------------------------------------------
!  Set master process.
#ifdef DISTRIBUTE
!  Set distribute-memory communicator and broadcast flag value to all
!  processes.
#endif
!-----------------------------------------------------------------------
!
#ifdef DISTRIBUTE
# ifdef DISJOINTED
#  ifdef CONCURRENT_KERNEL
      MyCOMM=TASK_COMM_WORLD
      MasterProcess=TaskRank.eq.0
#  else
      MyCOMM=FULL_COMM_WORLD
      MasterProcess=FullRank.eq.0
#  endif
# else
      MyCOMM=OCN_COMM_WORLD
      MasterProcess=MyRank.eq.0
# endif
#else
      MasterProcess=Master
#endif
#if defined DISTRIBUTE && defined DISJOINTED
!
      CALL mp_bcasti (ng, model, flag, MyCOMM)
#endif
!
!-----------------------------------------------------------------------
!  Scan array for requested string.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      IF (flag.ne.NoErr) THEN
        foundit=.TRUE.
        IF (MasterProcess) THEN
          WRITE (stdout,10) flag, line, TRIM(routine)
  10      FORMAT (' Found Error: ', i2.2, t20, 'Line: ', i0,            &
     &            t35, 'Source: ', a)
        END IF
        CALL my_flush (stdout)
      END IF

      RETURN
      END FUNCTION TaskError
!
      SUBROUTINE StandardName (Sname, variable, prefix, suffix)
!
!=======================================================================
!                                                                      !
!  This routine concatenates prefix and suffix strings to generate the !
!  'standard_name' attribute. Blank spaces in the prefix are replaced  !
!  with underscore.                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     variable     Standard name main variable (character)             !
!     prefix       Standard name prefix (OPTIONAL, character)          !
!     suffix       Standard name suffix (OPTIONAL, character)          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Sname        concatenated standard name (character)              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: variable
      character (len=*), optional, intent(in) :: prefix
      character (len=*), optional, intent(in) :: suffix
      character (len=*), intent(out) :: Sname
!
!  Local variable declarations.
!
      integer :: Icomma, Lstr, i
!
!-----------------------------------------------------------------------
!  Generate 'standard_name' attribute.
!-----------------------------------------------------------------------
!
      DO i=1,LEN(Sname)
        Sname(i:i)=CHAR(32)
      END DO
!
!  Replace blank space with underscore, CHAR(95).
!
      Lstr=LEN_TRIM(variable)
      Icomma=INDEX(variable,CHAR(44),BACK=.FALSE.)
      IF (Icomma.gt.0) THEN                               ! remove comma
        Sname=variable(1:Icomma-1) // variable(Icomma+1:Lstr)
        Lstr=LEN_TRIM(Sname)
      ELSE
        Sname(1:Lstr)=variable(1:Lstr)
      END IF
!
      DO i=1,Lstr
        IF (Sname(i:i).eq.CHAR(32)) THEN
          Sname(i:i)=CHAR(95)
        END IF
      END DO
!
!  Append prefix and or suffix string(s).

      IF (PRESENT(prefix)) THEN
        Sname=TRIM(ADJUSTL(prefix)) // TRIM(ADJUSTL(Sname))
      END IF
!
      IF (PRESENT(suffix)) THEN
        Sname=TRIM(ADJUSTL(Sname)) // TRIM(ADJUSTL(suffix))
      END IF
!
      RETURN
      END SUBROUTINE StandardName
!
      FUNCTION assign_string (A, string) RESULT (ErrFlag)
!
!=======================================================================
!                                                                      !
!  This routine assigns allocatable strings.  It allocates/reallocates !
!  output string variable.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     string       String to be assigned (character)                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A            Assigned allocatable string (character)             !
!     ErrFlag      Error flag (integer)                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=:), allocatable, intent(inout) :: A
      character (len=*), intent(in) :: string
!
!  Local variable declarations.
!
      integer :: lstr
      integer :: ErrFlag
!
!-----------------------------------------------------------------------
!  Allocate output string to the size of input string.
!-----------------------------------------------------------------------
!
      ErrFlag=-1
!
      lstr=LEN_TRIM(string)
      IF (.not.allocated(A)) THEN
        allocate ( character(LEN=lstr) :: A, STAT=ErrFlag)
      ELSE
        deallocate (A)
        allocate ( character(LEN=lstr) :: A, STAT=ErrFlag)
      END IF
!
!  Assign requested value.
!
      A=string
!
      RETURN
      END FUNCTION assign_string
!
      FUNCTION find_string (A, Asize, string, Aindex) RESULT (foundit)
!
!=======================================================================
!                                                                      !
!  This logical function scans an array of type character for an       !
!  specific string.                                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     A            Array of strings (character)                        !
!     Asize        Size of A (integer)                                 !
!     string       String to search (character)                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aindex       Array element containing the string (integer)       !
!     foundit      The value of the result is TRUE/FALSE if the        !
!                    string was found or not.                          !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Asize

      integer, intent(out) :: Aindex

      character (len=*), intent(in) :: A(Asize)
      character (len=*), intent(in) :: string
!
!  Local variable declarations.
!
      logical :: foundit

      integer :: i
!
!-----------------------------------------------------------------------
!  Scan array for requested string.
!-----------------------------------------------------------------------
!
      foundit=.FALSE.
      Aindex=0
      DO i=1,Asize
        IF (TRIM(A(i)).eq.TRIM(string)) THEN
          foundit=.TRUE.
          Aindex=i
          EXIT
        END IF
      END DO

      RETURN
      END FUNCTION find_string
!
      SUBROUTINE join_string (A, Asize, string, Lstring)
!
!=======================================================================
!                                                                      !
!  This routine concatenate a character array into a single string     !
!  with each element separated by commas.                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     A            Array of strings (character)                        !
!     Asize        Size of A (integer)                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     string       Concatenated string (character)                     !
!     Lstring      Length of concatenated string (integer)             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Asize

      integer, intent(out) :: Lstring

      character (len=*), intent(in) :: A(Asize)
      character (len=*), intent(out) :: string
!
!  Local variable declarations.
!
      integer :: i, ie, is, lstr
!
!-----------------------------------------------------------------------
!  Concatenate input character array.
!-----------------------------------------------------------------------
!
!  Initialize to blank string.
!
      lstr=LEN(string)
      DO i=1,lstr
        string(i:i)=' '
      END DO
!
!  Concatenate.
!
      is=1
      DO i=1,Asize
        lstr=LEN_TRIM(A(i))
        IF (lstr.gt.0) THEN
          ie=is+lstr-1
          string(is:ie)=TRIM(A(i))
          is=ie+1
          string(is:is)=','
          is=is+2
        END IF
      END DO
      Lstring=LEN_TRIM(string)-1

      RETURN
      END SUBROUTINE join_string
!
      FUNCTION lowercase (Sinp) RESULT (Sout)
!
!=======================================================================
!                                                                      !
!  This character function converts input string elements to           !
!  lowercase.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Sinp       String with uppercase elements (character)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Sout       Lowercase string (character)                          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Cooper Redwine, 1995: "Upgrading to Fortran 90", Springer-       !
!       Verlag, New York, pp 416.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (*), intent(in) :: Sinp
!
!  Local variable definitions.
!
      integer :: i, j, lstr

      character (LEN(Sinp)) :: Sout

      character (26), parameter :: Lcase = 'abcdefghijklmnopqrstuvwxyz'
      character (26), parameter :: Ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
!-----------------------------------------------------------------------
!  Convert input string to lowercase.
!-----------------------------------------------------------------------
!
      lstr=LEN(Sinp)
      Sout=Sinp
      DO i=1,lstr
        j=INDEX(Ucase, Sout(i:i))
        IF (j.ne.0) THEN
          Sout(i:i)=Lcase(j:j)
        END IF
      END DO

      RETURN
      END FUNCTION lowercase
!
      FUNCTION uppercase (Sinp) RESULT (Sout)
!
!=======================================================================
!                                                                      !
!  This character function converts input string elements to           !
!  uppercase.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Sinp       String with lowercase characters (character)          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Sout       Uppercase string (character)                          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Cooper Redwine, 1995: "Upgrading to Fortran 90", Springer-       !
!       Verlag, New York, pp 416.                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (*), intent(in) :: Sinp
!
!  Local variable definitions.
!
      integer :: i, j, lstr

      character (LEN(Sinp)) :: Sout

      character (26), parameter :: Lcase = 'abcdefghijklmnopqrstuvwxyz'
      character (26), parameter :: Ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
!-----------------------------------------------------------------------
!  Convert input string to uppercase.
!-----------------------------------------------------------------------
!
      lstr=LEN(Sinp)
      Sout=Sinp
      DO i=1,lstr
        j=INDEX(Lcase, Sout(i:i))
        IF (j.ne.0) THEN
          Sout(i:i)=Ucase(j:j)
        END IF
      END DO

      RETURN
      END FUNCTION uppercase

      END MODULE strings_mod
