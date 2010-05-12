#include "cppdefs.h"
      MODULE strings_mod
!
!svn $Id: strings.F 429 2009-12-20 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module contains several string manipulation functions:         !
!                                                                      !
!    find_string   scans a character array for a specified string      !
!    lowercase     converts input string characters to lowercase       !
!    uppercase     converts input string characters to uppercase       !
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

      PUBLIC :: find_string
      PUBLIC :: lowercase
      PUBLIC :: uppercase

      CONTAINS

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

      FUNCTION lowercase (Ainp) RESULT (Aout)
!
!=======================================================================
!                                                                      !
!  This character function converts input string elements to           !
!  lowercase.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ainp       String with uppercase elements (character)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aout       Lowercase string (character)                          !
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
      character (*), intent(in) :: Ainp
!
!  Local variable definitions.
!
      integer :: i, j, lstr

      character (LEN(Ainp)) :: Aout

      character (26), parameter :: Lcase = 'abcdefghijklmnopqrstuvwxyz'
      character (26), parameter :: Ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
!-----------------------------------------------------------------------
!  Convert input string to lowercase.
!-----------------------------------------------------------------------
!
      lstr=LEN(Ainp)
      Aout=Ainp
      DO i=1,lstr
        j=INDEX(Ucase, Aout(i:i))
        IF (j.ne.0) THEN
          Aout(i:i)=Lcase(j:j)
        END IF
      END DO

      RETURN
      END FUNCTION lowercase

      FUNCTION uppercase (Ainp) RESULT (Aout)
!
!=======================================================================
!                                                                      !
!  This character function converts input string elements to           !
!  uppercase.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ainp       String with lowercase characters (character)          !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aout       Uppercase string (character)                          !
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
      character (*), intent(in) :: Ainp
!
!  Local variable definitions.
!
      integer :: i, j, lstr

      character (LEN(Ainp)) :: Aout

      character (26), parameter :: Lcase = 'abcdefghijklmnopqrstuvwxyz'
      character (26), parameter :: Ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!
!-----------------------------------------------------------------------
!  Convert input string to uppercase.
!-----------------------------------------------------------------------
!
      lstr=LEN(Ainp)
      Aout=Ainp
      DO i=1,lstr
        j=INDEX(Lcase, Aout(i:i))
        IF (j.ne.0) THEN
          Aout(i:i)=Ucase(j:j)
        END IF
      END DO

      RETURN
      END FUNCTION uppercase

      END MODULE strings_mod
