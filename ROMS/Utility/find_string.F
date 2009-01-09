#include "cppdefs.h"
      FUNCTION find_string (A, Asize, string, Aindex) RESULT (foundit)
!
!svn $Id: find_string.F 195 2008-07-18 21:16:50Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This function scans an array of type character for an specific      !
!  string.                                                             !
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
