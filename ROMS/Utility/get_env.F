      MODULE get_env_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module retrieves and decodes system environmental variables    !
!  using Fortran 2003 GET_ENVIRONMENT_VARIABLE intrinsic function.     !
!                                                                      !
!=======================================================================
!
      USE mod_parallel, ONLY : Master
      USE mod_iounits,  ONLY : stdout
!
      implicit none
!
! Define public overloading function.
!
      INTERFACE get_env
        MODULE PROCEDURE get_env_i                ! gets integer value
        MODULE PROCEDURE get_env_l                ! gets logical value
        MODULE PROCEDURE get_env_s                ! gets string  value
      END INTERFACE get_env
!
      CONTAINS
!
      FUNCTION get_env_i (name, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  Reads and decodes environmental variable with an integer value.     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      integer,           intent(out) :: value
      character (len=*), intent(in ) :: name
!
!  Local variable declaration.
!
      integer                        :: Lstr, status
      character (len=40)             :: string
      character (len=512)            :: msg
!
!-----------------------------------------------------------------------
!  Read and decode environmental variable. Return "value=-1" if not
!  found.
!-----------------------------------------------------------------------
!
      CALL GET_ENVIRONMENT_VARIABLE (name, string, Lstr, status)
!
!  Convert to integer.
!
      IF ((Lstr.gt.0).and.(status.eq.0)) THEN
        READ (string, *, IOSTAT=status, IOMSG=msg) value
        IF (Master.and.(status.ne.0)) THEN
          WRITE (stdout,10) TRIM(name), TRIM(string), TRIM(msg)
        END IF
      ELSE
        status=1
        value=-1
        IF (Master) THEN
          WRITE (stdout,20) TRIM(name)
        END IF
      END IF
!
  10  FORMAT (/,' GET_ENV_I - Error while converting string to',        &
     &        ' integer, name = ',a,', value = ',a,/,13x,'ErrMsg: ',a)
  20  FORMAT (/,' GET_ENV_I - Cannot find environmental variable',      &
     &        ', name = ',a)
!
      RETURN
      END FUNCTION get_env_i
!
      FUNCTION get_env_l (name, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  Reads and decodes environmental variable with an logical value.     !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical,           intent(out) :: value
      character (len=*), intent(in ) :: name
!
!  Local variable declaration.
!
      integer                        :: Lstr, status
      character (len=40)             :: string
!
!-----------------------------------------------------------------------
!  Read and decode environmental variable. Return value=.FALSE. if not
!  found.
!-----------------------------------------------------------------------
!
      CALL GET_ENVIRONMENT_VARIABLE (name, string, Lstr, status)
!
!  Convert to logical.
!
      value=.FALSE.
!
      IF ((Lstr.gt.0).and.(status.eq.0)) THEN
        IF ((string(1:1).eq.'0') .or.                                   &
     &      (string(1:1).eq.'F') .or.                                   &
     &      (string(1:1).eq.'f')) THEN
          value=.FALSE.
        ELSE
          value=.TRUE.
        END IF
      END IF
!
      RETURN
      END FUNCTION get_env_l
!
      FUNCTION get_env_s (name, value) RESULT (status)
!
!***********************************************************************
!                                                                      !
!  Reads string environmental variable.                                !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      character (len=*), intent(in)  :: name
      character (len=*), intent(out) :: value
!
!  Local variable declaration.
!
      integer                        :: Lstr1, Lstr2, status
      character (len=1024)           :: string
!
!-----------------------------------------------------------------------
!  Read environmental variable. Return a blank space, value=CHAR(32),
!  if not found.
!-----------------------------------------------------------------------
!
      CALL GET_ENVIRONMENT_VARIABLE (name, string, Lstr1, status)
!
!  Load string value.
!
      Lstr2=LEN(value)
!
      IF (Lstr1.gt.Lstr2) THEN
        value=string(1:Lstr2)
        IF (Master) THEN
          WRITE (stdout,10) TRIM(name), TRIM(string), Lstr2, Lstr1
        END IF
      ELSE
        value=string(1:Lstr1)
      END IF
!
  10  FORMAT (/,' GET_ENV_S - Error while retrieving enviromental ',    &
     &        'variable, name = ',a,/,13x,"string = '",a,"'",/,13x,     &
     &        'value variable length = ',i0,                            &
     &        ' is less than the required length = ',i0)
!
      RETURN
      END FUNCTION get_env_s
!
      END MODULE get_env_mod
