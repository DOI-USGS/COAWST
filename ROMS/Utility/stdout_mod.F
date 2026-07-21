#include "cppdefs.h"
      MODULE stdout_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  It sets the ROMS standard output unit to write verbose execution    !
!  information. Notice that the default standard out unit in Fortran   !
!  is 6.                                                               !
!                                                                      !
!  In some applications like coupling or disjointed mpi-communications,!
!  it is advantageous to write standard output to a specific filename  !
!  instead of the default Fortran standard output unit 6. If that is   !
!  the case, it opens such formatted file for writing.                 !
!                                                                      !
!=======================================================================
!
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod,    ONLY : FoundError
!
      implicit none

      PUBLIC  :: stdout_unit
!
!-----------------------------------------------------------------------
!  Module parameters.
!-----------------------------------------------------------------------
!
!  The following switch tells if the standard output unit/file has
!  been specified. It must be set up only once to avoid errors and
!  is usually called at the beginning of ROMS_initialize. However,
!  it is required earlier during ESM coupling configurations.
!
      logical, save :: Set_StdOutUnit = .TRUE.
!
!-----------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------
!
      FUNCTION stdout_unit (MyMaster) RESULT (StdOutUnit)
!
!***********************************************************************
!                                                                      !
!  This function determines ROMS standard output unit to write its     !
!  running verbose information.                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     MyMaster    Switch indicating Master process (logical)           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     StdOutUnit  Assigned standard output unit (integer; default=6)   !
!                                                                      !
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, intent(in)  :: MyMaster
!
!  Local variable declararions
!
      integer :: io_err
      integer :: StdOutUnit
!
      character (len=10 )          :: stdout_file
      character (len=256)          :: io_errmsg
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", stdout_unit"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set ROMS standard input unit. If requested, set and open standard
!  output file for writing.
!-----------------------------------------------------------------------

#ifdef ROMS_STDOUT
!
!  Set the Fortran standard output unit (default unit 6) to direct the
!  ROMS run information to a specific filename. It is advantageous in
!  coupling or disjointed applications that need ROMS information
!  separated from other components.
!
# if defined DISTRIBUTE && defined DISJOINTED
      StdOutUnit=60+ForkColor
      WRITE (stdout_file,'(a,i2.2,a)') 'log', ForkColor+1, '.roms'
# else
      StdOutUnit=60
      stdout_file='log.roms'
# endif
#else
!
! Set default standard output unit in Fortran.
!
      StdOutUnit=6
#endif

#ifdef ROMS_STDOUT
!
!-----------------------------------------------------------------------
!  Open ROMS standard output file.
!-----------------------------------------------------------------------
!
      io_err=0
      IF (Lappend) THEN
        OPEN (StdOutUnit, FILE=TRIM(stdout_file),                       &
     &        FORM='formatted', STATUS='old', POSITION='append',        &
     &        ACTION='write', IOSTAT=io_err, IOMSG=io_errmsg)
      ELSE
        OPEN (StdOutUnit, FILE=TRIM(stdout_file),                       &
     &        FORM='formatted', STATUS='replace',                       &
     &        IOSTAT=io_err, IOMSG=io_errmsg)
      END IF
      IF (io_err.ne.0) THEN
        IF (MyMaster) WRITE (StdOutUnit,10) TRIM(stdout_file),          &
     &                                      TRIM(io_errmsg)
        exit_flag=5
 10     FORMAT (/,' STDOUT_UNIT - Cannot open standard output file: ',  &
     &          a,/,15x,'ERROR: ',a)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
#endif
!
      END FUNCTION stdout_unit
!
      END MODULE stdout_mod
