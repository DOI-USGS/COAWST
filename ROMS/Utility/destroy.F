#include "cppdefs.h"
      MODULE destroy_mod
!
!git $Id$
!svn $Id: destroy.F 1151 2023-02-09 03:08:53Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  It releases the space allocated for pointer variable in ROMS        !
!  kernel structures. After a variable has been deallocated, it        !
!  cannot be defined or referenced until it is allocated or            !
!  assigned again.                                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     Varray     Pointer variable to deallocate (real)                 !
!     routine    Calling routine (string)                              !
!     line       Calling routine line (integer)                        !
!     Vstring    Variable name (string)                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Lsuccess   Deallocation error switch (logical)                   !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_parallel, ONLY : Master
      USE mod_iounits,  ONLY : stdout
!
      implicit none
!
      INTERFACE destroy
        MODULE PROCEDURE destroy_1d_i    ! 1D integer array
        MODULE PROCEDURE destroy_1d_l    ! 1D logical array
#ifdef SINGLE_PRECISION
        MODULE PROCEDURE destroy_1d_dp   ! 1D real(dp) array
        MODULE PROCEDURE destroy_2d_dp   ! 2D real(dp) array
        MODULE PROCEDURE destroy_3d_dp   ! 3D real(dp) array
        MODULE PROCEDURE destroy_4d_dp   ! 4D real(dp) array
        MODULE PROCEDURE destroy_5d_dp   ! 5D real(dp) array
#endif
        MODULE PROCEDURE destroy_1d_r8   ! 1D real(r8) array
        MODULE PROCEDURE destroy_2d_r8   ! 2D real(r8) array
        MODULE PROCEDURE destroy_3d_r8   ! 3D real(r8) array
        MODULE PROCEDURE destroy_4d_r8   ! 4D real(r8) array
        MODULE PROCEDURE destroy_5d_r8   ! 5D real(r8) array
      END INTERFACE destroy
!
      integer, parameter :: Avar = 1
      integer, parameter :: Pvar = 2
!
      CONTAINS
!
!***********************************************************************
      FUNCTION destroy_1d_i (ng, Varray, routine, line, Vstring)        &
     &               RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line

      integer, pointer, intent(inout) :: Varray(:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 1D integer array.
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_1D_I - Grid ',i2.2,                           &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,16x,a)
!
      RETURN
      END FUNCTION destroy_1d_i
!
!***********************************************************************
      FUNCTION destroy_1d_l (ng, Varray, routine, line, Vstring)        &
     &               RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      logical, pointer, intent(inout) :: Varray(:)
!
      integer, intent(in) :: ng, line
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 1D logical array.
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_1D_L - Grid ',i2.2,                           &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,16x,a)
!
      RETURN
      END FUNCTION destroy_1d_l
!
!***********************************************************************
      FUNCTION destroy_1d_r8 (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(r8), pointer, intent(inout) :: Varray(:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 1D floating-point array (KIND=r8).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_1D_R8 - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_1d_r8
!
!***********************************************************************
      FUNCTION destroy_2d_r8 (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(r8), pointer, intent(inout) :: Varray(:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 2D floating-point array (KIND=r8).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_2D_R8 - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_2d_r8
!
!***********************************************************************
      FUNCTION destroy_3d_r8 (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(r8), pointer, intent(inout) :: Varray(:,:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 3D floating-point array (KIND=r8).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_3D_R8 - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_3d_r8
!
!***********************************************************************
      FUNCTION destroy_4d_r8 (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(r8), pointer, intent(inout) :: Varray(:,:,:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 4D floating-point array (KIND=r8).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_4D_R8 - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_4d_r8
!
!***********************************************************************
      FUNCTION destroy_5d_r8 (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(r8), pointer, intent(inout) :: Varray(:,:,:,:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 5D floating-point array (KIND=r8).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
#ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
#else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
#endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_5D_R8 - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_5d_r8

#ifdef SINGLE_PRECISION
!
!***********************************************************************
      FUNCTION destroy_1d_dp (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(dp), pointer, intent(inout) :: Varray(:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 1D double precision array (KIND=dp).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
# ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
# else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
# endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_1D_DP - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_1d_dp
!
!***********************************************************************
      FUNCTION destroy_2d_dp (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(dp), pointer, intent(inout) :: Varray(:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 2D double precision array (KIND=dp).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
# ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
# else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
# endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_2D_DP - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_2d_dp
!
!***********************************************************************
      FUNCTION destroy_3d_dp (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(dp), pointer, intent(inout) :: Varray(:,:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 3D double precision array (KIND=dp).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
# ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
# else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
# endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_3D_DP - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_3d_dp
!
!***********************************************************************
      FUNCTION destroy_4d_dp (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(dp), pointer, intent(inout) :: Varray(:,:,:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 4D double precision array (KIND=dp).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
# ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
# else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
# endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_4D_DP - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_4d_dp
!
!***********************************************************************
      FUNCTION destroy_5d_dp (ng, Varray, routine, line, Vstring)       &
     &                RESULT (Lsuccess)
!***********************************************************************
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, line
!
      real(dp), pointer, intent(inout) :: Varray(:,:,:,:,:)
!
      character (len=*) :: Vstring, routine
!
!  Local variable declarations.
!
      logical :: Lsuccess
!
      integer :: Derror
!
      character (len=:), allocatable :: Dmsg
!
!-----------------------------------------------------------------------
!  Deallocate 5D double precision array (KIND=dp).
!-----------------------------------------------------------------------
!
      Lsuccess=.TRUE.
      Derror=0
!
# ifdef DISTRIBUTE
      IF (associated(Varray))                                           &
     &  deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
# else
      IF (Master) THEN
        IF (associated(Varray))                                         &
     &    deallocate ( Varray, ERRMSG = Dmsg, STAT = Derror )
      END IF
# endif
!
!  Report if unsuccessful deallocation.
!
      IF (Derror.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,10) ng, Vstring, routine, line, Dmsg
        END IF
        Lsuccess=.FALSE.
      END IF
!
 10   FORMAT (/,' DESTROY_5D_DP - Grid ',i2.2,                          &
     &        ', error while deallocating: ''',a,''' in routine ''',a,  &
     &        ''' at line = ',i0,/,17x,a)
!
      RETURN
      END FUNCTION destroy_5d_dp
#endif
!
      END MODULE destroy_mod
