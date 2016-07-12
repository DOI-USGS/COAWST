#include "cppdefs.h"
      MODULE distribute_mod
#ifdef DISTRIBUTE

# define BOUNDARY_ALLREDUCE /* use mpi_allreduce in mp_boundary */
# undef  COLLECT_ALLGATHER  /* use mpi_allgather in mp_collect  */
# define COLLECT_ALLREDUCE  /* use mpi_allreduce in mp_collect  */
# define REDUCE_ALLGATHER   /* use mpi_allgather in mp_reduce   */
# undef  REDUCE_ALLREDUCE   /* use mpi_allreduce in mp_reduce   */
!
!svn $Id: distribute.F 795 2016-05-11 01:42:43Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines are used for distrubuted-memory communications       !
!  between parallel nodes:                                             !
!                                                                      !
!  mp_aggregate2d    aggregates 2D tiled data into a 2D global array   !
!  mp_aggregate3d    aggregates 3D tiled data into a 3D global array   !
!  mp_barrier        barrier sychronization                            !
!  mp_bcastf         broadcasts floating point variables               !
!  mp_bcasti         broadcasts integer variables                      !
!  mp_bcastl         broadcasts logical variables                      !
!  mp_bcasts         broadcasts character variables                    !
!  mp_boundary       exchanges boundary data between tiles             !
!  mp_assemblef_1d   assembles 1D floating point array from tiles      !
!  mp_assemblef_2d   assembles 2D floating point array from tiles      !
!  mp_assemblef_3d   assembles 3D floating point array from tiles      !
!  mp_assemblei_1d   assembles 1D integer array from tiles             !
!  mp_assemblei_2d   assembles 2D integer array from tiles             !
!  mp_collect_f      collects 1D floating point array from tiles       !
!  mp_collect_i      collects 1D integer array from tiles              !
!  mp_collect_l      collects 1D logical vector data from tiles        !
!  mp_dump           writes 2D and 3D tiles arrays for debugging       !
!  mp_gather2d       collects a 2D tiled array for output purposes     !
!  mp_gather3d       collects a 3D tiled array for output purposes     !
!  mp_gather_state   collects state vector for unpacking of variables  !
!  mp_ncread1d       reads  in  1D state array from NetCDF file        !
!  mp_ncread2d       reads  in  2D state array from NetCDF file        !
!  mp_ncwrite1d      writes out 1D state array into NetCDF file        !
!  mp_ncwrite2d      writes out 2D state array into NetCDF file        !
!  mp_reduce         global reduction operations                       !
!  mp_reduce2        global reduction operations (MINLOC, MAXLOC)      !
!  mp_scatter2d      scatters input data to a 2D tiled array           !
!  mp_scatter3d      scatters input data to a 3D tiled array           !
!  mp_scatter_state  scatters global data for packing of state vector  !
!                                                                      !
!  Notice that the tile halo exchange can be found in "mp_exchange.F"  !
!                                                                      !
!=======================================================================
!
      implicit none
!
      INTERFACE mp_assemble
        MODULE PROCEDURE mp_assemblef_1d
        MODULE PROCEDURE mp_assemblef_2d
        MODULE PROCEDURE mp_assemblef_3d
        MODULE PROCEDURE mp_assemblei_1d
        MODULE PROCEDURE mp_assemblei_2d
      END INTERFACE mp_assemble
!
      INTERFACE mp_bcastf
        MODULE PROCEDURE mp_bcastf_0d
        MODULE PROCEDURE mp_bcastf_1d
        MODULE PROCEDURE mp_bcastf_2d
        MODULE PROCEDURE mp_bcastf_3d
        MODULE PROCEDURE mp_bcastf_4d
      END INTERFACE mp_bcastf
!
      INTERFACE mp_bcastl
        MODULE PROCEDURE mp_bcastl_0d
        MODULE PROCEDURE mp_bcastl_1d
        MODULE PROCEDURE mp_bcastl_2d
      END INTERFACE mp_bcastl
!
      INTERFACE mp_bcasti
        MODULE PROCEDURE mp_bcasti_0d
        MODULE PROCEDURE mp_bcasti_1d
        MODULE PROCEDURE mp_bcasti_2d
      END INTERFACE mp_bcasti
!
      INTERFACE mp_bcasts
        MODULE PROCEDURE mp_bcasts_0d
        MODULE PROCEDURE mp_bcasts_1d
      END INTERFACE mp_bcasts
!
      INTERFACE mp_collect
        MODULE PROCEDURE mp_collect_f
        MODULE PROCEDURE mp_collect_i
        MODULE PROCEDURE mp_collect_l
      END INTERFACE mp_collect
!
      INTERFACE mp_reduce
        MODULE PROCEDURE mp_reduce_0d
        MODULE PROCEDURE mp_reduce_1d
      END INTERFACE mp_reduce
!
      CONTAINS
!
      SUBROUTINE mp_barrier (ng)
!
!***********************************************************************
!                                                                      !
!  This routine blocks the caller until all group members have called  !
!  it.                                                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: MyError
!
!-----------------------------------------------------------------------
!  Synchronize all distribute-memory nodes in the group.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_barrier (OCN_COMM_WORLD, MyError)
# endif

      RETURN
      END SUBROUTINE mp_barrier
!
      SUBROUTINE mp_bcastf_0d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a floating-point scalar variable to all     !
!  processors the in group. It is called by all the members in the     !
!  group.                                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      real(r8), intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, 1, MP_FLOAT, MyMaster, OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastf_0d
!
      SUBROUTINE mp_bcastf_1d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D floating-point, nontiled, array to     !
!  all processors processors in the group. It is called by all the     !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      real(r8), intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=UBOUND(A, DIM=1)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastf_1d
!
      SUBROUTINE mp_bcastf_2d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D floating-point, nontiled, array to     !
!  all processors processors in the group. It is called by all the     !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      real(r8), intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror

      integer :: Asize(2)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Npts=Asize(1)*Asize(2)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastf_2d
!
      SUBROUTINE mp_bcastf_3d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 3D floating-point, nontiled, array to     !
!  all processors processors in the group. It is called by all the     !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          3D array to broadcast (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 3D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      real(r8), intent(inout) :: A(:,:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror

      integer :: Asize(3)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      Npts=Asize(1)*Asize(2)*Asize(3)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_3D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastf_3d
!
      SUBROUTINE mp_bcastf_4d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 4D floating-point, nontiled, array to     !
!  all processors processors in the group. It is called by all the     !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          4D array to broadcast (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 4D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      real(r8), intent(inout) :: A(:,:,:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror

      integer :: Asize(4)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      Asize(4)=UBOUND(A, DIM=4)
      Npts=Asize(1)*Asize(2)*Asize(3)*Asize(4)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTF_4D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastf_4d
!
      SUBROUTINE mp_bcasti_0d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts an integer scalar variable to all           !
!  processors the in group.  It is called by all the members           !
!  in the group.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (integer).                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, 1, MPI_INTEGER, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTI_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcasti_0d
!
      SUBROUTINE mp_bcasti_1d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D nontiled, integer array to all         1
!  processors processors in the group. It is called by all the         !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (integer).                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=UBOUND(A, DIM=1)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MPI_INTEGER, MyMaster, OCN_COMM_WORLD,   &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTI_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcasti_1d
!
      SUBROUTINE mp_bcasti_2d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D nontiled, integer array to all         1
!  processors processors in the group. It is called by all the         !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (integer).                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: Asize(2)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Npts=Asize(1)*Asize(2)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MPI_INTEGER, MyMaster, OCN_COMM_WORLD,   &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTI_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcasti_2d
!
      SUBROUTINE mp_bcastl_0d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a logical scalar variable to all            !
!  processors the in group. It is called by all the members            !
!  in the group.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (logical).                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      logical, intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
# ifdef MPI
      CALL mpi_bcast (A, 1, MPI_LOGICAL, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTL_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastl_0d
!
      SUBROUTINE mp_bcastl_1d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D nontiled, logical array to all         !
!  processors processors in the group. It is called by all the         !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (logical).                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      logical, intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Npts=UBOUND(A, DIM=1)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MPI_LOGICAL, MyMaster, OCN_COMM_WORLD,   &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTL_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastl_1d
!
      SUBROUTINE mp_bcastl_2d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 2D nontiled, logical array to all         !
!  processors processors in the group. It is called by all the         !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          2D array to broadcast (logical).                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 2D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      logical, intent(inout) :: A(:,:)
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: Asize(2)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Npts=Asize(1)*Asize(2)

# ifdef MPI
      CALL mpi_bcast (A, Npts, MPI_LOGICAL, MyMaster, OCN_COMM_WORLD,   &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTL_2D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcastl_2d
!
      SUBROUTINE mp_bcasts_0d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a string scalar variable to all processors  !
!  in the group. It is called by all the members in the group.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          Variable to broadcast (string).                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted variable.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      character (len=*), intent(inout) :: A
!
!  Local variable declarations
!
      integer :: Lstr, MyError, Nchars, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      IF (Lwclock) THEN
        CALL wclock_on (ng, model, 44)
      END IF
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Nchars=LEN(A)
# ifdef MPI
      CALL mpi_bcast (A, Nchars, MPI_BYTE, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTS_0D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      IF (Lwclock) THEN
        CALL wclock_off (ng, model, 44)
      END IF
# endif

      RETURN
      END SUBROUTINE mp_bcasts_0d
!
      SUBROUTINE mp_bcasts_1d (ng, model, A)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts a 1D nontiled,  string array to all         !
!  processors processors in the group. It is called by all the         !
!  members in the group.                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     A          1D array to broadcast (string).                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Broadcasted 1D array.                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      character (len=*), intent(inout) :: A(:)
!
!  Local variable declarations
!
      integer :: Asize, Lstr, MyError, Nchars, Serror

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 44)
# endif
!
!-----------------------------------------------------------------------
!  Broadcast requested variable.
!-----------------------------------------------------------------------
!
      Asize=UBOUND(A, DIM=1)
      Nchars=LEN(A(1))*Asize

# ifdef MPI
      CALL mpi_bcast (A, Nchars, MPI_BYTE, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_BCASTS_1D - error during ',a,' call, Node = ',   &
     &          i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 44)
# endif

      RETURN
      END SUBROUTINE mp_bcasts_1d
!
      SUBROUTINE mp_boundary (ng, model, Imin, Imax,                    &
     &                        LBi, UBi, LBk, UBk,                       &
     &                        update, A)
!
!***********************************************************************
!                                                                      !
!  This routine exchanges boundary arrays between tiles.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Imin       Starting tile index.                                  !
!     Imax       Ending   tile index.                                  !
!     Jstr       Starting tile index in the J-direction.               !
!     Jend       Ending   tile index in the J-direction.               !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     LBk        K-dimension Upper bound, if any. Otherwise, a value   !
!                  of one is expected.                                 !
!     UBk        K-dimension Upper bound.                              !
!     update     Switch activated by the node that updated the         !
!                  boundary data.                                      !
!     A          Boundary array (1D or 2D) to process.                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Updated boundary array (1D or 2D).                    !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: update

      integer, intent(in) :: ng, model, Imin, Imax
      integer, intent(in) :: LBi, UBi, LBk, UBk

      real(r8), intent(inout) :: A(LBi:UBi,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Ilen, Ioff, Lstr, MyError, Nnodes, Npts, Serror
      integer :: i, ik, k, kc, rank

      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1)) :: Asend

# ifdef BOUNDARY_ALLREDUCE
      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1)) :: Arecv
# else
      real(r8), dimension((UBi-LBi+1)*(UBk-LBk+1),                      &
     &                    0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
# endif
      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 48)
# endif
!
!-----------------------------------------------------------------------
!  Pack boundary data.  Zero-out boundary array except points updated
!  by the appropriate node, so sum reduction can be perfomed during
!  unpacking.
!-----------------------------------------------------------------------
!
!  Initialize buffer to the full range so unpacking is correct with
!  summation.  This also allows even exchange of segments with
!  communication routine "mpi_allgather".
!
      Ilen=UBi-LBi+1
      Ioff=1-LBi
      Npts=Ilen*(UBk-LBk+1)
      DO i=1,Npts
        Asend(i)=0.0_r8
      END DO
!
!  If a boundary tile, load boundary data.
!
      IF (update) THEN
        DO k=LBk,UBk
          kc=(k-LBk)*Ilen
          DO i=Imin,Imax
            ik=i+Ioff+kc
            Asend(ik)=A(i,k)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
# ifdef MPI
#  ifdef BOUNDARY_ALLREDUCE
      CALL mpi_allreduce (Asend, Arecv, Npts, MP_FLOAT, MPI_SUM,        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
 10     FORMAT (/,' MP_BOUNDARY - error during ',a,' call, Node = ',    &
     &          i3.3,' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
#  else
      CALL mpi_allgather (Asend, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT, &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 10     FORMAT (/,' MP_BOUNDARY - error during ',a,' call, Node = ',    &
     &          i3.3,' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Unpack data: reduction sum.
!-----------------------------------------------------------------------
!
# ifdef BOUNDARY_ALLREDUCE
      ik=0
      DO k=LBk,UBk
        DO i=LBi,UBi
          ik=ik+1
          A(i,k)=Arecv(ik)
        END DO
      END DO
# else
      Nnodes=NtileI(ng)*NtileJ(ng)-1
      ik=0
      DO k=LBk,UBk
        DO i=LBi,UBi
          A(i,k)=0.0_r8
          ik=ik+1
          DO rank=0,Nnodes
            A(i,k)=A(i,k)+Arecv(ik,rank)
          END DO
        END DO
      END DO
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 48)
# endif

      RETURN
      END SUBROUTINE mp_boundary
!
      SUBROUTINE mp_assemblef_1d (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 1D floating-point array from all members   !
!  in the group.  The collection of data from all nodes is achieved    !
!  as a reduction sum.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          1D array to collect.                                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 1D array.                                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      real(r8), intent(in) :: Aspv

      real(r8), intent(inout) :: A(:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyNpts, Nnodes, Serror, i

      integer, dimension(MPI_STATUS_SIZE) :: status

      real(r8), dimension(Npts) :: Asend

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
      MyNpts=UBOUND(A, DIM=1)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0.0_r8) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Set input data to process.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MP_FLOAT, MPI_SUM,            &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif
!
 10   FORMAT (/,' MP_ASSEMBLEF_1D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEF_1D - illegal special value, Aspv = ',    &
     &        1p,e17.10,/,19x,'a zero value is needed for global ',     &
     &        'reduction.')
 30   FORMAT (/,' MP_ASSEMBLEF_1D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)

      RETURN
      END SUBROUTINE mp_assemblef_1d
!
      SUBROUTINE mp_assemblef_2d (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 2D floating-point array from all members   !
!  in the group.  The collection of data from all nodes is achieved    !
!  as a reduction sum.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          2D array to collect.                                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 2D array.                                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      real(r8), intent(in) :: Aspv

      real(r8), intent(inout) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyNpts, Nnodes, Serror

      integer :: Asize(2)

      integer, dimension(MPI_STATUS_SIZE) :: status

      real(r8), dimension(Npts) :: Asend, Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      MyNpts=Asize(1)*Asize(2)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0.0_r8) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Set input data to process.
!
      Asend=RESHAPE(A, (/Npts/))
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MP_FLOAT, MPI_SUM,        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load collected data.
!
      A=RESHAPE(Arecv, Asize)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif
!
 10   FORMAT (/,' MP_ASSEMBLEF_2D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEF_2D - illegal special value, Aspv = ',    &
     &        1p,e17.10,/,19x,'a zero value is needed for global ',     &
     &        'reduction.')
 30   FORMAT (/,' MP_ASSEMBLEF_2D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)

      RETURN
      END SUBROUTINE mp_assemblef_2d
!
      SUBROUTINE mp_assemblef_3d (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 3D floating-point array from all members   !
!  in the group.  The collection of data from all nodes is achieved    !
!  as a reduction sum.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          3D array to collect.                                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 3D array.                                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      real(r8), intent(in) :: Aspv

      real(r8), intent(inout) :: A(:,:,:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyNpts, Nnodes, Serror

      integer :: Asize(3)

      integer, dimension(MPI_STATUS_SIZE) :: status

      real(r8), dimension(Npts) :: Asend, Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      MyNpts=Asize(1)*Asize(2)*Asize(3)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0.0_r8) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Set input data to process.
!
      Asend=RESHAPE(A, (/Npts/))
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MP_FLOAT, MPI_SUM,        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load collected data.
!
      A=RESHAPE(Arecv, Asize)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif
!
 10   FORMAT (/,' MP_ASSEMBLEF_3D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEF_3D - illegal special value, Aspv = ',    &
     &        1p,e17.10,/,19x,'a zero value is needed for global ',     &
     &        'reduction.')
 30   FORMAT (/,' MP_ASSEMBLEF_3D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)

      RETURN
      END SUBROUTINE mp_assemblef_3d
!
!
      SUBROUTINE mp_assemblei_1d (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 1D integer array from all members in the   !
!  group.  The collection of data from all nodes is achieved as a      !
!  reduction sum.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          1D array to collect.                                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 1D array.                                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      integer, intent(in) :: Aspv

      integer, intent(inout) :: A(:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyNpts, Nnodes, Serror, i

      integer, dimension(MPI_STATUS_SIZE) :: status

      integer, dimension(Npts) :: Asend

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
      MyNpts=UBOUND(A, DIM=1)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Set input data to process.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MPI_INTEGER, MPI_SUM,         &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif
!
 10   FORMAT (/,' MP_ASSEMBLEI_1D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEI_1D - illegal special value, Aspv = ',i4, &
     &        /,19x,'a zero value is needed for global reduction.')
 30   FORMAT (/,' MP_ASSEMBLEI_1D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)

      RETURN
      END SUBROUTINE mp_assemblei_1d
!
      SUBROUTINE mp_assemblei_2d (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine assembles a 2D integer array from all members in the   !
!  group.  The collection of data from all nodes is achieved as a      !
!  reduction sum.                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points, PROD(SIZE(A)).       !
!     Aspv       Special value indicating that an array element is     !
!                  not operated by the current parallel node. It must  !
!                  be zero to collect data by a global reduction sum.  !
!     A          2D array to collect.                                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Assembled 2D array.                                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      integer, intent(in) :: Aspv

      integer, intent(inout) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyNpts, Nnodes, Serror

      integer :: Asize(2)

      integer, dimension(MPI_STATUS_SIZE) :: status

      integer, dimension(Npts) :: Asend, Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Check input parameters.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      MyNpts=Asize(1)*Asize(2)
      IF (Npts.ne.MyNpts) THEN
        IF (Master) THEN
          WRITE (stdout,10) Npts, MyNpts
        END IF
        exit_flag=7
      END IF
!
      IF (Aspv.ne.0) THEN
        IF (Master) THEN
          WRITE (stdout,20) Aspv
        END IF
        exit_flag=7
      END IF
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Set input data to process.
!
      Asend=RESHAPE(A, (/Npts/))
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, Arecv, Npts, MPI_INTEGER, MPI_SUM,     &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,30) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load collected data.
!
      A=RESHAPE(Arecv, Asize)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif
!
 10   FORMAT (/,' MP_ASSEMBLEI_2D - inconsistent array size, Npts = ',  &
     &        i10,2x,i10,/,19x,'number of addressed array elements ',   &
     &        'is incorrect.')
 20   FORMAT (/,' MP_ASSEMBLEI_2D - illegal special value, Aspv = ',i4, &
     &        /,19x,'a zero value is needed for global reduction.')
 30   FORMAT (/,' MP_ASSEMBLEI_2D - error during ',a,' call, Node = ',  &
     &        i3.3,' Error = ',i3,/,19x,a)

      RETURN
      END SUBROUTINE mp_assemblei_2d
!
      SUBROUTINE mp_collect_f (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 1D floating-point array from all members    !
!  in the group. Then, it packs distributed data by removing the       !
!  special values. This routine is used when extracting station        !
!  data from tiled arrays.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points.                      !
!     Aspv       Special value indicating no data.  This implies that  !
!                  desired data is tile unbouded.                      !
!     A          Collected data.                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data.                                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      real(r8), intent(in) :: Aspv

      real(r8), intent(inout) :: A(Npts)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Nnodes, Serror
      integer :: i, rank, request

      integer, dimension(MPI_STATUS_SIZE) :: status

# if defined COLLECT_ALLGATHER
      real(r8), dimension(Npts,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
# elif defined COLLECT_ALLREDUCE
      real(r8), dimension(Npts) :: Asend
# else
      real(r8), allocatable :: Arecv(:)
# endif

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
# if defined COLLECT_ALLGATHER
      CALL mpi_allgather (A, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT,     &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Pack data according to special values: sum or ignore.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
      IF (Aspv.eq.0.0_r8) THEN
        DO i=1,Npts
          A(i)=0.0_r8
          DO rank=0,Nnodes
            A(i)=A(i)+Arecv(i,rank)
          END DO
        END DO
      ELSE
        DO i=1,Npts
          DO rank=0,Nnodes
            IF (Arecv(i,rank).ne.Aspv) THEN
              A(i)=Arecv(i,rank)
            END IF
          END DO
        END DO
      END IF
# elif defined COLLECT_ALLREDUCE
!
!  Copy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MP_FLOAT, MPI_SUM,            &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# else
      IF (MyRank.eq.MyMaster) THEN
!
!  If master node, allocate and receive buffer.
!
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  If master node, loop over other nodes to receive and accumulate the
!  data.
!
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv, Npts, MP_FLOAT, rank, rank+5,          &
     &                    OCN_COMM_WORLD, request, MyError)
          CALL mpi_wait (request, status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
          DO i=1,Npts
            A(i)=A(i)+Arecv(i)
          END DO
        END DO
        deallocate (Arecv)
!
!  Otherwise, send data to master node.
!
      ELSE
        CALL mpi_isend (A, Npts, MP_FLOAT, MyMaster, MyRank+5,          &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast accumulated (full) data to all nodes.
!
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# endif
 10   FORMAT (/,' MP_COLLECT_F - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,14x,a)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif

      RETURN
      END SUBROUTINE mp_collect_f
!
      SUBROUTINE mp_collect_i (ng, model, Npts, Aspv, A)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 1D integer array from all members in        !
!  the group. Then, it packs distributed data by removing the          !
!  special values. This routine is used when extracting station        !
!  data from tiled arrays.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of collected data points.                      !
!     Aspv       Special value indicating no data.  This implies that  !
!                  desired data is tile unbouded.                      !
!     A          Collected data.                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data.                                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      integer, intent(in) :: Aspv

      integer, intent(inout) :: A(Npts)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Nnodes, Serror
      integer :: i, rank, request

      integer, dimension(MPI_STATUS_SIZE) :: status

# if defined COLLECT_ALLGATHER
      integer, dimension(Npts,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
# elif defined COLLECT_ALLREDUCE
      integer, dimension(Npts) :: Asend
# else
      integer, allocatable :: Arecv(:)
# endif

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 49)
# endif
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
# if defined COLLECT_ALLGATHER
      CALL mpi_allgather (A, Npts, MPI_INTEGER, Arecv, Npts,            &
     &                    MPI_INTEGER, OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Pack data according to special values: sum or ignore.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
      IF (Aspv.eq.0) THEN
        DO i=1,Npts
          A(i)=0
          DO rank=0,Nnodes
            A(i)=A(i)+Arecv(i,rank)
          END DO
        END DO
      ELSE
        DO i=1,Npts
          DO rank=0,Nnodes
            IF (Arecv(i,rank).ne.Aspv) THEN
              A(i)=Arecv(i,rank)
            END IF
          END DO
        END DO
      END IF
# elif defined COLLECT_ALLREDUCE
!
!  Copy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MPI_INTEGER, MPI_SUM,         &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# else
      IF (MyRank.eq.MyMaster) THEN
!
!  If master node, allocate and receive buffer.
!
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  If master node, loop over other nodes to receive and accumulate the
!  data.
!
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv, Npts, MPI_INTEGER, rank, rank+5,       &
     &                    OCN_COMM_WORLD, request, MyError)
          CALL mpi_wait (request, status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
          DO i=1,Npts
            A(i)=A(i)+Arecv(i)
          END DO
        END DO
        deallocate (Arecv)
!
!  Otherwise, send data to master node.
!
      ELSE
        CALL mpi_isend (A, Npts, MPI_INTEGER, MyMaster, MyRank+5,       &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast accumulated (full) data to all nodes.
!
      CALL mpi_bcast (A, Npts, MPI_INTEGER, MyMaster, OCN_COMM_WORLD,   &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# endif
 10   FORMAT (/,' MP_COLLECT_I - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,14x,a)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 49)
# endif
      RETURN
      END SUBROUTINE mp_collect_i

      SUBROUTINE mp_collect_l (ng, model, Npts, A)
!
!***********************************************************************
!                                                                      !
!  This routine collects requested buffer from all members in the      !
!  group. Then, it packs distributed data by removing the special      !
!  values. This routine is used when extracting station data from      !
!  tiled arrays.                                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Npts       Number of extracted data points.                      !
!     A          Extracted data.                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data.                                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Npts

      logical, intent(inout) :: A(Npts)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Nnodes, Serror
      integer :: i, rank, request

      integer, dimension(MPI_STATUS_SIZE) :: status

# if defined COLLECT_ALLGATHER
      logical, dimension(Npts,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
# elif defined COLLECT_ALLREDUCE
      logical, dimension(Npts) :: Asend
# else
      logical, allocatable :: Arecv(:)
# endif

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 47)
# endif
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
# if defined COLLECT_ALLGATHER
      CALL mpi_allgather (A, Npts, MPI_LOGICAL, Arecv, Npts,            &
     &                    MPI_LOGICAL, OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Pack data according to special values: sum or ignore.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
      DO i=1,Npts
        DO rank=0,Nnodes
          A(i)=Arecv(i,rank)
        END DO
      END DO
# elif defined COLLECT_ALLREDUCE
!
!  Copy data to send.
!
      DO i=1,Npts
        Asend(i)=A(i)
      END DO
!
!  Collect data from all nodes as a reduced sum.
!
      CALL mpi_allreduce (Asend, A, Npts, MPI_LOGICAL, MPI_LOR,         &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# else
      IF (MyRank.eq.MyMaster) THEN
!
!  If master node, allocate and receive buffer.
!
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  If master node, loop over other nodes to receive and accumulate the
!  data.
!
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv, Npts, MPI_LOGICAL, rank, rank+5,       &
     &                    OCN_COMM_WORLD, request, MyError)
          CALL mpi_wait (request, status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
          DO i=1,Npts
            A(i)=A(i)+Arecv(i)
          END DO
        END DO
        deallocate (Arecv)
!
!  Otherwise, send data to master node.
!
      ELSE
        CALL mpi_isend (A, Npts, MPI_LOGICAL, MyMaster, MyRank+5,       &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast accumulated (full) data to all nodes.
!
      CALL mpi_bcast (A, Npts, MP_FLOAT, MyMaster, OCN_COMM_WORLD,      &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# endif
 10   FORMAT (/,' MP_COLLECT - error during ',a,' call, Node = ',       &
     &        i3.3,' Error = ',i3,/,14x,a)

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 47)
# endif

      RETURN
      END SUBROUTINE mp_collect_l
!
      SUBROUTINE mp_gather2d (ng, model, LBi, UBi, LBj, UBj,            &
     &                        tindex, gtype, Ascl,                      &
# ifdef MASKING
     &                        Amask,                                    &
# endif
     &                        A, Npts, Awrk, SetFillVal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 2D tiled, floating-point array from each    !
!  spawned node and stores it into one dimensional global array. It    !
!  is used to collect and  pack output data.                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     tindex     Time record index to process.                         !
!     gtype      C-grid type. If negative and Land-Sea is available,   !
!                  only water-points processed.                        !
!     Ascl       Factor to scale field before writing.                 !
!     Amask      Land/Sea mask, if any.                                !
!     A          2D tiled, floating-point array to process.            !
!     SetFillVal Logical switch to set fill value in land areas        !
!                  (optional).                                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts       Number of points processed in Awrk.                   !
!     Awrk       Collected data from each node packed into 1D array    !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(out) :: Npts

      real(r8), intent(in) :: Ascl

# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in)  :: A(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
# endif
      integer :: Cgrid, ghost, rank
      integer :: Io, Ie, Jo, Je, Ioff, Joff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen
      integer :: Lstr, MyError, MyType, Serror, Srequest
      integer :: i, ic, j, jc, np

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: MySize
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus

      real(r8), dimension(TileSize(ng)) :: Asend

      real(r8), dimension(TileSize(ng),                                 &
     &                    NtileI(ng)*NtileJ(ng)-1) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT

      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Npts=Ilen*Jlen
!
!  Set physical, non-overlapping (no ghost-points) ranges according to
!  tile rank.
!
      ghost=0

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT                              ! RHO-points
          Cgrid=2
      END SELECT

      Imin=BOUNDS(ng) % Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng) % Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,MyRank)
!
!  Compute size of distributed buffers.
!
      DO rank=0,NtileI(ng)*NtileJ(ng)-1
        MySize(rank)=(BOUNDS(ng) % Imax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Imin(Cgrid,ghost,rank)+1)*           &
     &               (BOUNDS(ng) % Jmax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Jmin(Cgrid,ghost,rank)+1)
      END DO
!
!  Initialize local arrays to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Collect requested array data.
!-----------------------------------------------------------------------
!
!  Pack and scale input data.
!
      np=0
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          np=np+1
          Asend(np)=A(i,j)*Ascl
        END DO
      END DO

# ifdef MASKING
!
!  If overwriting Land/Sea mask or processing water-points only, flag
!  land-points with special value.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
      IF (gtype.lt.0) THEN
        np=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            np=np+1
            IF (Amask(i,j).eq.0.0_r8) THEN
              Asend(np)=spval
            END IF
          END DO
        END DO
      ELSE IF (LandFill) THEN
        np=0
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            np=np+1
            IF (Amask(i,j).eq.0.0_r8) THEN
              Asend(np)=spval
            END IF
          END DO
        END DO
      END IF
# endif
!
!  If master processor, unpack the send buffer since there is not
!  need to distribute.
!
      IF (MyRank.eq.MyMaster) THEN
        np=0
        DO j=Jmin,Jmax
          jc=(j-Joff)*Ilen
          DO i=Imin,Imax
            np=np+1
            ic=i+Ioff+jc
            Awrk(ic)=Asend(np)
          END DO
        END DO
      END IF
!
!  Send, receive, and unpack data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), MySize(rank), MP_FLOAT, rank,  &
     &                    rank+5, OCN_COMM_WORLD, Rrequest(rank),       &
     &                    MyError)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
 10         FORMAT (/,' MP_GATHER2D - error during ',a,' call, Node = ',&
     &              i3.3,' Error = ',i3,/,13x,a)
            exit_flag=2
            RETURN
          END IF

          np=0
          Imin=BOUNDS(ng) % Imin(Cgrid,ghost,rank)
          Imax=BOUNDS(ng) % Imax(Cgrid,ghost,rank)
          Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,rank)
          Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,rank)

          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen
            DO i=Imin,Imax
              np=np+1
              ic=i+Ioff+jc
              Awrk(ic)=Arecv(np,rank)
            END DO
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, MySize(MyRank), MP_FLOAT, MyMaster,      &
     &                  MyRank+5, OCN_COMM_WORLD, Srequest, MyError)
        CALL mpi_wait (Srequest, Sstatus, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF

# ifdef MASKING
!
! If pocessing only water-points, remove land points and repack.
!
      IF ((MyRank.eq.MyMaster).and.(gtype.lt.0)) THEN
        ic=0
        np=Ilen*Jlen
        DO i=1,np
          IF (Awrk(i).lt.spval) THEN
            ic=ic+1
            Awrk(ic)=Awrk(i)
          END IF
        END DO
        Npts=ic
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

      RETURN
      END SUBROUTINE mp_gather2d
!
      SUBROUTINE mp_gather3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk,  &
     &                        tindex, gtype, Ascl,                      &
# ifdef MASKING
     &                        Amask,                                    &
# endif
     &                        A, Npts, Awrk, SetFillVal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 3D tiled, floating-point array from each    !
!  spawned node and stores it into one dimensional global array. It    !
!  is used to collect and  pack output data.                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     tindex     Time record index to process.                         !
!     gtype      C-grid type. If negative and Land-Sea is available,   !
!                  only water-points processed.                        !
!     Ascl       Factor to scale field before writing.                 !
!     Amask      Land/Sea mask, if any.                                !
!     A          3D tiled, floating-point array to process.            !
!     SetFillVal Logical switch to set fill value in land areas        !
!                  (optional).                                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Npts       Number of points processed in Awrk.                   !
!     Awrk       Collected data from each node packed into 1D array    !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, tindex, gtype
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(out) :: Npts

      real(r8), intent(in) :: Ascl

# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in)  :: A(LBi:UBi,LBj:UBj,LBk:UBk)
      real(r8), intent(out) :: Awrk(:)
!
!  Local variable declarations.
!
# ifdef MASKING
      logical :: LandFill
# endif
      integer :: Cgrid, ghost, rank
      integer :: Io, Ie, Jo, Je, Ioff, Joff, Koff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: Lstr, MyError, MyType, Serror, Srequest
      integer :: i, ic, j, jc, k, kc, np

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: MySize
      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus

      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)) :: Asend

      real(r8), dimension(TileSize(ng)*(UBk-LBk+1),                     &
     &                    NtileI(ng)*NtileJ(ng)-1) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT

      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF

      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
      Npts=IJlen*Klen
!
!  Set tile physical, non-overlapping (no ghost-points) ranges according
!  to tile rank.
!
      ghost=0

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Cgrid=1
        CASE (r2dvar, r3dvar)
          Cgrid=2
        CASE (u2dvar, u3dvar)
          Cgrid=3
        CASE (v2dvar, v3dvar)
          Cgrid=4
        CASE DEFAULT                              ! RHO-points
          Cgrid=2
      END SELECT

      Imin=BOUNDS(ng) % Imin(Cgrid,ghost,MyRank)
      Imax=BOUNDS(ng) % Imax(Cgrid,ghost,MyRank)
      Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,MyRank)
      Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,MyRank)
!
!  Compute size of distributed buffers.
!
      DO rank=0,NtileI(ng)*NtileJ(ng)-1
        MySize(rank)=(BOUNDS(ng) % Imax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Imin(Cgrid,ghost,rank)+1)*           &
     &               (BOUNDS(ng) % Jmax(Cgrid,ghost,rank)-              &
     &                BOUNDS(ng) % Jmin(Cgrid,ghost,rank)+1)*           &
     &               (UBk-LBk+1)
      END DO
!
!  Initialize local arrays to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Collect requested array data.
!-----------------------------------------------------------------------
!
!  Pack and scale input data.
!
      np=0
      DO k=LBk,UBk
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            np=np+1
            Asend(np)=A(i,j,k)*Ascl
          END DO
        END DO
      END DO

# ifdef MASKING
!
!  If overwriting Land/Sea mask or processing water-points only, flag
!  land-points with special value.
!
      IF (PRESENT(SetFillVal)) THEN
        LandFill=SetFillVal
      ELSE
        LandFill=tindex.gt.0
      END IF
      IF (gtype.lt.0) THEN
        np=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              np=np+1
              IF (Amask(i,j).eq.0.0_r8) THEN
                Asend(np)=spval
              END IF
            END DO
          END DO
        END DO
      ELSE IF (LandFill) THEN
        np=0
        DO k=LBk,UBk
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              np=np+1
              IF (Amask(i,j).eq.0.0_r8) THEN
                Asend(np)=spval
              END IF
            END DO
          END DO
        END DO
      END IF
# endif
!
!  If master processor, unpack the send buffer since there is not
!  need to distribute.
!
      IF (MyRank.eq.MyMaster) THEN
        np=0
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          DO j=Jmin,Jmax
            jc=(j-Joff)*Ilen+kc
            DO i=Imin,Imax
              np=np+1
              ic=i+Ioff+jc
              Awrk(ic)=Asend(np)
            END DO
          END DO
        END DO
      END IF
!
!  Send, receive, and unpack data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), MySize(rank), MP_FLOAT, rank,  &
     &                    rank+5, OCN_COMM_WORLD, Rrequest(rank),       &
     &                    MyError)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
 10         FORMAT (/,' MP_GATHER3D - error during ',a,' call, Node = ',&
     &              i3.3,' Error = ',i3,/,13x,a)
            exit_flag=2
            RETURN
          END IF

          np=0
          Imin=BOUNDS(ng) % Imin(Cgrid,ghost,rank)
          Imax=BOUNDS(ng) % Imax(Cgrid,ghost,rank)
          Jmin=BOUNDS(ng) % Jmin(Cgrid,ghost,rank)
          Jmax=BOUNDS(ng) % Jmax(Cgrid,ghost,rank)

          DO k=LBk,UBk
            kc=(k-Koff)*IJlen
            DO j=Jmin,Jmax
              jc=(j-Joff)*Ilen+kc
              DO i=Imin,Imax
                np=np+1
                ic=i+Ioff+jc
                Awrk(ic)=Arecv(np,rank)
              END DO
            END DO
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, MySize(MyRank), MP_FLOAT, MyMaster,      &
     &                  MyRank+5, OCN_COMM_WORLD, Srequest, MyError)
        CALL mpi_wait (Srequest, Sstatus, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF

# ifdef MASKING
!
! If pocessing only water-points, remove land points and repack.
!
      IF ((MyRank.eq.MyMaster).and.(gtype.lt.0)) THEN
        ic=0
        np=IJlen*Klen
        DO i=1,np
          IF (Awrk(i).lt.spval) THEN
            ic=ic+1
            Awrk(ic)=Awrk(i)
          END IF
        END DO
        Npts=ic
      END IF
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

      RETURN
      END SUBROUTINE mp_gather3d
!
      SUBROUTINE mp_gather_state (ng, model, Mstr, Mend, Asize,         &
     &                            A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine gathers (threaded to global) state data to all nodes   !
!  in the group. This  routine  is used to unpack the state data for   !
!  the GST analysis propagators.                                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Mstr       Threaded array lower bound.                           !
!     Mend       Threaded array upper bound.                           !
!     Asize      Size of the full state.                               !
!     A          Threaded 1D array process.                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Awrk       Collected data from each node packed into 1D full     !
!                  state array.                                        !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend, Asize

      real(r8), intent(in)  :: A(Mstr:Mend)
      real(r8), intent(out) :: Awrk(Asize)
!
!  Local variable declarations.
!
      integer :: LB, Lstr, MyError, Serror
      integer :: i, np, rank, request

      integer :: my_bounds(2)
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(2,0:NtileI(ng)*NtileJ(ng)-1) :: Abounds

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Collect data from all nodes.
!-----------------------------------------------------------------------
!
!  Collect data lower and upper bound dimensions.
!
      my_bounds(1)=Mstr
      my_bounds(2)=Mend
      CALL mpi_allgather (my_bounds, 2, MPI_INTEGER, Abounds, 2,        &
     &                    MPI_INTEGER, OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 10     FORMAT (/,' MP_GATHER_STATE - error during ',a,                 &
     &          ' call, Node = ',i3.3,' Error = ',i3,/,13x,a)
        exit_flag=2
        RETURN
      END IF
!
!  If master node, loop over other nodes and receive the data.
!
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          np=Abounds(2,rank)-Abounds(1,rank)+1
          LB=Abounds(1,rank)
          CALL mpi_irecv (Awrk(LB:), np, MP_FLOAT, rank, rank+5,        &
     &                    OCN_COMM_WORLD, request, MyError)
          CALL mpi_wait (request, status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
        END DO
!
!  Load master node contribution.
!
        DO i=Mstr,Mend
          Awrk(i)=A(i)
        END DO
!
!  Otherwise, send data to master node.
!
      ELSE
        np=Mend-Mstr+1
        CALL mpi_isend (A(Mstr:), np, MP_FLOAT, MyMaster, MyRank+5,     &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast collected data to all nodes.
!
      CALL mpi_bcast (Awrk, Asize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,  &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

      RETURN
      END SUBROUTINE mp_gather_state
!
      FUNCTION mp_ncread1d (ng, model, ncid, ncvname, ncname, ncrec,    &
     &                      LB1, UB1, Ascale, A)
!
!***********************************************************************
!                                                                      !
!  This function reads floating point 1D state array from specified    !
!  NetCDF file and scatters it to the other nodes.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvname      NetCDF variable name.                               !
!     ncname       NetCDF file name.                                   !
!     ncrec        NetCDF record index to write. If negative, it       !
!                    assumes that the variable is recordless.          !
!     LB1          First-dimension Lower bound.                        !
!     UB1          First-dimension Upper bound.                        !
!     Ascale       Factor to scale field after reading (real).         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A            Field to read in (real).                            !
!     mp_ncread1d  Error flag (integer).                               !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1

      real(r8), intent(in) :: Ascale

      real(r8), intent(out) :: A(LB1:UB1)

      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: ibuffer(2), my_bounds(2), start(1), total(1)

      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: Asize

      integer :: mp_ncread1d

      real(r8), allocatable :: Asend(:)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 47)
# endif
!
!-----------------------------------------------------------------------
!  Read requested NetCDF file and scatter it to all nodes.
!-----------------------------------------------------------------------
!
      mp_ncread1d=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      CALL mpi_allgather (my_bounds, 2, MPI_INTEGER,                    &
     &                    Asize, 2, MPI_INTEGER,                        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If not master node, receive data from master node.
!
      IF (MyRank.ne.MyMaster) THEN
        np=UB1-LB1+1
        CALL mpi_irecv (A(LB1:), np, MP_FLOAT, MyMaster, MyRank+5,      &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_IRECV', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
!
!  Scale recieved (read) data.
!
        DO i=LB1,UB1
          A(i)=A(i)*Ascale
        END DO
!
!  Otherwise, if master node allocate the send buffer.
!
      ELSE
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=Asize(2,rank)-Asize(1,rank)+1
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Asend)) THEN
          allocate (Asend(Npts))
        END IF
!
!  If master node, loop over all nodes and read buffers to send.
!
        mp_ncread1d=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (mp_ncread1d.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
          exit_flag=2
          ioerror=mp_ncread1d
        END IF
        IF (exit_flag.eq.NoError) THEN
          DO rank=0,NtileI(ng)*NtileJ(ng)-1
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            mp_ncread1d=nf90_get_var(ncid, varid,  Asend, start, total)
            IF (mp_ncread1d.ne.nf90_noerr) THEN
              WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
              exit_flag=2
              ioerror=mp_ncread1d
              EXIT
            END IF
!
!  Send buffer to all nodes, except itself.
!
            IF (rank.eq.MyMaster) THEN
              np=0
              DO i=LB1,UB1
                np=np+1
                A(i)=Asend(np)*Ascale
              END DO
            ELSE
              np=Asize(2,rank)-Asize(1,rank)+1
              CALL mpi_isend (Asend, np, MP_FLOAT, rank, rank+5,        &
     &                        OCN_COMM_WORLD, request, MyError)
              CALL mpi_wait (request, status, MyError)
              IF (MyError.ne.MPI_SUCCESS) THEN
                CALL mpi_error_string (MyError, string, Lstr, Serror)
                Lstr=LEN_TRIM(string)
                WRITE (stdout,10) 'MPI_ISEND', rank, MyError,           &
     &                            string(1:Lstr)
                exit_flag=2
                RETURN
              END IF
            END IF
          END DO
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Deallocate send buffer.
!
      IF (allocated(Asend).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Asend)
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 47)
# endif

 10   FORMAT (/,' MP_NCREAD1D- error during ',a,' call, Node = ',       &
     &          i3.3,' Error = ',i3,/,13x,a)
 20   FORMAT (/,' MP_NCREAD1D - error while inquiring ID for',          &
     &          ' variable: ',a,/,13x,'in file: ',a)
 30   FORMAT (/,' MP_NCREAD1D - error while reading variable: ',        &
     &        a,/,13x,'in file: ',a)

      RETURN
      END FUNCTION mp_ncread1d

      FUNCTION mp_ncread2d (ng, model, ncid, ncvname, ncname, ncrec,    &
     &                      LB1, UB1, LB2, UB2, Ascale, A)
!
!***********************************************************************
!                                                                      !
!  This function reads floating point 2D state array from specified    !
!  NetCDF file and scatters it to the other nodes.                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvname      NetCDF variable name.                               !
!     ncname       NetCDF file name.                                   !
!     ncrec        NetCDF record index to write. If negative, it       !
!                    assumes that the variable is recordless.          !
!     LB1          First-dimension Lower bound.                        !
!     UB1          First-dimension Upper bound.                        !
!     LB2          Second-dimension Lower bound.                       !
!     UB2          Second-dimension Upper bound.                       !
!     Ascale       Factor to scale field after reading (real).         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A            Field to read in (real).                            !
!     mp_ncread2d  Error flag (integer).                               !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1, LB2, UB2

      real(r8), intent(in) :: Ascale

      real(r8), intent(out) :: A(LB1:UB1,LB2:UB2)

      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: ibuffer(2), my_bounds(4), start(2), total(2)

      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: Asize

      integer :: mp_ncread2d

      real(r8), allocatable :: Asend(:)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 47)
# endif
!
!-----------------------------------------------------------------------
!  Read requested NetCDF file and scatter it to all nodes.
!-----------------------------------------------------------------------
!
      mp_ncread2d=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      my_bounds(3)=LB2
      my_bounds(4)=UB2
      CALL mpi_allgather (my_bounds, 4, MPI_INTEGER,                    &
     &                    Asize, 4, MPI_INTEGER,                        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If not master node, receive data from master node.
!
      IF (MyRank.ne.MyMaster) THEN
        np=(UB1-LB1+1)*(UB2-LB2+1)
        CALL mpi_irecv (A(LB1,LB2), np, MP_FLOAT, MyMaster, MyRank+5,   &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_IRECV', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
!
!  Scale recieved (read) data.
!
        DO j=LB2,UB2
          DO i=LB1,UB1
            A(i,j)=A(i,j)*Ascale
          END DO
        END DO
!
!  Otherwise, if master node allocate the send buffer.
!
      ELSE
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=(Asize(2,rank)-Asize(1,rank)+1)*                           &
     &       (Asize(4,rank)-Asize(3,rank)+1)
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Asend)) THEN
          allocate (Asend(Npts))
        END IF
!
!  If master node, loop over all nodes and read buffers to send.
!
        mp_ncread2d=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (mp_ncread2d.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
          exit_flag=2
          ioerror=mp_ncread2d
        END IF
        IF (exit_flag.eq.NoError) THEN
          DO rank=0,NtileI(ng)*NtileJ(ng)-1
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            start(2)=Asize(3,rank)
            total(2)=Asize(4,rank)-Asize(3,rank)+1
            mp_ncread2d=nf90_get_var(ncid, varid,  Asend, start, total)
            IF (mp_ncread2d.ne.nf90_noerr) THEN
              WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
              exit_flag=2
              ioerror=mp_ncread2d
              EXIT
            END IF
!
!  Send buffer to all nodes, except itself.
!
            IF (rank.eq.MyMaster) THEN
              np=0
              DO j=LB2,UB2
                DO i=LB1,UB1
                  np=np+1
                  A(i,j)=Asend(np)*Ascale
                END DO
              END DO
            ELSE
              np=(Asize(2,rank)-Asize(1,rank)+1)*                       &
     &           (Asize(4,rank)-Asize(3,rank)+1)
              CALL mpi_isend (Asend, np, MP_FLOAT, rank, rank+5,        &
     &                        OCN_COMM_WORLD, request, MyError)
              CALL mpi_wait (request, status, MyError)
              IF (MyError.ne.MPI_SUCCESS) THEN
                CALL mpi_error_string (MyError, string, Lstr, Serror)
                Lstr=LEN_TRIM(string)
                WRITE (stdout,10) 'MPI_ISEND', rank, MyError,           &
     &                            string(1:Lstr)
                exit_flag=2
                RETURN
              END IF
            END IF
          END DO
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Deallocate send buffer.
!
      IF (allocated(Asend).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Asend)
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 47)
# endif

 10   FORMAT (/,' MP_NCREAD2D - error during ',a,' call, Node = ',      &
     &          i3.3,' Error = ',i3,/,13x,a)
 20   FORMAT (/,' MP_NCREAD2D - error while inquiring ID for',          &
     &          ' variable: ',a,/,13x,'in file: ',a)
 30   FORMAT (/,' MP_NCREAD2D - error while reading variable: ',        &
     &        a,/,13x,'in file: ',a)

      RETURN
      END FUNCTION mp_ncread2d

      FUNCTION mp_ncwrite1d (ng, model, ncid, ncvname, ncname, ncrec,   &
     &                       LB1, UB1, Ascale, A)
!
!***********************************************************************
!                                                                      !
!  This function collects floating point 1D state array data from the  !
!  other nodes and writes it into specified NetCDF file.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number.                                !
!     model         Calling model identifier.                          !
!     ncid          NetCDF file ID.                                    !
!     ncvname       NetCDF variable name.                              !
!     ncname        NetCDF file name.                                  !
!     ncrec         NetCDF record index to write. If negative, it      !
!                     assumes that the variable is recordless.         !
!     LB1           First-dimension Lower bound.                       !
!     UB1           First-dimension Upper bound.                       !
!     Ascale        Factor to scale field before writing (real).       !
!     A             Field to write out (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     mp_ncwrite1d  Error flag (integer).                              !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1

      real(r8), intent(in) :: Ascale

      real(r8), intent(in) :: A(LB1:UB1)

      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: ibuffer(2), my_bounds(2), start(1), total(1)

      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(2,0:NtileI(ng)*NtileJ(ng)-1) :: Asize

      integer :: mp_ncwrite1d

      real(r8), allocatable :: Arecv(:)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Collect and write data into requested NetCDF file.
!-----------------------------------------------------------------------
!
      mp_ncwrite1d=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      CALL mpi_allgather (my_bounds, 2, MPI_INTEGER,                    &
     &                    Asize, 2, MPI_INTEGER,                        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If master node, allocate the receive buffer.
!
      IF (MyRank.eq.MyMaster) THEN
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=(Asize(2,rank)-Asize(1,rank)+1)
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  Write out master node contribution.
!
        start(1)=LB1
        total(1)=UB1-LB1+1
        np=0
        DO i=LB1,UB1
          np=np+1
          Arecv(np)=A(i)
        END DO
        mp_ncwrite1d=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (mp_ncwrite1d.eq.nf90_noerr) THEN
          mp_ncwrite1d=nf90_put_var(ncid, varid, Arecv, start, total)
          IF (mp_ncwrite1d.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
            exit_flag=3
            ioerror=mp_ncwrite1d
          END IF
        ELSE
          WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
          exit_flag=3
          ioerror=mp_ncwrite1d
        END IF
!
!  If master node, loop over other nodes and receive the data.
!
        IF (exit_flag.eq.NoError) THEN
          DO rank=1,NtileI(ng)*NtileJ(ng)-1
            np=Asize(2,rank)-Asize(1,rank)+1
            CALL mpi_irecv (Arecv, np, MP_FLOAT, rank, rank+5,          &
     &                      OCN_COMM_WORLD, request, MyError)
            CALL mpi_wait (request, status, MyError)
            IF (MyError.ne.MPI_SUCCESS) THEN
              CALL mpi_error_string (MyError, string, Lstr, Serror)
              Lstr=LEN_TRIM(string)
              WRITE (stdout,10) 'MPI_IRECV', rank, MyError,             &
     &                          string(1:Lstr)
              exit_flag=3
              RETURN
            END IF
!
!  Write out data into NetCDF file.
!
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            DO i=1,np
              Arecv(i)=Arecv(i)*Ascale
            END DO
            mp_ncwrite1d=nf90_put_var(ncid, varid, Arecv, start, total)
            IF (mp_ncwrite1d.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
              exit_flag=3
              ioerror=mp_ncwrite1d
              EXIT
            END IF
          END DO
        END IF
!
!  Otherwise, send data to master node.
!
      ELSE
        np=UB1-LB1+1
        CALL mpi_isend (A(LB1:), np, MP_FLOAT, MyMaster, MyRank+5,      &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Deallocate receive buffer.
!
      IF (allocated(Arecv).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Arecv)
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

 10   FORMAT (/,' MP_NCWRITE1D - error during ',a,' call, Node = ',     &
     &          i3.3,' Error = ',i3,/,13x,a)
 20   FORMAT (/,' MP_NCWRITE1D - error while writing variable: ',       &
     &        a,/,13x,'in file: ',a)
 30   FORMAT (/,' MP_NCWRITE1D - error while inquiring ID for',         &
     &        ' variable: ',a,/,13x,'in file: ',a)

      RETURN
      END FUNCTION mp_ncwrite1d

      FUNCTION mp_ncwrite2d (ng, model, ncid, ncvname, ncname, ncrec,   &
     &                       LB1, UB1, LB2, UB2, Ascale, A)
!
!***********************************************************************
!                                                                      !
!  This function collects floating point 2D state array data from the  !
!  other nodes and writes it into specified NetCDF file.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng            Nested grid number.                                !
!     model         Calling model identifier.                          !
!     ncid          NetCDF file ID.                                    !
!     ncvname       NetCDF variable name.                              !
!     ncname        NetCDF file name.                                  !
!     ncrec         NetCDF record index to write. If negative, it      !
!                     assumes that the variable is recordless.         !
!     LB1           First-dimension Lower bound.                       !
!     UB1           First-dimension Upper bound.                       !
!     LB2           Second-dimension Lower bound.                      !
!     UB2           Second-dimension Upper bound.                      !
!     Ascale        Factor to scale field before writing (real).       !
!     A             Field to write out (real).                         !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     mp_ncwrite2d  Error flag (integer).                              !
!                                                                      !
!  Note: We cannot include "USE mod_netcdf" here because of cyclic     !
!        dependency. Instead we need original NetCDF library module    !
!        "USE netcdf".                                                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncrec
      integer, intent(in) :: LB1, UB1, LB2, UB2

      real(r8), intent(in) :: Ascale

      real(r8), intent(in) :: A(LB1:UB1,LB2:UB2)

      character (len=*), intent(in) :: ncvname
      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Npts, Serror
      integer :: i, j, np, rank, request, varid
      integer :: ibuffer(2), my_bounds(4), start(2), total(2)

      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: Asize

      integer :: mp_ncwrite2d

      real(r8), allocatable :: Arecv(:)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Collect and write data into requested NetCDF file.
!-----------------------------------------------------------------------
!
      mp_ncwrite2d=nf90_noerr
!
!  Collect data lower and upper bounds dimensions.
!
      my_bounds(1)=LB1
      my_bounds(2)=UB1
      my_bounds(3)=LB2
      my_bounds(4)=UB2
      CALL mpi_allgather (my_bounds, 4, MPI_INTEGER,                    &
     &                    Asize, 4, MPI_INTEGER,                        &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  If master node, allocate the receive buffer.
!
      IF (MyRank.eq.MyMaster) THEN
        Npts=0
        DO rank=0,NtileI(ng)*NtileJ(ng)-1
          np=(Asize(2,rank)-Asize(1,rank)+1)*                           &
     &       (Asize(4,rank)-Asize(3,rank)+1)
          Npts=MAX(Npts, np)
        END DO
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Npts))
        END IF
!
!  Write out master node contribution.
!
        start(1)=LB1
        total(1)=UB1-LB1+1
        start(2)=LB2
        total(2)=UB2-LB2+1
        np=0
        DO j=LB2,UB2
          DO i=LB1,UB1
            np=np+1
            Arecv(np)=A(i,j)
          END DO
        END DO
        mp_ncwrite2d=nf90_inq_varid(ncid, TRIM(ncvname), varid)
        IF (mp_ncwrite2d.eq.nf90_noerr) THEN
          mp_ncwrite2d=nf90_put_var(ncid, varid, Arecv, start, total)
          IF (mp_ncwrite2d.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
            exit_flag=3
            ioerror=mp_ncwrite2d
          END IF
        ELSE
          WRITE (stdout,30) TRIM(ncvname), TRIM(ncname)
          exit_flag=3
          ioerror=mp_ncwrite2d
        END IF
!
!  If master node, loop over other nodes and receive the data.
!
        IF (exit_flag.eq.NoError) THEN
          DO rank=1,NtileI(ng)*NtileJ(ng)-1
            np=(Asize(2,rank)-Asize(1,rank)+1)*                         &
     &         (Asize(4,rank)-Asize(3,rank)+1)
            CALL mpi_irecv (Arecv, np, MP_FLOAT, rank, rank+5,          &
     &                      OCN_COMM_WORLD, request, MyError)
            CALL mpi_wait (request, status, MyError)
            IF (MyError.ne.MPI_SUCCESS) THEN
              CALL mpi_error_string (MyError, string, Lstr, Serror)
              Lstr=LEN_TRIM(string)
              WRITE (stdout,10) 'MPI_IRECV', rank, MyError,             &
     &                          string(1:Lstr)
              exit_flag=3
              RETURN
            END IF
!
!  Write out data into NetCDF file.
!
            start(1)=Asize(1,rank)
            total(1)=Asize(2,rank)-Asize(1,rank)+1
            start(2)=Asize(3,rank)
            total(2)=Asize(4,rank)-Asize(3,rank)+1
            DO i=1,np
              Arecv(i)=Arecv(i)*Ascale
            END DO
            mp_ncwrite2d=nf90_put_var(ncid, varid, Arecv, start, total)
            IF (mp_ncwrite2d.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(ncvname), TRIM(ncname)
              exit_flag=3
              ioerror=mp_ncwrite2d
              EXIT
            END IF
          END DO
        END IF
!
!  Otherwise, send data to master node.
!
      ELSE
        np=(UB1-LB1+1)*(UB2-LB2+1)
        CALL mpi_isend (A(LB1:,LB2:), np, MP_FLOAT, MyMaster, MyRank+5, &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast error flags to all nodes.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
!
!  Deallocate receive buffer.
!
      IF (allocated(Arecv).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Arecv)
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

 10   FORMAT (/,' MP_NCWRITE2D - error during ',a,' call, Node = ',     &
     &          i3.3,' Error = ',i3,/,13x,a)
 20   FORMAT (/,' MP_NCWRITE2D - error while writing variable: ',       &
     &        a,/,13x,'in file: ',a)
 30   FORMAT (/,' MP_NCWRITE2D - error while inquiring ID for',         &
     &        ' variable: ',a,/,13x,'in file: ',a)

      RETURN
      END FUNCTION mp_ncwrite2d

      SUBROUTINE mp_reduce_0d (ng, model, Asize, A, op_handle)
!
!***********************************************************************
!                                                                      !
!  This routine collects and reduces requested variables from all      !
!  nodes in the group.  Then,  it broadcasts reduced variables to      !
!  all nodes in the group.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Asize      Number of scalar variables to reduce.                 !
!     A          Vector of scalar variables to reduce.                 !
!     op_handle  Reduction operation handle (string).  The following   !
!                  reduction operations are supported:                 !
!                  'MIN', 'MAX', 'SUM'                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Vector of reduced scalar variables.                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Asize

      character (len=*), intent(in) :: op_handle

      real(r8), intent(inout) :: A
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Serror
      integer :: handle, i, rank, request

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus

      real(r8) :: Areduce
      real(r8) :: Asend
      real(r8), dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 45)
# endif
!
!-----------------------------------------------------------------------
!  Collect and reduce requested scalar variables.
!-----------------------------------------------------------------------
!
!  Pack data to reduce.
!
      Asend=A
!
!  Collect and reduce.
!
# if defined REDUCE_ALLREDUCE
      IF (op_handle(1:3).eq.'MIN') THEN
        handle=MPI_MIN
      ELSE IF (op_handle(1:3).eq.'MAX') THEN
        handle=MPI_MAX
      ELSE IF (op_handle(1:3).eq.'SUM') THEN
        handle=MPI_SUM
      END IF
      CALL mpi_allreduce (Asend, Areduce, 1, MP_FLOAT, handle,          &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# elif defined REDUCE_ALLGATHER
      CALL mpi_allgather (Asend, 1, MP_FLOAT,                           &
     &                    Arecv, 1, MP_FLOAT,                           &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
      Areduce=Arecv(0)
      DO rank=1,NtileI(ng)*NtileJ(ng)-1
        IF (op_handle(1:3).eq.'MIN') THEN
          Areduce=MIN(Areduce,Arecv(rank))
        ELSE IF (op_handle(1:3).eq.'MAX') THEN
          Areduce=MAX(Areduce,Arecv(rank))
        ELSE IF (op_handle(1:3).eq.'SUM') THEN
          Areduce=Areduce+Arecv(rank)
        END IF
      END DO
# else
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(rank), 1, MP_FLOAT, rank,               &
     &                    rank+500, OCN_COMM_WORLD, Rrequest(rank),     &
     &                    MyError)
        END DO
        Areduce=Asend
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, Rerror, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
          IF (op_handle(1:3).eq.'MIN') THEN
            Areduce=MIN(Areduce,Arecv(rank))
          ELSE IF (op_handle(1:3).eq.'MAX') THEN
            Areduce=MAX(Areduce,Arecv(rank))
          ELSE IF (op_handle(1:3).eq.'SUM') THEN
            Areduce=Areduce+Arecv(rank)
          END IF
        END DO
      ELSE
        CALL mpi_isend (Asend, 1, MP_FLOAT, MyMaster, MyRank+500,       &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, Sstatus, MyError)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast reduced variables from process to all processes in the
!  group.
!
      CALL mpi_bcast (Areduce, 1, MP_FLOAT, MyMaster,                   &
     &                OCN_COMM_WORLD, MyError)
      IF (Serror.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# endif
 10   FORMAT (/,' MP_REDUCE_0D - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,16x,a)
!
!  Unpack.
!
      A=Areduce
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 45)
# endif

      RETURN
      END SUBROUTINE mp_reduce_0d
!
      SUBROUTINE mp_reduce_1d (ng, model, Asize, A, op_handle)
!
!***********************************************************************
!                                                                      !
!  This routine collects and reduces requested variables from all      !
!  nodes in the group.  Then,  it broadcasts reduced variables to      !
!  all nodes in the group.                                             !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Asize      Number of scalar variables to reduce.                 !
!     A          Vector of scalar variables to reduce.                 !
!     op_handle  Reduction operation handle (string).  The following   !
!                  reduction operations are supported:                 !
!                  'MIN', 'MAX', 'SUM'                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Vector of reduced scalar variables.                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Asize

      character (len=*), intent(in) :: op_handle(Asize)

      real(r8), intent(inout) :: A(Asize)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Serror
      integer :: handle, i, rank, request

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: Rstatus
      integer, dimension(MPI_STATUS_SIZE) :: Sstatus

      real(r8), dimension(Asize,0:NtileI(ng)*NtileJ(ng)-1) :: Arecv
      real(r8), dimension(Asize) :: Areduce
      real(r8), dimension(Asize) :: Asend

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 45)
# endif
!
!-----------------------------------------------------------------------
!  Collect and reduce requested scalar variables.
!-----------------------------------------------------------------------
!
!  Pack data to reduce.
!
      DO i=1,Asize
        Asend(i)=A(i)
      END DO
!
!  Collect and reduce.
!
# if defined REDUCE_ALLREDUCE
      DO i=1,Asize
        IF (op_handle(i)(1:3).eq.'MIN') THEN
          handle=MPI_MIN
        ELSE IF (op_handle(i)(1:3).eq.'MAX') THEN
          handle=MPI_MAX
        ELSE IF (op_handle(i)(1:3).eq.'SUM') THEN
          handle=MPI_SUM
        END IF
        CALL mpi_allreduce (Asend(i), Areduce(i), 1, MP_FLOAT, handle,  &
     &                      OCN_COMM_WORLD, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,           &
     &                      string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END DO
# elif defined REDUCE_ALLGATHER
      CALL mpi_allgather (Asend, Asize, MP_FLOAT,                       &
     &                    Arecv, Asize, MP_FLOAT,                       &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
      DO i=1,Asize
        Areduce(i)=Arecv(i,0)
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          IF (op_handle(i)(1:3).eq.'MIN') THEN
            Areduce(i)=MIN(Areduce(i),Arecv(i,rank))
          ELSE IF (op_handle(i)(1:3).eq.'MAX') THEN
            Areduce(i)=MAX(Areduce(i),Arecv(i,rank))
          ELSE IF (op_handle(i)(1:3).eq.'SUM') THEN
            Areduce(i)=Areduce(i)+Arecv(i,rank)
          END IF
        END DO
      END DO
# else
      IF (MyRank.eq.MyMaster) THEN
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv(1,rank), Asize, MP_FLOAT, rank,         &
     &                    rank+500, OCN_COMM_WORLD, Rrequest(rank),     &
     &                    MyError)
        END DO
        DO i=1,Asize
          Areduce(i)=Asend(i)
        END DO
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_wait (Rrequest(rank), Rstatus, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, Rerror, string(1:Lstr)
            exit_flag=2
            RETURN
          END IF
          DO i=1,Asize
            IF (op_handle(i)(1:3).eq.'MIN') THEN
              Areduce(i)=MIN(Areduce(i),Arecv(i,rank))
            ELSE IF (op_handle(i)(1:3).eq.'MAX') THEN
              Areduce(i)=MAX(Areduce(i),Arecv(i,rank))
            ELSE IF (op_handle(i)(1:3).eq.'SUM') THEN
              Areduce(i)=Areduce(i)+Arecv(i,rank)
            END IF
          END DO
        END DO
      ELSE
        CALL mpi_isend (Asend, Asize, MP_FLOAT, MyMaster, MyRank+500,   &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, Sstatus, MyError)
        IF (Serror.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, Serror, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast reduced variables from process to all processes in the
!  group.
!
      CALL mpi_bcast (Areduce, Asize, MP_FLOAT, MyMaster,               &
     &                OCN_COMM_WORLD, MyError)
      IF (Serror.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
# endif
 10   FORMAT (/,' MP_REDUCE_1D - error during ',a,' call, Node = ',     &
     &        i3.3,' Error = ',i3,/,16x,a)
!
!  Unpack.
!
      DO i=1,Asize
        A(i)=Areduce(i)
      END DO
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 45)
# endif

      RETURN
      END SUBROUTINE mp_reduce_1d
!
      SUBROUTINE mp_reduce2 (ng, model, Isize, Jsize, A, op_handle)
!
!***********************************************************************
!                                                                      !
!  This routine computes the global minimum/maximum and its associated !
!  qualifiers like: location and/or other scalar components.  Then, it !
!  it broadcasts reduced variables to all nodes in the group.          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Isize      Size of I-dimension: the minimum/maximum to reduce    !
!                  is in location A(1,:) and qualifiers A(2:Isize,:).  !
!     Jsize      Size of J-dimension: number of different sets of      !
!                  minimum and/or maximum to process.                  !
!     A          Matrix of variables and qualifiers to reduce.         !
!     op_handle  Reduction operation handle (string) of size Jsize.    !
!                  The following  reduction operations are supported:  !
!                  'MINLOC', 'MAXLOC'                                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Matrix of reduced variables and qualifiers.           !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Isize, Jsize

      character (len=*), intent(in) :: op_handle(Jsize)

      real(r8), intent(inout) :: A(Isize,Jsize)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Serror
      integer :: handle, i, j

      real(r8), dimension(2,Isize) :: Areduce
      real(r8), dimension(2,Isize) :: Asend

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 45)
# endif
!
!-----------------------------------------------------------------------
!  Reduce requested variables and qualifiers.
!-----------------------------------------------------------------------
!
!  Pack and reduce.
!
      DO j=1,Jsize
        DO i=1,Isize
          Asend(1,i)=A(1,j)
          Asend(2,i)=A(i,j)
        END DO
        IF (op_handle(j)(1:6).eq.'MINLOC') THEN
          handle=MPI_MINLOC
        ELSE IF (op_handle(j)(1:6).eq.'MAXLOC') THEN
          handle=MPI_MAXLOC
        END IF
        CALL mpi_allreduce (Asend, Areduce, Isize,                      &
# ifdef DOUBLE_PRECISION
     &                      MPI_2DOUBLE_PRECISION,                      &
# else
     &                      MPI_2REAL,                                  &
# endif
     &                      handle, OCN_COMM_WORLD, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ALLREDUCE', MyRank, MyError,           &
     &                      string(1:Lstr)
 10       FORMAT (/,' MP_REDUCE2 - error during ',a,' call, Node = ',   &
     &            i3.3,' Error = ',i3,/,16x,a)
          exit_flag=2
          RETURN
        END IF
!
!  Unpack.
!
        A(1,j)=Areduce(1,1)
        DO i=2,Isize
          A(i,j)=Areduce(2,i)
        END DO
      END DO

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 45)
# endif

      RETURN
      END SUBROUTINE mp_reduce2
!
      SUBROUTINE mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,           &
     &                         Nghost, gtype, Amin, Amax,               &
# if defined READ_WATER && defined MASKING
     &                         NWpts, IJ_water,                         &
# endif
     &                         Npts, A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts input global data, packed as 1D real array, !
!  to each spawned MPI node.  Because this routine is also used by the !
!  adjoint model,  the ghost-points in the halo region are NOT updated !
!  in the ouput tile array (Awrk).  It is used by the  master node  to !
!  scatter input global data to each tiled node.                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     Nghost     Number of ghost-points in the halo region.            !
!     gtype      C-grid type. If negative and Land-Sea mask is         !
!                  available, only water-points are processed.         !
!     Amin       Input array minimum value.                            !
!     Amax       Input array maximum value.                            !
!     NWpts      Number of water points.                               !
!     IJ_water   IJ-indices for water points.                          !
!     Npts       Number of points to processes in A.                   !
!     A          Input global data from each node packed into 1D array !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!     Npts       Number of points to processes in A.                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Awrk       2D tiled, floating-point array.                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Nghost, gtype, Npts

# if defined READ_WATER && defined MASKING
      integer, intent(in) :: NWpts
      integer, intent(in) :: IJ_water(NWpts)
# endif
      real(r8), intent(inout) :: Amin, Amax
      real(r8), intent(inout) :: A(Npts+2)
      real(r8), intent(out) :: Awrk(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer :: Io, Ie, Jo, Je, Ioff, Joff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, IJlen
      integer :: Lstr, MyError, MySize, MyType, Serror, ghost
      integer :: i, ic, ij, j, jc, mc, nc

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 47)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT

      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      IJlen=Ilen*Jlen
!
!  Set physical, non-overlapping (Nghost=0) or overlapping (Nghost>0)
!  ranges according to tile rank.
!
      IF (Nghost.eq.0) THEN
        ghost=0                                   ! non-overlapping
      ELSE
        ghost=1                                   ! overlapping
      END IF

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Imin=BOUNDS(ng) % Imin(1,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(1,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(1,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(1,ghost,MyRank)
        CASE (r2dvar, r3dvar)
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
        CASE (u2dvar, u3dvar)
          Imin=BOUNDS(ng) % Imin(3,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(3,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(3,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(3,ghost,MyRank)
        CASE (v2dvar, v3dvar)
          Imin=BOUNDS(ng) % Imin(4,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(4,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(4,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(4,ghost,MyRank)
        CASE DEFAULT                              ! RHO-points
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
      END SELECT
!
!  Size of broadcast buffer.
!
      IF (gtype.gt.0) THEN
        MySize=IJlen
      ELSE
        MySize=Npts
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Scatter requested array data.
!-----------------------------------------------------------------------
!
!  If master processor, append minimum and maximum values to the end of
!  the buffer.
!
      IF (MyRank.eq.MyMaster) Then
        A(MySize+1)=Amin
        A(MySize+2)=Amax
      END IF
      MySize=MySize+2
!
!  Broadcast data to all processes in the group, itself included.
!
      CALL mpi_bcast (A, MySize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
         Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_SCATTER2D - error during ',a,' call, Node = ',   &
     &          i3.3, ' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
!
!  If water points only, fill land points.
!
      IF (gtype.gt.0) THEN
        DO nc=1,MySize-2
          Arecv(nc)=A(nc)
        END DO
# if defined READ_WATER && defined MASKING
      ELSE
        ij=0
        mc=0
        nc=0
        DO j=Jo,Je
          jc=(j-Joff)*Ilen
          DO i=Io,Ie
            ij=ij+1
            ic=i+Ioff+jc
            IF (IJ_water(mc+1).eq.ij) THEN
              mc=mc+1
              nc=nc+1
              Arecv(ic)=A(nc)
            ELSE
              Arecv(ic)=0.0_r8
            ENDIF
          END DO
        END DO
# endif
      END IF
!
!  Unpack data buffer.
!
      DO j=Jmin,Jmax
        jc=(j-Joff)*Ilen
        DO i=Imin,Imax
          ic=i+Ioff+jc
          Awrk(i,j)=Arecv(ic)
        END DO
      END DO
      Amin=A(MySize-1)
      Amax=A(MySize)
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 47)
# endif

      RETURN
      END SUBROUTINE mp_scatter2d
!
      SUBROUTINE mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, LBk, UBk, &
     &                         Nghost, gtype, Amin, Amax,               &
# if defined READ_WATER && defined MASKING
     &                         NWpts, IJ_water,                         &
# endif
     &                         Npts, A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine broadcasts input global data, packed as 1D real array, !
!  to each spawned MPI node.  Because this routine is also used by the !
!  adjoint model,  the ghost-points in the halo region are NOT updated !
!  in the ouput tile array (Awrk).  It is used by the  master node  to !
!  scatter input global data to each tiled node.                       !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Nghost     Number of ghost-points in the halo region.            !
!     gtype      C-grid type. If negative and Land-Sea mask is         !
!                  available, only water-points are processed.         !
!     Amin       Input array minimum value.                            !
!     Amax       Input array maximum value.                            !
!     NWpts      Number of water points.                               !
!     IJ_water   IJ-indices for water points.                          !
!     Npts       Number of points to processes in A.                   !
!     A          Input global data from each node packed into 1D array !
!                  in column-major order. That is, in the same way     !
!                  that Fortran multi-dimensional arrays are stored    !
!                  in memory.                                          !
!     Npts       Number of points to processes in A.                   !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Awrk       3D tiled, floating-point array.                       !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Nghost, gtype, Npts

# if defined READ_WATER && defined MASKING
      integer, intent(in) :: NWpts
      integer, intent(in) :: IJ_water(NWpts)
# endif
      real(r8), intent(inout) :: Amin, Amax
      real(r8), intent(inout) :: A(Npts+2)
      real(r8), intent(out) :: Awrk(LBi:UBi,LBj:UBj,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Io, Ie, Jo, Je, Ioff, Joff, Koff
      integer :: Imin, Imax, Jmin, Jmax
      integer :: Ilen, Jlen, Klen, IJlen
      integer :: Lstr, MyError, MySize, MyType, Serror, ghost
      integer :: i, ic, ij, j, jc, k, kc, mc, nc

      real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*(UBk-LBk+1)) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 47)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Set full grid first and last point according to staggered C-grid
!  classification. Notice that the offsets are for the private array
!  counter.
!
      MyType=ABS(gtype)

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Io=IOBOUNDS(ng) % ILB_psi
          Ie=IOBOUNDS(ng) % IUB_psi
          Jo=IOBOUNDS(ng) % JLB_psi
          Je=IOBOUNDS(ng) % JUB_psi
          Ioff=0
          Joff=1
        CASE (r2dvar, r3dvar)
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
        CASE (u2dvar, u3dvar)
          Io=IOBOUNDS(ng) % ILB_u
          Ie=IOBOUNDS(ng) % IUB_u
          Jo=IOBOUNDS(ng) % JLB_u
          Je=IOBOUNDS(ng) % JUB_u
          Ioff=0
          Joff=0
        CASE (v2dvar, v3dvar)
          Io=IOBOUNDS(ng) % ILB_v
          Ie=IOBOUNDS(ng) % IUB_v
          Jo=IOBOUNDS(ng) % JLB_v
          Je=IOBOUNDS(ng) % JUB_v
          Ioff=1
          Joff=1
        CASE DEFAULT                              ! RHO-points
          Io=IOBOUNDS(ng) % ILB_rho
          Ie=IOBOUNDS(ng) % IUB_rho
          Jo=IOBOUNDS(ng) % JLB_rho
          Je=IOBOUNDS(ng) % JUB_rho
          Ioff=1
          Joff=0
      END SELECT

      IF (LBk.eq.0) THEN
        Koff=0
      ELSE
        Koff=1
      END IF

      Ilen=Ie-Io+1
      Jlen=Je-Jo+1
      Klen=UBk-LBk+1
      IJlen=Ilen*Jlen
!
!  Set physical, non-overlapping (Nghost=0) or overlapping (Nghost>0)
!  ranges according to tile rank.
!
      IF (Nghost.eq.0) THEN
        ghost=0                                   ! non-overlapping
      ELSE
        ghost=1                                   ! overlapping
      END IF

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          Imin=BOUNDS(ng) % Imin(1,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(1,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(1,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(1,ghost,MyRank)
        CASE (r2dvar, r3dvar)
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
        CASE (u2dvar, u3dvar)
          Imin=BOUNDS(ng) % Imin(3,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(3,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(3,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(3,ghost,MyRank)
        CASE (v2dvar, v3dvar)
          Imin=BOUNDS(ng) % Imin(4,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(4,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(4,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(4,ghost,MyRank)
        CASE DEFAULT                              ! RHO-points
          Imin=BOUNDS(ng) % Imin(2,ghost,MyRank)
          Imax=BOUNDS(ng) % Imax(2,ghost,MyRank)
          Jmin=BOUNDS(ng) % Jmin(2,ghost,MyRank)
          Jmax=BOUNDS(ng) % Jmax(2,ghost,MyRank)
      END SELECT
!
!  Size of broadcast buffer.
!
      IF (gtype.gt.0) THEN
        MySize=IJlen*Klen
      ELSE
        MySize=Npts
      END IF
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Scatter requested array data.
!-----------------------------------------------------------------------
!
!  If master processor, append minimum and maximum values to the end of
!  the buffer.
!
      IF (MyRank.eq.MyMaster) Then
        A(MySize+1)=Amin
        A(MySize+2)=Amax
      END IF
      MySize=MySize+2
!
!  Broadcast data to all processes in the group, itself included.
!
      CALL mpi_bcast (A, MySize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,    &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
         Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
 10     FORMAT (/,' MP_SCATTER3D - error during ',a,' call, Node = ',   &
     &          i3.3, ' Error = ',i3,/,15x,a)
        exit_flag=2
        RETURN
      END IF
!
!  If water points only, fill land points.
!
      IF (gtype.gt.0) THEN
        DO nc=1,MySize-2
          Arecv(nc)=A(nc)
        END DO
# if defined READ_WATER && defined MASKING
      ELSE
        nc=0
        DO k=LBk,UBk
          kc=(k-Koff)*IJlen
          ij=0
          mc=0
          DO j=Jo,Je
            jc=(j-Joff)*Ilen+kc
            DO i=Io,Ie
              ij=ij+1
              ic=i+Ioff+jc
              IF (IJ_water(mc+1).eq.ij) THEN
                mc=mc+1
                nc=nc+1
                Arecv(ic)=A(nc)
              ELSE
                Arecv(ic)=0.0_r8
              ENDIF
            END DO
          END DO
        END DO
# endif
      END IF
!
!  Unpack data buffer.
!
      DO k=LBk,UBk
        kc=(k-Koff)*IJlen
        DO j=Jmin,Jmax
          jc=(j-Joff)*Ilen+kc
          DO i=Imin,Imax
            ic=i+Ioff+jc
            Awrk(i,j,k)=Arecv(ic)
          END DO
        END DO
      END DO
      Amin=A(MySize-1)
      Amax=A(MySize)
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 47)
# endif

      RETURN
      END SUBROUTINE mp_scatter3d
!
      SUBROUTINE mp_scatter_state (ng, model, Mstr, Mend, Asize,        &
     &                             A, Awrk)
!
!***********************************************************************
!                                                                      !
!  This routine scatters (global to threaded) state data to all nodes  !
!  in the group. Before this can be done, the global data needs to be  !
!  collected from all the  nodes  by the master.  This is achieved by  !
!  summing the input values at each point.  This  routine  is used to  !
!  pack the state data for the GST analysis propagators.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     Mstr       Threaded array lower bound.                           !
!     Mend       Threaded array upper bound.                           !
!     Asize      Size of the .                                         !
!     A          Threaded 1D array process.                            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Collected data from all nodes.                        !
!     Awrk       Threaded block of data.                               !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend, Asize

      real(r8), intent(inout)  :: A(Asize)

      real(r8), intent(out) :: Awrk(Mstr:Mend)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, Serror
      integer :: i, rank, request

      integer, dimension(0:NtileI(ng)*NtileJ(ng)-1) :: Rrequest

      integer, dimension(MPI_STATUS_SIZE) :: status

      real(r8), allocatable :: Arecv(:)

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 47)
# endif
!
!-----------------------------------------------------------------------
!  Collect data blocks from all nodes and scatter the data to all nodes.
!-----------------------------------------------------------------------
!
!  All nodes have distinct pieces of the data and zero everywhere else.
!  So the strategy here is for the master node to receive the data from
!  the other nodes (excluding itself) and accumulate the sum at each
!  point. Then, the master node broadcast (itself included) its copy of
!  the accumlated data to other the nodes in the group. After this, each
!  node loads only the required block of the data into output array.
!
!  Notice that only the master node allocates the recieving buffer
!  (Arecv). It also receives only buffer at the time to avoid having
!  a very large communication array.  So here memory is more important
!  than time.
!
      IF (MyRank.eq.MyMaster) THEN
!
!  If master node, allocate and receive buffer.
!
        IF (.not.allocated(Arecv)) THEN
          allocate (Arecv(Asize))
        END IF
!
!  If master node, loop over other nodes to receive and accumulate the
!  data.
!
        DO rank=1,NtileI(ng)*NtileJ(ng)-1
          CALL mpi_irecv (Arecv, Asize, MP_FLOAT, rank, rank+5,         &
     &                    OCN_COMM_WORLD, Rrequest(rank), MyError)
          CALL mpi_wait (Rrequest(rank), status, MyError)
          IF (MyError.ne.MPI_SUCCESS) THEN
            CALL mpi_error_string (MyError, string, Lstr, Serror)
            Lstr=LEN_TRIM(string)
            WRITE (stdout,10) 'MPI_IRECV', rank, MyError, string(1:Lstr)
 10         FORMAT (/,' MP_SCATTER_STATE - error during ',a,            &
     &              ' call, Node = ', i3.3,' Error = ',i3,/,13x,a)
            exit_flag=2
            RETURN
          END IF
          DO i=1,Asize
            A(i)=A(i)+Arecv(i)
          END DO
        END DO
!
!  Otherwise, send data to master node.
!
      ELSE
        CALL mpi_isend (A, Asize, MP_FLOAT, MyMaster, MyRank+5,         &
     &                  OCN_COMM_WORLD, request, MyError)
        CALL mpi_wait (request, status, MyError)
        IF (MyError.ne.MPI_SUCCESS) THEN
          CALL mpi_error_string (MyError, string, Lstr, Serror)
          Lstr=LEN_TRIM(string)
          WRITE (stdout,10) 'MPI_ISEND', MyRank, MyError, string(1:Lstr)
          exit_flag=2
          RETURN
        END IF
      END IF
!
!  Broadcast accumulated (full) data to all nodes.
!
      CALL mpi_bcast (A, Asize, MP_FLOAT, MyMaster, OCN_COMM_WORLD,     &
     &                MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,10) 'MPI_BCAST', MyRank, MyError, string(1:Lstr)
        exit_flag=2
        RETURN
      END IF
!
!  Load appropriate data block into output array.
!
      DO i=Mstr,Mend
        Awrk(i)=A(i)
      END DO
!
!  Deallocate receive buffer.
!
      IF (allocated(Arecv).and.(MyRank.eq.MyMaster)) THEN
        deallocate (Arecv)
      END IF

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 47)
# endif

      RETURN
      END SUBROUTINE mp_scatter_state
!
      SUBROUTINE mp_dump (ng, tile, gtype,                              &
     &                    ILB, IUB, JLB, JUB, KLB, KUB, A, name)
!
!***********************************************************************
!                                                                      !
!  This routine is used to debug distributed-memory communications.    !
!  It writes field into an ASCII file for further post-processing.     !
!                                                                      !
!***********************************************************************
!

      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, gtype
      integer, intent(in) :: ILB, IUB, JLB, JUB, KLB, KUB

      real(r8), intent(in) :: A(ILB:IUB,JLB:JUB,KLB:KUB)

      character (len=*) :: name
!
!  Local variable declarations.
!
      common /counter/ nc
      integer :: nc

      logical, save :: first = .TRUE.

      integer :: Imin, Imax, Ioff, Jmin, Jmax, Joff
      integer :: unit

# include "set_bounds.h"
!
!------------------------------------------------------------------------
!  Write out requested field.
!------------------------------------------------------------------------
!
      IF (first) THEN
        nc=0
        first=.FALSE.
      END IF
      nc=nc+1
      IF (Master) THEN
        WRITE (10,'(a,i3.3,a,a)') 'file ', nc, ': ', TRIM(name)
        CALL my_flush (10)
      END IF
!
!  Write out field including ghost-points.
!
      Imin=0
      Imax=Lm(ng)+1
      IF (EWperiodic(ng)) THEN
        Ioff=3
      ELSE
        Ioff=1
      END IF

      Jmin=0
      Jmax=Mm(ng)+1
      IF (NSperiodic(ng)) THEN
        Joff=3
      ELSE
        Joff=1
      END IF

      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
        Imin=1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
        Jmin=1
      END IF

      unit=(MyRank+1)*1000+nc
      WRITE (unit,*) ILB, IUB, JLB, JUB, KLB, KUB,                      &
     &               Ioff, Joff, Imin, Imax, Jmin, Jmax,                &
     &               A(ILB:IUB,JLB:JUB,KLB:KUB)
      CALL my_flush (unit)
!
!  Write out non-overlapping field.
!
      Imin=IstrR
      Imax=IendR
      IF (EWperiodic(ng)) THEN
        Ioff=2
      ELSE
        Ioff=1
      END IF

      Jmin=JstrR
      Jmax=JendR
      IF (NSperiodic(ng)) THEN
        Joff=2
      ELSE
        Joff=1
      END IF

      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.u2dvar).or.(gtype.eq.u3dvar)) THEN
        Imin=Istr
        Ioff=Ioff-1
      END IF
      IF ((gtype.eq.p2dvar).or.(gtype.eq.p3dvar).or.                    &
     &    (gtype.eq.v2dvar).or.(gtype.eq.v3dvar)) THEN
        Jmin=Jstr
        Joff=Joff-1
      END IF

      unit=(MyRank+1)*10000+nc
      WRITE (unit,*) Imin, Imax, Jmin, Jmax, KLB, KUB,                  &
     &               Ioff, Joff, Imin, Imax, Jmin, Jmax,                &
     &               A(Imin:Imax,Jmin:Jmax,KLB:KUB)
      CALL my_flush (unit)

      RETURN
      END SUBROUTINE mp_dump
!
      SUBROUTINE mp_aggregate2d (ng, model, gtype,                      &
     &                           LBiT, UBiT, LBjT, UBjT,                &
     &                           LBiG, UBiG, LBjG, UBjG,                &
     &                           Atiled, Aglobal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 2D tiled, floating-point array from each    !
!  spawned node and stores it into 2D global array. If nesting, the    !
!  global array contains the contact points data.                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     gtype      C-grid type.                                          !
!     LBiT       Tiled  array, I-dimension Lower bound.                !
!     UBiT       Tiled  array, I-dimension Upper bound.                !
!     LBjT       Tiled  array, J-dimension Lower bound.                !
!     UBjT       Tiled  array, J-dimension Upper bound.                !
!     LBiG       Global array, I-dimension Lower bound.                !
!     UBiG       Global array, I-dimension Upper bound.                !
!     LBjG       Global array, J-dimension Lower bound.                !
!     UBjG       Global array, J-dimension Upper bound.                !
!     Atiled     2D tiled, floating-point array to process.            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aglobal    2D global array, all tiles are aggregated.            !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype
      integer, intent(in) :: LBiT, UBiT, LBjT, UBjT
      integer, intent(in) :: LBiG, UBiG, LBjG, UBjG

      real(r8), intent(in)  :: Atiled(LBiT:UBiT,LBjT:UBjT)
      real(r8), intent(out) :: Aglobal(LBiG:UBiG,LBjG:UBjG)
!
!  Local variable declarations.
!
      integer :: Lstr, MyError, MyType, Nnodes, Npts, Serror
      integer :: i, j, np, rank

      integer,  dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: my_bounds

      real(r8), dimension(TileSize(ng)) :: Asend
      real(r8), dimension(TileSize(ng)*                                 &
     &                    NtileI(ng)*NtileJ(ng)) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Number of nodes in the group.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
!
!  Set starting and ending indices to process including contact points
!  (if nesting) according to the staggered C-grid classification.
!
      MyType=ABS(gtype)

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
        CASE (r2dvar, r3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (u2dvar, u3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (v2dvar, v3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
      END SELECT
!
!  Determine the maximum number of points to process between all tiles.
!  In collective communications, the amount of data sent must be equal
!  to the amount of data received.
!
      Npts=0
      DO rank=0,Nnodes
        np=(my_bounds(2,rank)-my_bounds(1,rank)+1)*                     &
     &     (my_bounds(4,rank)-my_bounds(3,rank)+1)
        Npts=MAX(Npts, np)
      END DO

      IF (Npts.gt.TileSize(ng)) THEN
        IF (Master) THEN
          WRITE (stdout,10) ' TileSize = ', TileSize(ng), Npts
 10       FORMAT (/,' MP_AGGREGATE2D - communication buffer to small,', &
     &            a, 2i8)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Initialize local arrays to facilitate collective communicatios.
!  This also avoid denormalized values, which facilitates debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Pack tile data.
!-----------------------------------------------------------------------
!
      np=0
      DO j=my_bounds(3,MyRank),my_bounds(4,MyRank)
        DO i=my_bounds(1,MyRank),my_bounds(2,MyRank)
          np=np+1
          Asend(np)=Atiled(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Aggregate data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mpi_allgather (Asend, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT, &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,20) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 20     FORMAT (/,' MP_AGGREGATE2D - error during ',a,' call, Node = ', &
     &          i3.3,' Error = ',i3,/,18x,a)
        exit_flag=5
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Unpack data into a global 2D array.
!-----------------------------------------------------------------------
!
      DO rank=0,Nnodes
        np=rank*Npts
        DO j=my_bounds(3,rank),my_bounds(4,rank)
          DO i=my_bounds(1,rank),my_bounds(2,rank)
            np=np+1
            Aglobal(i,j)=Arecv(np)
          END DO
        END DO
      END DO

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

      RETURN
      END SUBROUTINE mp_aggregate2d
!
      SUBROUTINE mp_aggregate3d (ng, model, gtype,                      &
     &                           LBiT, UBiT, LBjT, UBjT,                &
     &                           LBiG, UBiG, LBjG, UBjG,                &
     &                           LBk,  UBk,                             &
     &                           Atiled, Aglobal)
!
!***********************************************************************
!                                                                      !
!  This routine collects a 3D tiled, floating-point array from each    !
!  spawned node and stores it into 3D global array. If nesting, the    !
!  global array contains the contact points data.                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     gtype      C-grid type.                                          !
!     LBiT       Tiled  array, I-dimension Lower bound.                !
!     UBiT       Tiled  array, I-dimension Upper bound.                !
!     LBjT       Tiled  array, J-dimension Lower bound.                !
!     UBjT       Tiled  array, J-dimension Upper bound.                !
!     LBkT       Tiled  array, K-dimension Lower bound.                !
!     UBkT       Tiled  array, K-dimension Upper bound.                !
!     LBiG       Global array, I-dimension Lower bound.                !
!     UBiG       Global array, I-dimension Upper bound.                !
!     LBjG       Global array, J-dimension Lower bound.                !
!     UBjG       Global array, J-dimension Upper bound.                !
!     LBkG       Global array, K-dimension Lower bound.                !
!     UBkG       Global array, K-dimension Upper bound.                !
!     Atiled     3D tiled, floating-point array to process.            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Aglobal    3D global array, all tiles are aggregated.            !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, gtype
      integer, intent(in) :: LBiT, UBiT, LBjT, UBjT
      integer, intent(in) :: LBiG, UBiG, LBjG, UBjG
      integer, intent(in) :: LBk,  UBk

      real(r8), intent(in)  :: Atiled(LBiT:UBiT,LBjT:UBjT,LBk:UBk)
      real(r8), intent(out) :: Aglobal(LBiG:UBiG,LBjG:UBjG,LBk:UBk)
!
!  Local variable declarations.
!
      integer :: Klen, Lstr, MyError, MyType, Nnodes, Npts, Serror
      integer :: i, j, k, np, rank

      integer,  dimension(4,0:NtileI(ng)*NtileJ(ng)-1) :: my_bounds

      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)) :: Asend

      real(r8), dimension(TileSize(ng)*(UBk-LBk+1)*                     &
     &                    NtileI(ng)*NtileJ(ng)) :: Arecv

      character (len=MPI_MAX_ERROR_STRING) :: string

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, model, 46)
# endif
!
!-----------------------------------------------------------------------
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!-----------------------------------------------------------------------
!
!  Number of nodes in the group.
!
      Nnodes=NtileI(ng)*NtileJ(ng)-1
!
!  Set starting and ending indices to process including contact points
!  (if nesting) according to the staggered C-grid classification.
!
      MyType=ABS(gtype)

      SELECT CASE (MyType)
        CASE (p2dvar, p3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
        CASE (r2dvar, r3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (u2dvar, u3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrP(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendP(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrT(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendT(rank)
          END DO
        CASE (v2dvar, v3dvar)
          DO rank=0,Nnodes
            my_bounds(1,rank)=BOUNDS(ng) % IstrT(rank)
            my_bounds(2,rank)=BOUNDS(ng) % IendT(rank)
            my_bounds(3,rank)=BOUNDS(ng) % JstrP(rank)
            my_bounds(4,rank)=BOUNDS(ng) % JendP(rank)
          END DO
      END SELECT
      Klen=UBk-LBk+1
!
!  Determine the maximum number of points to process between all tiles.
!  In collective communications, the amount of data sent must be equal
!  to the amount of data received.
!
      Npts=0
      DO rank=0,Nnodes
        np=(my_bounds(2,rank)-my_bounds(1,rank)+1)*                     &
     &     (my_bounds(4,rank)-my_bounds(3,rank)+1)*                     &
     &     Klen
        Npts=MAX(Npts, np)
      END DO

      IF (Npts.gt.TileSize(ng)*Klen) THEN
        IF (Master) THEN
          WRITE (stdout,10) ' TileSize = ', TileSize(ng)*Klen, Npts
 10       FORMAT (/,' MP_AGGREGATE3D - communication buffer to small,', &
     &            a, 2i8)
        END IF
        exit_flag=5
        RETURN
      END IF
!
!  Initialize local arrays to facilitate collective communicatios.
!  This also avoid denormalized values, which facilitates debugging.
!
      Asend=0.0_r8
      Arecv=0.0_r8
!
!-----------------------------------------------------------------------
!  Pack tile data.
!-----------------------------------------------------------------------
!
      np=0
      DO k=LBk,UBk
        DO j=my_bounds(3,MyRank),my_bounds(4,MyRank)
          DO i=my_bounds(1,MyRank),my_bounds(2,MyRank)
            np=np+1
            Asend(np)=Atiled(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Aggregate data from all nodes.
!-----------------------------------------------------------------------
!
      CALL mpi_allgather (Asend, Npts, MP_FLOAT, Arecv, Npts, MP_FLOAT, &
     &                    OCN_COMM_WORLD, MyError)
      IF (MyError.ne.MPI_SUCCESS) THEN
        CALL mpi_error_string (MyError, string, Lstr, Serror)
        Lstr=LEN_TRIM(string)
        WRITE (stdout,20) 'MPI_ALLGATHER', MyRank, MyError,             &
     &                    string(1:Lstr)
 20     FORMAT (/,' MP_AGGREGATE3D - error during ',a,' call, Node = ', &
     &          i3.3,' Error = ',i3,/,18x,a)
        exit_flag=5
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Unpack data into a global 2D array.
!-----------------------------------------------------------------------
!
      DO rank=0,Nnodes
        np=rank*Npts
        DO k=LBk,UBk
          DO j=my_bounds(3,rank),my_bounds(4,rank)
            DO i=my_bounds(1,rank),my_bounds(2,rank)
              np=np+1
              Aglobal(i,j,k)=Arecv(np)
            END DO
          END DO
        END DO
      END DO

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off time clocks.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, model, 46)
# endif

      RETURN
      END SUBROUTINE mp_aggregate3d
#endif
      END MODULE distribute_mod
