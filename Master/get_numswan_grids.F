#include "cppdefs.h"

      SUBROUTINE get_numswan_grids(tempname)
#if defined SWAN_MODEL
!
!=======================================================================
!                                                                      !
!  This routine reads the swan model input and gets the number         !
!  of swan grids.                                                      !
!                                                                      !
!=======================================================================
!
      USE swan_iounits
      USE M_COUPLING
      USE mod_coupler_kinds

      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
      character (len=160), intent(in) :: tempname
!
!  Local variable declarations.
!
      integer :: i, inp, out, status
      integer :: MyRank, MyError, MyMaster

      character (len=160) :: line
      character (len=160) :: aline
!
      inp=2
      out=6
      MyMaster=0
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Read SWAN input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      OPEN (inp, FILE=TRIM(tempname), FORM='formatted', STATUS='old',  &
     &      ERR=110)
      GO TO 130
 110  WRITE (out,120) tempname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(tempname)
!     exit_flag=4
      RETURN
 120  FORMAT (/,' get num swan grids - Unable to open swan input file.', &
     &        /,a80)
 130  CONTINUE
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=116,END=140) line
        aline=ADJUSTL(line)
        IF(aline(1:7).eq.'NSGRIDS') THEN
          read(aline(8:12),'(i5)') NUM_SGRIDS
# ifdef NESTING
          IF (NUM_SGRIDS.eq.1) THEN
            WRITE(out,'(a40)') 'Num swan grids should be > 1 for nesting'
          END IF
# else
          IF (NUM_SGRIDS.ne.1) THEN
            WRITE(out,'(a40)') 'Num swan grids should = 1 for no nesting'
          END IF
# endif
        END IF
      END DO
 116  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 140  CLOSE (inp)
!

  60  FORMAT (/,'read model inputs - Error while processing line: ',/,a)

#endif
      RETURN
      END SUBROUTINE get_numswan_grids
