#include "cppdefs.h"

      MODULE mct_coupler_utils_mod
#ifdef COAWST_COUPLING

      implicit none

      PRIVATE
      PUBLIC  :: mct_getarg

      CONTAINS

!-----------------------------------------------------------------------
      SUBROUTINE mct_getarg (Iarg, Carg)
!-----------------------------------------------------------------------
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Iarg
      character (len=*), intent(inout) :: Carg
!
!  Local variable declarations.
!
      integer :: Lstr, ierror

# ifdef CRAY
      CALL pxfgetarg (Iarg, Carg, Lstr, ierror)
# else
      CALL getarg (Iarg, Carg)
# endif

      RETURN
      END SUBROUTINE mct_getarg
#endif

      END MODULE mct_coupler_utils_mod
