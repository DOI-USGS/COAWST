#include "cppdefs.h"
      MODULE mod_inwave_params
!
!svn $Id: mod_inwave.F 790 2008-10-14 19:32:07Z jcwarner $
!================================================== John C. Warner =====
!                                                                      !
!  ND        number of directional dimensions                          !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE mod_param

      implicit none

      integer, dimension(Ngrids) :: ND

      END MODULE mod_inwave_params
