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

!     These indices are for the boundaries.
      integer  :: isAC3d                    ! 3D Action density
      integer  :: isCT3d                    ! 3D wave theta celerity
      integer  :: isCX3d                    ! 3D wave x-dir celerity
      integer  :: isCY3d                    ! 3D wave y-dir celerity

!     These indices are for the I/O.
      integer  :: idACen                    ! wave action 
      integer  :: idACcx                    ! wave action celerity in xi  
      integer  :: idACcy                    ! wave action celerity in eta
      integer  :: idACct                    ! wave action celerity in theta dir
      integer  :: idACtp                    ! absolute peak period

!  These integers are assigned values in mod_ncparam.F

      integer :: ND

      CONTAINS

      SUBROUTINE allocate_inwave_params
!
!=======================================================================
!                                                                      !
!  This routine allocates several variables in the module that depend  !
!  on the number of nested grids.                                      !
!                                                                      !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Allocate dimension parameters.
!-----------------------------------------------------------------------
!
!mai      allocate ( ND )

      RETURN
      END SUBROUTINE allocate_inwave_params


      END MODULE mod_inwave_params
