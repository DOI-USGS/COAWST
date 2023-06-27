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
      integer  :: idACtp                    ! absolute peak period
      integer  :: idACkc                    ! wave number
      integer  :: idACag                    ! energy angle direction
      integer  :: idACac                    ! energy angle direction centered
      integer  :: idACen                    ! wave action 
      integer  :: idACkl                    ! wave 3D wave length
      integer  :: idACcx                    ! wave action celerity in xi  
      integer  :: idACcy                    ! wave action celerity in eta
      integer  :: idACct                    ! wave action celerity in theta dir
      integer  :: idACbz                    ! wave action bound wave
      integer  :: NDInw                     ! number of inwave diags

#ifdef DIAGNOSTICS_INWAVE
!
!  InWave diagnostics.
!
      integer  :: NDInw                     ! number of inwave diags
      integer  :: iIrate                    ! temporal rate
      integer  :: iIxadv                    ! x advection
      integer  :: iIyadv                    ! y advection
      integer  :: iItadv                    ! theta advection
      integer  :: iIbrak                    ! breaking
      integer  :: iIfric                    ! friction

      integer  :: idIrate                    ! temporal rate
      integer  :: idIxadv                    ! x advection
      integer  :: idIyadv                    ! y advection
      integer  :: idItadv                    ! theta advection
      integer  :: idIbrak                    ! breaking
      integer  :: idIfric                    ! friction

#endif


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
      integer :: ic

#ifdef DIAGNOSTICS_INWAVE
!
!  Indices for inwave diagnostic variables in the Dia_ terms.
!
      iIrate=1
      iIxadv=2
      iIyadv=3
      iItadv=4
      iIbrak=5
      ic=5
# if defined WDISS_FRIC
      ic=ic+1
      iIfric=ic
# endif
      NDInw=ic

!  Indices for inwave diagnostic vars for varinfo.dat
!
      idIrate=1
      idIxadv=2
      idIyadv=3
      idItadv=4
      idIbrak=5
      ic=5
# if defined WDISS_FRIC
      ic=ic+1
      idIfric=ic
# endif

#else
      NDInw=0          ! No inwave diagnostics
#endif

      RETURN
      END SUBROUTINE allocate_inwave_params

      END MODULE mod_inwave_params
