#include "cppdefs.h"
      MODULE mod_coupler_iounits
#if defined COAWST_COUPLING && defined MCT_LIB
!
!svn $Id: mod_coupler_iounits.F jcwarner $
!=======================================================================
!                                                                      !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  Aname       Atmosphere model stadard input file name.               !
!  IWBNDname   Input boundary data file name for InWave model          !
!  IWSWNname   Input spectral SWAN data file name for InWave model     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
# ifdef SWAN_MODEL
      USE swan_iounits
# endif
# ifdef ROMS_MODEL
      USE mod_param
# endif
      implicit none

# ifndef ROMS_MODEL
      integer, parameter :: stdinp = 5
      integer, parameter :: stdout = 6
      integer :: ioerror
# endif

!# if defined MCT_INTERP_OC2WV || defined MCT_INTERP_OC2AT || \
!     defined MCT_INTERP_WV2AT
! this flag is temporary for SCRIP option.
      integer :: scrip_opt
      character (len=80) :: SCRIPname
!# endif
# ifdef ROMS_COUPLING
#  ifdef MCT_INTERP_OC2WV
      character (len=80), dimension(:,:), pointer :: O2Wname
      character (len=80), dimension(:,:), pointer :: W2Oname
#  endif
#  ifdef MCT_INTERP_OC2AT
      character (len=80), dimension(:,:), pointer :: O2Aname
      character (len=80), dimension(:,:), pointer :: A2Oname
#  endif
# endif
# ifdef SWAN_COUPLING
#  ifdef MCT_INTERP_WV2AT
      character (len=80), dimension(:,:), pointer :: W2Aname
      character (len=80), dimension(:,:), pointer :: A2Wname
#  endif
# endif
# ifdef WRF_MODEL
      character (len=80) :: Aname
# endif
# ifdef INWAVE_MODEL
!     character (len=80), dimension(Ngridss) :: IWINIname
!     character (len=80), dimension(Ngridss) :: IWBNDname
!     character (len=80), dimension(Ngridss) :: IWSWNname
# endif

      CONTAINS

      SUBROUTINE allocate_coupler_iounits
!=======================================================================
!                                                                      !
!  This routine initialize all the coupler io vars.                    !
!                                                                      !
!=======================================================================

      character (len=1 ), parameter :: blank = ' '
      integer :: i, io, ia, iw

# ifdef ROMS_COUPLING
#  ifdef MCT_INTERP_OC2WV
      allocate (O2Wname(Nocn_grids,Nwav_grids))
      allocate (W2Oname(Nwav_grids,Nocn_grids))
      DO io=1,Nocn_grids
        DO iw=1,Nwav_grids
          DO i=1,LEN(O2Wname(io,iw))
            O2Wname(io,iw)(i:i)=blank
            W2Oname(iw,io)(i:i)=blank
          END DO
        END DO
      END DO
#  endif

#  ifdef MCT_INTERP_OC2AT
      allocate (O2Aname(Nocn_grids,Natm_grids))
      allocate (A2Oname(Natm_grids,Nocn_grids))
      DO ia=1,Natm_grids
        DO io=1,Nocn_grids
          DO i=1,LEN(O2Aname(io,ia))
            O2Aname(io,ia)(i:i)=blank
            A2Oname(ia,io)(i:i)=blank
          END DO
        END DO
      END DO
#  endif

# endif
# ifdef SWAN_COUPLING
#  ifdef MCT_INTERP_WV2AT
      allocate (W2Aname(Nwav_grids,Natm_grids))
      allocate (A2Wname(Natm_grids,Nwav_grids))
      DO ia=1,Natm_grids
        DO iw=1,Nwav_grids
          DO i=1,LEN(A2Wname(iw,ia))
            A2Wname(ia,iw)(i:i)=blank
            W2Aname(iw,ia)(i:i)=blank
          END DO
        END DO
      END DO
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_coupler_iounits
#endif

      END MODULE mod_coupler_iounits
