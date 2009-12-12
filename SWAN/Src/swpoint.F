#include "cppdefs.h"
      MODULE SWPOINT_MOD

#if defined MODEL_COUPLING
!
!svn $Id: swpoint_mod.F 756 2008-09-14 20:18:28Z jcwarner $
!==================================================== John C. Warner ===
!                                                                      !
!                                                                      !
!=======================================================================
!
!
      USE TIMECOMM
      USE OCPCOMM2
      USE OCPCOMM4
      USE SWCOMM1
      USE SWCOMM2
      USE SWCOMM3
      USE SWCOMM4
      USE OUTP_DATA
      USE M_GENARR
      USE M_PARALL
      USE M_PARALL2
      USE M_BNDSPEC
      USE M_CVMESH

      implicit none
!
      PRIVATE

      PUBLIC :: ALLOCATE_SWAN_ARRAYS
      PUBLIC :: INIT_POINTERS
      PUBLIC :: DEALLOCATE_SWAN_ARRAYS

      include 'mpif.h'
!
!  Declarations.
!
      CONTAINS

      SUBROUTINE ALLOCATE_SWAN_ARRAYS (ng, Numgrids)
!
!=======================================================================
!                                                                      !
!  Initialize global arrays.                                           !
!                                                                      !
!=======================================================================
!
      INTEGER, INTENT(IN) :: ng, Numgrids
!
!
!  Allocate Type arrays
!
      IF (ng.eq.1) THEN
        ALLOCATE (M_GENARR_MOD(Numgrids))
        ALLOCATE (PARALL_MOD(Numgrids))
        ALLOCATE (PARALL2_MOD(Numgrids))
        ALLOCATE (BGPDATZ_MOD(Numgrids))
        ALLOCATE (OPSDATZ_MOD(Numgrids))
        ALLOCATE (ORQDATZ_MOD(Numgrids))
!
!  Allocate pointers
!
!  M_BNDSPEC
!
        ALLOCATE (LBGP_G(Numgrids))
!       LBGP_G(ng)=.FALSE.
!
!  M_PARALL
!
        IF (.not.allocated(MCGRDGL_G)) ALLOCATE (MCGRDGL_G(Numgrids))
        IF (.not.allocated(NGRBGL_G)) ALLOCATE (NGRBGL_G(Numgrids))
        IF (.not.allocated(MXCGL_G)) ALLOCATE (MXCGL_G(Numgrids))
        IF (.not.allocated(MYCGL_G)) ALLOCATE (MYCGL_G(Numgrids))
        IF (.not.allocated(MXF_G)) ALLOCATE (MXF_G(Numgrids))
        IF (.not.allocated(MXL_G)) ALLOCATE (MXL_G(Numgrids))
        IF (.not.allocated(MYF_G)) ALLOCATE (MYF_G(Numgrids))
        IF (.not.allocated(MYL_G)) ALLOCATE (MYL_G(Numgrids))
        IF (.not.allocated(NBGGL_G)) ALLOCATE (NBGGL_G(Numgrids))
        IF (.not.allocated(XCLMAX_G)) ALLOCATE (XCLMAX_G(Numgrids))
        IF (.not.allocated(YCLMAX_G)) ALLOCATE (YCLMAX_G(Numgrids))
        IF (.not.allocated(XCLMIN_G)) ALLOCATE (XCLMIN_G(Numgrids))
        IF (.not.allocated(YCLMIN_G)) ALLOCATE (YCLMIN_G(Numgrids))
        IF (.not.allocated(LMXF_G)) ALLOCATE (LMXF_G(Numgrids))
        IF (.not.allocated(LMXL_G)) ALLOCATE (LMXL_G(Numgrids))
        IF (.not.allocated(LMYF_G)) ALLOCATE (LMYF_G(Numgrids))
        IF (.not.allocated(LMYL_G)) ALLOCATE (LMYL_G(Numgrids))
!
!  OCPCOMM4
!
        IF (.not.allocated(FUNHI_G)) ALLOCATE (FUNHI_G(Numgrids))
        IF (.not.allocated(FUNLO_G)) ALLOCATE (FUNLO_G(Numgrids))
        IF (.not.allocated(HIOPEN_G)) ALLOCATE (HIOPEN_G(Numgrids))
        IF (.not.allocated(PRINTF_G)) ALLOCATE (PRINTF_G(Numgrids))
!
!  SWANCOMM3
!
        IF (.not.allocated(MXC_G)) ALLOCATE (MXC_G(Numgrids))
        IF (.not.allocated(MYC_G)) ALLOCATE (MYC_G(Numgrids))
        IF (.not.allocated(MCGRD_G)) ALLOCATE (MCGRD_G(Numgrids))
        IF (.not.allocated(NGRBND_G)) ALLOCATE (NGRBND_G(Numgrids))
        IF (.not.allocated(MSC_G)) ALLOCATE (MSC_G(Numgrids))
        IF (.not.allocated(MDC_G)) ALLOCATE (MDC_G(Numgrids))
        IF (.not.allocated(MTC_G)) ALLOCATE (MTC_G(Numgrids))
        IF (.not.allocated(XCGMIN_G)) ALLOCATE (XCGMIN_G(Numgrids))
        IF (.not.allocated(XCGMAX_G)) ALLOCATE (XCGMAX_G(Numgrids))
        IF (.not.allocated(YCGMIN_G)) ALLOCATE (YCGMIN_G(Numgrids))
        IF (.not.allocated(YCGMAX_G)) ALLOCATE (YCGMAX_G(Numgrids))
!
!  TIMECOMM
!
        IF (.not.allocated(TINIC_G)) ALLOCATE (TINIC_G(Numgrids))
        IF (.not.allocated(DT_G)) ALLOCATE (DT_G(Numgrids))
        IF (.not.allocated(TFINC_G)) ALLOCATE (TFINC_G(Numgrids))
        IF (.not.allocated(TIMCO_G)) ALLOCATE (TIMCO_G(Numgrids))
!
!  SWCOMM1
!
        IF (.not.allocated(CHTIME_G)) ALLOCATE (CHTIME_G(Numgrids))
!
! OUTP_DATA
!
        ALLOCATE (NREOQ_G(Numgrids))
!       NREOQ_G(ng)=0
        ALLOCATE (LOPS_G(Numgrids))
!       LOPS_G(ng)=.FALSE.
        ALLOCATE (LORQ_G(Numgrids))
!       LORQ_G(ng)=.FALSE.
        ALLOCATE (OUTP_FILES_G(Numgrids,1:MAX_OUTP_REQ))
!
!       ALLOCATE (ORQDATZ_MOD(ng)%FORQ)
!
!  SWANCOMM2
!
        IF (.not.allocated(CVLEFT_G)) ALLOCATE (CVLEFT_G(Numgrids))
        IF (.not.allocated(LXOFFS_G)) ALLOCATE (LXOFFS_G(Numgrids))
        IF (.not.allocated(VARFR_G)) ALLOCATE (VARFR_G(Numgrids))
        IF (.not.allocated(VARWI_G)) ALLOCATE (VARWI_G(Numgrids))
        IF (.not.allocated(VARWLV_G)) ALLOCATE (VARWLV_G(Numgrids))
        IF (.not.allocated(DYNDEP_G)) ALLOCATE (DYNDEP_G(Numgrids))
        IF (.not.allocated(VARAST_G)) ALLOCATE (VARAST_G(Numgrids))
        IF (.not.allocated(ICOND_G)) ALLOCATE (ICOND_G(Numgrids))
        IF (.not.allocated(NWAMN_G)) ALLOCATE (NWAMN_G(Numgrids))
        IF (.not.allocated(OPTG_G)) ALLOCATE (OPTG_G(Numgrids))
        IF (.not.allocated(NBFILS_G)) ALLOCATE (NBFILS_G(Numgrids))
        IF (.not.allocated(NBSPEC_G)) ALLOCATE (NBSPEC_G(Numgrids))
        IF (.not.allocated(NBGRPT_G)) ALLOCATE (NBGRPT_G(Numgrids))

        IF (.not.allocated(MXG_G)) ALLOCATE (MXG_G(Numgrids,NUMGRD))
        IF (.not.allocated(MYG_G)) ALLOCATE (MYG_G(Numgrids,NUMGRD))
        IF (.not.allocated(IGTYPE_G))                                   &
     &           ALLOCATE (IGTYPE_G(Numgrids,NUMGRD))
        IF (.not.allocated(LEDS_G)) ALLOCATE (LEDS_G(Numgrids,NUMGRD))

        IF (.not.allocated(COSVC_G)) ALLOCATE (COSVC_G(Numgrids))
        IF (.not.allocated(COSWC_G)) ALLOCATE (COSWC_G(Numgrids))
        IF (.not.allocated(RDTIM_G)) ALLOCATE (RDTIM_G(Numgrids))
        IF (.not.allocated(SINVC_G)) ALLOCATE (SINVC_G(Numgrids))
        IF (.not.allocated(SINWC_G)) ALLOCATE (SINWC_G(Numgrids))
        IF (.not.allocated(XOFFS_G)) ALLOCATE (XOFFS_G(Numgrids))
        IF (.not.allocated(YOFFS_G)) ALLOCATE (YOFFS_G(Numgrids))

        IF (.not.allocated(ALPG_G))  ALLOCATE (ALPG_G(Numgrids,NUMGRD))
        IF (.not.allocated(COSPG_G)) ALLOCATE (COSPG_G(Numgrids,NUMGRD))
        IF (.not.allocated(DXG_G))   ALLOCATE (DXG_G(Numgrids,NUMGRD))
        IF (.not.allocated(DYG_G))   ALLOCATE (DYG_G(Numgrids,NUMGRD))
        IF (.not.allocated(EXCFLD_G))                                   &
     &           ALLOCATE (EXCFLD_G(Numgrids,NUMGRD))
        IF (.not.allocated(SINPG_G)) ALLOCATE (SINPG_G(Numgrids,NUMGRD))
        IF (.not.allocated(STAGX_G)) ALLOCATE (STAGX_G(Numgrids,NUMGRD))
        IF (.not.allocated(STAGY_G)) ALLOCATE (STAGY_G(Numgrids,NUMGRD))
        IF (.not.allocated(XPG_G))   ALLOCATE (XPG_G(Numgrids,NUMGRD))
        IF (.not.allocated(YPG_G))   ALLOCATE (YPG_G(Numgrids,NUMGRD))

        IF (.not.allocated(IFLDYN_G))                                   &
     &           ALLOCATE (IFLDYN_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLIDL_G))                                   &
     &           ALLOCATE (IFLIDL_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLIFM_G))                                   &
     &           ALLOCATE (IFLIFM_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLNHF_G))                                   &
     &           ALLOCATE (IFLNHF_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLNHD_G))                                   &
     &           ALLOCATE (IFLNHD_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLNDS_G))                                   &
     &           ALLOCATE (IFLNDS_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLNDF_G))                                   &
     &           ALLOCATE (IFLNDF_G(Numgrids,NUMGRD))

        IF (.not.allocated(IFLBEG_G))                                   &
     &           ALLOCATE (IFLBEG_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLINT_G))                                   &
     &           ALLOCATE (IFLINT_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLEND_G))                                   &
     &           ALLOCATE (IFLEND_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLTIM_G))                                   &
     &           ALLOCATE (IFLTIM_G(Numgrids,NUMGRD))
        IF (.not.allocated(IFLFAC_G))                                   &
     &           ALLOCATE (IFLFAC_G(Numgrids,NUMGRD))
!
!     MODULE M_CVMESH
!
        ALLOCATE (XCSAVE_G(Numgrids))
        ALLOCATE (YCSAVE_G(Numgrids))
        ALLOCATE (MXITNR_G(Numgrids))
      END IF

      END SUBROUTINE ALLOCATE_SWAN_ARRAYS

      SUBROUTINE INIT_POINTERS (ng)
!
!=======================================================================
!                                                                      !
!  Establish swan pointers to the global arrays.                       !
!                                                                      !
!=======================================================================
!
      INTEGER, INTENT(IN) :: ng
!
!  Allocate pointers
!
!  M_BNDSPEC
!
!     LBGP=>LBGP_G(ng)
!
!  M_PARALL
!
      MCGRDGL=>MCGRDGL_G(ng)
      NGRBGL=>NGRBGL_G(ng)
      MXCGL=>MXCGL_G(ng)
      MYCGL=>MYCGL_G(ng)
      MXF=>MXF_G(ng)
      MXL=>MXL_G(ng)
      MYF=>MYF_G(ng)
      MYL=>MYL_G(ng)
      NBGGL=>NBGGL_G(ng)
      XCLMAX=>XCLMAX_G(ng)
      YCLMAX=>YCLMAX_G(ng)
      XCLMIN=>XCLMIN_G(ng)
      YCLMIN=>YCLMIN_G(ng)
      LMXF=>LMXF_G(ng)
      LMXL=>LMXL_G(ng)
      LMYF=>LMYF_G(ng)
      LMYL=>LMYL_G(ng)
!
!  OCPCOMM4
!
      FUNHI=>FUNHI_G(ng)
      FUNLO=>FUNLO_G(ng)
      HIOPEN=>HIOPEN_G(ng)
      PRINTF=>PRINTF_G(ng)
!
!  SWANCOMM3
!
      MXC=>MXC_G(ng)
      MYC=>MYC_G(ng)
      MCGRD=>MCGRD_G(ng)
      NGRBND=>NGRBND_G(ng)
      MSC=>MSC_G(ng)
      MDC=>MDC_G(ng)
      MTC=>MTC_G(ng)
      XCGMIN=>XCGMIN_G(ng)
      XCGMAX=>XCGMAX_G(ng)
      YCGMIN=>YCGMIN_G(ng)
      YCGMAX=>YCGMAX_G(ng)
!
!  TIMECOMM
!
      TINIC=>TINIC_G(ng)
      DT=>DT_G(ng)
      TFINC=>TFINC_G(ng)
      TIMCO=>TIMCO_G(ng)
!
!  SWCOMM1
!
      CHTIME=>CHTIME_G(ng)
!
! OUTP_DATA
!
      NREOQ=>NREOQ_G(ng)
      LOPS=>LOPS_G(ng)
!     LORQ=>LORQ_G(ng)
      OUTP_FILES=>OUTP_FILES_G(ng,:)
!
!  SWANCOMM2
!
      CVLEFT=>CVLEFT_G(ng)
      LXOFFS=>LXOFFS_G(ng)
      VARFR=>VARFR_G(ng)
      VARWI=>VARWI_G(ng)
      VARWLV=>VARWLV_G(ng)
      DYNDEP=>DYNDEP_G(ng)
      VARAST=>VARAST_G(ng)
      ICOND=>ICOND_G(ng)
      NWAMN=>NWAMN_G(ng)
      OPTG=>OPTG_G(ng)
      NBFILS=>NBFILS_G(ng)
      NBSPEC=>NBSPEC_G(ng)
      NBGRPT=>NBGRPT_G(ng)

      MXG=>MXG_G(ng,:)
      MYG=>MYG_G(ng,:)
      IGTYPE=>IGTYPE_G(ng,:)
      LEDS=>LEDS_G(ng,:)

      COSVC=>COSVC_G(ng)
      COSWC=>COSWC_G(ng)
      RDTIM=>RDTIM_G(ng)
      SINVC=>SINVC_G(ng)
      SINWC=>SINWC_G(ng)
      XOFFS=>XOFFS_G(ng)
      YOFFS=>YOFFS_G(ng)

      ALPG=>ALPG_G(ng,:)
      COSPG=>COSPG_G(ng,:)
      DXG=>DXG_G(ng,:)
      DYG=>DYG_G(ng,:)
      EXCFLD=>EXCFLD_G(ng,:)
      SINPG=>SINPG_G(ng,:)
      STAGX=>STAGX_G(ng,:)
      STAGY=>STAGY_G(ng,:)
      XPG=>XPG_G(ng,:)
      YPG=>YPG_G(ng,:)

      IFLDYN=>IFLDYN_G(ng,:)
      IFLIDL=>IFLIDL_G(ng,:)
      IFLIFM=>IFLIFM_G(ng,:)
      IFLNHF=>IFLNHF_G(ng,:)
      IFLNHD=>IFLNHD_G(ng,:)
      IFLNDS=>IFLNDS_G(ng,:)
      IFLNDF=>IFLNDF_G(ng,:)

      IFLBEG=>IFLBEG_G(ng,:)
      IFLINT=>IFLINT_G(ng,:)
      IFLEND=>IFLEND_G(ng,:)
      IFLTIM=>IFLTIM_G(ng,:)
      IFLFAC=>IFLFAC_G(ng,:)

      END SUBROUTINE INIT_POINTERS

      SUBROUTINE DEALLOCATE_SWAN_ARRAYS (ng, Numgrids)
!
!=======================================================================
!                                                                      !
!  Initialize global arrays.                                           !
!                                                                      !
!=======================================================================
!
      INTEGER, INTENT(IN) :: ng, Numgrids
!
!
!  DEAllocate Type arrays
!
      IF (ng.eq.Numgrids) THEN
        DEALLOCATE (M_GENARR_MOD)
        DEALLOCATE (PARALL_MOD)
        DEALLOCATE (PARALL2_MOD)
        DEALLOCATE (BGPDATZ_MOD)
        DEALLOCATE (OPSDATZ_MOD)
        DEALLOCATE (ORQDATZ_MOD)
!
!  DEAllocate pointers
!
!  M_BNDSPEC
!
        DEALLOCATE (LBGP_G)
!       LBGP_G(ng)=.FALSE.
!
!  M_PARALL
!
        DEALLOCATE (MCGRDGL_G)
        DEALLOCATE (NGRBGL_G)
        DEALLOCATE (MXCGL_G)
        DEALLOCATE (MYCGL_G)
        DEALLOCATE (MXF_G)
        DEALLOCATE (MXL_G)
        DEALLOCATE (MYF_G)
        DEALLOCATE (MYL_G)
        DEALLOCATE (NBGGL_G)
        DEALLOCATE (XCLMAX_G)
        DEALLOCATE (YCLMAX_G)
        DEALLOCATE (XCLMIN_G)
        DEALLOCATE (YCLMIN_G)
        DEALLOCATE (LMXF_G)
        DEALLOCATE (LMXL_G)
        DEALLOCATE (LMYF_G)
        DEALLOCATE (LMYL_G)
!
!  OCPCOMM4
!
        DEALLOCATE (FUNHI_G)
        DEALLOCATE (FUNLO_G)
        DEALLOCATE (HIOPEN_G)
        DEALLOCATE (PRINTF_G)
!
!  SWANCOMM3
!
        DEALLOCATE (MXC_G)
        DEALLOCATE (MYC_G)
        DEALLOCATE (MCGRD_G)
        DEALLOCATE (NGRBND_G)
        DEALLOCATE (MSC_G)
        DEALLOCATE (MDC_G)
        DEALLOCATE (MTC_G)
        DEALLOCATE (XCGMIN_G)
        DEALLOCATE (XCGMAX_G)
        DEALLOCATE (YCGMIN_G)
        DEALLOCATE (YCGMAX_G)
!
!  TIMECOMM
!
        DEALLOCATE (TINIC_G)
        DEALLOCATE (DT_G)
        DEALLOCATE (TFINC_G)
        DEALLOCATE (TIMCO_G)
!
!  SWCOMM1
!
        DEALLOCATE (CHTIME_G)
!
! OUTP_DATA
!
        DEALLOCATE (NREOQ_G)
!       NREOQ_G(ng)=0
        DEALLOCATE (LOPS_G)
!       LOPS_G(ng)=.FALSE.
        DEALLOCATE (LORQ_G)
!       LORQ_G(ng)=.FALSE.
        DEALLOCATE (OUTP_FILES_G)
!
!       ALLOCATE (ORQDATZ_MOD(ng)%FORQ)
!
!  SWANCOMM2
!
        DEALLOCATE (CVLEFT_G)
        DEALLOCATE (LXOFFS_G)
        DEALLOCATE (VARFR_G)
        DEALLOCATE (VARWI_G)
        DEALLOCATE (VARWLV_G)
        DEALLOCATE (DYNDEP_G)
        DEALLOCATE (VARAST_G)
        DEALLOCATE (ICOND_G)
        DEALLOCATE (NWAMN_G)
        DEALLOCATE (OPTG_G)
        DEALLOCATE (NBFILS_G)
        DEALLOCATE (NBSPEC_G)
        DEALLOCATE (NBGRPT_G)

        DEALLOCATE (MXG_G)
        DEALLOCATE (MYG_G)
        DEALLOCATE (IGTYPE_G)
        DEALLOCATE (LEDS_G)

        DEALLOCATE (COSVC_G)
        DEALLOCATE (COSWC_G)
        DEALLOCATE (RDTIM_G)
        DEALLOCATE (SINVC_G)
        DEALLOCATE (SINWC_G)
        DEALLOCATE (XOFFS_G)
        DEALLOCATE (YOFFS_G)

        DEALLOCATE (ALPG_G)
        DEALLOCATE (COSPG_G)
        DEALLOCATE (DXG_G)
        DEALLOCATE (DYG_G)
        DEALLOCATE (EXCFLD_G)
        DEALLOCATE (SINPG_G)
        DEALLOCATE (STAGX_G)
        DEALLOCATE (STAGY_G)
        DEALLOCATE (XPG_G)
        DEALLOCATE (YPG_G)
!
!     MODULE M_CVMESH
!
        DEALLOCATE (XCSAVE_G)
        DEALLOCATE (YCSAVE_G)
        DEALLOCATE (MXITNR_G)
      END IF

      END SUBROUTINE DEALLOCATE_SWAN_ARRAYS

#endif
      END MODULE SWPOINT_MOD
