! ###########################################################################################
! CI test driver for MYNN-EDMF scheme
! ###########################################################################################
program driver

use module_bl_mynnedmf_driver
use module_bl_mynnedmf_wrf_tests

implicit none

call wrf_test(case='clr',bl_mynn_closure=2.6,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2,               &
      bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1,   &
      bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,&
      tke_budget=1)

call wrf_test(case='clr',bl_mynn_closure=3.0,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2,               &
      bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1,   &
      bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,&
      tke_budget=1)

call wrf_test(case='marine_StCu',bl_mynn_closure=3.0,bl_mynn_cloudpdf=2,bl_mynn_mixlength=2,               &
      bl_mynn_edmf=1,bl_mynn_edmf_dd=1,bl_mynn_edmf_mom=1,bl_mynn_edmf_tke=1,bl_mynn_cloudmix=1,           &
      bl_mynn_mixqt=0,bl_mynn_mixscalars=0,bl_mynn_mixaerosols=0,bl_mynn_mixnumcon=0,bl_mynn_ess=1,        &
      tke_budget=1)

end program driver
