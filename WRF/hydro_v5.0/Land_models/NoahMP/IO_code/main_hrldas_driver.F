#include "wrfhydro_cpp.h"
!  Program Name:
!  Author(s)/Contact(s):
!  Abstract:
!  History Log:
! 
!  Usage:
!  Parameters: <Specify typical arguments passed>
!  Input Files:
!        <list file names and briefly describe the data they include>
!  Output Files:
!        <list file names and briefly describe the information they include>
! 
!  Condition codes:
!        <list exit condition or error codes returned >
!        If appropriate, descriptive troubleshooting instructions or
!        likely causes for failures could be mentioned here with the
!        appropriate error code
! 
!  User controllable options: <if applicable>

#if !defined COAWST_COUPLING
program Noah_hrldas_driver
! this is the main program to drive HRLDAS-Noah, HRLDAS-NoahMP, and other Land models.

# ifdef Noah1d
! this is used to drive Noah1d
   use module_noah1d_hrldas_driver, only: land_driver_ini, land_driver_exe
# else
! this is used to drive NoahMP
   use module_noahmp_hrldas_driver, only: land_driver_ini, land_driver_exe
# endif

   implicit none
   integer :: ITIME, NTIME

   call land_driver_ini(NTIME)
#ifdef WRF_HYDRO   
   do ITIME = 1, NTIME
#else
   do ITIME = 0, NTIME
#endif
       call land_driver_exe(ITIME)
   end do

# ifdef WRF_HYDRO
   write(6,*) "NTIME = ",NTIME, " ITIME =", ITIME
   call hydro_finish() 
# endif

END 

#else

      module main_hrldas_driver

! this is used to drive NoahMP
      use module_noahmp_hrldas_driver, only: land_driver_ini, land_driver_exe
      use mct_coupler_params
      use hydro_coupler

      implicit none
      include "mpif.h"
      integer :: ITIME, NTIME

      CONTAINS

      subroutine coawst_hydro_init(MyComm)
      integer, intent(in) :: MyComm
      write(*,*) 'calling land_driver_ini my comm is ', MyComm
      HYDRO_COMM_WORLD=MyComm
      call land_driver_ini(NTIME)
      call initialize_hyd_coupling(1)
      call initialize_hyd_routers
      end subroutine coawst_hydro_init
   
      subroutine coawst_hydro_run
      integer :: ih, io
      write(*,*) 'calling land driver run ', ntime

      do ITIME = 1, NTIME
!       write(*,*) 'looking to call hyd2ocn coupling ', itime, nhydfocn(1,1)
        IF (MOD(ITIME-1,nHYDFOCN(1,1)).eq.0) THEN
          DO io=1,Nocn_grids
            DO ih=1,Nhyd_grids
              call hyd2ocn_coupling(ih,io)
            END DO
          END DO
          DO io=1,Nocn_grids
            DO ih=1,Nhyd_grids
              call hydfocn_coupling(io,ih)
            END DO
          END DO
        END IF
        call land_driver_exe(ITIME)
      end do
      end subroutine coawst_hydro_run

      subroutine coawst_hydro_finalize
      write(6,*) "NTIME = ",NTIME, " ITIME =", ITIME
      write(*,*) 'calling land driver end ', ntime
      call finalize_hyd_coupling(1)
!     call hydro_finish()
      end subroutine coawst_hydro_finalize

      end module main_hrldas_driver

#endif

